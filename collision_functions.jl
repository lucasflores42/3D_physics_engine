# -----------------------------------------------------------------------------
#                           Parameters
# -----------------------------------------------------------------------------
const restitution = [0.5, 0.5, 0.3]
const restitution_angular = 0.5
const collision_min_distance = grid_size #* sqrt(2)
const max_velocity = 50.0
const max_angular_velocity = 20.0
const friction_coef = 0.3

include("rigidbody_functions.jl")
include("material_functions.jl")

# -----------------------------------------------------------------------------
#                           Calculate all collisions
# -----------------------------------------------------------------------------
function collision_physics!(particles, rigidbodies, powder, liquid, gas, id_grid, cell_of_particle)

    n_particles = length(particles)
    pos_correction = [@SVector zeros(3) for _ in 1:n_particles]
    vel_correction = [@SVector zeros(3) for _ in 1:n_particles]
    contact_count = zeros(Int, n_particles)

    n_rb = length(rigidbodies)
    cm_correction = [@SVector zeros(3) for _ in 1:n_rb]
    V_correction  = [@SVector zeros(3) for _ in 1:n_rb]
    ω_correction  = zeros(n_rb)
    rb_contact_count = zeros(Int, n_rb)

    pending_breaks = Tuple{rigidbody_struct, Tuple{Int,Int}}[]
    cells = keys(id_grid)

    for cell in cells

        i, j, k = cell
        cell_particles = id_grid[cell]

        for a in 1:length(cell_particles)
            for b in a+1:length(cell_particles)
                resolve_pair!(particles, rigidbodies, powder, liquid, gas, id_grid, cell_of_particle, cell_particles[a], cell_particles[b], pos_correction, vel_correction, contact_count, cm_correction, V_correction, ω_correction, rb_contact_count, pending_breaks)
            end
        end

        for di in -1:1
            for dj in -1:1
                if di == 0 && dj == 0
                    continue
                end
                if di < 0
                    continue
                end
                if di == 0 && dj < 0
                    continue
                end

                ni = i + di
                nj = j + dj

                if ni >= 1 && ni <= pixel_size_x && nj >= 1 && nj <= pixel_size_y
                    if haskey(id_grid, (ni, nj))
                        neighbor_particles = id_grid[(ni, nj)]
                        for a in cell_particles
                            for b in neighbor_particles
                                resolve_pair!(particles, rigidbodies, powder, liquid, gas, id_grid, cell_of_particle, a, b, pos_correction, vel_correction, contact_count, cm_correction, V_correction, ω_correction, rb_contact_count, pending_breaks)
                            end
                        end
                    end
                end
            end
        end
    end

    for i in 1:n_particles
        p = particles[i]
        if p.active == 1 && p.rigidbody == 0 && contact_count[i] > 0
            p.position = p.position + pos_correction[i]
            p.velocity = clamp_velocity(p.velocity + vel_correction[i], max_velocity)
        end
    end

    for rb in rigidbodies
        nc = max(rb_contact_count[rb.id], 1)

        rb.cm = rb.cm + cm_correction[rb.id] / nc
        rb.V  = clamp_velocity(rb.V + V_correction[rb.id] / nc, max_velocity)
        #rb.ω  = SVector(rb.ω[1], rb.ω[2], rb.ω[3] + ω_correction[rb.id] / nc)
        rb.ω  = SVector(rb.ω[1], rb.ω[2], clamp_angular_velocity(rb.ω[3] + ω_correction[rb.id] / nc, max_angular_velocity))

        for idx in rb.particle_indices
            particles[idx].position = particles[idx].position + cm_correction[rb.id] / nc
        end
    end

    # only one break per rigidbody per frame
    processed_rigidbody_ids = Set{Int}()
    for (rb, broken_bond) in pending_breaks
        if rb.id in processed_rigidbody_ids
            continue
        end
        push!(processed_rigidbody_ids, rb.id)
        split_rigidbody!(particles, rigidbodies, rb, broken_bond)
    end
end

# -----------------------------------------------------------------------------
#                           Single resolve function — handles every pair type
# -----------------------------------------------------------------------------
function resolve_pair!(particles, rigidbodies, powder, liquid, gas, id_grid, cell_of_particle, i, j, pos_correction, vel_correction, contact_count, cm_correction, V_correction, ω_correction, rb_contact_count, pending_breaks)

    n_particles = length(pos_correction)  
    if i > n_particles || j > n_particles
        return   
    end

    p1 = particles[i]
    p2 = particles[j]

    if p1.active == 0 && p2.active == 0
        return
    end
    if p1.collision == 0 || p2.collision == 0
        return
    end 
    if p1.rigidbody != 0 && p1.rigidbody == p2.rigidbody
        return   # same rigidbody, never self-collide
    end
    if p1.softbody != 0 && p1.softbody == p2.softbody
        return   # same softbody, handled by its own constraints
    end

    material_transformation(p1,p2)

    if p1.material == "liquid" && p2.material == "gas"
        return
    elseif p1.material == "gas" && p2.material == "liquid"
        return
    elseif p1.material == "powder" && p2.material == "gas"
        return
    elseif p1.material == "gas" && p2.material == "powder"
        return
    end

    r_vec = p1.position - p2.position
    r = norm(r_vec)

    if r >= collision_min_distance || r < 0.0001
        return
    end

    overlap = collision_min_distance - r
    x1, x2, x3 = p1.position, p2.position, p3.position
    v1, v2, v3 = p1.velocity, p2.velocity, p3.velocity
    normal = (x1 - x2) / r
    r_sq = r^2

    # ---- Case 1: both are rigidbody particles ----
    if p1.rigidbody != 0 && p2.rigidbody != 0

        rb1 = rigidbodies[p1.rigidbody]
        rb2 = rigidbodies[p2.rigidbody]
        m1, m2 = rb1.M, rb2.M
        total_mass = m1 + m2

        rb_contact_count[rb1.id] += 1
        rb_contact_count[rb2.id] += 1
   
        dv1 = (1 + restitution) * m2 / (m1 + m2) * dot(v1 - v2, x1 - x2) * (x1 - x2) / r_sq
        dv2 = - (1 + restitution) * m1 / (m1 + m2) * dot(v2 - v1, x2 - x1) * (x2 - x1) / r_sq

        shift1 = overlap * normal * (m2 / total_mass)
        shift2 = overlap * normal * (m1 / total_mass)

        Δp1 = m1 * dv1
        Δp2 = m2 * dv2

        # friction correction
        a = abs(normal[1]) < 0.9 ? SVector(1.0, 0.0, 0.0) : SVector(0.0, 1.0, 0.0)  # choose a direction not parallel to n

        tangent = normalize(cross(a, normal))  # tangent direction in the plane perpendicular to n

        v_rel = v1 - v2
        vt = dot(v_rel, tangent)

        jn = norm(Δp1)
        jt = -vt * (m1 * m2 / (m1 + m2))
        jt = clamp(jt, -friction_coef * jn, friction_coef * jn)

        friction_impulse = jt * tangent

        Δp1 = Δp1 + friction_impulse
        Δp2 = Δp2 - friction_impulse    

        r1_rel = p1.position - rb1.cm
        r2_rel = p2.position - rb2.cm

        I1 = calculate_inertia_tensor(particles, rb1)
        I2 = calculate_inertia_tensor(particles, rb2)
        invI1 = inv(I1)
        invI2 = inv(I2)
        tau1 = cross(r1_rel, Δp1)
        tau2 = cross(r2_rel, Δp2)

        cm_correction[rb1.id] = cm_correction[rb1.id] + shift1
        cm_correction[rb2.id] = cm_correction[rb2.id] - shift2

        V_correction[rb1.id] = V_correction[rb1.id] + Δp1 / m1
        V_correction[rb2.id] = V_correction[rb2.id] + Δp2 / m2

        # ω = I⁻¹ * (r × Δp)
        ω_correction[rb1.id] += restitution_angular * (invI1 * tau1)
        ω_correction[rb2.id] += restitution_angular * (invI2 * tau2)

    # ---- Case 2: only p1 is a rigidbody particle ----
    elseif p1.rigidbody != 0 && p2.rigidbody == 0

        rb1 = rigidbodies[p1.rigidbody]
        m1, m2 = rb1.M, p2.mass
        total_mass = m1 + m2

        rb_contact_count[rb1.id] += 1

        dv1 = (1 + restitution) * m2 / (m1 + m2) * dot(v1 - v2, x1 - x2) * (x1 - x2) / r_sq
        dv2 = - (1 + restitution) * m1 / (m1 + m2) * dot(v2 - v1, x2 - x1) * (x2 - x1) / r_sq

        shift = overlap * normal * (m2 / total_mass)
        Δp1 = m1 * dv1

        # friction correction
        a = abs(normal[1]) < 0.9 ? SVector(1.0, 0.0, 0.0) : SVector(0.0, 1.0, 0.0)  
        tangent = normalize(cross(a, normal))
        v_rel = v1 - v2
        vt = dot(v_rel, tangent)

        jn = norm(Δp1)
        jt = -vt * (m1 * m2 / (m1 + m2))
        jt = clamp(jt, -friction_coef * jn, friction_coef * jn)

        friction_impulse = jt * tangent
        Δp1 = Δp1 + friction_impulse

        r1_rel = p1.position - rb1.cm
        I1 = calculate_inertia(particles, rb1)
        invI1 = inv(I1)
        tau1 = cross(r1_rel, Δp1)

        cm_correction[rb1.id] = cm_correction[rb1.id] + shift
        V_correction[rb1.id] = V_correction[rb1.id] + Δp1 / m1
        ω_correction[rb1.id] += restitution_angular * (invI1 * tau1)

        if p2.active == 1
            p2.position = p2.position - overlap * normal * (m1 / total_mass)
            p2.velocity = p2.velocity + dv2
        end

        if norm(p2.velocity) > rb1.break_threshold && p2.material == "powder"

            # returns the positions within the rb1.bonds array where the predicate is true
            idxs = findall(bond_tuple -> i in bond_tuple, rb1.bonds)
   
            if !isempty(idxs)
                # deleta os bonds quebrados e adiciona à lista de pending_breaks
                for k in sort(idxs, rev=true)
                    the_bond1 = rb1.bonds[k]
                    deleteat!(rb1.bonds, k)
                    push!(pending_breaks, (rb1, the_bond1))
                    erase_particle!(particles[i], id_grid, cell_of_particle)
                    erase_particle!(particles[j], id_grid, cell_of_particle)
                end
                println("bonds broken: $(idxs)")
            end
        end

    # ---- Case 3: only p2 is a rigidbody particle ----
    elseif p1.rigidbody == 0 && p2.rigidbody != 0

        rb2 = rigidbodies[p2.rigidbody]
        m1, m2 = p1.mass, rb2.M
        total_mass = m1 + m2

        rb_contact_count[rb2.id] += 1

        dv1 = (1 + restitution) * m2 / (m1 + m2) * dot(v1 - v2, x1 - x2) * (x1 - x2) / r_sq
        dv2 = - (1 + restitution) * m1 / (m1 + m2) * dot(v2 - v1, x2 - x1) * (x2 - x1) / r_sq

        if p1.active == 1
            p1.position = p1.position + overlap * normal * (m2 / total_mass)
            p1.velocity = p1.velocity + dv1
        end

        shift = overlap * normal * (m1 / total_mass)
        Δp2 = m2 * dv2

        # friction correction
        a = abs(normal[1]) < 0.9 ? SVector(1.0, 0.0, 0.0) : SVector(0.0, 1.0, 0.0)  
        tangent = normalize(cross(a, normal)) 
        v_rel = v1 - v2
        vt = dot(v_rel, tangent)

        jn = norm(Δp2)
        jt = -vt * (m1 * m2 / (m1 + m2))
        jt = clamp(jt, -friction_coef * jn, friction_coef * jn)

        friction_impulse = jt * tangent
        Δp2 = Δp2 - friction_impulse

        r2_rel = p2.position - rb2.cm
        I2 = calculate_inertia(particles, rb2)
        invI2 = inv(I2)
        tau2 = cross(r2_rel, Δp2)

        cm_correction[rb2.id] = cm_correction[rb2.id] - shift
        V_correction[rb2.id] = V_correction[rb2.id] + Δp2 / m2
        ω_correction[rb2.id] += restitution_angular * (invI2 * tau2)

        if norm(p1.velocity) > rb2.break_threshold && p1.material == "powder"
            # j is the rb particle index
            # array of indices of all bonds that include j
            idxs = findall(bond_tuple -> j in bond_tuple, rb2.bonds)
            if !isempty(idxs)
                for k in sort(idxs, rev=true)
                    the_bond2 = rb2.bonds[k]
                    deleteat!(rb2.bonds, k)
                    push!(pending_breaks, (rb2, the_bond2))
                    erase_particle!(particles[j], id_grid, cell_of_particle)
                    erase_particle!(particles[i], id_grid, cell_of_particle)
                end
                println("bonds broken: $(idxs)")
            end
        end

    # ---- Case 4: neither is a rigidbody — covers free particles AND softbody particles ----
    else

        m1, m2 = p1.mass, p2.mass
        total_mass = m1 + m2

        dv1 = (1 + restitution) * m2 / (m1 + m2) * dot(v1 - v2, x1 - x2) * (x1 - x2) / r_sq
        dv2 = - (1 + restitution) * m1 / (m1 + m2) * dot(v2 - v1, x2 - x1) * (x2 - x1) / r_sq

        if p1.active == 1
            pos_correction[i] = pos_correction[i] + overlap * normal * (m2 / total_mass)
            vel_correction[i] = vel_correction[i] + dv1
            contact_count[i] += 1
        end

        if p2.active == 1
            pos_correction[j] = pos_correction[j] + (-overlap * normal * (m1 / total_mass))
            vel_correction[j] = vel_correction[j] + dv2
            contact_count[j] += 1
        end
    end
end

function clamp_velocity(v, max_speed)
    speed = norm(v)
    if speed > max_speed
        return v * (max_speed / speed)
    end
    return v
end

function clamp_angular_velocity(ω, max_ω)
    if abs(ω) > max_ω
        return sign(ω) * max_ω
    end
    return ω
end

function transform_particle!(particles, target_array, id_grid, cell_of_particle, p, p2, new_particle)

    if p.active == 0 || p2.active == 0
        return   # already transformed earlier this same scan
    end

    px = Int(floor(p.position[1] / grid_size)) + 1
    py = Int(floor(p.position[2] / grid_size)) + 1
    pz = Int(floor(p.position[3] / grid_size)) + 1

    # remove p from its grid cell 
    erase_particle!(p, id_grid, cell_of_particle)

    # remove p2 from its grid cell 
    erase_particle!(p2, id_grid, cell_of_particle)

    # spawn the new gas particle
    push!(target_array, new_particle)
    push!(particles, new_particle)
    push!(cell_of_particle, (px, py, pz))

    if !haskey(id_grid, (px, py, pz))
        id_grid[(px, py, pz)] = Int[]
    end
    push!(id_grid[(px, py, pz)], new_particle.id)
end