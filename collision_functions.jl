# -----------------------------------------------------------------------------
#                           Parameters
# -----------------------------------------------------------------------------
const restitution_x = 0.5
const restitution_y = 0.0
const collision_min_distance = grid_size #* sqrt(2)

include("rigidbody_functions.jl")

# -----------------------------------------------------------------------------
#                           Calculate all collisions
# -----------------------------------------------------------------------------
function collision_physics!(particles, rigidbodies, id_grid)

    n_particles = length(particles)
    pos_correction = [@SVector zeros(2) for _ in 1:n_particles]
    vel_correction = [@SVector zeros(2) for _ in 1:n_particles]
    contact_count = zeros(Int, n_particles)

    n_rb = length(rigidbodies)
    cm_correction = [@SVector zeros(2) for _ in 1:n_rb]
    V_correction  = [@SVector zeros(2) for _ in 1:n_rb]
    ω_correction  = zeros(n_rb)
    rb_contact_count = zeros(Int, n_rb)

    cells = keys(id_grid)

    for cell in cells

        i, j = cell
        cell_particles = id_grid[cell]

        for a in 1:length(cell_particles)
            for b in a+1:length(cell_particles)
                resolve_pair!(particles, rigidbodies, cell_particles[a], cell_particles[b], pos_correction, vel_correction, contact_count, cm_correction, V_correction, ω_correction, rb_contact_count)
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
                                resolve_pair!(particles, rigidbodies, a, b, pos_correction, vel_correction, contact_count, cm_correction, V_correction, ω_correction, rb_contact_count)
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
            p.velocity = p.velocity + vel_correction[i]
        end
    end

    for rb in rigidbodies
        nc = max(rb_contact_count[rb.id], 1)

        rb.cm = rb.cm + cm_correction[rb.id] / nc
        rb.V  = rb.V  + V_correction[rb.id] / nc
        rb.ω  = SVector(rb.ω[1], rb.ω[2], rb.ω[3] + ω_correction[rb.id] / nc)

        for idx in rb.particle_indices
            particles[idx].position = particles[idx].position + cm_correction[rb.id] / nc
        end
    end
end

# -----------------------------------------------------------------------------
#                           Single resolve function — handles every pair type
# -----------------------------------------------------------------------------
function resolve_pair!(particles, rigidbodies, i, j, pos_correction, vel_correction, contact_count, cm_correction, V_correction, ω_correction, rb_contact_count)

    p1 = particles[i]
    p2 = particles[j]

    if p1.active == 0 && p2.active == 0
        return
    end
    if p1.rigidbody != 0 && p1.rigidbody == p2.rigidbody
        return   # same rigidbody, never self-collide
    end
    if p1.softbody != 0 && p1.softbody == p2.softbody
        return   # same softbody, handled by its own constraints
    end
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
    x1, x2 = p1.position, p2.position
    v1, v2 = p1.velocity, p2.velocity
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

        #dv1 = - (1 + colision_restitution_coefficient) * m2 / (m1 + m2) * dot(v1 - v2, x1 - x2) * (x1 - x2) / r_sq
        #dv2 = - (1 + colision_restitution_coefficient) * m1 / (m1 + m2) * dot(v2 - v1, x2 - x1) * (x2 - x1) / r_sq

        dv1_x = - (1 + restitution_x) * m2 / (m1 + m2) * dot(v1 - v2, x1 - x2) * (x1[1] - x2[1]) / r_sq
        dv1_y = - (1 + restitution_y) * m2 / (m1 + m2) * dot(v1 - v2, x1 - x2) * (x1[2] - x2[2]) / r_sq
        dv1 = [dv1_x, dv1_y]
        
        dv2_x = - (1 + restitution_x) * m1 / (m1 + m2) * dot(v2 - v1, x2 - x1) * (x2[1] - x1[1]) / r_sq
        dv2_y = - (1 + restitution_y) * m1 / (m1 + m2) * dot(v2 - v1, x2 - x1) * (x2[2] - x1[2]) / r_sq
        dv2 = [dv2_x, dv2_y]

        shift1 = overlap * normal * (m2 / total_mass)
        shift2 = overlap * normal * (m1 / total_mass)

        Δp1 = m1 * dv1
        Δp2 = m2 * dv2

        r1_rel = p1.position - rb1.cm
        r2_rel = p2.position - rb2.cm

        I1 = calculate_inertia(particles, rb1)
        I2 = calculate_inertia(particles, rb2)

        cm_correction[rb1.id] = cm_correction[rb1.id] + shift1
        cm_correction[rb2.id] = cm_correction[rb2.id] - shift2

        V_correction[rb1.id] = V_correction[rb1.id] + Δp1 / m1
        V_correction[rb2.id] = V_correction[rb2.id] + Δp2 / m2

        ω_correction[rb1.id] += (r1_rel[1]*Δp1[2] - r1_rel[2]*Δp1[1]) / I1
        ω_correction[rb2.id] += (r2_rel[1]*Δp2[2] - r2_rel[2]*Δp2[1]) / I2

    # ---- Case 2: only p1 is a rigidbody particle ----
    elseif p1.rigidbody != 0 && p2.rigidbody == 0

        rb1 = rigidbodies[p1.rigidbody]
        m1, m2 = rb1.M, p2.mass
        total_mass = m1 + m2

        rb_contact_count[rb1.id] += 1

        dv1_x = - (1 + restitution_x) * m2 / (m1 + m2) * dot(v1 - v2, x1 - x2) * (x1[1] - x2[1]) / r_sq
        dv1_y = - (1 + restitution_y) * m2 / (m1 + m2) * dot(v1 - v2, x1 - x2) * (x1[2] - x2[2]) / r_sq
        dv1 = [dv1_x, dv1_y]
        
        dv2_x = - (1 + restitution_x) * m1 / (m1 + m2) * dot(v2 - v1, x2 - x1) * (x2[1] - x1[1]) / r_sq
        dv2_y = - (1 + restitution_y) * m1 / (m1 + m2) * dot(v2 - v1, x2 - x1) * (x2[2] - x1[2]) / r_sq
        dv2 = [dv2_x, dv2_y]

        shift = overlap * normal * (m2 / total_mass)
        Δp1 = m1 * dv1
        r1_rel = p1.position - rb1.cm
        I1 = calculate_inertia(particles, rb1)

        cm_correction[rb1.id] = cm_correction[rb1.id] + shift
        V_correction[rb1.id] = V_correction[rb1.id] + Δp1 / m1
        ω_correction[rb1.id] += (r1_rel[1]*Δp1[2] - r1_rel[2]*Δp1[1]) / I1

        if p2.active == 1
            p2.position = p2.position - overlap * normal * (m1 / total_mass)
            p2.velocity = p2.velocity + dv2
        end

    # ---- Case 3: only p2 is a rigidbody particle ----
    elseif p1.rigidbody == 0 && p2.rigidbody != 0

        rb2 = rigidbodies[p2.rigidbody]
        m1, m2 = p1.mass, rb2.M
        total_mass = m1 + m2

        rb_contact_count[rb2.id] += 1

        dv1_x = - (1 + restitution_x) * m2 / (m1 + m2) * dot(v1 - v2, x1 - x2) * (x1[1] - x2[1]) / r_sq
        dv1_y = - (1 + restitution_y) * m2 / (m1 + m2) * dot(v1 - v2, x1 - x2) * (x1[2] - x2[2]) / r_sq
        dv1 = [dv1_x, dv1_y]
        
        dv2_x = - (1 + restitution_x) * m1 / (m1 + m2) * dot(v2 - v1, x2 - x1) * (x2[1] - x1[1]) / r_sq
        dv2_y = - (1 + restitution_y) * m1 / (m1 + m2) * dot(v2 - v1, x2 - x1) * (x2[2] - x1[2]) / r_sq
        dv2 = [dv2_x, dv2_y]

        if p1.active == 1
            p1.position = p1.position + overlap * normal * (m2 / total_mass)
            p1.velocity = p1.velocity + dv1
        end

        shift = overlap * normal * (m1 / total_mass)
        Δp2 = m2 * dv2
        r2_rel = p2.position - rb2.cm
        I2 = calculate_inertia(particles, rb2)

        cm_correction[rb2.id] = cm_correction[rb2.id] - shift
        V_correction[rb2.id] = V_correction[rb2.id] + Δp2 / m2
        ω_correction[rb2.id] += (r2_rel[1]*Δp2[2] - r2_rel[2]*Δp2[1]) / I2

    # ---- Case 4: neither is a rigidbody — covers free particles AND softbody particles ----
    else

        m1, m2 = p1.mass, p2.mass
        total_mass = m1 + m2

        dv1_x = - (1 + restitution_x) * m2 / (m1 + m2) * dot(v1 - v2, x1 - x2) * (x1[1] - x2[1]) / r_sq
        dv1_y = - (1 + restitution_y) * m2 / (m1 + m2) * dot(v1 - v2, x1 - x2) * (x1[2] - x2[2]) / r_sq
        dv1 = [dv1_x, dv1_y]
        
        dv2_x = - (1 + restitution_x) * m1 / (m1 + m2) * dot(v2 - v1, x2 - x1) * (x2[1] - x1[1]) / r_sq
        dv2_y = - (1 + restitution_y) * m1 / (m1 + m2) * dot(v2 - v1, x2 - x1) * (x2[2] - x1[2]) / r_sq
        dv2 = [dv2_x, dv2_y]

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