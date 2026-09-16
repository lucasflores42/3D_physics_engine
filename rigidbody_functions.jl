# -----------------------------------------------------------------------------
#                           RigidBody physics
# ----------------------------------------------------------------------------- 
mutable struct rigidbody_struct
    id::Int
    particle_indices::Vector{Int}
    cm::SVector{2, Float64}
    V::SVector{2, Float64}
    ω::SVector{3, Float64}
    M::Float64
    bonds::Vector{Tuple{Int,Int}}   
    break_threshold::Float64
end

function rigidbody_physics(particles, rigidbodies)

    for rb in rigidbodies

        # translation
        #F_gravity = @SVector zeros(2)
        F_gravity = calculate_gravity(rb.cm, rb.M, 0, nothing)

        rb.V += (F_gravity / rb.M) * dt      
        translation = rb.V * dt
        new_cm = rb.cm + translation         

        # rotation
        angle = rb.ω[3] * dt
        cos_a = cos(angle)
        sin_a = sin(angle)

        for idx in rb.particle_indices
            p = particles[idx]

            r = p.position - rb.cm

            r_rot = SVector(cos_a*r[1] - sin_a*r[2], sin_a*r[1] + cos_a*r[2])

            p.position = new_cm + r_rot

            r_new = p.position - new_cm
            p.velocity = rb.V + SVector(-rb.ω[3]*r_new[2], rb.ω[3]*r_new[1])
        end

        rb.cm = new_cm
    end
end

function calculate_inertia(particles, rb)
    inertia = 0.0
    for i in rb.particle_indices
        p = particles[i]
        r = p.position - rb.cm
        inertia += p.mass * (r[1]^2 + r[2]^2)
    end
    return max(inertia, 1.0)
end

function calculate_inertia_tensor(particles, rb)

    I_tensor = zeros(2,2)
    I3 = Matrix{Float64}(I, 2, 2)  # identity

    for i in rb.particle_indices
        
        p = particles[i]
        r = p.position .- rb.cm   

        r2 = dot(r, r)
        rrT = r * transpose(r)
        I_tensor .+= p.mass .* (r2 .* I3 .- rrT)
    end

    return I_tensor
end

function calculate_center_of_mass(particles)
    total_mass = 0.0
    cm_x = 0.0
    cm_y = 0.0
    
    for p in particles
        total_mass += p.mass
        cm_x += p.mass * p.position[1]
        cm_y += p.mass * p.position[2]
    end
    
    return SVector(cm_x / total_mass, cm_y / total_mass), total_mass
end

function create_cube!(particles, rigidbodies, id, offset, v_init, ω_init, m, n)
    particle_radius = grid_size/2
    particle_diam = 2 * particle_radius

    positions = SVector{2,Float64}[]
    for row in 0:(m-1)
        for col in 0:(n-1)
            push!(positions, SVector(col * particle_diam, row * particle_diam))
        end
    end

    indices = Int[]

    for pos in positions
        p = solid_struct(
            length(particles)+1,
            offset .+ pos,
            @SVector(zeros(2)),
            @SVector(zeros(2)),
            particle_radius,
            10.0,            # mass
            id,
            0,
            1,              # active
            1,              # collision
            1,              # gravity
            "solid"
        )
        push!(particles, p)
        push!(indices, length(particles))
    end

    # Calculate center of mass
    cube_particles = [particles[i] for i in indices]
    cm, total_mass = calculate_center_of_mass(cube_particles)

    # Set initial velocities
    for i in indices
        r = particles[i].position - cm
        particles[i].velocity = v_init + SVector(-ω_init[1]*r[2], ω_init[1]*r[1])
    end

    bonds = build_grid_bonds(positions, particle_diam, indices)

    rb = rigidbody_struct(
        id,
        indices, #global indices of particles in the rigidbody
        cm,
        SVector(v_init[1], v_init[2]),
        SVector(0.0, 0.0, ω_init[1]),
        total_mass,
        bonds,
        10
    )
    push!(rigidbodies, rb)
end

function create_sphere!(particles, rigidbodies, id, offset, v_init, ω_init, r)
    particle_radius = grid_size/2
    particle_diam = 2 * particle_radius

    positions = SVector{2,Float64}[]
    for row in 0:(2*r)
        for col in 0:(2*r)
            if (row - r)^2 + (col - r)^2 >= r^2 - 10 && (row - r)^2 + (col - r)^2 <= r^2 + 10
                push!(positions, SVector(col * particle_diam, row * particle_diam))
            end
        end
    end

    indices = Int[]

    for pos in positions
        p = solid_struct(
            length(particles)+1,
            offset .+ pos,
            @SVector(zeros(2)),
            @SVector(zeros(2)),
            particle_radius,
            10.0,            # mass
            id,
            0,
            1,              # active
            1,              # collision
            1,              # gravity
            "solid"
        )
        push!(particles, p)
        push!(indices, length(particles))
    end

    # Calculate center of mass
    cube_particles = [particles[i] for i in indices]
    cm, total_mass = calculate_center_of_mass(cube_particles)

    # Set initial velocities
    for i in indices
        r = particles[i].position - cm
        particles[i].velocity = v_init + SVector(-ω_init[1]*r[2], ω_init[1]*r[1])
    end

    bonds = build_grid_bonds(positions, particle_diam, indices)

    rb = rigidbody_struct(
        id,
        indices, #global indices of particles in the rigidbody
        cm,
        SVector(v_init[1], v_init[2]),
        SVector(0.0, 0.0, ω_init[1]),
        total_mass,
        bonds,
        10
    )
    push!(rigidbodies, rb)
end


function build_grid_bonds(local_positions, spacing, indices)
    bonds = Tuple{Int,Int}[]
    n = length(local_positions)
    for a in 1:n
        for b in a+1:n
            d = norm(local_positions[a] - local_positions[b])
            if d < spacing * 1.1   # adjacent (no diagonals)
                push!(bonds, (indices[a], indices[b]))   # global ids
            end
        end
    end
    return bonds
end

# check the two groups of bounds
# subtract and recalculate the info of the first
# create a new body with the rest
function split_rigidbody!(particles, rigidbodies, rb, broken_bond)

    # drop any particle indices that have been destroyed
    rb.particle_indices = filter(idx -> particles[idx].active == 1, rb.particle_indices)

    if length(rb.particle_indices) <= 1
        return   # nothing left to check connectivity on
    end

    particle_set = Set(rb.particle_indices)

    # also drop any bonds referencing a now-dead particle
    rb.bonds = filter(bond -> bond[1] in particle_set && bond[2] in particle_set, rb.bonds)

    # build adjacency list for the rigidbody's particles
    adj = Dict{Int, Vector{Int}}()
    for pid in rb.particle_indices
        adj[pid] = Int[]
    end

    # add the bonds in both directions
    for (x, y) in rb.bonds
        push!(adj[x], y)
        push!(adj[y], x)
    end

    # BFS from whichever live particle we have, to find what's connected to it
    start = first(rb.particle_indices)
    connected_particles = Set([start])
    particles_to_visit = [start]
    next_particle = 1

    while next_particle <= length(particles_to_visit)
        current_particle = particles_to_visit[next_particle]
        next_particle += 1

        for neighboring_particle in adj[current_particle]
            #if neighboring_particle ∉ connected_particles
            if !(neighboring_particle in connected_particles)
                push!(connected_particles, neighboring_particle)
                push!(particles_to_visit, neighboring_particle)
            end
        end
    end

    # if everything is still reachable, nothing actually split
    if length(connected_particles) == length(rb.particle_indices)
        return
    end

    group_a = collect(connected_particles)
    group_b = [pid for pid in rb.particle_indices if pid ∉ connected_particles]

    # Keep the original rigidbody on the larger component
    if length(group_a) == 1 && length(group_b) > 1
        group_a, group_b = group_b, group_a
        connected_particles = Set(group_a)
    end

    original_bonds = rb.bonds

    # ---- group_a: either stays as rb, or becomes a free particle if alone ----
    if length(group_a) == 1
        lone_index = group_a[1]
        particles[lone_index].rigidbody = 0
    else
        rb.particle_indices = group_a
        new_bonds_a = Tuple{Int,Int}[]
        for bond in original_bonds
            x, y = bond
            if (x in connected_particles) && (y in connected_particles)
                push!(new_bonds_a, (x, y))
            end
        end
        rb.bonds = new_bonds_a

        piece_particles_a = [particles[i] for i in rb.particle_indices]
        cm_a, mass_a = calculate_center_of_mass(piece_particles_a)
        rb.cm = cm_a
        rb.M = mass_a
        println("group_a particles: ", [p.position for p in piece_particles_a], " -> cm=", cm_a)   # <-- HERE


        for i in rb.particle_indices
            particles[i].rigidbody = rb.id
        end
    end

    # ---- group_b: either becomes a free particle, or a brand new rigidbody ----
    if length(group_b) == 1
        lone_index = group_b[1]
        particles[lone_index].rigidbody = 0
    else
        new_particle_indices_b = group_b

        new_bonds_b = Tuple{Int,Int}[]
        for bond in original_bonds
            x, y = bond
            if (x ∉ connected_particles) && (y ∉ connected_particles)
                push!(new_bonds_b, (x, y))
            end
        end

        piece_particles_b = [particles[i] for i in new_particle_indices_b]
        cm_b, mass_b = calculate_center_of_mass(piece_particles_b)
        println("group_b particles: ", [p.position for p in piece_particles_b], " -> cm=", cm_b)   

        new_id = length(rigidbodies) + 1
        for i in new_particle_indices_b
            particles[i].rigidbody = new_id
        end

        new_rigidbody = rigidbody_struct(new_id, new_particle_indices_b, cm_b, rb.V, rb.ω, mass_b, new_bonds_b, rb.break_threshold)
        push!(rigidbodies, new_rigidbody)
    end
end