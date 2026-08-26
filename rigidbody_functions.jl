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
    return inertia 
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
            offset .+ pos,
            @SVector(zeros(2)),
            @SVector(zeros(2)),
            particle_radius,
            3.0,
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

    rb = rigidbody_struct(
        id,
        indices,
        cm,
        SVector(v_init[1], v_init[2]),
        SVector(0.0, 0.0, ω_init[1]),
        total_mass
    )
    push!(rigidbodies, rb)
end