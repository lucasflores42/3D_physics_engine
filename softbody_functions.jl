# -----------------------------------------------------------------------------
#                           SoftBody physics
# ----------------------------------------------------------------------------- 
mutable struct softbody_struct
    particle_indices::Vector{Int}
    constraints::Vector{Tuple{Int,Int,Float64}}   # (local_i, local_j, rest_length)
    stiffness::Float64
    pinned::Vector{Bool}
end

function softbody_physics(particles, softbodies)
    for sb in softbodies

        # integrate: gravity + velocity, skip pinned particles
        for k in 1:length(sb.particle_indices)
            if sb.pinned[k]
                continue
            end
            idx = sb.particle_indices[k]
            p = particles[idx]

            F_gravity = SVector(0.0, -10.0) * p.mass
            p.velocity = p.velocity + (F_gravity / p.mass) * dt
            p.position = p.position + p.velocity * dt
        end


        # enforce constraints: pull each connected pair back to rest_length
        for (li, lj, rest_length) in sb.constraints
            i, j = sb.particle_indices[li], sb.particle_indices[lj]
            p1, p2 = particles[i], particles[j]

            delta = p2.position - p1.position
            dist = norm(delta)
            if dist > 0.0001
                diff = (dist - rest_length) / dist
                correction = delta * diff * 0.5 * sb.stiffness

                if !sb.pinned[li]
                    p1.position = p1.position + correction
                    p1.velocity = p1.velocity + correction / dt
                end
                if !sb.pinned[lj]
                    p2.position = p2.position - correction
                    p2.velocity = p2.velocity - correction / dt
                end
            end
        end
    end
end

function create_rope!(particles, softbodies, id, offset, n, stiffness, rest_length)
    indices = Int[]

    for k in 0:(n-1)
        p = solid_struct(
            offset .+ SVector(-k*rest_length, 0),
            @SVector(zeros(2)),
            @SVector(zeros(2)),
            grid_size/2,
            1.0,
            0,
            id,       # softbody
            1, 1, 1,
            "solid"
            )
        push!(particles, p)
        push!(indices, length(particles))
    end

    constraints = Tuple{Int,Int,Float64}[]
    for k in 1:(n-1)
        push!(constraints, (k, k+1, rest_length))
    end

    pinned = fill(false, n)
    pinned[1] = true

    push!(softbodies, softbody_struct(indices, constraints, stiffness, pinned))
end