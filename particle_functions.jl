# -----------------------------------------------------------------------------
#                           Parameters
# ----------------------------------------------------------------------------- 
const gravity_coef = 0.1
const colision_restitution_coefficient = 0.5
const collision_min_distance = grid_size #* sqrt(2)

# -----------------------------------------------------------------------------
#                           Particles physics
# ----------------------------------------------------------------------------- 
mutable struct powder_struct
    id::Int64
    position::SVector{2, Float64}
    velocity::SVector{2, Float64}
    acceleration::SVector{2, Float64}
    radius::Float64
    mass::Float64

    rigidbody::Int64
    softbody::Int64

    active::Int64
    collision::Int64
    gravity::Int64

    material::String
end

mutable struct gas_struct
    id::Int64
    position::SVector{2, Float64}
    velocity::SVector{2, Float64}
    acceleration::SVector{2, Float64}
    radius::Float64 # half the grid size
    mass::Float64

    rigidbody::Int64
    softbody::Int64

    active::Int64
    collision::Int64
    gravity::Int64
    sph::Int64

    time_active::Int64
    lifetime::Int64

    material::String
end

mutable struct liquid_struct
    id::Int64
    position::SVector{2, Float64}
    velocity::SVector{2, Float64}
    acceleration::SVector{2, Float64}
    radius::Float64
    mass::Float64

    rigidbody::Int64
    softbody::Int64

    density::Float64
    pressure::Float64
    target_density::Float64
    stiff_coef::Float64       
    viscosity_coef::Float64 

    active::Int64
    collision::Int64
    gravity::Int64
    sph::Int64

    color_id::Int64
    material::String
end

mutable struct solid_struct
    id::Int64
    position::SVector{2, Float64}
    velocity::SVector{2, Float64}
    acceleration::SVector{2, Float64}
    radius::Float64
    mass::Float64

    rigidbody::Int64
    softbody::Int64

    active::Int64
    collision::Int64
    gravity::Int64

    material::String
end

function particle_physics(particles, liquid, liquid2, gas, powder, solid, id_grid, cell_of_particle)


    all_liquids = [liquid; liquid2]
    
    for p in all_liquids

        if p.active == 0
            continue
        end

        if p.sph == 1

            grad_pressure = @SVector zeros(2)
            laplacian_velocity = @SVector zeros(2)
            F_surface = @SVector zeros(2)
            calculate_density_pressure!(p, particles, id_grid)

            px = Int(floor(p.position[1] / grid_size)) + 1
            py = Int(floor(p.position[2] / grid_size)) + 1

            for di in -sph_cell_range:sph_cell_range
                for dj in -sph_cell_range:sph_cell_range
                    ni, nj = px + di, py + dj
                    if ni < 1 || ni > pixel_size_x || nj < 1 || nj > pixel_size_y || !haskey(id_grid, (ni, nj))
                        continue
                    end
                    for j in id_grid[(ni, nj)]
                        p2 = particles[j]
                        if p2.material != "liquid" ||  p2 === p
                            continue
                        end

                        r_vec = p.position - p2.position
                        r = norm(r_vec)

                        if r > 2*smoothing_length || r == 0
                            continue
                        end

                        grad_pressure += pressure_gradient(p, p2, r, r_vec)
                        laplacian_velocity += viscosity_laplacian(p, p2, r, r_vec)
                        #F_surface += surface_tension_force(p, p2, r, r_vec)
                    end
                end
            end

            F_pressure = -grad_pressure
            F_viscosity = p.mass * p.viscosity_coef * laplacian_velocity
        else
            F_pressure = @SVector zeros(2)
            F_viscosity = @SVector zeros(2)
        end

        if p.gravity == 1
            F_gravity = calculate_gravity(p.position, p.mass, 0, solid)
        else
            F_gravity = @SVector zeros(2)
        end
        
        F_total = F_gravity + F_pressure + F_viscosity + F_surface
        p.acceleration = F_total / p.mass
        p.velocity += p.acceleration * dt
        p.position += p.velocity * dt
    end

    for i in 1:length(gas)
       
        p = gas[i]
        if p.active == 0
            continue
        end

        p.time_active += 1

        if p.time_active >= p.lifetime
            erase_particle!(p, id_grid, cell_of_particle)
            continue   # skip the rest of this particle's physics 
        end

        if p.gravity == 1
            F_gravity = -1 * calculate_gravity(p.position, p.mass, 0, solid)
        else
            F_gravity = @SVector zeros(2)
        end
        
        F_total = F_gravity
        p.acceleration = F_total / p.mass
        p.velocity += p.acceleration * dt
        p.position += p.velocity * dt

    end

    for i in 1:length(powder)

        p = powder[i]
        if p.active == 0
            continue
        end


        if p.gravity == 1
            F_gravity = calculate_gravity(p.position, p.mass, 0, solid)
        else
            F_gravity = @SVector zeros(2)
        end
        
        F_total = F_gravity
        p.acceleration = F_total / p.mass
        p.velocity += p.acceleration * dt
        p.position += p.velocity * dt        
    end

    for i in 1:length(solid)

        p = solid[i]
        if p.active == 0
            continue
        end

        if p.gravity == 1
            F_gravity = calculate_gravity(p.position, p.mass, 0, solid)
        else
            F_gravity = @SVector zeros(2)
        end
            
        F_total = F_gravity
        p.acceleration = F_total / p.mass
        p.velocity += p.acceleration * dt
        p.position += p.velocity * dt
    end
end

function calculate_gravity(position, mass, id, solid)

    F_gravity = @SVector zeros(2)
    
    #=
    for j in 1:length(solid)

        p = solid[j]
              
        if p.material == "solid" && p.rigidbody != id
            r_vec = p.position - position
            r = norm(r_vec)
            
            if r > 0.0001  
                F_gravity += gravity_coef * mass * p.mass * r_vec / (r ^ 3)
            end
        end
    end
    =#
    
    #return F_gravity
    return mass * SVector(0.0, -10.0)
end

function erase_particle!(p, id_grid, cell_of_particle)
    
    #px = Int(floor(p.position[1] / grid_size)) + 1
    #py = Int(floor(p.position[2] / grid_size)) + 1
    px, py = cell_of_particle[p.id]

    if haskey(id_grid, (px, py))
        cell_ids = id_grid[(px, py)]
        filter!(x -> x != p.id, cell_ids)

        if isempty(cell_ids)
            delete!(id_grid, (px, py))
        end
    end

    p.active = 0
    p.collision = 0
end