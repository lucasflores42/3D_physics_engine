function clamp_particles(particles)
    for p in particles
        if p.position[1] < grid_size/2
            p.position[1] = grid_size/2
            p.velocity[1] = abs(p.velocity[1])  # Bounce back
        elseif p.position[1] > box_size_x - grid_size/2
            p.position[1] = box_size_x - grid_size/2
            p.velocity[1] = -abs(p.velocity[1])  # Bounce back
        end
        
        if p.position[2] < grid_size/2
            p.position[2] = grid_size/2
            p.velocity[2] = abs(p.velocity[2])  # Bounce back
        elseif p.position[2] > box_size_y - grid_size/2
            p.position[2] = box_size_y - grid_size/2
            p.velocity[2] = -abs(p.velocity[2])  # Bounce back
        end

        if p.position[3] < grid_size/2
            p.position[3] = grid_size/2
            p.velocity[3] = abs(p.velocity[3])  # Bounce back
        elseif p.position[3] > box_size_z - grid_size/2
            p.position[3] = box_size_z - grid_size/2
            p.velocity[3] = -abs(p.velocity[3])  # Bounce back
        end
    end
end

function material_code(material, color_id=0)
    if material == "solid"
        return 1
    elseif material == "liquid"
        return color_id == 2 ? 5 : 2
    elseif material == "gas"
        return 3
    elseif material == "powder"
        return 4
    else
        return 0
    end
end
function init_grids(particles)
    id_grid = Dict{Tuple{Int,Int,Int}, Vector{Int}}()
    cell_of_particle = Vector{Tuple{Int,Int,Int}}(undef, length(particles))

    for i in 1:length(particles)
        p = particles[i]
        px = Int(floor(p.position[1] / grid_size)) + 1
        py = Int(floor(p.position[2] / grid_size)) + 1
        pz = Int(floor(p.position[3] / grid_size)) + 1
        cell_of_particle[i] = (px, py, pz)

        if 1 <= px <= voxel_size_x && 1 <= py <= voxel_size_y && 1 <= pz <= voxel_size_z && p.collision == 1
            if !haskey(id_grid, (px, py, pz))
                id_grid[(px, py, pz)] = Int[]
            end
            push!(id_grid[(px, py, pz)], i)
        end
    end

    return id_grid, cell_of_particle
end

function update_grids!(particles, id_grid, cell_of_particle)
    for i in 1:length(particles)
        p = particles[i]
        if p.active == 0
            continue
        end

        px = Int(floor(p.position[1] / grid_size)) + 1
        py = Int(floor(p.position[2] / grid_size)) + 1
        pz = Int(floor(p.position[3] / grid_size)) + 1
        new_cell = (px, py, pz)
        old_cell = cell_of_particle[i]

        if new_cell == old_cell
            continue
        end

        if p.collision == 1 && haskey(id_grid, old_cell)
            filter!(x -> x != i, id_grid[old_cell])
            if isempty(id_grid[old_cell])
                delete!(id_grid, old_cell)
            end
        end

        if p.collision == 1
            if !haskey(id_grid, new_cell) # If new_cell is not a key in id_grid...
                id_grid[new_cell] = Int[]
            end
            push!(id_grid[new_cell], i)
        end

        cell_of_particle[i] = new_cell
    end
end

const material_priority = Dict("solid" => 1, "powder" => 2, "liquid" => 3, "gas" => 4)

function build_material_grid(particles, id_grid)
    material_grid = zeros(Int, voxel_size_x, voxel_size_y, voxel_size_z)

    for (cell, ids) in id_grid
        px, py, pz = cell

        if px < 1 || px > voxel_size_x || py < 1 || py > voxel_size_y || pz < 1 || pz > voxel_size_z
            continue
        end

        best_material = nothing
        best_color_id = 0
        best_priority = typemax(Int)

        for i in ids
            p = particles[i]
            pr = get(material_priority, p.material, 99)
            if pr < best_priority
                best_priority = pr
                best_material = p.material
                best_color_id = p.material == "liquid" ? p.color_id : 0
            end
        end

        if best_material !== nothing
            material_grid[px, py, pz] = material_code(best_material, best_color_id)
        end
    end

    return material_grid
end
