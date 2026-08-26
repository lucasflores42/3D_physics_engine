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
            p.position[2] = box_size - grid_size/2
            p.velocity[2] = -abs(p.velocity[2])  # Bounce back
        end
    end
end

function material_code(material)
    if material == "solid"
        return 1
    elseif material == "liquid"
        return 2
    elseif material == "gas"
        return 3
    elseif material == "powder"
        return 4
    else
        return 0
    end
end

function init_grids(particles)
    id_grid = Dict{Tuple{Int,Int}, Vector{Int}}()
    cell_of_particle = Vector{Tuple{Int,Int}}(undef, length(particles))

    for i in 1:length(particles)
        p = particles[i]
        px = Int(floor(p.position[1] / grid_size)) + 1
        py = Int(floor(p.position[2] / grid_size)) + 1
        cell_of_particle[i] = (px, py)

        if 1 <= px <= pixel_size_x && 1 <= py <= pixel_size_y && p.collision == 1
            if !haskey(id_grid, (px, py))
                id_grid[(px, py)] = Int[]
            end
            push!(id_grid[(px, py)], i)
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
        new_cell = (px, py)
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
            if !haskey(id_grid, new_cell)
                id_grid[new_cell] = Int[]
            end
            push!(id_grid[new_cell], i)
        end

        cell_of_particle[i] = new_cell
    end
end

const material_priority = Dict("solid" => 1, "powder" => 2, "liquid" => 3, "gas" => 4)

function build_material_grid(particles, id_grid)
    material_grid = zeros(Int, pixel_size_x, pixel_size_y)

    for (cell, ids) in id_grid
        px, py = cell

        if px < 1 || px > pixel_size_x || py < 1 || py > pixel_size_y
            continue
        end

        best_material = nothing
        best_priority = typemax(Int)

        for i in ids
            p = particles[i]
            pr = get(material_priority, p.material, 99)
            if pr < best_priority
                best_priority = pr
                best_material = p.material
            end
        end

        if best_material !== nothing
            material_grid[px, py] = material_code(best_material)
        end
    end

    return material_grid
end 

function spawn_particle!(particles, liquid, gas, powder, solid, id_grid, cell_of_particle, px, py, material)

    x = (px - 1) * grid_size + grid_size/2
    y = (py - 1) * grid_size + grid_size/2

    if material == "powder"
        p = powder_struct(SVector(x,y), @SVector(zeros(2)), @SVector(zeros(2)),
                           grid_size/2, 10.0, 0, 0, 1, 1, 1, "powder")
        push!(powder, p)
    elseif material == "liquid"
        p = liquid_struct(SVector(x,y), SVector(1*rand(),0.0), @SVector(zeros(2)),
                           grid_size/2, 0.1, 0, 0, 1000, 0.0, 1, 1, 1, 1, "liquid")
        push!(liquid, p)
    elseif material == "gas"
        p = gas_struct(SVector(x,y), @SVector(zeros(2)), @SVector(zeros(2)),
                        grid_size/2, 0.1, 0, 0, 1, 1, 1, 1, "gas")
        push!(gas, p)
    elseif material == "solid"
        p = solid_struct(SVector(x,y), @SVector(zeros(2)), @SVector(zeros(2)),
                          grid_size/2, 1.0, 0, 0, 0, 1, 0, "solid")
        push!(solid, p)
    else
        return
    end

    push!(particles, p)
    i = length(particles)

    push!(cell_of_particle, (px, py))

    if p.collision == 1
        if !haskey(id_grid, (px, py))
            id_grid[(px, py)] = Int[]
        end
        push!(id_grid[(px, py)], i)
    end
end