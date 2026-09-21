function material_transformation(p1,p2)

    if p1.material == "powder" && p2.material == "liquid"

        new_gas = gas_struct(length(particles)+1, p1.position, SVector(rand(),0.0), @SVector(zeros(2)),
                        grid_size/2, 0.1, 
                        0, 0, 
                        1, 1, 1, 1, 
                        0, 300, "gas")
        transform_particle!(particles, gas, id_grid, cell_of_particle, p1, p2, new_gas)
        return
    elseif p1.material == "liquid" && p2.material == "powder"

        new_gas = gas_struct(length(particles)+1, p2.position, SVector(rand(),0.0), @SVector(zeros(2)),
                        grid_size/2, 0.1, 
                        0, 0, 
                        1, 1, 1, 1, 
                        0, 300, "gas")
        transform_particle!(particles, gas, id_grid, cell_of_particle, p1, p2, new_gas)
        return
    end
end