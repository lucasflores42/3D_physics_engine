# -----------------------------------------------------------------------------
#                           Parameters
# ----------------------------------------------------------------------------- 
const smoothing_length = 0.2
const liquid_target_density = 0.4
const liquid_stiff_coef = 0.1
const liquid_viscosity_coef = 0.1
const sph_cell_range = Int(ceil(3 * smoothing_length / grid_size))
 
 
 function kernel(r)
    q = r / smoothing_length
    if q <= 1.0
        return (1.0 - 1.5*q*q + 0.75*q*q*q) / (π * smoothing_length^2)
    elseif q <= 2.0
        return 0.25 * (2.0 - q)^3 / (π * smoothing_length^2)
    else
        return 0.0
    end
end

function kernel_gradient(r, r_vec)
     q = r / smoothing_length

    if q <= 1.0
        factor = (-3.0 + 2.25*q) / (π * smoothing_length^4)
        return factor * r_vec
    elseif q <= 2.0
        factor = -0.75 * (2.0 - q)^2 / (π * smoothing_length^4 * q)
        return factor * r_vec
    else
        return @SVector zeros(2)
    end
end

function calculate_density_pressure!(p, particles, id_grid)

    p.density = 0.0

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
                if p2.material != "liquid"
                    continue
                end
                r_vec = p.position - p2.position
                r = norm(r_vec)
                p.density += p2.mass * kernel(r)
            end
        end
    end

    p.pressure = liquid_stiff_coef * ((p.density/liquid_target_density)^7 - 1)
end

function pressure_gradient(p, p2, r, r_vec)

    kernel_grad = kernel_gradient(r, r_vec)
    grad_pressure = @SVector zeros(2)

    d1 = max(p.density, 1.0)
    d2 = max(p2.density, 1.0)

    pressure_term = (p.pressure / d1^2 + p2.pressure / d2^2)
    grad_pressure += p2.mass * pressure_term * kernel_grad

    return grad_pressure
end

function viscosity_laplacian(p, p2, r, r_vec)

    kernel_grad = kernel_gradient(r, r_vec)

    laplacian_velocity = @SVector zeros(2)
    v_ij = p.velocity - p2.velocity
    dot_r_grad = dot(r_vec, kernel_grad)
    denominator = dot(r_vec, r_vec) + 0.01 * smoothing_length^2
    
    if denominator != 0
        laplacian_velocity += 2.0 * (p2.mass / p2.density) * 
                            v_ij * (dot_r_grad / denominator)
    end

    return laplacian_velocity
end

    


