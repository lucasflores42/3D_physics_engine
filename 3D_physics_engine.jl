#import Pkg
#Pkg.add(["StaticArrays", "Plots", "LinearAlgebra", "GLMakie])
using Plots, LinearAlgebra, StaticArrays #, GLMakie

# -----------------------------------------------------------------------------
#                           Parameters
# ----------------------------------------------------------------------------- 
# 256x256x256, total 16,777,216 voxels.
const grid_size = 1.0
const voxel_size_x = 100
const voxel_size_y = 100 
const voxel_size_z = 100 
const box_size_x = voxel_size_x * grid_size
const box_size_y = voxel_size_y * grid_size
const box_size_z = voxel_size_z * grid_size

const tmax = 1000.0
const dt = 0.01

include("sph_functions.jl")
include("rigidbody_functions.jl")
include("softbody_functions.jl")
include("particle_functions.jl")
include("collision_functions.jl")
include("other_functions.jl")
include("material_functions.jl")


# -----------------------------------------------------------------------------
#                           Create Scene
# ----------------------------------------------------------------------------- 
function create_scene()

    particles = Union{liquid_struct, solid_struct, gas_struct, powder_struct}[]
    liquid = liquid_struct[]
    solid = solid_struct[]
    gas = gas_struct[]
    powder = powder_struct[]
    rigidbodies = rigidbody_struct[]
    softbodies = softbody_struct[]

    # boundary of world
    for i in 1:voxel_size_x
        for j in 1:voxel_size_y
            for k in 1:voxel_size_z

                if k == 1

                    p = solid_struct(
                        length(particles) + 1,
                        SVector{3, Float64}(i, j, k),
                        SVector{3, Float64}(0.0, 0.0, 0.0),
                        SVector{3, Float64}(0.0, 0.0, 0.0),
                        grid_size / 2,
                        10000.0,

                        0,              # rigidbody
                        0,

                        0,              # active
                        1,              # collision
                        0,              # gravity

                        "solid"         # material
                    )
                    push!(solid, p)
                    push!(particles, p)

                end
            end
        end
    end

    # some liquid
    for i in 1:10

        x = 325
        y = 55 + 10*rand()
        z = 55 + 10*rand()

        p = liquid_struct(length(particles)+1, SVector(x,y,z), @SVector(zeros(3)), @SVector(zeros(3)),
                       grid_size/2, 0.1, 0, 0,
                       0.4, 0.0, 0.4, 0.1, 0.1,   # density, pressure, target_density, stiff_coef, viscosity_coef
                       1, 1, 1, 1, 
                       1, "liquid")
        push!(liquid, p)
        push!(particles, p)
    end

    create_sphere!(particles, rigidbodies, length(rigidbodies) + 1, [40.0, 100.0, 100.0], [0.0, 10.0, 0.0], [0.0, 0.0, 0.0], 10)
    create_sphere!(particles, rigidbodies, length(rigidbodies) + 1, [160.0, 100.0, 100.0], [0.0, -10.0, 0.0], [0.0, 0.0, 0.0], 10)

    return particles, liquid, gas, powder, solid, rigidbodies, softbodies
end

# -----------------------------------------------------------------------------
#                           Simulation step
# ----------------------------------------------------------------------------- 
function simulation_step(particles, liquid, gas, powder, solid, rigidbodies, softbodies, id_grid, cell_of_particle)

    particle_physics(particles, liquid, gas, powder, solid, id_grid, cell_of_particle)
    rigidbody_physics(particles, rigidbodies)
    softbody_physics(particles, softbodies)

    collision_physics!(particles, rigidbodies, powder, liquid, gas, id_grid, cell_of_particle)

    update_grids!(particles, id_grid, cell_of_particle)
end

# -----------------------------------------------------------------------------
#                           Visualization
# ----------------------------------------------------------------------------- 
function visualization(particles, id_grid, step)

    material_grid = build_material_grid(particles, id_grid)

    xs = Int[]
    ys = Int[]
    zs = Int[]
    cs = Symbol[]

    for x in 1:voxel_size_x
        for y in 1:voxel_size_y
            for z in 1:voxel_size_z
                v = material_grid[x, y, z]
                if v == 0
                    continue
                end

                push!(xs, x)
                push!(ys, y)
                push!(zs, z)

                color = if v == 1
                    :gray
                elseif v == 2
                    :blue
                elseif v == 3
                    :orange
                elseif v == 4
                    :brown
                elseif v == 5
                    :lightblue
                else
                    :white
                end

                push!(cs, color)
            end
        end
    end

    plt = scatter3d(
        xs, ys, zs;
        markercolor = cs,
        markersize = 1.5,
        xlabel = "X", ylabel = "Y", zlabel = "Z",
        title = "Time $(round(step, digits=2))s",
        xlim = (0, voxel_size_x),
        ylim = (0, voxel_size_y),
        zlim = (0, voxel_size_z),
        aspect_ratio = :equal,
        camera = (30, 35),
        legend = false
    )

    return plt
end

# -----------------------------------------------------------------------------
#                           Main Simulation
# ----------------------------------------------------------------------------- 
function main()
    t = 0.0
    step = 0

    particles, liquid, gas, powder, solid, rigidbodies, softbodies = create_scene()
    id_grid, cell_of_particle = init_grids(particles)

    while t < tmax
        step += 1

        step_time = @elapsed simulation_step(
            particles, liquid, gas, powder, solid, rigidbodies, softbodies,
            id_grid, cell_of_particle
        )

        if step % 10 == 0
            render_time = @elapsed begin
                plt = visualization(particles, id_grid, t)
                display(plt)
            end
            println("t = $(round(t, digits=2))s | step: $(round(step_time * 1000, digits=2))ms | render: $(round(render_time * 1000, digits=2))ms")
        else
            println("t = $(round(t, digits=2))s | step: $(round(step_time * 1000, digits=2))ms")
        end

        t += dt
    end
end

main()
