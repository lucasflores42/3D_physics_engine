#import Pkg
#Pkg.add(["StaticArrays", "Plots", "LinearAlgebra", "GLMakie])
using Plots, LinearAlgebra, StaticArrays #, GLMakie

# -----------------------------------------------------------------------------
#                           Parameters
# ----------------------------------------------------------------------------- 
# 270 height and 480 width, total 129,600 pixels.
const grid_size = 1.0
const pixel_size_x = 480
const pixel_size_y = 270 
const box_size_x = pixel_size_x * grid_size
const box_size_y = pixel_size_y * grid_size

const tmax = 1000.0
const dt = 0.01

include("sph_functions.jl")
include("rigidbody_functions.jl")
include("softbody_functions.jl")
include("particle_functions.jl")
include("collision_functions.jl")
include("other_functions.jl")


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
    for i in 1:pixel_size_x
        for j in 1:pixel_size_y

            if (i >= 1 && i <=5) || (i>=pixel_size_x-4 && i<=pixel_size_x) || (j>=1 && j<=5) || (j>=pixel_size_y-4 && j<=pixel_size_y)

                pos_x = (i - 1) * grid_size + grid_size/2
                pos_y = (j - 1) * grid_size + grid_size/2

                p = solid_struct(
                    [pos_x, pos_y],           
                    [0.0, 0.0],     # velocity
                    [0.0, 0.0],     # acceleration
                    grid_size/2,    # radius
                    10000.0,      # mass

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

            u_left = 300
            u_right = 350
            u_bottom = 50
            u_top = 150
            u_thickness = 3
            
            if (i >= u_left && i <= u_left + u_thickness && j >= u_bottom && j <= u_top) ||      # Left wall
            (i >= u_right - u_thickness && i <= u_right && j >= u_bottom && j <= u_top) ||       # Right wall
            (i >= u_left && i <= u_right && j >= u_bottom && j <= u_bottom + u_thickness)      # Bottom wall
            #(i >= u_left && i <= u_right && j >= u_top && j <= u_top + u_thickness)              # top wall

                
                pos_x = (i - 1) * grid_size + grid_size/2
                pos_y = (j - 1) * grid_size + grid_size/2

                p = solid_struct(
                    [pos_x, pos_y],           
                    [0.0, 0.0],     # velocity
                    [0.0, 0.0],     # acceleration
                    grid_size/2,    # radius
                    1000000.0,      # mass

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

    # some liquid
    for i in 1:0

        x = 305 + 40*rand()
        y = 60 + 200*rand()

        p = liquid_struct(
            [x, y],           
            [0.0, 0.0],     # velocity
            [0.0, 0.0],     # acceleration
            grid_size/2,    # radius
            0.1,            # mass

            0,              # rigidbody
            0,

            1000.0,         # density
            0.0,            # pressure

            1,              # active 
            1,              # collision
            1,              # gravity
            0,              # sph

            "liquid"        
        )
        push!(liquid, p)
        push!(particles, p)
    end

    # some gas
    for i in 1:0

        x = 325
        y = 55 + 80*rand()

        p = gas_struct(
            [x, y],           
            [0.0, 0.0],     # velocity
            [0.0, 0.0],     # acceleration
            grid_size/2,    # radius
            0.1,            # mass

            0,              # rigidbody
            0,

            1,              # active 
            1,              # collision
            1,              # gravity
            0,              # sph

            "gas"        
        )
        push!(gas, p)
        push!(particles, p)
    end

    # some powder
    for i in 1:0

        x = (1/3)*box_size_x + rand() * (1/6)*(box_size_x - 2 * grid_size)
        y = grid_size + rand() * (box_size_y - 2 * grid_size)

        p = powder_struct(
            [95,y],           
            [0.0, 0.0],     # velocity
            [0.0, 0.0],     # acceleration
            grid_size/2,    # radius
            10.0,            # mass

            0,              # rigidbody
            0,

            1,              # active 
            1,              # collision
            1,              # gravity

            "powder"        
        )
        push!(powder, p)
        push!(particles, p)
    end

    create_cube!(particles, rigidbodies, 1, [100.0, 8.0], [0.0, 0.0], [0.0],15, 3)
    create_cube!(particles, rigidbodies, 2, [100-6, 25.0], [0.0, 0.0], [0.0],2, 15)
    
    create_rope!(particles, softbodies, 1, [250.0, 180.0], 10, 0.1, grid_size)

    return particles, liquid, gas, powder, solid, rigidbodies, softbodies
end

# -----------------------------------------------------------------------------
#                           Simulation step
# ----------------------------------------------------------------------------- 
function simulation_step(particles, liquid, gas, powder, solid, rigidbodies, softbodies, id_grid, cell_of_particle)

    particle_physics(particles, liquid, gas, powder, solid, id_grid)
    rigidbody_physics(particles, rigidbodies)
    softbody_physics(particles, softbodies)

    collision_physics!(particles, rigidbodies, id_grid)

    update_grids!(particles, id_grid, cell_of_particle)
end

# -----------------------------------------------------------------------------
#                           Visualization
# ----------------------------------------------------------------------------- 
function visualization(particles, id_grid, step)

    material_grid = build_material_grid(particles, id_grid)

    colors = cgrad([:white, :brown, :blue, :gray, :orange], 5, categorical=true)

    plt = heatmap(material_grid', color=colors, clims=(0,4),
                  xlim=(0, box_size_x), ylim=(0, box_size_y),
                  title="Time $(round(step, digits=2))s",
                  xlabel="X", ylabel="Y",
                  size=(1920, 1080), aspect_ratio=:equal, legend=false)

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

        if step % 10 == 0
            print("time of plot:")
            plt = @time visualization(particles, id_grid, t)
            print("time of display:")
            @time display(plt)
        end
        print("time of step:")
        @time simulation_step(particles, liquid, gas, powder, solid, rigidbodies, softbodies, id_grid, cell_of_particle)
        t += dt

        println()
    end
end

#main()
