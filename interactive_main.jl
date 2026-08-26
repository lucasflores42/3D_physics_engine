# ============================================================================
# Interactive particle-placement sandbox for the 2D physics engine.
# Uses GLMakie for the window/UI (buttons, mouse events, live-updating plot)
# and reuses the physics engine defined in 2D_physics_engine.jl.
# ============================================================================

using Plots, LinearAlgebra, StaticArrays, GLMakie
# NOTE: Plots and GLMakie both export some of the same function names
# (heatmap!, xlims!, ylims!). Whenever we call one of those, we prefix it
# with "GLMakie." to make sure Julia uses the right one and doesn't error
# with "UndefVarError".

include("2D_physics_engine.jl")   # brings in create_scene, simulation_step,
                                   # build_material_grid, spawn_particle!, etc.

Makie.inline!(false)   # forces GLMakie to open a real desktop window instead
                        # of trying to embed the plot inline (e.g. in a notebook)

# ----------------------------------------------------------------------------
# Build the initial simulation state, exactly like the non-interactive main()
# ----------------------------------------------------------------------------
particles, liquid, gas, powder, solid, rigidbodies, softbodies = create_scene()
id_grid, cell_of_particle = init_grids(particles)

# ----------------------------------------------------------------------------
# Figure = the whole window. Everything else (axis, buttons) is placed
# inside it using a grid layout, addressed like fig[row, column].
# ----------------------------------------------------------------------------
fig = Figure(size = (1400, 850))

# A GridLayout is a "sub-grid" we can pack multiple buttons into, stacked
# vertically, all living in column 1 of the outer figure.
layout_botoes = fig[1, 1] = GridLayout(tellwidth = true, tellheight = false)

# The Axis is where the actual simulation gets drawn (column 2). Its title
# doubles as the "which material is currently selected" indicator - we
# update ax.title directly from the button callbacks below.
ax = Axis(fig[1, 2])

# Turns off GLMakie's default "drag a rectangle to zoom" interaction, since
# we want plain left-click/drag to mean "place a particle", not "zoom".
deregister_interaction!(ax, :rectanglezoom)

# Set the plot's coordinate range to match the simulation's actual size
# (in simulation units, not pixels) - so clicking at a given screen spot
# maps to a sensible simulation coordinate.
GLMakie.xlims!(ax, 0, box_size_x)
GLMakie.ylims!(ax, 0, box_size_y)

# material_grid_obs holds the current "what material is in each grid cell"
# matrix. It's an Observable - reassigning material_grid_obs[] later
# automatically refreshes the heatmap on screen, no manual redraw needed.
material_grid_obs = Observable(build_material_grid(particles, id_grid))

# Categorical color scale: index 0=white(empty), 1=brown(solid), 2=blue(liquid),
# 3=gray(gas), 4=orange(powder) - matches material_code(...) in the engine.
colors = cgrad([:white, :brown, :blue, :gray, :orange], 5, categorical=true)
GLMakie.heatmap!(ax, material_grid_obs, colorrange = (0,4), colormap = colors)

# ----------------------------------------------------------------------------
# Buttons - one per material, stacked in the left column (layout_botoes).
# ----------------------------------------------------------------------------
btn1 = Button(layout_botoes[1, 1], label = "Liquid", buttoncolor = :blue, labelcolor = :white)
btn2 = Button(layout_botoes[2, 1], label = "Solid", buttoncolor = :brown, labelcolor = :white)
btn3 = Button(layout_botoes[3, 1], label = "Gas", buttoncolor = :gray, labelcolor = :white)
btn4 = Button(layout_botoes[4, 1], label = "Powder", buttoncolor = :orange, labelcolor = :white)

# tipo_atual ("current type") tracks which material is selected right now.
# 0 = nothing selected yet.
tipo_atual = Observable(0)

# `on(btn1.clicks) do _ ... end` registers a callback: whenever btn1 is
# clicked, this code runs. The `_` just means "we don't care about the
# click event's details, only that it happened". Each callback also updates
# ax.title so you can see which material is currently active.
on(btn1.clicks) do _; tipo_atual[] = 1; end
on(btn2.clicks) do _; tipo_atual[] = 2; end
on(btn3.clicks) do _; tipo_atual[] = 3; end
on(btn4.clicks) do _; tipo_atual[] = 4; end

# Maps tipo_atual's numeric code to the actual material name string that
# spawn_particle! expects.
material_names = Dict(1 => "liquid", 2 => "solid", 3 => "gas", 4 => "powder")

# Tracks whether the left mouse button is currently held down, so we can
# support "click and drag to paint" instead of just single clicks.
mouse_pressionado = Observable(false)

# Called whenever we might want to place a particle (on click, or on mouse
# move while the button is held).
function adicionar_particula_se_ativo()
    if tipo_atual[] != 0   # only place something if a material is selected

        # mouseposition(ax.scene) converts the mouse's screen position into
        # the axis's data coordinates (i.e. simulation x/y), automatically.
        posicao = mouseposition(ax.scene)

        # Only place if the click is actually inside the simulation bounds.
        if 0 <= posicao[1] <= box_size_x && 0 <= posicao[2] <= box_size_y

            # Convert continuous simulation coordinates to the engine's
            # integer grid cell indices (same convention used everywhere
            # else in the engine).
            px = Int(floor(posicao[1]/grid_size)) + 1
            py = Int(floor(posicao[2]/grid_size)) + 1

            if 1 <= px <= pixel_size_x && 1 <= py <= pixel_size_y
                spawn_particle!(particles, liquid, gas, powder, solid, id_grid, cell_of_particle,
                                 px, py, material_names[tipo_atual[]])
            end
        end
    end
end

# Mouse button press/release handling.
on(ax.scene.events.mousebutton) do event
    if event.button == Mouse.left
        if event.action == Mouse.press
            mouse_pressionado[] = true
            adicionar_particula_se_ativo()   # place one immediately on click
        elseif event.action == Mouse.release
            mouse_pressionado[] = false
        end
    end
end

# Mouse movement handling - if the button is still held, keep placing
# particles as the mouse moves (like painting).
on(ax.scene.events.mouseposition) do _
    if mouse_pressionado[]
        adicionar_particula_se_ativo()
    end
end

# Actually opens the window and returns a handle to it (`screen`), which we
# use below to check whether the window is still open.
screen = display(fig)

# ----------------------------------------------------------------------------
# Physics loop. GLMakie's own event loop (handling clicks, redraws, etc.)
# runs continuously in the background while the window is open. To step our
# simulation forward WITHOUT blocking that event loop, we run our own loop
# inside @async - this lets both "loops" interleave, as long as we
# periodically yield control back with sleep()/yield().
# ----------------------------------------------------------------------------
t = 0.0
step = 0
@async begin
    while t < tmax && screen.window_open[]
        step_time = @elapsed simulation_step(particles, liquid, gas, powder, solid, rigidbodies, softbodies, id_grid, cell_of_particle)

        global step += 1

        render_time = 0.0
        if step % 3 == 0
            render_time = @elapsed (material_grid_obs[] = build_material_grid(particles, id_grid))
        end

        global t += dt
        println("t = $(round(t, digits=2))s | step: $(round(step_time*1000, digits=2))ms | render: $(round(render_time*1000, digits=2))ms")

        sleep(0.001)
        yield()
    end
end
# If running this file as a script (not from the REPL), keep the process
# alive until the window is closed - otherwise Julia would exit immediately
# after starting the async loop.
if !isinteractive()
    wait(screen)
end