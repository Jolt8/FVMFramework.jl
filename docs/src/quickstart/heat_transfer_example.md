

# Getting Started with FVMFramework

This tutorial will teach you how to use the tools provided by FVMFramework to solve a 3D heat transfer problem

In this tutorial, we will be solving the following heat transfer equation:

```math
Q_{a, in} = -k A \frac{T_{b} - T_{a}}{\Delta d}
```
Where:
- ``Q_{a,in}`` is the heat flux going into cell a in [W]
- ``k`` is the thermal conductivity of the material in [W/(m*K)]
    - Note that k is taken as the harmonic mean of the thermal conductivities of the two cells sharing the face
- ``A`` is the area of the face shared by both cells in [m^2]
- ``T_a`` and ``T_b`` are the temperatures of cell a and cell b respectively in [K]
- ``\Delta d`` is the distance between the cell centers of the neighboring cells in [m]

The general workflow of this package is to:
- Define a simple mesh using Ferrite.jl or import a mesh using FerriteGmsh.jl
- Create a config struct
- Define the physics types (Solid, Fluid, etc.) that will be used to map what flux physics are applied where
- Add cached fields, special caches, and optimized parameters
- Add regions (groups of cells) and define:
    - Their initial conditions for the state variables (u)
    - How cached variables will be updated 
    - What boundary conditions that region has
    - Other internal physics that happens inside the cell
- Add patches (groups of faces) and define:
    - What boundary conditions that patch has
    - Other unique physics that happen at the faces included in the patch
- Create flux physics functions that describe how fluxes are calculated between the faces of the different cell types
- Create a function that maps each combination of cell types to its respective flux physics function (using the previously defined physics types)
- And finally, call `finish_fvm_config` to create the internal data structures needed to solve the system

Moving on...

First, we must import all the packages necessary for this example to work

```julia
using Revise

using FVMFramework

using Ferrite
using FerriteGmsh
using OrdinaryDiffEq
using SparseArrays
using ComponentArrays
using NonlinearSolve
using Sparspak
import SparseConnectivityTracer, ADTypes
using ILUZero
using StaticArrays

using Unitful
```

To get started, we must first generate a grid for which the simulation will be performed on
Do note that with FerriteGmsh, custom meshes can be imported with either Hexahedral or Tetrahedral cells
```julia
grid_dimensions = (100, 10, 10)

#since base SI is used, both of these are in meters
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0)) 

grid = generate_grid(Hexahedron, grid_dimensions, left, right)
```

Next, we define our cell sets which will later form our regions and patches
```julia
cell_half_dist = (right[1] / grid_dimensions[1]) / 2
addcellset!(grid, "copper", x -> x[1] <= (left[1] + (right[1] / 2) + cell_half_dist))
addcellset!(grid, "steel", x -> x[1] >= left[1] + (right[1] / 2))
addfacetset!(grid, "bad_weld", x -> x[1] == (right[1] / 2))
```

Now we define a ComponentVector holding a vector of zeros with length equal to the amount of cells in our grid
Note that this does not define the initial conditions of the system, but rather acts as a template for what the initial conditions will look like.
Initial conditions will be defined later in `add_region!()`
```julia
u_proto = ComponentVector(
    temp = zeros(length(grid.cells))u"K",
)
```


Now we create the config to store all of our setup data and additionally define our physics types which will later determine which flux functions are used for each cell-type combination
```julia
config = create_fvm_config(grid, u_proto)

struct Solid <: AbstractPhysics end
```


For those coming fom other finite volume packages, this is likely one of the hardest things to get used to
Now, since FVMFramework prefers to cache certain values that could technically be done in a single step, we must define some cached values that will be used for physics
In this simple example, we will define a cached value for the heat in each cells
Important to note, when a cached value is define here, it ends up in both the du (state variable derivative) and the u (state variable) vectors
This means that we will have access to both `du.heat` (change in heat over time in [W]) and `u.heat` (heat present in a cell in [J])
Note that the derivatives of values are over time so `u.heat[cell_id]` is in Joules while `du.heat[cell_id]` is in Joules / second, i.e. Watts.
Also, if you're coming from another finite volume method package, you may be used to writing flux functions that update both cell_a and cell_b.
In FVMFramework, however, we instead apply the flux to each cell for each face it has instead of having to track which cell is cell_a and which is cell_b. 

#TODO: the explanation of why we use single sided flux functions needs to be improved

We define heat here because it means that we can cache the heat moving in to each cell through flux calculations and then cap the rate of heat change to a rate of temperature change (see `cap_heat_flux_to_temp_change!` below)
```julia
add_setup_syms!(
    config;
    cache_syms_and_units = (heat = u"J",),
    special_caches = ComponentVector(),
    optimized_parameters = ComponentVector()
)

#cap_heat_flux_to_temp_change! is a function already included in this framework,
#we're just showing it here to show the purpose of cap functions

function cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    #K/s = J/s / (m^3 * kg/m^3 * J/(kg*K))
    du.temp[cell_id] += du.heat[cell_id] / (vol * u.rho[cell_id] * u.cp[cell_id])
end
```

If you're wondering what the `special_caches` are for they are for more complicated caches such as but not limited to:
- values stored per face of a cell (ex. `du.mass_face` (mass flow rate across a cell's face in [kg/s]))
- caches that are stored per species (ex. `du.molar_concentrations` (change in molar concentrations for each species over time in [mol/s]))
- anything else that requires more work to calculate



Now we will add the properties for the "copper" region.

Its state variables is set to a temp of 270°C
And its properties are set to:
- A thermal conductivity (`k`) of 237 W/(m*K)
- A density (`rho`) of 2700 kg/m^3
- A specific heat capacity (`cp`) of 921 J/(kg*K)

The region is given the type `Solid()`, so will later need to define what flux physics happens between cells of this type
Furthermore, we will also add a boundary condition that continuously adds heat to the "copper" region. This is another area where FVMFramework differs from other finite volume method packages. Instead of explicitly naming boundary conditions (Dirichlet, Neumann, Robin, etc.) we define boundary conditions as functions that are applied to regions of cells or patches of faces.

Note, since these properties are defined in the properties area of a region, it means that they cannot be modified during the solve process

```julia
n_copper_cells = length(getcellset(grid, "copper"))
wattage_applied_to_region = ustrip(upreferred(10000.0u"W"))
#since this will be used in the solver, it cannot have units and must be in base SI
#if you want to ensure this you can do ustrip(upreferred(100u"W"))

per_cell_heating = wattage_applied_to_region / n_copper_cells

add_region!(
    config, "copper";
    type = Solid(),
    initial_conditions = ComponentVector(
        temp = 270.0u"°C",
    ),
    properties = ComponentVector(
        k = 237.0u"W/(m*K)", 
        rho = 2700.0u"kg/m^3",
        cp = 921.0u"J/(kg*K)",
        #if you want, you can also make per_cell_heating a property
        #per_cell_heating = (1000.0u"W" / n_copper_cells)
    ),
    property_update_function = function copper_property_update!(du, u, p, t, cell_id, vol, system)

    end,
    region_function =
    function heat_transfer!(du, u, cell_id, vol)
        du.heat[cell_id] += per_cell_heating

        #the alternative is this:
        #du.heat[cell_id] += u.per_cell_heating[cell_id]

        cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    end
)
```

Now we will add the properties for the "steel" region
Its state variables is set to a temp of 350°C
And its properties are set to:
- A thermal conductivity (`k`) of 123 W/(m*K)
- A density (`rho`) of 7800 kg/m^3
- A specific heat capacity (`cp`) of 450 J/(kg*K)
```julia
add_region!(
    config, "steel";
    type = Solid(),
    initial_conditions = ComponentVector(
        temp = 350.0u"°C",
    ),
    properties = ComponentVector(
        k = 123.0u"W/(m*K)",
        rho = 7800.0u"kg/m^3",
        cp = 450.0u"J/(kg*K)",
    ),
    property_update_function = 
    function steel_property_update!(du, u, p, t, cell_id, vol, system)
        
    end,
    region_function =
    function heat_transfer!(du, u, p, t, cell_id, vol)
        cap_heat_flux_to_temp_change!(du, u, cell_id, vol)
    end
)
```

Now, we define a flux function for heat diffusion between two cells
Notice that the syntax very closely resembles the original equation:
```math
Q_{a, in} = -k A \frac{T_{b} - T_{a}}{\Delta d}
```

Note that the arguments used by this function are expected to be used by all flux functions
**IMPORTANT:** The norm here (accesed by `norm[cell_id][face_idx]`) points from the cell center of the `idx_a` cell to the cell center of the `idx_b` cell
```julia
function heat_diffusion!(
    du, u, p, t,
    idx_a, idx_b, face_idx,
    area, norm, dist,
    vol_a, vol_b
)
    k_effective = 2 * u.k[idx_a] * u.k[idx_b] / (u.k[idx_a] + u.k[idx_b])

    grad_T = (u.temp[idx_b] - u.temp[idx_a]) / dist

    du.heat[idx_a] -= -k_effective * grad_T * area
end
```

While you could make heat_diffusion! use cell_volumes instead and ues it instead of this equation, it's preferrable to always wrap them in a function that looks like the one below
This allows you to add additional flux physics to the same region-region interface in the future
```julia
function solid_solid_flux!(
    du, u, p, t, 
    idx_a, idx_b, face_idx,
    cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances,
    cell_volumes
)
    heat_diffusion!(
        du, u, p, t,
        idx_a, idx_b, face_idx,
        cell_neighbor_areas[idx_a][face_idx], cell_neighbor_normals[idx_a][face_idx], cell_neighbor_distances[idx_a][face_idx],
        cell_volumes[idx_a], cell_volumes[idx_b]
    )
end
```

Now we define a function that will determine which flux functions are used for each region-region interface determined by its assigned type in `add_region!()`
```julia
function connection_map_function(type_a, type_b)
    typeof(type_a) <: Solid && typeof(type_b) <: Solid && return solid_solid_flux!
end
```

To show how patch functions work, we'll simulate a bad weld between the copper and steel regions.
Note that properties applied by the `add_patch!()` function are applied to cells that have a face contained in the facetset used to define the patch
Typically these would be defined before you've defined your flux physics, but I did this afterwards so you could see the regular heat transfer function first

```julia
add_patch!(
    config, "bad_weld";
    properties = ComponentVector(
        weld_k = 10.0u"W/(m*K)"
    ),
    patch_function = function bad_weld!(du, u, p, t, idx_a, idx_b, face_idx, cell_neighbor_areas, cell_neighbor_normals, cell_neighbor_distances, cell_volumes)
        grad_T = (u.temp[idx_b] - u.temp[idx_a]) / cell_neighbor_distances[idx_a][face_idx]

        du.heat[idx_a] -= -u.weld_k[idx_a] * grad_T * cell_neighbor_areas[idx_a][face_idx]
    end
)
```

Now we compile all the information we've provided to create the internal data structures needed to solve the system
```julia
du0_vec, u0_vec, geo, system = finish_fvm_config(config, connection_map_function, check_units = false);
```

Now we define a function that will be passed into the ODE problem
It determines the order of operations that will be performed during each time step

Importantly, if you wanted to overwrite a flux that would normally happen between two cell types,
you could create an `add_patch!()` that **overwrites** something like `du.mass_face` by doing something like `du.mass_face = 0.0`.
However, if you wanted to apply a flux to a patch of faces that **modifies** a normally happening flux you could write: `du.mass_face -= 0.1`

```julia
function solve_system!(du, u, p, t, geo, system)
    solve_connection_groups!(du, u, p, t, geo, system)
    solve_patch_groups!(du, u, p, t, geo, system)
    solve_region_groups!(du, u, p, t, geo, system)
end
```

Create a closure self-containing solve_system!, geo, and system into a single function to be used by the ODE solver
```julia
f_closure = (du, u, p, t) -> fvm_operator!(du, u, p, t, solve_system!, geo, system)
```


Also define a p_guess (must be a vector)
```julia
p_guess = [0.0]
```

Now that we have a fully defined ODE problem, we can use any ODE solver from `OrdinaryDiffEq.jl` to solve the system!

Let's start with an explicit solve:
```julia
test_prob = ODEProblem(f_closure, u0_vec, (0.0, 100.0), p_guess)
@time sol = solve(test_prob, Tsit5(), tspan = (0.0, 100.0), callback = approximate_time_to_finish_cb)
```

Since this framework has compatibility for ForwardDiff, compatibility with SparseConnectivityTracer is also included
This allows us to create an accurate sparsity map for the jacobian with minimal effort
Furthermore, compared to symbolic sparsity detection, it's much faster for larger systems
For smaller problems, you can skip this, but explicitly providng the sparsity pattern to the solver makes this siframework scale to hundreds of thousands of cells
```julia
detector = SparseConnectivityTracer.TracerLocalSparsityDetector()

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure(du, u, p_guess, 0.0), du0_vec, u0_vec, detector
)
```

Now we'll create the implicit problem
```julia
f_closure_implicit = ODEFunction(f_closure, jac_prototype = float.(jac_sparsity))

implicit_prob = ODEProblem(f_closure_implicit, u0_vec, (0.0, 100000.0), p_guess)
```

**IMPORTANT:** Here's a list of solvers that I've had great success with over the time I've spent working on this framework
The list is roughly ordered from best performing on large systems to best performing on small systems, but this can change depending on your specific problem
Usually, it's just best to test them out and see which one works best for you!

However, I've found that `Tsit5()` is best for explicit problems and `FBDF()` is the best for implicit problems

```julia
desired_steps = 300
save_interval = implicit_prob.tspan[2] / desired_steps

@time sol = solve(
    implicit_prob, 
    FBDF(linsolve = KrylovJL_GMRES(), precs = algebraicmultigrid, concrete_jac = true), 
    callback = approximate_time_to_finish_cb, 
    saveat = save_interval
) #for 1 000 000+ cells

@time sol = solve(
    implicit_prob, 
    FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), 
    callback = approximate_time_to_finish_cb, 
    saveat = save_interval
) #for 100 000 - 1 000 000 cells

@time sol = solve(
    implicit_prob, 
    FBDF(linsolve = SparspakFactorization()), 
    callback = approximate_time_to_finish_cb, 
    saveat = save_interval
) #for 1000 - 100 000 cells

@time sol = solve(
    implicit_prob, 
    FBDF(linsolve = KLUFactorization(), precs = iluzero, concrete_jac = true), 
    callback = approximate_time_to_finish_cb, 
    saveat = save_interval
) #usually slower than `SparspakFactorization()` but sometimes it's faster, 
#it really depends on the problem
```
**Note:** `SparspakFactorization()` is the only solver that works with `ForwardDiff`. **However**, for gradient-free optimization purposes, `KrylovJL_GMRES()` works just fine on larger systems. Then, a very slow `SparspakFactorization()` solve can be used to converge on the true optimal parameters using `ForwardDiff`.

Also, note the use of the `approximate_time_to_finish_cb` callback which displays time the simulation is at and an estimated time to finish. Trust me when I tell you that this will save you many headaches and will be invaluable for your sanity in the future! Furthermore, since it forces a console print every time the solver finishes a step, it makes Ctrl + C (stop program) almost always work without issue.

In addition, we can also use hybrid methods which switch between explicit and implicit methods:
These are usually the fastest
```julia
@time sol = solve(implicit_prob, AutoTsit5(FBDF(linsolve = SparspakFactorization())), callback = approximate_time_to_finish_cb, saveat = save_interval)
```

Now, we can finally save our solution to a vtk file to be viewed in Paraview!
```julia
record_sol = true

sim_file = @__FILE__

du_named, u_named = regenerate_fvm_state(sol, system, solve_system!, config.geo, p_guess)

output_dir = "C:\\Users\\wille\\OneDrive\\Desktop\\Julia_cfd_output_files"

#output_dir = joinpath(@__DIR__, "vtk_output")

if record_sol == true
    sol_to_vtk(sol, du_named, u_named, grid, config.geo, sim_file, output_dir; include_zeros_fields = false)
end
```

One more thing, it is highly recommended to create a new file on your desktop and make it your new `output_dir`. 
This makes it extremely easy to view solver results by just going to your desktop, going to that directory, going to the newest folder and viewing it in paraview

<video src="video.mp4" autoplay loop muted></video>