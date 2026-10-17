
# Here's some other guidelines you should follow

In general, try to stay pretty consistent with my solver architecture. If there's some data that cannot be easily accessed, first try to put it into the system.additional_data NamedTuple since that can store any arbitrary data in a safe manner. If something truly can't be added through already existing methods, data structures, etc. etc., just make the modification and make a note about it and then I'll review it later.

Do note that once you put something in the additional_data NamedTuple, you there's a function overload that allows system.field to work for any field you add to the additional_data NamedTuple. And try to do this because I want this to be convention.

The following repository-specific instructions override generic implementation preferences elsewhere in this prompt. Do not refactor existing FVMFramework conventions merely because another implementation would be more conventional. However, if there's something that you truly believe would make the code more understandable, more efficient, or easier to work with, definitely suggest and implement the change while making a note of it and then I'll review it later.

Also, before creating a new data structure, helper function, storage mechanism, or geometry-access pattern, search the repository for an existing equivalent and use it where practical.

Do not perform speculative micro-optimizations outside hot paths unless profiling or benchmarks justify them.

Prefer readable code, but try not to sacrifice too much performance for readability. That being said, don't be too strict about performance requirements unless you notice a big issue.

Also, optimization outside of the actual solver hot loop such as allocations in the main script do not matter at all. 

One thing that you should definitely know about is this very important piece of syntax for getting good performance iterating through nested fields of ComponentArrays.

Basically, we always try to use this looping method when we need to iterate through the fields of a ComponentArray:

"""
    foreach_field_at!(f, cell_id::Int, groups::Vararg{Any, N}) where {N}

Iterates over the elements of each field at a specific cell_id

# Example:
```julia
foreach_field_at!(cell_id, du.mass_fractions, du.molar_concentrations) do species, du_mass_fractions, du_molar_concentrations
    du_mass_fractions[species] += 1.0
    du_molar_concentrations[species] += 1.0
end
```

The reason this function exists is that this causes a ton of allocations and dynamic dispatch:
```julia
for name in propertynames(du.mass_fractions[1])
    getproperty(du.mass_fractions, name)[cell_id] += 1.0 
end
```

However, with how this function works, there are sometimes cases where it does not behave as you would expect.
For example, let's say you have these two vectors and wanted to iterate through both of them:
```julia
mass_fractions = ComponentVector(
    methane = 0.5,
    water = 0.5
)

elemental_compositions = ComponentVector(
    methane = (
        C = 1,
        H = 4,
        O = 0
    ),
    water = (
        C = 0,
        H = 2,
        O = 1
    )
)
```
\n
If you just do:
```julia
foreach_field_at!(cell_id, u.mass_fractions, u.elemental_compositions) do species, u_mass_fractions, u_elemental_compositions
    u_mass_fractions[species[cell_id]] += 1.0
    u_elemental_compositions[species[cell_id]] += 1.0
end
```
**species** will just be [1], [2] etc.
\n
This will not index [elemental_compositions] properly because it needs be indexed as [1:3], [4:6]
thus, **elemental_compositions[species]** will just return the value for **elemental_compositions.methane.C** and nothing else

To make this this doesn't happen, do this:
```julia
mass_fractions_species_idx = 1

foreach_field_at!(cell_id, u.elemental_compositions) do species, u_elemental_compositions
    view(u.mass_fractions, cell_id)[mass_fractions_species_idx] += 1.0
    u_elemental_compositions[species[cell_id]] += 1.0
    mass_fractions_species_idx += 1
end
```
While this is not the most ideal, it's pretty much the best option avaliable without introducing even more bloated looping functions
"""

Also, do note that my current VirtualFVMArray datastructure does not allow for completely arbitrary nested vectors. I'd say the most you'd be able to get away with is this:
```julia
grad_density = zeros(n_cells, 3)u"kg/m^4",
grad_vel_u = zeros(n_cells, 3)u"m/(s*m)",
grad_vel_v = zeros(n_cells, 3)u"m/(s*m)",
grad_vel_w = zeros(n_cells, 3)u"m/(s*m)",
grad_temperature = zeros(n_cells, 3)u"K/m",


vel_u_face = zeros(n_cells, n_faces)u"m/s",
vel_v_face = zeros(n_cells, n_faces)u"m/s",
vel_w_face = zeros(n_cells, n_faces)u"m/s",
temperature_face = zeros(n_cells, n_faces)u"J/m^3"
```

These will all be accessed like matricies like u.grad_density[cell_id, 3]
Right now, stuff like u.grad_density[1][3] is not supported and will lead to performance issues so avoid doing that.

However you can have arbitrarily nested ComponentVectors like:
```julia
species_atom_properties = ComponentVector(
    name_1 = (
        nested_name_2 = (
            nested_nested_name_3 = (
                nested_nested_nested_name_4 = 1,
                nested_nested_nested_name_5 = 12.011
            ),
            nested_nested_name_6 = (
                nested_nested_nested_name_7 = 4,
                nested_nested_nested_name_8 = 1.008
            ),
        )
    )
)
```

Also, just try to write code that's pretty human readable, in general I like more descriptive names. For example, I recently switched to u.density instead of u.rho, I like u.temperature instead of u.T, I like u.pressure instead of u.p. Although I usually use u.mu, I think u.dynamic_viscosity is more readable. However there are some exceptions, I like u.vel_u instead of u.velocity_u, etc. Try to be consistent with how I name things throughout the code.

Also, I like whitespace between kwargs, I like `K = 3.0` instead of `K=3.0`. Also `test <: Zoop` instead of `test<:Zoop`. But no whitespace for brackets or function arguments: `function foo(bar, baz)` instead of `function foo( bar, baz )`.

Also, don't always make physics calculations a multi-line thing. 

Example of times when multi-lines were good:
```julia 
normal_velocity_a =
        vel_u_a * cell_face_normal[1] +
        vel_v_a * cell_face_normal[2] +
        vel_w_a * cell_face_normal[3]

F_density, 
F_momentum_density_u,
F_momentum_density_v,
F_momentum_density_w,
F_volumetric_energy = 
hllc_flux(
    density_a,
    momentum_density_u_a,
    momentum_density_v_a,
    momentum_density_w_a,
    volumetric_energy_a,
    (u.cp[idx_a] / u.cv[idx_a]),

    density_b,
    momentum_density_u_b,
    momentum_density_v_b,
    momentum_density_w_b,
    volumetric_energy_b,
    (u.cp[idx_b] / u.cv[idx_b]),

    face_normal_a,
    low_mach_correction,
)

#I like this one because it clearly distinguishes that internal_energy and kinetic_energy are similar in form and units
volumetric_energy_face = density_face * (
    specific_internal_energy_face + specific_kinetic_energy_face
)
```

Example of a time when a multi-line was not needed:
```julia
corrected_normal_velocity_a =
        average_normal_velocity + velocity_jump_scaling * half_normal_velocity_jump
corrected_normal_velocity_b =
    average_normal_velocity - velocity_jump_scaling * half_normal_velocity_jump
```

OH, another thing:
Never ever use single line conditionals, always use explicit if statements.

Don't do this:
```julia
weight_power >= 0 || throw(ArgumentError("weight_power must be non-negative"))
```

Do this instead:
```julia
if weight_power < 0
    throw(ArgumentError("weight_power must be non-negative"))
end
```

Also, I know I do this a lot in my code, but I've reverted my decision. Before I always used zero(value) to make something zero, but I think just putting = 0.0 is better.

Example:
```julia
gradient_x = zero(center_value) #bad

gradient_x *= 0.0
#or
gradient_x = 0.0
```

Use this function to get geometry info on internal faces between two cells:
```julia
(
    dist,
    face_area_a, face_normal_a, face_distance_a, vol_a,
    face_area_b, face_normal_b, face_distance_b, vol_b
) = interface_geometry(geo, idx_a, face_a, idx_b, face_b)
```

Use this function to get geometry info on boundaries where there is only one cell on that face:
```julia
face_area, face_normal, face_distance, cell_volume = boundary_geometry(geo, idx, face_idx)
```

Use this function to get geometry info for a cell:
```julia
vol = cell_geometry(geo, cell_id)
```

Also, please keep in mind that units are never transferred to any variables within the code. Every time finish_fvm_config() is called, it upreferres.() and ustrips.() everything into base SI.

Speaking of which...
Here are the units used in the simulation:
- pressure = Pa
- temperature = K
- mass = kg
- mole = mol
- distance = m
- time = s
- energy = J

THAT BEING SAID!!!:
Whenever you're entering units into the simulation, you do not have to obey these conventions at all because unitful automatically converts them to the right units.
For example, you can totally use 1.0u"kJ/mol", 3.0u"inch", they will all be converted to base SI guaranteed so that should never be something you should be that concerned about.

Also, when defining mass fractions, just use 1.0u"kg/kg" for the unit instead of just 1.0.

Avoid julia -e for ad-hoc scratch/debug scripts because it has been unreliable in this development environment. The explicit Pkg.instantiate() command below is an exception.

If you need to run small scripts to test things, do not write scratch files to the project directory. 
Instead:
1. Create a scratch file in your designated artifact scratch directory
2. Use the `write_to_file` tool to create the scratch file at that absolute path.
3. Run the code in the terminal by providing the absolute path to Julia, like this: 
   `julia C:\Users\...\scratch\scratch.jl`

If you wanted to run a script in the project directory, use this:
julia --project=. examples\new_navier_stokes_test\new_navier_stokes_solver.jl

If you get a "Package not installed" error despite using `--project=.`, run this first to download missing dependencies:
julia --project=. -e "using Pkg; Pkg.instantiate()"


You can also profile the solver like this:

```julia
using Profile
println("Including solver script (this will compile and run it once)...")
# We include the script to load all the definitions and run the first (compilation) pass
include(normpath(joinpath(@__DIR__, "../../../../../OneDrive/Desktop/FVMFramework/examples/new_navier_stokes_test/new_navier_stokes_solver.jl")))
println("Running with @profile...")
Profile.clear()
@profile solve(
    implicit_prob,
    callback = callbacks,
)
println("Writing profile output to files...")
open(joinpath(@__DIR__, "prof_tree.txt"), "w") do f
    Profile.print(IOContext(f, :displaysize => (500, 500)), format=:tree)
end
open(joinpath(@__DIR__, "prof_flat.txt"), "w") do f
    Profile.print(IOContext(f, :displaysize => (500, 500)), format=:flat, sortedby=:count)
end
println("Done!")
```