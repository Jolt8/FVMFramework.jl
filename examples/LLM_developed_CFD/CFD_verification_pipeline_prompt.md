You are working on my Julia finite-volume simulation package, **FVMFramework**. I want you to build an automated verification and numerical-diagnostics framework for the compressible Navier–Stokes portion of the package.

The purpose of this framework is to make it practical for LLMs to implement CFD algorithms while relying on deterministic tests, analytical solutions, convergence studies, conservation checks, Jacobian diagnostics, and solver-performance diagnostics to determine whether an implementation is actually correct.

I do **not** want the verification system to merely check that a simulation finishes or that plots "look reasonable."

I want it to answer questions such as:

- Is the implementation mathematically correct?
- Does it conserve the quantities it should conserve?
- Does it preserve uniform flow?
- Does it achieve its expected spatial order of accuracy?
- Does it produce physically admissible states?
- Did a new reconstruction method change the residual/Jacobian stencil?
- Does the supplied Jacobian or automatic differentiation agree with finite differences?
- Does an explicit integrator work when the implicit integrator fails?
- Did a change make OrdinaryDiffEq dramatically harder to solve even though the final answer is still correct?
- Can we localize failures to a particular numerical component instead of debugging an entire CFD simulation?

The eventual goal is to be able to give a future coding agent a task such as:

"Implement MUSCL reconstruction."

and then have that agent repeatedly run the verification suite until the implementation is demonstrably correct and does not introduce unacceptable solver-performance regressions.

# IMPORTANT DEVELOPMENT PHILOSOPHY

The solver implementation and the verification oracle must be as independent as practical.

Do not write tests whose only criterion is that the new implementation agrees with itself.

Good verification targets include:

- analytical solutions,
- exact conservation identities,
- known mathematical limiting behavior,
- independently evaluated finite differences,
- manufactured solutions,
- published/canonical CFD benchmark solutions,
- expected convergence order,
- fixed reference datasets.

The implementation agent may write test code, but a failed test must NOT be "fixed" by weakening tolerances, changing reference values, or modifying analytical expectations unless there is a demonstrated error in the verifier itself.

Keep expected results, convergence requirements, tolerances, and benchmark data logically separate from implementation code.

# EXISTING SOFTWARE CONTEXT

The project is written in Julia.

The solver uses a finite-volume formulation and integrates the semi-discrete PDE system through SciML / OrdinaryDiffEq.

Relevant solver infrastructure may include:

- OrdinaryDiffEq implicit methods such as FBDF or other stiff integrators,
- Krylov linear solvers such as KrylovJL_GMRES,
- ILU0 or related preconditioners,
- automatic differentiation and/or sparsity detection,
- conservative compressible variables,
- face-based numerical fluxes,
- structured or unstructured finite-volume grids,
- custom FVMFramework system/state representations.

The existing compressible solver is known to work in a relatively simple configuration.

Current or planned numerical/physical components include:

- compressible Euler/Navier–Stokes,
- HLLC,
- Thornber low-Mach reconstruction correction,
- MUSCL,
- slope limiters,
- molecular viscosity,
- heat conduction,
- species advection,
- species diffusion,
- SST k-omega turbulence,
- reacting-flow models later,
- moving meshes later.

Do not assume exact function names or data structures before inspecting the repository.

First inspect the existing architecture and adapt the verification framework to the existing conventions rather than imposing an entirely new solver architecture.


# Here's some other guidelines you should follow

In general, try to stay pretty consistent with my solver architecture. If there's some data that cannot be easily accessed, first try to put it into the system.additional_data NamedTuple since that can store any arbitrary data in a safe manner. If something truly can't be added through already existing methods, data structures, etc. etc., just make the modification and make a note about it and then I'll review it later.

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

# PRIMARY DESIGN REQUIREMENT

Create a verification system with approximately the following user-facing concept:

```julia
report = validate(
    some_component_or_configuration;
    level = :unit
)
```

or an equivalent API that fits the existing package architecture.

I eventually want to be able to do things conceptually like:

```julia
validate(HLLC())
validate(MUSCL())
validate(ThornberCorrection())
validate(SpeciesTransport())
validate(CompressibleFlow(...))
```

If these exact constructors do not exist, do not invent a major API rewrite merely to match this syntax. Create a clean verification interface appropriate to the current architecture.

Support at least three levels of validation:

```julia
:unit
:integration
:full
```

Suggested intent:

- `:unit` should run very quickly, ideally seconds or less.
- `:integration` should exercise the semi-discrete PDE plus time integration.
- `:full` may run convergence studies and larger benchmark cases.

# OUTPUT REQUIREMENTS

The verifier must generate:

1. A concise human-readable terminal report.
2. A machine-readable report suitable for feeding directly back to an LLM.
3. Detailed diagnostics for failed tests.

JSON would be a good machine-readable format unless another existing project format is more appropriate.

Example conceptual output:

```text
FVMFramework Verification Report

COMPONENT TESTS
PASS  Uniform-state preservation
PASS  Face conservation
PASS  MUSCL constant reconstruction
PASS  MUSCL linear reconstruction
PASS  Thornber high-Mach recovery

RESIDUAL TESTS
PASS  Global mass conservation
PASS  Global energy conservation
PASS  Constant-state RHS
PASS  NaN / Inf detection
PASS  State admissibility

JACOBIAN TESTS
PASS  AD vs finite-difference Jv
FAIL  Declared sparsity
      Expected stencil radius: 1
      Observed stencil radius: 2

INTEGRATION
PASS  Explicit diagnostic solve
FAIL  FBDF + GMRES

PERFORMANCE
Baseline linear iterations: 420
Current linear iterations: 8700
Regression: +1971%

RESULT: FAIL
```

Machine-readable output should identify:

- test name,
- component,
- pass/fail,
- measured value,
- expected value or tolerance,
- relative/absolute error,
- useful diagnostic context,
- solver statistics where relevant.

Do not make an LLM parse thousands of lines of raw OrdinaryDiffEq output to determine what went wrong.

# PART 1 — COMPONENT-LEVEL TESTS

Create small deterministic tests that bypass OrdinaryDiffEq wherever possible.

These should test individual numerical operations directly.

## Uniform / identical-state flux test

For any conservative numerical flux, if the left and right states are identical:

```julia
U_L == U_R
```

then the numerical flux should reduce to the physical flux to numerical precision.

Test this for representative physically valid states.

## Face conservation

For an internal face, verify that the contribution applied to one finite volume is exactly equal and opposite to the contribution applied to its neighbor, with correct geometric conventions.

This should include all conserved quantities relevant to the currently enabled model.

At minimum:

- mass,
- momentum,
- total energy,
- species when present.

## Reconstruction tests

For reconstruction schemes such as MUSCL:

### Constant field

A constant field must reconstruct exactly to the same constant value at every face.

### Linear field

Where the limiter should not activate, a linear field should reconstruct the analytically correct face values to numerical precision or expected discretization accuracy.

### Smooth-field convergence

For a smooth periodic function, measure reconstruction error on multiple grids and estimate spatial order.

For a second-order MUSCL formulation, the observed convergence order should approach approximately second order in the asymptotic regime.

Do not hard-code unrealistically strict finite-grid expectations; compute measured order and use a justified acceptance threshold.

### Discontinuity robustness

Test reconstruction around a discontinuity and verify:

- no NaN/Inf,
- no catastrophic overshoot,
- physical-state constraints where applicable.

Do not require mathematically impossible monotonic behavior from a limiter that does not promise it.

## Thornber low-Mach correction

Add direct tests for expected limiting behavior.

At minimum:

- identical states remain unchanged,
- at sufficiently high/supersonic Mach number the correction recovers ordinary reconstruction,
- at low Mach number the reconstructed left-right velocity jump is reduced according to the implemented Thornber formulation,
- the correction does not alter unrelated variables unless intended by the chosen formulation.

Derive expected values independently from the published formula or an independently coded reference implementation.

Do not use the production Thornber routine itself to generate expected values.

## Species transport

For species advection/diffusion:

- uniform species fraction gives zero diffusive flux,
- zero concentration gradient gives zero diffusive flux,
- a prescribed linear composition gradient produces the expected diffusion direction and magnitude for the implemented constitutive model,
- species fluxes obey the conservation constraints of the chosen diffusion formulation,
- transported species remain consistent with any sum-to-one constraints when the numerical method is intended to enforce them.

# PART 2 — SEMI-DISCRETE RESIDUAL TESTING

The most important object to verify is the finite-volume residual before involving a time integrator.

Create an interface that can evaluate conceptually:

```julia
R = residual(u, p, t)
```

or invoke the existing `rhs!` directly in a controlled test.

## Uniform-flow / free-stream preservation

Construct a uniform physically valid flow field.

For a domain and boundary conditions that should preserve a uniform state:

```text
density = constant
velocity = constant
pressure = constant
temperature = constant
species = constant where applicable
```

the residual should be approximately zero.

Report separate norms for:

- mass,
- momentum components,
- energy,
- species,
- turbulence variables if applicable.

## Global conservation

For a periodic or closed system with no relevant source terms, calculate the volume-integrated rate of change of conserved quantities.

Verify that internal-face contributions cancel globally.

For example:

```julia
sum(V .* dρdt) ≈ 0
sum(V .* dρEdt) ≈ 0
```

and analogous momentum/species checks where appropriate.

Use absolute and normalized residuals so that the diagnostic remains meaningful across state magnitudes.

## State validity

Every residual and integration test should automatically check for:

- NaN,
- Inf,
- nonphysical density,
- nonphysical pressure,
- nonphysical temperature,
- invalid mass fractions,
- invalid turbulence quantities when applicable.

The verifier should report:

- variable,
- cell index,
- time,
- offending value,
- nearby state if practical.

Do not allow the solver to simply crash later without reporting where state admissibility first failed.

# PART 3 — EMPIRICAL JACOBIAN DEPENDENCY / STENCIL DETECTION

This is a high-priority feature.

I use implicit OrdinaryDiffEq solvers and Jacobian sparsity information.

A numerical change such as MUSCL can expand the residual dependency stencil.

I want the verification suite to automatically detect this.

Implement a diagnostic that:

1. Evaluates a baseline residual.
2. Perturbs one cell/state degree of freedom at a time, or uses an efficient grouping strategy if needed.
3. Reevaluates the residual.
4. Determines which residual entries changed above a carefully chosen numerical threshold.
5. Constructs an empirical dependency graph / sparsity pattern.
6. Compares it with the sparsity pattern currently declared or inferred by FVMFramework.

Example failure report:

```text
Jacobian sparsity mismatch

Variable perturbed:
    density, cell 44

Observed residual dependencies:
    cells 42:46

Declared residual dependencies:
    cells 43:45

Unexpected dependencies:
    cells 42, 46
```

The implementation should be practical for small verification meshes, not necessarily production-size meshes.

It is acceptable for this diagnostic to be expensive because it will operate on intentionally tiny test problems.

The goal is correctness and interpretability.

# PART 4 — JACOBIAN-VECTOR PRODUCT VERIFICATION

Create a generic diagnostic comparing the derivative used by the implicit solver / AD system to an independent finite-difference directional derivative.

For a random normalized perturbation vector `v`, compare the framework Jacobian-vector product to something conceptually like:

```julia
Jv_fd =
    (R(u + ε*v) - R(u - ε*v)) / (2ε)
```

Use central finite differences where physically admissible.

Carefully choose or adapt epsilon based on state scaling and floating-point precision rather than blindly using one hard-coded number.

Report:

```text
||Jv_AD - Jv_FD||
relative error
maximum component error
index of maximum error
```

Where feasible, repeat for several random directions.

The test should distinguish:

- good agreement,
- probable finite-difference truncation/roundoff issues,
- definite derivative mismatch.

If the current solver uses a full Jacobian rather than explicit Jv operations, adapt this test accordingly.

# PART 5 — EXPLICIT VS IMPLICIT DIAGNOSTIC SOLVES

For every important integration-level benchmark, provide a small explicit reference/smoke integration in addition to the production implicit solve.

The explicit solver is NOT intended to be efficient for real CFD simulations.

Its purpose is diagnostic classification.

The verification report should clearly distinguish:

```text
Explicit integration: PASS
Implicit integration: FAIL
```

from:

```text
Explicit integration: FAIL
Implicit integration: FAIL
```

Interpretation should be included in the machine-readable diagnostics.

If explicit passes and implicit fails, likely investigation areas include:

- Jacobian,
- sparsity,
- nonlinear solver,
- preconditioner,
- Krylov configuration,
- nonsmooth residual behavior.

If both fail, likely investigation areas include:

- spatial discretization,
- state admissibility,
- boundary conditions,
- instability,
- implementation error.

Do not claim these diagnoses are mathematically certain; label them as likely failure categories.

Use conservative stable timesteps for the explicit diagnostic cases.

# PART 6 — ORDINARYDIFFEQ SOLVER-HEALTH METRICS

This is another high-priority feature.

A numerical change may produce the correct final field but dramatically worsen implicit solver behavior.

Capture as many reliable OrdinaryDiffEq statistics as reasonably available, for example:

- RHS/function evaluations,
- Jacobian evaluations,
- linear-solver iterations if exposed,
- nonlinear iterations if exposed,
- accepted steps,
- rejected steps,
- minimum accepted timestep,
- maximum timestep,
- solver return code,
- runtime,
- allocation information where useful and stable.

If certain metrics are unavailable for some solver configurations, handle this gracefully.

Create a baseline mechanism.

A test should be able to compare current metrics against a stored/reference baseline and report significant regressions.

For example:

```text
Correctness: PASS

Performance regression:
runtime              1.8 s -> 4.9 s
rejected steps       3 -> 87
linear iterations    420 -> 6100

Performance status: FAIL
```

Do not make runtime alone an extremely strict criterion because runtime is noisy.

Prefer deterministic or relatively stable quantities such as:

- iteration counts,
- rejected steps,
- function evaluations,

and use generous thresholds for wall-clock runtime.

Design baselines so they can be intentionally updated when a justified algorithmic change occurs.

# PART 7 — CANONICAL VERIFICATION PROBLEMS

Build the framework so benchmark cases are modular and easy to add.

Do NOT necessarily implement every benchmark below in the first commit if that would make the task too large.

Prioritize infrastructure plus the tests relevant to currently existing solver features.

Desired benchmark library eventually includes:

## Euler / inviscid compressible flow

- uniform/free-stream preservation,
- contact discontinuity,
- Sod shock tube,
- smooth periodic advection where appropriate,
- isentropic vortex,
- converging-diverging nozzle.

## MUSCL

- smooth convergence test,
- discontinuous transport,
- shock-tube comparison.

## Low-Mach correction

- smooth low-Mach vortex or advection case,
- pressure-fluctuation scaling if practical,
- recovery of normal compressible scheme at high Mach.

## Viscous flow

- Couette flow,
- Poiseuille flow,
- analytical or high-quality reference viscous test.

## Heat conduction

- one-dimensional conduction with analytical solution.

## Species

- scalar advection,
- one-dimensional diffusion,
- advection-diffusion analytical/reference case.

## SST k-omega, when implemented

Use an accepted turbulent benchmark such as:

- turbulent flat plate,
- turbulent channel,
- or another standard benchmark compatible with the solver geometry.

Keep benchmark datasets separate from production implementations.

# PART 8 — METHOD OF MANUFACTURED SOLUTIONS

Design an extensible Method of Manufactured Solutions framework, even if only a simple example is implemented initially.

The intended workflow is:

1. Define smooth analytical primitive fields such as:

```julia
rho_exact(x, y, t)
u_exact(x, y, t)
v_exact(x, y, t)
T_exact(x, y, t)
```

2. Independently calculate the continuous PDE derivatives and corresponding forcing/source terms.

3. Run the finite-volume solver with those forcing terms.

4. Compare numerical and exact solutions.

5. Repeat across mesh refinement levels.

6. Calculate measured convergence order.

Where practical, use automatic differentiation or symbolic-independent differentiation inside the verification code to evaluate derivatives of the manufactured functions.

Do not reuse the production finite-volume gradient/reconstruction routines to generate the manufactured forcing, since that would destroy independence.

Initially, focus on a manufactured case simple enough to debug.

The eventual goal is to support combinations such as:

- inviscid compressible flow,
- viscous compressible flow,
- species transport,
- possibly turbulence equations.

# PART 9 — CONVERGENCE STUDIES

Create reusable utilities that run a problem on progressively refined meshes and calculate convergence order.

For a scalar error metric:

```text
E_h
E_h/2
E_h/4
```

calculate observed order using the usual logarithmic ratio.

Support norms such as:

- L1,
- L2,
- Linf.

Reports should include a table similar to:

```text
Grid       L2 error       observed order
32         2.34e-3        -
64         6.21e-4        1.91
128        1.57e-4        1.98
```

Tests should check a justified minimum observed order, not demand exact theoretical order on coarse meshes.

# PART 10 — DIFFERENTIAL FEATURE TESTING

Numerical features should be easy to toggle so the verifier can isolate interactions.

Where architecture permits, support comparisons conceptually equivalent to:

```text
first order
first order + Thornber
MUSCL
MUSCL + Thornber
MUSCL + viscosity
MUSCL + species
```

The verification framework should make it easy to identify:

"Everything works until component X is enabled."

If enabling a feature requires intrusive changes to the existing architecture, do not redesign the entire package solely for this. Use the cleanest practical mechanism available.

# PART 11 — MUTATION TESTING OF THE VERIFIER

The verifier itself must be tested.

Implement at least a small mutation-testing mechanism or a documented set of intentional faults that can be enabled in test-only code.

Examples:

- reverse an internal face flux sign,
- omit a neighbor contribution,
- perturb an energy flux,
- use an intentionally incorrect velocity-gradient component,
- disable a dependency in the declared sparsity pattern,
- alter a reconstruction coefficient.

Verify that the expected tests fail.

Generate a report such as:

```text
Mutation: reverse_internal_face_flux

Expected detectors:
    face conservation
    global conservation
    shock benchmark

Detected:
    face conservation   PASS
    global conservation PASS
    shock benchmark     PASS

Mutation detection status: PASS
```

Here "PASS" means the verifier successfully caught the deliberately introduced bug.

Do not put dangerous mutation hooks into normal production execution paths.

# PART 12 — TEST ORGANIZATION

Create a clean directory structure.

Something conceptually like:

```text
test/
    verification/
        components/
        residual/
        jacobian/
        integration/
        benchmarks/
        manufactured/
        performance/
        mutations/
        utilities/
```

Adapt this to Julia package conventions and the existing project layout.

Each test should state:

- what physical/numerical property it verifies,
- why that property should hold,
- what independent oracle is used,
- what tolerance is justified.

Avoid tests that contain unexplained magic numbers.

# PART 13 — MACHINE-READABLE FAILURE CONTEXT FOR LLMS

Every failed test should expose enough structured information for another coding agent to act on it.

Examples:

For convergence:

```json
{
  "test": "muscl_smooth_convergence",
  "status": "fail",
  "expected_min_order": 1.8,
  "measured_order": 1.04,
  "errors": [0.01, 0.0049, 0.0024]
}
```

For sparsity:

```json
{
  "test": "jacobian_dependency",
  "status": "fail",
  "perturbed_cell": 44,
  "perturbed_variable": "rho",
  "declared_dependencies": [43, 44, 45],
  "observed_dependencies": [42, 43, 44, 45, 46]
}
```

For implicit integration:

```json
{
  "test": "implicit_smoke",
  "status": "fail",
  "explicit_status": "pass",
  "implicit_status": "fail",
  "return_code": "...",
  "likely_categories": [
    "jacobian",
    "sparsity",
    "nonlinear solver",
    "preconditioner"
  ]
}
```

Do not over-diagnose from insufficient evidence.

# PART 14 — REPRODUCIBILITY

All randomized tests must use deterministic seeds.

Record:

- Julia version,
- relevant package versions if practical,
- solver algorithm,
- tolerances,
- test mesh,
- random seed,
- relevant numerical-method configuration.

A failure should be reproducible from the report wherever practical.

# PART 15 — PERFORMANCE OF THE VERIFIER

The verification suite itself should remain practical.

Target rough categories:

```text
unit:
seconds or less where possible

integration:
seconds to perhaps tens of seconds

full:
allowed to be considerably longer
```

Use tiny meshes for Jacobian and structural diagnostics.

Do not run production-size CFD meshes for tests that only need 10–50 cells.

# PART 16 — DO NOT OVER-ENGINEER THE FIRST VERSION

Implement this in stages.

## Stage 1 — highest priority

Build these first:

1. verification/reporting framework,
2. constant/free-stream preservation,
3. face/global conservation,
4. NaN/Inf and physical-state checks,
5. empirical Jacobian dependency detection,
6. AD vs finite-difference Jv,
7. explicit vs implicit smoke solve,
8. OrdinaryDiffEq solver statistics,
9. at least one canonical benchmark for the currently working compressible solver,
10. machine-readable report generation.

Make Stage 1 solid before adding large numbers of benchmark problems.

## Stage 2

Add:

- MUSCL-specific reconstruction and convergence tests,
- Thornber low-Mach tests,
- species advection/diffusion tests,
- general convergence-study utilities,
- additional compressible benchmarks.

## Stage 3

Add:

- Method of Manufactured Solutions,
- viscous analytical benchmarks,
- turbulence benchmarks,
- mutation testing,
- performance baseline management.

If some Stage 2/3 infrastructure is trivial to add while implementing Stage 1, that is fine.

# PART 17 — HOW TO WORK ON THIS REPOSITORY

Before modifying anything:

1. Inspect the repository structure.
2. Locate:
   - state representation,
   - grid representation,
   - residual/RHS evaluation,
   - face flux routines,
   - Jacobian/sparsity infrastructure,
   - OrdinaryDiffEq problem construction,
   - existing tests.
3. Summarize the relevant architecture.
4. Identify the smallest set of interfaces needed by the verifier.
5. Avoid unnecessary production-code changes.

Then implement incrementally.

After each substantial change:

1. run existing tests,
2. run new unit verification,
3. run integration verification,
4. investigate failures rather than bypassing them.

Do not rewrite working solver components merely because you prefer a different architecture.

# PART 18 — RULES WHEN A TEST FAILS

You are explicitly prohibited from solving test failures by casually:

- increasing tolerances,
- removing tests,
- changing analytical reference values,
- weakening expected convergence order,
- deleting difficult benchmark cases,
- switching off physical admissibility checks.

If you believe a test is genuinely wrong:

1. explain why,
2. derive the correct expectation independently,
3. demonstrate the problem,
4. then modify the verifier.

Otherwise fix the implementation.

# PART 19 — COMMENTS AND DOCUMENTATION

I do not want an enormous amount of explanatory boilerplate.

Document:

- why a verification test exists,
- what property it verifies,
- where the expected result comes from,
- why a tolerance was chosen,
- unusual numerical considerations.

Avoid comments that merely restate obvious Julia syntax.

# PART 20 — DELIVERABLES

At the end, provide:

1. the implemented verification framework,
2. all added tests,
3. any minimal production-code changes required to expose testable interfaces,
4. example human-readable report,
5. example machine-readable report,
6. instructions for running:
   - unit verification,
   - integration verification,
   - full verification,
7. a short architecture summary,
8. a list of verification capabilities that remain unimplemented,
9. any solver weaknesses discovered while creating the verifier.

Most importantly, leave the project in a state where another coding agent can implement a CFD feature, run the verifier, receive useful structured failures, and iterate without requiring me to manually inspect every detail.

The end goal is not to prove the entire CFD solver correct mathematically.

The end goal is to make numerical mistakes, integration problems, Jacobian mistakes, sparsity mistakes, conservation errors, convergence failures, physical-state violations, and major solver-performance regressions **difficult to introduce silently**.
