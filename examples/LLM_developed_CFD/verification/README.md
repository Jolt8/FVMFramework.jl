# Compressible-flow verification

The verifier exercises the existing conservative-state, FVMFramework group, HLLC,
and OrdinaryDiffEq paths on intentionally small hexahedral meshes. It provides:

- identical-state HLLC and pairwise face-conservation checks;
- free-stream and closed-domain global-conservation residual checks;
- localized NaN, Inf, density, pressure, and temperature diagnostics;
- finite-difference dependency detection against tracer-declared sparsity;
- ForwardDiff versus centered-finite-difference Jacobian-vector products;
- fixed-CFL SSPRK43 and implicit FBDF Sod shock-tube solves;
- exact Sod Riemann-solution errors and OrdinaryDiffEq health statistics;
- MUSCL constant, linear, discontinuity, free-stream, and smooth-order tests;
- general L1, L2, Linf, and observed-convergence-order utilities;
- empirical confirmation that MUSCL expands the cell stencil from radius 1 to 2;
- independent Thornber identical/high-Mach/low-Mach/tangential checks;
- a finite zero-velocity Thornber HLLC Jacobian and full-residual AD/finite-
  difference comparison;
- thermodynamic consistency after Thornber velocity reconstruction;
- a nonuniform, initially stagnant Thornber case through explicit and implicit
  integration;
- first-order, MUSCL, Thornber, and MUSCL-plus-Thornber feature combinations;
- conservative `rho*Y_k` species state, HLLC advection, and mixture-corrected
  Fickian diffusion checks;
- coupled explicit species-contact advection with an independent boundary-flux
  mass balance;
- direct weighted-least-squares and oblique corrected-face-gradient tests;
- analytical Couette, Poiseuille, and one-dimensional heat-conduction cases;
- viscous compressible and conservative species manufactured solutions;
- persistent solver-work baselines and long-time state/residual trajectories;
- gradient and transport-sign mutation checks;
- conservative SST `rho*k`/`rho*omega`, blending, source, flux, channel, and
  flat-plate checks;
- stationary-contact and first-order/MUSCL Sod benchmark comparisons;
- concise text and structured JSON reports.

Run from the repository root:

```powershell
julia --project=. examples/LLM_developed_CFD/verification/run_verification.jl unit
julia --project=. examples/LLM_developed_CFD/verification/run_verification.jl integration
julia --project=. examples/LLM_developed_CFD/verification/run_verification.jl full
```

An optional second argument selects the report directory. The Julia API is:

```julia
include("examples/LLM_developed_CFD/verification/FVMVerification.jl")
using .FVMVerification

report = validate(CompressibleFlowVerification(); level = :unit)
write_report(report, "path/to/report/directory")
```

`:unit` runs component, residual, admissibility, reconstruction, species,
low-Mach, viscous/thermal, MMS, mutation, sparsity, derivative, contact, and
convergence checks.
`:integration` adds 32-cell first-order and MUSCL Sod benchmarks with explicit
and implicit integrators, performance comparisons, and long-time cases.
`:full` uses 64-cell Sod cases and extends the smooth MUSCL convergence study
through 192 cells. Reports
record the mesh, tolerances, algorithms, package versions, and random seed.

## Stage 3 transport verification

The production internal-face correction now projects the complete interpolated
gradient onto the face normal and replaces only that normal component. Corrected
`v`, `w`, and temperature gradients are no longer overwritten by arithmetic
averages. Direct oblique-face tests reject both historical defects.

Couette and heat-conduction profiles are reproduced to roundoff. The
cell-centred Poiseuille benchmark converges at second order; its 64-cell Linf
velocity error is approximately `2.29e-4`. These compact one-dimensional cases
isolate the production gradient and transport signs from Riemann-solver
dissipation.

The viscous compressible MMS uses constant density/pressure with sinusoidal
tangential velocity. It verifies momentum diffusion and viscous work in total
energy, with asymptotic order about `1.98`. The conservative two-species MMS
obtains order `2.01` and preserves the species sum to roundoff. MMS norms omit
the two boundary cells; boundary flux behavior is covered by the analytical
transport cases.

Long-time JSON results include sampled minimum-density, minimum-pressure, and
scaled-residual histories. A nonuniform stationary contact remains unchanged
through `t=1`, with explicit/implicit Linf disagreement below `5e-15`. A smooth
acoustic/composition startup evolves by about `3.34e-2`, stays admissible, and
has explicit/implicit relative Linf disagreement near `1.30e-6` at `t=0.5`.

## Performance baselines

Versioned references live in
[`performance/performance_baselines.toml`](performance/performance_baselines.toml).
Every 32- and 64-cell first-order, MUSCL, and Thornber Sod path records accepted
and rejected steps, RHS calls, linear solves, nonlinear iterations, minimum
timestep, and runtime. Work counts and minimum timestep are gating; runtime is a
loose, non-gating warning because compilation and host load are variable.

At 32 cells, MUSCL FBDF currently uses 517 accepted / 94 rejected steps and
1,087 RHS calls, versus fixed-CFL SSPRK43's 23 / 0 and 92. The stored 64-cell
references are 854 / 163 / 1,972 versus 46 / 0 / 184. To refresh a reference,
run the corresponding tier, inspect the solution and old/new JSON statistics,
then deliberately update the TOML values and date.

## SST k-omega

The production SST implementation is in
[`turbulence/sst_k_omega.jl`](../navier_stokes_testing_grounds/turbulence/sst_k_omega.jl).
It provides the standard Menter constants, `F1`/`F2` blending, strain invariant,
eddy-viscosity and production limiters, destruction and cross-diffusion terms,
conservative HLLC advection, blended diffusion, and the wall `omega` treatment.
The existing viscous flux consumes effective molecular-plus-turbulent viscosity
and conductivity and includes the isotropic `-2*rho*k/3` Reynolds stress.

The species-capable Navier--Stokes example now carries conservative `rho*k` and
`rho*omega` through its region, inlet, outlet, internal-face, and wall paths. A
non-integrating setup smoke check is available:

```powershell
julia --project=. examples/LLM_developed_CFD/navier_stokes_testing_grounds/new_navier_stokes_solver.jl --setup-only
```

Verification covers blending limits, analytical homogeneous decay, conservative
face fluxes, the complete FVM operator, positive-state rejection, an equilibrium
turbulent-channel log layer, and a separate smooth-flat-plate skin-friction
dataset. The channel and flat-plate checks validate closure asymptotes and wall
stress; they are not claims of mesh-converged RANS solutions. A resolved
flat-plate or channel computation still needs an adequately refined wall-normal
mesh and documented inlet/wall-resolution study—the current example has only
one cell in each transverse direction.

## No-species solver diagnostics

The 100-cell no-species driver is intentionally configurable separately from
the small verification cases:

```powershell
# Verified default: ordinary HLLC, slip walls, linearly implicit transient
julia --project=. examples/LLM_developed_CFD/navier_stokes_testing_grounds/new_navier_stokes_solver_no_species.jl --tmax=0.1

# Also run the scaled nonlinear least-squares steady solve
julia --project=. examples/LLM_developed_CFD/navier_stokes_testing_grounds/new_navier_stokes_solver_no_species.jl --tmax=0.01 --steady

# Diagnostic modes
julia --project=. examples/LLM_developed_CFD/navier_stokes_testing_grounds/new_navier_stokes_solver_no_species.jl --tmax=0.01 --thornber
julia --project=. examples/LLM_developed_CFD/navier_stokes_testing_grounds/new_navier_stokes_solver_no_species.jl --no-slip-walls
```

`--tmax=`, `--abstol=`, and `--reltol=` set the transient controls, while
`--progress` enables throttled progress output. Ordinary HLLC is the default
because the prescribed inlet is approximately Mach 1.7. Thornber is retained
behind `--thornber` for low-Mach and diagnostic work rather than being applied
silently during the mixed-Mach startup.

The steady-state version works for the verified default configuration. It
scales the conservative variables, uses nonlinear least squares with
Levenberg-Marquardt and sparse finite differences, and starts from the known
uniform inlet state for the slip-wall problem. Ordinary and Thornber-uniform
runs both returned `Success`, an admissible state, and a scaled Linf residual
of approximately `2.06e-11`.

This is not evidence of global Newton convergence. The original stagnant guess
is far from the root, the old zero-speed Thornber norm generated non-finite AD
Jacobians, and HLLC wave-selection branches remain difficult for direct Newton.
The optional no-slip case uses a transient warm start and has not yet been
shown to converge to a steady root.

The shock-limited MUSCL case uses `FBDF(AutoFiniteDiff)` because the limiter
changes branches at discontinuities. Smooth-state MUSCL still undergoes the
normal ForwardDiff-versus-centered-finite-difference Jv verification. This
backend choice is recorded in the solver statistics.

The finite-difference sparsity checker treats extra declared entries as
inactive-at-this-state rather than errors. It fails only when a numerically
observed dependency is absent from the declared pattern.

The compressible state now evolves species densities `rho*Y_k`; mass fractions
are derived caches. Both explicit and implicit Sod integrations include these
extra degrees of freedom, and a separate transported-contact case verifies the
coupled species path. The present species model is passive: mixture `cp` and
`cv` remain prescribed fluid properties, species enthalpy diffusion and
chemistry are not yet coupled, and species face values remain first order when
the Euler variables use MUSCL.

See [ARCHITECTURE.md](ARCHITECTURE.md) for the adapter design, verification
oracles, current limitations, and solver weaknesses found during Stages 1–3.
