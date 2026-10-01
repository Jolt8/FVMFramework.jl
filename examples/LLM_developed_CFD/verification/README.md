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
low-Mach, sparsity, derivative, contact, and convergence checks.
`:integration` adds 32-cell first-order and MUSCL Sod benchmarks with explicit
and implicit integrators. `:full` uses 64-cell Sod cases and extends the smooth
MUSCL convergence study through 192 cells. Reports
record the mesh, tolerances, algorithms, package versions, and random seed.

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

Stage 3 items intentionally not included here are manufactured solutions,
viscous analytical cases, turbulence benchmarks, mutation testing, and
persistent performance-baseline management.

See [ARCHITECTURE.md](ARCHITECTURE.md) for the adapter design, verification
oracles, current limitations, and solver weaknesses found during Stages 1–2.
