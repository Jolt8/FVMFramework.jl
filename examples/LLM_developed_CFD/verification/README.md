# Compressible-flow verification

Stage 1 exercises the existing conservative-state, FVMFramework group, HLLC,
and OrdinaryDiffEq paths on intentionally small hexahedral meshes. It provides:

- identical-state HLLC and pairwise face-conservation checks;
- free-stream and closed-domain global-conservation residual checks;
- localized NaN, Inf, density, pressure, and temperature diagnostics;
- finite-difference dependency detection against tracer-declared sparsity;
- ForwardDiff versus centered-finite-difference Jacobian-vector products;
- fixed-CFL SSPRK43 and implicit FBDF Sod shock-tube solves;
- exact Sod Riemann-solution errors and OrdinaryDiffEq health statistics;
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

`:unit` runs component, residual, admissibility, sparsity, and derivative
checks. `:integration` adds a 32-cell Sod benchmark with both integrators.
`:full` runs the same Stage 1 properties with a 64-cell benchmark. Reports
record the mesh, tolerances, algorithms, package versions, and random seed.

The finite-difference sparsity checker treats extra declared entries as
inactive-at-this-state rather than errors. It fails only when a numerically
observed dependency is absent from the declared pattern.

Stage 2/3 items intentionally not included here are MUSCL order studies,
Thornber limiting tests, species transport, manufactured solutions, viscous
analytical cases, turbulence benchmarks, mutation testing, and persistent
performance-baseline management.

See [ARCHITECTURE.md](ARCHITECTURE.md) for the adapter design, verification
oracles, current limitations, and solver weaknesses found during Stage 1.
