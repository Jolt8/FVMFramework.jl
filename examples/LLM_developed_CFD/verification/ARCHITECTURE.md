# Stage 1 architecture and findings

## Architecture

`finish_fvm_config` flattens the conservative `ComponentVector` state and
retains axes that map it back to named cell fields. It also builds region,
patch, and internal-connection groups plus independent cache/property storage.
`fvm_operator!` creates `VirtualFVMArray` views, runs property updates, evaluates
each internal face once, evaluates boundary patches, and caps integrated face
flows by cell volume.

The verifier follows that production path on small generated hexahedral grids.
`CompressibleCase` contains the finished `FVMSystem`, geometry, conservative
initial state, and residual closure. `CompressibleFlowVerification` selects
the unit and benchmark mesh sizes. The reporting layer is independent of the
compressible adapter and stores every result as measured metrics, expectations,
and actionable diagnostics.

The empirical sparsity oracle uses centered finite differences. It compares
observed dependencies with `SparseConnectivityTracer` output and fails on
observed-but-undeclared entries. Declared-but-inactive entries are reported but
do not fail because a single physical state cannot activate every HLLC branch.
The derivative oracle compares ForwardDiff with independently evaluated,
scaled centered differences in three deterministic random directions.

The canonical case is the ideal-gas Sod shock tube. SSPRK43 uses a fixed
acoustic CFL of 0.2; FBDF uses tracer sparsity and its normal AD path. The
benchmark values come from a separately implemented exact Riemann solution.

## Current limitations

- Stage 1 covers the inviscid first-order HLLC path. MUSCL, Thornber, species,
  manufactured solutions, viscous analytical cases, turbulence, mutations,
  convergence studies, and persistent performance baselines remain later-stage
  work.
- State validation covers conservative ideal-gas density, pressure, and
  temperature plus optional mass fractions. Turbulence-specific admissibility
  rules will be needed when those state fields are introduced.
- The sparsity diagnostic intentionally uses tiny meshes and one smooth state;
  it reports inactive declared entries rather than claiming they are wrong.
- Solver statistics are captured but not yet compared with a versioned baseline.

## Solver weaknesses found

- `primitive_from_conservative` and `get_temperature_ideal` clamp density and
  internal energy to `1e-10`. That keeps downstream arithmetic finite but can
  mask where an inadmissible conservative state first appeared. The verifier
  therefore checks the unclamped conservative state before interpreting it.
- In `fluid_viscous_and_diffusive_flux!`, corrected face gradients are computed
  for all velocity components and temperature, but the `v`, `w`, and
  temperature results are immediately overwritten by arithmetic averages.
  The scalar-component form used inside `corrected_face_gradient` also deserves
  a dedicated analytical test. This routine is outside the inviscid Stage 1
  benchmark and was not modified here.
- On the 32-cell Sod run used during development, FBDF needed substantially
  more steps and allocations than fixed-CFL SSPRK43 and began with a very small
  accepted timestep. Both methods were correct, but these metrics are useful
  candidates for the future baseline mechanism.
- The existing example is monolithic and executes a large solve when included.
  The verifier therefore constructs a small equivalent configuration rather
  than importing the top-level example script.
