# Stage 1–2 architecture and findings

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

The canonical cases are a stationary contact discontinuity and the ideal-gas
Sod shock tube. SSPRK43 uses a fixed acoustic CFL of 0.2. First-order FBDF uses
ForwardDiff; shock-limited MUSCL uses `AutoFiniteDiff` because limiter branch
changes at the discontinuity prevent useful Newton derivatives from ForwardDiff.
The benchmark values come from a separately implemented exact Riemann solution.

MUSCL reconstructs density, pressure, and velocity. Reconstructing pressure
rather than density and temperature independently preserves constant pressure
across a material contact. A smooth algebraic intersection combines the
Venkatakrishnan and strict local-bound limiters, retaining the tighter limit
without a hard `min` between them. The general convergence utility supplies
weighted L1/L2/Linf norms and nonuniform-refinement order calculations.

The compressible state contains one conservative density `rho*Y_k` per species.
Region property updates recover `Y_k` into nested caches, HLLC uses its density
flux to advect the upwind composition, and the region cap converts integrated
species flow to `d(rho*Y_k)/dt`. Transmissive and prescribed-inlet boundaries
use the same density-flux coupling. The verifier checks the pointwise identities
`sum(rho*Y_k) = rho` and `sum(d(rho*Y_k)/dt) = d(rho)/dt` through the complete
semi-discrete operator.

Diffusion uses the production Fickian flux followed by the mixture correction
`J_k <- J_k - Y_k*sum(J)`. This makes the total diffusive mass flux zero even
for unequal species diffusivities. The original oriented mass-fraction helper
is retained and tested for compatibility, while the coupled solver writes
equal-and-opposite fluxes directly into conservative species-flow caches.

## Thornber low-Mach reconstruction

The original Thornber implementation evaluated
`sqrt(u^2 + v^2 + w^2)` at exactly stagnant states. Its value is finite at the
origin, but its derivative is undefined, and ForwardDiff produced non-finite
HLLC Jacobian entries. This was the cause of the immediate Newton failure when
the no-species solver started from a stationary domain.

The production correction now uses a smooth velocity norm with a Mach
regularization of `1e-3`. The local maximum and the cap at Mach one use compact
C1 transitions. These changes define finite derivatives at zero velocity,
equal neighboring Mach numbers, and the Mach-one transition. The verifier now
checks the complete Thornber HLLC ForwardDiff Jacobian at a nonuniform,
zero-velocity state and compares the full semi-discrete Jacobian-vector product
with centered finite differences.

After Thornber changes reconstructed velocity, HLLC rebuilds total energy by
adding the corrected-minus-original kinetic-energy density. This preserves the
original internal energy and pressure. Previously, corrected momentum was used
with the uncorrected total energy while HLLC continued using the original
pressure, so the reconstructed conservative and primitive states were not
thermodynamically consistent. A dedicated regression now verifies recovered
pressure after the correction.

The integration suite includes a nonuniform, initially stagnant Thornber Sod
case through both SSPRK43 and FBDF. This is intentionally separate from the
large supersonic-inlet example: passing a low-Mach component or Sod test does
not establish that applying the correction during a mixed-Mach startup is a
good modelling choice.

## No-species driver and steady-state scope

`new_navier_stokes_solver_no_species.jl` now makes the important numerical
choices explicit:

- ordinary HLLC is the default for its Mach-1.7 prescribed inlet;
- `--thornber` enables the low-Mach correction as a diagnostic mode;
- slip, adiabatic side walls are the default verified configuration;
- `--no-slip-walls` restores the viscous wall treatment, which needs adequate
  transverse resolution and a boundary-condition study;
- `Rosenbrock23(AutoFiniteDiff, Sparspak)` replaces the unspecified automatic
  ODE algorithm and avoids nested nonlinear failures during the tested startup;
- density and internal-energy admissibility are checked before accepting ODE
  states, saving is bounded, progress output is optional, and final time and
  tolerances are command-line configurable.

The steady path scales every conservative field, formulates the problem as
nonlinear least squares, and uses Levenberg-Marquardt with sparse finite
differences. For the default slip-wall case, the uniform prescribed-inlet state
is an exact boundary-informed steady guess. Both ordinary HLLC and the
high-Mach-uniform Thornber mode return `Success`, remain admissible, and have a
scaled Linf residual of approximately `2.06e-11` in the tested configuration.

That result has a deliberately narrow interpretation. It proves that the
configured steady residual recognizes and accepts the known uniform root. It
does not prove global convergence from the original stagnant state. Direct
Newton from the stagnant or merely transient-warmed state still encounters
HLLC branch sensitivity; trust-region, pseudo-transient, and least-squares
experiments reduced the residual but could stall. The optional no-slip case
uses the transient result as its steady initial guess and has not been shown to
converge to a steady root.

## Current limitations

- Species currently behave as passive scalars. Composition does not yet update
  mixture thermodynamics or transport properties, species enthalpy diffusion
  is absent, and chemistry source terms are not connected to `rho*Y_k`.
- Species advection is first-order upwind even when density, pressure, and
  velocity use MUSCL reconstruction. A bounded second-order species
  reconstruction needs dedicated gradient caches and convergence tests.
- Manufactured solutions, viscous analytical cases, turbulence, mutations,
  and persistent performance baselines remain Stage 3 work.
- State validation covers conservative ideal-gas density, pressure,
  temperature, non-negative species densities, and the species-density sum.
  Turbulence-specific admissibility rules will be needed when those state
  fields are introduced.
- The sparsity diagnostic intentionally uses tiny meshes and one smooth state;
  it reports inactive declared entries rather than claiming they are wrong.
- Solver statistics are captured but not yet compared with a versioned baseline.
- Long-time convergence of the 100-cell supersonic-inlet driver with Thornber
  is not verified. Its mixed-Mach startup develops a materially different flow
  from ordinary HLLC, including backflow by 1 s in the tested configuration.
- The steady solve is verified only for the known uniform slip-wall root. The
  no-slip configuration and convergence from a generic distant initial guess
  remain open problems.

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
- On the MUSCL Sod runs, FBDF with ForwardDiff could not accept its first step
  because the shock limiter changes branches. `AutoFiniteDiff` completes the
  solve, but on the latest 64-cell case it required 906 accepted and 167 rejected
  steps versus SSPRK43's 46 accepted fixed-CFL steps. This is a strong candidate
  for a future performance baseline and nonlinear-solver investigation.
- The existing example is monolithic and executes a large solve when included.
  The verifier therefore constructs a small equivalent configuration rather
  than importing the top-level example script.
- The former Thornber zero-speed norm generated non-finite ForwardDiff
  Jacobians and caused immediate Newton failure. Smooth Mach regularization
  fixes that derivative defect, but it does not make Thornber appropriate for
  the supersonic-inlet startup.
- Leaving the ODE algorithm unspecified obscured whether the driver was using
  an explicit or stiff path and led to misleading `MaxIters` failures. The
  no-species driver now selects a linearly implicit method explicitly and
  reports accepted/rejected steps and nonlinear failures.
