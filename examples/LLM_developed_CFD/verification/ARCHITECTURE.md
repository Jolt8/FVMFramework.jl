# Stage 1–3 architecture and findings

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

## Viscous and thermal verification

The internal-face `corrected_face_gradient` previously made three independent
scalar products instead of projecting the complete gradient. It now evaluates
`dot(grad_avg, normal)` once and adds the resulting correction along the normal.
The flux routine also no longer discards its corrected `v`, `w`, and temperature
gradients. The operation preserves the interpolated tangential gradient and
replaces only its normal component.

Three direct checks protect this code: weighted least squares must recover a
resolved linear field, an oblique-face oracle checks the full projection, and
distinct `u`, `v`, `w`, and temperature inputs verify every transported field.
The Couette, Poiseuille, and conduction cases use the same production correction
inside compact one-dimensional finite-volume operators. This exposes transport
signs and orders without conflating them with HLLC dissipation. Poiseuille gives
second-order convergence; Couette and conduction reproduce their linear steady
solutions to roundoff.

## Manufactured solutions

The first MMS is a special case of compressible Navier--Stokes: constant
density, pressure, and temperature with tangential
`u(y)=u0+A*sin(2*pi*y)`. Independently evaluated continuous sources cancel
Newtonian momentum diffusion and the divergence of viscous work in total
energy. Density remains exactly steady and the measured interior asymptotic
order is about 1.98.

The second MMS uses `Y_a=0.4+0.1*sin(2*pi*y)` and `Y_b=1-Y_a` with constant
density and equal diffusivity. It calls the production conservative,
mixture-corrected diffusion flux, obtains order 2.01, and keeps the summed
species residual at roundoff. Interior norms are used because the MMS supplies
analytical boundary fluxes; boundary behavior is independently covered by the
analytical transport cases.

## Persistent performance and long-time evidence

`performance/performance_baselines.toml` stores solver-work expectations for
32- and 64-cell first-order, MUSCL, and Thornber Sod runs. Accepted/rejected
steps, RHS evaluations, linear solves, nonlinear iterations, and minimum
timestep are gating. Runtime is only a loose, non-gating warning. Missing
baseline entries fail their checks.

The main cost is now quantified: 32-cell MUSCL FBDF requires 517 accepted / 94
rejected steps and 1,087 RHS calls versus SSPRK43's 23 / 0 and 92. The stored
64-cell reference is 854 / 163 and 1,972 versus 46 / 0 and 184.

Long-time integration records sampled minimum-density, minimum-pressure, and
scaled-residual trajectories. A nonuniform stationary contact is the known
steady oracle and remains unchanged through `t=1` with both integrators. A
smooth acoustic/composition perturbation provides a nontrivial startup; both
paths stay admissible and agree to about `1.30e-6` relative Linf at `t=0.5`.

## Mutation and SST k-omega

Mutation checks show that the suite rejects the historical componentwise
gradient projection and reversed viscous/species-diffusion signs. At 64 cells,
the sign mutations amplify their MMS errors by roughly 2,490 times.

The production SST implementation evolves conservative `rho*k` and `rho*omega`.
HLLC mass flux upwinds both quantities; internal faces add the molecular-plus-
SST diffusion coefficients; cell sources contain limited production, destruction,
and cross diffusion. `F1` blends the inner/outer constants, `F2` enters the
eddy-viscosity limiter, and the wall treatment imposes zero `k` with the standard
near-wall `omega` scale. State validation requires both conservative fields to
remain finite and strictly positive.

The momentum/energy viscous flux uses `mu + mu_t`, turbulent conductivity with
`Pr_t=0.9`, and the isotropic `-2*rho*k/3` Reynolds stress. The species-capable
Navier--Stokes example connects the same implementation through initialization,
inlet/outlet advection, internal diffusion, wall fluxes, source capping, and its
admissibility check.

Tests cover inner/outer blending limits, homogeneous decay source values,
equal-and-opposite face transport, a complete FVM residual, and time integration
against the analytical decay solution. Canonical closure checks verify the
channel log-layer identities `mu_t=rho*kappa*y*u_tau` and `P_k=beta*rho*k*omega`,
plus skin friction against a separate smooth-flat-plate reference table.
These are equilibrium/asymptotic closure validations, not mesh-converged RANS
solutions. The current example's one-cell transverse mesh is sufficient for a
plumbing smoke check but not a credible wall-bounded-flow calculation.

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
- A resolved SST flat-plate/channel solve still needs a wall-normal mesh study,
  inlet turbulence sensitivity, and comparison with profile/reference data.
- State validation covers conservative ideal-gas density, pressure,
  temperature, non-negative species densities, and the species-density sum.
  Conventional SST `k` and `omega` field names are also checked for strict
  positivity when present.
- The sparsity diagnostic intentionally uses tiny meshes and one smooth state;
  it reports inactive declared entries rather than claiming they are wrong.
- Performance references are specific to the 32- and 64-cell Sod cases; a new
  mesh or algorithm must add a reviewed baseline entry.
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
- The internal viscous face projection and overwritten `v`, `w`, and temperature
  gradients were real defects. They are fixed and protected by direct oblique-
  gradient, analytical transport, MMS, and mutation checks.
- On the MUSCL Sod runs, FBDF with ForwardDiff could not accept its first step
  because the shock limiter changes branches. `AutoFiniteDiff` completes the
  solve. The current reviewed 64-cell baseline is 854 accepted and 163 rejected
  steps versus SSPRK43's 46 accepted fixed-CFL steps; this remains a nonlinear-
  solver investigation target.
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
