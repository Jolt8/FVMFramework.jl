# Clapeyron `vt_flash` compilation investigation

## Summary

For the pure nitrous-oxide feed model, Clapeyron's generic `vt_flash` can
spend minutes compiling before the ODE solver begins. The thermodynamic solve
itself is not intrinsically that slow, and `ForwardDiff` produces sensible
derivatives once the code can execute.

The behavior was reproduced on:

- Julia 1.12.6
- Clapeyron 0.6.26
- `PR(["nitrous oxide"])`
- Windows
- 2026-10-22

The expensive path is inside Clapeyron's pure-component
`tx_flash_pure` machinery. Calling `saturation_pressure` or the direct VT
`pressure` method compiles normally. Passing a scalar amount, an `SVector`, or
calling `tx_flash_pure` directly did not remove the pathological compilation.
Disabling optimization with `-O0` also did not help.

Using `--compile=min` allows `vt_flash` to execute quickly enough to diagnose
the behavior. This strongly indicates an inference/native-code-generation
specialization problem rather than a flash convergence problem.

## Observed timings

The following are representative measurements from fresh Julia processes.
Times vary with hardware and cache state.

| Operation | First call | Warm call |
|---|---:|---:|
| Construct `PR` model | about 4.0 s | not measured |
| ForwardDiff through direct `pressure` | 0.77 s | 30 μs |
| `saturation_pressure` | 2.54 s | 36 μs |
| Normal `vt_flash` | over 180 s; did not finish | unavailable |
| `vt_flash` with `-O0` | over 180 s; did not finish | unavailable |
| Direct `tx_flash_pure` | over 120 s; did not finish | unavailable |
| `vt_flash` with `--compile=min` | 0.18 s | 0.12 s |
| ForwardDiff `vt_flash` with `--compile=min` | 0.39 s | 0.17 s |
| ForwardDiff flash-property Jacobian with `--compile=min` | 0.36 s | 0.19 s |
| Specialized pure flash, normal compilation | 2.58 s | 66 μs |
| ForwardDiff specialized flash, normal compilation | 5.79 s | 116 μs |

`--compile=min` is useful as a diagnostic. It is not a practical CFD
workaround: two flashes per RHS evaluation at roughly 0.1–0.2 seconds each
would make time integration extremely slow.

## ForwardDiff findings

With `--compile=min`, Clapeyron's `vt_flash` differentiated successfully in
single-phase liquid, two-phase equilibrium, and single-phase vapor states.
The tested output vector contained pressure, bulk density, bulk specific
internal energy, and bulk specific enthalpy.

For the initial liquid tank state, the pressure gradient with respect to
`[volume, temperature, mass]` was:

```text
[-1.0664468436178148e10, 365698.27013094584, 1.5046869045754004e7]
```

This matched the direct EOS `pressure` gradient.

For a constructed two-phase state, the pressure gradient was:

```text
[-0.0, 124364.16261538563, -0.0]
```

That behavior is physically consistent for a pure fluid inside the
saturation dome: pressure follows saturation temperature and is independent
of bulk volume and mass until a phase boundary is crossed.

Therefore, the problem is not a fundamental inability to propagate
`ForwardDiff.Dual` through the flash. It is compilation cost in normal Julia
compilation mode.

## SparseConnectivityTracer overloads

The probe deliberately includes `clapeyron_tracer_overloads.jl`, matching the
CFD case.

Those overloads affect mixed `SparseConnectivityTracer.Dual` and
`ForwardDiff.Dual` calculations. They do not affect ordinary ForwardDiff-only
calculations.

With `--compile=min`, local sparsity tracing through `vt_flash` completed and
returned a conservative dense dependency pattern for the tested properties.

Two overloads intentionally preserve dependency propagation rather than
correct numerical values:

```julia
Base.power_by_squaring(x::SparseConnectivityTracer.Dual,
                       p::SparseConnectivityTracer.Dual) = x + p

Base.ldexp(x::SparseConnectivityTracer.Dual, i::Int) = x
```

These must only be used for structural sparsity detection. Their primal values
and numerical derivatives are not physically meaningful. Because local
sparsity tracing uses primal values for control-flow decisions, they may also
produce conservative or imperfect branch-specific patterns.

## Implemented workaround

The injector train now uses `pure_component_vt_flash` from
`fvm/injector_delivery_train_scripts/injector_delivery_train_scripts.jl`.

For a pure, subcritical fluid it performs the VT equilibrium calculation
directly:

1. Call `saturation_pressure(model, temperature)` to obtain saturation
   pressure and saturated liquid/vapor molar volumes.
2. Compute the bulk molar volume `V/n`.
3. Below the saturated-liquid volume, use direct EOS pressure and construct a
   single liquid phase.
4. Above the saturated-vapor volume, use direct EOS pressure and construct a
   single vapor phase.
5. Inside the saturation dome, use saturation pressure and the lever rule to
   construct liquid and vapor phase amounts.
6. Return a normal Clapeyron `FlashResult`, so existing property calls continue
   to work.

The specialized result was compared with Clapeyron `vt_flash` for liquid,
two-phase, and vapor states. Pressure, density, specific internal energy, and
specific enthalpy agreed to a relative tolerance of `1e-9`. ForwardDiff also
returned the same initial liquid-state pressure gradient.

### Limitations

- The helper is for exactly one component.
- It assumes a subcritical temperature where saturation pressure and saturated
  phase volumes exist.
- It is not a replacement for a multicomponent flash.
- Phase-boundary transitions are physically nonsmooth; an AD calculation at
  the exact boundary should not be expected to have a unique derivative.
- A future supercritical feed state needs an explicit supercritical branch.

## Reproducing the issue

Run commands from the repository root.

Safe baseline and workaround probes:

```powershell
julia --project=. examples/1D_hybrid_rocket_engine/design_scripts/overloads/vt_flash_probe.jl baseline
julia --project=. examples/1D_hybrid_rocket_engine/design_scripts/overloads/vt_flash_probe.jl workaround
julia --project=. examples/1D_hybrid_rocket_engine/design_scripts/overloads/vt_flash_probe.jl injector
```

Compare the specialized flash with Clapeyron in diagnostic compilation mode:

```powershell
julia --project=. --compile=min examples/1D_hybrid_rocket_engine/design_scripts/overloads/vt_flash_probe.jl compare
```

Exercise ForwardDiff and SparseConnectivityTracer:

```powershell
julia --project=. --compile=min examples/1D_hybrid_rocket_engine/design_scripts/overloads/vt_flash_probe.jl forwarddiff
julia --project=. --compile=min examples/1D_hybrid_rocket_engine/design_scripts/overloads/vt_flash_probe.jl sparsity
```

Deliberately reproduce the normal-compilation problem with an external
timeout:

```powershell
julia --project=. examples/1D_hybrid_rocket_engine/design_scripts/overloads/vt_flash_probe.jl vt
```

## How to recognize the same problem later

The likely signature is:

- direct `pressure(model, V, T, n)` works;
- `saturation_pressure(model, T)` works;
- the first `vt_flash` consumes CPU and substantial memory for a very long
  time;
- the ODE solver reports no steps because compilation occurs before solving;
- `-O0` does not materially help;
- `--compile=min` reaches the flash result quickly;
- subsequent calls in a successfully compiled process are much faster.

If this behavior changes after upgrading Julia or Clapeyron, rerun the probe
before retaining the workaround. The generic upstream flash may eventually
compile normally and provide broader phase support.

## Separate energy-reference warning

The comparison probe returned approximately `-97.2 kJ/kg` for the initial
liquid tank specific internal energy. This is valid: absolute internal energy
depends on the thermodynamic model's reference state and does not have to be
positive.

Consequently, boundary-condition validation should check that specific
internal energy is finite, not that it is greater than zero. At the time of
this investigation, `subsonic_inlet_flux.jl` still contains a positivity
check. That check will reject a physically valid liquid nitrous state if the
mid-section becomes liquid.
