using ADTypes
using Clapeyron
using ComponentArrays
using ForwardDiff
using SparseConnectivityTracer
using StaticArrays
using Unitful

# Keep this include in the probe: the CFD case loads these mixed
# SparseConnectivityTracer/ForwardDiff overloads before constructing sparsity.
include(joinpath(@__DIR__, "clapeyron_tracer_overloads.jl"))

# Load the exact specialized flash used by the injector train rather than a
# duplicated approximation that could drift away from production behavior.
include(joinpath(
    @__DIR__,
    "..",
    "..",
    "fvm",
    "injector_delivery_train_scripts",
    "injector_delivery_train_scripts.jl",
))

function timed_call(f, label)
    GC.gc()
    result = @timed f()
    println(
        label,
        ": time=", round(result.time; digits = 6),
        " s, compile=", round(result.compile_time; digits = 6),
        " s, recompile=", round(result.recompile_time; digits = 6),
        " s, bytes=", result.bytes,
    )
    flush(stdout)
    return result.value
end

function print_usage()
    println("""
    Clapeyron pure-fluid VT flash probe

    Run from the repository root:

      julia --project=. examples/1D_hybrid_rocket_engine/design_scripts/overloads/vt_flash_probe.jl baseline
      julia --project=. examples/1D_hybrid_rocket_engine/design_scripts/overloads/vt_flash_probe.jl workaround
      julia --project=. examples/1D_hybrid_rocket_engine/design_scripts/overloads/vt_flash_probe.jl injector
      julia --project=. --compile=min examples/1D_hybrid_rocket_engine/design_scripts/overloads/vt_flash_probe.jl compare
      julia --project=. --compile=min examples/1D_hybrid_rocket_engine/design_scripts/overloads/vt_flash_probe.jl forwarddiff
      julia --project=. --compile=min examples/1D_hybrid_rocket_engine/design_scripts/overloads/vt_flash_probe.jl sparsity

    The `vt` mode deliberately calls Clapeyron.vt_flash in normal compilation
    mode and may take minutes or appear hung. Run it with an external timeout.
    `--compile=min` is useful for diagnosis but is too slow for the CFD solve.
    """)
end

mode = isempty(ARGS) ? "help" : only(ARGS)
mode == "help" && (print_usage(); exit())

println(
    "mode=", mode,
    ", Julia=", VERSION,
    ", Clapeyron=", Base.pkgversion(Clapeyron),
)
flush(stdout)

model = timed_call("construct PR model") do
    PR(["nitrous oxide"])
end

tank_volume = 1.6e-3
tank_temperature = 298.0
tank_mass = 1.134
molecular_weight = 0.044013
moles(mass) = mass / molecular_weight

x0 = [tank_volume, tank_temperature, tank_mass]

direct_pressure(x) =
    pressure(model, x[1], x[2], SVector(moles(x[3])))
standard_flash(x) =
    vt_flash(model, x[1], x[2], moles(x[3]))
specialized_flash(x) =
    pure_component_vt_flash(model, x[1], x[2], moles(x[3]))

standard_pressure(x) = pressure(standard_flash(x))
specialized_pressure(x) = pressure(specialized_flash(x))

function flash_properties(flash_function, x)
    result = flash_function(x)
    return [
        pressure(result),
        mass_density(model, result),
        mass_internal_energy(model, result),
        mass_enthalpy(model, result),
    ]
end

function phase_test_states()
    _, liquid_molar_volume, vapor_molar_volume =
        saturation_pressure(model, tank_temperature)
    two_phase_moles =
        tank_volume /
        (0.5 * liquid_molar_volume + 0.5 * vapor_molar_volume)
    vapor_moles = 0.5 * tank_volume / vapor_molar_volume

    return (
        liquid = x0,
        two_phase = [
            tank_volume,
            tank_temperature,
            two_phase_moles * molecular_weight,
        ],
        vapor = [
            tank_volume,
            tank_temperature,
            vapor_moles * molecular_weight,
        ],
    )
end

if mode == "baseline"
    p1 = timed_call("direct pressure first") do
        direct_pressure(x0)
    end
    p2 = timed_call("direct pressure second") do
        direct_pressure(x0)
    end
    saturation1 = timed_call("saturation_pressure first") do
        saturation_pressure(model, tank_temperature)
    end
    saturation2 = timed_call("saturation_pressure second") do
        saturation_pressure(model, tank_temperature)
    end
    @show p1 p2 saturation1 saturation2
elseif mode == "vt"
    println("WARNING: normal compilation of this call is the reproduced issue.")
    flush(stdout)
    flash1 = timed_call("Clapeyron vt_flash first") do
        standard_flash(x0)
    end
    flash2 = timed_call("Clapeyron vt_flash second") do
        standard_flash(x0)
    end
    @show pressure(flash1) pressure(flash2)
elseif mode == "workaround"
    flash1 = timed_call("specialized pure flash first") do
        specialized_flash(x0)
    end
    flash2 = timed_call("specialized pure flash second") do
        specialized_flash(x0)
    end
    gradient1 = timed_call("ForwardDiff specialized pressure first") do
        ForwardDiff.gradient(specialized_pressure, x0)
    end
    gradient2 = timed_call("ForwardDiff specialized pressure second") do
        ForwardDiff.gradient(specialized_pressure, x0)
    end
    @show pressure(flash1) pressure(flash2) gradient1 gradient2
elseif mode == "injector"
    injector_u0, properties = get_injector_train_u0_and_properties(model)
    u = ComponentVector(
        tank = (
            mass = injector_u0.tank.mass,
            temperature = injector_u0.tank.temperature,
        ),
        mid_section = (
            mass = injector_u0.mid_section.mass,
            temperature = injector_u0.mid_section.temperature,
        ),
    )
    du = similar(u)
    fill!(du, 0.0)
    system = (
        injector_train_properties = properties,
        nitrous_oxide_model = model,
    )
    timed_call("update injector thermodynamic properties") do
        update_injector_train_properties!(
            du,
            u,
            nothing,
            0.0,
            system,
            nothing,
        )
    end
    @show properties.tank properties.mid_section
elseif mode == "compare"
    for (phase_name, state) in pairs(phase_test_states())
        standard = timed_call("$phase_name Clapeyron vt_flash") do
            standard_flash(state)
        end
        specialized = timed_call("$phase_name specialized flash") do
            specialized_flash(state)
        end
        standard_values = flash_properties(standard_flash, state)
        specialized_values = flash_properties(specialized_flash, state)
        @assert isapprox(standard_values, specialized_values; rtol = 1e-9)
        println(
            phase_name,
            ": standard phases=", length(standard.fractions),
            ", specialized phases=", length(specialized.fractions),
            ", pressure=", pressure(specialized),
        )
        @show standard_values specialized_values
        flush(stdout)
    end
elseif mode == "forwarddiff"
    for (phase_name, state) in pairs(phase_test_states())
        pressure_gradient = timed_call(
            "$phase_name ForwardDiff vt_flash pressure",
        ) do
            ForwardDiff.gradient(standard_pressure, state)
        end
        property_jacobian = timed_call(
            "$phase_name ForwardDiff vt_flash properties",
        ) do
            ForwardDiff.jacobian(
                x -> flash_properties(standard_flash, x),
                state,
            )
        end
        @show phase_name pressure_gradient property_jacobian
        flush(stdout)
    end
elseif mode == "sparsity"
    detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
    pattern = timed_call("SparseConnectivityTracer vt_flash properties") do
        ADTypes.jacobian_sparsity(
            x -> flash_properties(standard_flash, x),
            x0,
            detector,
        )
    end
    @show Matrix(pattern)
else
    error("unknown mode: $mode")
end
