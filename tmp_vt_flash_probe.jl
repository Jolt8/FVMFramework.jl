using Clapeyron
using ForwardDiff

function timed_call(f, label)
    GC.gc()
    result = @timed f()
    println(
        label,
        ": time=",
        round(result.time; digits = 6),
        " s, compile=",
        round(result.compile_time; digits = 6),
        " s, recompile=",
        round(result.recompile_time; digits = 6),
        " s, bytes=",
        result.bytes,
    )
    flush(stdout)
    return result.value
end

println("Julia ", VERSION, ", Clapeyron probe")
flush(stdout)

model = timed_call("construct PR model") do
    PR(["nitrous oxide"])
end

tank_volume = 1.6e-3
tank_temperature = 298.0
tank_moles = 1.134 / 0.044013

direct_pressure() =
    pressure(model, tank_volume, tank_temperature, [tank_moles])
flash_at(T, n) = vt_flash(model, tank_volume, T, [n])
flash_pressure(T, n) = pressure(flash_at(T, n))
flash_internal_energy(T, n) = mass_internal_energy(model, flash_at(T, n))

timed_call("direct pressure first") do
    direct_pressure()
end
timed_call("direct pressure second") do
    direct_pressure()
end

flash_1 = timed_call("vt_flash Float64 first") do
    flash_at(tank_temperature, tank_moles)
end
println(
    "  phases=",
    length(flash_1.fractions),
    ", vapour_idx=",
    flash_1.data.vapour_idx,
    ", pressure=",
    pressure(flash_1),
)
flush(stdout)

timed_call("vt_flash Float64 second") do
    flash_at(tank_temperature, tank_moles)
end
timed_call("vt_flash Float64 100 calls") do
    for _ in 1:100
        flash_at(tank_temperature, tank_moles)
    end
end

timed_call("ForwardDiff direct pressure d/dT first") do
    ForwardDiff.derivative(
        T -> pressure(model, tank_volume, T, [tank_moles]),
        tank_temperature,
    )
end
timed_call("ForwardDiff direct pressure d/dT second") do
    ForwardDiff.derivative(
        T -> pressure(model, tank_volume, T, [tank_moles]),
        tank_temperature,
    )
end

for (label, f, x) in (
    ("ForwardDiff vt_flash pressure d/dT", T -> flash_pressure(T, tank_moles), tank_temperature),
    ("ForwardDiff vt_flash pressure d/dn", n -> flash_pressure(tank_temperature, n), tank_moles),
    ("ForwardDiff vt_flash internal energy d/dT", T -> flash_internal_energy(T, tank_moles), tank_temperature),
    ("ForwardDiff vt_flash internal energy d/dn", n -> flash_internal_energy(tank_temperature, n), tank_moles),
)
    try
        value = timed_call(label * " first") do
            ForwardDiff.derivative(f, x)
        end
        println("  derivative=", value)
        flush(stdout)
        timed_call(label * " second") do
            ForwardDiff.derivative(f, x)
        end
    catch error
        println(label, " ERROR: ", typeof(error), ": ", sprint(showerror, error))
        flush(stdout)
    end
end
