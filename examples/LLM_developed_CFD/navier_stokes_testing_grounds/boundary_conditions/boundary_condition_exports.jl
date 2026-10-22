include(joinpath(@__DIR__, "bcs/characteristic_farfield_flux.jl"))
export characteristic_farfield_flux!

include(joinpath(@__DIR__, "bcs/no_slip_adiabatic_wall_flux.jl"))
export no_slip_adiabatic_wall_flux!

include(joinpath(@__DIR__, "bcs/supersonic_outflow_flux.jl"))
export supersonic_outlet_flux!

include(joinpath(@__DIR__, "bcs/slip_wall_flux.jl"))
export slip_wall_flux!, symmetry_flux!

include(joinpath(@__DIR__, "bcs/supersonic_inlet_flux.jl"))
export supersonic_inlet_flux!