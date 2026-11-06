function get_air_properties!(vol, room_temp)
    air_density = 1.225u"kg/m^3"

    air_mass = air_density * vol

    initial_conditions = ComponentVector(
        mass = air_mass,
        temperature = room_temp,
        mass_fractions = normalize!(
            ComponentVector(
                mass_fractions = (
                    nitrous_oxide = 0.0,
                    ethylene = 0.0,
                    nitrogen = 0.781,
                    oxygen = 0.211,
                    carbon_dioxide = 0.0003,
                    water = 0.0
                ),
            )
        ),
    )

    properties = (
        molecular_weights = (
            nitrous_oxide = 44.01u"g/mol", 
            ethylene = 28.05u"g/mol", 
            nitrogen = 28.01u"g/mol", 
            oxygen = 32.00u"g/mol", 
            carbon_dioxide = 44.01u"g/mol", 
            water = 18.02u"g/mol"
        ),
        
        R_gas = 8.314u"J/(mol*K)",
        
        turbulent_kinetic_energy = 0.1u"m^2/s^2",
        specific_dissipation_rate = 100.0u"1/s",
        #k = 0.026u"W/(m*K)",
        prandtl_number = 0.705,
        #=
        diffusion_coefficients_ambient = (
            nitrous_oxide = 1.4e-5u"m^2/s",
            ethylene = 1.4e-5u"m^2/s",
            nitrogen = 1.8e-5u"m^2/s",
            oxygen = 1.9e-5u"m^2/s",
            carbon_dioxide = 1.6e-5u"m^2/s",
            water = 2.4e-5u"m^2/s",
        ),
        diffusion_coefficients_ambient_temperature = 21.13u"°C",
        diffusion_coefficients_ambient_pressure = 1.0u"atm"
        =#
        diffusion_coefficients = (
            nitrous_oxide = 1.4e-5u"m^2/s",
            ethylene = 1.4e-5u"m^2/s",
            nitrogen = 1.8e-5u"m^2/s",
            oxygen = 1.9e-5u"m^2/s",
            carbon_dioxide = 1.6e-5u"m^2/s",
            water = 2.4e-5u"m^2/s",
        )
    )

    special_caches = (

    )

    return initial_conditions, properties, cache_syms_and_units, special_caches
end