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
        molar_masses = (
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
    )

    cache_syms_and_units = (
        #updated based on temperature, density, and composition:
        density = u"kg/m^3",
        pressure = u"Pa",
        speed_of_sound = u"m/s",

        #updated based on temp:
        k = u"W/(m*K)",
        mu = 1e-5u"Pa*s",
        
        #updated based on temperature, density, and composition:
        cp = 1.005e3u"J/(kg*K)",
        cv = 718.0u"J/(kg*K)",

        #updated based on cp and cv
        gamma = u"1",
        
        #updated based on composition:
        mw_avg = u"g/mol",

        molar_concentrations = (
            nitrous_oxide = u"mol/m^3",
            ethylene = u"mol/m^3",
            nitrogen = u"mol/m^3",
            oxygen = u"mol/m^3",
            carbon_dioxide = u"mol/m^3",
            water = u"mol/m^3",
        ),
        diffusion_coefficients = (
            nitrous_oxide = u"m^2/s",
            ethylene = u"m^2/s",
            nitrogen = u"m^2/s",
            oxygen = u"m^2/s",
            carbon_dioxide = u"m^2/s",
            water = u"m^2/s",
        )
    )

    special_caches = (

    )

    return initial_conditions, properties, cache_syms_and_units, special_caches
end