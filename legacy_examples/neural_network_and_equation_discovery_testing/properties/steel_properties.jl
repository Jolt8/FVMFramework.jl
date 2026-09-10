function get_steel_properties()
    steel_properties = ComponentVector(
        R_gas = 8.314u"J/(mol*K)",

        k = 50.0u"W/(m*K)",
        cp = 500.0u"J/(kg*K)",
        solid_rho = 7850.0u"kg/m^3",
    )

    return steel_properties
end