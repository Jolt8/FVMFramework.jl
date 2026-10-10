Ok, after doing a lot of optimization, fiddling around with different parameters manually and cross checking against historical hybrid data, I think I've found a good range of parameters that we could aim for 

# Optimized Parameters:

Valve flow capacity factor = (8 mm^2 - 12 mm^2) **Selected maximum = 10 mm^2**
- In terms of optimization, this one really doesn't matter that much, as long as the valve that we buy can adjust to a capacity factor within this range, we should be fine
- 10 mm^2 - 12 mm^2 seems to be a good max flow capacity factor for our purposes

Injector orifice area = (3.5 mm^2 - 5.0 mm^2) **Selected = 4.0 mm^2**
- This one seems like it's going to be hard to get right
- A common value from historical data on hybrids is around 4.0 mm^2
- However, I'm unsure if most teams using an injector area within this range had an adjustable valve before the injector like we're planning on using
- So I'd say something like 4.0 - 5.0 mm^2 would be good
- In any case, we should design this to be cheap and easy to replace if we get it wrong

Additional fuel grain void diameter = (8.2 - 10 mm) **Selected = 10 mm**
- This basically represents the ideal amount of the fuel grain burned during the burn. The port diameter is calculated based on the phenolic_liner_inner_diameter - 2 * desired_residual_fuel_web_thickness - additional_fuel_grain_void_diameter 
- From historical data, it seems like most go with something on the smaller end at around 8.2 mm
- Simulations tell me that it should be closer to 9 mm to get the right burn time of 5 seconds, but I'm going to trust the historical data because they actually have regression rate data while my regression rate prediction is pretty rough
- I think we'll just go for 10 mm because we want a longer burn time to consume all of our oxidizer
    - Note: this results in a burn time of around 5.85 seconds which I think is fine.
- NOTE: we always have to make sure that the area created by the initial port diameter is greater than the nozzle throat area (which it currently is)

Fuel grain length = (280 mm - 310 mm) **Selected = 310 mm**
- This still seems pretty uncertain to me
- Based on historical data, it seems like most teams go with something in the range of 240 mm - 280 mm
- However, based on simulations, it seems like something closer to 300 mm would be better
- Also, I'm extra unsure because Purdue's hybrid that seems very similar to ours used a fuel grain length of around 240 mm, which is way off of what the simulation is telling me
- I think this needs a little bit more research, but I'm leaning towards something around 310 mm because the burn time we're aiming for is longer than the historical data motors that were aiming for a burn time of around 5 seconds
- We can also vary the fuel_grain_length and additional_fuel_grain_void_diameter if we want a different burn time
- This might need to be longer to consume all the oxidizer (probably around 375 mm to have only 0.1 kg of oxidizer left assuming a additional_fuel_grain_void_diameter of 10 mm)
- Given this, I think we should leave a bit of extra room in the design to have the fuel grain length eat into the post-combustion section's lenght a little to hit that 375 mm target if we need to
- I think a max fuel grain length of around 360 mm is reasonable to leave space for
- 375 mm is probably a bit excessive

Nozzle throat diameter = (14.0 mm - 15.0 mm) **Selected = 14.0 mm**
- I'm pretty confident on this one both from historical data and simulations
- The Purdue one used a 14 mm nozzle diameter, and the nozzle diameter that gives the best results seems to be around 14 mm for the simulation as well
- I think we'll go with 14 mm

Expansion ratio = (4.0 - 4.5) **Selected = 4.0**
- Based on historical data, most teams go with something in the range of 4.0-4.5
- When the simulation is allowed to adjust this, it adjusts to something around 7.0 suprisingly
- I think this is due to the fact that the model uses the vacuum isp instead of the sea-level isp that depends on chamber pressure, O/F ratio, and the expansion ratio
- Thus, I locked in the expansion ratio to 4.0, so we'll go for 4.0
- We might be able to do something slightly higher


# Other parameters that are locked in:

Initial tank oxidizer mass = 1.134u"kg" 
- Based on a 2.5lbs COTS tank

Phenolic liner inner diameter = 33.32 mm
- Based on a 1.312 inch standard COTS size

I'm still a little bit skeptical of the simulation because almost all the runs with parameters similar to the ones above leave about 0.3 kg of unused nitrous oxide

This one however:
```julia
viewable_system_design_loss(
    ComponentVector(
        valve_flow_capacity_factor = 1.0e-5,
        injector_orifice_area = 4.0e-6,
        additional_fuel_grain_void_diameter = 0.012,
        fuel_grain_length = 0.31,
        nozzle_throat_diameter = 0.014,
    ), properties_unitless
)
#These parameters work very well, giv eus a burn time of about 6.5 seconds 
#and result in a final cummulative impulse of aorund 2900 N*s
#It also uses up much more of the avaliable oxidizer while the others don't
```
Only leaves behind 0.1 kg of nitrous oxide while achieiving a cummulative impulse of around 2900 N*s. This simulation including all the others usually have chamber pressures of around 1.8 MPa (peak of around 2.0 MPa), which is slightly lower than the commonly cited 2.0 MPa for most student hybrids. This one has an oxidizer to fuel ratio of around 7.8 which is in the right range. Honestly, I think this one might be one of the most trustworthy simulations.

This one:
```julia
viewable_system_design_loss(
    ComponentVector(
        valve_flow_capacity_factor = 1.0e-5,
        injector_orifice_area = 6.0e-6,
        additional_fuel_grain_void_diameter = 0.009,
        fuel_grain_length = 0.31,
        nozzle_throat_diameter = 0.014,
    ), properties_unitless
)
```
Has a chamber pressure that's much more inline with historical data (around 2.0 MPa) and has a cummulative impulse of 2617 N*s. It also has an injector velocity of around 50 m/s which seems perfect. This one leaves behind 0.2 kg of nitrous oxide. This one also has an oxidizer to fuel ratio of around 8.8 which is a little too high.

I'm staring to get a feeling that my ISP calculations are wrong, 2400 N*s of impulse from around 0.9 kg of nitrous oxide doesn't seem realistic, although under ideal circumstances you could theoretically achieve around 2600 Ns with this much nitrous oxide. I wonder if I could just multiply the get_isp() funciton by about (200 / 260) to roughly correct it. That being said, I think everything else in the simulation (pressure, mass, O/F ratio, injector velocity, etc.) is very realistic.

# Parameters Used in These Simulations:
```julia
properties = ComponentVector(
    #Overall rocket properties
    desired_oxidizer_to_fuel_ratio = 7.6,
    desired_burn_time = 5.0u"s", #Flexible, anything between 5 and 8 seconds is fine
    desired_impulse = 2500.0u"N*s", #this will be calculated based on desired_average thrust and desired_burn_time later
    target_injector_pressure_drop_to_adjustable_valve_pressure_drop_ratio = 2.0, #not an optimized parameter, but the ideal choice is hard to know
    
    #Tank
    u0_tank_oxidizer_mass = 1.134u"kg", #No longer optimized because we're using a COTS N2O Tank
    tank_temperature = 21.0u"°C",

    #adjustable valve
    valve_flow_capacity_factor = 10.0u"mm^2", #a result of the final valve we choose #we should optimize this to get a good ideal of what we want
    ```
    Note: this is used to calculate adjustable valve flow like this:

    p.adjusted_valve_flow_capacity_factor = valve_opening * p.valve_flow_capacity_factor
    
    volumetric_flow = p.adjusted_valve_flow_capacity_factor * sqrt((p.tank_pressure - p.mid_section_pressure) / p.tank_density)

    Thus, we're going to have to convert this to a flow coefficient / flow factor to pick a commercial valve
    Kv is approximately Kv = 0.036 * valve_flow_capacity_factor (only works if valve flow capacity factor is in mm^2)
    In our case, Kv is 0.36 m^3/h at full opening based on a valve_flow_capacity_factor of 10.0 mm^2
    ```julia
    adjusted_valve_flow_capacity_factor = 0.0u"mm^2",
    valve_opening = 1.0,

    #Mid section
    mid_section_pressure = 1.0u"atm", #this determines the initial kg of nitrous oxide in the mid_section
    mid_section_temperature = 21.0u"°C",
    mid_section_volume = 100.0u"cm^3", #Changed from 10u"cm" to 100u"cm" because I wanted to decrease solver stiffness for debugging purposes

    #Injector properties
    injector_discharge_coefficient = 0.7,
    #injector_number_of_orifices = 20, #optimized
    #injector_orifice_diameter = 1.5u"mm", #optimized
    injector_orifice_area = 12.0u"mm^2", #optimized
    target_injector_velocity = 30.0u"m/s", #Would be fine with anything from 30-60 m/s

    #Chamber
    u0_chamber_nitrous_oxide_mass_fraction = 1.0,
    u0_chamber_hdpe_mass_fraction = 0.0,
    chamber_pressure = 1.0u"atm", #this determines the initial kg of nitrous oxide in the chamber
    chamber_temperature = 21.0u"°C",
    chamber_volume = 212.0u"cm^3", #this is very imperfect, we should update this in the future
    target_chamber_pressure = 20.0u"bar",

    #Fuel grain
    #u0_fuel_grain_void_diameter = 3.0u"cm",
    fuel_mass = 0.0u"kg", #this will be derived by substracting the volume of the cylinder formed by the u0_fuel_grain_void_diameter by the final_fuel_grain_void_diameter and then multiplying by the fuel density
    fuel_density = 950.0u"kg/m^3",
    #fuel_regression_rate = 0.5u"mm/s",
    additional_fuel_grain_void_diameter = 0.8u"cm", #optimized
    phenolic_liner_inner_diameter = 33.32u"mm", #not optimized, we're going to be using a COTS phenolic liner, so we're just going to use the ID that the manufactuerer specifies
    final_fuel_grain_void_diameter = 0.0u"mm", #the actual OD of the fuel grain should be about 33.32 mm #Update: this is not going to be defined by the phenolic_liner_inner_diameter - 2 * desired_residual_fuel_web_thickness
    fuel_grain_length = 27.0u"cm", #optimized
    desired_residual_fuel_web_thickness = 2.0u"mm", 

    #Fuel Grain Empirical Parameters
    #These are for when fuel_regression is measured in mm/s and oxidizer_mass_flux is measured in g/(cm^2*s)
    fuel_regression_coeff_a = 0.248,
    fuel_regression_coeff_n = 0.331,

    #Fuel properties
    fuel_lowering_heating_value = 47.42u"MJ/kg",
    stoichiometric_oxidizer_fuel_ratio = 9.41,
    fuel_heat_of_pyrolysis = 2.3u"MJ/kg",
    fuel_surface_temperature = 1300.0u"K",
    combustion_efficiency = 0.9,

    #Nozzle
    nozzle_throat_diameter = 15.0u"mm", #optimized
    nozzle_discharge_coefficient = 0.9,
    expansion_ratio = 4.0
)
```






