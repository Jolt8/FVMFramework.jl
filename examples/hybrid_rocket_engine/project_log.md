


# Project log from simplified_hybrid_digital_twin.jl:



## Switching the solver to correctly use the T_comb temperature after cea is turned off:
- Instead of having the fluid_domain temperature be set to basically room temperature after CEA is turned off, it is instead set to the CEA temperature, and then allowed to cool on its own
- It appears as if making this change does not make a significant difference in the temperatures experienced by components inside the motor
- This is mainly due to the fact that the nitrous oxide coming in from the injector "washes out" the still hot combustion gases out of the chamber almost immediately, causing this modified version to basically resemble the version where the combustion gases are set to room temperature after the fact

- IN ADDITION (CHANGES), I also made it so that when the motor is considered to be burned out, the oxidizer mass flux into the chamber stops
- This change made the pressure in the chamber drop to ambient almost immediatley, so a lot of evidence at least from this simulation points to post motor burnout temperatures not mattering that much and it seems like using the peak temperatures is a safe assumption

- EDIT: I also changed it so that all density changes in the fluid and heat exchange with the fluid is turned off once the motor burns out so we can observe how the temperature of the solids change after the motor burns out

- After making this change, I've observed that heat distribution after the motor burns out is quite significant, since a significant amount of thermal energy is transfered to the aluminum housing near the nozzle, the other bits of the alumium housing exchange heat with the extremely hot section of the aluminum housing, causing their temperatures to peak around 15 seconds after the motor is done firing

- In conclusion, while almost all parts experience their peak temperature once combustion stops, if something like the aluminum housing has a connection of the nozzle where a severe thermal gradient is allowed to form compared to the aluminum housing elsewhere, they will reach their peak temperatures later