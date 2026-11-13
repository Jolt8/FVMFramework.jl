Changes made to produce this plot:
- Instead of setting the tempreature of the combustion gases to room temeprature after the motor is assumed to be burned out, the temperature of the combustion gases are set to what they were during the time when CEA is being used
- This had little effect on the temperature of the combustion gases because the mass flux from the injector basically forces the combustion gases out so quickly that their temperatures do not matter
- IN ADDITION (CHANGES), I also made it so that when the motor is considered to be burned out, the oxidizer mass flux into the chamber stops
- This change made the pressure in the chamber drop to ambient almost immediatley, so a lot of evidence at least from this simulation points to post motor burnout temperatures not mattering that much and it seems like using the peak temperatures is a safe assumption
- Thus, for thermal analysis, I believe its safe to simulate only the time from when the motor starts to when it finishes burning

- EDIT: I also changed it so that all density changes in the fluid and heat exchange with the fluid is turned off once the motor burns out so we can observe how the temperature of the solids change after the motor burns out

- After making this change, I've observed that heat distribution after the motor burns out is quite significant, since a significant amount of thermal energy is transfered to the aluminum housing near the nozzle, the other bits of the alumium housing exchange heat with the extremely hot section of the aluminum housing, causing their temperatures to peak around 15 seconds after the motor is done firing

Also y1 is the bottom of the aluminum housing, y30 is the top