# 1. Instructions

[unknown] - When explaining equations, define and name every variable immediately after introducing the equation.

[2026-10-14] - When reviewing academic writing, provide sufficiently critical feedback on grammar, sentence flow, clarity, and awkward phrasing. Do not overlook grammatical problems simply because the underlying ideas are sound.

[2026-10-20] - When discussing engineering decisions, focus on physical reasoning, practical engineering value, and how a method contributes to designing, testing, or improving real systems.

[2026-10-26] - Code should follow an "engineering-first, physically explicit" style.

[2026-10-26] - Keep the main physical calculation visible. Avoid hiding important physics behind generic abstractions, excessive helper functions, or unnecessary framework machinery.

[2026-10-26] - Use descriptive, domain-specific variable names. Long names are acceptable when they make physical meaning clearer.

[2026-10-26] - Prefer explicit intermediate variables and straightforward loops over compact programming idioms.

[2026-10-26] - Group code according to the sequence of physical operations being performed.

[2026-10-26] - Use comments to explain physical assumptions, geometry, and implementation reasoning.

[2026-10-26] - Major physics functions should be readable from top to bottom without requiring frequent navigation into other functions.

[2026-10-27] - Avoid aggressive indentation and unnecessary line wrapping in code.

[2026-10-27] - Keep assignments on a single line where practical instead of wrapping the assignment operator and expression onto separate lines.

[2026-10-27] - Include blank lines between logically distinct calculations and physical operations.

[2026-10-27] - Prefer code formatting resembling:

solid_fuel_specific_energy = u.volumetric_energy[fuel_cell_id] / solid_fuel_density

injected_fuel_specific_energy = solid_fuel_specific_energy - p.fuel_heat_of_pyrolysis

pyrolyzed_fuel_energy_flow = pyrolyzed_fuel_mass_flow_rate * injected_fuel_specific_energy

du.volumetric_energy_flow[fluid_cell_id] += pyrolyzed_fuel_energy_flow

[2026-10-27] - Avoid formatting resembling:

solid_fuel_specific_energy =
            u.volumetric_energy[fuel_cell_id] / solid_fuel_density
injected_fuel_specific_energy =
            solid_fuel_specific_energy - p.fuel_heat_of_pyrolysis

[2026-10-27] - Prioritize readability and explicit physical meaning over clever syntax or highly abstract software architecture.


# 2. Identity

[unknown] - Name: Will Martin.

[unknown] - Location: Canada. Previously associated with Regina, Saskatchewan; currently attending university in Saskatoon, Saskatchewan.

[unknown] - Age: Exact age not stored. Previously a Grade 12 International Baccalaureate student and currently a first-year university student.

[unknown] - Education: Previously completed or pursued International Baccalaureate coursework, including English HL, History IA, and Theory of Knowledge.

[unknown] - Personal interests include chemical engineering, reaction engineering, fluid mechanics, thermodynamics, propulsion, simulation, numerical modeling, experimental engineering, and building physical systems.

[unknown] - Particularly interested in reactor design, chemical kinetics, multiphase packed beds, Darcy flow, computational fluid dynamics, turbomachinery, rocket engines, and energy systems.

[2026-09-07] - First-year engineering student at the University of Saskatchewan in Saskatoon.

[2026-09-07] - Plans to pursue Chemical Engineering.

[2026-09-07] - Interested in laboratory-scale research and development, especially work combining physical equipment, experimental measurements, and physics-based mathematical models.

[2026-09-07] - Member of the University of Saskatchewan Space Design Team (USST).

[2026-10-20] - Identifies strongly with the role of a systems engineer who combines system architecture, mathematical modeling, uncertainty assessment, experimentation, and iterative design improvements.


# 3. Career

[unknown] - Long-term career goal: Work in chemical engineering research and development, preferably in a technically focused individual-contributor role rather than primarily in management.

[unknown] - Interested in eventually following a principal engineer or principal individual-contributor career path.

[unknown] - Plans to pursue a master's degree, with a PhD remaining a possibility depending on future interests and career requirements.

[unknown] - Prefers bench-scale hardware development and physics-based mathematical modeling.

[unknown] - Technical experience includes Julia programming, finite-volume numerical methods, differential equations, parameter estimation, thermal modeling, reaction kinetics, compressible flow, and experimental instrumentation.

[unknown] - Experience developing custom scientific software, particularly the Julia-based FVMFramework.

[unknown] - Experience with microcontrollers, ADCs, analog electronics, multiplexers, ultrasonic transducers, and data acquisition.

[unknown] - Experience with FreeCAD, mechanical design, pressure-system concepts, piping and instrumentation diagrams, and preliminary propulsion-system design.

[2026-09] - Joined the USST hybrid rocket engine project, where significant technical work and project organization were needed.

[2026-09] - Took a leading role in technical modeling and design decisions for the student hybrid propulsion project.

[2026-09] - Interested in eventually expanding USST propulsion development from hybrid rockets to liquid-fueled rockets and potentially gas turbines and ramjets.

[2026-10-20] - Expressed a preference for being a systems engineer with simulation capabilities and the ability to integrate experimental data into models.

[2026-10-20] - Particularly enjoys identifying uncertain aspects of a design, building models to investigate them, designing experiments to resolve uncertainty, and updating engineering decisions using experimental results.

[2026-10-27] - Increasingly interested in conventional chemical engineering applications that benefit from combining reduced-order models with targeted CFD simulations.


# 4. Projects

[unknown] - FVMFramework: Developed a custom finite-volume modeling framework in Julia. The framework provides a flexible foundation for discretizing and solving physical conservation equations. Developed the original framework personally without AI assistance, which contributed significantly to understanding how to use it as an engineering tool. Uses DifferentialEquations.jl for time integration and nonlinear solving. Continues to extend and apply the framework to thermal modeling, reaction engineering, compressible flow, and other engineering problems.

[unknown] - Methanol Steam Reforming Packed-Bed Reactor: Developed a numerical model for a catalytic packed-bed reactor involving methanol steam reforming, methanol decomposition, and the water-gas shift reaction. Hardware includes a catalytic packed bed, five axial thermocouples, and an internal heater wire. Developed and fitted a thermal model involving approximately 13 thermal parameters. Achieved good agreement between modeled and experimental temperature profiles. Chemical kinetics fitting was deferred while thermal behavior was being established. Experiments involved heater pulses and relatively warm flow conditions around 50 °C. Investigated insulation, heater construction, and thermocouple installation. Learned that Superwool insulation performed better than fiberglass, mineral-insulated heater wire was preferable, and thermocouple beads could puncture mica insulation.

[unknown] - Dialysis Diffusion Modeling: Used FVMFramework to model diffusion in a dialysis-related system, including temperature-dependent behavior described using Arrhenius relationships.

[unknown] - Esterification Kinetics Modeling: Used Julia and FVMFramework-related numerical methods to model esterification kinetics and Arrhenius temperature dependence.

[2025-11-24] - Numerical Solver Performance Investigation: Compared manually implemented ODE-based numerical calculations against ModelingToolkit-related approaches. Observed approximately 1,000,000 cells solved in 10 seconds using a manually implemented approach, compared with approximately 1,000 cells in 10 minutes using the ModelingToolkit approach under the tested conditions. This influenced preferences for direct numerical implementations.

[2026-09] - Ultrasonic Tomography for Packed Beds: Developed an experimental ultrasonic time-of-flight measurement system intended to map internal conditions in packed beds. Investigated ultrasonic excitation frequencies of approximately 5–40 kHz and transmitter pulses around 18–40 V. Hardware involved a Raspberry Pi Pico 2-family controller, ultrasonic transducers, approximately 1 nF transducer capacitance, TC4420 drivers, analog multiplexers, and an AD9226 12-bit ADC. Planned a custom PCB and improved signal conditioning. Intended to use measurements for ultrasonic tomography and spatial reconstruction within packed-bed systems.

[2026-09] - Ultrasonic Tomography Electronics: Investigated the HV2707, ADG5208F multiplexers, 74HC4067 multiplexers, 74HC595 shift registers, MCP4161 digital potentiometer, and AD9226 ADC. Accepted a -5 V supply for the receiver multiplexer and considered AC coupling with biasing. Planned shorter signal wires, improved ADC clock integrity, shielding, and reduced gain before unwanted transients.

[2026-09-13] - Ultrasonic Time-of-Flight Prototype: Achieved initial successful transmitter-to-receiver measurements. Example Tx1-to-Rx2 time-of-flight values included 148.072 µs, 129.845 µs, 140.149 µs, and 145.557 µs. Observed cooldown periods around 5.1–5.2 ms. Found that the AD9226 operated reliably near 2.5 MHz sampling, while some sampling frequencies above 4 MHz or below 2 MHz produced instability. Identified an ADC pipeline delay equivalent to seven sampling periods. Encountered mirrored upper ADC bitlines D11–D5 and investigated wiring and sampling issues.

[2026-09] - Compressible CFD Development in FVMFramework: Began implementing a density-based compressible Navier–Stokes solver. Initially selected HLLC numerical fluxes, with AUSM-family schemes considered for later development. Planned first-order reconstruction followed by MUSCL higher-order reconstruction, weighted least-squares gradients, viscous fluxes, species advection and diffusion, wall boundary conditions, and SST k–ω turbulence modeling. Investigated implicit monolithic solution approaches and low-Mach-number corrections.

[2026-09] - CFD Solver Verification: Developed a first verification stage for the inviscid, first-order HLLC solver. Tests included the Sod shock tube, uniform-state flux identity, internal-face conservation, and free-stream residual calculations. Used SSPRK43 with CFL 0.2 and FBDF with KrylovJL_GMRES and ILUZero. Tested using Julia 1.12.5, FVMFramework 0.1.0, ForwardDiff 1.3.2, OrdinaryDiffEq 6.108.0, and a 64-cell mesh. Confirmed agreement with analytical Euler fluxes, exact internal-face conservation, and near-machine-precision free-stream residuals. Found that FBDF required more steps and allocations in the tested configuration.

[2026-09] - CFD Solver Debugging: Identified issues involving viscous face gradients for velocity components and temperature being overwritten by averages. Also identified that density and internal-energy clamps around 1e-10 could conceal invalid physical states, and that an earlier MUSCL implementation generated nonphysical inlet behavior.

[2026-09] - USST Hybrid Rocket Prototype: Joined a small University of Saskatchewan Space Design Team project to develop a nitrous-oxide hybrid rocket engine. The first iteration was intended for static-fire testing rather than flight. Prioritized low cost, high probability of successful operation, useful experimental measurements, and lessons for later propulsion systems. The team had access to machining, welding, pressure testing, data acquisition, controls, and finite-element analysis capabilities. A hot-fire location was not yet established. Approximately five members were available. Funding discussions included an initial budget around $1,000 and potential funding increases through grants.

[2026-09] - Hybrid Rocket Propellant Selection: Considered N2O as oxidizer and HTPB, HDPE, or ABS as solid fuels. ABS offered 3D-printing convenience, while HDPE and HTPB remained attractive fuel candidates. Investigated grain regression, oxidizer-to-fuel ratio, pressure losses, and combustion performance.

[2026-09] - Hybrid Rocket System-Level Simulator: Developed a reduced-order ODE model for the hybrid rocket. Included oxidizer storage, valve flow, injector flow, combustion, fuel regression, nozzle performance, thrust, specific impulse, and propellant depletion. Investigated optimization of valve flow capacity, injector area, fuel-grain geometry, throat diameter, and oxidizer loading. Used numerical optimization objectives including thrust, oxidizer-to-fuel ratio, injector velocity, pressure-drop relationships, unburned propellants, combustion duration, and total impulse.

[2026-09-24] - Hybrid Rocket Preliminary Optimization: Manually identified a configuration with approximately 2.06 kg initial oxidizer mass, valve flow capacity factor 0.0001, injector area 0.0006 m², initial port diameter 0.010 m, additional port diameter 0.0132 m, grain length 0.60 m, and nozzle throat diameter 0.0233 m. The simulation predicted approximately 5.284 seconds of burn time, 895.6 N average thrust, an oxidizer-to-fuel ratio of 7.45, and specific impulse of 283.3 seconds. The desired O/F ratio was approximately 7.6. Earlier performance targets included approximately 1,000 N average thrust for five seconds.

[2026-09] - Hybrid Rocket Valve and Injector Investigation: Considered using pressure gauges before and after an adjustable valve to correlate pressure drop with oxidizer mass flow. Investigated electronically controlled valves costing more than $400, manually adjustable valves, injector-only flow restriction, and valve pressure losses. Considered future requirements for on/off control and possible re-ignition.

[2026-09] - Hybrid Rocket Detailed 1D/1.5D Model: Planned a spatially resolved model containing an oxidizer tank, post-valve region, pre-chamber, combustion chamber, fuel-grain port, solid-fuel region, post-combustion region, and nozzle. Considered density, total specific energy, axial velocity, and species concentrations as state variables. Investigated moving meshes for fuel regression, thermophysical properties through Clapeyron, and gravity orientation. Encountered negative cell masses caused by excessive time steps and investigated isoutofdomain, minimum time steps near 1e-5 seconds, and fallback thermodynamic calculations.

[2026-09] - Hybrid Rocket CAD and P&ID: Developed preliminary CAD and piping-and-instrumentation designs. Used FreeCAD and investigated separate material regions, assemblies, local coordinate systems, and component interfaces. Produced an initial P&ID in approximately two hours. Investigated machining and assembly of the injector plate, NPT threaded inlet connection, flanges, retaining components, threaded holes, through-holes, and bolted connections.

[2026-09-29] - Hybrid Rocket Scholarship Application: Began preparing a funding application under the Projects category for a Small-Scale Hybrid Rocket Prototype. Proposed a 400 N bench-scale prototype focused on learning through hands-on engineering, expanding USST propulsion capabilities, gathering experimental data, and informing a future flight-capable engine. Intended to accompany the application with a detailed budget spreadsheet.

[2026-10-03] - Hybrid Rocket Fuel-Web Modeling: Recognized that the fuel grain should retain a protective residual fuel web rather than regress completely to the liner. Investigated terminating the simulation based on a minimum remaining web thickness and corresponding final port diameter.

[2026-10-06] - Hybrid Rocket Injector Mechanical Design: Investigated whether the NPT inlet thread should be machined directly into the injector plate or into a separate component connected through a flange system. Explored screw retention, drilled and tapped holes, threaded inserts, through-bolts, and nuts. Expressed a preference for through-holes with nuts on the opposite side where practical.

[2026-10] - Ethylene Combustion Kinetics: Investigated the Westbrook–Dryer kinetic mechanism for ethylene combustion, including the reaction C2H4 + 2O2 → 2CO + 2H2O. Used a pre-exponential factor of 7.5895e16 in base SI concentration units and activation energy of 125,520 J/mol. Implemented a rate proportional to ethylene concentration raised to 0.1 and oxygen concentration raised to 1.65, with the reaction rate limited by the minimum of chemical and turbulence-controlled rates. Requested verification against the original published mechanism.

[2026-10] - Hybrid Rocket Fuel-Regression Module: Developed a physically detailed computational module involving axial planes, fuel regression, moving mesh geometry, ALE conservation, port-radius adjustments, and fuel injection into adjacent fluid cells. Wrote an initial version personally, recognized remaining errors, and considered delegating revisions to an LLM.

[2026-10-20] - CFD Method Selection and Validation: Experimented with SLAU2 numerical fluxes after an LLM recommendation. Used comparisons against an established SLAU2 implementation to validate behavior. Accepted that full understanding of every numerical implementation detail was not necessary provided the method was appropriately verified.

[2026-10-27] - Hybrid Rocket Full-CFD Feasibility Investigation: Combined compressible flow, SLAU2 fluxes, combustion, turbulence, and other physics into a more detailed hybrid rocket simulation. Found that the computational time required for full-system CFD was prohibitively large even after substantial optimization. Concluded that fully resolved CFD was impractical for routine full-motor design iteration. Shifted toward reduced-order system models for design and optimization, with CFD reserved for targeted component-level investigations.


# 5. Preferences

[unknown] - Strong preference for engineering problems involving physical systems, real equipment, experiments, and mathematical modeling.

[unknown] - Enjoys building original modeling tools because doing so develops intuition about how the numerical and physical systems work.

[unknown] - Prefers Julia as a primary scientific computing language.

[unknown] - Prefers custom FVMFramework implementations for many one-dimensional and two-dimensional engineering problems.

[unknown] - Prefers DifferentialEquations.jl for differential equation solving and numerical integration.

[unknown] - Generally avoids relying on ModelingToolkit, MethodOfLines, and Symbolics for primary modeling work.

[unknown] - Prefers finite-volume methods because conservation equations and fluxes have intuitive physical interpretations.

[unknown] - Finds parameter fitting and model calibration particularly useful and interesting.

[unknown] - Values physical interpretability and numerical verification when developing simulations.

[unknown] - Prefers models with explicit physical state variables, including density, momentum density, energy, species concentrations, and temperature.

[unknown] - Dislikes excessive mathematical abstraction when it obscures physical meaning.

[unknown] - Enjoys experimental apparatus design, instrumentation, and translating numerical predictions into measurable physical quantities.

[unknown] - Interested in developing systems that combine mathematical models with experimental measurements.

[unknown] - Prefers practical engineering tradeoffs over pursuing mathematical sophistication for its own sake.

[unknown] - Prefers a technically focused engineering career over a primarily managerial career.

[unknown] - Interested in research environments where simulations can directly influence laboratory experiments and hardware improvements.

[2026-09] - Initially interested in a multiyear student propulsion development progression involving a hybrid rocket in Year 1, liquid rocket in Year 2, jet engine or gas turbine in Year 3, and rocket-boosted ramjet in Year 4.

[2026-09] - For the hybrid rocket project, prioritizes achieving a successful first static fire, containing costs, and collecting useful experimental data over optimizing the first prototype for flight.

[2026-09] - Prefers simulation-driven preliminary design followed by physical testing and iterative improvement.

[2026-09] - Interested in using reduced-order models to optimize system-level design variables.

[2026-09] - Appreciates high-performance computational methods, especially when they permit faster design iteration.

[2026-09] - Interested in sliding-mesh methods, overset meshes, and harmonic balance approaches for turbomachinery simulations. Found the potential computational speedups of harmonic balance particularly interesting.

[2026-09] - Prefers using established physical relationships wherever possible rather than introducing unnecessary empirical complexity.

[2026-09] - Finds some highly specialized CFD numerical methods difficult to understand intuitively.

[2026-09] - Dislikes OpenFOAM's workflow and generally prefers working within a custom Julia-based modeling environment.

[2026-09] - Comfortable using LLMs to help implement numerical methods that are difficult or time-consuming to develop personally, provided their outputs can be verified.

[2026-09] - Experienced an afternoon energy slump, approximately between noon and 9 PM. Found that coffee around noon could improve work capacity for several hours.

[2026-09] - Attended mathematics lectures but found mechanics coursework less enjoyable. Preferred concentrated examination preparation and sometimes relied on strong International Baccalaureate fundamentals.

[2026-10-20] - Expressed uncertainty about how much time to invest in understanding specialized numerical methods such as SST k–ω turbulence modeling, SLAU2, second-order upwinding, and FEM.

[2026-10-20] - Enjoys understanding the physical meaning and interactions of simulation state variables more than investigating every detail of the underlying numerical machinery.

[2026-10-20] - Increasingly believes that LLMs can implement specialized numerical components more efficiently than personally developing everything from first principles.

[2026-10-20] - Values having developed FVMFramework personally because it provides a foundation for understanding and controlling subsequent simulations.

[2026-10-20] - Wants to distinguish between engineering work that deserves personal attention and implementation work that can reasonably be delegated to LLMs.

[2026-10-20] - Expressed strong agreement with the following description of an ideal engineering workflow:

"Here's the propulsion architecture I think makes sense. Here's why. Here's the model I built to investigate the uncertain parts. Here's what I don't trust in the model. Here's the experiment we're going to run to resolve that uncertainty. Here's what the data says. And here's how we're changing the design because of it."

[2026-10-20] - Reported that this systems-engineering workflow closely matched the work being done on the hybrid rocket project and that the work was highly enjoyable.

[2026-10-20] - Identified systems engineering with simulation and experimental-data integration as an especially appealing long-term career direction.

[2026-10-20] - Increasingly questions the practical value of investing heavily in highly sophisticated numerical machinery when existing methods can be implemented or verified with LLM assistance.

[2026-10-20] - Began reconsidering the importance of differentiable meshing and large-scale geometry optimization. Suspected that such approaches were impractical for near-term engineering applications.

[2026-10-20] - Shifted emphasis toward using FVMFramework to solve present engineering problems rather than expanding its numerical capabilities primarily for hypothetical future geometry-optimization workflows.

[2026-10-26] - Prefers engineering code that communicates physical operations clearly and directly.

[2026-10-26] - Prefers domain-specific variable names over short, generic identifiers.

[2026-10-26] - Prefers explicitly calculated intermediate quantities over densely nested expressions.

[2026-10-26] - Dislikes excessive helper-function decomposition when it forces the reader to jump between functions to understand one physical calculation.

[2026-10-26] - Prefers comments explaining why physical calculations are performed, especially when assumptions or geometric relationships are involved.

[2026-10-27] - Prefers moderate indentation and visually separated calculation blocks rather than aggressively indented, tightly packed code.

[2026-10-27] - Finds the ability to follow the sequence of physical calculations more important than minimizing code length.

[2026-10-27] - Comfortable accepting code generated by an LLM when its structure is sufficiently readable and the physical behavior can be validated.

[2026-10-27] - After observing the computational cost of full-system hybrid rocket CFD, concluded that detailed CFD should generally be reserved for smaller, constrained component-level problems.

[2026-10-27] - Became substantially more enthusiastic about reduced-order system models because their computational efficiency makes parameter optimization and repeated design iteration practical.

[2026-10-27] - Expressed excitement that the engineering workflow most personally enjoyable—system-level physical modeling, experimentation, and model improvement—also appeared to be a practical and effective way to solve engineering problems.

[2026-10-27] - Became less interested in ramjets and gas turbines after recognizing that full-system CFD with large-scale geometry optimization would often be computationally impractical.

[2026-10-27] - Became significantly more interested in conventional chemical engineering problems that benefit from combining targeted CFD with fast reduced-order models.

[2026-10-27] - Expressed concern about strongly coupled physical systems that cannot easily be decomposed into smaller subsystems or approximated sufficiently accurately using reduced-order models.

[2026-10-27] - Increasingly favors a modeling hierarchy in which inexpensive system-level models handle overall architecture and optimization, while detailed numerical models investigate localized phenomena or important uncertainties.

[2026-10-27] - Values uncertainty reduction, experimental validation, and practical design improvements as primary measures of a model's engineering usefulness.