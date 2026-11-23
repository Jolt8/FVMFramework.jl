# Why was axisymmetric FVM the wrong approach this for this project?

Here are a couple issues that I remember from working on this simulation script:

Problems and Improvements in the future:

# The model was way too slow
Problem:
- While cold flow runs and even runs where the igniter physics were turned on took about 2-5 minutes to run, the arrhenius chemistry made everything slow to a crawl with the simulation likely taking around 12 hours (!) to finish, which is just way too long to get any valuable engineering information out of this model.

Improvements:
- In general, I feel like for a lot of the problems that I work on in the future, Arrhenius chemistry is just not going to be practical for almost any kind of simulation4

- I mean, as I understand it, two major things that make ODEProblems slow are:
    1. Wildly different timescales between coupled physics 
        - In this model, heat transfer in solid materials was a very long timescale problem
        - SLAU2 advection and convection through the nozzle was about in the middle in terms of timescale
            - SLAU2 was also really not the choice for this preliminary model because I really don't think that we had the data to inform the model enough to have it actually predict any kind of unwanted transients during startup, since this is really the only thing that SLAU2 gets us. 
        - Arrhenius chemistry was extremely fast, slowing down the problem way too much compared to these other models
    2. Physics that requires conditionals and forms a negative feedback loop with itself
        - Arrhenius type reaction physics are exceptionally bad for this because they typically can't handle when concentrations dip into the negative territory briefly, thus requiring a conditional, and then they form a negative feedback loop with themselves where they react on one step, bring the concentrations down too low for the next step for them to still be active, causing the solver to have to reject tons of steps and take significantly smaller steps 
        - A smooth cutoff might be able to aid in this, but I haven't had that many great experienced even with infinitely differentiable ramps in the past as they typically perform worse than conditionals 

- I sometimes wonder if the method of doing a split ODEProblem where stiff/short timescale physics gets a stiff solver like Rodas5P and then non-stiff parts of the problem get a non-stiff solver like FBDF, but I think we should just avoid this problem in the future

- The solution to this is probably going to be more heavily relying on gibbs-free energy minimization chemistry because you can either make it instant (avoids reaction timescale issues almost entirely) or have it follow a relaxation rate that can be easily customized. 

# Not having multiphase SLAU2, led to some hard choices in terms of modeling 
Problem:
- Since simulating multiphase SLAU2 would take significantly more modeling effort, injector cooling was neglected with a smooth ramp once the pre-chamber temperture went below a specific threshold

Improvements:
- I don't think there's a lot to improve here, I'm rather happy with this because implementing multiphase SLAU2 would've been hell

# Using SLAU2 in general makes the problem much harder to work on for me as an engineer
Problem:
- Since: 
    1. I had an LLM implement the SLAU2 physics
    2. SLAU2 and compressible navier stokes code is genuinely incomprehensible without if you haven't spent a significant amount of time on it (more time than it's actually worth) 
    3. Even if the two factors above did not exist, it still does not take away from the fact that writing code that respects compresisble navier stokes takes absolutely forever and is never that intuitive
- Working on this project was a pain, and for a significant majority of the physics (especially the regression model) I relied on a LLMs to either: 
    1. edit my code so significantly that it was hardly recognizable 
    2. implement something for me so that I barely even knew how it worked. 
- This is terrible mostly because maintaining technical ownership and writing the code myself for these engineering projects specifically is extremely important because I ultimately want to be the person making the assumptions and writing down comments where necessary. Furthermore, LLMs just write in a style that really goes against my own coding style, making it usually pretty hard to read.

Improvement:
- This is kinda sad to say, but I think for overall systems models, FVM with any kind of compressible navier stokes is just never going to be the kind of model fidelity that I reach for. 
- Once more, I'm not even sure if SLAU2 is going to be right to use for almost any kind of future problem that I work on because even in the cases where I'm working on specific components, I feel like empirical correlations would often be more accurate and easier for me to work with as an engineer
- I think it's going to be quite hard in the future to find examples where SLAU2 or any form of compressible or incompressible navier stokes are going to be the best choice

# Models where there's too much uncertainty of whether or not it's going to work or not really encourages LLM usage
Problem:
- While working on this model, due to the large amount of time that it was taking to complete and because I was using multiple types of physics that I just really didn't have a lot of confidence in such as:
    1. The moving fuel grain (which moved the geometry for a significant number of cells other bodies)
    2. The injector physics and the fact that I was using a density based approach for incompressible fluids
    3. The igniter's physics and how it had to slowly transition from a solid to a fluid_domain cell

Solution:
- Again, I really think that we just can't commit to any kind of detailed simulation before we have actual test data 
- Thus, I think of a majority of my future projects are likely going to use a mix of: 
    - Overall systems models (like the transient design script) for determining rough values at the start of design and acting as a source of truth
        - Will not use FVMFramework and should likely be custom built every single time
        - Not sure if uncertainty quantification belongs here
    - Simplified reduced order models (like the simplified ROM model currently in this directory) for identifyign current weak points in the design
        - the simplified_hybrid_digital_twin showed me how bad heat slowly soaking back from the graphite nozzle into the aluminum was going to be 
        - These ROM models also provide a great way to allow an LLM to interface with your model and collect the data that it needs to do its own analysis and report its findings 
    - Component level design scripts that use already established empirical models with uncertainty quantification 
        - I have not complete the nozzle one so far, but with how reliable empirical correlations usually are, I'm very confident that it would produce very useful results especially with uncertainty quantification partly covering our asses 




