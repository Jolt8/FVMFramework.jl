# Overall Takeaways

## Typical amount and progression of model fidelity throughout a project:
- I think of a majority of my future projects are likely going to use a mix of: 
    - Overall systems models (like the transient design script) for determining rough values at the start of design and acting as a source of truth
        - **for initial sizing and values that other scripts can work off of**
        - Will not use FVMFramework and should likely be custom built every single time
        - Not sure if uncertainty quantification belongs here
    - Simplified reduced order models (like the simplified ROM model currently in this directory)
        - **for identifying current weak points in the design**
        - the simplified_hybrid_digital_twin showed me how bad heat slowly soaking back from the graphite nozzle into the aluminum was going to be 
        - These ROM models also provide a great way to allow an LLM to interface with your model and collect the data that it needs to do its own analysis and report its findings
    - Component level design scripts that use already established empirical models with uncertainty quantification 
        - **For more intensive design and making sure that something will work even with uncertainty**
        - I have not complete the nozzle one so far, but with how reliable empirical correlations usually are, I'm very confident that it would produce very useful results especially with uncertainty quantification partly covering our asses 

## Do not build for data you do not have yet
As an example:
- Due to my ideal that I would one day be able to link the regression observed in the fuel grain to what happened during the test fire, I constantly wanted to use a model that could explicitly track this moving fuel grain over time. 
- In hindsight, this really was not the right choice for a few reasons:
    1. It adds an insane amount of extra complexity that made me rely on an LLM to get it working properly
    2. While I haven't confirmed this for sure, I don't think that the moving fuel grain would've had a significant impact on the information that we were trying to extract out of the model

Improvements:
- Don't preemptively build a model that can interface with future collected data because it's almost always better to do something when it's actually needed rather than immediately because you think it might be useful in the future
- While it's hard to do this, always critically evaluate whehther or not a mechanism that is complicated actually adds that much to the problem and don't get blinded by the belief that (more accurate physics == always better)

## LLM Policy
- Use LLMs as much as you want for information retrieval and ideas/suggestions
    - Also, for more specific engineering questions that would take a long time to answer and are moderately valuable, use Claude to run its own code (in whatever language it wants) to verify some of your conclusions and potentially identify things that you missed 
- Do not use LLMs for writing the core systems engineering models, but they can still probably be useful for creating some abstractions if you know for an absolute fact (and be very sure about this) that learning it is not worth your time
    - For example, an LLM wrote the FreeCAD workbench that allowed me to make the ROM model very easily
    - Also, in hindsight, I don't think that an LLM writing most of the machinery behind SLAU2 caused that big of a loss of ownership when I use SLAU2 in other scripts, but whenver I use the LLM to write the physics for a specific problem because I don't currently understand SLAU2 right now that well, then I lose ownership
    - I think the same would go for moving geometries, if it created some higher level abstractions and then it was my job to use them where they fit, I don't feel like I would feel the loss in ownership
    - This is mostly because these are generic tools and they don't immediately mean that the LLM made a modeling decision
    - Summed up: LLMs can write code with no assumptions in it (tooling) or with literature-standard assumptions I can look up. They don't write code that makes new modelling decisions.

## Most Important Transient Event
- For all the physics that you choose to use, make sure that they're not faster than the fastest thing you actually care about and if they're faster make them instant
- For example, SLAU2 was not useful because we were not trying to capture millisecond chamber startup transients
- The mid_section having its own volume was not useful once again because we did not need to capture its millisecond transients

## Start every model with the question that it's answering
- Always make sure that before you build a model it can answer this question: This model exists to decide X; it needs accuracy Y on quantity Z
- Try not to develop general purpose models (except for maybe the ROM models) because they're typically less useful

## Models as logbooks
- Models are extremely valuable as means of informing teammates or LLMs about the problem that you're working on because they encode pretty much every engineering decision you've made

# Model 1: The Design Script

## Stats
- Runtime: less than a second (probably around 0.3 seconds for most runs)
    - Perfect runtime for initial optimization
- Development Time: About a week and a half
- Decisions informed: probably the most of all models (that not saying much for an initial design model though)

## Things I liked

### I was able to maintain complete ownership of this model and I learned a lot about the problem while doing so
- Self explanatory, since this model is pretty simple all things considered, I did not have to use an LLM that much so I retained complete ownership of the script 

### Comparing to historical data 
- Using historical data gathered using NotebookLM deep research, I was able to compare the results that my solver came up with to historical data and they matched pretty closely after a couple tweaks, so I really trusted the model after that and the parameters that it gave me.

### Data gathered by this model
- This model pretty much became the source of truth used for all other models so this was probably the most important model of this entire project

Overall, compared to the previous models that I built, I was certainly the most happy with this one. 

## Problems and Improvements

### Passing dual numbers into the solve was a pain in the ass, forcing me to use FiniteDiff
Problem:
- Since I didn't use the common FVMFramework infrastructure, I couldn't easily leverage the work I had put in before to make sure that duals cleanly make it into the solver

Improvements (two of them here):
- Be more conscious of them in the future and make sure that they can be passed in relatively easily
- Just make most design scripts AD incompatible 
    - This isn't that bad because I've found myself using Tsit5 and BBO optimization more frequenlty with this style of script anyways


### Clapeyron.jl's uv_flash just sucks in general and causes some compilation errors that makes compiling the overall system take ages
Problem:
- For some reason, while using uv_flash or ph_flash with Clapeyron.jl, I had to constantly use this kind of weird function closures to get maximum performance

```julia
f_closure = let
    om = oxidizer_model
    cm = chamber_model
    p_axes_local = p_axes
    u_axes_local = u_axes

    (du, u, p, t) -> system_ode!(du, u, p, t, om, cm, p_axes_local, u_axes_local)
end
```

- This was extremely obnoxious because it basically mean that whatever parameters I entered into the solver that were not part of the optimized parameter were fixed until I restarted my repl
- It also made changing functions quite difficult 

Solution:
- I really don't have a good solution here, the only thing I could think of would be to open an issue on Clapeyron.jl
- I did not solve this really and basically just had to tough it out while working on the script 
- Check VT_FLASH_COMPILATION_FINDINGS.md for more info
- YES I GOT IT! Just tabulate some of the thermodynamics (especially for pure systems)


### The optimization routine is pretty tedious
Problem:
- If I ever wanted to adjust parameters or try a different design, there was quite a lot of friction in doing that because it would basically mean another 3 hour work session on the model that was starting to feel pretty fragile near the end, the fragility and friction can generally be attributed to:
    - flash calculations failing randomly
    - The injector was pretty unstable since the mid section and pre-chamber were modeled as separate volumes again
    - Changing non-optimized parameters required a julia restart
    - So many safeguards were required to stop the solver from crashing randomly

Solutions:
- Really not sure, the only thing that I could really hope to get rid of would be the Clapeyron.jl bs, but other than that, I think this is just something that has to be accepted for non-deterministic/not well defined/unconstrained design problems
    - non-deterministic as in most of the physics do not flow one way and are tightly coupled, requiring the use of actual optimization methods rather than being able to just detemine that this flow rate of nitrous oxide will result in this thrust directly
- The only other solution would be to just not use non-ideal EOS, but in the future we're going to have to determine whether or not ideal EOS is good enough
- This gets into the problem of sizing (required performance determines geometry) and rating (geometry is checked if it meets required performance), sizing is typically very fast to generate data on and should be used to constrain the design space, then rating should be used to find the best candidate



## Things left to decide:

### Common Design Parameter Interface
- I should've probably provided a better means of interfacing with the results of this solver rather than just purely referencing the values stored in my head
Notes:
- This one is slightly debatable
- On one hand, being able to access this solver directly would make sure that what we're designing for would always stay consistent with this source of truth.
- On the other hand, this project's needs, requirements, design specs, etc. have slowly morphed over time into something slightly different
    - Unsure if the different is better or worse
    - Also, using values stored in my head is much quicker
    - That being said, it would not be hard at all to import a common ComponentVector storing all the design parameters for this problem similar to the way I've used the freecad table data with the cad_data ComponentVector
- EDIT: using a common ComponentVector is better

# Model 2: The Axisymmetric FVM Model

## Stats
- Runtime: 
    - 12 hours for the full hot fire
    - About 2-5 minutes for the cold flow with the igniter
- Development Time: About a week and a half
- Decisions informed: basically nothing other than that a less detailed model is needed

## Problems and Improvements in the future:

### The model was way too slow
Problem:
- While cold flow runs and even runs where the igniter physics were turned on took about 2-5 minutes to run, the arrhenius chemistry made everything slow to a crawl with the simulation likely taking around 12 hours (!) to finish, which is just way too long to get any valuable engineering information out of this model.

Improvements:
- In general, I feel like for a lot of the problems that I work on in the future, non-integer-order Arrhenius chemistry is just not going to be practical for almost any kind of simulation for the kind of timescales that most engineering problems require

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

- NOTE: check brainstorming.md in the new_reaction_mechanism_brainstorming directory because I think I have a solution to this

- Corrections
    - Ok, kinetics are not impractical, especially when mixing is the main limiter instead of the reaction itself (Da number ~<= 1)
    - Also, fractional reaction order were one of the most significant contributors to this pain

### Not having multiphase SLAU2, led to some hard choices in terms of modeling 
Problem:
- Since simulating multiphase SLAU2 would take significantly more modeling effort, injector cooling was neglected with a smooth ramp once the pre-chamber temperture went below a specific threshold

Improvements:
- I don't think there's a lot to improve here, I'm rather happy with this because implementing multiphase SLAU2 would've been hell

### Using SLAU2 in general makes the problem much harder to work on for me as an engineer
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

### Models where there's too much uncertainty of whether or not it's going to work or not really encourages LLM usage
Problem:
- While working on this model, due to the large amount of time that it was taking to complete and because I was using multiple types of physics that I just really didn't have a lot of confidence in such as:
    1. The moving fuel grain (which moved the geometry for a significant number of cells other bodies)
    2. The injector physics and the fact that I was using a density based approach for incompressible fluids
    3. The igniter's physics and how it had to slowly transition from a solid to a fluid_domain cell

Solution:
- Again, I really think that we just can't commit to any kind of detailed simulation before we have actual test data 
- See 'Typical amount and progression of model fidelity throughout a project:' in the overall takeaways above

# Model 3: Hybrid Digital Twin

## Stats
- Runtime: 
    - 20-30 minutes for the hot fire although for some reason I feel like it was longer
- Development Time: About a week
- Decisions informed: basically nothing other than that a less detailed model is needed

## Problems and Improvements in the future:

## A moving fuel grain was used
Problem:
- Due to my ideal that I would one day be able to link the regression observed in the fuel grain to what happened during the test fire, I constantly wanted to use a model that could explicitly track this moving fuel grain over time. 
- In hindsight, this really was not the right choice for a few reasons:
    1. It adds an insane amount of extra complexity that made me rely on an LLM to get it working properly
    2. While I haven't confirmed this for sure, I don't think that the moving fuel grain would've had a significant impact on the information that we were trying to extract out of the model
- Since I used an LLM to get this implemented faster, I really lost ownership of certain parts of this model and I lost interest in it

Improvements:
- Don't preemptively build a model that can interface iwth future collected data because it's almost always better to do something when it's actually needed rather than because you think it might be useful in the future
- While it's hard to do this, always critically evaluate whehther or not a mechanism that is complicated actually adds that much to the problem and don't get blinded by the belief that (more accurate physics == always better)
- If you ever feel tempted to use an LLM specifically for writing the code for a systems engineering model, just DON'T, it's not worth it, I don't know if it's just because I'm overly petty and I hate when I lose ownership of my code, but it genuinely feels that for systems models even a slight amount of lost information makes everything significantly harder to work on 
- Perhaps we should develop a more generally applicable moving mesh methodology (even if it's generated with the help of an LLM) so that, even though it's abstracted, I get to choose where, when, and why something is or isn't used 
    - This might be impossible because moving meshes typically require that momentum and other conservation of energy implications be dealt with (at least for SLAU2)

## Arrenhius chemistry was used (not going to write this out fully, just see my explanation for the FVM script used in this project)

# Model 4: Simplified Hybrid Rocket 

## Stats
- Runtime: 
    - 1-2 minutes
- Development Time: 3 days (since it took advantage of the things established in model 3)
- Decisions informed: where some weakpoints in the design were (such as the heat absorbed by the graphite soaking back into the aluminum)

## Things I liked about this project

### An LLM realized that the injector could be merged to no longer require the mid section
- This is what I actually absolutely love LLMs for, they are extremely creative when it comes to proposing solutions that really wouldn't of ever occured to me
- Because we did this, we no longer had to deal with the extremely sensitive mid_section transient, making it significantly easier to solve the problem 

### No moving mesh was used
- Since this script did not use a moving mesh, it was significantly easier to get a working model and I never really had to rely on an LLM to code anything for me, allowing me to take full ownership of the model

### NASA CEA was used
- Compared to Arrhenius chemistry which was extremely stiff for this problem, NASA CEA was extremely practical here and worked very well

## Problems and Imporvements

### Using rocketcea and being forced to use HP mode
Problem:
- Since we're using rocketcea and it's forced to be in HP mode, it does not map especially cleanly to FVM

Improvements:
- Instead of using rocketcea in the future, we should ideally create our own CEA engine in Julia that utilizes its own gibbs or helmholtz free energy minimization routine combined with NASA's 9 coefficient property polynomials 
- Another option would just be to create an application that takes in the data from NASA CEARUN, but I would like to develop my own to make it significantly more customizable and so I understand it better

# Sub Model 1: Nozzle Thermals
## Stats
- Runtime: TBD
- Development Time: 3 days
- Decisions informed: TBD

## Things I liked:

### Uncertainty Quantification
- The uncertainty quantification and how that general system can be generally used to design something that will have the highest chance of success in addition to giving you information on what parameters are going to be the most valuable to retrieve from experiemtal data is amazing
- margin ratio = (capacity − demand) / combined uncertainty of capacity and demand

### Use of Empirical Models
- Due to the use of many empirical models, I trust this model significantly more than basically all of the models that I've used so far for this specific problem
- I didn't really choose what empirical models I was going to use and I pretty much just left that to Claude, but I didn't feel a significant loss of ownership despite this because it's just information retrieval
    - I'm also fine with the LLM estimating the uncertainty for these models based on consensus online and maybe just latent knowledge it has retained in its weights

### Writing code myself is a great teacher of the problem's structure
- Forcing myself to write the code really helped me understand this specific problem very well 
- I think this script was the final nail in the coffin showing just how important it is that I do not allow LLMs to significantly contribute to my systems-engineering models and even component models 
- I feel like there's sometimes could be a better means of learning about these problems quicker, but this combination really is just too perfect to give up because it lets me get work done while at the same time understanding the problem better
- I think this is actually great news for me because it means that other people who primarily use AI without question will fall behind in this specific domain of deeply understanding problems, but that might just be cope on my end 
    - Also means that my future position is more defendable, but only time will tell if AI can automate this whole process, but I'm guessing not really because it's very useful having the person who's performing the tests also understand what kind of data they want to collect 
    - Also, by the time that AI can pretty much fully automate this, I'm sure other careers won't exactly be safe either

## Things I didn't like:

### Nothing so far!