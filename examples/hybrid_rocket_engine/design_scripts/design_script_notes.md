# Model 1: Design Script Notes

## Things I liked

### I was able to maintain complete ownership of this model and I learned a lot about the problem while doing so
- Self explanatory, since this model is pretty simple all things considered, I did not have to use an LLM that much so I retained complete ownership of the script 

### Comparing to historical data 
- Using historical data gathered using NotebookLM deep research, I was able to compare the results that my solver came up with to historical data and they matched pretty closely after a couple tweaks, so I really trusted the model after that and the parameters that it gave me.

### Data gathered by this model
- This model pretty much became the source of truth used for all other models so this was probably the most important model of this entire project

## Problems and Improvements

### Passing dual numbers into the solve was a pain in the ass, forcing me to use an explicit ODE Solver
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


## Things left to decide:

### I should've probably provided a better means of interfacing with the results of this solver rather than just purely referencing the values stored in my head
Notes:
- This one is slightly debatable
- On one hand, being able to access this solver directly would make sure that what we're designign for would always stay consistent with this source of truth.
- On the other hand, this project's needs, requirements, design specs, etc. have slowly morphed over time into something slightly different
    - Unsure if the different is better or worse
    - Also, using values stored in my head is much quicker
    - That being said, it would not be hard at all to import a common ComponentVector storing all the design parameters for this problem similar to the way I've used the freecad table data with the cad_data ComponentVector