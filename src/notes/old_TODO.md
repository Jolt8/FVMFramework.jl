TODO:
    - I've just realized that since this solver is primarily going to be used for doing parameter fitting, making the actual internals of the solver parallelizable is actually pretty useless. Instead, I think it's going to be much better to distribute different guesses of the optimized parameters across a cluster and then let each core solve the problem independently and pick the best guess. 
        - This prevents bottlenecks due to the lack of communication needed between cores
        - It's great because it kinda limits me to doing smaller scale problems that are doable on a single core for debugging and then easily scaling them up with multiple parameters 
        - and I don't need to put any of the effort into making the solver parallel which would probably make everything uglier 
        - I'm pretty sure that EnsembleDistributed or EnsembleGPU could handle this for me and this should be pretty trivial to implement
        - However, one thing that would be tricky is that we wouldn't know if the solves that are taking the longest are the ones that are the most true to life or if they're just bad guesses. 
        - Perhaps we could check on the fly if lower solver times corresponds to higher accuracy and if it does, we could just terminate early for bad guesses, however this could get stuck in local minima. 
        - However, perhaps this wouldn't be necessary as higher stiffness is usually an indicator of bad guesses 
        - Basically, it would check if higher stiffness in the runs that finish very early creates a better fit to the data. If it does, we could let the next longest time to finish with similar parameters (very important*) to finish. As soon as we see that a higher stiffness problem as the one before produces parameters that are less accurate than the stiffness of the ones before, we stop all currently running simulations and then use that as the best guess.
            - *If a run with completely different parameters from another just happened to be slightly more or less stiff from another, 
            - Thus, we would probably have to check this along each row and column of possible parameters 
            - Furthermore, we would still want to check the parameters that enclose (are to the left and right) to make sure we don't miss the true minima between values.
            - Although, now that I'm thinking about this, how does this differ to just stopping rows/columns early that return worse fits?
            - **IMPORTANT** I've just realized that even taking stiffness into consideration is not necessary. Simply terminating runs that result in a worse fit is way better. 
            - I wonder if SciMLSensitivity.jl does this already?
    
    - For second-order spatial discretization, I've been recommended to use Venkatakrishnan limiters, Van Albada Limiters, or maybe even Learned Neural Flux Limiters. DONE!

    - We have to choose between implementing second-order spatial discretization, second-order meshes, moving meshes, and adaptive mesh refinement as the next project

    - Implement a method to check units right after finish_fvm_config 
        - This would also mean that unitful values would have to be stripped when finish_fvm_config is called
        - Note that you can turn any quantity into base SI units by doing upreferred(1.0u"cm")
        - You can also set your preferred units with Unitful.preferunits(u"s", u"m", u"kg", u"K")
        - This has kinda always bugged me, but there's no way to simplify units like in 1.0u"kg*m^2/s^2"|> u"J", but I guess that's for the best because it might convert something like Pa * s to Poise (who the fuck uses Poise)
        - Also, did you know that m_fix = Unitful.FixedUnits(u"m") exists where you can return errors when something like 1.0m_fix + 1.0u"cm" is called, that would probably be useful for checking unit consistency
    
    - Implement GPU acceleration later

    - Split internal physics, sources, and capacity parameters instead of all of them being lumped in one big region
        - This is a maybe, it would be more modular but would require the redefnition of properties and would not be great for connections
        - I'm going to ignore this one, but it might be implemented later if there's a need for it

    - Implement FEM eventually for things like stress analysis (long way away from this one)

    - Implement a parameter optimization pipeline for fitting (method 1) and finding (method 2) best parameters (probably two separate files/folders/methods)
        - Partially done, I don't think it's going to be worth it to try to create methods that try to cover the many different ways optimization is done
        - instead, I think it's much better we just show how users can do something and then let them implement their own methods on a case-by-case basis. 

    - Implement VLE (Vapor-Liquid-Equilibrium)

    - Implement Navier-Stokes

    - Create methods for interfacing with other FVM or DG solvers like Trixi.jl or Oceanigans.jl

    - Create methods to more easily diagnose problems through logging. Graphing a certain variable over time would be very useful
        - this really needs to be worked on 
        - Just some brainstorming:
            - It would be nice if every single variable was tracked
            - This could probably be achieevd through a debug version of each flux and internal physics function where if problems emerge such as values like NaN, Inf, -Inf, 0.0, etc.
        - We should probably look into the package JuliaDebugger that provides breakpoints that would be insanely useful

    - Create a method for doing distributing computing by recruiting my home PC as a worker to prevent having to do simulation on my laptop in class

    - Implement MPI