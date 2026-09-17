Design goal: To create a program that automatically runs a given incompressible NavierStokes solver until it either completes or crashes

I think we should generally evaluate the stability of a problem by working through these branches and immediately terminating if one succeeds
1. First try to run a NonlinearProblem to see if the solver can converge to steady state
    If this succeeds, the problem is basically stable for all intents and purposes
2. If that fails, then run an implicit solve with SparspakFactorization() and see if it converges within the desired timeframe
3. Do the same with Krylov_GMRES()
    In general, it seems like Krylov_GMRES() is much more stable because it acts like a natural dampener compared to SparspakFactorization()



Evaluate the acoustic CFL limit for both convective and viscous terms (this should be the last resort, because it seems like sticking to this would make solving everything extremely slow)
- You find this by looking at the 



