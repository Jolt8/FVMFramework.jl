Conclusions from my initial testing of neural networks through Lux and DiffEqFlux and symbolic regression through DataDrivenDiffEq and SymbolicRegression.jl

Neural networks with Lux and DiffEqFlux:
    - From some initial testing, these have worked incredibly well and are actually really performant for meshes with < ~ 3000 cells (could probably push it to 30000 cells if you were willing to wait long enough)
    - Getting a loss pretty close to the actual solution on a pretty terrible laptop only took like 10 minutes on a 1100 cells grid to predict the heat transfer equation
    - this equation:
    '''julia
    function heat_diffusion!(
        du, u, 
        idx_a, idx_b, face_idx,
        area, norm, dist,
        vol_a, vol_b
    )
        k_effective = 2 * u.k[idx_a] * u.k[idx_b] / (u.k[idx_a] + u.k[idx_b])

        grad_T = (u.temp[idx_b] - u.temp[idx_a]) / dist

        du.heat[idx_a] -= -k_effective * grad_T * area
    end
    '''
    - I have not tried running the neural network on a GPU yet, and it was clear during the training that the slowness caused by running the neural network was the only thing holding things back
    - furthermore, until Enzyme.jl gets updated to julia 1.12, we're forced to use forward mode automatic differentiation which was pretty slow to find the solution 
    - Once Enzyme is fixed, I'm sure that the optimization of the model itself won't take that long
    - CONCLUSION: neural networks will definitely be useful in the future for extremely unknown physics 

Symbolic Regression with DataDrivenDiffEq and SymbolicRegression.jl:
    - I've tried to get this to work though multiple tests, and it has been largely terrible 
    - I think the biggest issue for symbolic regression is kinda a "chicken and egg" problem where most of the fundamental equations used have already been discovered and the equations that have not been discovered yet are likely too complex to use symbolic regression to find
    - In these situations where an already existing equation does not exist, the best solution is likely just to use a neural network to predict the missing equation(s) since they do a great job of minimizing the loss for very unknown physics like the pressure drop over a completely new and unique geometry
    - It's a little less ideal for research papers because you can't exactly write out a neural layer, but for a non-academic settings, it probably doesn't matter that you didn't spend the months to find an equation that matches you data
    - Even for finding equations like empirical correlations for the performance of different heat exchanger types, the best solution would likely be to take and already existing equation structure and then use paramater optimization (very fast) to optimize for the parameters that cause the model to match the regime you're working in
    - CONCLUSION: symbolic regression is not going to be useful in the future and I'll likely avoid it 