"""
    algebraicmultigrid(W, du, u, p, t, newW, Plprev, Prprev, solverdata)

This function takes in the Jacobian of the FVM operator, and the solution vector, and returns the Algebraic Multigrid preconditioner. 

Parameters
    - W: The Jacobian of the FVM operator. 
    - du: The solution vector. 
    - u: The solution vector. 
    - p: The parameters. 
    - t: The time. 
    - newW: Whether the Jacobian needs to be updated. 
    - Plprev: The previous preconditioner. 
    - Prprev: The previous preconditioner. 
    - solverdata: The solver data. 

Returns
    - Pl: The preconditioner. 
    - nothing
"""
function algebraicmultigrid(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
    if newW === nothing || newW
        Pl = AlgebraicMultigrid.aspreconditioner(AlgebraicMultigrid.ruge_stuben(convert(AbstractMatrix, W)))
    else
        Pl = Plprev
    end
    Pl, nothing
end

"""
    iluzero(W, du, u, p, t, newW, Plprev, Prprev, solverdata)

This function takes in the Jacobian of the FVM operator, and the solution vector, and returns the ILU(0) preconditioner. 

Parameters
    - W: The Jacobian of the FVM operator. 
    - du: The solution vector. 
    - u: The solution vector. 
    - p: The parameters. 
    - t: The time. 
    - newW: Whether the Jacobian needs to be updated. 
    - Plprev: The previous preconditioner. 
    - Prprev: The previous preconditioner. 
    - solverdata: The solver data. 

Returns
    - Pl: The preconditioner. 
    - nothing
"""
function iluzero(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
    if newW === nothing || newW
        Pl = ilu0(convert(AbstractMatrix, W))
    else
        Pl = Plprev
    end
    Pl, nothing
end