
function classify_variables(var_access_logs)
    state_vars = Dict{Vector{Symbol}, VarAccessLog}()
    cache_vars = Dict{Vector{Symbol}, VarAccessLog}()
    fixed_vars = Dict{Vector{Symbol}, VarAccessLog}()

    for (path, log) in var_access_logs
        accesses = log.accesses

        last_du_access = nothing
        for i in length(accesses):-1:1
            if accesses[i] in (:du_read, :du_write)
                last_du_access = accesses[i]
                break
            end
        end

        if last_du_access == :du_write
            #the du of state variables are always written to last
            state_vars[path] = log
        elseif log.n_u_writes > 0 || log.n_du_writes > 0 || log.n_du_reads > 0
            #caches are never written to last, they're always read last
            cache_vars[path] = log
        else
            #u is the only thing read from fixed variables
            fixed_vars[path] = log
        end
    end

    return state_vars, cache_vars, fixed_vars
end

#state_vars, cache_vars, fixed_vars = classify_variables(var_access_logs)