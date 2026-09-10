
mutable struct VarAccessLog
    path::Vector{Symbol}
    accesses::Vector{Symbol}
    n_du_reads::Int
    n_du_writes::Int
    n_u_reads::Int
    n_u_writes::Int
end

function merge_trace_results(region_contexts)
    var_access_logs = Dict{Vector{Symbol}, VarAccessLog}()

    encountered_paths = Vector{Symbol}[]

    for var_access in region_contexts
        path = var_access.path

        access_type = var_access.access_type

        if !(path.path in encountered_paths) && path.root != :_expr
            push!(encountered_paths, path.path)
            push!(var_access_logs,
                path.path => VarAccessLog(
                    path.path,
                    Symbol[],
                    0, 0,
                    0, 0
                )
            )
        end

        if path.root == :du
            if access_type == :read
                var_access_logs[path.path].n_du_reads += 1
                push!(var_access_logs[path.path].accesses, :du_read)
            elseif access_type == :write
                var_access_logs[path.path].n_du_writes += 1
                push!(var_access_logs[path.path].accesses, :du_write)
            end
        elseif path.root == :u
            if access_type == :read
                var_access_logs[path.path].n_u_reads += 1
                push!(var_access_logs[path.path].accesses, :u_read)
            elseif access_type == :write
                var_access_logs[path.path].n_u_writes += 1
                push!(var_access_logs[path.path].accesses, :u_write)
            end
        end
    end

    return var_access_logs, encountered_paths
end

#var_access_logs, encountered_paths = merge_trace_results(ctx_du.access_logs)