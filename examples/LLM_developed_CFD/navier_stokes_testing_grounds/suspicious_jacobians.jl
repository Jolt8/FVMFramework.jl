using SparseArrays
using ComponentArrays

function check_suspicious_jacobian(jac_sparsity, du, u, p, t, system, geo, f_closure_implicit)
    jac_numeric = zeros(size(jac_sparsity)...)
    # 2. Evaluate the exact numeric Jacobian at the initial state
    # (This might take a moment, but it computes the actual derivative values)
    ForwardDiff.jacobian!(
        jac_numeric, 
        u -> (du = similar(u); f_closure_implicit(du, u, p_guess, 0.0); du), 
        u0_vec
    )
    # 3. Hunt for suspiciously tiny couplings (e.g., between 0 and 1e-12)
    tiny_threshold = 1e-12

    comp_u = ComponentVector(u, system.state_axes)
    u_labels = ComponentArrays.labels(comp_u)
    
    rows, cols, vals = findnz(sparse(jac_numeric))

    suspicious_couplings = 0
    coupling_counts = Dict{Tuple{String, String}, Int}()

    for i in eachindex(vals)
        val = vals[i]
        if val != 0.0 && abs(val) < tiny_threshold
            suspicious_couplings += 1
            
            # Strip the "[cell_index]" from the labels to group them
            row_var = replace(u_labels[rows[i]], r"\[\d+\]" => "")
            col_var = replace(u_labels[cols[i]], r"\[\d+\]" => "")
            
            key = (row_var, col_var)
            coupling_counts[key] = get(coupling_counts, key, 0) + 1
        end
    end
    
    println("Found $suspicious_couplings suspiciously small couplings in the Jacobian!")
    println("Unique coupling types:")
    for ((row_var, col_var), count) in sort(collect(coupling_counts), by=x->x[2], rev=true)
        println("  $row_var depends on $col_var ($count occurrences)")
    end
end