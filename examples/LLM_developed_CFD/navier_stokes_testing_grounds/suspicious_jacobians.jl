using SparseArrays
using LinearAlgebra
using ComponentArrays
using ForwardDiff
using Printf

"""
    check_phantom_dependencies(jac_sparsity, du, u, p, t, system, geo, f_closure_implicit;
                               phantom_threshold = 1e-15, max_display = 25, verbose = true)
    check_phantom_dependencies(jac_sparsity, jac_numeric, u, system;
                               phantom_threshold = 1e-15, max_display = 25, verbose = true)

Check for "phantom" dependencies in a sparse Jacobian: entries where `jac_sparsity`
declares a non-zero dependency (e.g. detected by SparseConnectivityTracer), but the actual
numerical Jacobian has a derivative of zero (or <= `phantom_threshold`).

Phantom dependencies bloat the sparsity pattern, increase the graph coloring number
(causing more ForwardDiff sweeps per time step), and add unnecessary fill-in to linear solvers.

Returns a NamedTuple containing diagnostics and summary metrics.
"""
function check_phantom_dependencies(
    jac_sparsity,
    u,
    system;
    phantom_threshold = 1e-15,
    max_display = 25,
    verbose = true
)
    jac_numeric = zeros(size(jac_sparsity)...)
    ForwardDiff.jacobian!(
        jac_numeric,
        u_in -> begin
            du_local = similar(u_in)
            f_closure_implicit(du_local, u_in, p_guess, 0.0)
            du_local
        end,
        u
    )
    sp_matrix = sparse(jac_sparsity)
    sp_rows, sp_cols, _ = findnz(sp_matrix)
    total_declared = length(sp_rows)

    comp_u = ComponentVector(u, system.state_axes)
    u_labels = ComponentArrays.labels(comp_u)

    phantom_count = 0
    exact_zero_count = 0
    near_zero_count = 0
    intracell_phantoms = 0
    intercell_phantoms = 0

    # Key: (row_var, col_var) => Dict with breakdown metrics
    coupling_data = Dict{Tuple{String, String}, Dict{Symbol, Any}}()

    for k in 1:total_declared
        r = sp_rows[k]
        c = sp_cols[k]
        val = jac_numeric[r, c]
        abs_val = abs(val)

        if abs_val <= phantom_threshold
            phantom_count += 1
            is_exact = (val == 0.0)
            if is_exact
                exact_zero_count += 1
            else
                near_zero_count += 1
            end

            # Extract variable name and cell index if present (e.g. "density[4]" -> "density", 4)
            r_match = match(r"^(.*)\[(\d+)\]$", u_labels[r])
            c_match = match(r"^(.*)\[(\d+)\]$", u_labels[c])

            r_var = r_match !== nothing ? r_match.captures[1] : replace(u_labels[r], r"\[\d+\]" => "")
            c_var = c_match !== nothing ? c_match.captures[1] : replace(u_labels[c], r"\[\d+\]" => "")
            r_cell = r_match !== nothing ? parse(Int, r_match.captures[2]) : 0
            c_cell = c_match !== nothing ? parse(Int, c_match.captures[2]) : 0

            is_intra = (r_cell != 0 && r_cell == c_cell)
            if is_intra
                intracell_phantoms += 1
            else
                intercell_phantoms += 1
            end

            key = (r_var, c_var)
            if !haskey(coupling_data, key)
                coupling_data[key] = Dict{Symbol, Any}(
                    :count => 0,
                    :exact_zeros => 0,
                    :near_zeros => 0,
                    :intra => 0,
                    :inter => 0,
                    :max_abs => 0.0,
                )
            end
            entry = coupling_data[key]
            entry[:count] += 1
            if is_exact
                entry[:exact_zeros] += 1
            else
                entry[:near_zeros] += 1
            end
            if is_intra
                entry[:intra] += 1
            else
                entry[:inter] += 1
            end
            entry[:max_abs] = max(entry[:max_abs], abs_val)
        end
    end

    # Check for missing dependencies: entries where numeric Jacobian is nonzero, but jac_sparsity is 0
    num_sp = sparse(jac_numeric)
    num_rows, num_cols, num_vals = findnz(num_sp)
    missing_couplings = 0
    missing_data = Dict{Tuple{String, String}, Int}()
    for k in eachindex(num_vals)
        r = num_rows[k]
        c = num_cols[k]
        if abs(num_vals[k]) > phantom_threshold && sp_matrix[r, c] == 0
            missing_couplings += 1
            r_var = replace(u_labels[r], r"\[\d+\]" => "")
            c_var = replace(u_labels[c], r"\[\d+\]" => "")
            key = (r_var, c_var)
            missing_data[key] = get(missing_data, key, 0) + 1
        end
    end

    phantom_pct = total_declared > 0 ? (phantom_count / total_declared) * 100.0 : 0.0
    active_count = total_declared - phantom_count

    if verbose
        println("="^78)
        println("                    PHANTOM DEPENDENCY AUDIT")
        println("="^78)
        @printf("Total Declared Sparsity Nonzeros: %d\n", total_declared)
        @printf("Active Couplings (|val| > %.1e): %d (%.1f%%)\n", phantom_threshold, active_count, 100.0 - phantom_pct)
        @printf("Phantom Couplings (|val| <= %.1e): %d (%.1f%%)\n", phantom_threshold, phantom_count, phantom_pct)
        @printf("  ├── Exact Zeros (val == 0.0)  : %d\n", exact_zero_count)
        @printf("  ├── Near Zeros (<= %.1e)     : %d\n", phantom_threshold, near_zero_count)
        @printf("  ├── Intra-Cell Phantoms       : %d\n", intracell_phantoms)
        @printf("  └── Inter-Cell Phantoms       : %d\n", intercell_phantoms)
        println("-"^78)

        if phantom_count > 0
            println("Top Phantom Couplings (Declared in jac_sparsity but zero numerically):")
            sorted_couplings = sort(collect(coupling_data), by=x->x[2][:count], rev=true)
            for (idx, ((row_var, col_var), data)) in enumerate(sorted_couplings)
                if idx > max_display
                    println("  ... and $(length(sorted_couplings) - max_display) more coupling types.")
                    break
                end
                @printf(
                    "  [%2d] %s depends on %s\n       Count: %d | Exact 0: %d | Intra: %d | Inter: %d | Max |val|: %.2e\n",
                    idx, row_var, col_var, data[:count], data[:exact_zeros], data[:intra], data[:inter], data[:max_abs]
                )
            end
        else
            println("No phantom dependencies detected. Sparsity pattern perfectly matches active derivatives.")
        end

        if missing_couplings > 0
            println("!"^78)
            @printf("CRITICAL WARNING: %d active derivatives are MISSING from jac_sparsity!\n", missing_couplings)
            for ((row_var, col_var), count) in sort(collect(missing_data), by=x->x[2], rev=true)
                println("  Missing: $row_var depends on $col_var ($count occurrences)")
            end
            println("!"^78)
        end
        println("="^78)
    end

    return (
        total_declared = total_declared,
        active_count = active_count,
        phantom_count = phantom_count,
        phantom_percentage = phantom_pct,
        exact_zeros = exact_zero_count,
        near_zeros = near_zero_count,
        intracell_phantoms = intracell_phantoms,
        intercell_phantoms = intercell_phantoms,
        phantom_couplings = coupling_data,
        missing_count = missing_couplings,
        missing_couplings = missing_data,
        jac_numeric = jac_numeric
    )
end

function check_phantom_dependencies(
    jac_sparsity,
    du,
    u,
    p,
    t,
    system,
    geo,
    f_closure_implicit;
    phantom_threshold = 1e-15,
    max_display = 25,
    verbose = true
)
    println("Evaluating exact numeric Jacobian at given state via ForwardDiff...")
    jac_numeric = zeros(size(jac_sparsity)...)
    ForwardDiff.jacobian!(
        jac_numeric,
        u_in -> begin
            du_local = similar(u_in)
            f_closure_implicit(du_local, u_in, p, t)
            du_local
        end,
        u
    )
    return check_phantom_dependencies(
        jac_sparsity,
        jac_numeric,
        u,
        system;
        phantom_threshold = phantom_threshold,
        max_display = max_display,
        verbose = verbose
    )
end

"""
    check_suspicious_jacobian(jac_sparsity, du, u, p, t, system, geo, f_closure_implicit;
                              tiny_threshold = 1e-12)

Evaluates numeric Jacobian and looks for non-zero couplings with magnitude below `tiny_threshold`.
Also accepts a precomputed `jac_numeric` matrix as the first argument.
"""
function check_suspicious_jacobian(
    jac_sparsity,
    du,
    u,
    p,
    t,
    system,
    geo,
    f_closure_implicit;
    tiny_threshold = 1e-12
)
    jac_numeric = zeros(size(jac_sparsity)...)
    ForwardDiff.jacobian!(
        jac_numeric, 
        u_in -> begin
            du_local = similar(u_in)
            f_closure_implicit(du_local, u_in, p, t)
            du_local
        end, 
        u
    )

    return check_suspicious_jacobian(jac_numeric, u, system; tiny_threshold = tiny_threshold)
end

function check_suspicious_jacobian(
    jac_numeric::AbstractMatrix,
    u,
    system;
    tiny_threshold = 1e-12
)
    comp_u = ComponentVector(u, system.state_axes)
    u_labels = ComponentArrays.labels(comp_u)
    
    rows, cols, vals = findnz(sparse(jac_numeric))

    suspicious_couplings = 0
    coupling_counts = Dict{Tuple{String, String}, Int}()

    for i in eachindex(vals)
        val = vals[i]
        if val != 0.0 && abs(val) < tiny_threshold
            suspicious_couplings += 1
            
            row_var = replace(u_labels[rows[i]], r"\[\d+\]" => "")
            col_var = replace(u_labels[cols[i]], r"\[\d+\]" => "")
            
            key = (row_var, col_var)
            coupling_counts[key] = get(coupling_counts, key, 0) + 1
        end
    end
    
    println("Found $suspicious_couplings suspiciously small non-zero couplings in the Jacobian (|val| < $tiny_threshold)!")
    println("Unique coupling types:")
    for ((row_var, col_var), count) in sort(collect(coupling_counts), by=x->x[2], rev=true)
        println("  $row_var depends on $col_var ($count occurrences)")
    end

    return coupling_counts
end