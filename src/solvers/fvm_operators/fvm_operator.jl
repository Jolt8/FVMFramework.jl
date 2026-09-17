
function fvm_operator!(du_vec, u_vec, p, t, system::FVMSystem, geo::FVMGeometry, solve_groups!)
    du, u = unpack_fvm_state(du_vec, u_vec, p, t, system)

    solve_groups!(du, u, p, t, system, geo)
end