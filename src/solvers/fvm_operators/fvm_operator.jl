
function fvm_operator!(du_vec, u_vec, p, t, solve_groups!, geo::FVMGeometry, system::FVMSystem)
    du, u = unpack_fvm_state(du_vec, u_vec, p, t, system)

    solve_groups!(du, u, p, t, geo, system)
end