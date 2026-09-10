function unpack_fvm_state(du_vec, u_vec, p, t, system)
    du_vec .= 0.0 #zero out all derivatives

    #if you're wondering why we zero out everything, it's because for some reason calling get_tmp() returns all undef values, so we have to zero it out ourselves. 
    #I think this is intended behaviour
    
    if (first(u_vec) + first(p)) isa SparseConnectivityTracer.Dual{Float64}
        du = VirtualFVMArray((du_vec, (get_tmp(system.du_diff_cache, first(u_vec) + first(p)) .= 0.0)), system.du_virtual_axes)
        u = VirtualFVMArray((u_vec, (get_tmp(system.u_diff_cache, first(u_vec) + first(p)) .= 0.0), system.properties_vec, p), system.u_virtual_axes)
    else
        du = VirtualFVMArray((du_vec, (get_tmp(system.du_diff_cache, first(u_vec) + first(p)) .= 0.0)), system.du_virtual_axes)
        u = VirtualFVMArray((u_vec, (get_tmp(system.u_diff_cache, first(u_vec) + first(p)) .= 0.0), system.properties_vec, p), system.u_virtual_axes)
    end
    return du, u
end