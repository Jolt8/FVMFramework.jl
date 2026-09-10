using ComponentArrays

function test_access(ca, ca_2, sym_access_list)
    for sym in sym_access_list
        #ca_at_sym = getproperty(ca, sym)
        getproperty(ca, sym) .+= 1.0 
        view(ca, sym) .+= 1.0
        getproperty(ca, sym) .+= getproperty(ca_2, sym)
        ca[sym] .+= ca_2[sym]
        #all four of the methods above cause a bunch of GC and runtime dispatch
    end
end

ca = ComponentVector(
    test = (
        haha = [1.0, 1.0],
        hehe = [1.0, 1.0],
    ),
    zoop = [1.0, 1.0]
)

ca_2 = deepcopy(ca)

using BenchmarkTools

#@btime test_access($(ca), $(ca_2))

sym_access_list = [:test, :zoop]

VSCodeServer.@profview @btime test_access($(ca), $(ca_2), $(sym_access_list))

function append_fixed_properties_to_u!(u, system)
    
end

function append_p_to_u!(u, p, system)
    
end

function append_fixed_properties_and_p_to_u!(u, p, system)
    fixed_properties = ComponentVector(system.fixed_properties_vec, system.fixed_properties_axes)

    append_fixed_properties_and_p_to_u!(u, system)

    append_p_to_u!(u, p, system)
end