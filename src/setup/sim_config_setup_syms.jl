function add_setup_syms!(
    config;
    cache_syms_and_units,
    special_caches,
    second_order_syms = [],
    optimized_parameters
)
    config.cache_syms_and_units = cache_syms_and_units
    config.special_caches = special_caches
    config.second_order_syms = second_order_syms
    config.optimized_parameters = optimized_parameters

    for cache_sym in propertynames(config.special_caches)
        config.cache_syms_and_units = (; config.cache_syms_and_units..., cache_sym => unit(first(getproperty(config.special_caches, cache_sym))))
    end

    return 
end