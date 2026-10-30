object_names = [
    "injector_bulkhead_top",
    "injector_bulkhead_bottom",
    "mid_section",
    "pre_chamber",
    "phenolic_liner_top",
    "aluminum_housing_top",

    #(1-n cells)
    "fluid_domain",
    "fuel_grain",
    "phenolic_liner",
    "aluminum_housing",
    #(1-n cells) end

    "recirculation_zone", #abstraction; not a physical volume

    "nozzle_top",
    "nozzle_bottom",
    "aluminum_housing_nozzle",
    "aluminum_housing_bottom",
    "nozzle_retainer"
]

room_temp = 294.13u"K"

#we would then add all the region functions (this part is very simple)

#we're still going to need a script in the middle here that gets the geometry, data strucutres, and topology we need to complete this
#this is going to be very difficult!

#let's do an example:
#this one would be applied to every single fluid_domain cell

fluid_domain_volume = geo.cell_volumes[fluid_domain[1]]

air_initial_conditions, air_properties, air_cache_syms_and_units, air_special_caches = get_air_properties(fluid_domain_volume, room_temp)

function overall_combustion_model!(du, u, p, t, system, geo, cell_id)
    decompose_nitrous_oxide!(du, u, p, t, system, geo, cell_id)
    combust_ethylene_and_oxygen_westbrook_dryer!(du, u, p, t, system, geo, cell_id) 
    #combust_ethylene_and_oxygen_simple!(du, u, p, t, system, geo, cell_id) #is also an option
end

function update_fluid_properties!(du, u, p, t, system, geo, cell_id)
    overall_fluid_property_update!(du, u, p, t, system, geo, cell_id)
    update_species_molar_concentrations!(du, u, p, t, system, geo, cell_id)
    update_fluid_diffusion_coefficients!(du, u, p, t, system, geo, cell_id)
end

function fluid_domain_physics!(du, u, p, t, system, geo, cell_id)
    overall_combustion_model!(du, u, p, t, system, geo, cell_id)
    cap_fluid_flow!(du, u, p, t, system, geo, cell_id)
end

add_region!(
    config, "fluid_domain";
    type = Solid(), #Not sure if we even need a type, my gut feeling says no
    initial_conditions = air_initial_conditions,
    properties = air_properties,
    cache_syms_and_units = air_cache_syms_and_units,
    special_caches = air_special_caches,
    property_update_function = 
    function update_fluid_domain_region!(du, u, p, t, system, geo, cell_id)
        update_fluid_properties!(du, u, p, t, system, geo, cell_id)
    end,
    region_function = 
    function fluid_domain_region!(du, u, p, t, system, geo, cell_id)
        fluid_domain_physics!(du, u, p, t, system, geo, cell_id)
    end
)

add_patch!(
    owners = ["mid_section"],
    connections = ["injector_bulkhead_top", "injector_bulkhead_bottom"],
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo,
        mid_section, injector_bulkhead_top, injector_bulkhead_bottom
    )
        mid_section_to_injector_bulkhead_empirical_heat_flux!(
            mid_section,
            injector_bulkhead_top,
            injector_bulkhead_top,
            geo
        )

        #this would do something like:
        #= 
        
        du.heat[mid_section] += (u.temperature[connection_idxs.injector_bulkhead_top] - u.temperature[owner_idxs.mid_section]) * 
            geo.interface_areas[mid_section, injector_bulkhead_top] *
            geo.interface_distances[mid_section, injector_bulkhead_top]

        etc. for the heats of other regions
        =#
    end
)

add_patch!(
    owner = "injector_bulkhead_top",
    connection = "injector_bulkhead_bottom",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo,
        injector_bulkhead_top, injector_bulkhead_bottom #I'm thinking that these could either store an Int or 1:Int
    )
        heat_diffusion!(du, u, p, t, system, geo, [injector_bulkhead_top], [injector_bulkhead_bottom])
    end
)

add_patch!(
    owner = "mid_section",
    connection = "pre_chamber",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo,
        mid_section, pre_chamber
    )
        mid_section_to_pre_chamber_flux!(du, u, p, t, system, geo, mid_section, pre_chamber) 
        #this would be the valve equation determining flow between the mid section and first chamber cell
        #purely algebraic
    end
)

#Axial fluxes from top cells to injector_bulkhead_bottom:
#start
add_patch!(
    owner = "pre_chamber",
    connection = "injector_bulkhead_bottom"
)

add_patch!(
    owner = "phenolic_liner_top",
    connection = "injector_bulkhead_bottom"
)

add_patch!(
    owner = "aluminum_housing_top",
    connection = "injector_bulkhead_bottom"
)
#end



#Radial fluxes across top cells
#start
add_patch!(
    owner = "pre_chamber",
    connection = "phenolic_liner_top"
)

add_patch!(
    owner = "phenolic_liner_top",
    connection = "aluminum_housing_top"
)
#end


add_patch!(
    owner = "pre_chamber",
    connection = "fluid_domain", #1-n cells,
)

#We could add a flux between the top fluid domain and the fuel grain, but I'm going to ignore it

#Axial fluxes along n cells
#start
add_patch!(
    owner = "fluid_domain",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, cells
    )
        for (upstream_cell, downstream_cell) in upstream_downstream(cells)

            fuel_port_id = system.fuel_port_ids[upstream_cell, downstream_cell]

            interface_area = pi * u.fuel_port_radius[fuel_port_id]^2

            axial_fluid_domain_flux!(du, u, p, t, system, interface_area, upstream_cell, downstream_cell)
        end
    end
)

add_patch!(
    owner = "fuel_grain",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, cells
    )
        for (upstream_cell, downstream_cell) in upstream_downstream(cells)
            fuel_port_id = system.fuel_port_ids[upstream_cell, downstream_cell]

            fluid_port_area = pi * u.fuel_port_radius[fuel_port_id]^2
            fuel_grain_area = pi * (geo.max_bounding_box_x[upstream_cell] - geo.centroids[upstream_cell][1])
            interface_area = fuel_grain_area - fluid_port_area

            custom_geo_heat_diffusion!(
                du, u, p, t, system, 
                interface_area, #we would also need cell centroids, cell volumes, and other geo information
                upstream_cell, downstream_cell
            )
        end
    end
)

add_patch!(
    owner = "phenolic_liner",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, cells
    )
        for (upstream_cell, downstream_cell) in upstream_downstream(cells)
            heat_diffusion!(du, u, p, t, system, geo, upstream_cell, downstream_cell)
            #we can use their default geometries because they don't change
        end
    end
)

add_patch!(
    owner = "aluminum_housing",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, cells
    )
        for (upstream_cell, downstream_cell) in upstream_downstream(cells)
            heat_diffusion!(du, u, p, t, system, geo, upstream_cell, downstream_cell)
            #we can use their default geometries because they don't change
        end
    end
)
#end

#Radial fluxes along cells
#start
add_patch!(
    owner = "fluid_domain",
    connection = "fuel_grain",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, fluid_domain_cells, fuel_grain_cells
    )
        for fluid_cell in fluid_domain_cells
            fuel_grain_cell = system.fluid_fuel_grain_neighbors[fluid_cell]

            #we get custom geometry here

            #do fluxes
        end
    end
)

add_patch!(
    owner = "fuel_grain",
    connection = "phenolic_liner",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, fuel_grain_cells, phenolic_liner_cells
    )
        for fuel_grain_cell in fuel_grain_cells
            phenolic_liner_cell = system.fuel_grain_phenolic_liner_neighbors[fluid_cell]

            #we get semi-custom geoemtry here
            #face areas between these two don't change, but the fuel grain centroid gets closer

            port_id = system.port_id_based_on_fuel_cell[fuel_grain_cell]

            fuel_grain_centroid = (geo.max_bounding_boxes_x[fuel_grain_cell] - u.fuel_port_radius[port_id])

            custom_geo_heat_diffusion!(
                du, u, p, t, system, 
                fuel_grain_centroid, #we would also need cell centroids, cell volumes, and other geo information
                upstream_cell, downstream_cell
            )

            #do fluxes
        end
    end
)

add_patch!(
    owner = "phenolic_liner",
    connection = "aluminum_housing",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, neighbors
    )
        for (idx_a, idx_b) in neighbors
            heat_diffusion!(du, u, p, t, system, geo, idx_a, idx_b)
        end
    end
)
#end

add_patch!(
    owner = "fluid_domain",
    connection = "nozzle_top",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, fluid_domain, nozzle_top
    )
        nozzle_outlet!(du, u, p, t, system, geo, fluid_domain[end], nozzle_top)
        heat_transfer_from_fluid_to_nozzle_top!(du, u, p, t, system, geo, fluid_domain[end], nozzle_top)
    end
)

add_patch!(
    owner = "fluid_domain",
    connection = "recirculation_zone",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, fluid_domain, recirculation_zone
    )
        recirculation_zone!(du, u, p, t, system, geo, fluid_domain[end], recirculation_zone)
    end
)

#we could also add a flux between the nozzle top and phenolic liner but I don't think it's necessary


#This is a radial flux between the nozzle graphite and aluminum housing
add_patch!(
    owner = "nozzle_top",
    connection = "aluminum_housing_nozzle",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, nozzle_top, recirculation_zone
    )
        nozzle_top_to_aluminum_housing_nozzle_flux!(du, u, p, t, system, geo, nozzle_top, recirculation_zone)
    end
)

add_patch!(
    owner = "nozzle_top",
    connection = "nozzle_bottom",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, nozzle_top, nozzle_bottom
    )
        heat_diffusion!(du, u, p, t, system, geo, nozzle_top, nozzle_bottom)
    end
)

add_patch!(
    owner = "nozzle_bottom",
    connection = "aluminum_housing_bottom",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, nozzle_bottom, aluminum_housing_bottom
    )
        heat_diffusion!(du, u, p, t, system, geo, nozzle_bottom, aluminum_housing_bottom)
    end
)

add_patch!(
    owner = "nozzle_bottom",
    connection = "nozzle_retainer",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, nozzle_bottom, nozzle_retainer
    )
        heat_diffusion!(du, u, p, t, system, geo, nozzle_bottom, nozzle_retainer)
    end
)

add_patch!(
    owner = "aluminum_housing_bottom",
    connection = "nozzle_retainer",
    patch_function = function fuel_wall_patch_flux!(
        du, u, p, t, system, geo, aluminum_housing_bottom, nozzle_retainer
    )
        heat_diffusion!(du, u, p, t, system, geo, aluminum_housing_bottom, nozzle_retainer)
    end
)