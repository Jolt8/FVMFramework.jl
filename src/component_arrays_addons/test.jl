using ComponentArrays
using FVMFramework
using BenchmarkTools
using Polyester

n_cells = 1_000_000

test_component_array = ComponentVector(
    mass_fractions = (
        oxygen = ones(n_cells).*0.5,
        nitrogen = ones(n_cells).*0.5,
    ),
    elemental_compositions = (
        oxygen = (
            O = ones(n_cells).*2.0,
            N = ones(n_cells).*0.0,
        ),
        nitrogen = (
            O = ones(n_cells).*0.0,
            N = ones(n_cells).*2.0,
        ),
    )
)

for species in Val.(propertynames(test_component_array.mass_fractions))
    getproperty(test_component_array.mass_fractions, species)
end

system = (species = Val.(propertynames(test_component_array.mass_fractions)), elements = Val.(propertynames(test_component_array.elemental_compositions.oxygen)))

function check_speed_with_system(u, system)
    for cell_id in eachindex(u.mass_fractions.oxygen)
        for species in system.species
            #view(u.mass_fractions, species)[cell_id] += 1.0
            getproperty(u.mass_fractions, species)[cell_id] += 1.0

            for element in system.elements
                getproperty(getproperty(u.elemental_compositions, species), element)[cell_id] += 1.0
            end
        end
    end
end

function check_speed_with_system_modified(u, system)
    for species in system.species
        u_mf = getproperty(u.mass_fractions, species)
        @inbounds @simd for cell_id in eachindex(u_mf)
            u_mf[cell_id] += 1.0
        end
        for element in system.elements
            u_elem = getproperty(getproperty(u.elemental_compositions, species), element)
            @inbounds @simd for cell_id in eachindex(u_elem)
                u_elem[cell_id] += 1.0
            end
        end
    end
end

function check_speed_with_forfields(u, system)
    for cell_id in eachindex(u.mass_fractions.oxygen)
        for_fields!(u.mass_fractions) do species, u_mass_fractions
            u_mass_fractions[species[cell_id]] += 1.0
        end
    end
end

function check_speed_with_system_batched(u, system)
    @batch for cell_id in eachindex(u.mass_fractions.oxygen)
        for species in system.species
            #view(u.mass_fractions, species)[cell_id] += 1.0
            getproperty(u.mass_fractions, species)[cell_id] += 1.0

            for element in system.elements
                getproperty(getproperty(u.elemental_compositions, species), element)[cell_id] += 1.0
            end
        end
    end
end

function check_speed_with_forfields_batched(u, system)
    @batch for cell_id in eachindex(u.mass_fractions.oxygen)
        for_fields!(u.mass_fractions) do species, u_mass_fractions
            u_mass_fractions[species[cell_id]] += 1.0
        end
    end
end

check_speed_with_forfields(test_component_array, system)

@btime check_speed_with_system($test_component_array, $system)
@btime check_speed_with_forfields($test_component_array, $system)

@btime check_speed_with_system_modified($test_component_array, $system)

@btime check_speed_with_system_batched($test_component_array, $system)
@btime check_speed_with_forfields_batched($test_component_array, $system)
