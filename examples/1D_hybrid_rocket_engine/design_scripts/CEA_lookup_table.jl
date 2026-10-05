using CSV
using DataFrames
using Interpolations
using Unitful

# 1. Load the data
df = CSV.read(joinpath(@__DIR__, "cea_tables", "cea_table.csv"), DataFrame)

# 2. Extract unique grid points and sort them
pressures = sort(unique(df.P_chamber_pascals)) #Already in Pa
of_ratios = sort(unique(df.OF_ratio))
expansion_ratios = sort(unique(df.expansion_ratio))

# 3. Reshape the flat CSV columns into 3D arrays
grid_size = (length(pressures), length(of_ratios), length(expansion_ratios))
isp_grid = fill(NaN, grid_size)
cstar_grid = fill(NaN, grid_size)

pressure_indices = Dict(value => index for (index, value) in enumerate(pressures))
of_ratio_indices = Dict(value => index for (index, value) in enumerate(of_ratios))
expansion_ratio_indices = Dict(value => index for (index, value) in enumerate(expansion_ratios))

for row in eachrow(df)
    p_idx = pressure_indices[row.P_chamber_pascals]
    of_idx = of_ratio_indices[row.OF_ratio]
    expansion_ratio_idx = expansion_ratio_indices[row.expansion_ratio]
    
    isp_grid[p_idx, of_idx, expansion_ratio_idx] = row.Isp_s
    cstar_grid[p_idx, of_idx, expansion_ratio_idx] = row.Cstar_m_s
end

any(isnan, isp_grid) && error("The CEA Isp table is not a complete rectangular grid")
any(isnan, cstar_grid) && error("The CEA C-star table is not a complete rectangular grid")

# 4. Create the 3D linear interpolators
# Gridded(Linear()) is fast and ODE-friendly
isp_interpolator_Pa_unclamped = interpolate(
    (pressures, of_ratios, expansion_ratios),
    isp_grid,
    Gridded(Linear()),
)

function isp_interpolator_Pa(P, OF, expansion_ratio)
    P = clamp(P, minimum(pressures), maximum(pressures))
    OF = clamp(OF, minimum(of_ratios), maximum(of_ratios))
    expansion_ratio = clamp(
        expansion_ratio,
        minimum(expansion_ratios),
        maximum(expansion_ratios),
    )

    return isp_interpolator_Pa_unclamped(P, OF, expansion_ratio)
end

cstar_interpolator_Pa_unclamped = interpolate(
    (pressures, of_ratios, expansion_ratios),
    cstar_grid,
    Gridded(Linear()),
)

function cstar_interpolator_Pa(P, OF, expansion_ratio)
    P = clamp(P, minimum(pressures), maximum(pressures))
    OF = clamp(OF, minimum(of_ratios), maximum(of_ratios))
    expansion_ratio = clamp(
        expansion_ratio,
        minimum(expansion_ratios),
        maximum(expansion_ratios),
    )

    return cstar_interpolator_Pa_unclamped(P, OF, expansion_ratio)
end
