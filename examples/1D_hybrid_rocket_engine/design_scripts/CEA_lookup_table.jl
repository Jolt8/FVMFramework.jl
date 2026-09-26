using CSV
using DataFrames
using Interpolations
using Unitful

# 1. Load the data
df = CSV.read(joinpath(@__DIR__, "cea_tables", "cea_table.csv"), DataFrame)

# 2. Extract unique grid points and sort them
pressures = sort(unique(df.P_chamber_pascals)) #Already in Pa
of_ratios = sort(unique(df.OF_ratio))

# 3. Reshape the flat CSV columns into 2D matrices
isp_grid = zeros(length(pressures), length(of_ratios))
cstar_grid = zeros(length(pressures), length(of_ratios))

for row in eachrow(df)
    p_idx = findfirst(x -> x == row.P_chamber_pascals, pressures)
    of_idx = findfirst(x -> x == row.OF_ratio, of_ratios)
    
    isp_grid[p_idx, of_idx] = row.Isp_s
    cstar_grid[p_idx, of_idx] = row.Cstar_m_s
end

# 4. Create the 2D linear interpolators
# Gridded(Linear()) is fast and ODE-friendly
isp_interpolator_Pa_unclamped = interpolate((pressures, of_ratios), isp_grid, Gridded(Linear()))

function isp_interpolator_Pa(P, OF)
    P = clamp(P, minimum(pressures), maximum(pressures))
    OF = clamp(OF, minimum(of_ratios), maximum(of_ratios))

    #=
    @show "ISP"
    @show P
    @show OF

    @show cstar_interpolator_Pa_unclamped(P, OF)
    =#

    return isp_interpolator_Pa_unclamped(P, OF)
end

cstar_interpolator_Pa_unclamped = interpolate((pressures, of_ratios), cstar_grid, Gridded(Linear()))

function cstar_interpolator_Pa(P, OF)
    P = clamp(P, minimum(pressures), maximum(pressures))
    OF = clamp(OF, minimum(of_ratios), maximum(of_ratios))

    #=
    @show "C_STAR"

    @show P
    @show OF

    @show cstar_interpolator_Pa_unclamped(P, OF)
    =#

    return cstar_interpolator_Pa_unclamped(P, OF)
end

#=
# Create an easy wrapper function to call
function get_cea_properties(P_chamber_pa, OF_ratio)
    # Be sure to handle out-of-bounds queries so your ODE solver doesn't crash!
    # (e.g., clamp the OF ratio between your min and max values)
    P_psia = P_chamber_pa * 0.000145038 
    
    clamped_P = clamp(P_psia, minimum(pressures), maximum(pressures))
    clamped_OF = clamp(OF_ratio, minimum(of_ratios), maximum(of_ratios))
    
    return isp_interpolator(clamped_P, clamped_OF), cstar_interpolator(clamped_P, clamped_OF)
end
=#
