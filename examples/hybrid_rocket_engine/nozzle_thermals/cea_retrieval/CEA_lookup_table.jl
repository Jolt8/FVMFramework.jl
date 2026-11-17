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
frozen_flags = sort(unique(df.frozen)) # 0 = shifting equilibrium, 1 = frozen

# 3. Reshape the flat CSV columns into 4D arrays (the last axis is the frozen flag)
grid_size = (length(pressures), length(of_ratios), length(expansion_ratios), length(frozen_flags))

chamber_cp_grid = fill(NaN, grid_size)
chamber_mu_grid = fill(NaN, grid_size)
chamber_k_grid = fill(NaN, grid_size)
chamber_prandtl_grid = fill(NaN, grid_size)
chamber_isp_vac_grid = fill(NaN, grid_size)
chamber_cstar_grid = fill(NaN, grid_size)
chamber_T_comb_grid = fill(NaN, grid_size)
chamber_mw_avg_grid = fill(NaN, grid_size)
chamber_gamma_grid = fill(NaN, grid_size)

N2O_grid = fill(NaN, grid_size)
C2H4_grid = fill(NaN, grid_size)
N2_grid = fill(NaN, grid_size)
O2_grid = fill(NaN, grid_size)
CO2_grid = fill(NaN, grid_size)
CO_grid = fill(NaN, grid_size)
H2O_grid = fill(NaN, grid_size)
H2_grid = fill(NaN, grid_size)

pressure_indices = Dict(value => index for (index, value) in enumerate(pressures))
of_ratio_indices = Dict(value => index for (index, value) in enumerate(of_ratios))
expansion_ratio_indices = Dict(value => index for (index, value) in enumerate(expansion_ratios))
frozen_indices = Dict(value => index for (index, value) in enumerate(frozen_flags))

for row in eachrow(df)
    p_idx = pressure_indices[row.P_chamber_pascals]
    of_idx = of_ratio_indices[row.OF_ratio]
    expansion_ratio_idx = expansion_ratio_indices[row.expansion_ratio]
    frozen_idx = frozen_indices[row.frozen]
    
    chamber_cp_grid[p_idx, of_idx, expansion_ratio_idx, frozen_idx] = row.chamber_cp
    chamber_mu_grid[p_idx, of_idx, expansion_ratio_idx, frozen_idx] = row.chamber_mu
    chamber_k_grid[p_idx, of_idx, expansion_ratio_idx, frozen_idx] = row.chamber_k
    chamber_prandtl_grid[p_idx, of_idx, expansion_ratio_idx, frozen_idx] = row.chamber_prandtl
    chamber_isp_vac_grid[p_idx, of_idx, expansion_ratio_idx, frozen_idx] = row.chamber_isp_vac
    chamber_cstar_grid[p_idx, of_idx, expansion_ratio_idx, frozen_idx] = row.chamber_cstar
    chamber_T_comb_grid[p_idx, of_idx, expansion_ratio_idx, frozen_idx] = row.chamber_T_comb
    chamber_mw_avg_grid[p_idx, of_idx, expansion_ratio_idx, frozen_idx] = row.chamber_mw_avg
    chamber_gamma_grid[p_idx, of_idx, expansion_ratio_idx, frozen_idx] = row.chamber_gamma

    N2O_grid[p_idx, of_idx, expansion_ratio_idx, frozen_idx] = row.N2O
    C2H4_grid[p_idx, of_idx, expansion_ratio_idx, frozen_idx] = row.C2H4
    N2_grid[p_idx, of_idx, expansion_ratio_idx, frozen_idx] = row.N2
    O2_grid[p_idx, of_idx, expansion_ratio_idx, frozen_idx] = row.O2
    CO2_grid[p_idx, of_idx, expansion_ratio_idx, frozen_idx] = row.CO2
    CO_grid[p_idx, of_idx, expansion_ratio_idx, frozen_idx] = row.CO
    H2O_grid[p_idx, of_idx, expansion_ratio_idx, frozen_idx] = row.H2O
    H2_grid[p_idx, of_idx, expansion_ratio_idx, frozen_idx] = row.H2
end

any(isnan, chamber_cp_grid) && error("The CEA chamber_cp table is not a complete rectangular grid")
any(isnan, chamber_mu_grid) && error("The CEA chamber_mu table is not a complete rectangular grid")
any(isnan, chamber_k_grid) && error("The CEA chamber_k table is not a complete rectangular grid")
any(isnan, chamber_prandtl_grid) && error("The CEA chamber_prandtl table is not a complete rectangular grid")
any(isnan, chamber_isp_vac_grid) && error("The CEA chamber_isp_vac table is not a complete rectangular grid")
any(isnan, chamber_cstar_grid) && error("The CEA chamber_cstar table is not a complete rectangular grid")
any(isnan, chamber_T_comb_grid) && error("The CEA chamber_T_comb table is not a complete rectangular grid")
any(isnan, chamber_mw_avg_grid) && error("The CEA chamber_mw_avg table is not a complete rectangular grid")
any(isnan, chamber_gamma_grid) && error("The CEA chamber_gamma table is not a complete rectangular grid")

any(isnan, N2O_grid) && error("The CEA N2O table is not a complete rectangular grid")
any(isnan, C2H4_grid) && error("The CEA C2H4 table is not a complete rectangular grid")
any(isnan, N2_grid) && error("The CEA N2 table is not a complete rectangular grid")
any(isnan, O2_grid) && error("The CEA O2 table is not a complete rectangular grid")
any(isnan, CO2_grid) && error("The CEA CO2 table is not a complete rectangular grid")
any(isnan, CO_grid) && error("The CEA CO table is not a complete rectangular grid")
any(isnan, H2O_grid) && error("The CEA H2O table is not a complete rectangular grid")
any(isnan, H2_grid) && error("The CEA H2 table is not a complete rectangular grid")

# 4. Create the 4D linear interpolators
# Gridded(Linear()) is fast and ODE-friendly
grid_axes = (pressures, of_ratios, expansion_ratios, frozen_flags)

# Wraps a 4D interp so it is still called as interp(P, OF, eps) and the frozen flag is a keyword
# frozen = false/true lands exactly on the 0/1 knot, so nothing is interpolated along that axis
with_frozen_flag(interp) = (P, OF, eps; frozen = false) -> interp(P, OF, eps, Int(frozen))

chamber_cp_interp = with_frozen_flag(interpolate(grid_axes, chamber_cp_grid, Gridded(Linear())))
chamber_mu_interp = with_frozen_flag(interpolate(grid_axes, chamber_mu_grid, Gridded(Linear())))
chamber_k_interp = with_frozen_flag(interpolate(grid_axes, chamber_k_grid, Gridded(Linear())))
chamber_prandtl_interp = with_frozen_flag(interpolate(grid_axes, chamber_prandtl_grid, Gridded(Linear())))
chamber_isp_vac_interp = with_frozen_flag(interpolate(grid_axes, chamber_isp_vac_grid, Gridded(Linear())))
chamber_cstar_interp = with_frozen_flag(interpolate(grid_axes, chamber_cstar_grid, Gridded(Linear())))
chamber_T_comb_interp = with_frozen_flag(interpolate(grid_axes, chamber_T_comb_grid, Gridded(Linear())))
chamber_mw_avg_interp = with_frozen_flag(interpolate(grid_axes, chamber_mw_avg_grid, Gridded(Linear())))
chamber_gamma_interp = with_frozen_flag(interpolate(grid_axes, chamber_gamma_grid, Gridded(Linear())))

N2O_interp = with_frozen_flag(interpolate(grid_axes, N2O_grid, Gridded(Linear())))
C2H4_interp = with_frozen_flag(interpolate(grid_axes, C2H4_grid, Gridded(Linear())))
N2_interp = with_frozen_flag(interpolate(grid_axes, N2_grid, Gridded(Linear())))
O2_interp = with_frozen_flag(interpolate(grid_axes, O2_grid, Gridded(Linear())))
CO2_interp = with_frozen_flag(interpolate(grid_axes, CO2_grid, Gridded(Linear())))
CO_interp = with_frozen_flag(interpolate(grid_axes, CO_grid, Gridded(Linear())))
H2O_interp = with_frozen_flag(interpolate(grid_axes, H2O_grid, Gridded(Linear())))
H2_interp = with_frozen_flag(interpolate(grid_axes, H2_grid, Gridded(Linear())))

mass_fractions_interps_vec = [
    N2O_interp,
    C2H4_interp,
    N2_interp,
    O2_interp,
    CO2_interp,
    CO_interp,
    H2O_interp,
    H2_interp,
]


