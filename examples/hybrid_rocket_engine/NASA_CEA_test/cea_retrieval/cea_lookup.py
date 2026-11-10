# pip install rocketcea
#from rocketcea.cea_obj import fuelCards, add_new_fuel
from rocketcea.cea_obj import fuelCards, add_new_fuel
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.units import add_user_units
import csv
import numpy as np
import os

hdpe_card = """
fuel Polyethylene C 2.0 H 4.0 wt%=100.00 h,cal=-13860.0 t(k)=298.15 rho=0.95
"""

add_user_units('millipoise', 'Pa-s', 1e-4)
add_user_units('mcal/cm-K-s', 'W/m-K', 0.4184)

add_new_fuel('HDPE', hdpe_card)

# Initialize CEA with N2O and Polyethylene
cea = CEA_Obj(
    oxName='N2O', 
    fuelName='HDPE',
    pressure_units = 'Pa',
    temperature_units='K',
    density_units='kg/m^3',
    specific_heat_units='J/kg-K',
    enthalpy_units='J/kg',
    cstar_units='m/s',
    viscosity_units='Pa-s',
    thermal_cond_units='W/m-K'
)

pressures = np.linspace((100_000), (10_000_000), 50) # Sweep chamber pressures
of_ratios = np.linspace(1.0, 15.0, 100)     # Sweep O/F ratios
expansion_ratios = np.linspace(3.0, 10.0, 15) # Match the current design bounds in 0.5 increments

species_used = [
    "N2O",
    "C2H4",
    "N2",
    "O2",
    "CO2",
    "CO",
    "H2O",
    "H2"
]

alternative_species_used = [
    "*N2O",
    "*C2H4",
    "*N2",
    "*O2",
    "*CO2",
    "*CO",
    "*H2O",
    "*H2"
]

# Ensure the directory exists relative to this script
csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cea_tables", "cea_table.csv")
os.makedirs(os.path.dirname(csv_path), exist_ok=True)

with open(csv_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["P_chamber_pascals", "OF_ratio", "expansion_ratio", "chamber_cp", "chamber_mu", "chamber_k", "chamber_prandtl", "chamber_isp_vac", "chamber_cstar", "chamber_T_comb", "chamber_mw_avg", "chamber_gamma", "untracked_mass_fraction"] + species_used)
    
    for P in pressures:
        for OF in of_ratios:
            chamber_cp, chamber_mu, chamber_k, chamber_prandtl = cea.get_Chamber_Transport(Pc=P, MR=OF)

            molecular_weights, mass_fractions = cea.get_SpeciesMassFractions(
                Pc=P, MR=OF, min_fraction=0.0,
                frozen=0, frozenAtThroat=0
            )

            chamber_mass_fractions = {
                species: fractions[1]
                for species, fractions in mass_fractions.items()
            }

            final_species_mass_fractions = [
                chamber_mass_fractions.get(sp, chamber_mass_fractions.get(alt_sp, 0.0))
                for sp, alt_sp in zip(species_used, alternative_species_used)
            ]

            selected_mass_fraction_sum = sum(final_species_mass_fractions)
            if selected_mass_fraction_sum <= np.finfo(float).eps:
                raise ValueError(
                    "Selected CEA species have zero total mass fraction at "
                    f"P={P} Pa and O/F={OF}."
                )

            untracked_mass_fraction = 1.0 - selected_mass_fraction_sum
            final_species_mass_fractions = [
                mass_fraction / selected_mass_fraction_sum
                for mass_fraction in final_species_mass_fractions
            ]
            
            for eps in expansion_ratios:
                chamber_isp_vac, chamber_cstar, chamber_T_comb, chamber_mw_avg, chamber_gamma = cea.get_IvacCstrTc_ChmMwGam(Pc=P, MR=OF, eps=eps)

                writer.writerow([P, OF, eps, chamber_cp, chamber_mu, chamber_k, chamber_prandtl, chamber_isp_vac, chamber_cstar, chamber_T_comb, chamber_mw_avg, chamber_gamma, untracked_mass_fraction] + final_species_mass_fractions)
