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

#value=[massfrac_injface, massfrac_chm, massfrac_tht, massfrac_exit]
molecular_weights, mass_fractions = cea.get_SpeciesMassFractions(
    Pc=2_000_000, MR=5.0, eps=10.0, min_fraction=0,
    frozen=0, frozenAtThroat=0
)

chamber_mass_fractions = dict(
    sorted(
        ((species, fractions[1]) for species, fractions in mass_fractions.items()),
        key=lambda item: item[1],
        reverse=True,
    )
)

print(len(chamber_mass_fractions))
print(chamber_mass_fractions)