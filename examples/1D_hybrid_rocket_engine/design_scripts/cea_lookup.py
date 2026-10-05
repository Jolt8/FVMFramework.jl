# pip install rocketcea
from rocketcea.cea_obj import CEA_Obj, fuelCards, add_new_fuel
import csv
import numpy as np
import os

hdpe_card = """
fuel Polyethylene C 2.0 H 4.0 wt%=100.00 h,cal=-13860.0 t(k)=298.15 rho=0.95
"""

add_new_fuel('HDPE', hdpe_card)

# Initialize CEA with N2O and Polyethylene
cea = CEA_Obj(oxName='N2O', fuelName='HDPE')

pressures_psia = np.linspace((100_000 / 6894.76), (10_000_000 / 6894.76), 50) # Sweep chamber pressures
of_ratios = np.linspace(1.0, 15.0, 100)     # Sweep O/F ratios
expansion_ratios = np.linspace(3.0, 10.0, 15) # Match the current design bounds in 0.5 increments

# Ensure the directory exists relative to this script
csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cea_tables", "cea_table.csv")
os.makedirs(os.path.dirname(csv_path), exist_ok=True)

with open(csv_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["P_chamber_pascals", "OF_ratio", "expansion_ratio", "Isp_s", "Cstar_m_s"])
    
    for P in pressures_psia:
        for OF in of_ratios:
            # C-star is a chamber property and does not depend on nozzle expansion
            # ratio, so calculate it once and repeat it across the epsilon grid.
            cstar_ft_s = cea.get_Cstar(Pc=P, MR=OF)
            cstar_m_s = cstar_ft_s * 0.3048

            for eps in expansion_ratios:
                isp_obj = cea.get_Isp(Pc=P, MR=OF, eps=eps)
                writer.writerow([P * 6894.76, OF, eps, isp_obj, cstar_m_s])
