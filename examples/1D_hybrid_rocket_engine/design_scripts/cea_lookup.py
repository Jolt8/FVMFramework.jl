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
eps = 10.0                                 # Constant expansion ratio

# Ensure the directory exists relative to this script
csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cea_tables", "cea_table.csv")
os.makedirs(os.path.dirname(csv_path), exist_ok=True)

with open(csv_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["P_chamber_pascals", "OF_ratio", "Isp_s", "Cstar_m_s"])
    
    for P in pressures_psia:
        for OF in of_ratios:
            # Get Isp and Cstar from CEA
            isp_obj = cea.get_Isp(Pc=P, MR=OF, eps=eps)
            cstar_ft_s = cea.get_Cstar(Pc=P, MR=OF)
            
            # Convert Cstar to m/s
            cstar_m_s = cstar_ft_s * 0.3048
            writer.writerow([P * 6894.76, OF, isp_obj, cstar_m_s])