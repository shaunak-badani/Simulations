"""
Run from root
python scripts/calculate_pressure.py
This script is incorrect as the product is taken for F_ija * r_ijb, not F_ia * r_ijb
"""

import MDAnalysis as mda

import warnings
# suppress some MDAnalysis warnings about PSF files
warnings.filterwarnings('ignore')
import numpy as np

from matplotlib import pyplot as plt

GRO = "./em/argon.gro"
TRR = "./nvt/argon_nvt.trr"


u = mda.Universe(GRO, TRR)
print(len(u.trajectory))
i = 0

# print("first few positions ", positions[:2])
kJ_mol_to_J = 1.6605391 * 1e-21
m_s_to_nm_ps = 1e-3
amu_to_kg = 1.6605391 * 1e-27
kJ_mol_nm_to_amu_nm_ps = kJ_mol_to_J * (m_s_to_nm_ps)**2 / (amu_to_kg)
# print("scaling factor  : ", kJ_mol_nm_to_amu_nm_ps)

positions = 0
velocities = 0
forces = 0


m = 1 # in amu
box = u.dimensions[:3] / 10 # box dimensions in nm now
volume = np.prod(box)
print("volume : ", volume)

nm_to_m = 1e-9
ps_to_s = 1e-12

amu_nm_ps_2_to_kg_m_s_2 = amu_to_kg / (nm_to_m * (ps_to_s)**2)
amu_nm_ps_2_to_bar = amu_nm_ps_2_to_kg_m_s_2 / 1e5
# print("pressure conversion factor : ", amu_nm_ps_2_to_bar)

def calculate_pressure(alpha, beta):
    assert alpha >= 0 and alpha < 3 and beta >= 0 and beta < 3, "invalid alpha and beta values!"
    first_term = m * np.sum(velocities[:, alpha] * velocities[:, beta])
    second_term = np.sum(forces[:, alpha] * positions[:, beta])
    total_term = first_term + second_term

    pressure_value = total_term / volume # this value is in amu / (ps^2 . nm)

    pressure_in_bars = pressure_value * amu_nm_ps_2_to_bar

    return pressure_in_bars

# N = len(u.trajectory)
N = 25
for i in range(0, N, 10):
    positions = (u.trajectory[i].positions / 10) # positions are now in nm
    velocities = (u.trajectory[i].velocities / 10) # velocities are now in nm / ps
    forces = u.trajectory[i].forces * 10 # forces are now in kJ / (mol . nm)
    forces *= kJ_mol_nm_to_amu_nm_ps
    print("pressure_value xx (in bars) : ", calculate_pressure(0, 0))