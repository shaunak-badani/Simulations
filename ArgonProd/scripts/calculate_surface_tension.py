import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# scaling factor
m_s_to_nm_ps = 1e-3
amu_to_kg = 1.6605391 * 1e-27
nm_to_m = 1e-9
ps_to_s = 1e-12

amu_nm_ps_2_to_kg_m_s_2 = amu_to_kg / (nm_to_m * (ps_to_s)**2)
amu_nm_ps_2_to_bar = amu_nm_ps_2_to_kg_m_s_2 / 1e5
bar_to_amu_nm_ps_2 = 1 / amu_nm_ps_2_to_bar
print("Conversion factor : ", bar_to_amu_nm_ps_2)

def get_pressure_values(filename):
    data_values = []    
    with open(filename, "r") as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith('#') or line.startswith('@'):
                continue
            # values = (line.strip(' ').strip('\n').split('  '))
            x, y = map(float, line.split())
            data_values.append([x, y])
    data_values = pd.DataFrame(data_values, columns = ['Time', 'Pressure'])
    # data_values['Pressure'] *= bar_to_amu_nm_ps_2
    return data_values

pxx = get_pressure_values("./prod/pxx.xvg")
pyy = get_pressure_values("./prod/pyy.xvg")
pzz = get_pressure_values("./prod/pzz.xvg")

# time offset to consider averages from
time_offset = 100
indices = pxx['Time'] >= time_offset
pxx_values = pxx[indices]['Pressure']
pyy_values = pyy[indices]['Pressure']
pzz_values = pzz[indices]['Pressure']
pxx_avg = pxx_values.mean()
pyy_avg = pyy_values.mean()
pzz_avg = pzz_values.mean()

print(pxx_avg)
print(pyy_avg)
print(pzz_avg)
Lz = 176.16

surface_tension =  (Lz) * (pzz_values - 0.5 * (pxx_values + pyy_values))
print("Surface tension value : \n", surface_tension)
plt.plot(pzz['Time'][indices], surface_tension)
plt.show()
