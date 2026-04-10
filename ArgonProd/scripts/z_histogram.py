"""
Run from root
python scripts/distance_histograms.py
"""

import MDAnalysis as mda

import warnings
# suppress some MDAnalysis warnings about PSF files
warnings.filterwarnings('ignore')
import numpy as np

from matplotlib import pyplot as plt

GRO = "./em/argon.gro"
u = mda.Universe(GRO)

# convert to nm
positions = u.trajectory[0].positions / 10

N = positions.shape[0]

z = positions[:, 2]

p, be = np.histogram(z, density=True)

coords = (be[1:] + be[:-1])/2
plt.plot(coords, p)
plt.scatter(coords, p)

ax = plt.gca()
current_ticks = list(ax.get_xticks())
x_min = np.min(z)
x_max = np.max(z)
ax.set_xticks(current_ticks + [x_min, x_max])
plt.show()