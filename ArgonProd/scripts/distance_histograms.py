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

GRO = "./prep_files/argon.gro"
u = mda.Universe(GRO)

# convert to nm if gro file
positions = u.trajectory[0].positions / 10

N = positions.shape[0]
print(positions[:1])
box = u.dimensions[:3]
distances = []
for i in range(N):
    for j in range(i + 1, N):
        dist = positions[i] - positions[j]
        for m in range(3):
            if dist[m] >= box[m] / 2:
                dist[m] -= box[m]
            if dist[m] <= -box[m] / 2:
                dist[m] += box[m]
        distance = np.linalg.norm(dist)
        distances.append(distance)

distances = np.array(distances)

p, be = np.histogram(distances, bins = 30)

coords = (be[1:] + be[:-1])/2
# plt.plot(coords, p)
# plt.scatter(coords, p)


pot_contribution = 4 * 1 * ((1 / coords)**12 - (1 / coords)**6)

plt.plot(coords, pot_contribution)
plt.scatter(coords, pot_contribution)

ax = plt.gca()
current_ticks = list(ax.get_xticks())
x_min = coords[0]
ax.set_xticks(current_ticks + [x_min])
plt.show()