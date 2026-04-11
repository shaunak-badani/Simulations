import numpy as np
import matplotlib.pyplot as plt

rdf = np.loadtxt("argon_94K/rdf.xvg", comments=['#', '@'], unpack=True)
t = rdf[0]
rdf= rdf[1]

fig = plt.figure(figsize=(6,5))

plt.plot(t, rdf, 'r-', lw=2, label="rdf")
plt.ylabel("rdf")
plt.xlabel("r (nm)")
plt.axvline(x = 0.37, linewidth = 2, color = 'g', label = '3.7 Angstrom')
plt.axvline(x = 0.7, linewidth = 2, color = 'blue', label = '7.0 Angstrom')
plt.axvline(x = 1, linewidth = 2, color = 'orange', label = '10.0 Angstrom')
plt.axhline(y = 1, linewidth = 1.2, color = 'red', label = '1 rdf line')
plt.legend(loc = 'best')
plt.savefig('rdf.png')
