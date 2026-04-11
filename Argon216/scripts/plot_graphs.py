import numpy as np
import matplotlib.pyplot as plt

energy = np.loadtxt("argon_94K/energy.xvg", comments=['#', '@'], unpack=True)
t = energy[0]
U, KE, Etot, T = energy[1:5]


# layout with 3 graphs in 2x2
fig = plt.figure(figsize=(6,5))
axE = fig.add_subplot(221)
axT = fig.add_subplot(223)
axhT = fig.add_subplot(224)

# energies
axE.plot(t, Etot, 'k-', lw=2, label="total")
axE.plot(t, U, 'b-', lw=2, label="potential")
axE.plot(t, KE, 'r-', lw=2, label="kinetic")
axE.set_ylabel("energy (kJ/mol)")
axE.legend(loc="best")
# plot T
# timeseries

axT.plot(t, T, 'r-', lw=2, label="temperature")
axT.set_ylabel("temperature (K)")
axT.set_xlabel("time (ps)")

# histogram
h,edges = np.histogram(T, bins=20, density=True)
midpoints = 0.5*(edges[1:] + edges[:-1])
# plot the histogram as "marginal" over the time series
axhT.plot(h, midpoints, 'r-')

axhT.set_ylim(axT.get_ylim())
axhT.get_xaxis().set_visible(False) # hide whole axis (counts/density)
# save figure as high resolution PNG
fig.savefig("energies.png", dpi=300)
# mean and fluctuations
print(Etot.mean(), Etot.std())