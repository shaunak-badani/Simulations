import numpy as np
import matplotlib.pyplot as plt

msd = np.loadtxt("argon_94K/msd.xvg", unpack=True)
t = msd[0]
msd = msd[1]

fig = plt.figure(figsize=(6,5))

plt.scatter(t, msd, 7,  label="rdf")
plt.ylabel("Mean Square displacement (nm)")
plt.xlabel("time (ps)")

# plt.savefig('msd.png')

slope, intercept = np.polyfit(t, msd, 1)

print(f"Equation: y = {slope:.7f}x + {intercept:.7f}")
line_fit = slope * t + intercept
plt.plot(t, line_fit, 'g-', lw=2, label = "line_fit")
plt.legend(loc = 'best')
plt.show()

# slope value is in (nm^2) / ps
nm_to_cm = 1e-7
ps_to_s = 1e-12
diffusion_coefficient = slope * nm_to_cm**2 / (6 * ps_to_s)


print("Obtained diffusion coefficitent : ", diffusion_coefficient)
