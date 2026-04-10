import numpy as np
from scipy.special import erfc


def phi6(beta, cutoff):
    a = beta * cutoff
    return (1 + a**2 + 0.5 * a**4) * np.exp(-a**2)

def elec_tol(beta, cutoff):
    a = beta * cutoff
    return erfc(a)

ewald_rtol_lj = phi6(10, 0.28)
ewald_rtol = elec_tol(10, 0.255)
print("ewald-rtol-lj value : ", ewald_rtol_lj)
print("ewald-rtol value : ", ewald_rtol)