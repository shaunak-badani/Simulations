#!/usr/bin/env python
# generate a lattice of atoms for a given density
# OB for ASU PHY494/598
# Placed into the Public Domain

import numpy as np
import mdIO as io

# see for instance share/gromacs/top/oplsaa.ff/atomtypes.atp
molecular_weight = {   # atomic units u == g/mol
    'He': 4.0026,
    'Ne': 20.1797,
    'Ar': 39.948,
    'Kr': 83.798,
    'Xe': 131.293,
    }

N_A = 6.022e23   # in mol**(-1)

def rho2specificvolume(rho, MW):
    """
    rho = m/V     --- rho in g.cm**(-3)
    n = N/N_A
    M = m/n       --- M in g/mol

    rho = N M / NA V

    N=1

    1/V = rho*NA/M

    cm**3 = (1e8 A)**3
    """
    return rho * 1e-24 * N_A / MW

# different ways to compute the lattice vectors
# timing see comments

# for testing
# import timeit
# timeit.repeat("generate_lattice.lattice_vector_3(generate_lattice.n,generate_lattice.a)", "import generate_lattice", number=100000)
#n = np.array([ 1.,  0.,  3.])
#a = np.array([[ 10.,   0., -10.], [  0.,   1.,   0.], [ -2.,   3.,   5.]])


def lattice_vector_0(n, a):
    """Generate Bravais lattice vector n*a.

    r = n1*a1 + n2*a2 + n3*a3

    where a1 = (a1x, a1y, a1z) etc are the unit cell vectors.

    :Arguments:
       *n*
          array([n1, n2, n3])
       *a*
          3x3 array, where a[0] = first unit cell vector, a[1] second
          etc
    """
    # 2.23
    # traditional element-by-element loop
    r = np.zeros(3)
    for i in range(3):
        # loop over basis vectors
        for j in range(3):
            # loop over coordinates x, y, z
            r[j] += n[i] * a[i,j]
    return r

# eliminate the x-y-z loop by making use of array's capabilities

def lattice_vector_1(n, a):
    # 2.99
    r = np.zeros(3)
    for i in range(3):
        # loop over basis vectors
        r += n[i] * a[i]
    return r

def lattice_vector_2(n, a):
    # 5.43
    r = np.sum([n[i] * a[i] for i in range(3)], axis=0)
    return r

def lattice_vector_3(n, a):
    # 1.86
    # 2.25 with asarray
    #n = np.asarray(n)
    #a = np.asarray(a)
    return np.sum(n[:,np.newaxis]*a, axis=0)  # uses broadcasting!

def lattice_vector_4(n, a):
    # 1.70
    return np.sum((n*a.T).T, axis=0)  # uses broadcasting!

def lattice_vector_5(n, a):
    # 2.59
    return sum((n*a.T).T)  # uses broadcasting and simple addition



# lattice vectors
bravais_lattice = {
    'cubic': np.array([[1,0,0], [0,1,0], [0,0,1]]),
    }

# to be used with a cubic bravais lattice: will generate
# the corresponding Bravais lattice (multiply with lattice constant a!)
basis = {
    'cubic': np.array([[0,0,0]]),
    'bcc': np.array([[0,0,0], [0.5,0.5,0.5]]),
    'fcc': np.array([[0,0,0], [0.5,0.5,0], [0.5,0,0.5], [0,0.5,0.5]]),
}

def simulation_cell_parameters(rho, N, lattice="fcc"):
    """Return unitcell length a, number of unitcells in each direction.

    V = N/rho    (rho in 1/Length = volume per particle)

    N_cells = ceil(N/N_u)  (actually, adjusted to be a perfect cube)

    M = N_cells**(1/3)     (repetition of unitcell along each axis)

    V_u = V/N_cells        (volume of a single unit cell)

    a = (V_u)**(1/3)       (length of unit cell, all on cubic lattice!)

    N_u: number of atoms per unitcell: cubic=1, bcc=2, fcc=3

    Note: N is adjusted to fill a cubic simulation cell

    :Returns: (a, M, N_u, N)
    """
    N_u = len(basis[lattice])
    N_cells = (round(np.ceil(N/float(N_u))**(1/3.)))**3
    N = int(N_u * N_cells)
    M = int(round(N_cells**(1./3)))
    V = N/rho
    V_u = V/N_cells
    a = V_u**(1./3.)

    print("--- Lattice parameters " + 36*"-")
    print("lattice: %s    rho = %f A**-3" % (lattice, rho))
    print("number of unitcells: %d" % N_cells)
    print("number of atoms:     %d" % N)
    print("a = %f A, M = %d, L=%g A, V_u = %g A**3" % (a, M, M*a, V_u))
    print(60*"-")

    return a, M, N_u, N


def generate_lattice(rho, N, atomname="Ar", lattice="fcc"):
    a, M, N_u, N = simulation_cell_parameters(rho, N, lattice=lattice)
    atoms = N * [atomname]
    box = np.array([M*a, M*a, M*a])
    coordinates = np.zeros((N,3))

    bvecs = bravais_lattice['cubic']
    b = basis[lattice]

    iatom = 0
    for i in range(M):
        for j in range(M):
            for k in range(M):
                n = np.array([i,j,k])
                v = lattice_vector_4(n, bvecs)
                x = a*(b + v)    # this is a N_u * 3 array!
                coordinates[iatom:iatom+N_u] = x
                iatom += N_u
    return atoms, coordinates, box


def make_lattice(rho, N, filename=None, fileformat="pdb", atomname="Ar", lattice="fcc"):
    atoms, coordinates, box = generate_lattice(rho, N, atomname=atomname, lattice=lattice)
    if filename is None:
        filename = "%s_rho=%g_N=%d_%s.%s" % (atomname, rho, len(atoms), lattice, fileformat)
    io.write_single(filename, atoms, coordinates, box)
    print("Wrote lattice to %r" % filename)
    return filename

def random_displace(coordinates, box, **kwargs):
    import numpy.random
    n = numpy.round(kwargs.pop('sigma', 3)*numpy.random.standard_normal(coordinates.shape))
    v = n * box     # works for cubic unitcells
    return coordinates + v

def make_displaced_lattice(atoms, coordinates, box, filename):
    newcoord = random_displace(coordinates, box)
    io.write_pdb(filename, atoms, newcoord, box)
    print("Wrote displaced lattice to %r" % filename)
    return filename


usage = """usage: %%prog [options] N

Generate a PDB or XYZ file with approximately N atoms, forming a lattice with
the given density rho.

lattice can be one of %r
""" % basis.keys()


if __name__ == "__main__":
    import sys
    import os
    from optparse import OptionParser

    parser = OptionParser(usage=usage)
    parser.add_option("-d", "--density", dest="density", default="1.374",
                      metavar="FLOAT", type="float",
                      help="density in g/cm**3 [%default]")
    parser.add_option("-n", "--atom-name", dest="atom", default="Ar",
                      metavar="STRING",
                      help="atom name [%default]")
    parser.add_option("-l", "--lattice", dest="lattice", default="fcc",
                      choices=list(basis.keys()),
                      metavar="LATTICE",
                      help="generate regular Bravais lattice [%default]")
    parser.add_option("-o", "--output", dest="output", default=None,
                      metavar="FILE",
                      help="write output to PDB file FILE; default is "
                      "to generate a filename")
    parser.add_option("-x", "--xyz", dest="xyz", action="store_true", default=False,
                      help="write output to XYZ file instead of PDB")
    parser.add_option("-w", "--molecular-weight", dest="mw", default=None,
                      metavar="FLOAT", type="float",
                      help="molecular weight of atom in g/mol; if none "
                      "supplied an internal lookup table is tried.")

    opts, args = parser.parse_args()

    try:
        N = int(args[0])
    except:
        print("Argument must be the approximate number of particles N. See --help")
        sys.exit(1)

    if not opts.mw:
        try:
            opts.mw = molecular_weight[opts.atom]
        except KeyError:
            print("No molecular weight found (known: %r); supply with --molecular-weight" % molecular_weight.keys())
            sys.exit(1)

    fileformat = "pdb"
    output = opts.output
    if opts.xyz:
        fileformat = "xyz"
        if opts.output is not None:
            root, ext = os.path.splitext(opts.output)
            output = root + ".xyz"

    v0 = rho2specificvolume(opts.density, opts.mw)
    make_lattice(v0, N, filename=output, fileformat=fileformat,
                 atomname=opts.atom, lattice=opts.lattice)

