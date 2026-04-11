# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# input/output module (r/w XYZ, w PDB)
#
# Oliver Beckstein 2013, for ASU Simulation Approaches to Bio- and Nanophysics
# PHY598/PHY494/CHM598
#
# Placed into the Public Domain

"""
Functions to read and write XYZ files.

Resources: the XYZ format was taken from
http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/xyzplugin.html and
is therefore compatible with VMD.

Format implemented here (elements in square brackets are optional)::

   N
   [box Lx Ly Lz   [frame N]]   molecule_name
   atom1 x y z
   atom2 x y z
   ...
   atomN x y z
"""
# from __future__ import with_statement   # for python 2.5

import numpy as np
import os

def read_xyz(filename):
    """Read XYZ file *filename*.

    Returns atoms (as array) and coordinates as Nx3 numpy array,
    and box [Lx, Ly, Lz].
    """
    with open(filename) as xyz:
        n_atoms = int(xyz.readline())
        title = xyz.readline().strip()
        if title.startswith("box"):
            # parses:
            # box 12.4 20.03 33 other stuff that is ignored
            # --> box = array([12.4, 20.03, 33.])
            fields = title.split()
            box = np.array([float(x) for x in fields[1:4]])
        else:
            box = np.zeros(3)
        iatom = 0
        # pre-allocate arrays and fill
        #atoms = np.array(n_atoms * [None])  # object array so that we can use strings of arbitrary length
        atoms = np.zeros(n_atoms, dtype="S6") # strings of maximum length 6
        coordinates = np.zeros((n_atoms, 3))
        for line in xyz:
            atom, x, y, z = line.split()
            try:
                atoms[iatom] = atom
                coordinates[iatom, :] = float(x), float(y), float(z)
            except IndexError:
                raise ValueError("There are more coordinates to be read than indicated in the header.")
            iatom += 1

    if iatom != n_atoms:
        raise ValueError("number of coordinates read %d does not agree with number of atoms stated in file %d" \
                             % (iatom, n_atoms))

    return atoms, coordinates, box

def write_xyz(filename, atoms, coordinates, box, title="simulation", mode="w"):
    """Write atoms and coordinates to XYZ file *filename*.

    :Arguments:
       *filename*
           name of the output file
       *atoms*
           list of the N atom names
       *coordinates*
           coordinates as Nx3 array
       *box*
           box lengths (Lx Ly Lz)
    """
    with open(filename, mode) as xyz:
        xyz.write("%d\n" % len(atoms))
        xyz.write("box %g %g %g  %s\n" % (box[0], box[1], box[2], title))
        for i in range(len(atoms)):
            x, y, z = coordinates[i]
            xyz.write("%8s  %10.5f %10.5f %10.5f\n" % (atoms[i], x, y, z))


def write_xyz_frame(xyz, atoms, coordinates, box, frame, title="simulation"):
    """Write frame in XYZ format to open file object *xyz*.

    :Arguments:
       *filename*
           name of the output file
       *atoms*
           list of the N atom names
       *coordinates*
           coordinates as Nx3 array
       *box*
           box lengths (Lx Ly Lz)
       *frame*
           frame number (or time)
    """
    xyz.write("%d\n" % len(atoms))
    xyz.write("box %g %g %g  frame %d  %s\n" % (box[0], box[1], box[2], frame, title))
    for i in range(len(atoms)):
        x, y, z = coordinates[i]
        xyz.write("%8s  %10.5f %10.5f %10.5f\n" % (atoms[i], x, y, z))


def write_single(*args, **kwargs):
    """Write coordinate file depending on file extension.

    Arguments must be appropriate for the underlying writer. Typically you need ::

       write_single(filename, atoms, coordinates, box, title="title text")


    Supported formats:

    - filename.pdb: simple PDB file format
    - filename.xyz: xyz format

    """
    writers = {'.xyz': write_xyz,
               '.pdb': write_pdb,
               }
    filename = args[0]
    root, ext = os.path.splitext(filename)
    return writers[ext](*args, **kwargs)

def write_pdb(filename, atoms, coordinates, box, mode="w"):
    """Very primitive PDB writer.

    :Arguments:
       *filename*
           name of the output file
       *atoms*
           list of the N atom names
       *coordinates*
           coordinates as Nx3 array (must be in Angstroem)
       *box*
           box lengths (Lx Ly Lz) (must be in Angstroem)

    See http://www.wwpdb.org/documentation/format32/sect9.html
    """
    with open(filename, mode) as xyz:
        xyz.write("HEADER    simple PDB file with %d atoms\n" % len(atoms))
        xyz.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n" %
                  (box[0], box[1], box[2], 90., 90., 90., "P 1", 1))
        for i in range(len(atoms)):
            serial = i+1
            name = resName = atoms[i]
            chainID = 'A'
            resSeq = i+1
            iCode = ''
            occupancy = 1.0
            tempFactor = 0.0
            x, y, z = coordinates[i]
            xyz.write("ATOM  %(serial)5d %(name)-4s %(resName)-4s%(chainID)1s%(resSeq)4d%(iCode)1s   %(x)8.3f%(y)8.3f%(z)8.3f%(occupancy)6.2f%(tempFactor)6.2f\n" % vars())


if __name__ == "__main__":

    testfiles = ["test_06_0.xyz", "test_06_1.xyz"]
    for f in testfiles:
        print("Testing %r" % f)
        atoms, coordinates, box = read_xyz(f)
        print(atoms[3:8])
        print(coordinates[3:8])
