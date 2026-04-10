```bash
packmol -i pack.inp
# The values in pack.inp need to be scaled by 10
# editconf scales down by 10 assuming the positions are in angstrom values
# in the pdb file
gmx_mpi editconf -f prep_files/argon.pdb -o prep_files/argon.gro -box 11.011 11.011 176.16 -noc
mkdir em
gmx_mpi grompp -f minim.mdp -c prep_files/argon.gro -p topol.top -o em/argon.tpr
mkdir nvt
gmx_mpi grompp -c em/argon.gro -p topol.top -f nvt.mdp -o nvt/argon_nvt.tpr
mkdir prod
gmx_mpi grompp -c nvt/argon_nvt.gro -p topol.top -f prod.mdp -o prod/argon_prod.tpr
```