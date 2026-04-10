# Simulation of water 

> Simulation details have been borrowed from the paper "Development and application of a particle-particle particle-mesh Ewald method for dispersion interactions"

```bash
mkdir prep_files
# packmol -i pack.inp
# gmx_mpi editconf -f prep_files/water.pdb -o prep_files/water.gro -box 5 5 15 -noc
gmx_mpi insert-molecules -box 5 5 5 -nmol 5000 -ci spc216.gro -o prep_files/water_norm.gro
gmx_mpi editconf -f prep_files/water_norm.gro -o prep_files/water.gro -box 5 5 15
mkdir em
gmx_mpi grompp -f em.mdp -c prep_files/water.gro -p topol.top -o em/water_min.tpr # with flexible
gmx_mpi grompp -f em2.mdp -c em/water_min.gro -p topol.top -o em/water_min2.tpr # without flexible
mkdir nvt
gmx_mpi grompp -f nvt.mdp -c em/water_min2.gro -p topol.top -o nvt/water_eq.tpr # equilibration
mkdir prod
gmx_mpi grompp -f nvt.mdp -c nvt/water_eq.gro -p topol.top -o prod/water_prod.tpr # production
```