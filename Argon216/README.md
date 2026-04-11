# Argon 216

> Following steps exactly as per [this tutorial](https://becksteinlab.physics.asu.edu/pages/courses/2013/SimBioNano/11/SimulatingliquidArgonwithGromacs.pdf).

## Steps

```bash
Simulating liquid Argon with Gromacs Documentation, Release 1.0
mkdir -p NAME/p11 && cd NAME/p11
curl -O http://becksteinlab.physics.asu.edu/pages/courses/2013/SimBioNano/11/Argon_GMX.tar.gz
tar -zxvf Argon_GMX.tar.gz
python scripts/generate_lattice.py -l cubic 512
# rename pdb generated to Ar512.pdb
mkdir argon_94K
gmx_mpi grompp -p Argon_GMX/Argon/argon.top -c Ar512.pdb -f Argon_GMX/Argon/LJ.mdp -o argon_94K/md.tpr

# to calculate rdf
gmx_mpi rdf -s md.tpr -f md.xtc -o rdf -xvg none
# select 0 then 0 then Ctrl-D

# to calculate msd
gmx_mpi msd -s md.tpr -f md.xtc -o msd -xvg none
# select 0 then Ctrl-D
```

