# Simulating Alanine dipeptide in gromacs

> Alanine dipeptide molecule serves as the "hello world" of MD simulations. It is the simplest molecule where we can perform MD simulations and analyze trajectories.

Follow [this](https://cbp-unitn.gitlab.io/QCB/tutorial2_gromacs.html) tutorial. I'm noting the condensed set of commands I will be using in this tutorial.

Commands:

Without solvent

```bash
gmx pdb2gmx -f tutorial_gromacs/alanine_dip/dip.pdb -ff amber99sb-ildn  -water none -o dip_processed.gro
gmx editconf -f dip_processed.gro -o dip_newbox.gro -c -d 1.0 -bt cubic
gmx grompp -f tutorial_gromacs/alanine_dip/minim.mdp -c dip_newbox.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
gmx grompp -f tutorial_gromacs/alanine_dip/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt
gmx grompp -f tutorial_gromacs/alanine_dip/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt

# Running NVE simulations
gmx grompp -f tutorial_gromacs/alanine_dip/nve.mdp -c em.gro -r em.gro -p topol.top -o nve.tpr
gmx mdrun -deffnm nve
```

- Optional: Convert tpr to human readable format:
```bash
gmx dump -s nvt.tpr > nvt_verbose_vacuum.tpr
```
