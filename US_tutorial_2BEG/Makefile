# Tutorial borrowed from http://www.mdtutorials.com/gmx/umbrella/01_pdb2gmx.html

DOWNLOADS=0.downloaded_files
PDB_NAME=2BEG_model1_capped.pdb
STRUCTURES=1.generated_files
GMX=gmx_mpi

gen: ${DOWNLOADS}/${PDB_NAME}
	mkdir -p ${STRUCTURES}
	${GMX} pdb2gmx -f ${DOWNLOADS}/${PDB_NAME} -ignh -ter -ff gromos53a6 -water spc -o ${STRUCTURES}/complex.gro
# 	Select "None" for the N-termini, and "COO-" for the C-termini for each chain.
# 	This command generates posre_*.itps and topol_*.itps

# Run after gen
add-position-restraints: posre_Protein_chain_B.itp
	@printf "\n#ifdef POSRES_B\n#include \"posre_Protein_chain_B.itp\"\n#endif\n" >> topol_Protein_chain_B.itp
	
add-box: ${STRUCTURES}/complex.gro
# 	We will be pulling a total distance of 5.0 nm in a 12.0-nm box, to avoid the complications described above. 
#   The center of mass of the protofibril will be placed at (3.280, 2.181, 2.4775) in a box of dimensions 6.560 x 4.362 x 12. 
# 	Use editconf to place the protofibril at this location:
	${GMX} editconf -f ${STRUCTURES}/complex.gro -o ${STRUCTURES}/newbox.gro -center 3.280 2.181 2.4775 -box 6.560 4.362 12


solvate: ${STRUCTURES}/newbox.gro topol.top
	${GMX} solvate -cp ${STRUCTURES}/newbox.gro -cs spc216.gro -o ${STRUCTURES}/solv.gro -p topol.top

add-ions: ${DOWNLOADS}/ions.mdp ${STRUCTURES}/solv.gro topol.top
	${GMX} grompp -f ${DOWNLOADS}/ions.mdp -c ${STRUCTURES}/solv.gro -p topol.top -o ${STRUCTURES}/ions.tpr
	${GMX} genion -s ${STRUCTURES}/ions.tpr -o ${STRUCTURES}/solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.1