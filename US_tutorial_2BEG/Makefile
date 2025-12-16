# Tutorial borrowed from http://www.mdtutorials.com/gmx/umbrella/01_pdb2gmx.html

DOWNLOADS=0.downloaded_files
PDB_NAME=2BEG_model1_capped.pdb
STRUCTURES=1.generated_files

gen:
	mkdir -p ${STRUCTURES}
	gmx_mpi pdb2gmx -f ${DOWNLOADS}/${PDB_NAME} -ignh -ter -ff gromos53a6 -water spc -o ${STRUCTURES}/complex.gro
# 	Select "None" for the N-termini, and "COO-" for the C-termini for each chain.