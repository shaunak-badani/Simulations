# Tutorial borrowed from http://www.mdtutorials.com/gmx/umbrella/01_pdb2gmx.html
SHELL := /bin/bash
DOWNLOADS=0.downloaded_files
PDB_NAME=2BEG_model1_capped.pdb
STRUCTURES=1.generated_files
SYSTEM_SUFFIX=ab_42
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
	${GMX} grompp -maxwarn 2 -f ${DOWNLOADS}/ions.mdp -c ${STRUCTURES}/solv.gro -p topol.top -o ${STRUCTURES}/ions.tpr
	${GMX} genion -s ${STRUCTURES}/ions.tpr -o ${STRUCTURES}/solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.1
# 	Select group 13 (SOL)

generate-em: ${DOWNLOADS}/minim.mdp ${STRUCTURES}/solv_ions.gro topol.top
	mkdir -p em
	${GMX} grompp -maxwarn 2 -f ${DOWNLOADS}/minim.mdp -c ${STRUCTURES}/solv_ions.gro -p topol.top -o em/em_${SYSTEM_SUFFIX}.tpr

run-em: em/em_${SYSTEM_SUFFIX}.tpr
	cd em && ${GMX} mdrun -v -deffnm em_${SYSTEM_SUFFIX}

generate-npt: ${DOWNLOADS}/npt.mdp em/em_${SYSTEM_SUFFIX}.gro topol.top
	mkdir -p npt
	${GMX} grompp -maxwarn 3 -f ${DOWNLOADS}/npt.mdp -c em/em_${SYSTEM_SUFFIX}.gro -p topol.top -r em/em_${SYSTEM_SUFFIX}.gro -o npt/npt_${SYSTEM_SUFFIX}.tpr

run-npt: npt/npt_${SYSTEM_SUFFIX}.tpr
	cd npt && ${GMX} mdrun -v -deffnm npt_${SYSTEM_SUFFIX}


def-custom-index-groups: npt/npt_${SYSTEM_SUFFIX}.gro
	${GMX} make_ndx -f npt/npt_${SYSTEM_SUFFIX}.gro
	sed -i -e 's/r_1-27/Chain_A/g' -e 's/r_28-54/Chain_B/g' index.ndx
	mkdir -p steered_md && mv index.ndx steered_md/
# 	...
#  > r 1-27
#  > name 19 Chain_A -> This doesn't work anymore, replaced with sed command above after make_ndx is run
#  > r 28-54
#  > name 20 Chain_B -> This doesn't work anymore, replaced  with sed command above after make_ndx is run
#  > q

generate-pull-tpr: ${DOWNLOADS}/md_pull.mdp npt/npt_${SYSTEM_SUFFIX}.gro topol.top steered_md/index.ndx npt/npt_${SYSTEM_SUFFIX}.cpt
	${GMX} grompp -maxwarn 3 -f ${DOWNLOADS}/md_pull.mdp -c npt/npt_${SYSTEM_SUFFIX}.gro -p topol.top -r npt/npt_${SYSTEM_SUFFIX}.gro -n steered_md/index.ndx -t npt/npt_${SYSTEM_SUFFIX}.cpt -o steered_md/pull_${SYSTEM_SUFFIX}.tpr

run-pull-tpr: steered_md/pull_${SYSTEM_SUFFIX}.tpr
	cd steered_md && ${GMX} mdrun -deffnm pull_${SYSTEM_SUFFIX} -pf pullf.xvg -px pullx.xvg

extract-frames: steered_md/pull_${SYSTEM_SUFFIX}.tpr steered_md/pull_${SYSTEM_SUFFIX}.xtc
	mkdir -p umbrella_sampling
	${GMX} trjconv -s steered_md/pull_${SYSTEM_SUFFIX}.tpr -f steered_md/pull_${SYSTEM_SUFFIX}.xtc -o umbrella_sampling/conf.gro -sep
# 	save the whole system, group 0, when prompted

get-distances: umbrella_sampling/conf0.gro steered_md/pull_${SYSTEM_SUFFIX}.tpr steered_md/index.ndx
	mkdir -p distances
	i=1 ; while [[ $$i -l 501 ]] ; do \
 		${GMX} distance -s steered_md/pull_${SYSTEM_SUFFIX}.tpr -f umbrella_sampling/conf$${i}.gro -n steered_md/index.ndx -select 'com of group "Chain_A" plus com of group "Chain_B"' -oall distances/dist$${i}.xvg; \
        ((i = i + 1)) ; \
    done
	touch summary_distances.dat
	i=1 ; while [[ $$i -lt 501 ]] ; do \
		d=`tail -n 1 distances/dist$${i}.xvg | awk '{print $$2}'`; \
		echo "$${i} $${d}" >> distances/summary_distances.dat;  \
        ((i = i + 1)) ; \
	done

setup-pre-umbrella-sampling-npt: ${DOWNLOADS}/npt_umbrella.mdp umbrella_sampling/conf0.gro topol.top 
# 	Identified 23 configurations.
# 	These configurations have distance 0.5-5nm, with a spacing of 0.2nm
# 	Handpicked by looking at summary_distances file previously generated.
# 	6 175 189 195 203 211 216 223 228 235 250 260 274 289 309 322 347 364 376 399 414 430 451
	mkdir -p umbrella_sampling_simulations
	CONFIGURATIONS=(6 175 189 195 203 211 216 223 228 235 250 260 274 289 309 322 347 364 376 399 414 430 451); \
	for configuration_id in "$${CONFIGURATIONS[@]}"; do \
		echo "Processing configuration : $$configuration_id"; \
		mkdir -p umbrella_sampling_simulations/$$configuration_id; \
		${GMX} grompp -maxwarn 3 -f ${DOWNLOADS}/npt_umbrella.mdp -c umbrella_sampling/conf$$configuration_id.gro -p topol.top \
			-r umbrella_sampling/conf$$configuration_id.gro -n steered_md/index.ndx -o umbrella_sampling_simulations/$$configuration_id/npt$$configuration_id.tpr; \
	done

delete-not-necessary-configurations:
	CONFIGURATIONS=(6 175 189 195 203 211 216 223 228 235 250 260 274 289 309 322 347 364 376 399 414 430 451); \
	keep_regex=$$(printf '|conf%s.gro' "$${CONFIGURATIONS[@]}"); \
	keep_regex=$${keep_regex:1}; \
	echo $${keep_regex}; \
	ls umbrella_sampling | grep -Ev "$$keep_regex" | xargs -I{} rm umbrella_sampling/{};

run-pre-us-npt:
	CONFIGURATIONS=(6 175 189 195 203 211 216 223 228 235 250 260 274 289 309 322 347 364 376 399 414 430 451); \
	for configuration_id in "$${CONFIGURATIONS[@]}"; do \
		echo "Running simulation id : $$configuration_id"; \
		(cd umbrella_sampling_simulations/$$configuration_id && ${GMX} mdrun -deffnm npt$$configuration_id) \
	done

setup-umbrella-sampling:
	CONFIGURATIONS=(6 175 189 195 203 211 216 223 228 235 250 260 274 289 309 322 347 364 376 399 414 430 451); \
	for configuration_id in "$${CONFIGURATIONS[@]}"; do \
		echo "Processing configuration : $$configuration_id"; \
		mkdir -p umbrella_sampling_simulations/$$configuration_id; \
		${GMX} grompp -maxwarn 3 -f ${DOWNLOADS}/md_umbrella.mdp -c umbrella_sampling_simulations/$$configuration_id/npt$$configuration_id.gro -p topol.top \
			-r umbrella_sampling_simulations/$$configuration_id/npt$$configuration_id.gro -n steered_md/index.ndx -o umbrella_sampling_simulations/$$configuration_id/umbrella$$configuration_id.tpr; \
	done 

run-us:
	CONFIGURATIONS=(6 175 189 195 203 211 216 223 228 235 250 260 274 289 309 322 347 364 376 399 414 430 451); \
	for configuration_id in "$${CONFIGURATIONS[@]}"; do \
		echo "Running simulation id : $$configuration_id"; \
		(cd umbrella_sampling_simulations/$$configuration_id && ${GMX} mdrun -deffnm umbrella$$configuration_id) \
	done

delete-not-xvg:
# 	You only need xvg files for final wham analysis, delete the rest from umbrella_sampling_simulations
	find umbrella_sampling_simulations -maxdepth 2 -type f ! -name '*.xvg' -delete


delete-trajectories:
	find . -type f -name *.trr -delete
	find . -type f -name *.xtc -delete