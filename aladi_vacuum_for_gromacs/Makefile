CLEAN_PATTERNS = \
    "*.log" "*.edr" "*.trr" "*.cpt" "*.tprls" "*.out" \
    "*.prev" "*.bak" "\#*\#" \
    "traj.trr"

.PHONY: clean

clean:
	@echo "Removing GROMACS-generated files..."
	@for p in $(CLEAN_PATTERNS); do \
	    rm -rf $$p 2>/dev/null || true; \
	done
	@echo "Done."

# make organize dir=<name>
organize:
	@if [ -z "$(dir)" ]; then \
	    echo "Usage: make organize dir=<name>"; exit 1; \
	fi
	mkdir -p tmp
	mv $(dir).edr $(dir).log $(dir).trr tmp/
	mv tmp $(dir)

move-forces:
	mkdir forces
	mv angle_forces.dat bonded_forces.dat energies.dat \
		periodic_improper_dihedral_forces.dat proper_dihedral_forces.dat forces/

delete-trajectories:
	find . -type f -name *.trr -delete
