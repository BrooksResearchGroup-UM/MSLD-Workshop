#! /bin/bash

# Load anaconda and slurm modules
module load slurm
module load anaconda/3.5.3.0

# Load conda environment with pyCHARMM installed
source activate pycharmm47

# Postprocess trajectories
for dir in `ls physical_ligands`; do
    echo $dir
    cd physical_ligands/$dir

    python ../../get_free_energy.py
    
    cd ../../
done
