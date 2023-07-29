#! /bin/bash

# Load anaconda and slurm modules
module load slurm
module load anaconda/3.5.3.0

# Load conda environment with pyCHARMM installed
source activate pycharmm47

# Add SLURM options 
EXCLUDE="--exclude=gollum003"
RE_Q="--no-requeue"

# For each ligand, run a simulation with both charge sets
for dir in `ls physical_ligands`; do
    echo $dir
    cd physical_ligands/$dir
    sbatch --time=7-00:00:00 --ntasks=1 --tasks-per-node=1 --cpus-per-task=4 -p gpu --gres=gpu:1 --export=ALL $EXCLUDE ../../md_job.sh ff
    sbatch --time=7-00:00:00 --ntasks=1 --tasks-per-node=1 --cpus-per-task=4 -p gpu --gres=gpu:1 --export=ALL $EXCLUDE ../../md_job.sh crn
    cd ../../ 
done
