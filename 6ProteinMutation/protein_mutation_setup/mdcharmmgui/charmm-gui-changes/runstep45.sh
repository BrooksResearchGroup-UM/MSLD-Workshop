#!/bin/bash
#SBATCH -A rhayes1_lab
#SBATCH -p free-gpu
#SBATCH --time=240
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1

export CHARMMDIR=~rhayes1/CHARMM_EXE
source $CHARMMDIR/modules
export CHARMMEXEC=$CHARMMDIR/gnu/charmm

$CHARMMEXEC -i step4_equilibration.inp -o step4_equilibration.out
$CHARMMEXEC cnt=1 -i step5_production.inp -o step5_production_1.out
