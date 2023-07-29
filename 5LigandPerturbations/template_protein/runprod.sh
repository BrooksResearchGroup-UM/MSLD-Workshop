#!/bin/bash

source ~/.bashrc
source activate pychm_msld_wksh

iend=$(( $SLURM_ARRAY_TASK_ID * $nitt ))
ibeg=$(( $iend - $nitt ))
python -c "import alf; alf.runprod($step,'$p',$ibeg,$iend,engine='pycharmm')"
