#! /bin/bash

source ~/.bashrc
source activate pychm_msld_wksh

python -c "import alf; alf.postprocess($i,$eqS,$S,$N,$skipE,True,engine='pycharmm')"
