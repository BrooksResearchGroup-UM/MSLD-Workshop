#!/bin/bash

source ~/.bashrc
source activate pychm_msld_wksh

python -c "import alf; alf.initialize(engine='pycharmm')"
python -c "import alf; alf.runflat(1,100,13000,39000,engine='pycharmm')"
python -c "import alf; alf.runflat(101,113,125000,375000,engine='pycharmm')"
