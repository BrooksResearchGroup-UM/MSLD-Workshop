#! /bin/bash

module load anaconda/3.5.3.0
source activate pycharmm47

python ../../charge_correction.py $1
