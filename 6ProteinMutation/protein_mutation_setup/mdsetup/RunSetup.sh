#! /bin/bash

CHARMMDIR=~rhayes1/CHARMM_EXE
source $CHARMMDIR/modules
CHARMMEXEC=$CHARMMDIR/gnu/charmm

name=`cat name`

./0Preprocess.sh
mpirun -np 1 $CHARMMEXEC -i 1Setup.inp
./2Solvate.sh
mpirun -np 1 $CHARMMEXEC sysname=\"$name -i 3Setup.inp
