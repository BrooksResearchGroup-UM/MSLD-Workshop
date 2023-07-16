#! /bin/bash

CHARMMDIR=~rhayes1/CHARMM_EXE
source $CHARMMDIR/modules
CHARMMEXEC=$CHARMMDIR/gnu/charmm

name=T16
source /dfs8/rhayes1_lab/rhayes1/75_MSLD/06_alfdebug/alf-20230707/setupenv

mkdir prep
cp alf_info.py nsubs prep/
python -c "import alf; alf.initialize()"

cp ../mdsetup/scratch3/annealed.crd prep/system.crd
sed "s/BOXSUB/`cat ../mdsetup/scratch3/box.dat`/g" $name.inp > prep/$name.inp
cp -r ../mdsetup/toppar aa_stream prep/
cp alchemical_definitions.inp prep/
cp ../mdsetup/prep/nbond* prep/

mkdir scratch4
mpirun -np 1 $CHARMMEXEC sysname=\"$name -i 4SetupBlock.inp

# util/rotate_c2b.py prep/minimized_charmm.crd > prep/minimized.crd
