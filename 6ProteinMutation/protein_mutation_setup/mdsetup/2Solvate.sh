#! /bin/bash

# module load mmtsb

INDIR=scratch1
OUTDIR=scratch2

name=`cat name`

mkdir $OUTDIR

util/solvate.py $INDIR/system.crd rhdo > $OUTDIR/systemS.crd
util/ionizeBymM.py $OUTDIR/systemS.crd `cat scratch1/charge.dat` SOD:100.0=CLA:100.0 > $OUTDIR/systemI.crd
util/rotate_b2c.py $OUTDIR/systemI.crd > $OUTDIR/system.crd

sed "s/BOXSUB/`head -n 1 abcabc.dat`/g" ${name}0.inp > scratch2/${name}0.inp

mkdir prep
cp -r toppar $OUTDIR/system.crd prep/
util/makenbondchm.py `head -n 1 abcabc.dat` > prep/nbond.str
util/makenbondbld.py `head -n 1 abcabc.dat` > prep/nbond_blade.str
