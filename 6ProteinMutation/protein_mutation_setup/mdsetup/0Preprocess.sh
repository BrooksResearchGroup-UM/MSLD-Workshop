
PDB=1pga.pdb
OUTDIR=scratch0

mkdir $OUTDIR

grep "^ATOM\|^TER\|^HETATM\|^END" $PDB | ./FixHis.sh > $OUTDIR/clean.pdb 

util/ConvPDB.py $OUTDIR/clean.pdb > $OUTDIR/system.crd

CHARMMDIR=/share/crsp/lab/rhayes1/share/git/CHARMM/dev-release
mkdir toppar
cp $CHARMMDIR/toppar/par_all36m_prot.prm toppar/
cp $CHARMMDIR/toppar/par_all36_cgenff.prm toppar/
cp $CHARMMDIR/toppar/top_all36_prot.rtf toppar/
cp $CHARMMDIR/toppar/top_all36_cgenff.rtf toppar/
cp $CHARMMDIR/toppar/toppar_water_ions.str toppar/

