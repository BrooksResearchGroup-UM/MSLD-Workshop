grep -i $1 step2_solvator.pdb > $1.pdb
echo TER >> $1.pdb
echo END >> $1.pdb
change=`echo $1 | tr '[a-z]' '[A-Z]'`
sed -i'.pdb' "s/IONS/${change}/" $1.pdb

