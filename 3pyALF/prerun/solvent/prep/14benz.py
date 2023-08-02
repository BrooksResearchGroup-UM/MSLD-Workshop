##
## pyCHARMM script drafted from many examples
##
##   ** MSLD with BLADE (no OMM support) **
##

##############################################
# Set up global parameters

# variables
box = 32.964000
pmegrid = 32

# msld variables
fnex = 5.5

nsites=len(nsubs)

##############################################
# Read in toppar files, coordinate files, etc.

# toppar files
read.rtf('prep/top_all36_msld.rtf')
read.rtf('prep/full_ligand.rtf',append=True)
read.prm('prep/par_all36_msld.prm',flex=True)
read.prm('prep/full_ligand.prm',flex=True,append=True)

ligseg = 'LIG' 
resnum = '1'
read.sequence_pdb('prep/full_ligand.pdb')
gen.new_segment(ligseg,setup_ic=True)
read.pdb('prep/full_ligand.pdb',resid=True)

solvated=True
if solvated==True:
  read.sequence_pdb('prep/solvent.pdb')
  gen.new_segment('WT00',setup_ic=True,angle=False, dihedral=False)
  read.pdb('prep/solvent.pdb',resid=True)

# bomblev -1  ! JZV

#  Hybrid Ligand Setup
#  (1) read in patch files
for i in range(len(nsubs)):
  for j in range(nsubs[i]):
    read.rtf(f'prep/site{i+1}_sub{j+1}_pres.rtf',append=True)

# (3) add atoms for each substituent
pycharmm.lingo.charmm_script('ic generate')

pycharmm.lingo.charmm_script('autogen nopatch') # Don't autogen every step
for i in range(len(nsubs)):
  for j in range(nsubs[i]):
    pycharmm.lingo.charmm_script(f'patch p{i+1}_{j+1} {ligseg} {resnum} setup')
    read.pdb(f'prep/site{i+1}_sub{j+1}_frag.pdb',resid=True)
    ic.prm_fill(replace_all=False) # ic param
    ic.build() # ic build
select.store_selection('togenerate',pycharmm.SelectAtoms(seg_id=ligseg,res_id=resnum))
pycharmm.lingo.charmm_script('auto angle dihe sele .bonded. .bonded. .bonded. togenerate end')

# (2) delete atoms in common core ligand 
#    atoms taken from site1_sub1.txt site2_sub1.txt
select.store_selection('todelete',pycharmm.SelectAtoms().by_res_and_type(ligseg,resnum,'C4 C5 H4 H5'))
pycharmm.lingo.charmm_script('delete atom select todelete end')

# Hybrid Ligand Block
# Substituent definitions
for site in range(len(nsubs)):
  for sub in range(nsubs[site]):
    selname='site'+str(site+1)+'sub'+str(sub+1)
    # extract alchem patch atoms from patch file
    sub_atoms=[]
    rtffile = 'prep/site'+str(site+1)+'_sub'+str(sub+1)+'_pres.rtf'
    for line in open(rtffile,'r'):
      if line[0:4] == 'ATOM': sub_atoms.append(line.split()[1].upper())
    atoms_in_sub = pycharmm.SelectAtoms().by_res_and_type(ligseg,resnum,' '.join(sub_atoms))
    # saves the identified atoms in the charmm variable
    select.store_selection(selname,atoms_in_sub)

# delete angles and dihedrals between alchem groups
# pycharmm.lingo.charmm_script('auto angle dihe')
settings.set_bomb_level(-1)
for site in range(nsites):
    for sub1 in range(nsubs[site]):
        for sub2 in range(sub1+1,nsubs[site]):
            pycharmm.lingo.charmm_script('dele connectivity sele {} show end sele {} show end'
            .format(f'site{site+1}sub{sub1+1}',f'site{site+1}sub{sub2+1}'))


# # write out psf, crd, pdb files
# write.psf_card('patch.psf')
# write.coor_card('patch.crd')
# write.coor_pdb('patch.pdb')


##############################################
# Create water box & periodic images

# MODIFY to set up periodic images or SBC
coor.stat()
crystal.define_cubic(box)
# crystal.build(14.0)
pycharmm.lingo.charmm_script('open read card unit 14 name prep/cubic.xtl\ncrystal read card unit 14\nclose unit 14')
image.setup_segment(0.0, 0.0, 0.0, ligseg )
if solvated==True:
  image.setup_residue(0.0, 0.0, 0.0, 'WT00')

##############################################
# Set up BLOCK module for MSLD

# check that the system's net charge is 0
pycharmm.lingo.charmm_script('set charge = ?cgtot')
netQ = pycharmm.lingo.get_charmm_variable('CHARGE')
tol=1e-8
if (netQ > tol) or (netQ < (-1*tol)):
    print("ERROR: system net charge not equal to zero!! Exiting...")
    #pycharmm.lingo.charmm_script('stop')


# MSLD BLOCK module
# ** the block module passed by pycharmm.lingo CANNOT be divided into parts
#    it must be passed as one complete unit
# ** Therefore, multiple strings are created and passed at once to lingo
blockplusone = nblocks + 1
knoe = 118.4  # for newest version of CATS, use 118.4 for everything
# initialize block
block_init='''
!! BLOCK setup
BLOCK {}
   clear
END
BLOCK {}
'''.format(blockplusone,blockplusone)
# load blocks
block_call=''
ii=2
sub0=np.sum(nsubs)-nsubs
for site in range(nsites):
    for sub in range(nsubs[site]):
        block_call+='Call {} sele {} show end\n'.format(ii,f'site{site+1}sub{sub+1}')
        ii+=1
# for softcore atoms (not sure how to loop-create CATS atoms...?)
block_parm='''
! scat on
! scat k {}
! cats sele atom ?segid ?resid ?atomname .or. [list of atom names to cat]

qldm theta
lang temp {}
soft on
pmel ex

ldin 1 1.0  0.0  5.0  0.0  5.0'''.format(knoe,temp)
# ldin lines
block_ldin=''
sitestr=''
sub0=np.cumsum(nsubs)-nsubs
for site in range(nsites):
    for sub in range(nsubs[site]):
        if sub == 0:
            tmplmb=1.0-(0.01*(nsubs[site]-1))
        else:
            tmplmb=0.01
        iii=sub0[site]+sub
        bii=iii+2
        block_ldin+='ldin {} {:.4f} 0.0 5.0 {} 5.0\n'.format(bii,tmplmb,bias['b'][0,iii])
        sitestr+=str(site+1)+'  '
# add in exclusions with adex
block_adex=''
for site in range(nsites):
    for sub1 in range(nsubs[site]):
        for sub2 in range(sub1+1,nsubs[site]):
            block_adex+='adex {} {}\n'.format(sub0[site]+sub1+2,sub0[site]+sub2+2)

# msld parameters
block_msld='''
!!rmla bond thet dihe impr
rmla bond thet impr
msld 0  {} fnex {}
msma
'''.format(sitestr,fnex)

# msld variable biases
block_varb='ldbi {}\n'.format((5*nblocks*(nblocks-1))//2)
sub0=np.cumsum(nsubs)-nsubs
ibias=0
for si in range(nsites):
    for sj in range(si,nsites):
        for ii in range(nsubs[si]):
            for jj in range(nsubs[sj]):
                if (si != sj) or (jj > ii):
                    # bii should be iii+2
                    iii=sub0[si]+ii
                    jjj=sub0[sj]+jj
                    bii=iii+2
                    bjj=jjj+2
                    # !vbrex! these lines needed in vb.inp, not here, add to biases below
                    # !vbrex! if si==sj:
                    # !vbrex!     c_shift=2.0*(myrep-ncentral)
                    # !vbrex!     s_shift=0.5*(myrep-ncentral)
                    block_varb+='ldbv {} {} {}  6  0.00 {} 0\n'.format(ibias+1,bii,bjj,-bias['c'][iii,jjj])
                    block_varb+='ldbv {} {} {} 10 -5.56 {} 0\n'.format(ibias+2,bii,bjj,-bias['x'][iii,jjj])
                    block_varb+='ldbv {} {} {}  8 0.017 {} 0\n'.format(ibias+3,bii,bjj,-bias['s'][iii,jjj])
                    block_varb+='ldbv {} {} {} 10 -5.56 {} 0\n'.format(ibias+4,bjj,bii,-bias['x'][jjj,iii])
                    block_varb+='ldbv {} {} {}  8 0.017 {} 0\n'.format(ibias+5,bjj,bii,-bias['s'][jjj,iii])
                    ibias+=5
block_varb+='END\n'

#print('''** TEST **
pycharmm.lingo.charmm_script('''
{}
{}
{}
{}
{}
{}
{}'''.format(block_init,block_call,block_parm,block_ldin,block_adex,block_msld,block_varb))
#pycharmm.lingo.charmm_script('stop')


