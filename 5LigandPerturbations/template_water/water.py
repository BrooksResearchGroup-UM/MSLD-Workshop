##
## pyCHARMM script drafted from many examples
##
##   ** MSLD with BLADE (no OMM support) **
##

from prep.msld_patch import *

##############################################
# Set up global parameters

# system variables
sysname = 'water'
box = 40.988588
pmegrid = 40.0
# temp = 298.15 # Set by variables.py
builddir = './prep'

# msld variables
fnex = 5.5

# nblocks=alf_info['nblocks'] # Set by variables.py
# nsubs=alf_info['nsubs'] # Set by variables.py
nsites=len(nsubs)

# dictionary to define perturbations
# types = 'ligand','side_chain'
pert={}                    # dict of nested dicts; index = site # (zero indexed)

# ligand perturbations
pert[0]={'subs':['1','2','3','4','5','6'],  # site1_sub# (0-index = site 1)
         'segid':'LIG',
         'resid':'1',
         'type':'ligand'}

## protein side chain perturbations
#pert[0]={'subs':['nat','ile','leu','val'],
#         'segid':'PROC',
#         'resid':'3',
#         'type':'side_chain'}



##############################################
# Read in toppar files, coordinate files, etc.

# toppar files
settings.set_bomb_level(-2)
#pycharmm.lingo.charmm_script('stream ./variables'+str(nrun)+'.inp')
pycharmm.lingo.charmm_script('stream '+builddir+'/toppar.str')

# load pert toppar files
load_alchem_toppar(nsites,nsubs,pert,builddir)
settings.set_bomb_level(-2)

# # read in protein
# read.sequence_pdb(builddir+'/proc.pdb')
# gen.new_segment('PROC','NTER','CTER', setup_ic=True)   # (SEGID, FIRST, LAST, [options])
# read.pdb(builddir+'/proc.pdb',resid=True)
# 
# pycharmm.lingo.charmm_script('''
# ic generate
# ic param
# ic build
# hbuild sele hydrogen end
# auto angle dihe''')

# read in ligand core (written by msld-py-prep)
read.sequence_pdb(builddir+'/core.pdb')
gen.new_segment(pert[0]['segid'], setup_ic=True)   # (SEGID, FIRST, LAST, [options])
read.pdb(builddir+'/core.pdb',resid=True)


# load in alchemical patches 
load_alchem_patches(nsites,nsubs,pert,builddir)


# delete angles and dihedrals between alchem groups
pycharmm.lingo.charmm_script('auto angle dihe')
settings.set_bomb_level(-1)
for site in range(nsites):
    for sub1 in range(nsubs[site]):
        for sub2 in range(sub1+1,nsubs[site]):
            pycharmm.lingo.charmm_script('dele connectivity sele {} show end sele {} show end'
            .format(pert[site]['select'][sub1],pert[site]['select'][sub2]))


# read in solvent and ions
read.sequence_pdb(builddir+'/solvent.pdb')
gen.new_segment('WT00', setup_ic=True, angle=False, dihedral=False)
read.pdb(builddir+'/solvent.pdb',resid=True)

#read.sequence_pdb(builddir+'/ions.pdb')
#gen.new_segment('IONS', setup_ic=True, angle=False, dihedral=False)
#read.pdb(builddir+'/ions.pdb',resid=True)

pycharmm.lingo.charmm_script('print coor sele .not. init end')
settings.set_bomb_level(0)

# # write out psf, crd, pdb files
# write.psf_card('patch.psf')
# write.coor_card('patch.crd')
# write.coor_pdb('patch.pdb')


##############################################
# Create water box & periodic images

crystal.define_cubic(box)
crystal.build(14.0)

# center at 0.0 for blade and regular charmm; center at boxhalf for omm!  # way to automate this section?
image.setup_segment(0.0, 0.0, 0.0, 'LIG' )  # list out one by one
image.setup_residue(0.0, 0.0, 0.0, 'WT00')
image.setup_residue(0.0, 0.0, 0.0, 'IONS')

pycharmm.lingo.charmm_script('coor copy comp') # for truncated system


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
pert[site]['buffer']=[]  # could also use a dict
ii=2
for site in range(nsites):
    for sub in range(nsubs[site]):
        block_call+='Call {} sele {} show end\n'.format(ii,pert[site]['select'][sub])
        pert[site]['buffer'].append(ii)
        ii+=1
# for softcore atoms (not sure how to loop-create CATS atoms...?)
block_parm='''
! scat on
! scat k {}
! cats sele atom ?segid ?resid ?atomname .or. [list of atom names to cat]

qldm theta
lang temp {}
soft w14
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
        # block_ldin+='ldin {} {:.4f} 0.0 5.0 @lams{}s{} 5.0\n'.format(pert[site]['buffer'][sub],tmplmb,site+1,sub+1)
        bii=pert[site]['buffer'][sub]
        iii=sub0[site]+sub
        block_ldin+='ldin {} {:.4f} 0.0 5.0 {} 5.0\n'.format(bii,tmplmb,bias['b'][0,iii])
        sitestr+=str(site+1)+'  '
# add in exclusions with adex
block_adex=''
for site in range(nsites):
    for sub1 in range(nsubs[site]):
        for sub2 in range(sub1+1,nsubs[site]):
            block_adex+='adex {} {}\n'.format(pert[site]['buffer'][sub1],pert[site]['buffer'][sub2])

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
                    bii=pert[si]['buffer'][ii]
                    bjj=pert[sj]['buffer'][jj]
                    iii=sub0[si]+ii
                    jjj=sub0[sj]+jj
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


