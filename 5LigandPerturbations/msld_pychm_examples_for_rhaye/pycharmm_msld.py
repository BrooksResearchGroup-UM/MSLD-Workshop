##
## pyCHARMM script drafted from many examples
##
##   ** MSLD with BLADE (no OMM support) **
##

import os
import sys
import numpy as np
import pandas
from msld_patch import *

##############################################
# Load pyCHARMM libraries

import pycharmm
import pycharmm.generate as gen
import pycharmm.ic as ic
import pycharmm.coor as coor
import pycharmm.energy as energy
import pycharmm.dynamics as dyn
import pycharmm.nbonds as nbonds
import pycharmm.minimize as minimize
import pycharmm.crystal as crystal
import pycharmm.image as image
import pycharmm.psf as psf
import pycharmm.read as read
import pycharmm.write as write
import pycharmm.settings as settings
import pycharmm.cons_harm as cons_harm
import pycharmm.cons_fix as cons_fix
import pycharmm.select as select
import pycharmm.shake as shake
import pycharmm.scalar as scalar
from pycharmm.lib import charmm as libcharmm


##############################################
# Set up global parameters

# variables
sysname = 'water'
box = 43.0
pmegrid = 40.0
temp = 298.15
builddir = './prep'

# nonbonded conditions
nb_fswitch = False          # normal fswitching functions
nb_pme = True               # normal PME

# dynamics conditions
blade = True

# dynamics variables
cpt_on   = True             # run with CPT for NPT?
timestep = 0.002            # ps
ns1      = 500000           # number of MD steps per 1 ns
total_ns = 1                # total number of production ns sampling
nequil = int(ns1*(1/10))    # equil for 100 ps
nprod  = ns1*total_ns       # prod sampling for 5 ns
nsavc  = 1000               # dcd save frequency

# msld variables
fnex = 5.5
nrun = 116                  # msld run number (could also be cmd ln argu)

nblocks=np.loadtxt('./nblocks',dtype='int')
nsubs=np.loadtxt('./nsubs',dtype='int',ndmin=1)
nsites=len(nsubs)

# dictionary to define perturbations
# types = 'ligand','side_chain'
pert={}                    # dict of nested dicts; index = site # (zero indexed)

# ligand perturbations
pert[0]={'subs':['1','3','4','5','6','7'],  # site1_sub# (0-index = site 1)
         'segid':'LIG',
         'resid':'1',
         'type':'ligand'}

# protein side chain perturbations
pert[0]={'subs':['nat','ile','leu','val'],
         'segid':'PROC',
         'resid':'3',
         'type':'side_chain'}



##############################################
# Read in toppar files, coordinate files, etc.

# toppar files
settings.set_bomb_level(-2)
pycharmm.lingo.charmm_script('stream ./variables'+str(nrun)+'.inp')
pycharmm.lingo.charmm_script('stream '+builddir+'/toppar.str')

# load pert toppar files
load_alchem_toppar(nsite,nsubs,pert,builddir)
settings.set_bomb_level(-2)

# read in protein
read.sequence_pdb(builddir+'/proc.pdb')
gen.new_segment('PROC','NTER','CTER', setup_ic=True)   # (SEGID, FIRST, LAST, [options])
read.pdb(builddir+'/proc.pdb',resid=True)

pycharmm.lingo.charmm_script('''
ic generate
ic param
ic build
hbuild sele hydrogen end
auto angle dihe''')

# read in ligand core (written by msld-py-prep)
read.sequence_pdb(builddir+'/core.pdb')
gen.new_segment(lpert[0]['segid'], setup_ic=True)   # (SEGID, FIRST, LAST, [options])
read.pdb(builddir+'/core.pdb',resid=True)


# load in alchemical patches 
load_alchem_patches(nsite,nsubs,pert,builddir)


# delete angles and dihedrals between alchem groups
pycharmm.lingo.charmm_script('auto angle dihe')
settings.set_bomb_level(-1)
for site in range(nsites):
    for sub1 in range(nsubs[site]):
        for sub2 in range(sub1+1,nsubs[site]):
            pycharmm.lingo.charmm_script('dele connectivity sele {} show end sele {} show end'
            .format(pert[site]['select'][sub1],pert[site]['select'][sub2]))


# read in solvent and ions
read.sequence_pdb(builddir+'/solv.pdb')
gen.new_segment('SOLV', setup_ic=True, angle=False, dihedral=False)
read.pdb(builddir+'/solv.pdb',resid=True)

read.sequence_pdb(builddir+'/ions.pdb')
gen.new_segment('IONS', setup_ic=True, angle=False, dihedral=False)
read.pdb(builddir+'/ions.pdb',resid=True)

pycharmm.lingo.charmm_script('print coor sele .not. init end')
settings.set_bomb_level(0)

# write out psf, crd, pdb files
write.psf_card('patch.psf')
write.coor_card('patch.crd')
write.coor_pdb('patch.pdb')


##############################################
# Create water box & periodic images

crystal.define_cubic(box)
crystal.build(14.0)

# center at 0.0 for blade and regular charmm; center at boxhalf for omm!  # way to automate this section?
image.setup_segment(0.0, 0.0, 0.0, 'PROC')  # list out one by one
image.setup_segment(0.0, 0.0, 0.0, 'LIG' )  # list out one by one
image.setup_residue(0.0, 0.0, 0.0, 'SOLV')
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
    pycharmm.lingo.charmm_script('stop')


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
for site in range(nsites):
    for sub in range(nsubs[site]):
        if sub == 0:
            tmplmb=1.0-(0.01*(nsubs[site]-1))
        else:
            tmplmb=0.01
        block_ldin+='ldin {} {:.4f} 0.0 5.0 @lams{}s{} 5.0\n'.format(pert[site]['buffer'][sub],tmplmb,site+1,sub+1)
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
block_varb='''
! Check block.doc for functional form of these biasing potentials
calc nbiaspot = 5 * ( @nblocks * ( @nblocks - 1 ) ) / 2
ldbi @nbiaspot
set ibias = 1
set iblock = 0
set si = 1
label loop5_vb
if @si .le. @nsites then
   set jblock = @iblock
   set sj = @si
   label loop5b_vb
   if @sj .le. @nsites then
      set ii = 1
      label loop6_vb
      if @ii .le. @nsubs@@{si} then
         calc ip1 = @ii + 1 + @iblock
         set jj = 1
         if @si .eq. @sj then
            calc jj = @ii + 1
         endif
         label loop7_vb
         if @jj .le. @nsubs@@{sj} then
            calc jp1 = @jj + 1 + @jblock

            !! !! adds in an pairwise exclusion
            !! if @si .eq. @sj then
            !!    adex @ip1 @jp1
            !! endif

            set c_shift = 0.0
            set x_shift = 0.0
            set s_shift = 0.0
            !vbrex! these lines needed in vb.inp, not here
            !vbrex! if @si .eq. @sj then
            !vbrex!    calc c_shift = 2.0 * (@myrep - @ncentral)
            !vbrex!    calc s_shift = 0.5 * (@myrep - @ncentral)
            !vbrex! endif

            calc coeff = @cs@@{si}s@@{ii}s@@{sj}s@@{jj} + @{c_shift}
            ldbv @ibias @ip1 @jp1 6 0.0 @coeff 0
            calc ibias = @ibias + 1
            calc coeff = @xs@@{si}s@@{ii}s@@{sj}s@@{jj} + @{x_shift}
            ldbv @ibias @ip1 @jp1 10 -5.56 @coeff 0
            calc ibias = @ibias + 1
            calc coeff = @ss@@{si}s@@{ii}s@@{sj}s@@{jj} + @{s_shift}
            ldbv @ibias @ip1 @jp1 8 0.017 @coeff 0
            calc ibias = @ibias + 1
            calc coeff = @xs@@{sj}s@@{jj}s@@{si}s@@{ii} + @{x_shift}
            ldbv @ibias @jp1 @ip1 10 -5.56 @coeff 0
            calc ibias = @ibias + 1
            calc coeff = @ss@@{sj}s@@{jj}s@@{si}s@@{ii} + @{s_shift}
            ldbv @ibias @jp1 @ip1 8 0.017 @coeff 0
            calc ibias = @ibias + 1
            calc jj = @jj + 1
            goto loop7_vb
         endif
         calc ii = @ii + 1
         goto loop6_vb
      endif
      calc jblock = @jblock + @nsubs@@{sj}
      calc sj = @sj + 1
      goto loop5b_vb
   endif
   calc iblock = @iblock + @nsubs@@{si}
   calc si = @si + 1
   goto loop5_vb
endif
END
'''

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


##############################################
# Set NonBonded settings & SP energy calc

cutnb = 14.0
cutim = cutnb
ctofnb = 12.0
ctonnb = 10.0

## nbond switching
## use a dictionary so that it becomes easy to switch between w/ vs w/o PME
nbonds_dict = {'cutnb':cutnb,'cutim':cutim,
           'ctonnb':ctonnb,'ctofnb':ctofnb,
           'atom':True,'vatom':True,
           'cdie':True,'eps':1.0,
           'inbfrq':-1, 'imgfrq':-1}

if nb_pme:
    nbonds_dict['switch']=True
    nbonds_dict['vfswitch']=True
    nbonds_dict['ewald']=True
    nbonds_dict['pmewald']=True
    nbonds_dict['kappa']=0.32
    nbonds_dict['fftx']=pmegrid
    nbonds_dict['ffty']=pmegrid
    nbonds_dict['fftz']=pmegrid
    nbonds_dict['order']=6

elif nb_fswitch:
    nbonds_dict['fswitch']=True
    nbonds_dict['vfswitch']=True
    nbonds_dict['ewald']=False
    nbonds_dict['pmewald']=False

else: 
    print("NonBonded Parameter Error - both pme and switch are false")
    pycharmm.lingo.charmm_script('stop')

nbonds=pycharmm.NonBondedScript(**nbonds_dict)
nbonds.run()
energy.show()

##############################################
# Minimize the system

minimize.run_sd(nstep=250,nprint=50,step=0.005,tolenr=1e-3,tolgrd=1e-3)
energy.show()
#minimize.run_abnr(nstep=250,nprint=50,tolenr=1e-3,tolgrd=1e-3)
#energy.show()

# write out psf, crd, pdb files
write.psf_card('minimized.psf')
write.coor_card('minimized.crd')
write.coor_pdb('minimized.pdb')


##############################################
# Set up and run Dynamics

# dynamics conditions
if blade:
    useblade = 'prmc pref 1 iprs 100 prdv 100'
    gscal = 0.1
    ntrfrq=0
    leap = True
    openmm = False
else: 
    print("MSLD can only be run with BLADE - exiting...")
    pycharmm.lingo.charmm_script('stop')

# set shake
shake.on(bonh=True,fast=True,tol=1e-7)
dyn.set_fbetas(np.full((psf.get_natom()),gscal,dtype=float))

# initialize blade
pycharmm.lingo.charmm_script('energy blade')

# set up output directories
if not os.path.isdir('res'): os.system('mkdir res')
if not os.path.isdir('dcd'): os.system('mkdir dcd')

# set up dynamics dictionary of parameters
dynamics_dict = {'cpt':cpt_on,'leap':True,'langevin':True,
    'timestep':timestep,
    'nsavc':nsavc,
    'nsavl':10,  # frequency for saving lambda values in lamda-dynamics
    'nprint': 1000, # Frequency to write to output
    'iprfrq': 10000, # Frequency to calculate averages
    'isvfrq': 10000, # Frequency to save restart file
    'ntrfrq':ntrfrq,
    'firstt':temp,'finalt':temp,'tstruct':temp,'tbath':temp,
    'iasors': 1,'iasvel':1,'iscvel': 0,'iscale': 0,
    'ihtfrq':0,'ieqfrq':0,'ichecw': 0,
    'inbfrq':-1,'imgfrq':-1,'ihbfrq':0,'ilbfrq':0,
    'echeck': -1}


if cpt_on:
    dynamics_dict['pconstant'] = True
    dynamics_dict['pmass'] = psf.get_natom()*0.12
    dynamics_dict['pref'] = 1.0
    dynamics_dict['pgamma'] = 20.0
    dynamics_dict['hoover'] = True
    dynamics_dict['reft'] = temp
    dynamics_dict['tmass'] = 1000


if blade:
    dynamics_dict['omm'] = False
    dynamics_dict['blade'] = useblade


# MD equilibration
dcd_file = pycharmm.CharmmFile(file_name='dcd/{}_{}.dcd'.format(sysname,'heat'), 
               file_unit=1,formatted=False,read_only=False)
res_file = pycharmm.CharmmFile(file_name='res/{}_{}.res'.format(sysname,'heat'), 
               file_unit=2,formatted=True,read_only=False)
lam_file = pycharmm.CharmmFile(file_name='res/{}_{}.lam'.format(sysname,'heat'), 
               file_unit=3,formatted=False,read_only=False)

dynamics_dict['start']  = True
dynamics_dict['nstep']  = nequil
dynamics_dict['iunrea'] = -1
dynamics_dict['iunwri'] = res_file.file_unit
dynamics_dict['iuncrd'] = dcd_file.file_unit
dynamics_dict['iunldm'] = lam_file.file_unit

equil_dyn = pycharmm.DynamicsScript(**dynamics_dict)
equil_dyn.run()

dcd_file.close()
res_file.close()
lam_file.close()

write.coor_pdb('dcd/{}_fframe.{}.pdb'.format(sysname,'equil')) # write out final frame


# MD production
dcd_file = pycharmm.CharmmFile(file_name='dcd/{}_{}.dcd'.format(sysname,'prod'), 
               file_unit=1,formatted=False,read_only=False)
res_file = pycharmm.CharmmFile(file_name='res/{}_{}.res'.format(sysname,'prod'), 
               file_unit=2,formatted=True,read_only=False)
prv_rest = pycharmm.CharmmFile(file_name='res/{}_{}.res'.format(sysname,'heat'), 
               file_unit=4,formatted=True,read_only=False)
lam_file = pycharmm.CharmmFile(file_name='res/{}_{}.lam'.format(sysname,'prod'), 
               file_unit=3,formatted=False,read_only=False)

dynamics_dict['start']  = False
dynamics_dict['restart']= True
dynamics_dict['nstep']  = nprod
dynamics_dict['iunrea'] = prv_rest.file_unit
dynamics_dict['iunwri'] = res_file.file_unit
dynamics_dict['iuncrd'] = dcd_file.file_unit
dynamics_dict['iunldm'] = lam_file.file_unit

prod_dyn = pycharmm.DynamicsScript(**dynamics_dict)
prod_dyn.run()

dcd_file.close()
res_file.close()
lam_file.close()

write.coor_pdb('dcd/{}_fframe.{}.pdb'.format(sysname,'prod')) # write out final frame

#if openmm: pycharmm.lingo.charmm_script('omm clear')
#if blade: pycharmm.lingo.charmm_script('blade off')

# collect lambda statistics
proc_lam = pycharmm.CharmmFile(file_name='res/{}_{}.lam'.format(sysname,'prod'), 
           file_unit=33,formatted=False,read_only=False)
pycharmm.lingo.charmm_script('traj lamb print ctlo 0.95 cthi 0.99 first {} nunit {}'.format(proc_lam.file_unit,1))


##############################################
# FINISHED


#pycharmm.lingo.charmm_script('pref')
pycharmm.lingo.charmm_script('stop')


