#! /usr/bin/env python

import os, sys, shutil, subprocess
import alf

gotcwd=os.getcwd()

chmguidir=sys.argv[1]
setupdir=gotcwd+'/charmm-gui-setup'
tooldir=gotcwd+'/charmm-gui-tools'

# md equilibration

if not os.path.exists(setupdir):
  os.mkdir(setupdir)
if os.path.exists(setupdir+'/prep'):
  shutil.rmtree(setupdir+'/prep')
os.mkdir(setupdir+'/prep')

shutil.copy(chmguidir+'/step3_pbcsetup.psf',setupdir+'/prep/')
shutil.copy(chmguidir+'/step3_pbcsetup.crd',setupdir+'/prep/')
shutil.copy(chmguidir+'/step3_pbcsetup.str',setupdir+'/prep/')
shutil.copy(chmguidir+'/crystal_image.str',setupdir+'/prep/')

shutil.copy(chmguidir+'/toppar.str',setupdir+'/prep/')
if os.path.exists(setupdir+'/prep/toppar'):
  shutil.rmtree(setupdir+'/prep/toppar')
shutil.copytree(chmguidir+'/toppar',setupdir+'/prep/toppar')

# Fix TIP3 rtf
fpin=open(chmguidir+'/toppar/toppar_water_ions.str','r')
fpout=open(setupdir+'/prep/toppar/toppar_water_ions.str','w')
for line in fpin:
  if len(line.split())>=2 and line.split()[0]=='RESI' and line.split()[1]=='TIP3':
    fpout.write(' '.join(line.split()[0:3]+['NOANG','NODIH']+line.split()[3:])+'\n')
  else:
    fpout.write(line)
fpin.close()
fpout.close()

if not os.path.exists(tooldir):
  print('Error missing tools directory')

shutil.copy(tooldir+'/nbond.str',setupdir+'/prep/')

shutil.copy(tooldir+'/eqchmgui.inp',setupdir+'/')

CHARMM=os.environ['CHARMMEXEC']

os.chdir(setupdir)
# Set OMP_NUM_THREADS=1 for BLaDE
subprocess.call(['mpirun','-np','1','-x','OMP_NUM_THREADS=1',CHARMM,'-i','eqchmgui.inp'])
os.chdir(gotcwd)

# msld setup

shutil.copy(tooldir+'/generic.inp',setupdir+'/prep/l33m.inp')

alf_info_str="""
import numpy as np
import os
alf_info={}
alf_info['name']='l33m'
alf_info['nsubs']=[2]
alf_info['nblocks']=np.sum(alf_info['nsubs'])
alf_info['ncentral']=0
alf_info['nreps']=1
alf_info['nnodes']=1
alf_info['enginepath']=os.environ['CHARMMEXEC']
alf_info['temp']=298.15
"""
fp=open(setupdir+'/prep/alf_info.py','w')
fp.write(alf_info_str)
fp.close()

alchemical_definitions_str="""
! List all mutation sites and mutants at each site
! j is hsp
set resid1 = 33
set s1seq1 = 0 ! l
set s1seq2 = m
set segid1 = PROA

! Set the terminal properties of any segid mutated above
set nterdel_proa = 0 ! 0 means don't do it
set nterres_proa = 31
set ntercap_proa = ace
set nterc_proa = 4 ! 2 nter, 3 cter, 4 ace, 5 ct3
set cterdel_proa = 0 ! 0 means don't do it
set cterres_proa = 35
set ctercap_proa = ct3
set cterc_proa = 5  ! 2 nter, 3 cter, 4 ace, 5 ct3

! Only modify aainitl and aafinal if you're mutating things besides proteins
set aainitl = 0
set aafinal = @nsites
"""
fp=open(setupdir+'/prep/alchemical_definitions.inp','w')
fp.write(alchemical_definitions_str)
fp.close()

if os.path.exists(setupdir+'/prep/aa_stream'):
  shutil.rmtree(setupdir+'/prep/aa_stream')
shutil.copytree('aa_stream',setupdir+'/prep/aa_stream')

os.chdir(setupdir)
sys.path.insert(0,'') # so alf can find prep after os.chdir
alf.initialize(engine='bladelib')
alf.runflat(1,2,13000,39000,engine='bladelib')
