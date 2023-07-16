#! /usr/bin/env python

import sys
import numpy as np
import copy

import crdio

# Charmm coor format
#         title
#         NATOM (I10)
#         ATOMNO RESNO   RES  TYPE  X     Y     Z   SEGID RESID Weighting
#           I10   I10 2X A8 2X A8       3F20.10     2X A8 2X A8 F20.10

# PDB format
# text IATOM  TYPE  RES  IRES      X  Y  Z    W
#  A6   I5  2X A4   A4    I5  4X     3F8.3 6X F6.2

if len(sys.argv) != 2:
  sys.stderr.write("Error, expects one argument\n")
  quit()

abc=np.loadtxt('abc.dat')

hth=np.matmul(abc,np.transpose(abc))
[hthw,hthv]=np.linalg.eigh(hth)
hthwrt=np.diag(np.sqrt(hthw))
abc_charmm=np.matmul(np.matmul(hthv,hthwrt),np.transpose(hthv))

kabc=np.zeros((3,3))
Vabc=np.dot(np.cross(abc[0,:],abc[1,:]),abc[2,:])
kabc[0,:]=np.cross(abc[1,:],abc[2,:])/Vabc
kabc[1,:]=np.cross(abc[2,:],abc[0,:])/Vabc
kabc[2,:]=np.cross(abc[0,:],abc[1,:])/Vabc

kabc_charmm=np.zeros((3,3))
kabc_charmm[0,:]=np.cross(abc_charmm[1,:],abc_charmm[2,:])/Vabc
kabc_charmm[1,:]=np.cross(abc_charmm[2,:],abc_charmm[0,:])/Vabc
kabc_charmm[2,:]=np.cross(abc_charmm[0,:],abc_charmm[1,:])/Vabc

Rtocharmm=np.matmul(np.transpose(abc_charmm),kabc)

fp=open(sys.argv[1],"r")
solute=crdio.readcrd(fp)
fp.close()

for atom in solute['data']:
  pos=atom['xyz']
  posnew=np.transpose(np.matmul(Rtocharmm,np.transpose(pos)))
  atom['xyz']=posnew

solute['title'].insert(-1,'* crd rotated to charmm frame by rotate.py\n')

crdio.writecrd(sys.stdout,solute)
