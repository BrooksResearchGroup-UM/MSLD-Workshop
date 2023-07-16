#! /usr/bin/env python

import sys
import numpy as np
import copy

import crdio
import maxdimension

# Charmm coor format
#         title
#         NATOM (I10)
#         ATOMNO RESNO   RES  TYPE  X     Y     Z   SEGID RESID Weighting
#           I10   I10 2X A8 2X A8       3F20.10     2X A8 2X A8 F20.10

# PDB format
# text IATOM  TYPE  RES  IRES      X  Y  Z    W
#  A6   I5  2X A4   A4    I5  4X     3F8.3 6X F6.2

solventSegid='WT01'
cutoff=10
solvcut=2.0

if len(sys.argv) < 3:
  sys.stderr.write("Error, expects at least two arguments\n")
  quit()

if len(sys.argv) < 4:
  rmax=maxdimension.main(sys.argv[1])
  boxL=rmax+2*cutoff
else:
  boxL=float(sys.argv[3])
  rmax=boxL-2*cutoff

if len(sys.argv) > 4:
  sys.stderr.write("Error, only uses three arguments\n")
  quit()

crystalType=sys.argv[2]
if not crystalType in ['rhdo','octa','cubi']:
  sys.stderr.write("Error, only rhdo, octa, and cubi are recognized box types for second argument")
  quit()

abcabc=np.zeros((6,))
abcabc[0:3]=boxL
if crystalType=='rhdo':
  abcabc[3]=np.arccos(1./2.)
  abcabc[4]=np.arccos(0.)
  abcabc[5]=np.arccos(1./2.)
elif crystalType=='octa':
  abcabc[3:6]=np.arccos(-1./3.)
elif crystalType=='cubi':
  abcabc[3:6]=np.arccos(0.)
np.savetxt('abcabc.dat',abcabc,fmt="%20.15f")

abc=np.zeros((3,3))
abc[0,0]=abcabc[0]
abc[1,0]=abcabc[1]*np.cos(abcabc[5])
abc[1,1]=abcabc[1]*np.sin(abcabc[5])
abc[2,0]=abcabc[2]*np.cos(abcabc[4])
abc[2,1]=abcabc[2]*(np.cos(abcabc[3])-np.cos(abcabc[4])*np.cos(abcabc[5]))/np.sin(abcabc[5])
abc[2,2]=np.sqrt(abcabc[2]*abcabc[2]-abc[2,0]*abc[2,0]-abc[2,1]*abc[2,1])
np.savetxt('abc.dat',abc,fmt="%20.15f")


# hth=np.matmul(abc,np.transpose(abc))
# [hthw,hthv]=np.linalg.eigh(hth)
# hthwrt=np.diag(np.sqrt(hthw))
# abc_charmm=np.matmul(np.matmul(hthv,hthwrt),np.transpose(hthv))
# 
# kabc=np.zeros((3,3))
# Vabc=np.dot(np.cross(abc[0,:],abc[1,:]),abc[2,:])
# kabc[0,:]=np.cross(abc[1,:],abc[2,:])/Vabc
# kabc[1,:]=np.cross(abc[2,:],abc[0,:])/Vabc
# kabc[2,:]=np.cross(abc[0,:],abc[1,:])/Vabc
# 
# kabc_charmm=np.zeros((3,3))
# kabc_charmm[0,:]=np.cross(abc_charmm[1,:],abc_charmm[2,:])/Vabc
# kabc_charmm[1,:]=np.cross(abc_charmm[2,:],abc_charmm[0,:])/Vabc
# kabc_charmm[2,:]=np.cross(abc_charmm[0,:],abc_charmm[1,:])/Vabc
# 
# Rtocharmm=np.matmul(np.transpose(abc_charmm),kabc)

abcd=np.diagonal(abc)

fp=open(sys.argv[1],"r")
solute=crdio.readcrd(fp)
fp.close()
fp=open('util/convpdb_waterbox/waterbox.crd',"r")
solventTile=crdio.readcrd(fp)
fp.close()
solventTileL=np.loadtxt('util/convpdb_waterbox/waterbox.L')

maximum=np.zeros((1,3))
minimum=np.zeros((1,3))
for i in range(3):
  maximum[0,i]=solute['data'][0]['xyz'][0,i]
  minimum[0,i]=solute['data'][0]['xyz'][0,i]
for atom in solute['data']:
  for i in range(3):
    if atom['xyz'][0,i]>maximum[0,i]:
      maximum[0,i]=atom['xyz'][0,i]
    if atom['xyz'][0,i]<minimum[0,i]:
      minimum[0,i]=atom['xyz'][0,i]

shift=-0.5*(maximum+minimum)
# print(solute['data'][0]['xyz'])
for atom in solute['data']:
  atom['xyz']+=shift
# print(solute['data'][0]['xyz'])

sys.stderr.write("Setting up grid search\n")
lookup={}
for atom in solute['data']:
  posz=atom['xyz']
  zimage0=-np.floor((posz[0,2]+solvcut)/abcd[2]+0.5)
  for zimagei in range(2):
    posy=posz+(zimagei+zimage0)*abc[2,:]
    yimage0=-np.floor((posy[0,1]+solvcut)/abcd[1]+0.5)
    for yimagei in range(2):
      posx=posy+(yimagei+yimage0)*abc[1,:]
      ximage0=-np.floor((posx[0,0]+solvcut)/abcd[0]+0.5)
      for ximagei in range(2):
        pos=posx+(ximagei+ximage0)*abc[0,:]
        if np.all((np.abs(pos)+solvcut)<(0.5*abcd)):
          id3=tuple(np.floor(pos/solvcut).squeeze().astype(int))
          if not id3 in lookup:
            lookup[id3]=[]
          lookup[id3].append(pos)
lookupNeighbor={}
for key in lookup:
  for i in range(3):
    for j in range(3):
      for k in range(3):
        keyNeighbor=(key[0]+i-1,key[1]+j-1,key[2]+k-1)
        if not keyNeighbor in lookupNeighbor:
          lookupNeighbor[keyNeighbor]=[]
        lookupNeighbor[keyNeighbor].extend(lookup[key])
sys.stderr.write("Done setting up grid search\n")

sys.stderr.write("Setting up tiling structure\n")
solventTileResidue=[]
for atom in solventTile['data']:
  resid=int(atom['resid'])
  if resid>len(solventTileResidue):
    solventTileResidue.append([])
  solventTileResidue[resid-1].append(atom)
sys.stderr.write("Done setting up tiling structure\n")

resid=1
solvcut2=solvcut*solvcut
dx=np.zeros((1,3))
halfTileCount=np.floor(0.5*abcd/solventTileL+0.5).astype(int) # (inclusive)
for i in range(-halfTileCount[0],halfTileCount[0]+1):
  dx[0,0]=solventTileL*i
  sys.stderr.write("Tile %d * *\n" % (i,))
  for j in range(-halfTileCount[1],halfTileCount[1]+1):
    dx[0,1]=solventTileL*j
    for k in range(-halfTileCount[2],halfTileCount[2]+1):
      dx[0,2]=solventTileL*k
      for residue in solventTileResidue:
        xSolvent=residue[0]['xyz']+dx
        id3=tuple(np.floor(xSolvent/solvcut).squeeze().astype(int))
        # Check if it's in the box
        include=np.all(np.abs(xSolvent)<0.5*abcd)
        # Check if it's too close to the protein
        if include and id3 in lookupNeighbor:
          for xSolute in lookupNeighbor[id3]:
            if np.sum((xSolute-xSolvent)*(xSolute-xSolvent))<solvcut2:
              include=False
              break
        if include:
          for atom in residue:
            newAtom=copy.deepcopy(atom)
            newAtom['resid']=("%d" % (resid,))
            newAtom['xyz']+=dx
            newAtom['segid']=solventSegid
            solute['data'].append(newAtom)
          resid+=1

solute['title'].insert(-1,'* crd solvated by solvate.py\n')

crdio.writecrd(sys.stdout,solute)
