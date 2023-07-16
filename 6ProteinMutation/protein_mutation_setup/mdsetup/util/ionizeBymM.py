#! /usr/bin/env python

from __future__ import print_function

import sys
import numpy as np
import copy
import random

import crdio

solventSegid='WT01'
ionSegid='IONS'
cutoff=10
solvcut=4.5
ionCharge={'SOD':1,'POT':1,'CLA':-1}
H2OmM=55509.

if len(sys.argv) != 4:
  sys.stderr.write("Error, expects one argument for crd file followed by one argument for charge followed by one argument with requested ion concentrations (mM), e.g. POT:100.0=CLA:100.0\n")
  quit()

fp=open(sys.argv[1],"r")
system=crdio.readcrd(fp)
fp.close()

charge=int(round(float(sys.argv[2])))

ionRequestmM={}
for token in sys.argv[3].split('='):
  ionRequestmM[token.split(':')[0]]=float(token.split(':')[1])
for key in ionRequestmM:
  sys.stderr.write("Request %f mM %s ions\n" % (ionRequestmM[key],key))

ionTile={}
fp=open('util/ions/sod.crd',"r")
ionTile['SOD']=crdio.readcrd(fp)
fp.close()
fp=open('util/ions/pot.crd',"r")
ionTile['POT']=crdio.readcrd(fp)
fp.close()
fp=open('util/ions/cla.crd',"r")
ionTile['CLA']=crdio.readcrd(fp)
fp.close()

abc=np.loadtxt('abc.dat')
abcd=np.diagonal(abc)

sys.stderr.write("Setting up grid search\n")
i=0
lookup={}
residSolvent=[]
lookupSolvent={}
for atom in system['data']:
  if atom['segid']!=solventSegid:
    # Consider other images, even though they probably don't matter...
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
            lookup[id3].append(i)
    # Or if you're SURE they don't matter:
    # id3=tuple(np.floor(atom['xyz']/solvcut).squeeze().astype(int))
    # if not id3 in lookup:
    #   lookup[id3]=[]
    # lookup[id3].append(i)
  else:
    if not atom['resid'] in residSolvent:
      residSolvent.append(atom['resid'])
      lookupSolvent[atom['resid']]=[]
    lookupSolvent[atom['resid']].append(i)
  i+=1
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

# Convert ionRequestmM to ionRequest
C=np.zeros((len(ionRequestmM),len(ionRequestmM)))
V=np.zeros((len(ionRequestmM),))
c=np.zeros((len(ionRequestmM),))
q=np.zeros((len(ionRequestmM),))
NH2O0=len(residSolvent)
gammaNorm=0
ionList=ionRequestmM.keys() # list in python2, view in python3, can't index view
i=0 # manually index ions
for ion in ionList: # use ion instead of ionList[i]
  c[i]=ionRequestmM[ion]/H2OmM
  q[i]=ionCharge[ion]
  gammaNorm+=q[i]*q[i]*c[i]
  i+=1
for i in range(len(ionList)):
  C[i,i]+=1
  V[i]+=NH2O0*c[i]
  V[i]+=-q[i]*c[i]*charge/gammaNorm
  for j in range(len(ionList)):
    C[i,j]+=c[i]
Nfloat=np.matmul(np.linalg.inv(C),V)
Napprox=np.floor(Nfloat+0.5).astype(int)
Qapprox=charge+np.dot(q,Napprox).astype(int)
print(Nfloat,file=sys.stderr)
print(Napprox,file=sys.stderr)
print(Qapprox,file=sys.stderr)
Nfinal=Napprox+0
if True: # not Qapprox==0:
  iFix=-1
  for i in range(len(ionList)):
    if abs(q[i])==1:
      if iFix==-1:
        iFix=i
      elif abs(Nfloat[i]-Napprox[i])>abs(Nfloat[i]-Napprox[i]):
        iFix=i
  Nfinal[iFix]-=Qapprox*q[iFix].astype(int)
print(Nfinal,file=sys.stderr)
ionRequest={}
i=0 # manually index ions again
for ion in ionList: # use ion instead of ionList[i]
  ionRequest[ion]=Nfinal[i]
  i+=1
for key in ionRequest:
  sys.stderr.write("Request %f %s ions\n" % (ionRequest[key],key))

resid=1
solvcut2=solvcut*solvcut
random.seed(2401)
residSolventReplace=[]
ions={'data':[]}
for ionResname in ionRequest:
  sys.stderr.write("Place %s ions\n" % (ionResname,))
  i=0
  attempt=0
  while i<ionRequest[ionResname]:
    if attempt==1000:
      sys.stderr.write("Error: could not place a single ion in 1000 attempts\n")
      quit()
    # else:
    #   sys.stderr.write("Try ion %d attempt %d\n" % (i,attempt))
    attempt+=1
    target=residSolvent[random.randrange(len(residSolvent))]
    # If we haven't already replaced it
    if not target in residSolventReplace:
        # See if it's too close to the protein
        xSolvent=system['data'][lookupSolvent[target][0]]['xyz']
        id3=tuple(np.floor(xSolvent/solvcut).squeeze().astype(int))
        include=True
        # Check if it's too close to the protein
        if include and id3 in lookupNeighbor:
          for soluteId in lookupNeighbor[id3]:
            xSolute=system['data'][soluteId]['xyz']
            if np.sum((xSolute-xSolvent)*(xSolute-xSolvent))<solvcut2:
              include=False
              break
        if include:
          for atom in ionTile[ionResname]['data']:
            newAtom=copy.deepcopy(atom)
            newAtom['resid']=("%d" % (resid,))
            newAtom['xyz']+=xSolvent
            newAtom['segid']=ionSegid
            ions['data'].append(newAtom)
          resid+=1
          i+=1
          attempt=0
          residSolventReplace.append(target)
sys.stderr.write("Finished replacing waters with ions\n")

Q=0
for ionResname in ionRequest:
  N=ionRequest[ionResname]
  C=(N*55.8*1000)/(len(residSolvent)-len(residSolventReplace))
  sys.stderr.write("%s ion concentration %f mM\n" % (ionResname,C))
  Q+=ionCharge[ionResname]*N
sys.stderr.write("Total ion charge is %d (should cancel protein charge)\n" % (Q,))

newSystem={}
newSystem['title']=copy.deepcopy(system['title'])
newSystem['title'].insert(-1,'* crd ions added by ionize.py\n')
newSystem['data']=[]
# Add solute atoms
for atom in system['data']:
  if atom['segid']!=solventSegid:
    newSystem['data'].append(copy.deepcopy(atom))
# Add unreplaced water ions
newResid=1
for resid in residSolvent:
  if not resid in residSolventReplace:
    for i in lookupSolvent[resid]:
      newAtom=copy.deepcopy(system['data'][i])
      newAtom['resid']=("%d" % (newResid,))
      newSystem['data'].append(newAtom)
    newResid+=1
# Add ions
for atom in ions['data']:
  newSystem['data'].append(atom)

crdio.writecrd(sys.stdout,newSystem)
