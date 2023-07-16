#! /usr/bin/env python

import sys
import numpy as np

# Charmm coor format
#         title
#         NATOM (I10)
#         ATOMNO RESNO   RES  TYPE  X     Y     Z   SEGID RESID Weighting
#           I10   I10 2X A8 2X A8       3F20.10     2X A8 2X A8 F20.10

# PDB format
# text IATOM  TYPE  RES  IRES      X  Y  Z    W
#  A6   I5  2X A4   A4    I5  4X     3F8.3 6X F6.2

def readcrd(fp):
  N=-1
  crd={}
  crd['title']=[]
  crd['data']=[]
  for line in fp:
    if len(line)>=1 and line[0]=='*':
      crd['title'].append(line)
    elif N<0:
      N=int(line[0:10])
    else:
      # id=int(line[0:10])
      # resno=int(line[10:20])
      resname=line[22:30].strip()
      name=line[32:40].strip()
      x=float(line[40:60])
      y=float(line[60:80])
      z=float(line[80:100])
      segid=line[102:110].strip()
      resid=line[112:120].strip()
      # crd['data'].append({'id':id,'resno':resno,'resname':resname,'name':name,'xyz':np.array([[x,y,z]]),'segid':segid,'resid':resid})
      crd['data'].append({'resname':resname,'name':name,'xyz':np.array([[x,y,z]]),'segid':segid,'resid':resid})
  if N!=len(crd['data']):
    sys.stderr.write("Error, wrong number of atoms in crd file\n")
    quit()
  return crd

def writecrd(fp,crd):
  for line in crd['title']:
    fp.write(line)
  fp.write("%10d  EXT\n" % (len(crd['data']),))

  seq={}
  resno=0
  for id in range(len(crd['data'])):
    atom=crd['data'][id]
    if not atom['segid'] in seq:
      seq[atom['segid']]={}
    if not atom['resid'] in seq[atom['segid']]:
      resno+=1
      seq[atom['segid']][atom['resid']]=resno
    elif seq[atom['segid']][atom['resid']] != resno:
      sys.stderr.write("Error: found more atoms from an old residue after printing another residue. Perhaps you have redundant resid? Current id %d, previous resno %d, current resno\n" % (id+1,seq[atom['segid']][atom['resid']],resno))
      quit()
    fp.write("%10d%10d  %-8s  %-8s%20.10f%20.10f%20.10f  %-8s  %-8s%20.10f\n" % (id+1,resno,atom['resname'],atom['name'],atom['xyz'][0,0],atom['xyz'][0,1],atom['xyz'][0,2],atom['segid'],atom['resid'],0))
