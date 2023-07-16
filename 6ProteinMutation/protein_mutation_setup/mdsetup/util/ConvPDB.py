#! /usr/bin/env python

import sys
import math
import numpy as np

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
  print("Error, expects one argument")
  quit()

def cndistance(c,n):
  if len(c)==0:
    return -1
  else:
    d2=0
    for i in range(3):
      d2+=(c[i]-n[i])**2
    return math.sqrt(d2)

allowatom=['ALA','CYS','ASP','GLU','PHE','GLY','HIS','HSD','HSE','HSP','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
allowhetatm=['TIP3']

renameres={'HIS':'HSD','HOH':'TIP3'}
renameatom={'TIP3':{'O':'OH2'},'ILE':{'CD1':'CD'}}

resmap={'TIP3':'WT00'}
chainmap={'A':'PROT'}
renumber=['WT00']

crd={}
crd['title']=[]
crd['title'].append("* CHARMM coor file converted from pdb file "+sys.argv[1]+"\n")
crd['title'].append("*\n")

crd['data']=[]
segid=""
resid=""
cxyz=()
sequence={}

fp=open(sys.argv[1],"r")
for line in fp:
  if len(line)>=6 and (line[0:4]=="ATOM" or line[0:6]=="HETATM"):
    name=line[13:17].strip()
    resname=line[17:21].strip()
    chain=line[21:22]
    prevresid=resid
    resid=line[22:27].strip()
    x=float(line[30:38])
    y=float(line[38:46])
    z=float(line[46:54])
    # Rename residues as requested
    if resname in renameres:
      resname=renameres[resname]
    if resname in renameatom:
      if name in renameatom[resname]:
        name=renameatom[resname][name]
    # Check for disallowed residues
    if line[0:4]=="ATOM":
      if not resname in allowatom:
        sys.stderr.write("WARNING: unknown resname\n")
        sys.stderr.write(line)
        continue
    if line[0:6]=="HETATM":
      if not resname in allowhetatm:
        sys.stderr.write("WARNING: unknown resname\n")
        sys.stderr.write(line)
        continue
    # Check chain connectivity
    if name=="C":
      cxyz=(x,y,z)
    if name=="N":
      nxyz=(x,y,z)
      d=cndistance(cxyz,nxyz)
      if d>=0 and (d<1.2 or d>1.5):
        sys.stderr.write("WARNING: chain connectivity likely broken with distance %f\n" % (d,))
        sys.stderr.write(line)
    # Assign SEGID
    prevsegid=segid
    if resname in resmap:
      segid=resmap[resname]
    elif chain in chainmap:
      segid=chainmap[chain]
    else:
      sys.stderr.write("ERROR: Cannot identify segid\n")
      sys.stderr.write(line)
      quit()
    # Increment counts
    if (not prevsegid==segid) or (not prevresid==resid):
      if not segid in sequence:
        sequence[segid]=[]
      sequence[segid].append(resname)
    # Renumber
    if (segid in renumber):
      residout=("%d" % (len(sequence[segid]),))
    else:
      residout=resid
    # Save entry
    crd['data'].append({'resname':resname,'name':name,'xyz':np.array([[x,y,z]]),'segid':segid,'resid':residout})
  # Reset previous segid and resid to help identify new residue
  elif len(line)>3 and line[0:3]=="TER":
    segid=""
    resid=""
    cxyz=()
sys.stderr.write('Found %d atoms\n' % (len(crd['data']),))


# Write what we got
crdio.writecrd(sys.stdout,crd)
