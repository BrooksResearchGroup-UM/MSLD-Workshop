#! /usr/bin/env python

import sys
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

def main(fnm):

  ignoreResname=['TIP3']

  fp=open(fnm,"r")
  crd=crdio.readcrd(fp)
  fp.close()

  N=0
  CoM=np.zeros((1,3))

  xyz=np.zeros((len(crd['data']),3))
  for atom in crd['data']:
    if not atom['resname'] in ignoreResname:
      xyz[N,:]=atom['xyz']
      CoM+=xyz[N,:]
      N+=1
  CoM/=N
  xyz=xyz[0:N,:]

  maximum=np.tile(CoM,[3,1])
  minimum=np.tile(CoM,[3,1])
  for i in range(N):
    for j in range(3):
      if xyz[i,j]>maximum[j,j]:
        maximum[j,:]=xyz[i,:]
      if xyz[i,j]<minimum[j,j]:
        minimum[j,:]=xyz[i,:]
  rmax=2*np.max(np.concatenate((np.diagonal(maximum)-CoM,CoM-np.diagonal(minimum)),axis=0))
  sys.stderr.write("Old bad estimate %f\n" % (rmax,))
  rmax=np.max(maximum[0::4]-minimum[0::4])
  sys.stderr.write("1st under estimate %f\n" % (rmax,))
  dr=maximum-minimum
  drmag=np.sqrt(np.sum(dr*dr,axis=1))
  rmax=np.max(drmag)
  sys.stderr.write("2nd under estimate %f\n" % (rmax,))

  for i in range(3):
    if drmag[i]==rmax:
      j=i
  u=dr[j,:]/drmag[j]

  umaximum=CoM
  uminimum=CoM
  dotmaximum=np.dot(umaximum,u)
  dotminimum=np.dot(uminimum,u)
  for i in range(N):
    d=np.dot(xyz[i,:],u)
    if d>dotmaximum:
      # print("max",i,xyz[i,:])
      dotmaximum=d
      umaximum=xyz[i,:]
    if d<dotminimum:
      # print("min",i,xyz[i,:])
      dotminimum=d
      uminimum=xyz[i,:]
  rmax=np.max(np.sqrt(np.sum((umaximum-uminimum)*(umaximum-uminimum))))
  sys.stderr.write("3rd under estimate %f\n" % (rmax,))

  # Any distance longer than the current rmax must have at least one end outside a sphere of radius rmax/2
  # Find those ends, and check them against all other atoms
  umid=0.5*(umaximum+uminimum)
  targetDist2=xyz[0:N,:]-np.tile(umid,[N,1])
  targetDist2=np.sum(targetDist2*targetDist2,axis=1)
  rout=0.5*rmax
  rin=rmax-np.sqrt(np.max(targetDist2))
  targetout=(targetDist2>0.999999*rout*rout)
  # print(np.argwhere(targetout))
  xyzout=xyz[targetout,:]
  targetin=(targetDist2>0.999999*rin*rin)
  # print(np.argwhere(targetin))
  xyzin=xyz[targetin,:]
  sys.stderr.write("Found comparing %d outer targets to %d inner targets\n" % (xyzout.shape[0],xyzin.shape[0]))
  rmax2=rmax*rmax
  for i in range(xyzout.shape[0]):
    for j in range(xyzin.shape[0]):
      dx=xyzout[i,:]-xyzin[j,:]
      if np.sum(dx*dx)>rmax2:
        rmax2=np.sum(dx*dx)
  rmax=np.sqrt(rmax2)
  sys.stderr.write("Exact answer %f\n" % (rmax,))

  return rmax

if __name__ == '__main__':
  if len(sys.argv) != 2:
    print("Error, expects one argument")
    quit()
  rmax=main(sys.argv[1])
  sys.stdout.write("%f\n" % (rmax,))
