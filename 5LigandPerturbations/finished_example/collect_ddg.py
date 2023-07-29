#! /usr/bin/env python

##
## Collect results from protein/analysis[###]/Result.txt and water/analysis[###]/Result.txt
## to compute a d.d.G(binding) for MSLD ligand perturbations
##
## Usage: python collect_ddg.py [###]
##

import sys
import numpy as np

# get run number from cmd ln
if len(sys.argv) != 2:
    print("Supply run number as 1st cmd ln argu")
    quit()
else:
    runnum=sys.argv[1]

# load protein data
p_ddg = np.loadtxt('./protein/analysis'+runnum+'/Result.txt', usecols=1,dtype=float)
p_err = np.loadtxt('./protein/analysis'+runnum+'/Result.txt', usecols=3,dtype=float)

# load water data
w_ddg = np.loadtxt('./water/analysis'+runnum+'/Result.txt', usecols=1,dtype=float)
w_err = np.loadtxt('./water/analysis'+runnum+'/Result.txt', usecols=3,dtype=float)

# calc ddG
ddG=p_ddg-w_ddg

# cald err
err=np.sqrt(np.square(p_err)+np.square(w_err))


## print everything out
print("%11s %11s %11s %11s %11s %11s" % ("Protein dG","Protein err","Water dG","Water err","ddG(bind)","err"))
for row in range(ddG.shape[0]):
    print("%11.3f %11.3f %11.3f %11.3f %11.3f %11.3f" % (p_ddg[row],p_err[row],w_ddg[row],w_err[row],ddG[row],err[row]))


