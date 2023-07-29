#! /usr/bin/env python

## list out all the segids in a pdb file (such as after 
## solvating a system with CHARMM-GUI)

import sys,os

if len(sys.argv) != 2:
    print("Supply name of file as first line argument")
    quit()
else:
    filename=sys.argv[1]

print("Split",filename,"into separate segids")

segids=[]
for line in open(filename,'r'):
    if line[0:4] == 'ATOM':
        tmp=line.split()
        if not (tmp[-1] in segids):
            segids.append(tmp[-1])

print(segids)

## extract segid lines into separate files
fp=open(filename,'r')
for s in segids:
    fp.seek(0) # allows the order of s in segids to be nonchronological
    line=fp.readline()
    output=open(s.lower()+'.pdb','w')
    while line:
        if line[0:4] == 'ATOM':
            tmp=line.split()
            if tmp[-1] == s:
                output.write(line)
            else:
                pass
        line=fp.readline()
    output.write('TER\n')
    output.write('END\n')
    output.close()
fp.close()
    

