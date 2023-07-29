#! /usr/bin/env python

####
#### Executable script to build MSLD ready ligand files
#### JZV 02/2019
####

import msld_chk
import msld_mcs
import msld_crn
import msld_prm
import msld_wrt
import glob

###
### You may want to execute this in steps to give you finer control 
### over how the final msld ready files are constructed
###
### All ligand structure files (mol2) and toppar files must be 
### available prior to running this script!
###
### Reflig is printed in mcsout & can be changed between MCS searching
### and performing charge renormalization
###


#####################################################################
## (1) Define System and File Variables

sysname = "BiCycle2"                           # name of future output files
molfile = "mol_list.txt"                  # list of mol2 file names

mcsout = 'MCS_for_MSLD.txt'               # MCS output filename
outdir = 'build.'+sysname                 # MSLD output directory

cgenff=True                               # Are CGenFF/ParamChem parameters being used?

inFrag=[[]]  # reflig core atoms to include in each fragment at each site (list of nsub lists)
inCore=[[]]  # reflig atoms (usually anchor atoms) to include in the core (list of nsub lists)

#####################################################################
if len(glob.glob(mcsout)) == 0:
    ## (2) Check molfile and toppar files before getting started
    msld_chk.MsldCHK(molfile)
    print("chk finished")
    
    #####################################################################
    ## (3) Maximum Common SubStruct Search with bonded-environments
    ## "mcsout" = results of the search - edit this file to manual edit the splicing
    ## cutoff = RMSD & distance cutoff to differentiate different atoms
    ## change verbose to True to get more data printed to stdout
    
    reflig = msld_mcs.MsldMCS(molfile,mcsout,cutoff=0.8,debug=False)
    print("MCS results printed to "+mcsout)
    print("Reference Ligand is "+reflig)
    quit()


#####################################################################
## (4) Perform Charge-Renormalization 
## To manually change what atoms are in the core/fragments, manually change mcsout
## or use "inFrag" and "inCore" nested lists above. "Anchor atoms" (and connected Hs)
## are automatically included in each fragment unless specifically stated to be "inCore"

msld_crn.MsldCRN(mcsout,outdir,inFrag,inCore,ChkQChange=True,verbose=True,debug=False)


#####################################################################
## (5) Write Ligand Parameters & the Charmm ALF input scripts
#msld_prm.MsldPRM(outdir,cgenff,verbose=True,debug=True)

msld_prm.MsldPRM(outdir,cgenff,verbose=False,debug=False)
msld_wrt.writeALF_Files(sysname,outdir,cgenff)


## Final Notes to the user
print("default TOPPAR parameters copied into build."+sysname+". Check to make sure these work for your system!")


## FINISHED

