# Set of functions to analyze and view MSLD trajectories
# Written by: Luis F Cervantes and Furyal Ahmed (2/23)

# from pymol import cmd
import numpy as np
import sys
import string
import os
import subprocess
from functools import reduce
import shutil
from itertools import product

import pycharmm
import pycharmm.lingo as lingo
import pycharmm.generate as gen
import pycharmm.ic as ic
import pycharmm.coor as coor
import pycharmm.energy as energy
import pycharmm.dynamics as dyn
import pycharmm.nbonds as nbonds
import pycharmm.minimize as mini
import pycharmm.crystal as crystal
import pycharmm.image as image
import pycharmm.psf as psf
import pycharmm.read as read
import pycharmm.write as write
import pycharmm.settings as settings
import pycharmm.cons_harm as cons_harm
import pycharmm.cons_fix as cons_fix
import pycharmm.select as select
import pycharmm.shake as shake
import pycharmm.charmm_file as charmm_file

from pycharmm.lib import charmm as libcharmm

sysname = 'SYSNAME'  # System name (str)
nsubs = [[],[]]       # [nsubs1, nsubs2, ...] List of number of subs at each site
psfPath = './prep'
sysInpPath = './prep'

psfFile=os.path.join(psfPath,'minimized.psf')
sysInpFile= os.path.join(sysInpPath, f'{sysname}.inp')


def take_overlap(*input):
    """
    Given an unpacked list `input` of np.arrays,
    take the overlap of the arrays given the first
    column of the arrays.

    In:
    `input` : unpacked list of np.arrays

    Out:
    `result`: list of np.arrays
    """
    n = len(input)
    maxIndex = max(array[:, 0].max() for array in input)
    indicator = np.zeros(maxIndex + 1, dtype=int)
    for array in input:
        indicator[array[:, 0]] += 1
    indicator = indicator == n

    result = []
    for array in input:
        # Look up each integer in the indicator array
        mask1 = indicator[array[:, 0]]
        # Use boolean indexing to get the sub array
        result.append(array[mask1])

    return result


def get_subs_on(LambdaFile, nsubs, cutoff=0.99, nsavc=10000, nsavl=10):
    """
    Create subs_on matrix of shape (nframes x nsites), where 
    A_{i,j} is the substituent index that is on at the ith frame 
    and jth site (zero-based indexing)

    In:
        LambdaFile           (str) : path to lambda trajectory file
        nsubs       [int, int,...] : list of nsubs per site
        cutoff             (float) : lambda cutoff
        nsavc                (int) : freq of saving frames
        nsavl                (int) : freq of saving lambdas

    Out:
        subs_on      2-D np.array  : subs_on matrix
    """
    assert (nsavc/nsavl).is_integer, f"Frequency of saving lambdas and frames are not multiples of each other"
    nsites = len(nsubs)
    skip = int(nsavc/nsavl)
    lams = np.loadtxt(LambdaFile)
    physical_subs = []
    lambdasPerFrame = []
    for site in range(nsites):
        if site == 0:
            index1 = 0
        else:
            index1 = np.cumsum(nsubs[:site+1])[site-1]
    
        index2 = np.cumsum(nsubs[:site+1])[-1] -1  
     
        lams_site = lams[skip-1::skip,index1:index2+1]
        lambdasPerFrame.append(lams_site)
        mask = lams_site >= cutoff
        subs_on = np.argwhere(mask)
        physical_subs.append(subs_on)

    lambdasPerFrame = np.stack(*lambdasPerFrame, 0) 
    subs_on = physical_subs[0]
    if nsites != 1:
        physical_subs = take_overlap(*physical_subs) # Only get fully physical end states
        subs_on = physical_subs[0]
        for arr in physical_subs[1:]: 
            subs_on = np.append(subs_on,np.reshape(arr[:,1],(subs_on.shape[0],1)),1)
    return subs_on, lambdasPerFrame 

def get_selections(sysInpFile):
    with open(sysInpFile,'r') as f:
        lines = f.readlines()

    # This method takes care of commented out !define lines and indentations.
    # Do not think this will break. If BLOCK atom definitions start with
    # something different than site{}_sub{}, then beginds should be changed
    # accordingly. 
    beginds = [i for i,l in enumerate(lines) if l.startswith('define site')]
    endinds = [i for i,l in enumerate(lines) if l.startswith('   none ) end')]
    atSels = []
    groupNames = []
    for _,(i,j) in enumerate(zip(beginds,endinds)):
        site = int(lines[i].split()[1].split('site')[1].split('_')[0]) - 1
        # while len(atSels) != site + 1:
        #     atSels.append({})
        sub = int(lines[i].split()[1].split('sub')[1]) - 1
        atSel = [l.split()[3] for l in lines[i+2:j] if l.split()[0] == 'atom']
        groupNames.append(f'site{site+1}_sub{sub+1}')
        atSels.append(atSel)
    return groupNames, atSels     



def visualize_pymol(lambdasPerFrame, psfFile, trajPath, atSels, groupNames, eqS=5, nunits=1,rep=0,include_bonds=False, centerLig=True):
    """
    Only works for trajectories of the same replica. Does not display first `eqS`
    nanoseconds since there is no lambda data for them. 
    """
    # Load psf for connectivity info
    cmd.load(psfFile,'traj')
    
    # Load `nunits` ns trajectories starting at `eqS+1`
    for i in range(eqS+1,eqS+nunits+1):
        # PyMol does not like additions to dcd extension suffix. So
        # we modify and create a temporary file to load into pymol
        trajFile=os.path.join(trajPath,f'{sysname}_prod{i}.dcd_{rep}')
        loadFile = trajFile.split('.dcd')[0]+'.dcd'
        shutil.copyfile(trajFile,loadFile) 
    
        # Load ns trajectory 
        cmd.load_traj(loadFile,'traj')
    
        # Remove temp file
        os.remove(loadFile)

    # Get trajectory information
    nframes = cmd.count_states(selection='traj')
    lambdas = lambdasPerFrame[0:nframes, :]
    cmd.show('stick','resname LIG')
    if centerLig:
        cmd.intra_fit('traj and resname LIG',1)
 
    # Define objects
    for isub, atSel in enumerate(atSels):
        pymolSel = 'resname LIG and (name '
        pymolSel += ' or name '.join(atSel)
        pymolSel += ' )'
        if include_bonds:
            for iframe in range(nframes):
                lam = lambdas[iframe,isub]
                cmd.set('stick_transparency',1-lam, pymolSel, iframe+1)
        else:
            cmd.extract(groupNames[isub], pymolSel)
            # cmd.create(groupNames[isub], pymolSel)
            # cmd.remove(f'traj and {pymolSel}')
 
    set_lambdas = []
    if not include_bonds:
        for isub, groupName in enumerate(groupNames):
            for iframe in range(nframes):
                lam = lambdas[iframe,isub]
                if isub == 0:
                    set_lambdas.append(lam)
                cmd.set(f'stick_transparency', 1-lam, groupName, iframe+1)

    # print(lambdas[:,0])
    # print(set_lambdas)
    return None

def streamFileToCHARMM(FilePath):
    """
    Read a stream file in CHARMM via pyCHARMM lingo
    In:
    FilePath  (str): path to file to stream

    Out:
    CHARMM output in stdout
    """
    # Open file
    with open(FilePath,'r') as f:
        stream = f.readlines()

    # Exclude title from stream
    stream = [line for line in stream if not line.startswith('*')]
    stream = ''.join(stream)

    # Stream
    lingo.charmm_script(stream)    

def getAtomCharges(rtfPath):
    """
    Given a path to an rtf/parameter str file, return a dictionary
    {atname:charge}
    """
    with open(rtfPath,'r') as f:
        lines = f.readlines()

    atoms = [l.split()[1] for l in lines if l.startswith('ATOM')]
    charges = [float(l.split()[3]) for l in lines if l.startswith('ATOM')]
    return dict(zip(atoms, charges)) 

def getCRNCharges(prepPath, comb):
    """
    Given a prep directory path containing *rtfs and the ligand combination,
    get charge dictionary
    """
    with open(os.path.join(prepPath,'core.rtf'),'r') as f:
        lines = f.readlines()

    for site, sub in enumerate(comb):
        with open(os.path.join(prepPath,f'site{site+1}_sub{sub+1}_pres.rtf'),'r') as f:
            lines.extend(f.readlines())

        
    atoms = [l.split()[1] for l in lines if l.startswith('ATOM')]
    charges = [float(l.split()[3]) for l in lines if l.startswith('ATOM')]

    return dict(zip(atoms,charges))


def createStream(ligPath, prepPath, comb):
    ffCharges  = getAtomCharges(os.path.join(ligPath,'lig.str'))
    crnCharges = getCRNCharges(prepPath, comb)
            
    with open(os.path.join(ligPath,'ff_charges.txt'),'w') as f:
        f.write('''* Stream files to convert to original force field charges
* (LFC 3/10)
*

''')
        for at in ffCharges.keys():
            f.write(f'scalar charge set {ffCharges[at]} sele atom LIG 1 {at} end\n')

    with open(os.path.join(ligPath,'crn_charges.txt'),'w') as f:
        f.write('''* Stream files to convert to renormalized charges
* (LFC 3/10)
*

''')
        for at in ffCharges.keys():
            f.write(f'scalar charge set {crnCharges[at]} sele atom LIG 1 {at} end\n')

    return crnCharges, ffCharges 

def calcDipole(segid='LIG'):
    lingo.charmm_script(f'coor dipole mass sele segid {segid} end') 
    dipoleVec = [] 
    for direction in ['XDIP','YDIP','ZDIP']: 
        dipoleVec.append(lingo.get_energy_value(direction)) 
   
    dipoleMag = lingo.get_energy_value('RDIP')
 
    return np.array(dipoleVec), dipoleMag


def calcDipAngle(dip1, dip2):
    unitVecs = []
    for vector in (dip1,dip2):
        unitVecs.append(vector / np.linalg.norm(vector))
    return np.degrees(np.arccos(np.clip(np.dot(*unitVecs), -1.0, 1.0))) 

def calcChargeDipDiff(ligPath, write=True):
    crnDipV, crnDipMag = calcDipole()
    
    streamFileToCHARMM(os.path.join(ligPath,'ff_charges.txt'))
    ffDipV, ffDipMag = calcDipole()
    
    dipAngle = calcDipAngle(crnDipV, ffDipV)      
    dipDiff = crnDipMag - ffDipMag

    if write:
        with open(os.path.join(ligPath,'dipole.txt'),'w') as f:
            f.write('dipdiff,dipangle\n') 
            f.write(f'{dipDiff},{dipAngle}') 

    else:
        return dipDiff, dipAngle  


def calcQRMSD(crnCharges, ffCharges, writePath=None):
    diff2 = []
    for at in ffCharges.keys():
        diff2.append( (crnCharges[at] - ffCharges[at])**2 )

    rmsd = np.sqrt(np.mean(diff2))

    if writePath:
        with open(os.path.join(writePath,'rmsd.txt'),'w') as f:
            f.write('rmsd\n') 
            f.write(str(rmsd))
 
    return rmsd


def performChecks(variablesFile):
    """
    Performs all checks necessary to run splitSystem
    1) Check for variables1.inp file
    2) Check for loaded cgenff and openbabel modules
    """
    assert os.path.exists(variablesFile), f"File '{variablesFile}' was not found."

    outp = subprocess.getoutput("which cgenff")
    outp1 = subprocess.getoutput("which openbabel")
    
    assert not "no cgenff" in outp, "Please load the CGenFF module."
    assert not "no openbabel" in outp1, "Please load the OpenBabel module."

 
 

def splitSystem(sysInpFile, psfFile, nsubs, minimize=True, writeDir='./physical_ligands',setupChargePert=False,variablesFile = f'./variables1.inp'):
    """
    Assumes MSLD setup system with hybrid ligands has the same name for its psf and 
    pdb files. Need to have a variables1.inp file. Need to have cgenff and openbabel
    loaded and available if setupChargePert option set to True. Assumes same force
    field was used for system setup (current CGenFF).
    """
    # Get sysInpFile dir location
    prepPath = os.path.split(sysInpFile)[0]

    # Create writeDir
    if not os.path.exists(writeDir):
        os.mkdir(writeDir)

    # Perform checks
    performChecks(variablesFile)
    quit()

    # Stream system input file to CHARMM  
    streamFileToCHARMM(variablesFile)
    streamFileToCHARMM(sysInpFile)

    # groupNames, atSels = get_selections(sysInpFile)
    # Groups = dict(zip(groupNames, atSels))

    # Get ligand combinations
    nsites = len(nsubs) 
    subb = [range(n) for n in nsubs] 
    ligCombs = list(product(*subb))

    # Iterate through each physical ligand comb
    for comb in ligCombs:
        # Remove hybrid ligand atoms and retain `comb` lig
        settings.set_bomb_level(-1)
        dontInclude = []
        for site, subRange in enumerate(subb):
            for sub in subRange:
                if sub != comb[site]:
                    # Delete site{site+1}_sub{sub+1} atoms
                    # # lingo.charmm_script(f'delete atom select site{site+1}_sub{sub+1} end ')
                    dontInclude.append(f'site{site+1}_sub{sub+1}')
        dontInclude = ' .or. -\n'.join(dontInclude)
        selection= f'sele ({dontInclude}) end'
        lingo.charmm_script(f'delete atom {selection}')
        settings.set_bomb_level(0)

        # Minimize
        if minimize:
            mini.run_sd(**{'nstep': 50,
                           'tolenr': 1e-6,
                           'tolgrd': 1e-3})

        # Make physical ligand dir
        ligPath = os.path.join(writeDir,'_'.join(list(map(lambda x: str(x+1),comb))))
        if not os.path.exists(ligPath):
            os.mkdir(ligPath)
        
        # Write to dir
        write.coor_pdb(os.path.join(ligPath,'system.pdb'))
        write.coor_card(os.path.join(ligPath,'system.crd'))
        write.psf_card(os.path.join(ligPath,'system.psf'))

        if setupChargePert:
            # Remove everything but ligand of interest
            settings.set_bomb_level(-1)
            lingo.charmm_script(f'delete atom sele .not. segid LIG end ')
            settings.set_bomb_level(0)

            # Write ligand
            write.coor_pdb(os.path.join(ligPath,'lig.pdb'))
            write.psf_card(os.path.join(ligPath,'lig.psf'))

            # Get ff parameters
            os.system(f"babel -ipdb {os.path.join(ligPath,'lig.pdb')} -omol2 {os.path.join(ligPath,'lig.mol2')}")
            os.system(f"cgenff {os.path.join(ligPath,'lig.mol2')} > {os.path.join(ligPath,'lig.str')}")

            # Create stream files to interconvert between charge sets
            ffCharges, crnCharges = createStream(ligPath, prepPath, comb)
            
            # Calculate dipole and charge differences between charge sets
            calcChargeDipDiff(ligPath)
            calcQRMSD(crnCharges, ffCharges,writePath=ligPath)


        # Reinitialize psf
        settings.set_bomb_level(-1)
        psf.delete_atoms(pycharmm.SelectAtoms().all_atoms())
        read.psf_card(psfFile)
        read.coor_card(psfFile.replace('.psf','.crd'))
    


if __name__ = '__main__':
    # Set up system for CRN->FF Correction
    splitSystem(sysInpFile, psfFile, nsubs, setupChargePert=True, minimize=False)
