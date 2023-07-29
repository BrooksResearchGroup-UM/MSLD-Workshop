#! /usr/bin/env python

from rdkit import Chem
from rdkit import ML
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import TemplateAlign
#from rdkit.Chem import MCS
from rdkit.Chem import rdFMCS
from rdkit.ML.Cluster import ClusterVis
from rdkit.ML.Cluster import ClusterUtils
from rdkit.Chem import rdRGroupDecomposition
from rdkit.Chem import PandasTools
from copy import deepcopy
import functools
from functools import reduce
import pandas as pd
import numpy as np
import sys,os
import csv
import re

####
#### Run with rdenv environment activated!
#### (source activate rdenv)
####

#############################################
# (1) Load in mol2 files 
#############################################

def load_sdf(fn):
    """
    Load sdf specified file fn into rdkit and returns list of 
    RDKit mol objects for each molecule in sdf. OpenBabel
    kekulizes molecules in the sdf file when converting from a
    mol2 file. This facilitates input into RDKit. RDKit reads
    in atoms in same order as specified in provided file.
    """
    suppl = Chem.SDMolSupplier(fn,removeHs=False)
    return [x for x in suppl]

def load_mol2(fn):
    """
    Load mol2 specified file fn into rdkit and returns mol object.
    If molecule is kekulized in mol2, this is the better way of
    loading into RDKit.
    """
    mol = Chem.MolFromMol2File(fn,removeHs=False)
    if mol == None:
        print("Could not load %s into RDKit successfully." % fn)
        quit()
    return mol

def fix2DDepiction(mol): 
    """
    If loaded into RDKit from an SDF, 2D projection from 3D coords
    sometimes results in structure with overlapping atoms. This
    function converts the mol object into a smiles and then uses 
    that smiles conversion to depict molecule without overlapping
    atoms. Atom ordering is assumed not to change. Function returns
    rdkit mol object created from smiles and indices of mol that
    correspond to the smiles molecule indices. Note that the input
    molecule's atoms' coordinates are now changed.
    """
    smi = Chem.MolToSmiles(mol)
    x = Chem.MolFromSmiles(smi,sanitize=False)
    AllChem.Compute2DCoords(x)
    AllChem.GenerateDepictionMatching2DStructure(mol,x)
    return x, mol.GetSubstructMatch(x)

def code_cgenff(fn):
    """
    Helper function takes in .txt file name fn containing all
    CGenFF atom types and returns dictionaries to encode types
    into numbers and back 
    """
    # Define CGenFF atom type dictionaries 
    f = open(fn,'r')
    cgenfftypes = f.readlines()
    f.close()
    cgenfftypes = [x.rstrip() for x in cgenfftypes]
    typenum = [i+1 for i in range(len(cgenfftypes))]
    CGenFF2Num = dict(zip(cgenfftypes,typenum))
    Num2CGenFF = dict(zip(typenum, cgenfftypes))
    return CGenFF2Num, Num2CGenFF 

def get_current_atom_types():
    """
    Get current atom types in "CGenFF_atomtypes.txt"
    """
    attypesPath = os.path.dirname(__file__)
    with open(f"{attypesPath}/CGenFF_atomtypes.txt",'r') as f:
        current_ats = f.readlines()
    current_ats = [at.strip() for at in current_ats]
    return current_ats

def update_atomtypes(attypes):
    current_ats = get_current_atom_types()
    
    # Update CGenFF_atomtypes for any new atom types 
    current_ats.extend(attypes)
    current_ats = list(set(current_ats))
    attypesPath = os.path.dirname(__file__)
    if os.path.exists(f"{attypesPath}/CGenFF_atomtypes.txt"):
        os.remove(f"{attypesPath}/CGenFF_atomtypes.txt")
    with open(f"{attypesPath}/CGenFF_atomtypes.txt",'w') as f:
        f.write("\n".join(current_ats))


def read_rtf(fn):
    """
    Function takes in rtf file name + extension and returns 3 lists*:
    1) Atom Types (where atom types are coded according to a dictionary
    defined herein - uses CGenFF_atomtypes.txt)
    2) Partial Charges
    3) Actual CGenFF Atom Types
    *Lists are ordered according to order in rtf.
    """
    # Read Atom Types and Partial Charges
    f = open(fn,'r')
    lines = f.readlines()
    f.close()
    lines = [line.split() for line in lines]
    attypes = []
    partial_charges = []
    atnames = []
    for line in lines:
        if line:
            if line[0] == 'ATOM':
                if 'LP' not in line[1]:
                    atnames.append(line[1])
                    attypes.append(line[2])
                    partial_charges.append(line[3])

    update_atomtypes(attypes)
    attypesPath = os.path.dirname(__file__)
    CGenFF2Num, Num2CGenFF = code_cgenff(f"{attypesPath}/CGenFF_atomtypes.txt")
    Attypes = list(map(lambda x: CGenFF2Num[x], attypes))
    return Attypes, partial_charges, attypes, atnames


def get_property(fn, prop):
    """
    Helper function for read_biovia_mol2. Takes in mol_property prop
    and returns as a list of lists for each molecule (index of outer list)
    and atom (index of inner lists).
    """
    f = open(fn, 'r')
    line = f.readline()
    MATCHTypes = []
    while line:
        if len(line.split()) == 4:
            if line.split()[1] == 'Creation':
                line = f.readline() # Indexed at a molecule start
                while line:
                    if line.split():
                        if line.split()[0] == '@<SCITEGIC>MOL_PROPERTY':
                            line = f.readline()
                            if line.split()[0] == prop:
                                line = f.readline()
                                line = f.readline() # Indexed at atom types
                                molattypes = [] # Create new list for molecule atom types
                                while line.split():
                                    molattypes.append(line.split()[0])
                                    line = f.readline() 
                                MATCHTypes.append(molattypes)


                    line = f.readline()
 
        line = f.readline()
    f.close()
    return MATCHTypes


def read_biovia_mol2(fn):
    """
    Gets Biovia mol2 atom types assigned by MATCH and associated
    partial charges. Returns 3 lists*:
    1) Atom Types (where atom types are coded according to dictionary)
    2) Partial Charges
    3) Actual CGenFF atom types
    * Since Biovia contains multiple molecules in a mol2, this function
    assumes that all molecules are in a single mol2 and in same order
    as mol001...mol*.
    """
    MATCHTypes = []
    Actualtypes = get_property(fn,'ForcefieldType')
    MATCHpcs = get_property(fn,'PartialCharge')
    attypesPath = os.path.dirname(__file__)
    CGenFF2Num, Num2CGenFF = code_cgenff(f"{attypesPath}/CGenFF_atomtypes.txt") 
    for i in Actualtypes:
        MATCHTypes.append(list(map(lambda x: CGenFF2Num[x],i)))
    return MATCHTypes, MATCHpcs, Actualtypes 


def get_cg_rtf_types(fn):
    """
    Helper function for read_cg_rtf. Finds property (prop) for each
    atom in a molecule rtf generated by charmm.
    """
    f = open(fn,'r')
    line = f.readline()
    attypes = []
    partial_charges = []
    lpidx = []
    while line:
        if line.split():
            if line.split()[-1] == 'LIG':
                while line:
                    if line.split():
                        if line.split()[0:2] == ['ATOM','TYPE']:
                            rawlist=line.split()[2:]
                            lpindex = list(filter(lambda x: rawlist[x][0:2] == 'LP', range(len(rawlist))))
                            lpidx.append(lpindex)
                            nolplist = list(filter(lambda x: x[0:2] != 'LP', rawlist))
                            attypes.extend(nolplist)               
                            line = f.readline()
                            if line.split():
                                while line.split()[0] != 'TYPE':
                                    lpind = list(filter(lambda x: line.split()[x][0:2] == 'LP', range(len(line.split()))))
                                    lpidx.append(lpind)
                                    nolp = list(filter(lambda x: x[0:2] != 'LP', line.split()))
                                    attypes.extend(nolp)               
                                    
                                    line = f.readline()

 
                    line = f.readline()             



        line = f.readline()
                
    lpidx = [list(map(lambda x: x - 1, lpidx[i])) for i in range(len(lpidx))]
    f.close()
    return attypes, lpidx


def get_cg_rtf_pcs(fn,lpidx):
    """
    Helper function for read_cg_rtf. Finds property (prop) for each
    atom in a molecule rtf generated by charmm. lpdix is the lone
    pair index for each patch iteration gotten from
    get_cg_rtf_types(). 
    """
    f = open(fn,'r')
    lpidx = [list(map(lambda x: x +1, i)) for i in lpidx]
    line = f.readline()
    attypes = []
    partial_charges = []
    count = 0
    while line:
        if line.split():
            if line.split()[-1] == 'LIG':
                while line:
                    if line.split():
                        if line.split()[0] == 'CHARGE':
                            rawlist=line.split()[1:]
                            nolplist = rawlist
                            if lpidx[count]:
                                nolplist = [j for i,j in enumerate(rawlist) if i not in lpidx[count]]
                            attypes.extend(nolplist)
                            
                            count+=1               
                            line = f.readline()
                            if line.split():
                                while line.split()[0] != 'ALPHA':
                                    rawlist = line.split()
                                    nolplist = rawlist
                                    if lpidx[count]:
                                        nolplist = [j for i,j in enumerate(rawlist) if i not in lpidx[count]]
                                    attypes.extend(nolplist) 
                                    count+=1
                                    line = f.readline()

 
                    line = f.readline()             



        line = f.readline()
                

    f.close()
    return attypes

def read_cg_rtf(fn):
    """
    Function reads charm generated (cg) rtf (fn) and returns:  
    1) Atom Types (where atom types are coded according to a dictionary
    defined herein - uses CGenFF_atomtypes.txt)
    2) Partial Charges
    3) Actual CGenFF Atom Types
    *Lists are ordered according to order in rtf.
    Since pieces of the molecules are patched, this function looks for
    core first (LIG) and then for the other patches after this.
    """
    attypesPath = os.path.dirname(__file__)
    CGenFF2Num, Num2CGenFF = code_cgenff(f"{attypesPath}/CGenFF_atomtypes.txt") 
    # Read Atom Types and Partial Charges
    Actualtypes,lpindex = get_cg_rtf_types(fn) 
    pcs = get_cg_rtf_pcs(fn,lpindex)
    qnormtypes = list(map(lambda x: CGenFF2Num[x],Actualtypes))

    return qnormtypes, pcs, Actualtypes 
    

def get_ligand_names(combos):
    """
    Helper function takes in combos ("combinations.txt") file with
    all ligand combinations specified and returns a list of rtf file names.
    """
    f = open(combos,'r')
    lines = f.readlines()
    f.close()
    lines = [x.rstrip().split() for x in lines]
    ss=""
    for i in range(len(lines[0])):
        ss = ss+"s"+str(i+1)+"s%s."
    ss += 'rtf'
    fnames = [ss % tuple(lines[i]) for i in range(len(lines))]
    return fnames

def get_og_lig_names(ogfile):
    """
    Helper function takes in combos ("combinations.txt") file with
    all ligand combinations specified and returns a list of rtf file names.
    """
    f = open(ogfile,'r')
    lines = f.readlines()
    f.close()
    lines = [x.rstrip() for x in lines]
    fnames = list(map(lambda x: x+".sdf",lines))
    return fnames

def read_cg_rtfs(combos):
    """
    Function reads combos file ("combinations.txt") and uses read_cg_rtf() to read
    in all molecule combinations in that file. Returns three lists of lists:
    1) Encoded atom types per molecule
    2) Atom partial charges per molecule
    3) Actual atom types per molecule
    """
    fnames = get_ligand_names(combos)
    CodedTypes = []
    PCs = []
    Types = []
    for rtf in fnames:
        codedtypes, pcs, actualtypes = read_cg_rtf(rtf)
        if len(codedtypes) != len(pcs) and len(pcs) != len(actualtypes):
            print("CHECK %s!! Mismatch in number of atoms for Atom Types and Charges." % rtf)
        CodedTypes.append(codedtypes)
        PCs.append(pcs)
        Types.append(actualtypes)

    return CodedTypes, PCs, Types


def read_og_rtfs(ognames):
    """
    Function reads combos file ("combinations.txt") and uses read_cg_rtf() to read
    in all molecule combinations in that file. Returns three lists of lists:
    1) Encoded atom types per molecule
    2) Atom partial charges per molecule
    3) Actual atom types per molecule
    """
    fnames = list(map(lambda x: x.replace(".sdf",".rtf"),ognames))
    CodedTypes = []
    PCs = []
    Types = []
    Names = []
    for rtf in fnames:
        codedtypes, pcs, actualtypes, atnames = read_rtf(rtf)
        if len(codedtypes) != len(pcs) and len(pcs) != len(actualtypes) and len(actualtypes) != len(atnames):
            print("CHECK %s!! Mismatch in number of atoms for Atom Types and Charges." % rtf)
        CodedTypes.append(codedtypes)
        PCs.append(pcs)
        Types.append(actualtypes)
        Names.append(atnames)

    return CodedTypes, PCs, Types, Names

def split_biovia_mol2(fn):
    """
    Helper fuunction splits Biovia mol2 into individual mol2 files.
    """
    f = open(fn,'r')
    line = f.readline()
    s = 0
    while line:
        if line.split():
            if line.split()[0] == "@<TRIPOS>MOLECULE":
                ff = open("biovia-%s.mol2" % (s+1), 'w')
                while line[0:13] != "@<TRIPOS>SUBS":
                    ff.write(line)
                    line = f.readline()
                ff.close()
                s+=1

        line = f.readline()
    f.close()



def match_atoms_lists(qnormmol,bioviamol, proplist1, proplist2):
    """
    Function rearranges MATCH molecules to have same atom indices as those for
    the CGenFF molecules based on their respective mol2 files and sdf files.
    """
    rearrange_idxs = bioviamol.GetSubstructMatch(qnormmol) # Get indices of bioviamol that correspond to qnormmol
    if len(rearrange_idxs) != len(proplist1):
        print("COULD NOT FIND SUBSTRUCTURE BETWEEN QNORM AND BIOVIA STRUCTURES")
        quit()
    print("\n")
    print(len(rearrange_idxs))
    print(len(proplist2))
    newlist = [proplist2[i] for i in list(rearrange_idxs)]

    return proplist1, newlist

def mol_with_atom_index( mol ):
    """
    Function returns mol object with labeled atoms based on index. Can be converted to 
    a png file for visualization.
    """
    atoms = mol.GetNumAtoms()
    for idx in range( atoms ):
        mol.GetAtomWithIdx( idx ).SetProp( 'molAtomMapNumber', str( mol.GetAtomWithIdx( idx ).GetIdx() ) )
    return mol

def get_differences(qnormmol, bioviamol, proplist1, proplist2):
    """
    Function takes in rdkit mols file and gets differences between:
    1) Atom Types
    2) Partial Chargesi
    proplist1 and proplist2 are ordered according to qnormfn and
    bioviafn, respectively. They need not be aligned.
    """
    qnormlist, biovialist = match_atoms_lists(qnormmol,bioviamol, proplist1, proplist2)
    diff_list = list(map(lambda x,y: float(x) - float(y), qnormlist, biovialist))
    
    return diff_list

def write_png(name, qnormmol, bioviamol, qnormprop, bioviaprop, differences):
    """
    Function takes in qnormmol and respective qnormprop (same order of atoms),
    bioviamol and bioviaprop (same order of atoms), and differences (same order
    of atoms as qnormmol) and writes a .png file containing three molecules (same)
    with qnorm properties, biovia properties, and differences properties. Variable
    name is to specify the name of the png file.
    """
    atoms1 = qnormmol.GetNumAtoms()
    atoms2 = bioviamol.GetNumAtoms()
    if atoms1 != atoms2:
        print("ATOM NUMBER MISMATCH IN %s" % name)
        quit()
    diffmol = Chem.Mol(qnormmol)

    for idx in range(atoms1):
        qnormmol.GetAtomWithIdx( idx ).SetProp( 'atomLabel', str( round(float(qnormprop[idx]), 3 ) ) )
        bioviamol.GetAtomWithIdx( idx ).SetProp( 'atomLabel', str( round(float(bioviaprop[idx]), 3) ) )
        diffmol.GetAtomWithIdx( idx ).SetProp( 'atomLabel', str( round(float(differences[idx]),3 )))
    mols = [qnormmol, bioviamol, diffmol]
    img=Draw.MolsToGridImage(mols,molsPerRow=3,subImgSize=(400,400),legends=("Charge Renormalized CGenFF", "Original CGenFF", "Difference"))

    img.save("%s.png" % name)


def check_atoms(mname, qnormmol, bioviamol, proplist1, proplist2):
    qnormlist, biovialist = match_atoms_lists(qnormmol, bioviamol, proplist1, proplist2)
    print(qnormlist)
    print(biovialist)
    if qnormlist == biovialist:
        print("Same atom types for molecule %s" % mname)


def get_all_differences(combos, bioviaf,prop):
    """
    Function takes in combos file and biovia file name and returns two lists of lists:
    1) Atom Types differences
    2) Partial Charges differences
    prop argument is 'attypes' or 'charges'
    """
    split_biovia_mol2(bioviaf)
    l,m,n = read_biovia_mol2(bioviaf)
    o,p,q = read_cg_rtfs(combos)
    if prop == 'attypes':
        prop = [o,l]
    if prop == 'charges':
        prop = [p,m]
    attypes = [q, n]
    qnormsdfs = list(map(lambda x: x.replace(".rtf", ".sdf"), get_ligand_names(combos))) 
    differences = []
    for i in range(len(o)):
        print("PROCESSING MOLECULE %s WITH INDEX %s" % (qnormsdfs[i].replace(".sdf",""), i))
        x = load_mol2("biovia-%s.mol2" % (i+1))
        y = load_sdf(qnormsdfs[i])
        y = y[0]
        smi1,_ = fix2DDepiction(y)
        smi2,_ = fix2DDepiction(x)
        if i == 0:
            print("HELLO")
            Draw.MolToFile(x, "Biovia-%s.png" % i)
            Draw.MolToFile(y,"Qnorm-%s.png" % i)
            lst = x.GetSubstructMatch(y)
            print(lst)
        check_atoms(i, y, x, attypes[0][i],attypes[1][i])
        diff_list = get_differences(y, x, prop[0][i], prop[1][i])
        differences.append(diff_list)
        write_png("molecule-%s" % i, y, x, prop[0][i], prop[1][i],diff_list)
    return differences 


def get_all_diffs_qrn_og(combos, ognames, prop):
    """
    Function takes in combos and ognames text files that have the substituents
    at a site per molecule and the corresponding mol0* molecule, respectively. Output
    is same as get_all_differences.
    """
    qnromrtfs = get_ligand_names(combos)
    qnormsdfs = list(map(lambda x: x.replace(".rtf", ".sdf"), qnromrtfs))
    ogrtfs = get_og_lig_names(ognames) 
    ogsdfs = list(map(lambda x: x.replace(".rtf",".sdf"), ogrtfs))

    l,m,n = read_og_rtfs(ogrtfs)
    o,p,q = read_cg_rtfs(combos)

    if prop == 'attypes':
        prop = [o,l]
    if prop == 'charges':
        prop = [p,m]
    attypes = [q, n]
    differences = []
    for i in range(len(ogrtfs)):
        print("PROCESSING MOLECULE %s WITH INDEX %s" % (qnormsdfs[i].replace(".sdf",""), i))
        x = load_sdf(ogsdfs[i])[0]
        y = load_sdf(qnormsdfs[i])[0]
        smi1,_ = fix2DDepiction(y)
        smi2,_ = fix2DDepiction(x)
        check_atoms(qnormsdfs[i].replace(".sdf",""), y, x, attypes[0][i],attypes[1][i])
        diff_list = get_differences(y, x, prop[0][i], prop[1][i])
        differences.append(diff_list)
    return differences 


def check_rdkit_loading(mols):
    """
    Checks if RDKit successfully loaded each molecule in a list. 
    Basically just checks if there is a None element in the provided
    list. If all molecules were not successfully loaded, it provides
    the index of the ligand in mol_list.txt that were not loaded.
    """
    notsuccessful = list(filter(lambda mol: mols[mol] == None, range(len(mols))))
    if notsuccessful:
        print("Ligands with indices %s were not loaded successfully into RDKit" % notsuccessful)
        quit()
    else:
        print("All ligands were successfully loaded into RDKit") 


def check_MCS_results(MCSIndices):
    """
    Checks if number of atoms in MCS result is the same for every
    molecule. MCSIndices is a list of tuples where the outer index
    corresponds to the molecule index and the tuple index corresponds
    to the atom within a molecule.
    """
    check = list(map(lambda x: len(MCSIndices[x])==len(MCSIndices[x+1]), range(len(MCSIndices)-1)))
    if False in check:
        print("MCS Substructure pattern could not be properly found across all molecules")
        quit()
    else:
        print("Number of atoms in MCS result is the same for every molecule\n")


def groups_to_df_html(htmlname, groups,mols,include_core=False,redraw_sidechains=False):
    """
    Create an html file containing results from the R Group decomposition. 
    """
    cols = ['Mol']+list(groups.keys())
    if redraw_sidechains:
        for k,vl in groups.items():
            if k=='Core':
                continue
            for i,v in enumerate(vl):
                vl[i] = Chem.RemoveHs(v)
                rdDepictor.Compute2DCoords(vl[i])

    
    if not include_core:
        cols.remove('Core')
        del groups['Core']
    groups['Mol'] = mols
    frame = pd.DataFrame(groups,columns=cols)
    PandasTools.ChangeMoleculeRendering(frame, renderer="PNG")
    dfhtml = open(htmlname, "w")
    dfhtml.write(frame.to_html())
    dfhtml.close()
    return frame


def fragment(mol, mcspatt, bond_ind):
    i, at_cut = bond_ind
    bond = mol.GetBondBetweenAtoms(i, at_cut)
    frag = Chem.FragmentOnBonds(mol, [bond.GetIdx()])
    fragsmiles = Chem.MolToSmiles(frag, isomericSmiles=True).split('.')
    # print(fragsmiles)
    
    #fragsmiles = [x.split('*]')[1] for x in fragsmiles]
    fragmols = [Chem.MolFromSmiles(smiles) for smiles in fragsmiles]
    
    for m in fragmols:     
        MaxComSubst = rdFMCS.FindMCS([m, mcspatt],
                  atomCompare=rdFMCS.AtomCompare.CompareIsotopes,
                  ringMatchesRingOnly=False,
                  completeRingsOnly=False,
                  timeout=5,
                  verbose=True
                  )
        if not MaxComSubst.smartsString:
            return m
    return None     

def r_group_decomposition(mol, core_indices, mcspatt):
    """
1) 
    Given a scaffold, iterate through each atom in the molecule
that belong in the scaffold. Get neighbors and if all neighboring
atoms are in the scaffold, then discard. Otherwise add to
anchor atom list [[],[]] list of lists (in case of multiple
atoms at one R group)


2) Identify number of R sites:
    Iterate through each anchor atom and look at their neighbors,
    if the atom has a neighbor in the anchor atom, add anchor
    atom to touple and eliminate neighbor anchor atom - these belong
    to the same site.

3) Shortest paths to anchor atoms:
    Make a list that looks like the anchor atoms list of lists
    For each atom not in the core, identify the minimum shortest path
    to each of the anchor atoms identified (if more than one per site, use
    the first one)
    Add the atom to the site with the shortest path
    """
    mcsatoms = re.findall('\d+', Chem.MolToSmiles(mcspatt))
    anchor_inds = []
    counter = 0
    # for i, at in enumerate(mol.GetAtoms()):
    #     if i in core_indices:
    for j in core_indices: 
        for i, at in enumerate(mol.GetAtoms()):
            if i == j:
                neighbor_idxs = list(map(lambda x: x.GetIdx(), at.GetNeighbors()))
                set_diff = set(neighbor_idxs) - set(core_indices)
                if set_diff:
                    anchor_inds.append(i)

    anchors = []
    all_frags = []
    # for i, at in enumerate(mol.GetAtoms()):
    #     if i in anchor_inds:
    for j in anchor_inds: 
        for i, at in enumerate(mol.GetAtoms()):
            if i == j:
                neighbor_idxs = list(map(lambda x: x.GetIdx(), at.GetNeighbors()))
                intersection = set(neighbor_idxs).intersection(set(anchor_inds))
                to_cut = list(set(neighbor_idxs) - intersection - set(neighbor_idxs).intersection(set(core_indices)))
                # frags = []
                # for at_cut in to_cut:
                #     m = fragment(mol, mcspatt, (i,at_cut))
                #     if m is not None:
                #         fragsmiles = Chem.MolToSmiles(m)
                #         frags.append(fragsmiles)

                # all_frags.append(frags) 
                anchors.append([i])
                if intersection:
                    anchors[-1].extend(list(intersection))
                    for at in intersection:
                        anchor_inds.remove(at)
                neighbor_idxs = set(neighbor_idxs) - intersection
                print(neighbor_idxs)
                # quit()
    R_groups = [[] for i in anchor_inds]
    shortest_paths = [[] for i in anchor_inds]
    for i, at in enumerate(mol.GetAtoms()):
        if i not in core_indices:
            path_lengths = [len(Chem.rdmolops.GetShortestPath(mol,i,j)) for j in anchor_inds]
            R_groups[np.argmin(path_lengths)].append(i)
            shortest_paths[np.argmin(path_lengths)].append(min(path_lengths))
    
    # Rearrange based on shortest paths
    for i, group in enumerate(R_groups):
        R_groups[i] = [x for _, x in sorted(zip(shortest_paths[i], group))]
        shortest_paths[i] = [y for y, _ in sorted(zip(shortest_paths[i], group))]


    return R_groups, shortest_paths, anchors 


def writeMCS(molnames, corenames, rgroupnames, anchor_names,mcsout="MCS_for_MSLD.txt"):
    nsites = len(rgroupnames)
    reflig = molnames[0]
    text = f"""# Maximum Common Substructure Search for MSLD (JV/LC 2023)
# {len(molnames)} molecules processed

NSUBS {' '.join([str(len(x)) for x in rgroupnames])}

REFLIG {reflig}

CORE
"""
    
    for mol, core_atoms in zip(molnames, corenames):
        text += f"{mol} {' '.join(core_atoms)}\n"

    text += "\nANCHOR ATOMS\n"
    for i, (mol,anchs) in enumerate(zip(molnames,anchor_names)):
        anchors = [x[0] if len(x) == 1 else "DUM" for x in anchs]
        text += f"{mol} {' '.join(anchors)}\n"

    text +="\n"

    for site in range(nsites):
        text += f"SITE {site+1} FRAGMENTS\n" 
        for mol,gnames in zip(molnames, rgroupnames[site]):
            text += f"{mol} {' '.join(gnames)}\n"
        text += '\n'

    text += 'END'

    with open(mcsout,'w') as fout:
        fout.write(text)

    return reflig

def MCSS_RDecomp(mol_list,mcsout="MCS_for_MSLD.txt"):
    """
    Function takes in a list of ligands, identifies the common core and
    shows charge distribution of each common core atom across all ligands
    in core_charges.csv. It also shows matching core atom names for each 
    molecule in core_names.csv
    """
    sdfnames = get_og_lig_names(mol_list)
    rtfnames = list(map(lambda x: x.replace(".sdf",".str"),sdfnames))
    molnames = list(map(lambda x: x.replace(".str",""), rtfnames))
     
    # Get partial charge, type, and name of every atom in each molecule 
    CodedTypes, PartialCharges, AtomTypes, AtomNames = read_og_rtfs(rtfnames)
    
    # Load ligands into RDKit
    mols = [load_sdf(mol)[0] for mol in sdfnames]
    
    # Check if load was successful
    check_rdkit_loading(mols) 
   
    # Assign atom types to each atom of each RDKit molecule object
    # RDKit uses "Isotope Labeling" for this
    for molidx, mol in enumerate(mols): 
        for atidx, atom in enumerate(mol.GetAtoms()):
            atom.SetIsotope(CodedTypes[molidx][atidx])

    # Fix 2D Depiction of molecules
    for mol in mols:
        _, _ = fix2DDepiction(mol)
 
    # Do MCS Search
    MaxComSubst = rdFMCS.FindMCS(mols,
                              atomCompare=rdFMCS.AtomCompare.CompareIsotopes,
                              ringMatchesRingOnly=False,
                              completeRingsOnly=False,
                              timeout=5,
                              verbose=True
                              )

    # Get MCS search results in a smarts string and load as mol object
    # into RDKit
    reflig = MaxComSubst.queryMol
    patt = MaxComSubst.smartsString
    print(patt)
    print("\nMCS search Smarts String Result is:\n%s\n" % patt)
    patt = Chem.MolFromSmarts(patt)
    # Get indices for each atom that match MCS pattern
    MatchIndices = [mol.GetSubstructMatch(patt) for mol in mols]

    # Check to see if number of indices is the same for each molecule
    check_MCS_results(MatchIndices)

    # Rearrange partial charges and atom names based on core atom indices
    CoreCharges = []
    CoreNames = []
    for molidx in range(len(MatchIndices)):
        core_partial_charges = [PartialCharges[molidx][atom] for atom in MatchIndices[molidx]]
        core_atom_names = [AtomNames[molidx][atom] for atom in MatchIndices[molidx]]
        CoreCharges.append(core_partial_charges)
        CoreNames.append(core_atom_names)

    # Write common core charges into dataframe
    charges_data = np.array(CoreCharges)
    df = pd.DataFrame(data=charges_data, index=molnames, columns=CoreNames[0])
    df.to_csv("core_charges.csv") 

    # Write common core atom names into dataframe
    names_data = np.array(CoreNames,dtype="<U50")
    df = pd.DataFrame(data=names_data, index=molnames, columns=CoreNames[0])
    df.to_csv("core_names.csv") 

    # Clear Isotopes
    for mol in mols:
        for at in mol.GetAtoms():
            at.SetAtomMapNum(at.GetIdx())

    # Perform R group decomposition based on MaxSubstruct search pattern patt
    print("Beginning R Group Decomposition")
    params = rdRGroupDecomposition.RGroupDecompositionParameters()
    # # params.alignment = rdRGroupDecomposition.RGroupCoreAlignment.MCS
    params.onlyMatchAtRGroups = True
    params.labels = rdRGroupDecomposition.RGroupLabels.AtomIndexLabels
    params.rgrouplabelling = rdRGroupDecomposition.RGroupLabelling.Isotope
    groups, unmatched = rdRGroupDecomposition.RGroupDecompose([patt],mols,asSmiles=False,asRows=False,options=params)
   
    print("R Group Decomposition complete") 
    # Check if any molecules were unmatched
    if not groups:
        print("R Group Decomposition Returned the following results:\n%s" % groups)
    if unmatched:
        print("One or more molecules were unmatched during R group decomposition step.\nDefaulting to custom R decomposition. No html or xlsx file will be created. This will most likely fail for highly symmetric cores\nbecause it assumes that the core atom indices are arranged such that the indices per molecule correspond to the same atom across all molecules")
        r_group_decomposition(mols[0],MatchIndices[0], patt)
        #quit()
        decomp_results = [r_group_decomposition(mol, indices, patt) for mol,indices in zip(mols, MatchIndices)]
        # fragment(mols, decomp_results[0], decomp_results[2])
        #quit()
        R_groups = []
        paths = []
        anchors = []
        nsites = len(decomp_results[0][0])
        for i in range(nsites):
            R_groups.append([decomp_results[j][0][i] for j,_ in enumerate(mols)])
            paths.append([decomp_results[j][1][i] for j,_ in enumerate(mols)])
        for mol,_ in enumerate(decomp_results):
            anchs = decomp_results[mol][2]
            to_append = []
            for siteanch in anchs:
                to_append.append([AtomNames[mol][at] for at in siteanch])
            anchors.append(to_append)
       
         
        R_group_names = [[] for i in range(nsites)]
        for site in range(nsites):
            for molidx in range(len(mols)):
                R_group_names[site].append([AtomNames[molidx][atom] for atom in R_groups[site][molidx]])
  
         
        reflig = writeMCS(molnames, CoreNames, R_group_names,anchors, mcsout=mcsout)

        return reflig

    # Fix 2D Depiction of every "core" fragment for every molecule
    for i in range(len(mols)):
        _,_ = fix2DDepiction(groups['Core'][i])
    
    # Write R group decomposition results to html and .xlsx
    img = Draw.MolToImage(groups['Core'][0], size=(800,800))
    img.save("core.png")
    frame = groups_to_df_html("r_group_decomposition.html", groups, mols, include_core=True) # Need to rewrite function from RDKit locally 
                                                                                             # to modify image size in dataframe
    PandasTools.SaveXlsxFromFrame(frame,"r_group_decomposition.xlsx", molCol='Core', size=(600,600)) # Need to fix by modifying function
                                                                                             # from RDKit locally to render images

    # Check if a decomposed R site is the same as another. Lets user know if two sites have the same substituents.
    nsites = sum(list(map(lambda x: x[0]=='R', groups.keys()))) # count number of sites
    todelete = []
    for site1 in range(nsites):
        for site2 in range(nsites):
            non_unique = []
            if site2 > site1:
                for frag1, frag2 in zip(groups["R%s" % str(site1+1)], groups["R%s" % str(site2+1)]):
                    print(Chem.CanonSmiles(Chem.MolToSmiles(frag1)))
                    frag1smiles = Chem.CanonSmiles(Chem.MolToSmiles(frag1,kekuleSmiles=True,allBondsExplicit=True,allHsExplicit=True))
                    frag2smiles = Chem.CanonSmiles(Chem.MolToSmiles(frag2,kekuleSmiles=True,allBondsExplicit=True,allHsExplicit=True))
                    if frag1smiles == frag2smiles:
                        non_unique.append(0) 
            if len(non_unique) == len(groups["R%s" % str(site1+1)]) and non_unique:
                print("Site %s has the same substituents as Site %s" % (str(site2+1), str(site1+1)))
                todelete.append("R%s" % str(site2+1))

    if todelete:
        nsites -= len(todelete)
        for key in todelete:
            del groups[key]
   
        sitekeys = list(filter(lambda x: x[0] == 'R',groups.keys())) 
        for idx,key in enumerate(sitekeys):
            groups["R%s" % str(idx+1)] = groups.pop(key)

    print("NEW NUMBER OF SITES IS %s" % nsites)
    print("NEW GROUPS IS: %s" % groups.keys())

    # Write R group decomposition results to html and .xlsx
    img = Draw.MolToImage(groups['Core'][0], size=(800,800))
    img.save("core.png")
    frame = groups_to_df_html("r_group_decomposition.html", groups, mols, include_core=True) # Need to rewrite function from RDKit locally 
                                                                                             # to modify image size in dataframe
    PandasTools.SaveXlsxFromFrame(frame,"r_group_decomposition.xlsx", molCol='Core', size=(600,600)) # Need to fix by modifying function
                                                                                             # from RDKit locally to render images
    # Check for unique fragments per site. Get dictionary Uniqueness
    Uniqueness = {}
    for site in range(nsites):
        non_unique = []
        non_unique_tuples = []
        Uniqueness["R%s" % str(site+1)] = []
        for frag1idx, frag1 in enumerate(groups["R%s" % str(site+1)]):
            if frag1idx not in non_unique:
                Uniqueness["R%s" % str(site+1)].append(frag1idx)
            else:
                Uniqueness["R%s" % str(site+1)].append([idx for idx in non_unique_tuples if idx[0] == frag1idx][0][1])
                 
            for frag2idx, frag2 in enumerate(groups["R%s" % str(site+1)]):
                if frag2idx > frag1idx and frag2idx <= len(groups["R%s" % str(site+1)])-1 and frag2idx not in non_unique:
                    frag1smiles = Chem.CanonSmiles(Chem.MolToSmiles(frag1))
                    frag2smiles = Chem.CanonSmiles(Chem.MolToSmiles(frag2))
                    # print(frag1smiles)
                    # print(frag2smiles)
                    if frag1smiles == frag2smiles:
                        print("For site %s, Molecule %s has the same substituent as Molecule %s" % (str(site+1),str(frag1idx),str(frag2idx)))
                        non_unique_tuples.append((frag2idx,frag1idx))
                        non_unique.append(frag2idx)
                        continue
            
    # Iterate through every molecule and every atom to match every R group to original mol files.
    RGroupNames = {} # Dictionary where every key is Ri, and value is list of lists for each molecule's frag names at site i. 
    RGroupCharges = {} # Dictionary where every key is Ri, and value is list of lists for each molecule's frag charges at site i.
    anchoratoms = [] 
    for site in range(nsites):
        RGroupNames["R%s" % str(site+1)] = []
        RGroupCharges["R%s" % str(site+1)] = []
        anchoratoms.append([])
        for fragidx, frag in enumerate(groups["R%s" % str(site+1)]):
            fragatlist = []
            pclist = []
            countanchors = 0
            for atidx, atom in enumerate(frag.GetAtoms()):
                if atom.GetIsotope() != 0: # Get rid of attachment dummy atom (not coded as isotope)
                    fragatlist.append(str(AtomNames[fragidx][atom.GetAtomMapNum()]))
                    pclist.append(round(float(PartialCharges[fragidx][atom.GetAtomMapNum()]),3))
                else:
                    countanchors += 1
            if countanchors > 1:
                anchoratoms[site].append("DUM")
            RGroupNames["R%s" % str(site+1)].append(fragatlist)
            RGroupCharges["R%s" % str(site+1)].append(pclist)

    print("anchoratoms %s" % anchoratoms)
    # Get Unique Fragments only for creation of MCS_for_MSLD.txt file 
    UniqueRGroupNames = deepcopy(RGroupNames)
    UniqueRGroupCharges = deepcopy(RGroupCharges)
    for site in range(nsites):
        keep = list(set(Uniqueness["R%s" % str(site+1)]))
        print("keep: %s" % keep)
        UniqueRGroupNames["R%s" % str(site+1)] = [UniqueRGroupNames["R%s" % str(site+1)][idx] for idx in keep]
        UniqueRGroupCharges["R%s" % str(site+1)] = [UniqueRGroupCharges["R%s" % str(site+1)][idx] for idx in keep]

    # To do: Get anchor atoms per unique substituent, get molecule names of unique substituent.
    # Get rid of duplicate sites (due to two sites of attachment)
    # finish getting anchor atom info
    for site in range(nsites):
        if "DUM" not in anchoratoms[site]:
            print("SITE %s" % str(site))
            attachment = []
            for frag1idx, frag1 in enumerate(groups["R%s" % str(site+1)]):
                print("FRAG %s" % frag1idx)
                fragatoms = [atom.GetAtomMapNum() for atom in frag1.GetAtoms() if atom.GetIsotope() != 0]
                for atidx, atom in enumerate(frag1.GetAtoms()):
                    if atom.GetIsotope() != 0:
                        print("ATIDX %s" % str(atidx))
                        fragneighboridxs = list(map(lambda x: x.GetAtomMapNum(), atom.GetNeighbors()))
                        molneighboridxs = list(map(lambda x: x.GetIdx(),mols[frag1idx].GetAtomWithIdx(atom.GetAtomMapNum()).GetNeighbors()))
                        #print(fragneighboridxs)
                        print(fragatoms)
                        print(molneighboridxs)
                        if not set(molneighboridxs).issubset(set(fragatoms)):
                            attachment.extend(set(molneighboridxs).difference(set(fragatoms)))
        else:
            continue
        anchoratoms[site] = [AtomNames[molidx][atomidx] for molidx, atomidx  in enumerate(attachment)]
    print(anchoratoms)
    
    RMolNames = deepcopy(Uniqueness)
    for site in Uniqueness.keys():
        uniqueindices = list(set(Uniqueness[site]))
        RMolNames[site] = [molnames[i] for i in uniqueindices]
    print(RMolNames)
    # Now that we have the Core and the R fragment names and anchor atoms, we can assemble the MCS_for_MSLD.txt file
    names_data = np.insert(names_data, 0, molnames,axis=1)
    nsubs = [len(set(Uniqueness[site])) for site in Uniqueness.keys()]
    substring = ' '.join([' '.join(str(i)) for i in nsubs])
    print(substring)
    fn = open("MCS_for_MSLD.txt", 'w', newline='')
    fn.write("""# MCS Search for MSLD (LFC/JZV 2020)
# %s molecules processed.

NSUBS %s

REFLIG %s

CORE
""" % (len(molnames),substring, molnames[0]))
    np.savetxt(fn,names_data,fmt='%s')
    fn.write("\nANCHOR ATOMS\n")
    anchor_data = np.array(anchoratoms,dtype="<U50")
    anchor_data = np.transpose(anchor_data)
    anchor_data = np.insert(anchor_data, 0, molnames, axis=1)
    np.savetxt(fn, anchor_data, fmt='%s')
    for i in range(nsites):
        fn.write("\nSITE %s FRAGMENTS\n" % str(i+1))
        fragnames = [molnames[i] for i in set(Uniqueness["R%s" % str(i+1)])]
        frags = UniqueRGroupNames['R%s' % str(i+1)]
        frags = [[fragnames[i]] + frags[i] for i in range(len(frags))]
        csv.writer(fn,delimiter=" ",lineterminator="\n").writerows(frags)
    fn.write("END")
    fn.close()

    return molnames[0]

if __name__ == '__main__':
    mol_list = 'mol_list.txt' 
    MCSS_RDecomp(mol_list)
