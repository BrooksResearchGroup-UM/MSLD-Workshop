'''
PyMOL MSLD Prep Plugin
Written by: Luis F Cervantes 05/05/2022

This plugin incorporates the msld-py-prep scripts
by Jonah Z Vilseck, Luis F Cervantes and 
Charles L Brooks III published in:

It automates most of the process for 
creating the CHARMM input files necessary to
MLSD simulations straight from molecules in a 
CSV file or loaded in PyMOL.
'''

from __future__ import absolute_import
from __future__ import print_function

# Avoid importing "expensive" modules here (e.g. scipy), since this code is
# executed on PyMOL's startup. Only import such modules inside functions.

import os

def __init_plugin__(app=None):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('MSLD Prep', run_plugin_gui)


# global reference to avoid garbage collection of our dialog
dialog = None


def run_plugin_gui():
    '''
    Open our custom dialog
    '''
    global dialog

    if dialog is None:
        dialog = make_dialog()

    dialog.show()

def show_warning(warning,info):
    from PyQt5.QtWidgets import QMessageBox
    msg = QMessageBox(dialog)
    msg.setIcon(QMessageBox.Critical)
    msg.setText(warning)
    msg.setInformativeText(info)
    msg.setWindowTitle("Error")
    msg.exec_()

# def run_askMCS(i):
#     text = i.text()
#     if text == "Yes":
#         print(f"A previous MCS has been found in: {savedir}")
#         vis_check.vis_check(mcsout)
#         with open(mcsout,'r') as f:
#             mcsinfo = f.readlines()
#             reflig = [line.split()[1] for line in mcsinfo if line.startswith('REFLIG')]
#             
#         print("Reference Ligand is "+reflig[0])
#         print(f"This previous MCS has been loaded.")
#         
#     if text == "No":
#         
# 
# def askMCS(mcsout):
#     from PyQt5.QtWidgets import QMessageBox
#     msg = QMessageBox(dialog)
#     msg.setIcon(QMessageBox.Question)
#     msg.setText("Do you want to use available MCS?")
#     msg.setInformativeText(f"A previous MCS has been found. If you do not wish to rerun the MCS search, press NO.")
#     msg.setWindowTitle("MCS File Available")
#     msg.setStandardButtons(QMessageBox.Yes|QMessageBox.No|QMessageBox.Cancel)
#     msg.setDefaultButton(QMessageBox.Yes)
#     msg.buttonClicked.connect(lambda: run_askMCS(mcsout))
#     msg.exec_()


def truncate_at_names(fname):
    """
    PyMOL outputs atom names with greater than 4 characters for
    halogens. Need to truncate for CHARMM.
    """

    import string 
    import os
    from copy import deepcopy
    
    with open(fname,'r') as f:
        lines = f.readlines()
    
    newlines = deepcopy(lines)
    idx = [i for i in range(len(lines)) if lines[i].startswith('@<TRIPOS>ATOM') or lines[i].startswith('@<TRIPOS>BOND')]
    
    atnames = [l.rstrip().split()[1] for l in lines[idx[0]+1:idx[1]]]
    atnames = [at for at in atnames]
    alphabet = list(string.ascii_uppercase)
    for i,line in enumerate(lines[idx[0]+1:idx[1]]):
        if any(patt in line for patt in ['Cl','CL','Br','BR']):
            atname = line.rstrip().split()[1]
            if len(atname) <= 4:
                continue
            newatname = atname[0:3]+atname[-1]
            atnames[i] == newatname
            it=0
            while newatname in atnames:
                newatname = atname[0:3]+alphabet[it]
                it+=1
            newline = line.replace(atname, newatname)
            newlines[i+idx[0]+1] = newline
        
    with open(fname,'w') as f:
        f.write("".join(newlines))

def get_bonded_matrix(fname):
    import numpy as np
    with open(fname, 'r') as f:
        x = f.readlines()
    
    indices = [x.index(y) for y in x if y.startswith('@<TRIPOS>ATOM') or y.startswith('@<TRIPOS>BOND')]
    atoms = [line.split()[1] for line in x[indices[0]+1:indices[1]] if line.split()]
    indices = [x.index(y) for y in x if y.startswith('@<TRIPOS>BOND') or y.startswith('@<TRIPOS>SUBSTRUCTURE')]
    bonds = [line.split()[1:-1] for line in x[indices[0]+1:indices[1]] if line.split()]

    n_atoms = len(atoms)
    bonded_matrix = np.zeros((n_atoms,n_atoms))
    
    for bond in bonds:
        bonded_matrix[int(bond[0])-1,int(bond[1])-1] = 1
        bonded_matrix[int(bond[1])-1,int(bond[0])-1] = 1

    return bonded_matrix, atoms


def get_rep_mol_atoms(lines,nsubs,mol_names):
    nsites = len(nsubs)
    rep_mol = mol_names[0]
    all_site_atoms = []
    all_site_names = []
    for site in range(nsites):
        if site == nsites - 1:
            indices = [lines.index(line) for line in lines if line.startswith(f'SITE {site+1}') or line.startswith(f'END')]
        else:
            indices = [lines.index(line) for line in lines if line.startswith(f'SITE {site+1}') or line.startswith(f'SITE {site+2}')]
        site_atoms = [line.split()[1:] for line in lines[indices[0]+1:indices[1]] if line.split()]
        site_names = [line.split()[0] for line in lines[indices[0]+1:indices[1]] if line.split()]
        site_atom_dict = dict(zip(site_names, site_atoms))
        all_site_atoms.append(site_atom_dict)
        all_site_names.append(set(site_names))

    rep_mol = list(set.intersection(*all_site_names))[0]
    site_atoms = []
    for site in range(nsites):
        site_atoms.append(all_site_atoms[site][rep_mol])

    return rep_mol, site_atoms




def write_new_mcsout(to_remove,new_core,site,mcsout='MCS_for_MSLD.txt'):
    import glob
    from copy import deepcopy

    sitee = deepcopy(site)
    with open(mcsout,'r') as f:
        lines = f.readlines()

    # Version control for undo button
    mcsoutfiles = glob.glob(mcsout.replace('.txt','*'))
    if len(mcsoutfiles) > 1:
        mcsoutfiles = [f for f in mcsoutfiles if f != mcsout]
        latest = max([int(fname.split('.txt')[0].split('__')[1]) for fname in mcsoutfiles])
        with open(mcsout.replace(".txt",f"__{latest+1}.txt"),'w') as fnew:
            for line in lines:
                fnew.write(line)
    else:
        with open(mcsout.replace(".txt","__1.txt"),'w') as fnew:
            for line in lines:
                fnew.write(line)

    
    # Modify lines
    indices = [lines.index(line) for line in lines if line.startswith('CORE') or line.startswith('ANCHOR')]
    mol_names = [line.split()[0] for line in lines[indices[0]+1:indices[1]] if line.split()]
    nsubs = [line.split()[1:] for line in lines if line.startswith('NSUBS')][0]
    nsubs = list(map(int,nsubs))
    reflig = [line.split()[1] for line in lines if line.startswith('REFLIG')][0]
    nsites = len(nsubs) 

    indices = [lines.index(line) for line in lines if line.startswith('ANCHOR') or line.startswith('SITE')]
    indices = [lines.index(line) for line in lines if line.startswith('ANCHOR') or line.startswith('SITE')]
    anchor_atoms = [line.split()[1:] for line in lines[indices[0]+1:indices[1]] if line.split()]

    # Find molecule that has unique fragments in all sites
    rep_mol, site_atoms = get_rep_mol_atoms(lines,nsubs,mol_names)
    bonded_matrix, atoms = get_bonded_matrix(rep_mol+'.mol2')
    rep_mol_idx = mol_names.index(rep_mol)

    # This is all to change anchor atoms based on core atoms modified
    rep_mol, site_atoms = get_rep_mol_atoms(lines,nsubs,mol_names)
    rep_mol_idx = mol_names.index(rep_mol)
    bonded_matrix, atoms = get_bonded_matrix(rep_mol+'.mol2')
    all_site_atoms = [inner for outer in site_atoms for inner in outer]
    rep_mol_core = [at for at in atoms if at not in all_site_atoms]
    for site in range(nsites):
        anchor_atom = anchor_atoms[rep_mol_idx][site]
        if anchor_atom in to_remove[rep_mol_idx]:
            # find new anchor atom based on connectivity
            anchor_atom_idx = atoms.index(anchor_atom)
            bonded = bonded_matrix[anchor_atom_idx,:]
            bonded = [atoms[i] for i in range(len(bonded)) if bonded[i] ==1]
            site_atoms_bonded = set(bonded).intersection(set(site_atoms[site]))
            core_atoms_bonded = list(set(bonded).intersection(set(rep_mol_core)))
            heavy_core_atoms_bonded = [at for at in core_atoms_bonded if not at.startswith('H')]
            print(heavy_core_atoms_bonded)
            # if there are two core atoms that connect to a site n atom, then new_anchor_atom is 'DUM'
            # else: the new_anchor_atom is the only core atom that connects to one or multiple site n atoms.
            if len(heavy_core_atoms_bonded) == 1:
                new_anchor_atom = heavy_core_atoms_bonded[0]
                new_anchor_atom_idx = new_core[rep_mol_idx].index(new_anchor_atom)
                for i in range(len(mol_names)):
                    anchor_atoms[i][site] = new_core[i][new_anchor_atom_idx]
            else:
                new_anchor_atom = 'DUM'
                for i in range(len(mol_names)):
                    anchor_atoms[i][site] = new_anchor_atom
        elif anchor_atom == 'DUM':
            # Check to see if only one core atom connects to any core atoms. If yes, then
            # This is the new anchor atom. Otherwise, still 'DUM'.
            core_ats = new_core[rep_mol_idx]
            site_ats = site_atoms[site]
            bonded_ats = [bonded_matrix[atoms.index(at)] for at in site_ats]
            bonded_ats = list(map(lambda bonded_ats: [atoms[i] for i in range(len(bonded_ats)) if bonded_ats[i] == 1],bonded_ats))
            bonded_cores = list(map(lambda bonded_ats: set(bonded_ats).intersection(set(core_ats)),bonded_ats))
            bonded_cores = [ats for ats in bonded_cores if ats]
            if all(ats == bonded_cores[0] for ats in bonded_cores) and len(list(bonded_cores[0])) == 1:
                new_anchor_atom = list(bonded_cores[0])[0]
                new_anchor_atom_idx = new_core[rep_mol_idx].index(new_anchor_atom)
                for i in range(len(mol_names)):
                    anchor_atoms[i][site] = new_core[i][new_anchor_atom_idx]
            # Else the anchor atom is still 'DUM'

    # Get site molnames and atoms
    all_site_atoms = []
    all_site_names = []
    for site in range(nsites):
        if site == nsites - 1:
            indices = [lines.index(line) for line in lines if line.startswith(f'SITE {site+1}') or line.startswith(f'END')]
        else:
            indices = [lines.index(line) for line in lines if line.startswith(f'SITE {site+1}') or line.startswith(f'SITE {site+2}')]
        site_atoms = [line.split()[1:] for line in lines[indices[0]+1:indices[1]] if line.split()]
        site_names = [line.split()[0] for line in lines[indices[0]+1:indices[1]] if line.split()]
        site_atom_dict = dict(zip(site_names, site_atoms))
        all_site_atoms.append(site_atom_dict)
        all_site_names.append(site_names)
    
    print(all_site_atoms)
    toRemove = dict(zip(mol_names,to_remove))
    mols = list(all_site_atoms[sitee].keys())
    for mol in toRemove.keys():
        if mol in all_site_atoms[sitee].keys():
            all_site_atoms[sitee][mol].extend(toRemove[mol])

    nsubs = list(map(str,nsubs))

    to_write =f"""# Maximum Common Substructure Search for Multisite Lambda Dynamics (JV 2022)
# {len(mols)} molecules processed

NSUBS {' '.join(nsubs)}

REFLIG {reflig}

"""

     
    CORESPECS = "CORE\n"
    for i in range(len(mol_names)):
        CORESPECS += f"{mol_names[i]} {' '.join(new_core[i])}\n"
    to_write += CORESPECS

    ANCHORSPECS = '\nANCHOR ATOMS\n'
    for i in range(len(mol_names)):
        ANCHORSPECS += f"{mol_names[i]} {' '.join(anchor_atoms[i])}\n"
    to_write += ANCHORSPECS

    for site in range(nsites):
        SITESPEC = f"\nSITE {site+1} FRAGMENTS\n"
        for mol in all_site_atoms[site].keys():
            SITESPEC += f"{mol} {' '.join(all_site_atoms[site][mol])}\n"
        to_write += SITESPEC

    to_write += '\nEND'

    with open(mcsout,'w') as f:
        f.write(to_write)

    return None
     


def remove_from_core(at_ind, mol_ind,mol_selection_core_atoms, name, core_atoms,to_site=1, mcsout="MCS_for_MSLD.txt"):
    from .msld_py_prep import vis_check
    
    # Remove specified atom
    to_remove = [core_atoms[mol_ind][at_ind]]
    mol_selection_core_atoms = [x for x in mol_selection_core_atoms if x != to_remove[0]]
    to_move = list(map(lambda mol: [mol[i] for i in range(len(mol)) if i == at_ind],core_atoms))
    new_core_atoms = list(map(lambda mol: [mol[i] for i in range(len(mol)) if i != at_ind],core_atoms))
    
    # Remove hydrogens bound to heavy atom, if any
    bonded_matrix, atoms = get_bonded_matrix(f"{name}.mol2")
    hydrogen_indices = [ind for ind in range(len(atoms)) if atoms[ind].startswith('H')]

    if not hydrogen_indices:
        write_new_mcsout(to_move,new_core_atoms, to_site,mcsout)
        vis_check.vis_check(mcsout)
        return None


    h_at = []
    for h in hydrogen_indices:
        heavy_bonded = bonded_matrix[h,:].nonzero()[0]
        heavy_bonded = [atoms[i] for i in heavy_bonded][0]
        if heavy_bonded in to_remove: 
            print(f'{atoms[h]} will also be removed since it is bound to {heavy_bonded}')
            h_at.append(atoms[h])

    ind = [at_ind for at_ind in range(len(mol_selection_core_atoms)) if mol_selection_core_atoms[at_ind] in h_at]
    to_hmove = list(map(lambda mol: [mol[i] for i in range(len(mol)) if i in ind],new_core_atoms))
    new_core_atoms = list(map(lambda mol: [mol[i] for i in range(len(mol)) if i not in ind],new_core_atoms))
    for i in range(len(to_move)):
        for j in range(len(to_hmove[i])):
            to_move[i].append(to_hmove[i][j])

    write_new_mcsout(to_move,new_core_atoms,to_site,mcsout)
    vis_check.vis_check(mcsout)

    return None
      
def make_dialog():
    # entry point to PyMOL's API
    from pymol import cmd

    # pymol.Qt provides the PyQt5 interface, but may support PyQt4
    # and/or PySide as well
    from pymol.Qt import QtWidgets
    from pymol.Qt.utils import loadUi
    from pymol.Qt.utils import getSaveFileNameWithExt
    from PyQt5.QtWidgets import QFileDialog, QMessageBox
    import subprocess

    # create a new Window
    dialog = QtWidgets.QDialog()

    # populate the Window from our *.ui file which was created with the Qt Designer
    uifile = os.path.join(os.path.dirname(__file__), 'crnwidget.ui')
    form = loadUi(uifile, dialog)

    # populate "Align to" dropdown with loadedPyMOL objects
    form.input_objecttoalign.addItem("")
    for pymolobj in cmd.get_names('objects'):
        form.input_objecttoalign.addItem(pymolobj)    

    # callback for the "Browse" button
    def browse_filename():
        filename,_ = QFileDialog.getOpenFileName(
            dialog, 'Open SMILES Dataset', filter='CSV File (*.csv)')
        if filename:
            form.input_dataset.setText(filename)

    # Align loaded molecules
    def alignMols(mols, wd, target):
        from . import generate_3D
        import pandas as pd
        from .msld_py_prep import rename_atoms
        from rdkit import Chem

        file_formats = ['sdf','mol','xyz']
        loadingFunctions = {'sdf': generate_3D.load_sdf, 'mol': generate_3D.load_mol,\
                            'mol2': generate_3D.load_mol2}
        molsext = list(map(lambda x: x+".sdf",mols))
        molObjList = []
        for fn,mol in zip(molsext,mols):
            try:
                molObjList.append(loadingFunctions['sdf'](fn)[0])
            except:
                cmd.save(os.path.join(wd,mol+'.sdf'),selection=mol,format='sdf')
                molObjList.append(loadingFunctions['sdf'](fn)[0])

        smileslist = list(map(Chem.MolToSmiles,molObjList))
        d = {'id': mols,'smiles': smileslist}
        df = pd.DataFrame(data=d)
        datasetName = 'generated_mol_dataset.csv' 
        df.to_csv(os.path.join(wd,datasetName),index=False)

        for obj in mols:
            cmd.delete(obj)

        generate3DMols(datasetName, wd,target=target+'.sdf')

    # 3D conformer generation helper functions
    def generate3DMols(fn, wd, smiles_col=1, names_col=0, delim=',',pattern=None,target=None):
        from . import generate_3D
        from .msld_py_prep import rename_atoms
 
        originalPath = os.getcwd()
        if wd:
            os.chdir(wd)
        mols = generate_3D.RDKit_Tools(fn,smiles_col=smiles_col,names_col=names_col,delim=delim)
        
        generated_mols = mols.generate_3D(ref=target)
        print(f"Molecules have been generated and saved in {wd}")
        # print(generated_mols) 
        for molname in generated_mols:
            # Load mol2 files into pymol
            cmd.load(f"{molname}.mol")
            # Convert to mol2 format
            cmd.save(f"{molname}.mol2",selection=molname,format="mol2")
            # Change atom naming to 4 characters
            truncate_at_names(f"{molname}.mol2")
            # Remove mol objects from PyMOL session
            cmd.delete(molname)
            # Uniquely rename atoms (required for PyPrep Scripts)
            rename_atoms.rename_atoms(molname)  
            # Load mol2 files into pymol
            cmd.load(f"{molname}.mol2")
            # Save as sdf format as well 
            cmd.save(f"{molname}.sdf",selection=molname,format="sdf")

        if wd:
            os.chdir(originalPath)
         



    # callback for the "Align" button
    def run_align():
        # get form data
        target = form.input_objecttoalign.currentText()
        pattern = form.input_patterntoalign.text().rstrip()
        dataset = form.input_dataset.text()
 
        # Verify alignment target is specified 
        if not (target or pattern):
            show_warning('No Molecule Selection!','Please specify the\
                        molecule you wish to align to.')
            return None

        if pattern and not (dataset or cmd.get_names('objects')):
            show_warning('No Molecules to Align!', 'Please load or specify\
                        the molecule(s) you wish to align.')
            return None

        savedir = QFileDialog.getExistingDirectory(
            dialog, caption='Save Aligned Molecules to...')

        if not savedir:
            return None
        # print(savedir)

        # Generate 3D conformers and align using RDKit if SMILES csv is
        # specified
        if dataset:
            if target:
                try:
                    os.remove(os.path.join(savedir,f"{target}.sdf"))
                    cmd.save(os.path.join(savedir,f"{target}.sdf"),selection=target,format='sdf')
                except:
                    cmd.save(os.path.join(savedir,f"{target}.sdf"),selection=target,format='sdf')

                generate3DMols(dataset,wd=savedir,target=os.path.join(savedir,f"{target}.sdf"))
                if pattern:
                    # Could modify scripts to add a substructure pattern to align to
                    pass
                else:
                    pass
            else:
                pass
 
        # Otherwise assume all molecule queries are loaded in PyMOL session
        else:
            if target:
                if not os.path.exists(os.path.join(savedir,f"{target}.sdf")):
                    cmd.save(os.path.join(savedir,f"{target}.sdf"),selection=target,format='sdf')

                molecules = cmd.get_names('objects')                
                for mol in molecules:
                    if not os.path.exists(os.path.join(savedir,f"{mol}.sdf")):
                        cmd.save(f"{mol}.sdf",selection=mol,format="sdf")
                alignMols(molecules, wd=savedir,target=target)  
    
    def displayTxt(mcsout):
        """

        """
        form.display_mcssoutput.clear()
        
        with open(mcsout,'r') as f:
            displaytxt = f.readlines()
 
        form.display_mcssoutput.append("".join(displaytxt))

    def run_mcss():
        from .msld_py_prep import msld_chk
        from .msld_py_prep import msld_mcs
        from .msld_py_prep import vis_check
        import glob
        
        # Get form data
        sysname = form.input_sysname.text().rstrip()
        forcefield = form.input_forcefield.currentText()

        # Verify system name is specified 
        if not sysname:
            show_warning('No System Name!','Please specify the\
                        name of your system.')
            return None

        # Verify force field is specified 
        if not forcefield:
            show_warning('No ForceField!','Please specify the\
                        force field with which you wish to parameterize your system.')
            return None

        # Have not enabled GAFF compatibility yet 
        if forcefield == 'GAFF':
            show_warning('GAFF Not Supported','GAFF is not yet supported.\
                        Please use another force field.')

            # Comment out warning and generate CHARMM compatible GAFF
            # parameter files here if this is needed.

            return None

        # Get pymol objects
        molObjs = cmd.get_names('objects', enabled_only=1)
        
        # Verify more than one molecule is loaded
        if len(molObjs) <= 1:
            show_warning('No Molecules Loaded!','Please load the\
                        molecules you wish to do an MCS Search for.')
            return None

        savedir = QFileDialog.getExistingDirectory(
            dialog, caption='Save Aligned Molecules to...')

        cwd = os.getcwd()
        os.chdir(savedir) 
        


        for ob in molObjs:
            # if there are no rtf files for pymol objects, then try to create them
            if not os.path.exists(ob+'.str'):
                return_code = subprocess.call(f"cgenff {ob+'.mol2'} > {ob+'.str'}",shell=True)

                if return_code == 127:
                    subprocess.run([f"rm {ob+'.str'}"],shell=True,capture_output=True).stdout.decode('utf-8')
                    load_cgenff_output = subprocess.run([f"module load cgenff"],shell=True,capture_output=True).stdout.decode('utf-8')
                    print(f"cgenff loading output is {load_cgenff_output}")
                    which_bash = subprocess.run([f"which bash"],shell=True,capture_output=True).stdout.decode('utf-8')
                    print(f"which bash is: {which_bash}")
                    return_code = subprocess.call(f"cgenff {ob+'.mol2'} > {ob+'.str'}",shell=True)

                    if return_code == 127:
                        show_warning("ParamChem was not found!", "Please install\
                                    the ParamChem module and alias the tool as 'cgenff'.\
                                    Otherwise create CHARMM stream files (.str) containing\
                                    both rtf and prm information.")
                        return None
 
            cmd.save(ob+'.mol2',selection=ob,format='mol2')
            truncate_at_names(ob+'.mol2')

        # Create molfile
        molfile = "mol_list.txt"                  # list of mol2 file names
        molfiled = "\n".join(molObjs)
        with open(molfile,'w') as f:
            f.write(molfiled)      

        # Specify output files 
        mcsout = 'MCS_for_MSLD.txt'               # MCS output filename
        outdir = 'build.'+sysname                 # MSLD output directory
      
        cgenff = forcefield == 'CGenFF'           # Are CGenFF/ParamChem parameters being used?

        if len(glob.glob(mcsout)) == 0:
            ## (2) Check molfile and toppar files before getting started
            msld_chk.MsldCHK(molfile)
            print("chk finished")
            os.chdir(savedir)
            reflig = msld_mcs.MsldMCS(molfile,mcsout,cutoff=0.8,debug=False)
            if not reflig:
                return None
            print("MCS results printed to "+os.path.join(savedir,mcsout))
            print("Reference Ligand is "+reflig)
            vis_check.vis_check(mcsout)

        else:
            print(f"A previous MCS has been found in: {savedir}")
            vis_check.vis_check(mcsout)
            with open(mcsout,'r') as f:
                mcsinfo = f.readlines()
                reflig = [line.split()[1] for line in mcsinfo if line.startswith('REFLIG')]
                
            print("Reference Ligand is "+reflig[0])
            print(f"This previous MCS has been loaded.")

        os.chdir(cwd)

        displayTxt(mcsout)
 
    def run_rdkitmcss():
        from .msld_py_prep import msld_chk
        from .msld_py_prep import vis_check
        from .msld_py_prep import msld_mcs_rdecomp
        import glob
        
        # Get form data
        sysname = form.input_sysname.text().rstrip()
        forcefield = form.input_forcefield.currentText()

        # Verify system name is specified 
        if not sysname:
            show_warning('No System Name!','Please specify the\
                        name of your system.')
            return None

        # Verify force field is specified 
        if not forcefield:
            show_warning('No ForceField!','Please specify the\
                        force field with which you wish to parameterize your system.')
            return None

        # Have not enabled GAFF compatibility yet 
        if forcefield == 'GAFF':
            show_warning('GAFF Not Supported','GAFF is not yet supported.\
                        Please use another force field.')

            # Comment out warning and generate CHARMM compatible GAFF
            # parameter files here if this is needed.

            return None

        # Get pymol objects
        molObjs = cmd.get_names('objects', enabled_only=1)
        
        # Verify more than one molecule is loaded
        if len(molObjs) <= 1:
            show_warning('No Molecules Loaded!','Please load the\
                        molecules you wish to do an MCS Search for.')
            return None

        savedir = QFileDialog.getExistingDirectory(
            dialog, caption='Save Aligned Molecules to...')

        cwd = os.getcwd()
        os.chdir(savedir) 
        


        for ob in molObjs:
            # if there are no rtf files for pymol objects, then try to create them
            if not os.path.exists(ob+'.str'):
                return_code = subprocess.call(f"cgenff {ob+'.mol2'} > {ob+'.str'}",shell=True)

                if return_code == 127:
                    subprocess.run([f"rm {ob+'.str'}"],shell=True,capture_output=True).stdout.decode('utf-8')
                    load_cgenff_output = subprocess.run([f"module load cgenff"],shell=True,capture_output=True).stdout.decode('utf-8')
                    print(f"cgenff loading output is {load_cgenff_output}")
                    which_bash = subprocess.run([f"which bash"],shell=True,capture_output=True).stdout.decode('utf-8')
                    print(f"which bash is: {which_bash}")
                    return_code = subprocess.call(f"cgenff {ob+'.mol2'} > {ob+'.str'}",shell=True)

                    if return_code == 127:
                        show_warning("ParamChem was not found!", "Please install\
                                    the ParamChem module and alias the tool as 'cgenff'.\
                                    Otherwise create CHARMM stream files (.str) containing\
                                    both rtf and prm information.")
                        return None
 
            cmd.save(ob+'.mol2',selection=ob,format='mol2')
            truncate_at_names(ob+'.mol2')

        # Create molfile
        molfile = "mol_list.txt"                  # list of mol2 file names
        molfiled = "\n".join(molObjs)
        with open(molfile,'w') as f:
            f.write(molfiled)      

        # Specify output files 
        mcsout = 'MCS_for_MSLD.txt'               # MCS output filename
        outdir = 'build.'+sysname                 # MSLD output directory
      
        cgenff = forcefield == 'CGenFF'           # Are CGenFF/ParamChem parameters being used?

        if len(glob.glob(mcsout)) == 0:
            ## (2) Check molfile and toppar files before getting started
            msld_chk.MsldCHK(molfile)
            print("chk finished")
            os.chdir(savedir)
            reflig = msld_mcs_rdecomp.MCSS_RDecomp(molfile,mcsout=mcsout)
            if not reflig:
                return None
            print("MCS results printed to "+os.path.join(savedir,mcsout))
            print("Reference Ligand is "+reflig)
            vis_check.vis_check(mcsout)

        else:
            print(f"A previous MCS has been found in: {savedir}")
            vis_check.vis_check(mcsout)
            with open(mcsout,'r') as f:
                mcsinfo = f.readlines()
                reflig = [line.split()[1] for line in mcsinfo if line.startswith('REFLIG')]
                
            print("Reference Ligand is "+reflig[0])
            print(f"This previous MCS has been loaded.")

        os.chdir(cwd)
        
        displayTxt(mcsout)

    def run_crn():
        from .msld_py_prep import msld_crn
        from .msld_py_prep import msld_prm
        from .msld_py_prep import msld_wrt
        import glob
        from pymol import cmd

        # Get form data
        sysname = form.input_sysname.text().rstrip()
        forcefield = form.input_forcefield.currentText()

        # Verify system name is specified 
        if not sysname:
            show_warning('No System Name!','Please specify the\
                        name of your system.')
            return None

        # Verify force field is specified 
        if not forcefield:
            show_warning('No ForceField!','Please specify the\
                        force field with which you wish to parameterize your system.')
            return None

        # Have not enabled GAFF compatibility yet 
        if forcefield == 'GAFF':
            show_warning('GAFF Not Supported','GAFF is not yet supported.\
                        Please use another force field.')

            # Comment out warning and generate CHARMM compatible GAFF
            # parameter files here if this is needed.

            return None

        # Ask for savedir
        savedir = QFileDialog.getExistingDirectory(
            dialog, caption='Select MCS directory...')

        if not savedir:
            return None


        # Check for MCS_for_MSLD.txt file in savedir
        if not os.path.exists(os.path.join(savedir,'MCS_for_MSLD.txt')):
            show_warning("MCS Not Found!","Please select a directory with an \
                        MCS file present.")
            return None
        
        outdir = 'build.'+sysname       # MSLD output directory
        mcsout = 'MCS_for_MSLD.txt'     # MCS output filename
        cgenff = forcefield == 'CGenFF' # Are CGenFF/ParamChem parameters being used?
        inFrag=[[]]                     # reflig core atoms to include in each frag/site (list of nsub lists)
        AnCore=[[]]                     # anchor atoms at each site to include in the core (list of nsub lists)

        cwd = os.getcwd()
        os.chdir(savedir) 

        msld_crn.MsldCRN(mcsout,outdir,inFrag,AnCore,ChkQChange=True,verbose=True,debug=False)
        
        
        #####################################################################
        ## (5) Write Ligand Parameters & the Charmm ALF input scripts
        msld_prm.MsldPRM(outdir,cgenff,verbose=False,debug=False)
        msld_wrt.writeALF_Files(sysname,outdir,cgenff)


        ## Final Notes to the user
        print("default TOPPAR parameters copied into build."+sysname+". \
              Check to make sure these work for your system!")
        
        
        # Unload current molecule files from PyMOL
        with open(mcsout,'r') as f:
            x = f.readlines()
        indices = [x.index(y) for y in x if y.startswith('CORE') or y.startswith('ANCHOR')]
        mol_names = [line.split()[0] for line in x[indices[0]+1:indices[1]] if line.split()]
        loaded = cmd.get_names('objects')
        for obj in mol_names:
            if obj in loaded:
                cmd.delete(obj)

        # Load files from build.{sysname} directory
        file_list = glob.glob(os.path.join(savedir,f'build.{sysname}','site*pdb'))
        file_list.append(os.path.join(savedir,f'build.{sysname}','core.pdb'))
        for name in file_list:
            cmd.load(name)

        os.chdir(cwd) 

    def run_moveselectedtosite():
        # Get form data
        site = form.input_siteselectedtomove.text().rstrip()
        if not site:
            show_warning("No Site Selected!","Please specify the site number you wish to move the selected atom to.")
            return None

        site = int(site)
        print(f"Will be moving to site {site}")
 
        mcsout = 'MCS_for_MSLD.txt'

        if not os.path.exists(mcsout):
            show_warning('MCS Not Found!','An MCS file has not been found in the current working directory. Please change directories or generate an MCS file.')
            return None

        atom_selection = cmd.get_model('pk1').atom[0].name
        mol_selection = cmd.get_names('objects',selection='pk1')[0]
        
        with open(mcsout,'r') as f:
            x = f.readlines()

        indices = [x.index(y) for y in x if y.startswith('CORE') or y.startswith('ANCHOR')]
        mol_names = [line.split()[0] for line in x[indices[0]+1:indices[1]] if line.split()]
        core_atoms = [line.split()[1:] for line in x[indices[0]+1:indices[1]] if line.split()]
        mol_selection_core_atoms = core_atoms[mol_names.index(mol_selection)]
        ind = [at_ind for at_ind in range(len(mol_selection_core_atoms)) if mol_selection_core_atoms[at_ind] == atom_selection]

        if not ind:
            show_warning("Core Atom Not Selected!","Please make sure the atom you selected is a core atom.")
            return None
        ind = ind[0]
        mol_ind = mol_names.index(mol_selection)
        remove_from_core(ind,mol_ind, mol_selection_core_atoms, mol_selection,core_atoms,to_site=site-1)
 
        displayTxt(mcsout)

        return None 


    def run_undomove():
        import glob
        from .msld_py_prep import vis_check

        mcsout = 'MCS_for_MSLD.txt'
        # Version control for undo button
        mcsoutfiles = glob.glob(mcsout.replace('.txt','*'))
        if len(mcsoutfiles) > 1:
            mcsoutfiles = [f for f in mcsoutfiles if f != mcsout]
            latest = max([int(fname.split('.txt')[0].split('__')[1]) for fname in mcsoutfiles])
            with open(mcsout.replace(".txt",f"__{latest}.txt"),'r') as ff:
                redolines = ff.readlines()
            with open(mcsout,'w') as f_redo:
                for line in redolines:
                    f_redo.write(line)
            os.remove(mcsout.replace(".txt",f"__{latest}.txt"))
            vis_check.vis_check(mcsout)
        else:
            show_warning("Already at oldest change.","Cannot undo any further.")

        displayTxt(mcsout)

    # hook up button callbacks
    form.button_align.clicked.connect(run_align)
    form.button_mcss.clicked.connect(run_mcss)
    form.button_rdkitmcss.clicked.connect(run_rdkitmcss)
    form.button_crn.clicked.connect(run_crn)
    form.button_moveselectedtosite.clicked.connect(run_moveselectedtosite)
    form.button_undomove.clicked.connect(run_undomove)
    form.button_browse.clicked.connect(browse_filename)
    form.button_close.clicked.connect(dialog.close)

    return dialog
