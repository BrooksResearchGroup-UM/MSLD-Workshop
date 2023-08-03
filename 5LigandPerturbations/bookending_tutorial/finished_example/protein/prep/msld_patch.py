#
# Load different kinds of alchemical patches for MSLD
# Alchemical perturbations are loaded by "type"
# Available types:
#	- ligand
#	- protein side chains
#
# Functions for each type are in:
#	- load_alchem_toppar
#	- load_alchem_patches
#

def load_alchem_toppar(nsites,nsubs,pert,builddir):
    """ Reads through pert, identifies perturbation types
        and loads the corresponding toppar files """

    import pycharmm
    import pycharmm.read as read
    for site in range(nsites):
        # TYPE: LIGAND 
        if pert[site]['type'] == 'ligand':
            # for ligands generated with msld-py-prep
            read.rtf(builddir+'/core.rtf', append=True)
            read.prm(builddir+'/full_ligand.prm', append=True, flex=True)

        # TYPE: PROTEIN SIDE CHAIN
        elif pert[site]['type'] == 'side_chain':
            # for protein side chains
            aastream=[]
            for site in range(nsites):
                for sub in range(1,nsubs[site]): # skip loading "nat" patch (for prot only)
                    if not (pert[site]['subs'][sub] in aastream):
                        aastream.append(pert[site]['subs'][sub])
                        pycharmm.lingo.charmm_script('stream '+builddir+'/aa_stream/ca_'+pert[site]['subs'][sub]+'.str')

        else:
            print("Perturbation Type NOT recognized - check and resubmit")
            pycharmm.lingo.charmm_script('stop')

    return


def load_alchem_patches(nsites,nsubs,pert,builddir):
    """ Load perturbation patches for MSLD """

    import pycharmm
    import pycharmm.read as read
    import pycharmm.select as select
    for site in range(nsites):
        # TYPE: LIGAND 
        if pert[site]['type'] == 'ligand':
            # load in patches for alchem ligand
            pycharmm.lingo.charmm_script('ic generate')
            for site in range(nsites):
                for sub in range(nsubs[site]):
                    read.rtf(builddir+'/site'+str(site+1)+'_sub'+pert[site]['subs'][sub]+'_pres.rtf', append=True)
                    pycharmm.lingo.charmm_script('''
patch p{}_{} {} {} setup
read coor pdb resid name {}
ic param
ic build
'''.format(str(site+1),
           pert[site]['subs'][sub],
           pert[site]['segid'],
           pert[site]['resid'],
           builddir+'/site'+str(site+1)+'_sub'+pert[site]['subs'][sub]+'_frag.pdb')
)

            # read in LonePair sites (if applicable)
            pycharmm.lingo.charmm_script('stream '+builddir+'/lpsites.inp')
            
            # define ligand substituent selections
            for site in range(nsites):
                pert[site]['select']=[]
                for sub in range(nsubs[site]):
                    # append charmm variable name for substituent selection
                    pert[site]['select'].append('site'+str(site+1)+'_sub'+pert[site]['subs'][sub])
                    # extract alchem patch atoms from patch file
                    sub_atoms=''
                    rtffile = builddir+'/site'+str(site+1)+'_sub'+pert[site]['subs'][sub]+'_pres.rtf'
                    for line in open(rtffile,'r'):
                        if line[0:4] == 'ATOM': sub_atoms=sub_atoms+line.split()[1].upper()+' '
                    sub_atoms=sub_atoms[:-1] # remove trailing space
                    atoms_in_sub = pycharmm.SelectAtoms().by_res_and_type(pert[site]['segid'],pert[site]['resid'],sub_atoms)
                    select.store_selection(pert[site]['select'][sub],atoms_in_sub) 
                    # saves the identified atoms in the charmm variable


        # TYPE: PROTEIN SIDE CHAIN
        elif pert[site]['type'] == 'side_chain':
            # load in patches for alchem (protein) side chains
            for site in range(nsites):
                for sub in range(1,nsubs[site]): # skip loading "nat" patch (for prot only)
                    pycharmm.lingo.charmm_script('''
patch ca_{} {} {} setup
! read coor pdb resid name pdb.pdb
ic param
ic build
'''.format(pert[site]['subs'][sub],
           pert[site]['segid'],
           pert[site]['resid'])
)

            # set alchem selection definitions (prot)
            backbone_atoms = 'N HN CA HA C O HT1 HT2 HT3 OXT'
            
            # figure out the "native" selection; also loads patch selection stream files
            for site in range(nsites):
            
                # native sub is all residue atoms minus backbone
                atoms_in_resid = pycharmm.SelectAtoms(seg_id=pert[site]['segid']) & pycharmm.SelectAtoms(res_id=pert[site]['resid'])
                atoms_in_bkbn  = pycharmm.SelectAtoms().by_res_and_type(pert[site]['segid'],pert[site]['resid'],backbone_atoms)
                atoms_in_nat   = atoms_in_resid & ~ atoms_in_bkbn
            
                #print('atoms_in_resid',atoms_in_resid.get_atom_types())   # get atom names of atoms selected
                #print('atoms_in_bkbn ',atoms_in_bkbn.get_atom_types())   # get atom names of atoms selected
                #print('atoms_in_nat  ',atoms_in_nat.get_atom_types())   # get atom names of atoms selected
            
                # create list of selection names
                patch_nat = 'seg'+pert[site]['segid']+'site'+pert[site]['resid']+'sub'+pert[site]['subs'][0]
                pert[site]['select']=[patch_nat]
            
                for sub in range(1,nsubs[site]): # skip loading "nat" patch (for prot only)
                    # load alchemical patch; subtract patch from native
                    patch_sub = 'seg'+pert[site]['segid']+'site'+pert[site]['resid']+'sub'+pert[site]['subs'][sub]
                    pert[site]['select'].append(patch_sub)
            
                    # extract alchem patch atoms from patch file
                    sub_atoms=''
                    for line in open(builddir+'/aa_stream/ca_'+pert[site]['subs'][sub]+'.str','r'):
                        if line[0:4] == 'ATOM': sub_atoms=sub_atoms+line.split()[1]+' '
                    sub_atoms=sub_atoms[:-1] # remove trailing space
                    atoms_in_sub = pycharmm.SelectAtoms().by_res_and_type(pert[site]['segid'],pert[site]['resid'],sub_atoms)
                    select.store_selection(pert[site]['select'][sub],atoms_in_sub)
                    atoms_in_nat = atoms_in_nat & ~ atoms_in_sub
            
                # save native selection
                select.store_selection(pert[site]['select'][0],atoms_in_nat)

        else:
            print("Perturbation Type NOT recognized - check and resubmit")
            pycharmm.lingo.charmm_script('stop')

    return






