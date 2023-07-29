import numpy as np
import sys
import string
import os
from functools import reduce
import shutil
from itertools import product
import pymbar

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


topparDir = '../../prep/toppar/'
nbondDir = '../../'
prepDir = '../../prep/'
box = 94 
trajname = sys.argv[1]


def is_factor(n):
    if (n % 2 != 0): return False  # favors even number
    while n:
        flag = False
        for x in (2,3,5):
            if n % x == 0:
               n = n / x
               flag = True
               break

        if flag: continue
        break

    if n == 1: return True
    return False

def checkfft(n, margin = 5):
    n = int(n) + margin
    while 1:
        if is_factor(n): break
        else: n += 1
    return n

def run_md(trajname,useomm=False,useblade=False,nequil=1000,nsteps=5000,nsavc=100,leap=True,lang=True,temp=298.15):
    dyn.set_fbetas(np.full((psf.get_natom()),10.0,dtype=float))
   
    res_file = pycharmm.CharmmFile(file_name=f'res/{trajname}.res', file_unit=2,
                                   formatted=True,read_only=False)
    my_dyn = pycharmm.DynamicsScript(leap=leap, lang=lang, start=True,
                                     nstep=nequil, timest=0.002,
                                     firstt=temp, finalt=temp, tbath=temp,
                                     tstruc=temp,
                                     teminc=0.0, twindh=0.0, twindl=0.0,
                                     iunwri=res_file.file_unit,
                                     inbfrq=-1, imgfrq=-1,
                                     iasors=0, iasvel=1, ichecw=0, iscale=0,
                                     iscvel=0,echeck=-1, nsavc=0, nsavv=0, nsavl=0, ntrfrq=0,
                                     isvfrq=nsavc,
                                     iprfrq=2*nsavc, nprint=nsavc, ihtfrq=0, ieqfrq=0,
                                     ilbfrq=0,ihbfrq=0,
                                     omm=useomm, blade=useblade)
    my_dyn.run()

    res_file.close()
    # open unit 2 write form name res/{}.res
    res_file = pycharmm.CharmmFile(file_name=f'res/{trajname}.res', file_unit=2,
                                   formatted=True,read_only=False)
    # open unit 1 write file name dcd/{}.dcd
    dcd_file = pycharmm.CharmmFile(file_name=f'dcd/{trajname}.dcd', file_unit=1,
                                   formatted=False,read_only=False)

    my_dyn = pycharmm.DynamicsScript(leap=leap, lang=lang, start=False, restart = True,
                                     nstep=nsteps, timest=0.002,
                                     firstt=temp, finalt=temp, tbath=temp,
                                     tstruc=temp,
                                     teminc=0.0, twindh=0.0, twindl=0.0,
                                     iunwri=res_file.file_unit,
                                     iunrea=res_file.file_unit,
                                     iuncrd=dcd_file.file_unit,
                                     inbfrq=-1, imgfrq=-1,
                                     iasors=0, iasvel=1, ichecw=0, iscale=0,
                                     iscvel=0,echeck=-1, nsavc=nsavc, nsavv=0, nsavl=0, ntrfrq=0,
                                     isvfrq=nsavc,
                                     iprfrq=2*nsavc, nprint=nsavc, ihtfrq=0, ieqfrq=0,
                                     ilbfrq=0,ihbfrq=0,
                                     omm=useomm, blade=useblade)
    my_dyn.run()

    res_file.close()
    dcd_file.close()
    return


def streamCGenFF(topparDir):
    """

    """

    stream = f"""! protein topology and parameter
open read card unit 10 name {os.path.join(topparDir,'top_all36_prot.rtf')}
read  rtf card unit 10

open read card unit 20 name {os.path.join(topparDir,'par_all36m_prot.prm')}
read para card unit 20 flex

! nucleic acids
open read card unit 10 name {os.path.join(topparDir,'top_all36_na.rtf')}
read  rtf card unit 10 append

open read card unit 20 name {os.path.join(topparDir,'par_all36_na.prm')}
read para card unit 20 append flex

! carbohydrates
open read card unit 10 name {os.path.join(topparDir,'top_all36_carb.rtf')}
read  rtf card unit 10 append

open read card unit 20 name {os.path.join(topparDir,'par_all36_carb.prm')}
read para card unit 20 append flex

! CGenFF
open read card unit 10 name {os.path.join(topparDir,'top_all36_cgenff.rtf')}
read  rtf card unit 10 append

bomblev -5
open read card unit 20 name {os.path.join(topparDir,'par_all36_cgenff.prm')}
read para card unit 20 append flex
bomblev 0

! Water
stream {os.path.join(topparDir,'toppar_water_ions.str')}"""
    lingo.charmm_script(stream)

# Set-up short dynamics
if not os.path.isdir('res'): os.system('mkdir res')
if not os.path.isdir('dcd'): os.system('mkdir dcd')

# Load Force Field
lingo.charmm_script('set temp = 298.15')
streamCGenFF(topparDir)
lingo.charmm_script(f"read para card flex append name {os.path.join(prepDir,'full_ligand.prm')}")

lingo.charmm_script("ioformat extended")

# Load system
read.psf_card('system.psf')
lingo.charmm_script("open read card unit 10 name system.crd\n\
read coor card unit 10")

# Mimic MSLD conditions
# Set up PBC
lingo.charmm_script(f"""coor stat
crystal define cubic {box} {box} {box} 90. 90. 90.
crystal build cutoff 14 nope 0
image byres xcen 0 ycen 0 zcen 0 sele resn tip3 .or. resn sod .or. resn cla end
image byseg xcen 0 ycen 0 zcen 0 sele .not. ( resn tip3 .or. resn sod .or. resn cla ) end""")

# Setup nonbonded
pmegrid = checkfft(box,margin=1)
lingo.charmm_script(f'set pmegrid = {pmegrid}')
lingo.charmm_script(f"stream {os.path.join(nbondDir,'nbond.str')}")

initialState = f"stream {trajname}_charges.txt"
if trajname == 'ff':
    lingo.charmm_script(initialState)
    finalState = f"stream crn_charges.txt" 

elif trajname == 'crn':
    lingo.charmm_script(initialState)
    finalState = f"stream ff_charges.txt" 

mini.run_sd(**{'nstep': 100,
               'tolenr': 1e-6,
               'tolgrd': 1e-3})

# nsteps = 2500000
# nsavc = 5000
# nequil = 500000
nsteps = 50
nsavc = 10
nequil = 10
nframes = int(nsteps/nsavc)
run_md(trajname, nsteps=nsteps, nsavc=nsavc, nequil=nequil)

# Open traj 
lingo.charmm_script(f'''open unit 51 read unform name dcd/{trajname}.dcd
    traj first 51 nunit 1 skip {nsavc}
    ''')

# Get Energies
e1 = []
e2 = []
for _ in range(nframes):
    # Get original trajectory energy
    lingo.charmm_script('''traj read
energy''')
    e1.append(lingo.get_energy_value('ENER'))

    # Swap charges and get energy
    lingo.charmm_script(finalState)
    lingo.charmm_script('energy')
    e2.append(lingo.get_energy_value('ENER'))

    # Swap back to original trajectory charges
    lingo.charmm_script(initialState)

e1 = np.array(e1,dtype=float)
e2 = np.array(e2,dtype=float)

np.savetxt(f"{trajname}_e1.csv",e1)
np.savetxt(f"{trajname}_e2.csv",e2)


