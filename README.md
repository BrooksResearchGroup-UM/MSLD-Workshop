# MSLD-Workshop
## Materials for five day workshop on using MSLD (including ALF and MSLD_py_prep) for drug discovery and protein engineering
### List of to dos before making public
- describe requirements for participation
- other stuff?

_Note, I have added a file with a template for building a conda environment from which CHARMM/pyCHARMM can be run and a number of other Python packages are included. This can be added to, but should certainly be tested._

_Also, I just added two examples from the pyCHARMM-Workshop github, the basic building of an alanine dipeptide then solvating it, and the set-up of a protein system (w/o paying attention to ionizable groups) and running dynamics with OpenMM or BLaDe. All visualization is now done with nglview._

![Workshop flyer](https://github.com/BrooksResearchGroup-UM/MSLD-Workshop/blob/main/flyer.jpg)

### Workshop schedule
| Date         | Program                                                                                                     |
|:------------:|:------------------------------------------------------------------------------------------------------------|
| Jul-31, 2023 | Installing CHARMM, pyCHARMM, and pyALF                                                                      |
| Aug-1, 2023  | Overview of the theories behind free energy methods<br/>Overview of theory behind MSLD and ALF<br/> Basic steps of writing CHARMM/pyCHARMM scripts |
| Aug-2, 2023  | Running ALF<br/> Preparing inputs for ALF using `msld_py_prep` tool.                                        |
| Aug-3, 2023  | Using ALF/MSLD to perform ligand perturbation and compute relative free energies                            |
| Aug-4, 2023  | Using ALF/MSLD to perform protein residue mutation<br/>Using ALF/MSLD to perform constant pH MD simulations |
