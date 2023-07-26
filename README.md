# MSLD-Workshop
## Materials for five day workshop on using MSLD (including ALF and MSLD_py_prep) for drug discovery and protein engineering
### List of to dos before making public
- Most of these scripts will be run using the Slurm Workload Manager.
- We highly encourage you work in a Linux environment, since CUDA is no longer supported using the Apple installation.

- While note explicitly required, please try to have a CHARMM-GUI account, which is free of charge at: charmm-gui.org

  
_Note, I have added a file with a template for building a conda environment from which CHARMM/pyCHARMM can be run and a number of other Python packages are included. This can be added to, but should certainly be tested._

_Also, I just added two examples from the pyCHARMM-Workshop github, the basic building of an alanine dipeptide then solvating it, and the set-up of a protein system (w/o paying attention to ionizable groups) and running dynamics with OpenMM or BLaDe. All visualization is now done with nglview._

### Workshop schedule
| Date         | Program                                                                                                     |
|:------------:|:------------------------------------------------------------------------------------------------------------|
| Mon, July 31, 2023 | Installing CHARMM, pyCHARMM, and pyALF                                                                      |
| Tue, August 1, 2023  | Overview of the theories behind free energy methods<br>Overview of theory behind MSLD and ALF<br> Basic steps of writing CHARMM/pyCHARMM scripts |
| Wed, August 2, 2023  | Running ALF<br> Preparing inputs for ALF using `msld_py_prep` tool.                                        |
| Thu, August 3, 2023  | Using ALF/MSLD to perform ligand perturbation and compute relative free binding energies                            |
| Fri, August 4, 2023  | Using ALF/MSLD to perform protein perturbation to <br>compute relative folding free energies<br>Using ALF/MSLD to <br>perform constant pH MD simulations |

![Workshop flyer](https://github.com/BrooksResearchGroup-UM/MSLD-Workshop/blob/main/flyer.jpg)


