## FAQ

Find solutions to some of the frequently encountered issues in installing and running the scripts presented in the workshop.

#### The list of questions is constantly evolving.

#### 1. When I installed cuda 12.0.0 using `mamba install -y -c 'nvidia/label/cuda-12.0.0' cuda`, and gcc 13.1.0 using `mamba install -y -c conda-forge gcc gxx gfortran`, I got an error while compiling pyCHARMM: 
`error -- unsupported GNU version! gcc versions later than 12 are not supported! The nvcc flag '-allow-unsupported-compiler' can be used to override this version check; however, using an unsupported host compiler may cause compilation failure or incorrect run time execution. Use at your own risk.`
- Use a lower version of gcc, like `mamba install -y -c conda-forge gcc==12.1 gxx==12.1 gfortran==12.1`

#### 2. I am using a HPC cluster to install CHARMM/pyCHARMM. While installing CHARMM, I get the following error while configuring using: `../configure --with-blade --with-fftdock -u  -D nvcc_ptx_target=52 -p <charmm_install_path>`

```
CMake Error at CMakeLists.txt:480 (message):
  MPI was not found and is required for DOMDEC_GPU.
```
- It is possible that the right version of gcc was not installed. You can check the version of gcc by typing `gcc --version`. Ask your HPCC admin to locate the right version of GCC. In satyr/gollum, one can explicitly load module gcc/12.1.0 using `module load gcc/12.1.0` with cuda/12.0 or try a bunch of other gcc versions that are available (and compatible with the cuda toolkit).

#### 3. When I build the CHARMM/pyCHARMM compatible environment with YAML file (Steps 1b), it takes extremely long. Is there an alternative?

- Take the steps listed under 1a instead of 1b to expedite installation. Installation of so many packages using conda is known to take a long time. That is why installing packages using mamba is recommended. 
