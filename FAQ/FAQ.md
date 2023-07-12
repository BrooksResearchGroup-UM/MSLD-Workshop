## FAQ

Find solutions to some of the frequently encountered issues in installing and running the scripts presented in the workshop.

#### The list of questions is constantly evolving.

### 1. When I installed cuda 12.0.0 using "mamba install -y -c 'nvidia/label/cuda-12.0.0' cuda", and gcc 13.1.0 using "mamba install -y -c conda-forge gcc gxx gfortran", I got an error while compiling pyCHARMM: "error -- unsupported GNU version! gcc versions later than 12 are not supported! The nvcc flag '-allow-unsupported-compiler' can be used to override this version check; however, using an unsupported host compiler may cause compilation failure or incorrect run time execution. Use at your own risk."
- ** Use a lower version of gcc, like "mamba install -y -c conda-forge gcc==12.1 gxx==12.1 gfortran==12.1"
