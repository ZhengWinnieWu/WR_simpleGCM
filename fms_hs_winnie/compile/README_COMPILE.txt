How to compile this code. In parentheses are examples.

1) Add compilers and libraries. Need compiler (ifort), MPI (mpich2,openmpi), and netCDF
2) Make sure the compiler variables FC and CC in intel.mk correspond to your compilers (mpif90,icc)
3) Add the library path to the environment. (Needed is the command 'nc-config', coming from netCDF). 
4) replace /PATH/TO in Makefile.fms_spectral_solo with the correct path where the src folder is located
5) Compile:
> make
> cd ../src/mppnccombine
> ./mppnccombine_compile
