#!/bin/csh
set echo
module load intel mvapich2 netcdf-c/4.3.3.1 netcdf-f/4.4.0.i
cp intel.mk Makefile.fms_spectral_solo Makefile obj
cd obj
make
cd ../../src/mppnccombine
./mppnccombine_compile
