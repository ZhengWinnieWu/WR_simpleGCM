# $Id: intel.mk,v 1.1.2.1 2011/01/25 01:10:50 afy Exp $
# template for the Intel fortran compiler
# typical use with mkmf
# mkmf -t template.ifc -c"-Duse_libMPI -Duse_netCDF" path_names /usr/local/include
############
# commands #
############
FC = ifort #single CPU 
FC = mpif90 #using MPI
CC = mpicc
LD = mpif90
LD = $(FC)
#########
# flags #
#########
DEBUG =
REPRO =
VERBOSE =
OPENMP =

# often, there's a pre-defined library and include variable
# we call them MPI_LIB and MPI_INC here
# otherwise, try the mpif90 command, but that doesn't always work
LIB_MPI := #-L$(MPI_LIB) $(shell mpif90 -showme:link)
INC_MPI := #-I$(MPI_INC) $(shell mpif90 -showme:compiler)

INC_NETCDF := -I$(NETCDF_INCCDIR) -I$(NETCDF_INCFDIR) #$(shell nc-config --fflags)
LIB_NETCDF :=  -L$(NETCDF_LIBCDIR) -Wl,-rpath=$(NETCDF_LIBCDIR) -L$(NETCDF_LIBFDIR) -Wl,-rpath=$(NETCDF_LIBFDIR) -lnetcdff -lnetcdf #$(shell nc-config --flibs)

FPPFLAGS := -fpp -Wp,-w
CPPFLAGS := $(INC_NETCDF) $(INC_MPI)

FFLAGS := -axCORE-AVX2,AVX,SSE4.2 -fno-alias -stack_temps -safe_cray_ptr -ftz -assume byterecl -i4 -r8 -nowarn -g -convert big_endian -fPIC -mcmodel=medium 
FFLAGS += $(INC_NETCDF) $(INC_MPI)
FFLAGS_OPT = -O2
FFLAGS_REPRO = -fltconsistency
FFLAGS_DEBUG = -O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -debug variable_locations -fpe0 -traceback -ftrapuv -convert big_endian -fPIC -mcmodel=medium 
#FFLAGS_DEBUG = -O0 -traceback -ftrapuv
FFLAGS_OPENMP = -openmp
FFLAGS_VERBOSE = -v -V -what


CFLAGS := -axCORE-AVX2,AVX,SSE4.2 -D__IFC -mcmodel=medium 
CFLAGS_OPT = -O2
CFLAGS_OPENMP = -openmp
CFLAGS_DEBUG = -O0 -g -ftrapuv -traceback

LDFLAGS := -axCORE-AVX2,AVX,SSE4.2 -mcmodel=medium -shared-intel
LDFLAGS_OPENMP := -openmp
LDFLAGS_VERBOSE := -Wl,-V,--verbose,-cref,-M

LIBS := $(shell nc-config --flibs)
LDFLAGS += $(LIBS)
#LDFLAGS += -lz

ifneq ($(REPRO),)
CFLAGS += $(CFLAGS_REPRO)
FFLAGS += $(FFLAGS_REPRO)
endif
ifneq ($(DEBUG),)
CFLAGS += $(CFLAGS_DEBUG)
FFLAGS += $(FFLAGS_DEBUG)
else
CFLAGS += $(CFLAGS_OPT)
FFLAGS += $(FFLAGS_OPT)
endif

ifneq ($(OPENMP),)
CFLAGS += $(CFLAGS_OPENMP)
FFLAGS += $(FFLAGS_OPENMP)
LDFLAGS += $(LDFLAGS_OPENMP)
endif

ifneq ($(VERBOSE),)
CFLAGS += $(CFLAGS_VERBOSE)
FFLAGS += $(FFLAGS_VERBOSE)
LDFLAGS += $(LDFLAGS_VERBOSE)
endif

ifeq ($(NETCDF),3)
  # add the use_LARGEFILE cppdef
  ifneq ($(findstring -Duse_netCDF,$(CPPDEFS)),)
    CPPDEFS += -Duse_LARGEFILE
  endif
endif

FFLAGS += $(CPPDEFS) $(FPPFLAGS)

LIBS := $(LIB_MPI) $(LIB_NETCDF)
LDFLAGS += $(LIBS)

#---------------------------------------------------------------------------
# you should never need to change any lines below.

# see the MIPSPro f90 manual for more details on some of the file extensions
# discussed here.
# this makefile template recognizes fortran sourcefiles with extensions
# .f, .f90, .F, .f90. Given a sourcefile <file>.<ext>, where <ext> is one of
# the above, this provides a number of default actions:

# make <file>.opt	create an optimization report
# make <file>.o		create an object file
# make <file>.s		create an assembly listing
# make <file>.x		create an executable file, assuming standalone
#			source
# make <file>.i		create a preprocessed file (for .F)
# make <file>.i90	create a preprocessed file (for .f90)

# The macro TMPFILES is provided to slate files like the above for removal.

RM = rm -f
SHELL = /bin/csh -f
TMPFILES = .*.m *.B *.L *.i *.i90 *.l *.s *.mod *.opt

.SUFFIXES: .F .f90 .H .L .T .f .f90 .h .i .i90 .l .o .s .opt .x

.f.L:
	$(FC) $(FFLAGS) -c -listing $*.f
.f.opt:
	$(FC) $(FFLAGS) -c -opt_report_level max -opt_report_phase all -opt_report_file $*.opt $*.f
.f.l:
	$(FC) $(FFLAGS) -c $(LIST) $*.f
.f.T:
	$(FC) $(FFLAGS) -c -cif $*.f
.f.o:
	$(FC) $(FFLAGS) -c $*.f
.f.s:
	$(FC) $(FFLAGS) -S $*.f
.f.x:
	$(FC) $(FFLAGS) -o $*.x $*.f *.o $(LDFLAGS)
.f90.L:
	$(FC) $(FFLAGS) -c -listing $*.f90
.f90.opt:
	$(FC) $(FFLAGS) -c -opt_report_level max -opt_report_phase all -opt_report_file $*.opt $*.f90
.f90.l:
	$(FC) $(FFLAGS) -c $(LIST) $*.f90
.f90.T:
	$(FC) $(FFLAGS) -c -cif $*.f90
.f90.o:
	$(FC) $(FFLAGS) -c $*.f90
.f90.s:
	$(FC) $(FFLAGS) -c -S $*.f90
.f90.x:
	$(FC) $(FFLAGS) -o $*.x $*.f90 *.o $(LDFLAGS)
.F.L:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -listing $*.F
.F.opt:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -opt_report_level max -opt_report_phase all -opt_report_file $*.opt $*.F
.F.l:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c $(LIST) $*.F
.F.T:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -cif $*.F
.F.f:
	$(FC) $(CPPDEFS) $(FPPFLAGS) -EP $*.F > $*.f
.F.i:
	$(FC) $(CPPDEFS) $(FPPFLAGS) -P $*.F
.F.o:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c $*.F
.F.s:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -S $*.F
.F.x:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -o $*.x $*.F *.o $(LDFLAGS)
.f90.L:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -listing $*.f90
.f90.opt:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -opt_report_level max -opt_report_phase all -opt_report_file $*.opt $*.f90
.f90.l:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c $(LIST) $*.f90
.f90.T:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -cif $*.f90
.f90.f90:
	$(FC) $(CPPDEFS) $(FPPFLAGS) -EP $*.f90 > $*.f90
.f90.i90:
	$(FC) $(CPPDEFS) $(FPPFLAGS) -P $*.f90
.f90.o:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c $*.f90
.f90.s:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -S $*.f90
.f90.x:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -o $*.x $*.f90 *.o $(LDFLAGS)
