NDEBUG := 
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

MICROPHYSICS_HOME := ../../..

# programs to be compiled
programs += testburn

EOS_DIR := helmholtz
NETWORK_DIR := rprox
INTEGRATOR_DIR := VODE

TEST_DIR = $(MICROPHYSICS_HOME)/networks/rprox/test

ALL:

include $(MICROPHYSICS_HOME)/GMicrophysics.mak

ALL: $(pnames)

# include the fParallel Makefile rules
include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak

