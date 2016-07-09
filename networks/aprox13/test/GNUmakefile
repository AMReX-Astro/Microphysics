NDEBUG := t
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

MICROPHYSICS_HOME := ../../..

# programs to be compiled
programs += testburn

EOS_DIR := helmholtz
NETWORK_DIR := aprox13

INTEGRATOR_DIR := BS

TEST_DIR = $(MICROPHYSICS_HOME)/networks/aprox13/test

ALL: 

include $(MICROPHYSICS_HOME)/GMicrophysics.mak

ALL: $(pnames)

# include the fParallel Makefile rules
include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak

