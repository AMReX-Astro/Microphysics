NDEBUG := 
MPI    :=
OMP    := 

SDC := t

COMP := gfortran

MKVERBOSE := t

programs += main

MICROPHYSICS_HOME := ../..

EOS_DIR := helmholtz

NETWORK_DIR := aprox13

INTEGRATOR_DIR := VODE

TEST_DIR := $(MICROPHYSICS_HOME)/unit_test/test_sdc

PROBIN_PARAMETER_DIRS := $(TEST_DIR)

ALL:


include $(MICROPHYSICS_HOME)/GMicrophysics.mak


ALL: $(pnames)

# include the fParallel Makefile rules
include $(FBOXLIB_HOME)/Tools/F_mk/GMakerules.mak





