NDEBUG := t
MPI    :=
OMP    := 

COMP := gfortran

MKVERBOSE := t

programs += main

MICROPHYSICS_HOME := ../..

EOS_DIR := helmholtz

NETWORK_DIR := aprox13

INTEGRATOR_DIR := VODE

TEST_DIR := $(MICROPHYSICS_HOME)/util/cj_detonation

PROBIN_PARAMETER_DIRS := $(TEST_DIR)

ALL:


include $(MICROPHYSICS_HOME)/GMicrophysics.mak


ALL: $(pnames)

# include the fParallel Makefile rules
include $(AMREX_HOME)/Tools/F_mk/GMakerules.mak




