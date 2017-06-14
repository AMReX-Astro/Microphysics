NDEBUG := t
MPI    :=
OMP    := t

COMP := gfortran

MKVERBOSE := t

programs += main

MICROPHYSICS_HOME := ../..

EOS_DIR := helmholtz

NETWORK_DIR := aprox19

TEST_DIR := $(MICROPHYSICS_HOME)/unit_test/test_eos

PROBIN_PARAMETER_DIRS := $(TEST_DIR)

ALL:

include $(MICROPHYSICS_HOME)/GMicrophysics.mak

ALL: $(pnames)

# include the fParallel Makefile rules
include $(AMREX_HOME)/Tools/F_mk/GMakerules.mak




