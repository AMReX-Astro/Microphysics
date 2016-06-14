NDEBUG := t
MPI    :=
OMP    := t

COMP := gfortran

MKVERBOSE := t

programs += main

MICROPHYSICS_DIR := ../..

EOS_DIR := helmholtz

NETWORK_DIR := aprox19

TEST_DIR := $(MICROPHYSICS_DIR)/unit_test/test_eos

PROBIN_PARAMETER_DIRS := $(TEST_DIR)

ALL:

include $(MICROPHYSICS_DIR)/GMicrophysics.mak

ALL: $(pnames)



