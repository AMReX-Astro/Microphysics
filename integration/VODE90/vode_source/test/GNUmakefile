NDEBUG := t
MPI    :=
OMP    := t

COMP := gfortran

MKVERBOSE := t

programs += tvode

MICROPHYSICS_HOME := ../../../../

EOS_DIR := gamma_law_general

NETWORK_DIR := ignition_reaclib/C-test

INTEGRATOR_DIR := VODE90

TEST_DIR := $(MICROPHYSICS_HOME)/integration/VODE90/vode_source

PROBIN_PARAMETER_DIRS := $(TEST_DIR)

ALL:


include $(MICROPHYSICS_HOME)/GMicrophysics.mak
include ../GPackage.mak
include GPackage.mak

ALL: $(pnames)

# include the fParallel Makefile rules
include $(AMREX_HOME)/Tools/F_mk/GMakerules.mak





