NDEBUG := t
MPI    :=
OMP    := 

COMP := gfortran

MKVERBOSE := t

programs += main

MICROPHYSICS_HOME := ../..

EOS_DIR := helmholtz

#NETWORK_DIR := rprox
#NETWORK_DIR := triple_alpha_plus_cago
#NETWORK_DIR := ignition_chamulak
NETWORK_DIR := aprox13
#NETWORK_DIR := ignition_reaclib/URCA-simple

INTEGRATOR_DIR := BS

TEST_DIR := $(MICROPHYSICS_HOME)/unit_test/burn_cell

PROBIN_PARAMETER_DIRS := $(TEST_DIR)

ALL:


include $(MICROPHYSICS_HOME)/GMicrophysics.mak


ALL: $(pnames)

# include the fParallel Makefile rules
include $(AMREX_HOME)/Tools/F_mk/GMakerules.mak




