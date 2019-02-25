NEBUG := t
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
#NETWORK_DIR := subch
#NETWORK_DIR := ignition_reaclib/URCA-simple
NETWORK_DIR := ECSN_weak

INTEGRATOR_DIR := VODE

TEST_DIR := $(MICROPHYSICS_HOME)/unit_test/burn_cell

PROBIN_PARAMETER_DIRS := $(TEST_DIR)

ALL:


include $(MICROPHYSICS_HOME)/GMicrophysics.mak


ALL: $(pnames)

# include the fParallel Makefile rules
include $(FBOXLIB_HOME)/Tools/F_mk/GMakerules.mak

F90FLAGS += -ffree-line-length-none
FFLAGS += -ffree-line-length-none
