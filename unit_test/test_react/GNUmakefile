NDEBUG := t
MPI    :=
OMP    := t

COMP := gfortran

MKVERBOSE := t

programs += main

MICROPHYSICS_HOME := ../..

EOS_DIR := helmholtz

#NETWORK_DIR := rprox
#NETWORK_DIR := triple_alpha
#NETWORK_DIR := ignition_chamulak
#NETWORK_DIR := aprox13
NETWORK_DIR := ignition_simple

INTEGRATOR_DIR := VODE

TEST_DIR := $(MICROPHYSICS_HOME)/unit_test/test_react

PROBIN_PARAMETER_DIRS := $(TEST_DIR)

ALL:


ALL: $(pnames)

include $(MICROPHYSICS_HOME)/GMicrophysics.mak



