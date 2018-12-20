NDEBUG := f
MPI    :=
OMP    := 

COMP := gfortran

MKVERBOSE := t

programs += evaluate_rhs

MICROPHYSICS_HOME := ../..

EOS_DIR := helmholtz

#NETWORK_DIR := rprox
#NETWORK_DIR := triple_alpha_plus_cago
#NETWORK_DIR := ignition_chamulak
#NETWORK_DIR := aprox13
#NETWORK_DIR := ignition_simple
NETWORK_DIR := ignition_reaclib/URCA-simple

TEST_DIR := $(MICROPHYSICS_HOME)/util/eval_rhs

PROBIN_PARAMETER_DIRS := $(TEST_DIR)

ALL:


include $(MICROPHYSICS_HOME)/GMicrophysics.mak


ALL: $(pnames)

# include the fParallel Makefile rules
include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak





