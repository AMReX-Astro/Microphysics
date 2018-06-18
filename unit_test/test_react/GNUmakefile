NDEBUG := t
MPI    :=
OMP    :=
CUDA   := 

COMP := PGI

MKVERBOSE := t

programs += main

MICROPHYSICS_HOME := ../..

EOS_DIR := helmholtz

#NETWORK_DIR := rprox
#NETWORK_DIR := triple_alpha_plus_cago
#NETWORK_DIR := ignition_chamulak
NETWORK_DIR := aprox13
#NETWORK_DIR := ignition_simple

INTEGRATOR_DIR := VODE90

TEST_DIR := $(MICROPHYSICS_HOME)/unit_test/test_react

PROBIN_PARAMETER_DIRS := $(TEST_DIR)

ALL:


include $(MICROPHYSICS_HOME)/GMicrophysics.mak


ALL: $(pnames)

# include the fParallel Makefile rules
include $(FBOXLIB_HOME)/Tools/F_mk/GMakerules.mak

# # CUDA-specific options ala StarLord
# # Limit the maximum number of registers available.
# # Increase the maximum number of threads to compensate.
# ifeq ($(COMP),PGI)
#   ifeq ($(CUDA),t)
#     FFLAGS += -Mcuda=maxregcount:64
#     F90FLAGS += -Mcuda=maxregcount:64
#   endif
# endif
