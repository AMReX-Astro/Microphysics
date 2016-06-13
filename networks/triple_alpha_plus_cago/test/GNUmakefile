NDEBUG := 
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

# programs to be compiled
programs += testburn
#programs += eval
#programs += testjacobian

EOS_DIR := helmholtz
NETWORK_DIR := triple_alpha_plus_cago

TEST_DIR = $(MICROPHYSICS_DIR)/networks/triple_alpha_plus_cago/test

ALL:

include $(MICROPHYSICS_DIR)/GMicrophysics.mak

ALL: $(pnames)

