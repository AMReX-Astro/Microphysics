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
NETWORK_DIR := ignition_simple

TEST_DIR = $(MICROPHYSICS_DIR)/networks/ignition_simple/test

ALL:

include $(MICROPHYSICS_DIR)/GMicrophysics.mak

ALL: $(pnames)

