NDEBUG := t
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

MICROPHYSICS_DIR := ../../..

# programs to be compiled
programs += testburn

EOS_DIR := helmholtz
NETWORK_DIR := xrb_simple


TEST_DIR = $(MICROPHYSICS_DIR)/networks/xrb_simple/test

ALL:

include $(MICROPHYSICS_DIR)/GMicrophysics.mak

ALL: $(pnames)
