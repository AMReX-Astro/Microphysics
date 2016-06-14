NDEBUG := t
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

MICROPHYSICS_DIR := ../../..

# programs to be compiled
programs += testburn

EOS_DIR := helmholtz
NETWORK_DIR := rprox


TEST_DIR = $(MICROPHYSICS_DIR)/networks/rprox/test

ALL:

include $(MICROPHYSICS_DIR)/GMicrophysics.mak

ALL: $(pnames)
