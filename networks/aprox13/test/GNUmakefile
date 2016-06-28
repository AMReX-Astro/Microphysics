NDEBUG := t
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

# programs to be compiled
programs += testburn

EOS_DIR := helmholtz
NETWORK_DIR := aprox13

INTEGRATOR_DIR := BS

TEST_DIR = $(MICROPHYSICS_HOME)/networks/aprox13/test

ALL: 

include $(MICROPHYSICS_HOME)/GMicrophysics.mak

ALL: $(pnames)

