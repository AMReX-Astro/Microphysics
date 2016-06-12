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

TEST_DIR = $(MICROPHYSICS_DIR)/networks/aprox13/test

ALL: 

include $(MICROPHYSICS_DIR)/GMicrophysics.mak

ALL: $(pnames)


%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o $@ $< $(objects) $(libraries)
else
	@echo "Linking $@ ... "
	@$(LINK.f90) -o $@ $< $(objects) $(libraries)
endif

