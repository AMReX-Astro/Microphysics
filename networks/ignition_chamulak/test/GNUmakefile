NDEBUG := 
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

MAESTRO_TOP_DIR := ../../../..

programs += testburn
programs += eval
programs += testjacobian

EOS_DIR := helmeos
CONDUCTIVITY_DIR := constant
NETWORK_DIR := ignition_chamulak

BOXLIB_ONLY := t
EXTRA_DIR = Microphysics/networks/ignition_chamulak/test

testAll: table

EXTRA_LOCATIONS += $(MAESTRO_TOP_DIR)/Source  # needed for constants.f90

include $(MAESTRO_TOP_DIR)/GMaestro.mak

testAll: $(pnames)

%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o $@ $< $(objects) $(libraries)
else
	@echo "Linking $@ ... "
	@$(LINK.f90) -o $@ $< $(objects) $(libraries)
endif



