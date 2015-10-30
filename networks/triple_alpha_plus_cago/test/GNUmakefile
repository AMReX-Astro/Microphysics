NDEBUG := 
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

# location of root directory
MAESTRO_TOP_DIR := $(MAESTRO_HOME)


# programs to be compiled
#programs += testburn
programs += eval
programs += testjacobian

EOS_DIR := helmeos
CONDUCTIVITY_DIR := constant
NETWORK_DIR := triple_alpha_plus_cago

BOXLIB_ONLY := t

EXTRA_DIR = Microphysics/networks/triple_alpha_plus_cago/test

testAll:

include $(MAESTRO_TOP_DIR)/GMaestro.mak

testAll: $(pnames) table


%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o $@ $< $(objects) $(libraries)
else
	@echo "Linking $@ ... "
	@$(LINK.f90) -o $@ $< $(objects) $(libraries)
endif

