NDEBUG := 
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

# location of root directory
MAESTRO_TOP_DIR := $(MAESTRO_HOME)


# programs to be compiled
programs += testburn

EOS_DIR := helmeos
CONDUCTIVITY_DIR := constant

NETWORK_TOP_DIR := $(ASTRODEV_DIR)/networks
NETWORK_DIR := aprox13

BOXLIB_ONLY := t

EXTRA_TOP_DIR := $(ASTRODEV_DIR)
EXTRA_DIR = networks/triple_alpha_plus_cago/test

testAll:

include $(MAESTRO_TOP_DIR)/GMaestro.mak

testAll: $(pnames)


%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o $@ $< $(objects) $(libraries)
else
	@echo "Linking $@ ... "
	@$(LINK.f90) -o $@ $< $(objects) $(libraries)
endif

