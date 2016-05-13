NDEBUG := t
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

# location of root directories
MAESTRO_TOP_DIR := $(MAESTRO_HOME)

# programs to be compiled
programs += testburn

EOS_DIR := helmeos
CONDUCTIVITY_DIR := constant
NETWORK_DIR := xrb_simple

BOXLIB_ONLY := t

EXTRA_DIR := Microphysics/networks/xrb_simple/test

testAll:

include $(MAESTRO_TOP_DIR)/GMaestro.mak

#F90FLAGS += -lm
CXXFLAGS += -lm

testAll: $(pnames) table


%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o $@ $< $(objects) $(libraries)
else
	@echo "Linking $@ ... "
	@$(LINK.f90) -o $@ $< $(objects) $(libraries)
endif

