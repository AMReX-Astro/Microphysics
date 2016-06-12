NDEBUG := t
MPI    :=
OMP    := t

COMP := gfortran

MKVERBOSE := t

programs += varden


EOS_DIR := helmholtz

#NETWORK_DIR := rprox
#NETWORK_DIR := triple_alpha
#NETWORK_DIR := ignition_chamulak
#NETWORK_DIR := aprox13
NETWORK_DIR := ignition_simple

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


