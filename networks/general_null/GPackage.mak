f90sources += general_null_network.f90
f90sources += general_null_burner.f90

# general_null_network.f90 is created at build time for this network
network.f90:   $(GENERAL_NET_INPUTS) $(MICROPHYSICS_DIR)/Microphysics/networks/general_null/network.template
	@echo " "
	@echo "---------------------------------------------------------------------------"
	@echo "${bold}WRITING general_null_network.f90${normal}"
	$(MICROPHYSICS_DIR)/Microphysics/networks/general_null/write_network.py \
            -t $(MAESTRO_TOP_DIR)/Microphysics/networks/general_null/network.template \
            -s $(GENERAL_NET_INPUTS) \
            -o general_null_network.f90
	@echo "---------------------------------------------------------------------------"
	@echo " "


# remove network.f90 for 'make clean' and therefore 'make realclean'
clean::
	$(RM) general_null_network.f90
