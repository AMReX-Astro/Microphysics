f90sources += actual_network.f90
f90sourcse += actual_network_data.f90
f90sources += actual_burner.f90
f90sources += actual_burner_data.f90

# actual_network.f90 is created at build time for this network
network.f90:   $(GENERAL_NET_INPUTS) $(MICROPHYSICS_DIR)/Microphysics/networks/general_null/network.template
	@echo " "
	@echo "---------------------------------------------------------------------------"
	@echo "${bold}WRITING actual_network.f90${normal}"
	$(MICROPHYSICS_DIR)/Microphysics/networks/general_null/write_network.py \
            -t $(MAESTRO_TOP_DIR)/Microphysics/networks/general_null/network.template \
            -s $(GENERAL_NET_INPUTS) \
            -o actual_network.f90
	@echo "---------------------------------------------------------------------------"
	@echo " "


# remove actual_network.f90 for 'make clean' and therefore 'make realclean'
clean::
	$(RM) actual_network.f90
