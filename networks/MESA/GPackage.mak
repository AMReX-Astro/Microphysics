f90sources += network.f90
f90sources += burner.f90
f90sources += net_utils.f90
f90sources += shutdown_mesa_net.f90
sf90sources += setup_mesa_net.f90
sf90sources += Do_One_Zone_Burn.f90
sf90sources += burn_solout.f90

xtr_libraries += -L$(MESA_DIR)/lib -lnet -leos -lscreen -lrates -lreaclib -lweak -lchem -linterp_2d -linterp_1d -lnum -lutils -lalert -lconst -lmtx -lmesaklu -lmesalapack -lmesablas

FINCLUDE_LOCATIONS += $(MESA_DIR)/include

# network.f90 is created at build time for this network
network.f90:   $(GENERAL_NET_INPUTS) $(ASTRODEV_DIR)/networks/MESA/network.template.mesa
	@echo " "
	@echo "${bold}WRITING network.f90${normal}"
	$(ASTRODEV_DIR)/networks/MESA/write_network.py \
            -t $(ASTRODEV_DIR)/networks/MESA/network.template.mesa \
            -s $(GENERAL_NET_INPUTS) \
            -o network.f90
	@echo " "


# remove network.f90 for 'make clean' and therefore 'make realclean'
clean::
	$(RM) network.f90
