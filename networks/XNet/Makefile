# Makefile for source/physics/sourceTerms/Burn/nuclearBurn/XNet

Burn += bn_initNetwork.o bn_xnetData.o bn_xnetInit.o bn_xnetFinalize.o \
	xnet_interface.o xnet_interface_mpi.o \
	xnet_abundances.o xnet_common.o xnet_conditions.o xnet_constants.o xnet_control.o \
	xnet_data.o xnet_data_distribute_mpi.o xnet_eos.o xnet_fd.o xnet_ffn.o xnet_flux.o \
	xnet_full_net.o xnet_jacobian_bcast_dense.o xnet_jacobian_dense.o xnet_match.o \
	xnet_net_preprocess.o xnet_output.o xnet_screening.o xnet_solve_be.o xnet_terminate.o \
	xnet_types.o xnet_util.o

bn_xnetData.o : xnet_control.o xnet_data.o xnet_eos.o

xnet_abundances.o : xnet_conditions.o xnet_control.o xnet_data.o xnet_eos.o xnet_interface.o xnet_types.o
xnet_conditions.o : xnet_control.o xnet_interface.o xnet_types.o
xnet_constants.o : xnet_types.o
xnet_control.o : xnet_interface.o xnet_types.o
xnet_data.o : xnet_control.o xnet_types.o
xnet_eos.o : Eos_interface.o Simulation_interface.o xnet_constants.o xnet_control.o xnet_data.o xnet_fd.o xnet_types.o
xnet_ffn.o : xnet_types.o
xnet_flux.o : xnet_types.o
xnet_interface.o : xnet_types.o
xnet_jacobian_dense.o : xnet_types.o
xnet_match.o : xnet_types.o
xnet_net_preprocess.o : xnet_types.o
