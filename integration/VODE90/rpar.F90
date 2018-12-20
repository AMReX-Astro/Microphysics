! rpar is a set of real quantities that we use to convert auxiliary data back
! and forth between the Microphysics data format and the specific data format
! used by the integrators.

! SDC TODO: we need to store the advective sources here.  Also some bits
! are not needed in the SDC implementation (like the cv and cp stuff)

module rpar_indices

#ifndef SDC
  use actual_network, only: nspec, nspec_evolve
  use burn_type_module, only: neqs
#else
  use sdc_type_module, only: SVAR, SVAR_EVOLVE
#endif

  implicit none


#ifndef SDC
  integer, parameter :: n_not_evolved = nspec - nspec_evolve

  integer, parameter :: irp_dens = 1
  integer, parameter :: irp_cv = irp_dens + 1
  integer, parameter :: irp_cp = irp_cv + 1
  integer, parameter :: irp_nspec = irp_cp + 1
  integer, parameter :: irp_abar = irp_nspec + n_not_evolved
  integer, parameter :: irp_zbar = irp_abar + 1
  integer, parameter :: irp_eta = irp_zbar + 1
  integer, parameter :: irp_ye = irp_eta + 1
  integer, parameter :: irp_cs = irp_ye + 1
  integer, parameter :: irp_dx = irp_cs + 1
  integer, parameter :: irp_t_sound = irp_dx + 1
  integer, parameter :: irp_y_init = irp_t_sound + 1
  integer, parameter :: irp_self_heat = irp_y_init + neqs
  integer, parameter :: irp_Told = irp_self_heat + 1
  integer, parameter :: irp_dcvdt = irp_Told + 1
  integer, parameter :: irp_dcpdt = irp_dcvdt + 1
  integer, parameter :: irp_t0 = irp_dcpdt + 1

  integer, parameter :: n_rpar_comps = irp_t0 + 1
#else
  ! Note: we require these components to be first, to allow for offset
  ! indexing with irp_ydot_a and irp_u_init
  integer, parameter :: irp_SRHO = 1
  integer, parameter :: irp_SMX  = 2
  integer, parameter :: irp_SMY  = 3
  integer, parameter :: irp_SMZ  = 4

  ! these are the SVAR advective terms (including those that we do
  ! not explicitly evolve).  Note: these need to be indexed using
  ! the indicies defined in sdc_type_module
  integer, parameter :: irp_ydot_a = 5

  ! these are the SVAR - SVAR_EVOLVE initial values of the unevolved
  ! components of the conserved state.  Note: these need to be indexed
  ! using the components defined above (irp_SRHO, irp_SMX, ...)
  integer, parameter :: irp_u_init = irp_ydot_a + SVAR

  integer, parameter :: irp_self_heat = irp_u_init + SVAR - SVAR_EVOLVE
  integer, parameter :: irp_T_from_eden = irp_self_heat + 1
  integer, parameter :: irp_t0 = irp_self_heat + 1

  integer, parameter :: n_rpar_comps = irp_t0
#endif
  

end module rpar_indices
