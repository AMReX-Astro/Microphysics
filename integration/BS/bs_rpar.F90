! rpar is a set of real quantities that we use to convert auxiliary data back
! and forth between the Microphysics data format and the specific data format
! used by the integrators.

module bs_rpar_indices

#ifndef SIMPLIFIED_SDC
  use actual_network, only: nspec, nspec_evolve
  use burn_type_module, only: neqs
#else
  use sdc_type_module, only: SVAR_EVOLVE
#endif

  use amrex_fort_module, only : rt => amrex_real
  implicit none

#ifndef SIMPLIFIED_SDC
  integer, parameter :: n_not_evolved = nspec - nspec_evolve

  integer, parameter :: irp_nspec = 1
  integer, parameter :: irp_y_init = irp_nspec + n_not_evolved
  integer, parameter :: irp_t_sound = irp_y_init + neqs
  integer, parameter :: irp_t0 = irp_t_sound + 1

  integer, parameter :: n_rpar_comps = irp_t0
#else

  ! note: we require these components to be first, to allow for use to
  ! index bs % u_init and bs % udot_a
  integer, parameter :: irp_SRHO = 1

#if defined(SDC_EVOLVE_ENERGY)

  integer, parameter :: irp_SMX  = 2
  integer, parameter :: irp_SMY  = 3
  integer, parameter :: irp_SMZ  = 4

  integer, parameter :: irp_t0 = irp_SMZ + 1

#elif defined(SDC_EVOLVE_ENTHALPY)

  integer, parameter :: irp_p0 = irp_SRHO + 1
  integer, parameter :: irp_t0 = irp_p0 + 1

#endif

  integer, parameter :: n_rpar_comps = irp_t0

#endif

end module bs_rpar_indices
