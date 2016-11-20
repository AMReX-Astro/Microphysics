! rpar is a set of real quantities that we use to convert auxiliary data back
! and forth between the Microphysics data format and the specific data format
! used by the integrators.

module rpar_indices

#ifndef SDC
  use actual_network, only: nspec, nspec_evolve
  use burn_type_module, only: neqs
#else
  use sdc_type_module, only: SVAR_EVOLVE
#endif

  implicit none

#ifndef SDC
  integer, parameter :: n_not_evolved = nspec - nspec_evolve

  integer, parameter :: irp_nspec = 1
  integer, parameter :: irp_y_init = irp_nspec + n_not_evolved
  integer, parameter :: irp_t_sound = irp_y_init + neqs
  integer, parameter :: irp_t0 = irp_t_sound + 1

  integer, parameter :: n_rpar_comps = irp_t0 + 1
#else
  integer, parameter :: irp_SRHO = 1
  integer, parameter :: irp_SMX  = 2
  integer, parameter :: irp_SMY  = 3
  integer, parameter :: irp_SMZ  = 4
  integer, parameter :: irp_t0 = irp_SMZ + 1

  integer, parameter :: n_rpar_comps = irp_t0
#endif

end module rpar_indices
