! rpar is a set of real quantities that we use to convert auxiliary data back
! and forth between the Microphysics data format and the specific data format
! used by the integrators.

module sdc_rpar_indices

  use actual_network, only: nspec, nspec_evolve
  use burn_type_module, only: neqs

  implicit none

  integer, parameter :: n_not_evolved = nspec - nspec_evolve

  integer, parameter :: irp_nspec = 1
  integer, parameter :: irp_y_init = irp_nspec + n_not_evolved
  integer, parameter :: irp_t_sound = irp_y_init + neqs
  integer, parameter :: irp_t0 = irp_t_sound + 1

  integer, parameter :: n_rpar_comps = irp_t0

end module sdc_rpar_indices
