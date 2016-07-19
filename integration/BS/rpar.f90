! rpar is a set of real quantities that we use to convert auxiliary data back
! and forth between the Microphysics data format and the specific data format
! used by the integrators.

module rpar_indices

  use actual_network, only: nspec, nspec_evolve, nrates
  use burn_type_module, only: neqs, num_rate_groups

  implicit none

  integer, parameter :: irp_t_sound = 1
  integer, parameter :: irp_y_init = irp_t_sound + 1

  integer, parameter :: n_rpar_comps = irp_y_init + nspec

end module rpar_indices
