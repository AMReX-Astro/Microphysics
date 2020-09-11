! rpar is a set of real quantities that we use to convert auxiliary data back
! and forth between the Microphysics data format and the specific data format
! used by the integrators.

module bs_rpar_indices

  use actual_network, only: nspec
  use burn_type_module, only: neqs

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, parameter :: irp_t0 = 1

  integer, parameter :: n_rpar_comps = irp_t0

end module bs_rpar_indices
