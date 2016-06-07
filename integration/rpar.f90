! rpar is a set of real quantities that we use to convert auxiliary data back
! and forth between the Microphysics data format and the specific data format
! used by the integrators.

module rpar_indices

  use network
  use burn_type_module

  implicit none

  integer, parameter :: irp_dens = 1
  integer, parameter :: irp_cv = irp_dens + 1
  integer, parameter :: irp_cp = irp_cv + 1
  integer, parameter :: irp_dedY = irp_cp + 1
  integer, parameter :: irp_dhdY = irp_dedY + nspec
  integer, parameter :: irp_nspec = irp_dhdY + nspec
  integer, parameter :: irp_abar = irp_nspec + nspec - nspec_evolve
  integer, parameter :: irp_zbar = irp_abar + 1
  integer, parameter :: irp_eta = irp_zbar + 1
  integer, parameter :: irp_ye = irp_eta + 1
  integer, parameter :: irp_self_heat = irp_ye + 1
  integer, parameter :: irp_have_rates = irp_self_heat + 1
  integer, parameter :: irp_rates = irp_have_rates + 1
  integer, parameter :: irp_Told = irp_rates + num_rate_groups * nrates
  integer, parameter :: irp_dcvdt = irp_Told + 1
  integer, parameter :: irp_dcpdt = irp_dcvdt + 1

  integer, parameter :: n_rpar_comps = irp_dcpdt

end module rpar_indices
