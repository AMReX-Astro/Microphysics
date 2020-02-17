! rpar is a set of real quantities that we use to convert auxiliary data back
! and forth between the Microphysics data format and the specific data format
! used by the integrators.

module vbdf_rpar_indices

  use actual_network, only: nspec, nspec_evolve
  use burn_type_module, only: neqs

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, parameter :: n_not_evolved = nspec - nspec_evolve

  integer, parameter :: irp_dens = 1
  integer, parameter :: irp_cx = irp_dens + 1
  integer, parameter :: irp_nspec = irp_cx + 1
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
  integer, parameter :: irp_dcxdt = irp_Told + 1
  integer, parameter :: irp_t0 = irp_dcxdt + 1

  integer, parameter :: n_rpar_comps = irp_t0

end module vbdf_rpar_indices
