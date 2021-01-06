! rpar is a set of real quantities that we use to convert auxiliary data back
! and forth between the Microphysics data format and the specific data format
! used by the integrators.

module vode_rpar_indices

#ifdef TRUE_SDC
  use actual_network, only: nspec
#else
  use actual_network, only: nspec, naux
  use burn_type_module, only: neqs
#endif

  use amrex_fort_module, only : rt => amrex_real
  implicit none

#ifdef TRUE_SDC
  ! f_source is function we are zeroing.  There are nspec + 2
  ! components (density and energy), since those are the unknowns for
  ! the nonlinear system
  integer, parameter :: irp_f_source = 1

  ! dt is the timestep (1 component)
  integer, parameter :: irp_dt = irp_f_source + nspec + 2

  ! mom is the momentum (3 components)
  integer, parameter :: irp_mom = irp_dt + 1

  ! evar is the other energy variable (rho e if we are solving for rho E, and vice versa)
  integer, parameter :: irp_evar = irp_mom + 3

  ! the temperature -- used as a guess in the EOS
  integer, parameter :: irp_temp = irp_evar + 1

  integer, parameter :: n_rpar_comps = irp_temp

#else
  integer, parameter :: irp_dens = 1
  integer, parameter :: irp_cv = irp_dens + 1
  integer, parameter :: irp_cp = irp_cv + 1
  integer, parameter :: irp_abar = irp_cp + 1
  integer, parameter :: irp_zbar = irp_abar + 1
  integer, parameter :: irp_eta = irp_zbar + 1
  integer, parameter :: irp_ye = irp_eta + 1
  integer, parameter :: irp_self_heat = irp_ye + 1
  integer, parameter :: irp_Told = irp_self_heat + 1
  integer, parameter :: irp_dcvdt = irp_Told + 1
  integer, parameter :: irp_dcpdt = irp_dcvdt + 1
  integer, parameter :: irp_t0 = irp_dcpdt + 1

  integer, parameter :: n_rpar_comps = irp_t0
#endif


end module vode_rpar_indices
