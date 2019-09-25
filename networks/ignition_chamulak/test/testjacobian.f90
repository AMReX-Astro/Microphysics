! a simple code to check the analytic Jacobian via numerical
! differencing

program testjacobian

  use amrex_constants_module
  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real
  use network
  use burn_type_module
  use microphysics_module
  use numerical_jac_module

  implicit none

  real(rt) :: dens, temp
  real(rt) :: Xin(nspec)

  type (burn_t) :: state

  call microphysics_init()

  dens = 2.6e9_rt
  temp = 6.e8_rt

  Xin(ic12)  = 0.5_rt
  Xin(io16)  = 0.5_rt
  Xin(iash)  = 0.0_rt

  state % rho = dens
  state % T = temp
  state % xn(:) = Xin(:)

  call test_numerical_jac(state)

  call microphysics_finalize()

end program testjacobian
