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

  dens = 1.0e6_rt
  temp = 2.e8_rt

  Xin = 0.0e0_rt
  Xin(ihe4) = 0.5e0_rt
  Xin(ic12) = 0.5e0_rt

  state % rho = dens
  state % T = temp
  state % xn(:) = Xin(:)

  call test_numerical_jac(state)

  call microphysics_finalize()

end program testjacobian
