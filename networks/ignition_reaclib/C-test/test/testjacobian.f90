! a simple code to check the analytic Jacobian via numerical 
! differencing

program testjacobian

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module
  use amrex_error_module
  use network
  use burn_type_module
  use microphysics_module
  use numerical_jac_module

  implicit none

  real(kind=rt) :: dens, temp
  real(kind=rt) :: Xin(nspec)

  type (burn_t) :: state

  call microphysics_init()

  dens = 3.e9_rt
  temp = 5.e8_rt

  Xin(jhe4)   = 0.4_rt
  Xin(jc12)   = 0.4_rt
  Xin(jo16)   = 0.2_rt

  state % rho = dens
  state % T = temp
  state % xn(:) = Xin(:)

  call test_numerical_jac(state)

  call microphysics_finalize()

end program testjacobian
