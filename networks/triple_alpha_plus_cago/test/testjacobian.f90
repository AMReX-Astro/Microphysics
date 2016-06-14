! a simple code to check the analytic Jacobian via numerical 
! differencing

program testjacobian

  use bl_types
  use bl_constants_module
  use bl_error_module
  use network
  use burn_type_module
  use microphysics_module
  use numerical_jac_module

  implicit none

  real(kind=dp_t) :: dens, temp
  real(kind=dp_t) :: Xin(nspec)

  type (burn_t) :: state

  call microphysics_init()

  dens = 1.0e6_dp_t
  temp = 2.e8_dp_t

  Xin = 0.0d0
  Xin(ihe4) = 0.5d0
  Xin(ic12) = 0.5d0

  state % rho = dens
  state % T = temp
  state % xn(:) = Xin(:)

  call test_numerical_jac(state)

  call microphysics_finalize()

end program testjacobian
