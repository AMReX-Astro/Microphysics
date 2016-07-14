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

  dens = 2.6e9_dp_t
  temp = 6.e8_dp_t

  Xin(ic12_)  = 0.5_dp_t
  Xin(io16_)  = 0.5_dp_t
  Xin(iash_)  = 0.0_dp_t

  state % rho = dens
  state % T = temp
  state % xn(:) = Xin(:)

  call test_numerical_jac(state)

  call microphysics_finalize()

end program testjacobian
