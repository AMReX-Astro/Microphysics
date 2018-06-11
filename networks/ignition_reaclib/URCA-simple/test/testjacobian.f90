! a simple code to check the analytic Jacobian via numerical 
! differencing

program testjacobian

  use bl_types
  use bl_constants_module
  use amrex_error_module
  use network
  use burn_type_module
  use microphysics_module
  use numerical_jac_module

  implicit none

  real(kind=dp_t) :: dens, temp
  real(kind=dp_t) :: Xin(nspec)

  type (burn_t) :: state

  call microphysics_init()

  dens = 3.e9_dp_t
  temp = 5.e8_dp_t

  Xin(jn)     = 0.0000_dp_t
  Xin(jp)     = 0.0000_dp_t
  Xin(jhe4)   = 0.0000_dp_t
  Xin(jc12)   = 0.4999_dp_t
  Xin(jo16)   = 0.4999_dp_t
  Xin(jne20)  = 0.0000_dp_t
  Xin(jne23)  = 0.0001_dp_t
  Xin(jna23)  = 0.0001_dp_t
  Xin(jmg23)  = 0.0000_dp_t

  state % rho = dens
  state % T = temp
  state % xn(:) = Xin(:)

  call test_numerical_jac(state)

  call microphysics_finalize()

end program testjacobian
