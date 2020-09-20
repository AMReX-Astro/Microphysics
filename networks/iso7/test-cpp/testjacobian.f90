! a simple code to check the analytic Jacobian via numerical 
! differencing

subroutine test_jacobian() bind(C)

  use network
  use eos_module
  use burner_module
  use actual_burner_module
  use actual_rhs_module
  use burn_type_module
  use numerical_jac_module
  
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  type (burn_t) :: state_ana, state_num

  type (eos_t) :: eos_state

  integer :: i, j

  state_ana % rho   = 2.0e7_rt
  state_ana % T     = 8.0e9_rt
  
  state_ana % xn(:) = ONE / nspec

  call test_numerical_jac(state_ana)
  
end subroutine test_jacobian
