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

  character (len=32) :: probin_file
  integer :: probin_pass(32)
  integer :: i, j

  probin_file = "probin"
  do i = 1, len(trim(probin_file))
     probin_pass(i) = ichar(probin_file(i:i))
  enddo

  call runtime_init(probin_pass(1:len(trim(probin_file))), len(trim(probin_file)))

  call network_init()
  call actual_rhs_init()
  call burner_init()
  call eos_init()

  state_ana % rho   = 2.0e7_rt
  state_ana % T     = 8.0e9_rt
  
  state_ana % xn(:) = ONE / nspec

  call test_numerical_jac(state_ana)
  
end subroutine test_jacobian
