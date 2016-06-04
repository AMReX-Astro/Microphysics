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
  
  implicit none

  type (burn_t) :: state_ana, state_num

  type (eos_t) :: eos_state

  character (len=32) :: probin_file
  integer :: probin_pass(32)
  integer :: i, j

  double precision, parameter :: delta = 1.d-6
  double precision, parameter :: SMALL = 1.d-12
  double precision :: num_jac

  character(len=16) :: namei,namej  
  
  probin_file = "probin"
  do i = 1, len(trim(probin_file))
     probin_pass(i) = ichar(probin_file(i:i))
  enddo

  call runtime_init(probin_pass(1:len(trim(probin_file))), len(trim(probin_file)))

  call network_init()
  call actual_rhs_init()
  call burner_init()
  call eos_init()

  ! Set up state

  state_ana % rho   = 2.0d7
  state_ana % T     = 8.0d9

  state_ana % xn(:) = ONE / nspec

  call burn_to_eos(state_ana, eos_state)
  call normalize_abundances(eos_state)
  call eos(eos_input_rt, eos_state)
  call eos_to_burn(eos_state, state_ana)

  state_ana % self_heat = .true.

  state_num = state_ana

  ! Evaluate the analytical Jacobian. Note that we
  ! need to call f_rhs first because that will fill
  ! the state with the rates that the Jacobian needs.

  call actual_rhs(state_ana)
  call actual_jac(state_ana)

  call actual_rhs(state_num)
  call numerical_jac(state_num)
  
888 format(a,"-derivatives that don't match:")
999 format(5x, "df(",a,")/dy(",a,")", g18.10, g18.10, g18.10)

  ! Now evaluate a numerical estimate of the Jacobian
  ! using the RHS.
  
  do j = 1, neqs
     
     if (j <= nspec_evolve) then
        namej = short_spec_names(j)
     else if (j==net_ienuc) then
        namej = "e"
     else if (j==net_itemp) then
        namej = "T"
     endif

     write(*,888) trim(namej)

     do i = 1, neqs

        if (i <= nspec_evolve) then
           namei = short_spec_names(i)
        else if (i==net_ienuc) then
           namei = "e"
        else if (i==net_itemp) then
           namei = "T"
        endif

        ! only dump the ones that don't match
        if (state_num % jac(i,j) /= state_ana % jac(i,j)) then
           if (state_num % jac(i,j) /= ZERO) then
              write (*,999) trim(namei), &
                   trim(namej), state_num % jac(i,j), state_ana % jac(i,j), &
                   (state_num % jac(i,j)-state_ana % jac(i,j))/state_num % jac(i,j)
           else
              write (*,999) trim(namei), &
                   trim(namej), state_num % jac(i,j), state_ana % jac(i,j)
           endif
        endif

     enddo
  enddo  
  
end subroutine test_jacobian
