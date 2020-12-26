! a simple code to check the analytic Jacobian via numerical
! differencing

subroutine test_jacobian() bind(C)

  use network
  use eos_module
  use burner_module
  use actual_burner_module
  use actual_rhs_module
  use burn_type_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  type (burn_t) :: state, statep, statem

  type (eos_t) :: eos_state

  integer :: i, j

  real(rt)        , parameter :: delta = 1.e-6_rt
  real(rt)        , parameter :: SMALL = 1.e-12_rt
  real(rt)         :: num_jac

  character(len=16) :: namei,namej
  real(rt) :: ydot(neqs), ydotm(neqs), ydotp(neqs)
  real(rt) :: jac(njrows, njcols)

  ! Set up state

  state % rho   = 2.0e7_rt
  state % T     = 8.0e9_rt

  state % xn(:) = ONE / nspec

  call burn_to_eos(state, eos_state)
  call normalize_abundances(eos_state)
  call eos(eos_input_rt, eos_state)
  call eos_to_burn(eos_state, state)

  state % self_heat = .true.

  ! Evaluate the analytical Jacobian. Note that we
  ! need to call f_rhs first because that will fill
  ! the state with the rates that the Jacobian needs.

  call actual_rhs(state, ydot)
  call actual_jac(state, jac)

888 format(a,"-derivatives that don't match:")
999 format(5x, "df(",a,")/dy(",a,")", g18.10, g18.10, g18.10)

  ! Now evaluate a numerical estimate of the Jacobian
  ! using the RHS.

  do j = 1, neqs

     statep = state
     ydotp = ydot
     statem = state
     ydotm = ydot

     if (j <= nspec) then
        statep % xn(j) = (ONE + delta) * state % xn(j)
     else if (j == net_itemp) then
        statep % T     = (ONE + delta) * state % T
     else if (j == net_ienuc) then
        statep % e     = (ONE + delta) * state % e
     endif

     call actual_rhs(statep, ydotp)

     if (j <= nspec) then
        statem % xn(j) = (ONE - delta) * state % xn(j)
     else if (j == net_itemp) then
        statem % T     = (ONE - delta) * state % T
     else if (j == net_ienuc) then
        statem % e     = (ONE - delta) * state % e
     endif

     call actual_rhs(statem, ydotm)

     if (j <= nspec) then
        namej = short_spec_names(j)
     else if (j==net_ienuc) then
        namej = "e"
     else if (j==net_itemp) then
        namej = "T"
     endif

     write(*,888) trim(namej)

     do i = 2, neqs

        if (j <= nspec) then
           num_jac = (ydotp(i) - ydotm(i))/(statep % xn(j) - statem % xn(j) + SMALL)
        else if (j == net_itemp) then
           num_jac = (ydotp(i) - ydotm(i))/(statep % T - statem % T + SMALL)
        else if (j == net_ienuc) then
           num_jac = (ydotp(i) - ydotm(i))/(statep % e - statem % e + SMALL)
        endif

        if (i <= nspec) then
           namei = short_spec_names(i)
        else if (i==net_ienuc) then
           namei = "e"
        else if (i==net_itemp) then
           namei = "T"
        endif

        ! only dump the ones that don't match
        if (num_jac /= jac(i,j)) then
           if (num_jac /= ZERO) then
              write (*,999) trim(namei), &
                   trim(namej), num_jac, jac(i,j), (num_jac - jac(i,j))/num_jac
           else
              write (*,999) trim(namei), &
                   trim(namej), num_jac, jac(i,j)
           endif
        endif

     enddo
  enddo

end subroutine test_jacobian
