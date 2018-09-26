module numerical_jac_module

  use amrex_constants_module, only: ZERO, HALF, ONE
  use amrex_fort_module, only : rt => amrex_real
  use network
  use burn_type_module

  implicit none

contains

  subroutine numerical_jac(state)

    !$acc routine seq

    use actual_rhs_module, only: actual_rhs
    use extern_probin_module, only : centered_diff_jac

    implicit none

    type (burn_t)    :: state

    integer          :: m, n

    type (burn_t)    :: state_del, state_delm

    ! the choice of eps should be ~ sqrt(eps), where eps is machine epsilon. 
    ! this balances truncation vs. roundoff error in the differencing
    real(rt), parameter :: eps = 1.d-8

    !$gpu

    state % jac(:,:) = ZERO

    ! default
    call actual_rhs(state)


    if (centered_diff_jac) then
       call copy_burn_t(state_del, state)
       call copy_burn_t(state_delm, state)

       ! species derivatives
       do n = 1, nspec_evolve
          ! perturb species
          state_del % xn = state % xn
          state_del % xn(n) = state % xn(n) * (ONE + eps)

          call actual_rhs(state_del)

          ! We integrate X, so convert from the Y we got back from the RHS

          state_del % ydot(1:nspec_evolve) = state_del % ydot(1:nspec_evolve) * aion(1:nspec_evolve)

          state_delm % xn = state % xn
          state_delm % xn(n) = state % xn(n) * (ONE - eps)

          call actual_rhs(state_delm)

          state_delm % ydot(1:nspec_evolve) = state_del % ydot(1:nspec_evolve) * aion(1:nspec_evolve)

          do m = 1, neqs
             state % jac(m,n) = HALF*(state_del % ydot(m) - state_delm % ydot(m)) / &
                  (eps * state % xn(n))
          enddo
       enddo

       ! temperature derivative
       state_del % xn = state % xn
       state_del % T  = state % T * (ONE + eps)

       call actual_rhs(state_del)

       state_delm % xn = state % xn
       state_delm % T  = state % T * (ONE - eps)

       call actual_rhs(state_delm)

       do m = 1, neqs
          state % jac(m,net_itemp) = HALF*(state_del % ydot(m) - state_delm % ydot(m)) / &
               (eps * state % T)
       enddo

       ! energy derivatives -- these are all 0! (yay!)
       do m = 1, neqs
          state % jac(m,net_ienuc) = ZERO
       enddo

    else
       call copy_burn_t(state_del, state)

       ! species derivatives
       do n = 1, nspec_evolve
          ! perturb species -- we send in X, but ydot is in terms of dY/dt, not dX/dt
          state_del % xn = state % xn
          state_del % xn(n) = state % xn(n) * (ONE + eps)

          call actual_rhs(state_del)

          do m = 1, neqs
             state % jac(m,n) = (state_del % ydot(m) - state % ydot(m)) / (eps * state % xn(n))
          enddo
       enddo

       ! temperature derivative
       state_del % xn = state % xn
       state_del % T  = state % T * (ONE + eps)

       call actual_rhs(state_del)

       do m = 1, neqs
          state % jac(m,net_itemp) = (state_del % ydot(m) - state % ydot(m)) / (eps * state % T)
       enddo

       ! energy derivatives -- these are all 0! (yay!)
       do m = 1, neqs
          state % jac(m,net_ienuc) = ZERO
       enddo

    endif

  end subroutine numerical_jac

#ifndef AMREX_USE_CUDA
  subroutine test_numerical_jac(state)
    ! compare the analytic Jacobian to the numerically differenced one

    use actual_rhs_module
    use eos_module, only : eos
    use eos_type_module, only : eos_t, eos_input_rt, normalize_abundances

    type (burn_t) :: state
    type (burn_t) :: state_num
    type (eos_t) :: eos_state

    integer :: i, j
    character(len=16) :: namei, namej  
      
    ! Set up state

    call burn_to_eos(state, eos_state)
    call normalize_abundances(eos_state)
    call eos(eos_input_rt, eos_state)
    call eos_to_burn(eos_state, state)

    state % self_heat = .true.

    state_num = state

    ! Evaluate the analytical Jacobian. Note that we need to call
    ! f_rhs first because that will fill the state with the rates that
    ! the Jacobian needs.
    call actual_rhs(state)
    call actual_jac(state)

    ! Now compute the numerical Jacobian.
    call actual_rhs(state_num)
    call numerical_jac(state_num)
  
888 format(a,"-derivatives that don't match:")
999 format(5x, "df(",a,")/dy(",a,")", g18.10, g18.10, g18.10)

    ! how'd we do?
    do j = 1, neqs
     
       if (j <= nspec_evolve) then
          namej = short_spec_names(j)
       else if (j == net_ienuc) then
          namej = "e"
       else if (j == net_itemp) then
          namej = "T"
       endif

       write(*,888) trim(namej)

       do i = 1, neqs

          if (i <= nspec_evolve) then
             namei = short_spec_names(i)
          else if (i == net_ienuc) then
             namei = "e"
          else if (i == net_itemp) then
             namei = "T"
          endif

          ! only dump the ones that don't match
          if (state_num % jac(i,j) /= state % jac(i,j)) then
             if (state_num % jac(i,j) /= ZERO) then
                write (*,999) trim(namei), &
                     trim(namej), state_num % jac(i,j), state % jac(i,j), &
                     (state_num % jac(i,j) - state % jac(i,j))/state_num % jac(i,j)
             else
                write (*,999) trim(namei), &
                     trim(namej), state_num % jac(i,j), state % jac(i,j)
             endif
          endif

       enddo
    enddo

  end subroutine test_numerical_jac
#endif
end module numerical_jac_module
