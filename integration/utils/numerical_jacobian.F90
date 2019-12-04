module numerical_jac_module

  use amrex_constants_module, only: ZERO, HALF, ONE
  use amrex_fort_module, only : rt => amrex_real
  use network
  use burn_type_module

  implicit none

contains

  subroutine numerical_jac(state, jac)

    !$acc routine seq

    use actual_rhs_module, only: actual_rhs
    use extern_probin_module, only : centered_diff_jac
    use jacobian_sparsity_module, only: set_jac_zero, set_jac_entry

    implicit none

    type (burn_t)    :: state
    real(rt) :: jac(neqs, neqs)
    integer          :: m, n

    real(rt) :: ydotp(neqs), ydotm(neqs)

    type (burn_t)    :: state_delp, state_delm

    ! the choice of eps should be ~ sqrt(eps), where eps is machine epsilon.
    ! this balances truncation vs. roundoff error in the differencing
    real(rt), parameter :: eps = 1.e-8_rt
    real(rt) :: scratch, h

    !$gpu

    call set_jac_zero(jac)


    if (centered_diff_jac) then
       call copy_burn_t(state_delp, state)
       call copy_burn_t(state_delm, state)

       ! species derivatives
       do n = 1, nspec_evolve
          ! perturb species
          state_delp % xn = state % xn
          state_delp % xn(n) = state % xn(n) * (ONE + eps)

          call actual_rhs(state_delp, ydotp)

          ! We integrate X, so convert from the Y we got back from the RHS

          ydotp(1:nspec_evolve) = ydotp(1:nspec_evolve) * aion(1:nspec_evolve)

          state_delm % xn = state % xn
          state_delm % xn(n) = state % xn(n) * (ONE - eps)

          call actual_rhs(state_delm, ydotm)

          ydotm(1:nspec_evolve) = ydotm(1:nspec_evolve) * aion(1:nspec_evolve)

          do m = 1, neqs
             scratch = HALF*(ydotp(m) - ydotm(m)) / (eps * state % xn(n))
             call set_jac_entry(jac, m, n, scratch)
          enddo
       enddo

       ! temperature derivative
       state_delp % xn = state % xn
       state_delp % T  = state % T * (ONE + eps)

       call actual_rhs(state_delp, ydotp)

       ydotp(1:nspec_evolve) = ydotp(1:nspec_evolve) * aion(1:nspec_evolve)

       state_delm % xn = state % xn
       state_delm % T  = state % T * (ONE - eps)

       call actual_rhs(state_delm, ydotm)

       ydotm(1:nspec_evolve) = ydotm(1:nspec_evolve) * aion(1:nspec_evolve)

       do m = 1, neqs
          scratch = HALF*(ydotp(m) - ydotm(m)) / (eps * state % T)
          call set_jac_entry(jac, m, net_itemp, scratch)
       enddo

       ! energy derivatives -- these are all 0! (yay!)
       scratch = ZERO
       do m = 1, neqs
          call set_jac_entry(jac, m, net_ienuc, scratch)
       enddo

    else
       call copy_burn_t(state_delp, state)
       call copy_burn_t(state_delm, state)   ! note: delm here is actually just the input

       ! default
       call actual_rhs(state_delm, ydotm)

       ydotm(1:nspec_evolve) = ydotm(1:nspec_evolve) * aion(1:nspec_evolve)

       ! species derivatives
       do n = 1, nspec_evolve
          ! perturb species -- we send in X, but ydot is in terms of dY/dt, not dX/dt
          state_delp % xn = state % xn

          h = eps * abs(state % xn(n))
          if (h == 0) then
             h = eps
          endif

          state_delp % xn(n) = state % xn(n) + h

          call actual_rhs(state_delp, ydotp)

          ! We integrate X, so convert from the Y we got back from the RHS

          ydotp(1:nspec_evolve) = ydotp(1:nspec_evolve) * aion(1:nspec_evolve)

          do m = 1, neqs
             scratch = (ydotp(m) - ydotm(m)) / h
             call set_jac_entry(jac, m, n, scratch)
          enddo
       enddo

       ! temperature derivative
       state_delp % xn = state % xn

       h = eps * abs(state % T)
       if (h == 0) then
          h = eps
       endif

       state_delp % T = state % T + h

       call actual_rhs(state_delp, ydotp)

       ydotp(1:nspec_evolve) = ydotp(1:nspec_evolve) * aion(1:nspec_evolve)

       do m = 1, neqs
          scratch = (ydotp(m) - ydotm(m)) / h
          call set_jac_entry(jac, m, net_itemp, scratch)
       enddo

       ! energy derivatives -- these are all 0! (yay!)
       scratch = ZERO
       do m = 1, neqs
          call set_jac_entry(jac, m, net_ienuc, scratch)
       enddo

    endif

  end subroutine numerical_jac

#ifndef AMREX_USE_CUDA
  subroutine test_numerical_jac(state)
    ! compare the analytic Jacobian to the numerically differenced one

    use actual_rhs_module
    use eos_module, only : eos
    use eos_type_module, only : eos_t, eos_input_rt, normalize_abundances
    use jacobian_sparsity_module, only: get_jac_entry, scale_jac_entry

    type (burn_t) :: state
    type (burn_t) :: state_num
    type (eos_t) :: eos_state

    real(rt) :: scratch, scratch_num
    real(rt) :: ydot(neqs)
    real(rt) :: jac_analytic(neqs, neqs), jac_numerical(neqs, neqs)

    integer :: i, j, n
    character(len=16) :: namei, namej

    ! Set up state

    call burn_to_eos(state, eos_state)
    call normalize_abundances(eos_state)
    call eos(eos_input_rt, eos_state)
    call eos_to_burn(eos_state, state)

    state % self_heat = .true.

    state_num = state

    ! Evaluate the analytical Jacobian.
    call actual_jac(state, jac_analytic)

    ! the analytic Jacobian is in terms of Y, since that's what the
    ! nets work with, so we convert it to derivatives with respect to
    ! X and of mass fraction creation rates
    do n = 1, nspec_evolve
       do j = 1, neqs
          call scale_jac_entry(jac_analytic, n, j, aion(n))
          call scale_jac_entry(jac_analytic, j, n, aion_inv(n))
       enddo
    enddo

    ! Now compute the numerical Jacobian.
    call actual_rhs(state_num, ydot)
    call numerical_jac(state_num, jac_numerical)

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
          call get_jac_entry(jac_numerical, i, j, scratch_num)
          call get_jac_entry(jac_analytic, i, j, scratch)
          if (scratch_num /= scratch) then
             if (scratch_num /= ZERO) then
                write (*,999) trim(namei), &
                     trim(namej), scratch_num, scratch, &
                     (scratch_num - scratch)/scratch_num
             else
                write (*,999) trim(namei), &
                     trim(namej), scratch_num, scratch
             endif
          endif

       enddo
    enddo

  end subroutine test_numerical_jac
#endif
end module numerical_jac_module
