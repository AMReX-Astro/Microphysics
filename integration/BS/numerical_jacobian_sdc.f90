module numerical_jac_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

  subroutine numerical_jac(state)

    !$acc routine seq

    use bs_rhs_module, only: f_bs_rhs
    use extern_probin_module, only: centered_diff_jac
    use sdc_type_module, only: sdc_t, SVAR_EVOLVE
    use bs_type_module, only: bs_t
    use amrex_constants_module, only: ZERO, HALF, ONE
    use amrex_fort_module, only : rt => amrex_real

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (bs_t)    :: state
    type (bs_t)    :: state_delp, state_delm

    integer        :: m, n

    ! the choice of eps should be ~ sqrt(eps), where eps is machine epsilon. 
    ! this balances truncation vs. roundoff error in the differencing
    real(rt), parameter :: eps = 1.e-8_rt, SMALL = 1.e-8_rt

    state % jac(:,:) = ZERO

    ! default
    call f_bs_rhs(state)

    if (centered_diff_jac) then

       state_delp = state
       state_delm = state

       do n = 1, SVAR_EVOLVE

          state_delp % y = state % y
          state_delp % y(n) = state % y(n) * (ONE + eps)

          call f_bs_rhs(state_delp)

          state_delm % y = state % y
          state_delm % y(n) = state % y(n) * (ONE - eps)

          call f_bs_rhs(state_delm)

          do m = 1, SVAR_EVOLVE

             state % jac(m,n) = (state_delp % ydot(m) - state_delm % ydot(m)) / ((state_delp % y(n) - state_delm % y(n)) + SMALL)

          enddo

       enddo

    else

       state_delp = state

       do n = 1, SVAR_EVOLVE

          state_delp % y = state % y
          state_delp % y(n) = state % y(n) * (ONE + eps)

          call f_bs_rhs(state_delp)

          do m = 1, SVAR_EVOLVE

             state % jac(m,n) = (state_delp % ydot(m) - state % ydot(m)) / ((state_delp % y(n) - state % y(n)) + SMALL)

          enddo

       enddo

    endif

  end subroutine numerical_jac

end module numerical_jac_module
