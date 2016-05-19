module numerical_jac_module

  use network
  use burn_type_module

  implicit none

contains

  subroutine numerical_jac(state)

    !$acc routine seq

    use bl_types
    use bl_constants_module, only: ZERO
    use eos_module
    use actual_rhs_module, only: actual_rhs

    implicit none

    type (burn_t)    :: state

    integer          :: m, n

    type (burn_t)    :: state_del

    real(dp_t), parameter :: eps = 1.d-8

    state % jac(:,:) = ZERO

    ! default
    call actual_rhs(state)

    state_del = state

    ! species derivatives
    do n = 1, nspec
       ! perturb species -- we send in X, but ydot is in terms of dY/dt, not dX/dt
       state_del % xn = state % xn
       state_del % xn(n) = state % xn(n) * (ONE + eps)

       call actual_rhs(state_del)

       do m = 1, neqs
          state % jac(m,n) = (state_del % ydot(m) - state % ydot(m)) / (eps * state % xn(n) / aion(n))
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

  end subroutine numerical_jac

end module numerical_jac_module
