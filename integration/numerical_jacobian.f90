module numerical_jac_module

  use network
  use burn_type_module

  implicit none

contains

  subroutine numerical_jac(state)

    use bl_types
    use bl_constants_module, only: ZERO
    use eos_module
    use actual_rhs_module, only: actual_rhs

    implicit none

    type (burn_t)    :: state

    integer          :: m, n

    double precision :: rho, temp, abar, zbar
    double precision :: y(nspec)

    type (burn_t)    :: state_del

    double precision, parameter :: eps = 1.d-8

    state % jac(:,:) = ZERO

    ! default
    call actual_rhs(state)

    rho  = state % rho
    temp = state % T

    abar = state % abar
    zbar = state % zbar

    y    = state % xn / aion

    state_del % rho = rho
    state_del % T = temp
    state_del % xn = y * aion
    state_del % abar = abar
    state_del % zbar = zbar

    state_del % T_old = state % T_old
    state_del % cv = state % cv
    state_del % dcvdT = state % dcvdT
    state_del % cp = state % cp
    state_del % dcpdT = state % dcpdT



    ! species derivatives
    do n = 1, nspec
       ! perturb species -- we send in X, but ydot is in terms of dY/dt, not dX/dt
       state_del % xn = y * aion
       state_del % xn(n) = y(n)*aion(n)*(ONE + eps)

       call actual_rhs(state_del)

       do m = 1, neqs
          state % jac(m,n) = (state_del % ydot(m) - state % ydot(m))/(eps*y(n))
       enddo
    enddo

    ! temperature derivative
    state_del % xn = y * aion
    state_del % T = temp*(ONE + eps)

    call actual_rhs(state_del)

    do m = 1, neqs
       state % jac(m,net_itemp) = (state_del % ydot(m) - state % ydot(m))/(eps*temp)
    enddo

    ! energy derivatives -- these are all 0! (yay!)
    do m = 1, neqs
       state % jac(m,net_ienuc) = ZERO
    enddo

  end subroutine numerical_jac

end module numerical_jac_module
