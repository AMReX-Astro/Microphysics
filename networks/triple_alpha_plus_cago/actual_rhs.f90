module actual_rhs_module

  use network
  use burn_type_module
  use screen_module
  use rates_module
  use dydt_module

  implicit none

  public

contains

  subroutine actual_rhs_init()

    implicit none

  end subroutine actual_rhs_init



  subroutine actual_rhs(state)

    use temperature_integration_module, only: temperature_rhs

    implicit none

    type (burn_t) :: state

    integer :: k

    double precision :: temp, dens
    double precision :: rates(nrates), dratesdt(nrates)
    double precision :: ymol(nspec)

    state % ydot = ZERO

    temp = state % T
    dens = state % rho
    ymol = state % xn / aion

    ! set up the species ODEs for the reaction network
    ! species inputs are in molar fractions but come out in mass fractions
    call make_rates(temp, dens, rates, dratesdt)

    call screen(temp, dens, ymol, rates, dratesdt)

    call dydt(ymol, rates, state % ydot(1:nspec_evolve))

    state % rates(1,:) = rates(:)
    state % rates(2,:) = dratesdt(:)

    ! Energy generation rate

    call ener_gener_rate(state % ydot(1:nspec_evolve), state % ydot(net_ienuc))

    call temperature_rhs(state)

  end subroutine actual_rhs



  subroutine actual_jac(state)

    use temperature_integration_module, only: temperature_jac

    implicit none

    type (burn_t) :: state

    double precision :: ymol(nspec)
    double precision :: rates(nrates), dratesdt(nrates)

    integer :: i, j

    rates(:)    = state % rates(1,:)
    dratesdt(:) = state % rates(2,:)

    ! initialize
    state % jac = ZERO

    ymol = state % xn / aion

    ! ======================================================================
    ! THESE ARE IN TERMS OF MOLAR FRACTIONS

    ! helium jacobian elements
    state % jac(ihe4,ihe4)  = - NINE * ymol(ihe4) * ymol(ihe4) * rates(ir3a) &
                              - ONE * ymol(ic12) * rates(ircago)
    state % jac(ihe4,ic12)  = - ONE * ymol(ihe4) * rates(ircago)

    ! carbon jacobian elements
    state % jac(ic12,ihe4) =   THREE * ymol(ihe4) * ymol(ihe4) * rates(ir3a) &
                             - ONE * ymol(ic12) * rates(ircago)
    state % jac(ic12,ic12) = - ONE * ymol(ihe4) * rates(ircago)

    ! oxygen jacobian elements
    state % jac(io16,ihe4) = ONE * ymol(ic12) * rates(ircago)
    state % jac(io16,ic12) = ONE * ymol(ihe4) * rates(ircago)

    ! ======================================================================

    ! Add the temperature derivatives: df(y_i) / dT

    call dydt(ymol, dratesdt, state % jac(1:nspec_evolve,net_itemp))

    ! Energy generation rate Jacobian elements with respect to species

    do j = 1, nspec_evolve
       call ener_gener_rate(state % jac(1:nspec_evolve,j), state % jac(net_ienuc,j))
    enddo

    ! Jacobian elements with respect to temperature

    call ener_gener_rate(state % jac(1:nspec_evolve,net_itemp), state % jac(net_ienuc,net_itemp))

    call temperature_jac(state)

  end subroutine actual_jac



  subroutine ener_gener_rate(dydt, enuc)

    use network

    implicit none

    double precision :: dydt(nspec_evolve), enuc

    enuc = -sum(dydt(:) * aion(1:nspec_evolve) * ebin(1:nspec_evolve))

  end subroutine ener_gener_rate

  subroutine update_unevolved_species(state)

    !$acc routine seq

    implicit none

    type (burn_t)    :: state

    ! although we nspec_evolve < nspec, we never change the Fe56
    ! abundance, so there is no algebraic relation we need to
    ! enforce here.

  end subroutine update_unevolved_species

end module actual_rhs_module
