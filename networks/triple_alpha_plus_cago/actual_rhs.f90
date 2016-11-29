module actual_rhs_module

  use network
  use burn_type_module
  use screen_module
  use rates_module
  use dydt_module
  use rate_type_module

  implicit none

  public

contains

  subroutine actual_rhs_init()

    implicit none

  end subroutine actual_rhs_init


  subroutine get_rates(state, rr)

    !$acc routine seq
    
    type (burn_t), intent(in) :: state
    type (rate_t), intent(out) :: rr

    double precision :: temp, dens
    double precision :: ymol(nspec)

    double precision :: rates(nrates), dratesdt(nrates)

    temp = state % T
    dens = state % rho
    ymol = state % xn * aion_inv

    call make_rates(temp, dens, rates, dratesdt)
    call screen(temp, dens, ymol, rates, dratesdt)
    
    rr % rates(1,:) = rates(:)
    rr % rates(2,:) = dratesdt(:)

    rr % T_eval = temp

  end subroutine get_rates
    

  subroutine actual_rhs(state)

    !$acc routine seq
    
    use temperature_integration_module, only: temperature_rhs

    implicit none

    type (burn_t) :: state
    type (rate_t) :: rr

    double precision :: ymol(nspec)
    integer :: k

    state % ydot = ZERO

    ymol = state % xn * aion_inv

    ! set up the species ODEs for the reaction network
    ! species inputs are in molar fractions but come out in mass fractions

    call get_rates(state, rr)

    call dydt(ymol, rr % rates(1,:), state % ydot(1:nspec_evolve))

    ! Energy generation rate

    call ener_gener_rate(state % ydot(1:nspec_evolve), state % ydot(net_ienuc))

    call temperature_rhs(state)

  end subroutine actual_rhs



  subroutine actual_jac(state)

    !$acc routine seq

    use burn_type_module, only : neqs
    use temperature_integration_module, only: temperature_jac

    implicit none

    type (burn_t) :: state
    type (rate_t) :: rr

    double precision :: ymol(nspec)
    double precision :: rates(nrates), dratesdt(nrates)

    integer :: i, j

    call get_rates(state, rr)

    rates(:)    = rr % rates(1,:)
    dratesdt(:) = rr % rates(2,:)

    ! initialize
    do j = 1, neqs
       state % jac(:,j) = ZERO
    enddo
    
    ymol = state % xn * aion_inv

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

    !$acc routine seq
    
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
