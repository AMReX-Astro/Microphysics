module actual_rhs_module

  use network
  use burn_type_module
  use temperature_integration_module, only: temperature_rhs, temperature_jac
  use microphysics_type_module, only: rt, ZERO

  implicit none

contains

  subroutine actual_rhs_init()

    implicit none

  end subroutine actual_rhs_init



  subroutine actual_rhs(state)

    use extern_probin_module, only: f_act, T_burn_ref, rho_burn_ref, rtilde, nu

    implicit none

    type (burn_t)    :: state

    real(rt) :: xfueltmp
    real(rt) :: dens, temp, rate, y(nspec)

    state % ydot = ZERO

    xfueltmp = max(state % xn(ifuel_), ZERO)
    dens     = state % rho
    temp     = state % T

    ! Rate is expressed in mass fraction form

    if (temp < f_act * T_burn_ref) then
       rate = ZERO
    else
       rate = rtilde * (dens/rho_burn_ref) * xfueltmp**2 * (temp/T_burn_ref)**nu
    endif

    state % ydot(ifuel_)  = -rate
    state % ydot(iash_)   =  rate

    ! Convert back to molar form

    state % ydot(1:nspec_evolve) = state % ydot(1:nspec_evolve) * aion_inv(1:nspec_evolve)

    call ener_gener_rate(state % ydot(1:nspec_evolve), state % ydot(net_ienuc))

    call temperature_rhs(state)

  end subroutine actual_rhs



  ! At present the analytical Jacobian is not implemented.

  subroutine actual_jac(state)

    implicit none

    type (burn_t) :: state

    state % jac(:,:) = ZERO

  end subroutine actual_jac



  ! Computes the instantaneous energy generation rate

  subroutine ener_gener_rate(dydt, enuc)

    use network

    implicit none

    real(rt) :: dydt(nspec_evolve), enuc

    ! This is basically e = m c**2

    enuc = sum(dydt(:) * aion(1:nspec_evolve) * ebin(1:nspec_evolve))

  end subroutine ener_gener_rate

  subroutine update_unevolved_species(state)

    !$acc routine seq

    implicit none

    type (burn_t)    :: state

  end subroutine update_unevolved_species

end module actual_rhs_module
