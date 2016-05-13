module actual_rhs_module

  use bl_types
  use bl_constants_module
  use network
  use burn_type_module
  use actual_burner_data, only: ener_gener_rate
  use temperature_integration_module, only: temperature_rhs, temperature_jac

  implicit none

contains

  subroutine actual_rhs(state)

    use extern_probin_module, only: f_act, T_burn_ref, rho_burn_ref, rtilde, nu

    implicit none

    type (burn_t)    :: state

    double precision :: xfueltmp
    double precision :: dens, temp, rate, y(nspec)

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

    state % ydot(1:nspec_evolve) = state % ydot(1:nspec_evolve) / aion(1:nspec_evolve)

    call ener_gener_rate(state % ydot(1:nspec_evolve), state % ydot(net_ienuc))

    call temperature_rhs(state)

  end subroutine actual_rhs



  ! At present the analytical Jacobian is not implemented.

  subroutine actual_jac(state)

    implicit none

    type (burn_t) :: state

    state % jac(:,:) = ZERO

  end subroutine actual_jac

end module actual_rhs_module
