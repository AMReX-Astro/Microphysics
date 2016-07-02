module rhs_module

contains

  ! The rhs routine provides the right-hand-side for the BS solver.
  ! This is a generic interface that calls the specific RHS routine in the
  ! network you're actually using.

  subroutine f_rhs(bs)

    !$acc routine seq

    use bl_types, only: dp_t
    use burn_type_module, only: burn_t, net_ienuc
    use bl_constants_module, only: ZERO
    use actual_rhs_module, only: actual_rhs
    use extern_probin_module, only: renormalize_abundances, burning_mode, burning_mode_factor
    use bs_type_module, only: bs_t, clean_state, renormalize_species, update_thermodynamics, &
                              burn_to_bs, bs_to_burn
    use integration_data, only: aionInv

    implicit none

    type (bs_t) :: bs

    type (burn_t) :: burn_state

    real(dp_t) :: limit_factor, t_sound, t_enuc

    ! We are integrating a system of
    !
    ! y(1:nspec)   = dX/dt
    ! y(net_itemp) = dT/dt
    ! y(net_ienuc) = denuc/dt

    bs % dydt(:) = ZERO

    ! Fix the state as necessary.

    call clean_state(bs)

    ! Renormalize abundances as necessary.

    if (renormalize_abundances) then
       call renormalize_species(bs)
    endif

    ! Update the thermodynamic quantities as necessary.

    call update_thermodynamics(bs)

    burn_state % have_rates = .false.

    ! Call the specific network routine to get the RHS.

    call bs_to_burn(bs, burn_state)
    call actual_rhs(burn_state)

    ! For burning_mode == 3, limit the rates.
    ! Note that this relies on burn_state % e being a relatively
    ! decent representation of the zone's current internal energy.

    if (burning_mode == 3) then

       t_enuc = burn_state % e / max(abs(burn_state % ydot(net_ienuc)), 1.d-50)
       t_sound = burn_state % dx / burn_state % cs

       limit_factor = min(1.0d0, burning_mode_factor * t_enuc / t_sound)

       burn_state % ydot = limit_factor * burn_state % ydot

    endif

    call burn_to_bs(burn_state, bs)

    ! Increment the evaluation counter.

    bs % n_rhs = bs % n_rhs + 1

  end subroutine f_rhs



  ! Analytical Jacobian

  subroutine jac(bs)

    !$acc routine seq

    use bl_types, only: dp_t
    use bl_constants_module, only: ZERO
    use actual_rhs_module, only: actual_jac
    use numerical_jac_module, only: numerical_jac
    use extern_probin_module, only: jacobian, burning_mode, burning_mode_factor
    use burn_type_module, only: burn_t, net_ienuc
    use bs_type_module, only: bs_t, bs_to_burn, burn_to_bs

    implicit none

    type (bs_t) :: bs

    type (burn_t) :: state

    real(dp_t) :: limit_factor, t_sound, t_enuc

    state % have_rates = .false.

    bs % jac(:,:) = ZERO

    ! Call the specific network routine to get the Jacobian.

    call bs_to_burn(bs, state)

    if (jacobian == 1) then
       call actual_jac(state)
    else
       call numerical_jac(state)
    endif

    ! For burning_mode == 3, limit the rates.
    ! Note that this relies on burn_state % e being a relatively
    ! decent representation of the zone's current internal energy.

    if (burning_mode == 3) then

       t_enuc = state % e / max(abs(state % ydot(net_ienuc)), 1.d-50)
       t_sound = state % dx / state % cs

       limit_factor = min(1.0d0, burning_mode_factor * t_enuc / t_sound)

       state % jac = limit_factor * state % jac

    endif

    call burn_to_bs(state, bs)

    ! Increment the evaluation counter.

    bs % n_jac = bs % n_jac + 1

  end subroutine jac

end module rhs_module
