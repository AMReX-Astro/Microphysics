  ! The f_rhs routine provides the right-hand-side for the DVODE solver.
  ! This is a generic interface that calls the specific RHS routine in the
  ! network you're actually using.

  subroutine f_rhs(neq, time, y, ydot, rpar, ipar)

    use actual_network, only: aion, nspec_evolve
    use bl_types, only: dp_t
    use burn_type_module, only: burn_t, net_ienuc, net_itemp
    use bl_constants_module, only: ZERO, ONE
    use actual_rhs_module, only: actual_rhs
    use extern_probin_module, only: call_eos_in_rhs, dT_crit, renormalize_abundances, &
                                    burning_mode, burning_mode_factor, &
                                    integrate_temperature, integrate_energy, integrate_molar_fraction
    use vode_type_module, only: clean_state, renormalize_species, update_thermodynamics, &
                                burn_to_vode, vode_to_burn
    use rpar_indices, only: n_rpar_comps, irp_have_rates, irp_y_init, irp_t_sound

    implicit none

    integer,    intent(IN   ) :: neq, ipar
    real(dp_t), intent(INOUT) :: time, y(neq)
    real(dp_t), intent(INOUT) :: rpar(n_rpar_comps)
    real(dp_t), intent(  OUT) :: ydot(neq)

    type (burn_t) :: burn_state

    real(dp_t) :: limit_factor, t_sound, t_enuc

    ! We are integrating a system of
    !
    ! y(1:nspec_evolve) = dX/dt
    ! y(net_itemp) = dT/dt
    ! y(net_ienuc) = denuc/dt

    ydot = ZERO

    ! Fix the state as necessary.

    call clean_state(y, rpar)

    ! Renormalize the abundances as necessary.

    if (renormalize_abundances) then
       call renormalize_species(y, rpar)
    endif

    ! Update the thermodynamics as necessary.

    call update_thermodynamics(y, rpar)

    ! Call the specific network routine to get the RHS.

    rpar(irp_have_rates) = -ONE

    call vode_to_burn(y, rpar, burn_state)
    call actual_rhs(burn_state)

    ! Allow integration of X instead of Y.

    if (.not. integrate_molar_fraction) then

       burn_state % ydot(1:nspec_evolve) = burn_state % ydot(1:nspec_evolve) * aion(1:nspec_evolve)

    endif

    ! Allow temperature and energy integration to be disabled.

    if (.not. integrate_temperature) then

       burn_state % ydot(net_itemp) = ZERO

    endif

    if (.not. integrate_energy) then

       burn_state % ydot(net_ienuc) = ZERO

    endif

    ! For burning_mode == 3, limit the rates.
    ! Note that we are limiting with respect to the initial zone energy.

    if (burning_mode == 3) then

       t_enuc = rpar(irp_y_init + net_ienuc - 1) / max(abs(burn_state % ydot(net_ienuc)), 1.d-50)
       t_sound = rpar(irp_t_sound)

       limit_factor = min(1.0d0, burning_mode_factor * t_enuc / t_sound)

       burn_state % ydot = limit_factor * burn_state % ydot

    endif

    call burn_to_vode(burn_state, y, rpar, ydot = ydot)

  end subroutine f_rhs



  ! Analytical Jacobian

  subroutine jac(neq, time, y, ml, mu, pd, nrpd, rpar, ipar)

    use actual_network, only: aion, nspec_evolve
    use integration_data, only: aionInv
    use bl_constants_module, only: ZERO
    use actual_rhs_module, only: actual_jac
    use burn_type_module, only: burn_t, net_ienuc, net_itemp
    use vode_type_module, only: vode_to_burn, burn_to_vode
    use rpar_indices, only: n_rpar_comps, irp_y_init, irp_t_sound
    use bl_types, only: dp_t
    use extern_probin_module, only: burning_mode, burning_mode_factor, &
                                    integrate_temperature, integrate_energy, integrate_molar_fraction

    implicit none

    integer   , intent(IN   ) :: neq, ml, mu, nrpd, ipar
    real(dp_t), intent(INOUT) :: y(neq), rpar(n_rpar_comps), time
    real(dp_t), intent(  OUT) :: pd(neq,neq)

    type (burn_t) :: state
    real(dp_t) :: limit_factor, t_sound, t_enuc
    integer :: n

    ! Call the specific network routine to get the Jacobian.

    call vode_to_burn(y, rpar, state)
    call actual_jac(state)

    ! Allow integration of X instead of Y.

    if (.not. integrate_molar_fraction) then

       do n = 1, nspec_evolve
          state % jac(n,:) = state % jac(n,:) * aion(n)
          state % jac(:,n) = state % jac(:,n) * aionInv(n)
       enddo

    endif

    ! Allow temperature and energy integration to be disabled.

    if (.not. integrate_temperature) then

       state % jac(net_itemp,:) = ZERO

    endif

    if (.not. integrate_energy) then

       state % jac(net_ienuc,:) = ZERO

    endif

    ! For burning_mode == 3, limit the rates.
    ! Note that we are limiting with respect to the initial zone energy.

    if (burning_mode == 3) then

       t_enuc = rpar(irp_y_init + net_ienuc - 1) / max(abs(state % ydot(net_ienuc)), 1.d-50)
       t_sound = rpar(irp_t_sound)

       limit_factor = min(1.0d0, burning_mode_factor * t_enuc / t_sound)

       state % jac = limit_factor * state % jac

    endif

    call burn_to_vode(state, y, rpar, jac = pd)

  end subroutine jac
