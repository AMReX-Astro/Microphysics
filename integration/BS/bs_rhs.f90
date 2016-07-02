module rhs_module

contains

  ! The rhs routine provides the right-hand-side for the BS solver.
  ! This is a generic interface that calls the specific RHS routine in the
  ! network you're actually using.

  subroutine f_rhs(bs)

    !$acc routine seq

    use eos_module
    use bl_types
    use bs_convert_module
    use burn_type_module
    use bl_constants_module, only: ZERO
    use actual_rhs_module, only: actual_rhs
    use extern_probin_module, only: call_eos_in_rhs, dT_crit, renormalize_abundances, &
                                    burning_mode, burning_mode_factor
    use rpar_indices
    use bs_type_module, only: bs_t, clean_state, renormalize_species
    use integration_data, only: aionInv

    implicit none

    type (bs_t) :: bs

    type (eos_t)  :: eos_state
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

    ! Several thermodynamic quantities come in via bs % upar -- note: these
    ! are evaluated at the start of the integration, so if things change
    ! dramatically, they will fall out of sync with the current
    ! thermodynamics.

    call bs_to_eos(eos_state, bs)

    ! If the temperature is smaller than the EOS can handle, allow it to
    ! reset the temperature accordingly.

    eos_state % reset = .true.

    ! Do not check the validity of inputs going into the EOS call.
    ! Sometimes we may stray into a meaningless state and we want
    ! to be able to get through the EOS call without a failure so
    ! that we can return to a meaningful state in the convergence.

    eos_state % check_inputs = .false.

    ! Evaluate the thermodynamics -- if desired. Note that
    ! even if this option is selected, we don't need to do it
    ! for non-self-heating integrations because the temperature
    ! isn't being updated. Also, if it is, we can optionally
    ! set a fraction dT_crit such that we don't call the EOS
    ! if the last temperature we evaluated the EOS at is relatively
    ! close to the current temperature.

    ! Otherwise just do the composition calculations since
    ! that's needed to construct dY/dt. Then make sure
    ! the abundances are safe.

    if (call_eos_in_rhs .and. bs % upar(irp_self_heat) > ZERO) then

       call eos(eos_input_rt, eos_state)

    else if (abs(eos_state % T - bs % upar(irp_Told)) > dT_crit * eos_state % T .and. bs % upar(irp_self_heat) > ZERO) then

       call eos(eos_input_rt, eos_state)

       bs % upar(irp_dcvdt) = (eos_state % cv - bs % upar(irp_cv)) / (eos_state % T - bs % upar(irp_Told))
       bs % upar(irp_dcpdt) = (eos_state % cp - bs % upar(irp_cp)) / (eos_state % T - bs % upar(irp_Told))
       bs % upar(irp_Told)  = eos_state % T

    else

       call composition(eos_state)

    endif

    call eos_to_bs(eos_state, bs)

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

    use actual_rhs_module, only: actual_jac
    use numerical_jac_module, only: numerical_jac
    use extern_probin_module, only: jacobian, burning_mode, burning_mode_factor
    use bs_convert_module
    use burn_type_module
    use bs_type_module

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
