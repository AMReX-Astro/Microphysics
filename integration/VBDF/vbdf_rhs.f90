  ! The rhs routine provides the right-hand-side for the VBDF solver.
  ! This is a generic interface that calls the specific RHS routine in the
  ! network you're actually using.

  subroutine rhs(ts)

    !$acc routine seq

    use eos_module
    use bl_types
    use vbdf_convert_module
    use burn_type_module
    use bl_constants_module, only: ZERO, ONE
    use actual_rhs_module, only: actual_rhs
    use extern_probin_module, only: call_eos_in_rhs, dT_crit, renormalize_abundances
    use rpar_indices
    use bdf_type_module

    implicit none

    type (bdf_ts) :: ts

    type (eos_t)  :: eos_state
    type (burn_t) :: burn_state

    real(dp_t) :: nspec_sum

    ! Ensure that mass fractions always stay positive.

    ts % y(1:nspec_evolve,1) = max(ts % y(1:nspec_evolve,1) * aion(1:nspec_evolve), 1.d-200) / aion(1:nspec_evolve)

    ! Optionally, renormalize them so they sum to unity.

    if (renormalize_abundances) then
       nspec_sum = sum(ts % y(1:nspec_evolve,1) * aion(1:nspec_evolve)) + &
                   sum(ts % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1,1) * aion(nspec_evolve+1:nspec))

       ts % y(1:nspec_evolve,1) = ts % y(1:nspec_evolve,1) / nspec_sum
       ts % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1,1) = ts % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1,1) / nspec_sum
    endif

    ! We are integrating a system of
    !
    ! y(1:nspec)   = dX/dt
    ! y(net_itemp) = dT/dt
    ! y(net_ienuc) = denuc/dt

    ts % yd(:,1) = ZERO

    ! Several thermodynamic quantities come in via ts % upar -- note: these
    ! are evaluated at the start of the integration, so if things change
    ! dramatically, they will fall out of sync with the current
    ! thermodynamics.

    call vbdf_to_eos(eos_state, ts)

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

    if (call_eos_in_rhs .and. ts % upar(irp_self_heat,1) > ZERO) then

       call eos(eos_input_rt, eos_state)

    else if (abs(eos_state % T - ts % upar(irp_Told,1)) > dT_crit * eos_state % T .and. ts % upar(irp_self_heat,1) > ZERO) then

       call eos(eos_input_rt, eos_state)

       ts % upar(irp_dcvdt,1) = (eos_state % cv - ts % upar(irp_cv,1)) / (eos_state % T - ts % upar(irp_Told,1))
       ts % upar(irp_dcpdt,1) = (eos_state % cp - ts % upar(irp_cp,1)) / (eos_state % T - ts % upar(irp_Told,1))
       ts % upar(irp_Told,1)  = eos_state % T

    else

       call composition(eos_state)

    endif

    call eos_to_vbdf(eos_state, ts)

    burn_state % have_rates = .false.

    ! Call the specific network routine to get the RHS.

    call vbdf_to_burn(ts, burn_state)
    call actual_rhs(burn_state)
    call burn_to_vbdf(burn_state, ts)

  end subroutine rhs



  ! Analytical Jacobian

  subroutine jac(ts)

    !$acc routine seq

    use actual_rhs_module, only: actual_jac
    use numerical_jac_module, only: numerical_jac
    use extern_probin_module, only: jacobian
    use vbdf_convert_module
    use burn_type_module
    use bdf_type_module
    use bl_error_module

    implicit none

    type (bdf_ts) :: ts

    type (burn_t) :: state

    state % have_rates = .false.

    ! Call the specific network routine to get the Jacobian.

    call vbdf_to_burn(ts, state)

    if (jacobian == 1) then
       call actual_jac(state)
    else
       call numerical_jac(state)
    endif

    call burn_to_vbdf(state, ts)

  end subroutine jac
