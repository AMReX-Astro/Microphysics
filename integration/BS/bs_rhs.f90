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
    use bl_constants_module, only: ZERO, ONE
    use actual_rhs_module, only: actual_rhs
    use extern_probin_module, only: call_eos_in_rhs, dT_crit, renormalize_abundances
    use rpar_indices
    use bs_type_module
    use integration_data, only: aionInv

    implicit none

    type (bs_t) :: bs

    type (eos_t)  :: eos_state
    type (burn_t) :: burn_state

    real(dp_t) :: nspec_sum

    ! Ensure that mass fractions always stay positive.

    bs % y(1:nspec_evolve) = max(bs % y(1:nspec_evolve) * aion(1:nspec_evolve), 1.d-30) * aionInv(1:nspec_evolve)
    bs % y(1:nspec_evolve) = min(bs % y(1:nspec_evolve) * aion(1:nspec_evolve), ONE) * aionInv(1:nspec_evolve)

    ! Ensure that the temperature always stays within reasonable limits.

    bs % y(net_itemp) = min(1.0d11 / temp_scale, max(bs % y(net_itemp), 1.0d4 / temp_scale))

    ! Optionally, renormalize them so they sum to unity.

    if (renormalize_abundances) then
       nspec_sum = sum(bs % y(1:nspec_evolve) * aion(1:nspec_evolve)) + &
                   sum(bs % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) * aion(nspec_evolve+1:nspec))

       bs % y(1:nspec_evolve) = bs % y(1:nspec_evolve) / nspec_sum
       bs % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) = bs % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) / nspec_sum
    endif

    ! We are integrating a system of
    !
    ! y(1:nspec)   = dX/dt
    ! y(net_itemp) = dT/dt
    ! y(net_ienuc) = denuc/dt

    bs % dydt(:) = ZERO

    ! Several thermodynamic quantities come in via bs % upar -- note: these
    ! are evaluated at the start of the integration, so if things change
    ! dramatically, they will fall out of sync with the current
    ! thermodynamics.

    call bs_to_eos(eos_state, bs)

    ! If the temperature is smaller than the EOS can handle, allow it to
    ! reset the temperature accordingly.

    eos_state % reset = .true.

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
    call burn_to_bs(burn_state, bs)

    ! Increment the evaluation counter.

    bs % n_rhs = bs % n_rhs + 1

  end subroutine f_rhs



  ! Analytical Jacobian

  subroutine jac(bs)

    !$acc routine seq

    use actual_rhs_module, only: actual_jac
    use numerical_jac_module, only: numerical_jac
    use extern_probin_module, only: jacobian
    use bs_convert_module
    use burn_type_module
    use bs_type_module

    implicit none

    type (bs_t) :: bs

    type (burn_t) :: state

    state % have_rates = .false.

    bs % jac(:,:) = ZERO

    ! Call the specific network routine to get the Jacobian.

    call bs_to_burn(bs, state)

    if (jacobian == 1) then
       call actual_jac(state)
    else
       call numerical_jac(state)
    endif

    call burn_to_bs(state, bs)

    ! Increment the evaluation counter.

    bs % n_jac = bs % n_jac + 1

  end subroutine jac

end module rhs_module
