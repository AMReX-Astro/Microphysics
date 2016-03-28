  ! The f_rhs routine provides the right-hand-side for the DVODE solver.
  ! This is a generic interface that calls the specific RHS routine in the
  ! network you're actually using.

  subroutine f_rhs(neq, time, y, ydot, rpar, ipar)

    use eos_module
    use bl_types
    use rpar_indices
    use vode_convert_module
    use burn_type_module
    use bl_constants_module, only: ZERO, ONE
    use actual_rhs_module, only: actual_rhs
    use extern_probin_module, only: call_eos_in_rhs, renormalize_abundances
    use rpar_indices

    implicit none

    integer,          intent(IN   ) :: neq, ipar
    double precision, intent(INOUT) :: time, y(neq)
    double precision, intent(INOUT) :: rpar(n_rpar_comps)
    double precision, intent(  OUT) :: ydot(neq)

    type (eos_t)  :: eos_state
    type (burn_t) :: burn_state

    ! Ensure that mass fractions always stay positive.

    y(1:nspec) = max(y(1:nspec) * aion, 1.d-200) / aion

    ! Optionally, renormalize them so they sum to unity.

    if (renormalize_abundances) then
       y(1:nspec) = y(1:nspec) / sum(y(1:nspec) * aion)
    endif

    ! We are integrating a system of
    !
    ! y(1:nspec)   = dX/dt
    ! y(net_itemp) = dT/dt
    ! y(net_ienuc) = denuc/dt

    ydot = ZERO

    ! Several thermodynamic quantities come in via rpar -- note: these
    ! are evaluated at the start of the integration, so if things change
    ! dramatically, they will fall out of sync with the current
    ! thermodynamics.

    call vode_to_eos(eos_state, y, rpar)

    ! Evaluate the thermodynamics -- if desired. Note that
    ! even if this option is selected, we don't need to do it
    ! for non-self-heating integrations because the temperature
    ! isn't being updated.

    ! Otherwise just do the composition calculations since
    ! that's needed to construct dY/dt. Then make sure
    ! the abundances are safe.

    if (call_eos_in_rhs .and. rpar(irp_self_heat) > ZERO) then
       call eos(eos_input_rt, eos_state)
    else
       call composition(eos_state)
    endif

    call eos_to_vode(eos_state, y, rpar)

    ! Call the specific network routine to get the RHS.

    call vode_to_burn(y, rpar, burn_state)
    call actual_rhs(burn_state)
    call burn_to_vode(burn_state, y, rpar, ydot = ydot)

  end subroutine f_rhs



  ! Analytical Jacobian

  subroutine jac(neq, time, y, ml, mu, pd, nrpd, rpar, ipar)

    use rpar_indices
    use bl_constants_module, only: ZERO
    use actual_rhs_module, only: actual_jac
    use vode_convert_module
    use burn_type_module
    use network, only: nspec
    use rpar_indices

    implicit none

    integer         , intent(IN   ) :: neq, ml, mu, nrpd, ipar
    double precision, intent(INOUT) :: y(neq), rpar(n_rpar_comps), time
    double precision, intent(  OUT) :: pd(neq,neq)

    type (burn_t) :: state

    ! Call the specific network routine to get the Jacobian.

    call vode_to_burn(y, rpar, state)
    call actual_jac(state)
    call burn_to_vode(state, y, rpar, jac = pd)

  end subroutine jac



