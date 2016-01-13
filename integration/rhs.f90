  ! The f_rhs routine provides the right-hand-side for the DVODE solver.  It 
  ! converts the mass fractions into molar abundances before calling the 
  ! make_rates, screen, and dydt routines.  It also checks to see if the 
  ! temperature has changed much since the last call - if so, it updates the 
  ! temperature to get a better estimate of the reaction rates.
  !
  ! The jac routine provides an explicit Jacobian to the DVODE solver.

  subroutine f_rhs(neq, time, y, ydot, rpar, ipar)

    use eos_module
    use bl_types
    use rpar_indices
    use vode_indices
    use bl_constants_module, only: ZERO, ONE
    use actual_rhs_module, only: actual_rhs
    use extern_probin_module, only: call_eos_in_rhs

    implicit none

    integer,          intent(IN   ) :: neq, ipar
    double precision, intent(IN   ) :: time, y(neq)
    double precision, intent(INOUT) :: rpar(n_rpar_comps)
    double precision, intent(  OUT) :: ydot(neq)

    double precision :: rates(nrates)

    integer :: k

    type (eos_t) :: state

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

    call vode_to_eos(state, y, rpar)

    ! Evaluate the thermodynamics -- if desired. Note that
    ! even if this option is selected, we don't need to do it
    ! for non-self-heating integrations because the temperature
    ! isn't being updated.

    ! Otherwise just do the composition calculations since
    ! that's needed to construct dY/dt. Then make sure
    ! the abundances are safe.

    if (call_eos_in_rhs .and. rpar(irp_self_heat) > ZERO) then
       call eos(eos_input_rt, state)
    else
       call composition(state)
    endif

    call eos_to_vode(state, y, rpar)

    ! Call the specific network routine to get dY/dt and de/dt.

    call actual_rhs(neq,time,y,ydot,rpar)

  end subroutine f_rhs



  ! Analytical Jacobian

  subroutine jac(neq, time, y, ml, mu, pd, nrpd, rpar, ipar)

    use rpar_indices, only: n_rpar_comps
    use actual_rhs_module, only: actual_jac

    implicit none

    integer         , intent(IN   ) :: neq, ml, mu, nrpd, ipar
    double precision, intent(IN   ) :: y(neq), rpar(n_rpar_comps), time
    double precision, intent(  OUT) :: pd(neq,neq)

    call actual_jac(neq, time, y, pd, rpar)

  end subroutine jac



  ! Given an rpar array and the integration state, set up an EOS state.
  ! We could fill the energy component by storing the initial energy in
  ! rpar if we wanted, but we won't do that because (1) if we do call the EOS,
  ! it is always in (rho, T) mode and (2) converting back would imply subtracting
  ! off the nuclear energy from the zone's internal energy, which could lead to
  ! issues from roundoff error if the energy released from burning is small.

  subroutine vode_to_eos(state, y, rpar)

    use network, only: nspec, aion
    use eos_type_module, only: eos_t
    use rpar_indices
    use vode_module, only: neq
    use vode_indices

    implicit none

    type (eos_t)     :: state
    double precision :: rpar(n_rpar_comps)
    double precision :: y(neq)

    state % rho     = rpar(irp_dens)
    state % T       = y(net_itemp)
    state % xn(:)   = y(1:nspec) * aion(:)
    state % cp      = rpar(irp_cp)
    state % cv      = rpar(irp_cv)
    state % abar    = rpar(irp_abar)
    state % zbar    = rpar(irp_zbar)
    state % dhdX(:) = rpar(irp_dhdY:irp_dhdY-1+nspec) / aion(:)
    state % dedX(:) = rpar(irp_dedY:irp_dedY-1+nspec) / aion(:)

  end subroutine vode_to_eos



  ! Given an EOS state, fill the rpar and integration state data.

  subroutine eos_to_vode(state, y, rpar)

    use network, only: nspec, aion
    use eos_type_module, only: eos_t
    use rpar_indices
    use vode_module, only: neq
    use vode_indices

    implicit none

    type (eos_t)     :: state
    double precision :: rpar(n_rpar_comps)
    double precision :: y(neq)

    rpar(irp_dens)                  = state % rho
    y(net_itemp)                    = state % T
    y(1:nspec)                      = state % xn(:) / aion(:)
    rpar(irp_cp)                    = state % cp
    rpar(irp_cv)                    = state % cv
    rpar(irp_abar)                  = state % abar
    rpar(irp_zbar)                  = state % zbar
    rpar(irp_dhdY:irp_dhdY+nspec-1) = state % dhdX(:) * aion(:)
    rpar(irp_dedY:irp_dedY+nspec-1) = state % dedX(:) * aion(:)

  end subroutine eos_to_vode
