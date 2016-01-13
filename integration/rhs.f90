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
    state % rho     = rpar(irp_dens)
    state % cp      = rpar(irp_cp)
    state % cv      = rpar(irp_cv)
    state % abar    = rpar(irp_abar)
    state % zbar    = rpar(irp_zbar)
    state % xn(:)   = y(1:nspec) * aion(:)
    state % dhdX(:) = rpar(irp_dhdX:irp_dhdX-1+nspec)
    state % dedX(:) = rpar(irp_dedX:irp_dedX-1+nspec)

    ! Temperature is one of the quantities that we are integrating --
    ! always use the current T.
    state % T       = y(net_itemp)

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

    rpar(irp_dhdX:irp_dhdX-1+nspec) = state % dhdX
    rpar(irp_dedX:irp_dedX-1+nspec) = state % dedX
    rpar(irp_cp) = state % cp
    rpar(irp_cv) = state % cv
    rpar(irp_abar) = state % abar
    rpar(irp_zbar) = state % zbar

    ! Call the specific network routine to get dY/dt and de/dt.

    call actual_rhs(neq,time,state,y,ydot,rpar)

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
