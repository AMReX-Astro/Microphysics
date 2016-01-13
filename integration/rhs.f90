  ! The f_rhs routine provides the right-hand-side for the DVODE solver.
  ! This is a generic interface that calls the specific RHS routine in the
  ! network you're actually using.

  subroutine f_rhs(neq, time, y, ydot, rpar, ipar)

    use eos_module
    use bl_types
    use rpar_indices
    use vode_data
    use bl_constants_module, only: ZERO, ONE
    use actual_rhs_module, only: actual_rhs
    use extern_probin_module, only: call_eos_in_rhs, do_constant_volume_burn
    use rpar_indices

    implicit none

    integer,          intent(IN   ) :: neq, ipar
    double precision, intent(INOUT) :: time, y(neq)
    double precision, intent(INOUT) :: rpar(n_rpar_comps)
    double precision, intent(  OUT) :: ydot(neq)

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

    ! Undo the scaling for the user-facing routines.

    y(net_itemp)   = y(net_itemp) * temp_scale
    rpar(irp_dens) = rpar(irp_dens) * dens_scale

    ! Call the specific network routine to get dY/dt and de/dt.

    call actual_rhs(neq,time,y,ydot,rpar)

    ! Re-apply the scaling.

    y(net_itemp)   = y(net_itemp) / temp_scale
    rpar(irp_dens) = rpar(irp_dens) / dens_scale

    ydot(net_itemp) = ydot(net_itemp) / temp_scale

  end subroutine f_rhs



  ! Sets up the temperature equation. This should be called from
  ! within the actual_rhs routine but is provided here as a convenience
  ! since most networks will use the same temperature ODE.

  subroutine temperature_rhs(neq, y, ydot, rpar)

    use network, only: nspec
    use rpar_indices

    implicit none

    integer          :: neq
    double precision :: y(neq), ydot(neq)
    double precision :: rpar(n_rpar_comps)

    double precision :: cv, cp, dhdY(nspec), dedY(nspec)

    dhdY = rpar(irp_dhdY:irp_dhdY+nspec-1)
    dedY = rpar(irp_dedY:irp_dedY+nspec-1)

    cv   = rpar(irp_cv)
    cp   = rpar(irp_cp)

    if (rpar(irp_self_heat) > ZERO) then

       ! Set up the temperature ODE.  For constant pressure, Dp/Dt = 0, we
       ! evolve :
       !    dT/dt = (1/c_p) [ -sum_i (xi_i omega_i) + Hnuc]
       !
       ! For constant volume, div{U} = 0, and we evolve:
       !    dT/dt = (1/c_v) [ -sum_i ( {e_x}_i omega_i) + Hnuc]
       !
       ! See paper III, including Eq. A3 for details.

       if (do_constant_volume_burn) then
          ydot(net_itemp) = ( ydot(net_ienuc) - sum( dedY(:) * ydot(1:nspec) ) ) / cv
       else
          ydot(net_itemp) = ( ydot(net_ienuc) - sum( dhdY(:) * ydot(1:nspec) ) ) / cp
       endif

    endif

  end subroutine temperature_rhs



  ! Analytical Jacobian

  subroutine jac(neq, time, y, ml, mu, pd, nrpd, rpar, ipar)

    use rpar_indices
    use bl_constants_module, only: ZERO
    use actual_rhs_module, only: actual_jac
    use extern_probin_module, only: do_constant_volume_burn
    use vode_data
    use network, only: nspec
    use rpar_indices

    implicit none

    integer         , intent(IN   ) :: neq, ml, mu, nrpd, ipar
    double precision, intent(INOUT) :: y(neq), rpar(n_rpar_comps), time
    double precision, intent(  OUT) :: pd(neq,neq)

    integer          :: j

    ! Undo the scaling for the user-facing routines.

    y(net_itemp)   = y(net_itemp) * temp_scale
    rpar(irp_dens) = rpar(irp_dens) * dens_scale

    ! Call the specific network routine to get the Jacobian.

    call actual_jac(neq, time, y, pd, rpar)

    ! Re-apply the scaling.

    y(net_itemp)   = y(net_itemp) / temp_scale
    rpar(irp_dens) = rpar(irp_dens) / dens_scale

    pd(net_itemp,:) = pd(net_itemp,:) / temp_scale

  end subroutine jac



  ! Sets up the temperature entries in the Jacobian. This should be called from
  ! within the actual_jac routine but is provided here as a convenience
  ! since most networks will use the same temperature ODE.

  subroutine temperature_jac(neq, y, pd, rpar)

    use network, only: nspec
    use rpar_indices

    implicit none

    integer          :: neq
    double precision :: y(neq), pd(neq, neq)
    double precision :: rpar(n_rpar_comps)

    double precision :: cv, cp, dhdY(nspec), dedY(nspec)

    dhdY = rpar(irp_dhdY:irp_dhdY+nspec-1)
    dedY = rpar(irp_dedY:irp_dedY+nspec-1)

    cv   = rpar(irp_cv)
    cp   = rpar(irp_cp)

    ! Temperature Jacobian elements

    if (rpar(irp_self_heat) > ZERO) then

       if (do_constant_volume_burn) then

          ! d(itemp)/d(yi)
          do j = 1, nspec
             pd(net_itemp,j) = ( pd(net_ienuc,j) - sum( dEdY(:) * pd(1:nspec,j) ) ) / cv
          enddo

          ! d(itemp)/d(temp)
          pd(net_itemp,net_itemp) = ( pd(net_ienuc,net_itemp) - sum( dEdY(:) * pd(1:nspec,net_itemp) ) ) / cv

       else

          ! d(itemp)/d(yi)
          do j = 1, nspec
             pd(net_itemp,j) = ( pd(net_ienuc,j) - sum( dhdY(:) * pd(1:nspec,j) ) ) / cp
          enddo

          ! d(itemp)/d(temp)
          pd(net_itemp,net_itemp) = ( pd(net_ienuc,net_itemp) - sum( dhdY(:) * pd(1:nspec,net_itemp) ) ) / cp

       endif

    endif

  end subroutine temperature_jac



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
    use vode_data

    implicit none

    type (eos_t)     :: state
    double precision :: rpar(n_rpar_comps)
    double precision :: y(neq)

    state % rho     = rpar(irp_dens) * dens_scale
    state % T       = y(net_itemp) * temp_scale
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
    use vode_data

    implicit none

    type (eos_t)     :: state
    double precision :: rpar(n_rpar_comps)
    double precision :: y(neq)

    rpar(irp_dens)                  = state % rho / dens_scale
    y(net_itemp)                    = state % T / temp_scale
    y(1:nspec)                      = state % xn(:) / aion(:)
    rpar(irp_cp)                    = state % cp
    rpar(irp_cv)                    = state % cv
    rpar(irp_abar)                  = state % abar
    rpar(irp_zbar)                  = state % zbar
    rpar(irp_dhdY:irp_dhdY+nspec-1) = state % dhdX(:) * aion(:)
    rpar(irp_dedY:irp_dedY+nspec-1) = state % dedX(:) * aion(:)

  end subroutine eos_to_vode
