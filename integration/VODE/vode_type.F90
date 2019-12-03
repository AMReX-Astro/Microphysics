module vode_type_module

  use cuvode_parameters_module, only: VODE_NEQS, VODE_LMAX, VODE_LENWM
  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains
  
  subroutine clean_state(y, rpar)

    !$acc routine seq

    use amrex_constants_module, only: ONE
    use actual_network, only: aion, nspec, nspec_evolve
    use burn_type_module, only: neqs, net_itemp
    use vode_rpar_indices, only: n_rpar_comps
    use eos_type_module, only : eos_get_small_temp
    use extern_probin_module, only: renormalize_abundances, SMALL_X_SAFE, MAX_TEMP

    implicit none

    real(rt) :: y(neqs), rpar(n_rpar_comps)

    real(rt) :: small_temp

    !$gpu

    ! Ensure that mass fractions always stay positive and less than or equal to 1.

    y(1:nspec_evolve) = max(min(y(1:nspec_evolve), ONE), SMALL_X_SAFE)

    ! Renormalize the abundances as necessary.

    if (renormalize_abundances) then
       call renormalize_species(y, rpar)
    endif

    ! Ensure that the temperature always stays within reasonable limits.

    call eos_get_small_temp(small_temp)

    y(net_itemp) = min(MAX_TEMP, max(y(net_itemp), small_temp))

  end subroutine clean_state


  subroutine renormalize_species(y, rpar)

    !$acc routine seq
    
    use network, only: aion, aion_inv, nspec, nspec_evolve
    use burn_type_module, only: neqs
    use vode_rpar_indices, only: n_rpar_comps, irp_nspec, n_not_evolved

    implicit none

    real(rt) :: y(neqs), rpar(n_rpar_comps)

    real(rt) :: nspec_sum

    !$gpu

    nspec_sum = &
         sum(y(1:nspec_evolve)) + &
         sum(rpar(irp_nspec:irp_nspec+n_not_evolved-1))

    y(1:nspec_evolve) = y(1:nspec_evolve) / nspec_sum
    rpar(irp_nspec:irp_nspec+n_not_evolved-1) = &
         rpar(irp_nspec:irp_nspec+n_not_evolved-1) / nspec_sum

  end subroutine renormalize_species


  subroutine update_thermodynamics(y, rpar)

    !$acc routine seq
    
    use amrex_constants_module, only: ZERO
    use extern_probin_module, only: call_eos_in_rhs, dT_crit
    use eos_type_module, only: eos_t, eos_input_rt
    use eos_composition_module, only : composition
    use eos_module, only: eos
    use vode_rpar_indices, only: n_rpar_comps, irp_self_heat, irp_cp, irp_cv, irp_Told, irp_dcpdt, irp_dcvdt
    use burn_type_module, only: neqs

    implicit none

    real(rt) :: y(neqs), rpar(n_rpar_comps)

    type (eos_t) :: eos_state

    !$gpu

    ! Several thermodynamic quantities come in via rpar -- note: these
    ! are evaluated at the start of the integration, so if things change
    ! dramatically, they will fall out of sync with the current
    ! thermodynamics.

    call vode_to_eos(eos_state, y, rpar)

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

    if (call_eos_in_rhs .and. rpar(irp_self_heat) > ZERO) then

       call eos(eos_input_rt, eos_state)

    else if (abs(eos_state % T - rpar(irp_Told)) > dT_crit * eos_state % T .and. rpar(irp_self_heat) > ZERO) then

       call eos(eos_input_rt, eos_state)

       rpar(irp_dcvdt) = (eos_state % cv - rpar(irp_cv)) / &
            (eos_state % T - rpar(irp_Told))
       rpar(irp_dcpdt) = (eos_state % cp - rpar(irp_cp)) / &
            (eos_state % T - rpar(irp_Told))
       rpar(irp_Told)  = eos_state % T

       ! note: the update to rpar(irp_cv) and irp_cp is done
       ! in the call to eos_to_bs that follows this block.

    else

       call composition(eos_state)

    endif

    call eos_to_vode(eos_state, y, rpar)

  end subroutine update_thermodynamics



  ! Given an rpar array and the integration state, set up an EOS state.
  ! We could fill the energy component by storing the initial energy in
  ! rpar if we wanted, but we won't do that because (1) if we do call the EOS,
  ! it is always in (rho, T) mode and (2) converting back would imply subtracting
  ! off the nuclear energy from the zone's internal energy, which could lead to
  ! issues from roundoff error if the energy released from burning is small.

  subroutine vode_to_eos(state, y, rpar)

    !$acc routine seq

    use integrator_scaling_module, only: dens_scale, temp_scale
    use network, only: nspec, nspec_evolve, aion, aion_inv
    use eos_type_module, only: eos_t
    use vode_rpar_indices, only: irp_dens, irp_nspec, irp_cp, irp_cv, irp_abar, irp_zbar, &
                            irp_eta, irp_ye, irp_cs, n_rpar_comps, n_not_evolved
    use burn_type_module, only: neqs, net_itemp

    implicit none

    type (eos_t) :: state
    real(rt)   :: rpar(n_rpar_comps)
    real(rt)   :: y(neqs)

    !$gpu

    state % rho     = rpar(irp_dens) * dens_scale
    state % T       = y(net_itemp) * temp_scale

    state % xn(1:nspec_evolve) = y(1:nspec_evolve)
    state % xn(nspec_evolve+1:nspec) = &
         rpar(irp_nspec:irp_nspec+n_not_evolved-1)

    state % cp      = rpar(irp_cp)
    state % cv      = rpar(irp_cv)
    state % abar    = rpar(irp_abar)
    state % zbar    = rpar(irp_zbar)
    state % eta     = rpar(irp_eta)
    state % y_e     = rpar(irp_ye)
    state % cs      = rpar(irp_cs)

  end subroutine vode_to_eos



  ! Given an EOS state, fill the rpar and integration state data.

  subroutine eos_to_vode(state, y, rpar)

    !$acc routine seq

    use integrator_scaling_module, only: inv_dens_scale, inv_temp_scale
    use network, only: nspec, nspec_evolve, aion, aion_inv
    use eos_type_module, only: eos_t
    use vode_rpar_indices, only: irp_dens, irp_nspec, irp_cp, irp_cv, irp_abar, irp_zbar, &
                            irp_eta, irp_ye, irp_cs, n_rpar_comps, n_not_evolved
    use burn_type_module, only: neqs, net_itemp

    implicit none

    type (eos_t) :: state
    real(rt)   :: rpar(n_rpar_comps)
    real(rt)   :: y(neqs)

    !$gpu

    rpar(irp_dens) = state % rho * inv_dens_scale
    y(net_itemp) = state % T * inv_temp_scale

    y(1:nspec_evolve) = state % xn(1:nspec_evolve)
    rpar(irp_nspec:irp_nspec+n_not_evolved-1) = &
         state % xn(nspec_evolve+1:nspec)

    rpar(irp_cp)                    = state % cp
    rpar(irp_cv)                    = state % cv
    rpar(irp_abar)                  = state % abar
    rpar(irp_zbar)                  = state % zbar
    rpar(irp_eta)                   = state % eta
    rpar(irp_ye)                    = state % y_e
    rpar(irp_cs)                    = state % cs

  end subroutine eos_to_vode



  ! Given a burn state, fill the rpar and integration state data.

  subroutine burn_to_vode(state, y, rpar, ydot)

    !$acc routine seq

    use integrator_scaling_module, only: inv_dens_scale, inv_temp_scale, inv_ener_scale, temp_scale, ener_scale
    use amrex_constants_module, only: ONE
    use network, only: nspec, nspec_evolve
    use vode_rpar_indices, only: irp_dens, irp_nspec, irp_cp, irp_cv, irp_abar, irp_zbar, &
                            irp_ye, irp_eta, irp_cs, irp_dx, &
                            irp_Told, irp_dcvdt, irp_dcpdt, irp_self_heat, &
                            n_rpar_comps, n_not_evolved
    use burn_type_module, only: neqs, burn_t, net_itemp, net_ienuc

    implicit none

    type (burn_t) :: state
    real(rt)    :: rpar(n_rpar_comps)
    real(rt)    :: y(neqs)
    real(rt), optional :: ydot(neqs)

    !$gpu

    rpar(irp_dens) = state % rho * inv_dens_scale
    y(net_itemp) = state % T * inv_temp_scale

    y(1:nspec_evolve) = state % xn(1:nspec_evolve)
    rpar(irp_nspec:irp_nspec+n_not_evolved-1) = state % xn(nspec_evolve+1:nspec)

    y(net_ienuc)                             = state % e * inv_ener_scale

    rpar(irp_cp)                             = state % cp
    rpar(irp_cv)                             = state % cv
    rpar(irp_abar)                           = state % abar
    rpar(irp_zbar)                           = state % zbar
    rpar(irp_ye)                             = state % y_e
    rpar(irp_eta)                            = state % eta
    rpar(irp_cs)                             = state % cs
    rpar(irp_dx)                             = state % dx

    rpar(irp_Told)                           = state % T_old
    rpar(irp_dcvdt)                          = state % dcvdt
    rpar(irp_dcpdt)                          = state % dcpdt

    if (present(ydot)) then
       ydot = state % ydot
       ydot(net_itemp) = ydot(net_itemp) * inv_temp_scale
       ydot(net_ienuc) = ydot(net_ienuc) * inv_ener_scale
    endif

<<<<<<< HEAD
=======
    if (present(jac)) then
       jac = state % jac
       jac(net_itemp,:) = jac(net_itemp,:) * inv_temp_scale
       jac(net_ienuc,:) = jac(net_ienuc,:) * inv_ener_scale
       jac(:,net_itemp) = jac(:,net_itemp) * temp_scale
       jac(:,net_ienuc) = jac(:,net_ienuc) * ener_scale
    endif

>>>>>>> bye_bye_vode
    if (state % self_heat) then
       rpar(irp_self_heat) = ONE
    else
       rpar(irp_self_heat) = -ONE
    endif

  end subroutine burn_to_vode



  ! Given an rpar array and the integration state, set up a burn state.

  subroutine vode_to_burn(y, rpar, state)

    !$acc routine seq

    use integrator_scaling_module, only: dens_scale, temp_scale, ener_scale
    use amrex_constants_module, only: ZERO
    use network, only: nspec, nspec_evolve, aion, aion_inv
    use vode_rpar_indices, only: irp_dens, irp_nspec, irp_cp, irp_cv, irp_abar, irp_zbar, &
                            irp_ye, irp_eta, irp_cs, irp_dx, &
                            irp_Told, irp_dcvdt, irp_dcpdt, irp_self_heat, &
                            n_rpar_comps, n_not_evolved
    use burn_type_module, only: neqs, burn_t, net_itemp, net_ienuc

    implicit none

    type (burn_t) :: state
    real(rt)    :: rpar(n_rpar_comps)
    real(rt)    :: y(neqs)

    integer :: n

    !$gpu

    state % rho      = rpar(irp_dens) * dens_scale
    state % T        = y(net_itemp) * temp_scale
    state % e        = y(net_ienuc) * ener_scale

    state % xn(1:nspec_evolve) = y(1:nspec_evolve)
    state % xn(nspec_evolve+1:nspec) = &
         rpar(irp_nspec:irp_nspec+n_not_evolved-1)

    state % cp       = rpar(irp_cp)
    state % cv       = rpar(irp_cv)
    state % abar     = rpar(irp_abar)
    state % zbar     = rpar(irp_zbar)
    state % y_e      = rpar(irp_ye)
    state % eta      = rpar(irp_eta)
    state % cs       = rpar(irp_cs)
    state % dx       = rpar(irp_dx)

    state % T_old    = rpar(irp_Told)
    state % dcvdt    = rpar(irp_dcvdt)
    state % dcpdt    = rpar(irp_dcpdt)

    if (rpar(irp_self_heat) > ZERO) then
       state % self_heat = .true.
    else
       state % self_heat = .false.
    endif

  end subroutine vode_to_burn

end module vode_type_module
