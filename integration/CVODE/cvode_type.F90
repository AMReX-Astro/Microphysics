module cvode_type_module

  use amrex_fort_module, only: rt => amrex_real
  use burn_type_module, only: neqs

  implicit none

  integer, parameter :: VODE_NEQS = neqs
  
contains


  subroutine sk_clean_state(y, rpar) bind(C, name="sk_clean_state")

    !$acc routine seq
    
    use actual_network, only: aion, nspec, nspec_evolve
    use burn_type_module, only: neqs
    use cvode_rpar_indices, only: n_rpar_comps

    implicit none

    real(rt) :: y(neqs), rpar(n_rpar_comps)

    !$gpu

    ! Ensure that mass fractions always stay positive.

    y(1:nspec_evolve) = max(y(1:nspec_evolve), 1.e-200_rt)

  end subroutine sk_clean_state


  subroutine sk_renormalize_species(y, rpar) bind(C, name="sk_renormalize_species")

    !$acc routine seq
    
    use network, only: aion, aion_inv, nspec, nspec_evolve
    use burn_type_module, only: neqs
    use cvode_rpar_indices, only: n_rpar_comps, irp_nspec, n_not_evolved
    use extern_probin_module, only: renormalize_abundances

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

  end subroutine sk_renormalize_species

  
  subroutine sk_update_thermodynamics(y, rpar) bind(C, name="sk_update_thermodynamics")

    !$acc routine seq
    
    use amrex_constants_module, only: ZERO
    use extern_probin_module, only: call_eos_in_rhs, dt_crit
    use eos_type_module, only: eos_t, eos_input_rt, composition
    use eos_module, only: eos
    use cvode_rpar_indices, only: n_rpar_comps, irp_self_heat, irp_cx, irp_Told, irp_dcxdt
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

    else if (abs(eos_state % T - rpar(irp_Told)) > dt_crit * eos_state % T .and. rpar(irp_self_heat) > ZERO) then

       call eos(eos_input_rt, eos_state)

       if (do_constant_volume_burn) then
          rpar(irp_dcxdt) = (eos_state % cv - rpar(irp_cx)) / &
               (eos_state % T - rpar(irp_Told))
       else
          rpar(irp_dcxdt) = (eos_state % cp - rpar(irp_cx)) / &
               (eos_state % T - rpar(irp_Told))
       end if
       rpar(irp_Told)  = eos_state % T

       ! note: the update to rpar(irp_cv) and irp_cp is done
       ! in the call to eos_to_bs that follows this block.

    else

       call composition(eos_state)

    endif

    call eos_to_vode(eos_state, y, rpar)

  end subroutine sk_update_thermodynamics



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
    use cvode_rpar_indices, only: irp_dens, irp_nspec, irp_cx, irp_abar, irp_zbar, &
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

    state % cx      = rpar(irp_cx)
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
    use cvode_rpar_indices, only: irp_dens, irp_nspec, irp_cx, irp_abar, irp_zbar, &
                            irp_eta, irp_ye, irp_cs, n_rpar_comps, n_not_evolved
    use burn_type_module, only: neqs, net_itemp, net_ienuc

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

    rpar(irp_cx)                    = state % cx
    rpar(irp_abar)                  = state % abar
    rpar(irp_zbar)                  = state % zbar
    rpar(irp_eta)                   = state % eta
    rpar(irp_ye)                    = state % y_e
    rpar(irp_cs)                    = state % cs

  end subroutine eos_to_vode



  ! Given a burn state, fill the rpar and integration state data.
  subroutine burn_to_vode(state, y, rpar, ydot, jac)

    !$acc routine seq

    use integrator_scaling_module, only: inv_dens_scale, inv_temp_scale, inv_ener_scale, temp_scale, ener_scale
    use amrex_constants_module, only: ONE
    use network, only: nspec, nspec_evolve, aion, aion_inv, NETWORK_SPARSE_JAC_NNZ
    use cvode_rpar_indices, only: irp_dens, irp_nspec, irp_cx, irp_abar, irp_zbar, &
                            irp_ye, irp_eta, irp_cs, irp_dx, &
                            irp_Told, irp_dcxdt, irp_self_heat, &
                            n_rpar_comps, n_not_evolved
    use burn_type_module, only: neqs, burn_t, net_itemp, net_ienuc
    use jacobian_sparsity_module, only: scale_csr_jac_entry

    implicit none

    type (burn_t) :: state
    real(rt)    :: rpar(n_rpar_comps)
    real(rt)    :: y(neqs)
    real(rt), optional :: ydot(neqs)
    
#ifdef REACT_SPARSE_JACOBIAN
    real(rt), optional :: jac(NETWORK_SPARSE_JAC_NNZ)
#else
    real(rt), optional :: jac(neqs, neqs)
#endif

    integer :: i

    !$gpu

    rpar(irp_dens) = state % rho * inv_dens_scale
    y(net_itemp) = state % T * inv_temp_scale

    y(1:nspec_evolve) = state % xn(1:nspec_evolve)
    rpar(irp_nspec:irp_nspec+n_not_evolved-1) = state % xn(nspec_evolve+1:nspec)

    y(net_ienuc)                             = state % e * inv_ener_scale

    rpar(irp_cx)                             = state % cx
    rpar(irp_abar)                           = state % abar
    rpar(irp_zbar)                           = state % zbar
    rpar(irp_ye)                             = state % y_e
    rpar(irp_eta)                            = state % eta
    rpar(irp_cs)                             = state % cs
    rpar(irp_dx)                             = state % dx

    rpar(irp_Told)                           = state % T_old
    rpar(irp_dcxdt)                          = state % dcxdt

    if (present(ydot)) then
       ydot = state % ydot
       ydot(net_itemp) = ydot(net_itemp) * inv_temp_scale
       ydot(net_ienuc) = ydot(net_ienuc) * inv_ener_scale
    endif

    if (present(jac)) then
#ifdef REACT_SPARSE_JACOBIAN
       jac = state % sparse_jac
       do i = 1, neqs
          call scale_csr_jac_entry(jac, net_itemp, i, inv_temp_scale)
          call scale_csr_jac_entry(jac, net_ienuc, i, inv_ener_scale)
          call scale_csr_jac_entry(jac, i, net_itemp, temp_scale)
          call scale_csr_jac_entry(jac, i, net_ienuc, ener_scale)
       enddo
#else
       jac = state % jac
       jac(net_itemp,:) = jac(net_itemp,:) * inv_temp_scale
       jac(net_ienuc,:) = jac(net_ienuc,:) * inv_ener_scale
       jac(:,net_itemp) = jac(:,net_itemp) * temp_scale
       jac(:,net_ienuc) = jac(:,net_ienuc) * ener_scale
#endif
    endif

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
    use cvode_rpar_indices, only: irp_dens, irp_nspec, irp_cx, irp_abar, irp_zbar, &
                            irp_ye, irp_eta, irp_cs, irp_dx, &
                            irp_Told, irp_dcxdt, irp_self_heat, &
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

    state % cx       = rpar(irp_cx)
    state % abar     = rpar(irp_abar)
    state % zbar     = rpar(irp_zbar)
    state % y_e      = rpar(irp_ye)
    state % eta      = rpar(irp_eta)
    state % cs       = rpar(irp_cs)
    state % dx       = rpar(irp_dx)

    state % T_old    = rpar(irp_Told)
    state % dcxdt    = rpar(irp_dcxdt)

    if (rpar(irp_self_heat) > ZERO) then
       state % self_heat = .true.
    else
       state % self_heat = .false.
    endif

  end subroutine vode_to_burn

  
  subroutine sk_initialize_cell(y, rpar) bind(C, name="sk_initialize_cell")

    use integrator_scaling_module, only: inv_dens_scale, inv_temp_scale, inv_ener_scale
    use amrex_constants_module, only: ONE, ZERO
    use extern_probin_module, only: dT_crit, do_constant_volume_burn
    use burn_type_module, only: neqs, burn_t, net_itemp, net_ienuc, burn_to_eos
    use burn_type_module, only: normalize_abundances_burn
    use temperature_integration_module, only: self_heat
    use eos_type_module, only: eos_t, eos_input_rt, copy_eos_t
    use eos_module, only: eos
    use network, only: nspec_evolve, aion_inv
    use cvode_rpar_indices, only: n_rpar_comps, irp_self_heat, irp_t_sound, irp_dx, irp_dens, &
                            irp_Told, irp_dcxdt, irp_y_init, irp_energy_offset

    implicit none

    real(rt), intent(inout) :: y(neqs), rpar(n_rpar_comps)
    
    type(burn_t) :: burn_state
    type(eos_t)  :: eos_state, eos_state_temp

    !$gpu

    ! Prescale the temperature, energy, density
    y(net_itemp) = y(net_itemp) * inv_temp_scale
    y(net_ienuc) = y(net_ienuc) * inv_ener_scale
    rpar(irp_dens) = rpar(irp_dens) * inv_dens_scale

    ! Do the EOS call
    call vode_to_burn(y, rpar, burn_state)

    ! Normalize the abundances
    call normalize_abundances_burn(burn_state)

    ! Call the EOS
    call burn_to_eos(burn_state, eos_state)
    call eos(eos_input_rt, eos_state)
    call eos_to_vode(eos_state, y, rpar)

    ! Initialize specific internal energy
    rpar(irp_energy_offset) = eos_state % e * inv_ener_scale
    y(net_ienuc) = eos_state % e * inv_ener_scale

    ! Detect burning mode
    if (self_heat) then
       rpar(irp_self_heat) = ONE
    else
       rpar(irp_self_heat) = -ONE
    endif

    ! Set sound crossing time
    rpar(irp_t_sound) = rpar(irp_dx)/eos_state % cs

    ! Set up for dT_crit functionality so we can use a linear
    ! interpolation of the specific heat in between EOS calls
    rpar(irp_Told) = eos_state % T

    if (dT_crit < 1.0e19_rt) then

       call copy_eos_t(eos_state_temp, eos_state)
       eos_state_temp % T = eos_state % T * (ONE + sqrt(epsilon(ONE)))

       call eos(eos_input_rt, eos_state_temp)

       if (do_constant_volume_burn) then
          rpar(irp_dcxdt) = (eos_state_temp % cv - eos_state % cv) / &
               (eos_state_temp % T - eos_state % T)
       else
          rpar(irp_dcxdt) = (eos_state_temp % cp - eos_state % cp) / &
               (eos_state_temp % T - eos_state % T)
       end if

    endif
    
    ! Save the initial state.
    rpar(irp_y_init:irp_y_init + neqs - 1) = y
    
  end subroutine sk_initialize_cell


  subroutine sk_finalize_cell(y, rpar) bind(C, name="sk_finalize_cell")

    use integrator_scaling_module, only: dens_scale, temp_scale, ener_scale
    use extern_probin_module, only: burning_mode
    use network, only: nspec, nspec_evolve
    use actual_rhs_module, only : update_unevolved_species
    use burn_type_module, only: neqs, burn_t, net_ienuc, net_itemp, normalize_abundances_burn
    use cvode_rpar_indices, only: n_rpar_comps, irp_energy_offset, irp_dens

    implicit none

    real(rt), intent(inout) :: y(neqs), rpar(n_rpar_comps)
    
    type(burn_t) :: burn_state

    !$gpu

    ! Subtract off the initial energy offset to return
    ! the total *generated* energy as y(net_ienuc)
    y(net_ienuc) = y(net_ienuc) - rpar(irp_energy_offset)

    ! Convert to burn state out
    call vode_to_burn(y, rpar, burn_state)
    
    ! Update unevolved species
    if (nspec_evolve < nspec) then
       call update_unevolved_species(burn_state)
    endif

    ! TODO: support burning mode 3 limiting
    
    ! Normalize abundances
    call normalize_abundances_burn(burn_state)

    ! Convert back to vode data
    call burn_to_vode(burn_state, y, rpar)

    ! Rescale the temperature, energy, density
    y(net_itemp) = y(net_itemp) * temp_scale
    y(net_ienuc) = y(net_ienuc) * ener_scale
    rpar(irp_dens) = rpar(irp_dens) * dens_scale

  end subroutine sk_finalize_cell
  
end module cvode_type_module
