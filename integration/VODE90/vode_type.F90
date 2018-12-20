module vode_type_module

  use cuvode_parameters_module, only: VODE_NEQS, VODE_LMAX, VODE_LENWM
  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains
  
  subroutine clean_state(y, rpar)

    !$acc routine seq
    
    use actual_network, only: aion, nspec, nspec_evolve
    use burn_type_module, only: neqs
    use rpar_indices, only: n_rpar_comps
    use extern_probin_module, only: integrate_relative_species

    implicit none

    real(rt) :: y(neqs), rpar(n_rpar_comps)

    !$gpu

    ! Ensure that mass fractions always stay positive.
    if (integrate_relative_species) then
       y(1:nspec_evolve) = max(y(1:nspec_evolve), -1.d0)
    else
       y(1:nspec_evolve) = max(y(1:nspec_evolve), 1.d-200)
    end if

  end subroutine clean_state

  subroutine renormalize_species(y, rpar)

    !$acc routine seq
    
    use network, only: aion, aion_inv, nspec, nspec_evolve
    use burn_type_module, only: neqs
    use rpar_indices, only: n_rpar_comps, irp_nspec, n_not_evolved, irp_y_init
    use extern_probin_module, only: integrate_relative_species

    implicit none

    real(rt) :: y(neqs), rpar(n_rpar_comps)

    real(rt) :: nspec_sum

    integer :: i

    !$gpu

    if (integrate_relative_species) then
       do i = 1, nspec_evolve
          y(i) = y(i) * rpar(irp_y_init - 1 + i) + rpar(irp_y_init - 1 + i)
       end do
    end if

    nspec_sum = &
         sum(y(1:nspec_evolve)) + &
         sum(rpar(irp_nspec:irp_nspec+n_not_evolved-1))

    y(1:nspec_evolve) = y(1:nspec_evolve) / nspec_sum
    rpar(irp_nspec:irp_nspec+n_not_evolved-1) = &
         rpar(irp_nspec:irp_nspec+n_not_evolved-1) / nspec_sum

    if (integrate_relative_species) then
       do i = 1, nspec_evolve
          y(i) = (y(i) - rpar(irp_y_init - 1 + i))/rpar(irp_y_init - 1 + i)
       end do
    end if

  end subroutine renormalize_species

  subroutine update_thermodynamics(y, rpar)

    !$acc routine seq
    
    use amrex_constants_module, only: ZERO
    use extern_probin_module, only: call_eos_in_rhs, dt_crit
    use eos_type_module, only: eos_t, eos_input_rt, composition
    use eos_module, only: eos
    use rpar_indices, only: n_rpar_comps, irp_self_heat, irp_cp, irp_cv, irp_Told, irp_dcpdt, irp_dcvdt
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


  subroutine initialize_relative_species(eos_state, rpar)

    use extern_probin_module, only: relatively_small_x
    use eos_type_module, only: eos_t
    use network, only: nspec_evolve
    use rpar_indices, only: n_rpar_comps, irp_y_init

    implicit none

    type (eos_t), intent(inout) :: eos_state
    real(rt), intent(inout) :: rpar(n_rpar_comps)    
    
    integer :: i

    !$gpu

    print *, 'in initialize_relative_species...'
    print *, 'got eos_state % xn = ', eos_state % xn
    print *, 'using relatively_small_x = ', relatively_small_x
    
    ! Use initial mass fractions no lower than relatively_small_x
    eos_state % xn(i) = max(eos_state % xn(i), relatively_small_x)

    ! Save the evolved species mass fractions for integrating relative X
    do i = 1, nspec_evolve
       rpar(irp_y_init - 1 + i) = eos_state % xn(i)
    end do

    print *, 'set y init to = ', rpar(irp_y_init:irp_y_init + nspec_evolve - 1)

  end subroutine initialize_relative_species


  subroutine store_initial_state(dvode_state)

    ! Save the initial integration vector
    !
    ! If we're integrating relative mass fractions
    ! then initialize_relative_species should have set
    ! the absolute initial mass fractions already and
    ! in this routine we save only the initial enuc and temp.

    use extern_probin_module, only: integrate_relative_species
    use network, only: nspec_evolve
    use rpar_indices, only: n_rpar_comps, irp_y_init
    use burn_type_module, only: net_ienuc, net_itemp
    use cuvode_types_module, only: dvode_t

    implicit none

    type (dvode_t), intent(inout) :: dvode_state
    
    integer :: i

    !$gpu

    if (.not. integrate_relative_species) then
       do i = 1, nspec_evolve
          dvode_state % rpar(irp_y_init - 1 + i) = dvode_state % y(i)
       end do
    end if

    dvode_state % rpar(irp_y_init - 1 + net_ienuc) = dvode_state % y(net_ienuc)
    dvode_state % rpar(irp_y_init - 1 + net_itemp) = dvode_state % y(net_itemp)

    ! print out the initial integration vector for debugging
    print *, 'integration vector at t=0: ', dvode_state % y
  end subroutine store_initial_state


  ! Given an rpar array and the integration state, set up an EOS state.
  ! We could fill the energy component by storing the initial energy in
  ! rpar if we wanted, but we won't do that because (1) if we do call the EOS,
  ! it is always in (rho, T) mode and (2) converting back would imply subtracting
  ! off the nuclear energy from the zone's internal energy, which could lead to
  ! issues from roundoff error if the energy released from burning is small.
  subroutine vode_to_eos(state, y, rpar)

    !$acc routine seq

    use integrator_scaling_module, only: dens_scale, temp_scale
    use extern_probin_module, only: integrate_relative_species
    use network, only: nspec, nspec_evolve, aion, aion_inv
    use eos_type_module, only: eos_t
    use rpar_indices, only: irp_dens, irp_nspec, irp_cp, irp_cv, irp_abar, irp_zbar, &
                            irp_eta, irp_ye, irp_cs, n_rpar_comps, n_not_evolved, irp_y_init
    use burn_type_module, only: neqs, net_itemp

    implicit none

    type (eos_t) :: state
    real(rt)   :: rpar(n_rpar_comps)
    real(rt)   :: y(neqs)

    integer :: i

    !$gpu

    state % rho     = rpar(irp_dens) * dens_scale
    state % T       = y(net_itemp) * temp_scale

    if (integrate_relative_species) then
       do i = 1, nspec_evolve
          state % xn(i) = y(i) * rpar(irp_y_init - 1 + i) + rpar(irp_y_init - 1 + i)
       end do
    else
       state % xn(1:nspec_evolve) = y(1:nspec_evolve)
    end if

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
    use extern_probin_module, only: integrate_relative_species
    use network, only: nspec, nspec_evolve, aion, aion_inv
    use eos_type_module, only: eos_t
    use rpar_indices, only: irp_dens, irp_nspec, irp_cp, irp_cv, irp_abar, irp_zbar, &
                            irp_eta, irp_ye, irp_cs, n_rpar_comps, n_not_evolved, irp_y_init
    use burn_type_module, only: neqs, net_itemp

    implicit none

    type (eos_t) :: state
    real(rt)   :: rpar(n_rpar_comps)
    real(rt)   :: y(neqs)

    integer :: i

    !$gpu

    rpar(irp_dens) = state % rho * inv_dens_scale
    y(net_itemp) = state % T * inv_temp_scale

    if (integrate_relative_species) then
       do i = 1, nspec_evolve
          print *, 'i = ', i
          y(i) = (state % xn(i) - rpar(irp_y_init - 1 + i))/rpar(irp_y_init - 1 + i)
          print *, 'xn(i) = ', state % xn(i)
          print *, 'y init (i) = ', rpar(irp_y_init - 1 + i)
       end do
    else
       y(1:nspec_evolve) = state % xn(1:nspec_evolve)
    end if
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
  subroutine burn_to_vode(state, y, rpar, ydot, jac)

    !$acc routine seq

    use integrator_scaling_module, only: inv_dens_scale, inv_temp_scale, inv_ener_scale, temp_scale, ener_scale
    use extern_probin_module, only: integrate_relative_species
    use amrex_constants_module, only: ONE
    use network, only: nspec, nspec_evolve, aion, aion_inv
    use rpar_indices, only: irp_dens, irp_nspec, irp_cp, irp_cv, irp_abar, irp_zbar, &
                            irp_ye, irp_eta, irp_cs, irp_dx, &
                            irp_Told, irp_dcvdt, irp_dcpdt, irp_self_heat, &
                            n_rpar_comps, n_not_evolved, irp_y_init
    use burn_type_module, only: neqs, burn_t, net_itemp, net_ienuc

    implicit none

    type (burn_t) :: state
    real(rt)    :: rpar(n_rpar_comps)
    real(rt)    :: y(neqs)
    real(rt), optional :: ydot(neqs), jac(neqs, neqs)

    integer :: i

    !$gpu

    rpar(irp_dens) = state % rho * inv_dens_scale
    y(net_itemp) = state % T * inv_temp_scale

    if (integrate_relative_species) then
       do i = 1, nspec_evolve
          y(i) = (state % xn(i) - rpar(irp_y_init - 1 + i))/rpar(irp_y_init - 1 + i)
       end do
    else       
       y(1:nspec_evolve) = state % xn(1:nspec_evolve)
    end if

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
       
       if (integrate_relative_species) then
          do i = 1, nspec_evolve
             ydot(i) = ydot(i)/rpar(irp_y_init - 1 + i)
          end do
       end if
       
       ydot(net_itemp) = ydot(net_itemp) * inv_temp_scale
       ydot(net_ienuc) = ydot(net_ienuc) * inv_ener_scale
    endif

    if (present(jac)) then
       jac = state % jac

       if (integrate_relative_species) then
          do i = 1, nspec_evolve
             jac(i,:) = jac(i,:) / rpar(irp_y_init - 1 + i)
             jac(:,i) = jac(:,i) * rpar(irp_y_init - 1 + i)
          end do
       end if

       jac(net_itemp,:) = jac(net_itemp,:) * inv_temp_scale
       jac(net_ienuc,:) = jac(net_ienuc,:) * inv_ener_scale
       jac(:,net_itemp) = jac(:,net_itemp) * temp_scale
       jac(:,net_ienuc) = jac(:,net_ienuc) * ener_scale
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
    use extern_probin_module, only: integrate_relative_species
    use amrex_constants_module, only: ZERO
    use network, only: nspec, nspec_evolve, aion, aion_inv
    use rpar_indices, only: irp_dens, irp_nspec, irp_cp, irp_cv, irp_abar, irp_zbar, &
                            irp_ye, irp_eta, irp_cs, irp_dx, &
                            irp_Told, irp_dcvdt, irp_dcpdt, irp_self_heat, &
                            n_rpar_comps, n_not_evolved, irp_y_init
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

    if (integrate_relative_species) then
       do n = 1, nspec_evolve
          state % xn(n) = y(n) * rpar(irp_y_init - 1 + n) + rpar(irp_y_init - 1 + n)
       end do
    else
       state % xn(1:nspec_evolve) = y(1:nspec_evolve)
    end if

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
