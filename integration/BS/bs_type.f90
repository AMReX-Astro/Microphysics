module bs_type_module

  use bl_types, only: dp_t
  use burn_type_module, only: neqs
  use rpar_indices, only: n_rpar_comps

  implicit none

  ! BS parameters -- see the discussion in 16.4
  integer, parameter :: KMAXX = 7
  integer :: nseq(KMAXX+1)

  ! error codes
  integer, parameter :: IERR_NONE = 0
  integer, parameter :: IERR_DT_TOO_SMALL = -100
  integer, parameter :: IERR_TOO_MANY_STEPS = -101
  integer, parameter :: IERR_DT_UNDERFLOW = -102
  integer, parameter :: IERR_NO_CONVERGENCE = -103

  integer, parameter :: IERR_LU_DECOMPOSITION_ERROR = -200

  real(kind=dp_t), parameter :: S1 = 0.25_dp_t
  real(kind=dp_t), parameter :: S2 = 0.7_dp_t

  real(kind=dp_t), parameter :: RED_BIG_FACTOR = 0.7_dp_t
  real(kind=dp_t), parameter :: RED_SMALL_FACTOR = 1.e-5_dp_t
  real(kind=dp_t), parameter :: SCALMX = 0.1_dp_t

  type bs_t
     logical :: first
     real(kind=dp_t) :: eps_old
     real(kind=dp_t) :: dt_did
     real(kind=dp_t) :: dt_next
     real(kind=dp_t) :: a(KMAXX+1)
     real(kind=dp_t) :: alpha(KMAXX, KMAXX)
     real(kind=dp_t) :: t_new
     integer :: kmax
     integer :: kopt

     real(kind=dp_t) :: y(neqs), dydt(neqs), jac(neqs, neqs)
     real(kind=dp_t) :: atol(neqs), rtol(neqs)
     real(kind=dp_t) :: upar(n_rpar_comps)
     real(kind=dp_t) :: t, dt, tmax
     integer         :: n
     integer         :: n_rhs, n_jac

     integer :: i, j, k

  end type bs_t

  !$acc declare create(nseq)

contains

  subroutine clean_state(state)

    !$acc routine seq

    use bl_constants_module, only: ONE
    use extern_probin_module, only: small_x, integrate_molar_fraction
    use actual_network, only: aion, nspec, nspec_evolve
    use integration_data, only: aionInv, temp_scale
    use burn_type_module, only: net_itemp
    use eos_module, only : eos_get_small_temp

    implicit none

    ! this should be larger than any reasonable temperature we will encounter
    real (kind=dp_t), parameter :: MAX_TEMP = 1.0d11

    ! this is the absolute cutoff for species -- note that this might
    ! be larger than small_x that the user set, but the issue is that
    ! we can have underflow issues if the integrator has to keep track
    ! of species mass fractions much smaller than this.
    real (kind=dp_t), parameter :: SMALL_X_SAFE = 1.0d-30
    real (kind=dp_t) :: small_temp

    type (bs_t) :: state

    ! Ensure that mass fractions always stay positive.
    if (integrate_molar_fraction) then
       state % y(1:nspec_evolve) = &
            max(min(state % y(1:nspec_evolve) * aion(1:nspec_evolve), ONE), &
                SMALL_X_SAFE) * aionInv(1:nspec_evolve)
    else
       state % y(1:nspec_evolve) = &
            max(min(state % y(1:nspec_evolve), ONE), SMALL_X_SAFE)
    endif


    ! Ensure that the temperature always stays within reasonable limits.
    call eos_get_small_temp(small_temp)

    state % y(net_itemp) = min(MAX_TEMP / temp_scale, &
                               max(state % y(net_itemp), small_temp / temp_scale))

  end subroutine clean_state


  subroutine renormalize_species(state)

    !$acc routine seq

    use actual_network, only: aion, nspec, nspec_evolve
    use rpar_indices, only: irp_nspec
    use extern_probin_module, only: integrate_molar_fraction

    implicit none

    type (bs_t) :: state

    real(dp_t) :: nspec_sum

    if (integrate_molar_fraction) then
       nspec_sum = &
            sum(state % y(1:nspec_evolve) * aion(1:nspec_evolve)) + &
            sum(state % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) * aion(nspec_evolve+1:nspec))
    else
       nspec_sum = &
            sum(state % y(1:nspec_evolve)) + &
            sum(state % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1))
    endif

    state % y(1:nspec_evolve) = state % y(1:nspec_evolve) / nspec_sum
    state % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) = &
         state % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) / nspec_sum

  end subroutine renormalize_species


  subroutine update_thermodynamics(state)

    !$acc routine seq

    use bl_constants_module, only: ZERO
    use eos_type_module, only: eos_t, composition
    use eos_module, only: eos_input_rt, eos
    use extern_probin_module, only: call_eos_in_rhs, dT_crit
    use rpar_indices, only: irp_self_heat, irp_cv, irp_cp, irp_dcvdt, irp_dcpdt, irp_Told

    implicit none

    type (bs_t) :: state
    type (eos_t) :: eos_state

    ! Several thermodynamic quantities come in via bs % upar -- note: these
    ! are evaluated at the start of the integration, so if things change
    ! dramatically, they will fall out of sync with the current
    ! thermodynamics.

    call bs_to_eos(eos_state, state)

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

    if (call_eos_in_rhs .and. state % upar(irp_self_heat) > ZERO) then

       call eos(eos_input_rt, eos_state)

    else if (abs(eos_state % T - state % upar(irp_Told)) > dT_crit * eos_state % T .and. state % upar(irp_self_heat) > ZERO) then

       call eos(eos_input_rt, eos_state)

       state % upar(irp_dcvdt) = (eos_state % cv - state % upar(irp_cv)) / &
            (eos_state % T - state % upar(irp_Told))
       state % upar(irp_dcpdt) = (eos_state % cp - state % upar(irp_cp)) / &
            (eos_state % T - state % upar(irp_Told))
       state % upar(irp_Told)  = eos_state % T

       ! note: the update to state % upar(irp_cv) and irp_cp is done
       ! in the call to eos_to_bs that follows this block.
    else

       call composition(eos_state)

    endif

    call eos_to_bs(eos_state, state)

  end subroutine update_thermodynamics



  ! Given a bs_t, set up an EOS state.
  ! We could fill the energy component by storing the initial energy in
  ! ts % upar if we wanted, but we won't do that because (1) if we do call the EOS,
  ! it is always in (rho, T) mode and (2) converting back would imply subtracting
  ! off the nuclear energy from the zone's internal energy, which could lead to
  ! issues from roundoff error if the energy released from burning is small.

  subroutine bs_to_eos(state, bs)

    !$acc routine seq

    use actual_network, only: nspec, nspec_evolve, aion
    use integration_data, only: aionInv, dens_scale, temp_scale
    use eos_type_module, only: eos_t
    use rpar_indices, only: irp_dens, irp_nspec, irp_cp, irp_cv, irp_abar, irp_zbar, &
                            irp_eta, irp_ye, irp_cs, irp_dhdY, irp_dedY
    use burn_type_module, only: net_itemp
    use extern_probin_module, only: integrate_molar_fraction

    implicit none

    type (eos_t) :: state
    type (bs_t) :: bs

    state % rho     = bs % upar(irp_dens) * dens_scale
    state % T       = bs % y(net_itemp) * temp_scale

    if (integrate_molar_fraction) then
       state % xn(1:nspec_evolve) = bs % y(1:nspec_evolve) * aion(1:nspec_evolve)
       state % xn(nspec_evolve+1:nspec) = &
            bs % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) * aion(nspec_evolve+1:nspec)
    else
       state % xn(1:nspec_evolve) = bs % y(1:nspec_evolve)
       state % xn(nspec_evolve+1:nspec) = &
            bs % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1)
    endif

    state % cp      = bs % upar(irp_cp)
    state % cv      = bs % upar(irp_cv)
    state % abar    = bs % upar(irp_abar)
    state % zbar    = bs % upar(irp_zbar)
    state % eta     = bs % upar(irp_eta)
    state % y_e     = bs % upar(irp_ye)
    state % cs      = bs % upar(irp_cs)

    if (integrate_molar_fraction) then
       state % dhdX(:) = bs % upar(irp_dhdY:irp_dhdY-1+nspec) * aionInv(:)
       state % dedX(:) = bs % upar(irp_dedY:irp_dedY-1+nspec) * aionInv(:)
    else
       state % dhdX(:) = bs % upar(irp_dhdY:irp_dhdY-1+nspec)
       state % dedX(:) = bs % upar(irp_dedY:irp_dedY-1+nspec)
    endif

  end subroutine bs_to_eos



  ! Given an EOS state, fill the rpar and integration state data.

  subroutine eos_to_bs(state, bs)

    !$acc routine seq

    use actual_network, only: nspec, nspec_evolve, aion
    use integration_data, only: aionInv, inv_dens_scale, inv_temp_scale
    use eos_type_module, only: eos_t
    use rpar_indices, only: irp_dens, irp_nspec, irp_cp, irp_cv, irp_abar, irp_zbar, &
                            irp_eta, irp_ye, irp_cs, irp_dhdY, irp_dedY
    use integration_data, only: temp_scale, dens_scale
    use burn_type_module, only: net_itemp
    use extern_probin_module, only: integrate_molar_fraction

    implicit none

    type (eos_t) :: state
    type (bs_t) :: bs

    bs % upar(irp_dens) = state % rho * inv_dens_scale
    bs % y(net_itemp) = state % T * inv_temp_scale

    if (integrate_molar_fraction) then
       bs % y(1:nspec_evolve) = state % xn(1:nspec_evolve) * aionInv(1:nspec_evolve)
       bs % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) = &
            state % xn(nspec_evolve+1:nspec) * aionInv(nspec_evolve+1:nspec)
    else
       bs % y(1:nspec_evolve) = state % xn(1:nspec_evolve)
       bs % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) = &
            state % xn(nspec_evolve+1:nspec) 
    endif

    bs % upar(irp_cp)                    = state % cp
    bs % upar(irp_cv)                    = state % cv
    bs % upar(irp_abar)                  = state % abar
    bs % upar(irp_zbar)                  = state % zbar
    bs % upar(irp_eta)                   = state % eta
    bs % upar(irp_ye)                    = state % y_e
    bs % upar(irp_cs)                    = state % cs

    if (integrate_molar_fraction) then
       bs % upar(irp_dhdY:irp_dhdY+nspec-1) = state % dhdX(:) * aion(:)
       bs % upar(irp_dedY:irp_dedY+nspec-1) = state % dedX(:) * aion(:)
    else
       bs % upar(irp_dhdY:irp_dhdY+nspec-1) = state % dhdX(:)
       bs % upar(irp_dedY:irp_dedY+nspec-1) = state % dedX(:)
    endif

  end subroutine eos_to_bs



  ! Given a burn state, fill the bs_t.

  subroutine burn_to_bs(state, bs)

    !$acc routine seq

    use actual_network, only: nspec, nspec_evolve, aion, nrates
    use integration_data, only: aionInv, inv_dens_scale, inv_temp_scale, inv_ener_scale
    use rpar_indices, only: irp_dens, irp_nspec, irp_cp, irp_cv, irp_abar, irp_zbar, &
                            irp_ye, irp_eta, irp_cs, irp_dx, irp_dhdY, irp_dedY, &
                            irp_self_heat, irp_have_rates, irp_rates, irp_Told, irp_dcvdt, irp_dcpdt
    use burn_type_module, only: burn_t, net_itemp, net_ienuc, num_rate_groups
    use bl_constants_module, only: ONE
    use extern_probin_module, only: integrate_molar_fraction

    implicit none

    type (burn_t) :: state
    type (bs_t) :: bs

    integer :: n

    bs % upar(irp_dens) = state % rho * inv_dens_scale
    bs % y(net_itemp) = state % T * inv_temp_scale

    if (integrate_molar_fraction) then
       bs % y(1:nspec_evolve) = state % xn(1:nspec_evolve) * aionInv(1:nspec_evolve)
       bs % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) = &
            state % xn(nspec_evolve+1:nspec) * aionInv(nspec_evolve+1:nspec)
    else
       bs % y(1:nspec_evolve) = state % xn(1:nspec_evolve)
       bs % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) = &
            state % xn(nspec_evolve+1:nspec) 
    endif

    bs % y(net_ienuc)                             = state % e * inv_ener_scale
    bs % upar(irp_cp)                             = state % cp
    bs % upar(irp_cv)                             = state % cv
    bs % upar(irp_abar)                           = state % abar
    bs % upar(irp_zbar)                           = state % zbar
    bs % upar(irp_ye)                             = state % y_e
    bs % upar(irp_eta)                            = state % eta
    bs % upar(irp_cs)                             = state % cs
    bs % upar(irp_dx)                             = state % dx

    if (integrate_molar_fraction) then
       bs % upar(irp_dhdY:irp_dhdY+nspec-1) = state % dhdX(:) * aion(:)
       bs % upar(irp_dedY:irp_dedY+nspec-1) = state % dedX(:) * aion(:)
    else
       bs % upar(irp_dhdY:irp_dhdY+nspec-1) = state % dhdX(:)
       bs % upar(irp_dedY:irp_dedY+nspec-1) = state % dedX(:)
    endif

    bs % upar(irp_Told)                           = state % T_old
    bs % upar(irp_dcvdt)                          = state % dcvdt
    bs % upar(irp_dcpdt)                          = state % dcpdt

    bs % dydt = state % ydot
    bs % dydt(net_itemp) = bs % dydt(net_itemp) * inv_temp_scale
    bs % dydt(net_ienuc) = bs % dydt(net_ienuc) * inv_ener_scale

    do n = 1, neqs
       bs % jac(:,n) = state % jac(:,n)
    enddo

    bs % jac(net_itemp,:) = bs % jac(net_itemp,:) * inv_temp_scale
    bs % jac(net_ienuc,:) = bs % jac(net_ienuc,:) * inv_ener_scale

    if (state % have_rates) then
       bs % upar(irp_have_rates) = ONE
    else
       bs % upar(irp_have_rates) = -ONE
    endif

    do n = 1, nrates
       bs % upar(irp_rates+(n-1)*num_rate_groups:irp_rates+n*num_rate_groups-1) = &
            state % rates(:,n)
    enddo

    if (state % self_heat) then
       bs % upar(irp_self_heat) = ONE
    else
       bs % upar(irp_self_heat) = -ONE
    endif

  end subroutine burn_to_bs



  ! Given a bs_t, set up a burn state.

  subroutine bs_to_burn(bs, state)

    !$acc routine seq

    use actual_network, only: nspec, nspec_evolve, aion, nrates
    use integration_data, only: aionInv, dens_scale, temp_scale, ener_scale
    use rpar_indices, only: irp_dens, irp_nspec, irp_cp, irp_cv, irp_abar, irp_zbar, &
                            irp_ye, irp_eta, irp_cs, irp_dx, irp_dhdY, irp_dedY, &
                            irp_self_heat, irp_have_rates, irp_rates, irp_Told, irp_dcvdt, irp_dcpdt
    use burn_type_module, only: burn_t, net_itemp, net_ienuc, num_rate_groups
    use bl_constants_module, only: ZERO, ONE
    use extern_probin_module, only: integrate_molar_fraction

    implicit none

    type (burn_t) :: state
    type (bs_t) :: bs

    integer :: n

    state % rho      = bs % upar(irp_dens) * dens_scale
    state % T        = bs % y(net_itemp) * temp_scale
    state % e        = bs % y(net_ienuc) * ener_scale

    if (integrate_molar_fraction) then
       state % xn(1:nspec_evolve) = bs % y(1:nspec_evolve) * aion(1:nspec_evolve)
       state % xn(nspec_evolve+1:nspec) = &
            bs % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) * aion(nspec_evolve+1:nspec)
    else
       state % xn(1:nspec_evolve) = bs % y(1:nspec_evolve)
       state % xn(nspec_evolve+1:nspec) = &
            bs % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1)
    endif

    state % cp       = bs % upar(irp_cp)
    state % cv       = bs % upar(irp_cv)
    state % abar     = bs % upar(irp_abar)
    state % zbar     = bs % upar(irp_zbar)
    state % y_e      = bs % upar(irp_ye)
    state % eta      = bs % upar(irp_eta)
    state % cs       = bs % upar(irp_cs)
    state % dx       = bs % upar(irp_dx)

    if (integrate_molar_fraction) then
       state % dhdX(:)  = bs % upar(irp_dhdY:irp_dhdY-1+nspec) * aionInv(:)
       state % dedX(:)  = bs % upar(irp_dedY:irp_dedY-1+nspec) * aionInv(:)
    else
       state % dhdX(:)  = bs % upar(irp_dhdY:irp_dhdY-1+nspec)
       state % dedX(:)  = bs % upar(irp_dedY:irp_dedY-1+nspec)
    endif

    state % T_old    = bs % upar(irp_Told)
    state % dcvdt    = bs % upar(irp_dcvdt)
    state % dcpdt    = bs % upar(irp_dcpdt)

    state % ydot = bs % dydt
    state % ydot(net_itemp) = state % ydot(net_itemp) * temp_scale
    state % ydot(net_ienuc) = state % ydot(net_ienuc) * ener_scale

    !do n = 1, neqs
    !   state % jac(:,n) = bs % jac(:,n)
    !enddo

    !state % jac(net_itemp,:) = state % jac(net_itemp,:) * temp_scale
    !state % jac(net_ienuc,:) = state % jac(net_ienuc,:) * ener_scale

    if (bs % upar(irp_have_rates) > ZERO) then
       state % have_rates = .true.
    else
       state % have_rates = .false.
    endif

    do n = 1, nrates
       state % rates(:,n) = bs % upar(irp_rates+(n-1)*num_rate_groups:irp_rates+n*num_rate_groups-1)
    enddo

    if (bs % upar(irp_self_heat) > ZERO) then
       state % self_heat = .true.
    else
       state % self_heat = .false.
    endif

    state % time = bs % t

    state % i = bs % i
    state % j = bs % j
    state % k = bs % k

  end subroutine bs_to_burn

end module bs_type_module
