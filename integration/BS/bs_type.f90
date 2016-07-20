module bs_type_module

  use bl_types, only: dp_t
  use burn_type_module, only: neqs, burn_t
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

     real(kind=dp_t) :: y(neqs)
     real(kind=dp_t) :: atol(neqs), rtol(neqs)
     real(kind=dp_t) :: upar(n_rpar_comps)
     real(kind=dp_t) :: t, dt, tmax
     integer         :: n
     integer         :: n_rhs, n_jac

     type(burn_t) :: burn_s

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
    use rpar_indices, only: irp_nspec, n_not_evolved
    use extern_probin_module, only: integrate_molar_fraction

    implicit none

    type (bs_t) :: state

    real(dp_t) :: nspec_sum

    if (integrate_molar_fraction) then
       nspec_sum = &
            sum(state % y(1:nspec_evolve) * aion(1:nspec_evolve)) + &
            sum(state % upar(irp_nspec:irp_nspec+n_not_evolved-1) * aion(nspec_evolve+1:nspec))
    else
       nspec_sum = &
            sum(state % y(1:nspec_evolve)) + &
            sum(state % upar(irp_nspec:irp_nspec+n_not_evolved-1))
    endif

    state % y(1:nspec_evolve) = state % y(1:nspec_evolve) / nspec_sum
    state % upar(irp_nspec:irp_nspec+n_not_evolved-1) = &
         state % upar(irp_nspec:irp_nspec+n_not_evolved-1) / nspec_sum

  end subroutine renormalize_species


  subroutine update_thermodynamics(state)

    !$acc routine seq

    use bl_constants_module, only: ZERO
    use eos_type_module, only: eos_t, composition
    use eos_module, only: eos_input_rt, eos
    use extern_probin_module, only: call_eos_in_rhs, dT_crit

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

    if (call_eos_in_rhs .and. state % burn_s % self_heat) then

       call eos(eos_input_rt, eos_state)
       call eos_to_bs(eos_state, state)

    else if (abs(eos_state % T - state % burn_s % T_old) > &
         dT_crit * eos_state % T .and. state % burn_s % self_heat) then

       call eos(eos_input_rt, eos_state)

       state % burn_s % dcvdt = (eos_state % cv - state % burn_s % cv) / &
            (eos_state % T - state % burn_s % T_old)
       state % burn_s % dcpdt = (eos_state % cp - state % burn_s % cp) / &
            (eos_state % T - state % burn_s % T_old)
       state % burn_s % T_old  = eos_state % T

       ! note: the update to state % upar(irp_cv) and irp_cp is done
       ! in the call to eos_to_bs that follows 
       call eos_to_bs(eos_state, state)

    else

       call composition(eos_state)

       ! just update what is needed here
       state % burn_s % y_e = eos_state % y_e
       state % burn_s % abar = eos_state % abar
       state % burn_s % zbar = eos_state % zbar

    endif



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
    use rpar_indices, only: irp_nspec, n_not_evolved
    use burn_type_module, only: net_itemp
    use extern_probin_module, only: integrate_molar_fraction

    implicit none

    type (eos_t) :: state
    type (bs_t) :: bs

    state % rho     = bs % burn_s % rho * dens_scale
    state % T       = bs % y(net_itemp) * temp_scale

    if (integrate_molar_fraction) then
       state % xn(1:nspec_evolve) = bs % y(1:nspec_evolve) * aion(1:nspec_evolve)
       state % xn(nspec_evolve+1:nspec) = &
            bs % upar(irp_nspec:irp_nspec+n_not_evolved-1) * aion(nspec_evolve+1:nspec)
    else
       state % xn(1:nspec_evolve) = bs % y(1:nspec_evolve)
       state % xn(nspec_evolve+1:nspec) = &
            bs % upar(irp_nspec:irp_nspec+n_not_evolved-1)
    endif

    ! we don't copy any of the other quantities, since we can always
    ! access them through the original bs type

  end subroutine bs_to_eos



  ! Given an EOS state, fill the rpar and integration state data.

  subroutine eos_to_bs(state, bs)

    !$acc routine seq

    use actual_network, only: nspec, nspec_evolve, aion
    use integration_data, only: aionInv, inv_dens_scale, inv_temp_scale
    use eos_type_module, only: eos_t
    use rpar_indices, only: irp_nspec, n_not_evolved
    use integration_data, only: temp_scale, dens_scale
    use burn_type_module, only: net_itemp
    use extern_probin_module, only: integrate_molar_fraction

    implicit none

    type (eos_t) :: state
    type (bs_t) :: bs

    bs % burn_s % rho = state % rho * inv_dens_scale
    bs % burn_s % T = state % T * inv_temp_scale

    if (integrate_molar_fraction) then
       bs % y(1:nspec_evolve) = state % xn(1:nspec_evolve) * aionInv(1:nspec_evolve)
       bs % upar(irp_nspec:irp_nspec+n_not_evolved-1) = &
            state % xn(nspec_evolve+1:nspec) * aionInv(nspec_evolve+1:nspec)
    else
       bs % y(1:nspec_evolve) = state % xn(1:nspec_evolve)
       bs % upar(irp_nspec:irp_nspec+n_not_evolved-1) = &
            state % xn(nspec_evolve+1:nspec) 
    endif

    bs % burn_s % cp = state % cp
    bs % burn_s % cv = state % cv
    bs % burn_s % abar = state % abar
    bs % burn_s % zbar = state % zbar
    bs % burn_s % eta = state % eta
    bs % burn_s % y_e = state % y_e
    bs % burn_s % cs = state % cs

  end subroutine eos_to_bs



  subroutine burn_to_bs(bs)

    ! Given a burn state, fill the bs_t.  For this implementation, we just
    ! modify the burn_t that is part of the bs_t, in place

    !$acc routine seq

    use actual_network, only: nspec, nspec_evolve, aion, nrates
    use integration_data, only: aionInv, inv_dens_scale, inv_temp_scale, inv_ener_scale
    use rpar_indices, only: irp_nspec, n_not_evolved
    use burn_type_module, only: burn_t, net_itemp, net_ienuc, num_rate_groups
    use bl_constants_module, only: ONE
    use extern_probin_module, only: integrate_molar_fraction

    implicit none

    type (burn_t) :: state
    type (bs_t) :: bs

    integer :: n

    bs % burn_s % rho = bs % burn_s % rho * inv_dens_scale
    bs % y(net_itemp) = bs % burn_s % T * inv_temp_scale

    if (integrate_molar_fraction) then
       bs % y(1:nspec_evolve) = bs % burn_s % xn(1:nspec_evolve) * aionInv(1:nspec_evolve)
       bs % upar(irp_nspec:irp_nspec+n_not_evolved-1) = &
            bs % burn_s % xn(nspec_evolve+1:nspec) * aionInv(nspec_evolve+1:nspec)
    else
       bs % y(1:nspec_evolve) = bs % burn_s % xn(1:nspec_evolve)
       bs % upar(irp_nspec:irp_nspec+n_not_evolved-1) = &
            bs % burn_s % xn(nspec_evolve+1:nspec) 
    endif

    bs % y(net_ienuc) = bs % burn_s % e * inv_ener_scale

    ! we don't need to do anything to thermodynamic quantities, cp, cv, ...

    bs % burn_s % ydot(net_itemp) = bs % burn_s % ydot(net_itemp) * inv_temp_scale
    bs % burn_s % ydot(net_ienuc) = bs % burn_s % ydot(net_ienuc) * inv_ener_scale

    bs % burn_s % jac(net_itemp,:) = bs % burn_s % jac(net_itemp,:) * inv_temp_scale
    bs % burn_s % jac(net_ienuc,:) = bs % burn_s % jac(net_ienuc,:) * inv_ener_scale

    ! fixme: net_itemp, net_itemp shouldn't be scaled

    ! we don't need to do anything to have_rates, rates, or self_heat

  end subroutine burn_to_bs



  subroutine bs_to_burn(bs)
    ! Given a bs_t, set up a burn state.  For this implementation, we just
    ! modify the burn_t that is part of the bs_t, in place
    
    !$acc routine seq

    use actual_network, only: nspec, nspec_evolve, aion, nrates
    use integration_data, only: aionInv, dens_scale, temp_scale, ener_scale
    use rpar_indices, only: irp_nspec, n_not_evolved
    use burn_type_module, only: burn_t, net_itemp, net_ienuc, num_rate_groups
    use bl_constants_module, only: ZERO, ONE
    use extern_probin_module, only: integrate_molar_fraction

    implicit none

    type (bs_t) :: bs

    integer :: n

    bs % burn_s % rho = bs % burn_s % rho * dens_scale
    bs % burn_s % T = bs % y(net_itemp) * temp_scale
    bs % burn_s % e = bs % y(net_ienuc) * ener_scale

    if (integrate_molar_fraction) then
       bs % burn_s % xn(1:nspec_evolve) = bs % y(1:nspec_evolve) * aion(1:nspec_evolve)
       bs % burn_s % xn(nspec_evolve+1:nspec) = &
            bs % upar(irp_nspec:irp_nspec+n_not_evolved-1) * aion(nspec_evolve+1:nspec)
    else
       bs % burn_s % xn(1:nspec_evolve) = bs % y(1:nspec_evolve)
       bs % burn_s % xn(nspec_evolve+1:nspec) = &
            bs % upar(irp_nspec:irp_nspec+n_not_evolved-1)
    endif

    ! all the other thermodynamic quantities (cp, cv, ...) are already
    ! in the burn_t

    ! fix the scaling of the derivative terms
    bs % burn_s % ydot(net_itemp) = bs % burn_s % ydot(net_itemp) * temp_scale
    bs % burn_s % ydot(net_ienuc) = bs % burn_s % ydot(net_ienuc) * ener_scale

    bs % burn_s % jac(net_itemp,:) = bs % burn_s % jac(net_itemp,:) * temp_scale
    bs % burn_s % jac(net_ienuc,:) = bs % burn_s % jac(net_ienuc,:) * ener_scale

    ! have_rates and the rates should already be set

    ! self_heat is already set

    bs % burn_s % time = bs % t

  end subroutine bs_to_burn

end module bs_type_module
