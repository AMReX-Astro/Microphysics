module bdf_type_module

  use bl_types, only: dp_t
  use burn_type_module, only: neqs
  use rpar_indices, only: n_rpar_comps

  implicit none

  integer, parameter :: bdf_npt = 1         ! number of points
  integer, parameter :: bdf_max_order = 5   ! maximum order (1 to 6)
  integer :: A(0:bdf_max_order, 0:bdf_max_order) ! pascal matrix, shared by all ts
  !$acc declare create(A)

  !
  ! bdf time-stepper
  !
  type :: bdf_ts

     integer  :: npt
     integer  :: neq
     integer  :: max_order
     integer  :: max_steps                  ! maximum allowable number of steps
     integer  :: max_iters                  ! maximum allowable number of newton iterations
     integer  :: verbose                    ! verbosity level
     real(dp_t) :: dt_min                   ! minimum allowable step-size
     real(dp_t) :: eta_min                  ! minimum allowable step-size shrink factor
     real(dp_t) :: eta_max                  ! maximum allowable step-size growth factor
     real(dp_t) :: eta_thresh               ! step-size growth threshold
     integer  :: max_j_age                  ! maximum age of Jacobian
     integer  :: max_p_age                  ! maximum age of newton iteration matrix

     logical  :: debug
     integer  :: dump_unit

     real(dp_t) :: rtol(neqs)               ! relative tolerances
     real(dp_t) :: atol(neqs)               ! absolute tolerances

     ! state
     real(dp_t) :: t                        ! current time
     real(dp_t) :: t1                       ! final time
     real(dp_t) :: dt                       ! current time step
     real(dp_t) :: dt_nwt                   ! dt used when building newton iteration matrix
     integer  :: k                          ! current order
     integer  :: n                          ! current step
     integer  :: j_age                      ! age of Jacobian
     integer  :: p_age                      ! age of newton iteration matrix
     integer  :: k_age                      ! number of steps taken at current order
     real(dp_t) :: tq(-1:2)                 ! error coefficients (test quality)
     real(dp_t) :: tq2save
     logical  :: refactor

     real(dp_t) :: J(neqs,neqs,bdf_npt)               ! Jacobian matrix
     real(dp_t) :: P(neqs,neqs,bdf_npt)               ! Newton iteration matrix
     real(dp_t) :: z(neqs,bdf_npt,0:bdf_max_order)    ! Nordsieck histroy array, indexed as (dof, p, n)
     real(dp_t) :: z0(neqs,bdf_npt,0:bdf_max_order)   ! Nordsieck predictor array
     real(dp_t) :: h(0:bdf_max_order)                 ! time steps, h = [ h_n, h_{n-1}, ..., h_{n-k} ]
     real(dp_t) :: l(0:bdf_max_order)                 ! predictor/corrector update coefficients
     real(dp_t) :: shift(0:bdf_max_order)             ! scratch array to hold shifted arrays
     real(dp_t) :: upar(n_rpar_comps,bdf_npt)         ! array of user parameters (passed to
                                                      !    user's Jacobian and f)
     real(dp_t) :: y(neqs,bdf_npt)                    ! current y
     real(dp_t) :: yd(neqs,bdf_npt)                   ! current \dot{y}
     real(dp_t) :: rhs(neqs,bdf_npt)                  ! solver rhs
     real(dp_t) :: e(neqs,bdf_npt)                    ! accumulated correction
     real(dp_t) :: e1(neqs,bdf_npt)                   ! accumulated correction, previous step
     real(dp_t) :: ewt(neqs,bdf_npt)                  ! cached error weights
     real(dp_t) :: b(neqs,bdf_npt)                    ! solver work space
     integer    :: ipvt(neqs,bdf_npt)                 ! pivots (neq,npts)

     ! counters
     integer :: nfe                         ! number of function evaluations
     integer :: nje                         ! number of Jacobian evaluations
     integer :: nlu                         ! number of factorizations
     integer :: nit                         ! number of non-linear solver iterations
     integer :: nse                         ! number of non-linear solver errors
     integer :: ncse                        ! number of consecutive non-linear solver errors
     integer :: ncit                        ! number of current non-linear solver iterations
     integer :: ncdtmin                     ! number of consecutive times we tried to shrink beyond the minimum time step

  end type bdf_ts

contains

  subroutine clean_state(state)

    use bl_constants_module, only: ONE
    use actual_network, only: nspec, nspec_evolve, aion
    use integration_data, only: aionInv, temp_scale

    implicit none

    type (bdf_ts) :: state

    ! Ensure that mass fractions always stay positive.

    state % y(1:nspec_evolve,1) = max(state % y(1:nspec_evolve,1) * aion(1:nspec_evolve), 1.d-200) * aionInv(1:nspec_evolve)

  end subroutine clean_state

  subroutine renormalize_species(state)

    use actual_network, only: nspec, nspec_evolve, aion
    use rpar_indices, only: irp_nspec

    implicit none

    type (bdf_ts) :: state

    real(dp_t) :: nspec_sum

    ! Optionally, renormalize them so they sum to unity.

    nspec_sum = sum(state % y(1:nspec_evolve,1) * aion(1:nspec_evolve)) + &
                sum(state % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1,1) * aion(nspec_evolve+1:nspec))

    state % y(1:nspec_evolve,1) = state % y(1:nspec_evolve,1) / nspec_sum
    state % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1,1) = state % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1,1) / nspec_sum

  end subroutine renormalize_species



  subroutine update_thermodynamics(state)

    use bl_constants_module, only: ZERO
    use eos_type_module, only: eos_t, composition
    use eos_module, only: eos_input_rt, eos
    use extern_probin_module, only: call_eos_in_rhs, dT_crit
    use rpar_indices, only: irp_cp, irp_cv, irp_Told, irp_dcpdt, irp_dcvdt, irp_self_heat

    implicit none

    type (bdf_ts) :: state

    type (eos_t) :: eos_state

    ! Several thermodynamic quantities come in via ts % upar -- note: these
    ! are evaluated at the start of the integration, so if things change
    ! dramatically, they will fall out of sync with the current
    ! thermodynamics.

    call vbdf_to_eos(eos_state, state)

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

    if (call_eos_in_rhs .and. state % upar(irp_self_heat,1) > ZERO) then

       call eos(eos_input_rt, eos_state)

    else if (abs(eos_state % T - state % upar(irp_Told,1)) > dT_crit * eos_state % T .and. &
             state % upar(irp_self_heat,1) > ZERO) then

       call eos(eos_input_rt, eos_state)

       state % upar(irp_dcvdt,1) = (eos_state % cv - state % upar(irp_cv,1)) / (eos_state % T - state % upar(irp_Told,1))
       state % upar(irp_dcpdt,1) = (eos_state % cp - state % upar(irp_cp,1)) / (eos_state % T - state % upar(irp_Told,1))
       state % upar(irp_Told,1)  = eos_state % T

    else

       call composition(eos_state)

    endif

    call eos_to_vbdf(eos_state, state)

  end subroutine update_thermodynamics



  ! Given a bdf_ts, set up an EOS state.
  ! We could fill the energy component by storing the initial energy in
  ! ts % upar if we wanted, but we won't do that because (1) if we do call the EOS,
  ! it is always in (rho, T) mode and (2) converting back would imply subtracting
  ! off the nuclear energy from the zone's internal energy, which could lead to
  ! issues from roundoff error if the energy released from burning is small.

  subroutine vbdf_to_eos(state, ts)

    !$acc routine seq

    use actual_network, only: nspec, nspec_evolve, aion
    use integration_data, only: aionInv, dens_scale, temp_scale
    use eos_type_module, only: eos_t
    use burn_type_module, only: net_itemp
    use rpar_indices, only: irp_dens, irp_nspec, irp_cp, irp_cv, irp_abar, irp_zbar, &
                            irp_eta, irp_ye, irp_cs, irp_dhdY, irp_dedY

    implicit none

    type (eos_t)  :: state
    type (bdf_ts) :: ts

    state % rho     = ts % upar(irp_dens,1) * dens_scale
    state % T       = ts % y(net_itemp,1) * temp_scale
    state % xn(1:nspec_evolve) = ts % y(1:nspec_evolve,1) * aion(1:nspec_evolve)
    state % xn(nspec_evolve+1:nspec) = ts % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1,1) * aion(nspec_evolve+1:nspec)
    state % cp      = ts % upar(irp_cp,1)
    state % cv      = ts % upar(irp_cv,1)
    state % abar    = ts % upar(irp_abar,1)
    state % zbar    = ts % upar(irp_zbar,1)
    state % eta     = ts % upar(irp_eta,1)
    state % y_e     = ts % upar(irp_ye,1)
    state % cs      = ts % upar(irp_cs,1)
    state % dhdX(:) = ts % upar(irp_dhdY:irp_dhdY-1+nspec,1) * aionInv(:)
    state % dedX(:) = ts % upar(irp_dedY:irp_dedY-1+nspec,1) * aionInv(:)

  end subroutine vbdf_to_eos



  ! Given an EOS state, fill the rpar and integration state data.

  subroutine eos_to_vbdf(state, ts)

    !$acc routine seq

    use actual_network, only: nspec, nspec_evolve, aion
    use integration_data, only: aionInv, inv_dens_scale, inv_temp_scale
    use eos_type_module, only: eos_t
    use rpar_indices, only: irp_dens, irp_nspec, irp_cp, irp_cv, irp_abar, irp_zbar, &
                            irp_eta, irp_ye, irp_cs, irp_dhdY, irp_dedY
    use integration_data, only: temp_scale, dens_scale
    use burn_type_module, only: net_itemp

    implicit none

    type (eos_t)  :: state
    type (bdf_ts) :: ts

    ts % upar(irp_dens,1)                  = state % rho * inv_dens_scale
    ts % y(net_itemp,1)                    = state % T * inv_temp_scale
    ts % y(1:nspec_evolve,1)               = state % xn(1:nspec_evolve) * aionInv(1:nspec_evolve)
    ts % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1,1) = state % xn(nspec_evolve+1:nspec) * aionInv(nspec_evolve+1:nspec)
    ts % upar(irp_cp,1)                    = state % cp
    ts % upar(irp_cv,1)                    = state % cv
    ts % upar(irp_abar,1)                  = state % abar
    ts % upar(irp_zbar,1)                  = state % zbar
    ts % upar(irp_eta,1)                   = state % eta
    ts % upar(irp_ye,1)                    = state % y_e
    ts % upar(irp_cs,1)                    = state % cs
    ts % upar(irp_dhdY:irp_dhdY+nspec-1,1) = state % dhdX(:) * aion(:)
    ts % upar(irp_dedY:irp_dedY+nspec-1,1) = state % dedX(:) * aion(:)

  end subroutine eos_to_vbdf



  ! Given a burn state, fill the bdf_ts.

  subroutine burn_to_vbdf(state, ts)

    !$acc routine seq

    use actual_network, only: nspec, nspec_evolve, aion, nrates
    use integration_data, only: aionInv, inv_dens_scale, inv_temp_scale, inv_ener_scale
    use rpar_indices, only: irp_dens, irp_nspec, irp_cp, irp_cv, irp_abar, irp_zbar, &
                            irp_ye, irp_eta, irp_cs, irp_dx, irp_dhdY, irp_dedY, &
                            irp_Told, irp_dcvdt, irp_dcpdt, irp_self_heat, irp_have_rates, irp_rates
    use burn_type_module, only: burn_t, net_itemp, net_ienuc, num_rate_groups
    use bl_constants_module, only: ONE

    implicit none

    type (burn_t) :: state
    type (bdf_ts) :: ts

    integer :: n

    ts % upar(irp_dens,1)                           = state % rho * inv_dens_scale
    ts % y(net_itemp,1)                             = state % T * inv_temp_scale
    ts % y(1:nspec_evolve,1)                        = state % xn(1:nspec_evolve) * aionInv(1:nspec_evolve)
    ts % y(net_ienuc,1)                             = state % e * inv_ener_scale
    ts % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1,1) = state % xn(nspec_evolve+1:nspec) * aionInv(nspec_evolve+1:nspec)
    ts % upar(irp_cp,1)                             = state % cp
    ts % upar(irp_cv,1)                             = state % cv
    ts % upar(irp_abar,1)                           = state % abar
    ts % upar(irp_zbar,1)                           = state % zbar
    ts % upar(irp_ye,1)                             = state % y_e
    ts % upar(irp_eta,1)                            = state % eta
    ts % upar(irp_cs,1)                             = state % cs
    ts % upar(irp_dx,1)                             = state % dx
    ts % upar(irp_dhdY:irp_dhdY+nspec-1,1)          = state % dhdX(:) * aion(:)
    ts % upar(irp_dedY:irp_dedY+nspec-1,1)          = state % dedX(:) * aion(:)
    ts % upar(irp_Told,1)                           = state % T_old
    ts % upar(irp_dcvdt,1)                          = state % dcvdt
    ts % upar(irp_dcpdt,1)                          = state % dcpdt

    ts % yd(:,1) = state % ydot
    ts % yd(net_itemp,1) = ts % yd(net_itemp,1) * inv_temp_scale
    ts % yd(net_ienuc,1) = ts % yd(net_ienuc,1) * inv_ener_scale

    ts % J(:,:,1) = state % jac
    ts % J(net_itemp,:,1) = ts % J(net_itemp,:,1) * inv_temp_scale
    ts % J(net_ienuc,:,1) = ts % J(net_ienuc,:,1) * inv_ener_scale

    if (state % have_rates) then
       ts % upar(irp_have_rates,1) = ONE
    else
       ts % upar(irp_have_rates,1) = -ONE
    endif

    do n = 1, nrates
       ts % upar(irp_rates+(n-1)*num_rate_groups:irp_rates+n*num_rate_groups-1,1) = state % rates(:,n)
    enddo

    if (state % self_heat) then
       ts % upar(irp_self_heat,1) = ONE
    else
       ts % upar(irp_self_heat,1) = -ONE
    endif

  end subroutine burn_to_vbdf


  ! Given a bdf_ts, set up a burn state.

  subroutine vbdf_to_burn(ts, state)

    !$acc routine seq

    use actual_network, only: nspec, nspec_evolve, aion, nrates
    use integration_data, only: aionInv, dens_scale, temp_scale, ener_scale
    use rpar_indices, only: irp_dens, irp_nspec, irp_cp, irp_cv, irp_abar, irp_zbar, &
                            irp_ye, irp_eta, irp_cs, irp_dx, irp_dhdY, irp_dedY, &
                            irp_Told, irp_dcvdt, irp_dcpdt, irp_self_heat, irp_have_rates, irp_rates
    use burn_type_module, only: burn_t, net_itemp, net_ienuc, num_rate_groups
    use bl_constants_module, only: ZERO

    implicit none

    type (burn_t) :: state
    type (bdf_ts) :: ts

    integer :: n

    state % rho      = ts % upar(irp_dens,1) * dens_scale
    state % T        = ts % y(net_itemp,1) * temp_scale
    state % e        = ts % y(net_ienuc,1) * ener_scale
    state % xn(1:nspec_evolve) = ts % y(1:nspec_evolve,1) * aion(1:nspec_evolve)
    state % xn(nspec_evolve+1:nspec) = ts % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1,1) * aion(nspec_evolve+1:nspec)
    state % cp       = ts % upar(irp_cp,1)
    state % cv       = ts % upar(irp_cv,1)
    state % abar     = ts % upar(irp_abar,1)
    state % zbar     = ts % upar(irp_zbar,1)
    state % y_e      = ts % upar(irp_ye,1)
    state % eta      = ts % upar(irp_eta,1)
    state % cs       = ts % upar(irp_cs,1)
    state % dx       = ts % upar(irp_dx,1)
    state % dhdX(:)  = ts % upar(irp_dhdY:irp_dhdY-1+nspec,1) * aionInv(:)
    state % dedX(:)  = ts % upar(irp_dedY:irp_dedY-1+nspec,1) * aionInv(:)
    state % T_old    = ts % upar(irp_Told,1)
    state % dcvdt    = ts % upar(irp_dcvdt,1)
    state % dcpdt    = ts % upar(irp_dcpdt,1)

    state % ydot = ts % yd(:,1)
    state % ydot(net_itemp) = state % ydot(net_itemp) * temp_scale
    state % ydot(net_ienuc) = state % ydot(net_ienuc) * ener_scale

    state % jac = ts % J(:,:,1)
    state % jac(net_itemp,:) = state % jac(net_itemp,:) * temp_scale
    state % jac(net_ienuc,:) = state % jac(net_ienuc,:) * ener_scale

    if (ts % upar(irp_have_rates,1) > ZERO) then
       state % have_rates = .true.
    else
       state % have_rates = .false.
    endif

    do n = 1, nrates
       state % rates(:,n) = ts % upar(irp_rates+(n-1)*num_rate_groups:irp_rates+n*num_rate_groups-1,1)
    enddo

    if (ts % upar(irp_self_heat,1) > ZERO) then
       state % self_heat = .true.
    else
       state % self_heat = .false.
    endif

  end subroutine vbdf_to_burn

end module bdf_type_module
