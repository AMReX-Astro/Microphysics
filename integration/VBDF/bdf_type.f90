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

    !$acc routine seq

    use bl_constants_module, only: ONE
    use extern_probin_module, only: integrate_molar_fraction
    use network, only: nspec, nspec_evolve, aion, aion_inv
    use burn_type_module, only: net_itemp
    use eos_module, only : eos_get_small_temp

    implicit none

    ! this should be larger than any reasonable temperature we will encounter
    real (kind=dp_t), parameter :: MAX_TEMP = 1.0d11

    ! this is the absolute cutoff for species -- note that this might
    ! be larger than small_x that the user set, but the issue is that
    ! we can have underflow issues if the integrator has to keep track
    ! of species mass fractions much smaller than this.
    real (kind=dp_t), parameter :: SMALL_X_SAFE = 1.0d-200
    real (kind=dp_t) :: small_temp

    type (bdf_ts) :: state

    ! Ensure that mass fractions always stay positive.
    if (integrate_molar_fraction) then
       state % y(1:nspec_evolve,1) = &
            max(state % y(1:nspec_evolve,1) * aion(1:nspec_evolve), &
                SMALL_X_SAFE) * aion_inv(1:nspec_evolve)
    else
       state % y(1:nspec_evolve,1) = max(state % y(1:nspec_evolve,1), SMALL_X_SAFE)
    endif

    ! Ensure that the temperature always stays within reasonable limits.
    call eos_get_small_temp(small_temp)

    state % y(net_itemp,1) = min(MAX_TEMP, max(state % y(net_itemp,1), small_temp))


  end subroutine clean_state


  subroutine renormalize_species(state)

    !$acc routine seq

    use actual_network, only: nspec, nspec_evolve, aion
    use rpar_indices, only: irp_nspec, n_not_evolved
    use extern_probin_module, only: integrate_molar_fraction

    implicit none

    type (bdf_ts) :: state

    real(dp_t) :: nspec_sum

    if (integrate_molar_fraction) then
       nspec_sum = &
            sum(state % y(1:nspec_evolve,1) * aion(1:nspec_evolve)) + &
            sum(state % upar(irp_nspec:irp_nspec+n_not_evolved-1,1) * aion(nspec_evolve+1:nspec))
    else
       nspec_sum = &
            sum(state % y(1:nspec_evolve,1)) + &
            sum(state % upar(irp_nspec:irp_nspec+n_not_evolved-1,1))
    endif

    state % y(1:nspec_evolve,1) = state % y(1:nspec_evolve,1) / nspec_sum
    state % upar(irp_nspec:irp_nspec+n_not_evolved-1,1) = &
         state % upar(irp_nspec:irp_nspec+n_not_evolved-1,1) / nspec_sum

  end subroutine renormalize_species


  subroutine update_thermodynamics(state)

    !$acc routine seq

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

       state % upar(irp_dcvdt,1) = (eos_state % cv - state % upar(irp_cv,1)) / &
            (eos_state % T - state % upar(irp_Told,1))
       state % upar(irp_dcpdt,1) = (eos_state % cp - state % upar(irp_cp,1)) / &
            (eos_state % T - state % upar(irp_Told,1))
       state % upar(irp_Told,1)  = eos_state % T

       ! note: the update to state % upar(irp_cv) and irp_cp is done
       ! in the call to eos_to_bs that follows this block.

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

    use network, only: nspec, nspec_evolve, aion, aion_inv
    use eos_type_module, only: eos_t
    use rpar_indices, only: irp_dens, irp_nspec, irp_cp, irp_cv, irp_abar, irp_zbar, &
                            irp_eta, irp_ye, irp_cs, n_not_evolved
    use burn_type_module, only: net_itemp
    use extern_probin_module, only: integrate_molar_fraction

    implicit none

    type (eos_t)  :: state
    type (bdf_ts) :: ts

    state % rho     = ts % upar(irp_dens,1)
    state % T       = ts % y(net_itemp,1)

    if (integrate_molar_fraction) then
       state % xn(1:nspec_evolve) = ts % y(1:nspec_evolve,1) * aion(1:nspec_evolve)
       state % xn(nspec_evolve+1:nspec) = &
            ts % upar(irp_nspec:irp_nspec+n_not_evolved-1,1) * aion(nspec_evolve+1:nspec)
    else
       state % xn(1:nspec_evolve) = ts % y(1:nspec_evolve,1)
       state % xn(nspec_evolve+1:nspec) = &
            ts % upar(irp_nspec:irp_nspec+n_not_evolved-1,1)
    endif

    state % cp      = ts % upar(irp_cp,1)
    state % cv      = ts % upar(irp_cv,1)
    state % abar    = ts % upar(irp_abar,1)
    state % zbar    = ts % upar(irp_zbar,1)
    state % eta     = ts % upar(irp_eta,1)
    state % y_e     = ts % upar(irp_ye,1)
    state % cs      = ts % upar(irp_cs,1)

  end subroutine vbdf_to_eos



  ! Given an EOS state, fill the rpar and integration state data.

  subroutine eos_to_vbdf(state, ts)

    !$acc routine seq

    use network, only: nspec, nspec_evolve, aion, aion_inv
    use eos_type_module, only: eos_t
    use rpar_indices, only: irp_dens, irp_nspec, irp_cp, irp_cv, irp_abar, irp_zbar, &
                            irp_eta, irp_ye, irp_cs, n_not_evolved
    use burn_type_module, only: net_itemp
    use extern_probin_module, only: integrate_molar_fraction

    implicit none

    type (eos_t)  :: state
    type (bdf_ts) :: ts

    ts % upar(irp_dens,1)                  = state % rho
    ts % y(net_itemp,1)                    = state % T

    if (integrate_molar_fraction) then
       ts % y(1:nspec_evolve,1) = state % xn(1:nspec_evolve) * aion_inv(1:nspec_evolve)
       ts % upar(irp_nspec:irp_nspec+n_not_evolved-1,1) = &
            state % xn(nspec_evolve+1:nspec) * aion_inv(nspec_evolve+1:nspec)
    else
       ts % y(1:nspec_evolve,1) = state % xn(1:nspec_evolve)
       ts % upar(irp_nspec:irp_nspec+n_not_evolved-1,1) = &
            state % xn(nspec_evolve+1:nspec)
    endif

    ts % upar(irp_cp,1)                    = state % cp
    ts % upar(irp_cv,1)                    = state % cv
    ts % upar(irp_abar,1)                  = state % abar
    ts % upar(irp_zbar,1)                  = state % zbar
    ts % upar(irp_eta,1)                   = state % eta
    ts % upar(irp_ye,1)                    = state % y_e
    ts % upar(irp_cs,1)                    = state % cs

  end subroutine eos_to_vbdf



  ! Given a burn state, fill the bdf_ts.

  subroutine burn_to_vbdf(state, ts)

    !$acc routine seq

    use network, only: nspec, nspec_evolve, aion, aion_inv, nrates
    use rpar_indices, only: irp_dens, irp_nspec, irp_cp, irp_cv, irp_abar, irp_zbar, &
                            irp_ye, irp_eta, irp_cs, irp_dx, &
                            irp_Told, irp_dcvdt, irp_dcpdt, irp_self_heat, &
                            n_not_evolved
    use burn_type_module, only: burn_t, net_itemp, net_ienuc
    use bl_constants_module, only: ONE
    use extern_probin_module, only: integrate_molar_fraction

    implicit none

    type (burn_t) :: state
    type (bdf_ts) :: ts

    integer :: n

    ts % upar(irp_dens,1)                           = state % rho
    ts % y(net_itemp,1)                             = state % T

    if (integrate_molar_fraction) then
       ts % y(1:nspec_evolve,1) = state % xn(1:nspec_evolve) * aion_inv(1:nspec_evolve)
       ts % upar(irp_nspec:irp_nspec+n_not_evolved-1,1) = &
         state % xn(nspec_evolve+1:nspec) * aion_inv(nspec_evolve+1:nspec)
    else
       ts % y(1:nspec_evolve,1) = state % xn(1:nspec_evolve)
       ts % upar(irp_nspec:irp_nspec+n_not_evolved-1,1) = &
         state % xn(nspec_evolve+1:nspec)
    endif

    ts % y(net_ienuc,1)                             = state % e
    ts % upar(irp_cp,1)                             = state % cp
    ts % upar(irp_cv,1)                             = state % cv
    ts % upar(irp_abar,1)                           = state % abar
    ts % upar(irp_zbar,1)                           = state % zbar
    ts % upar(irp_ye,1)                             = state % y_e
    ts % upar(irp_eta,1)                            = state % eta
    ts % upar(irp_cs,1)                             = state % cs
    ts % upar(irp_dx,1)                             = state % dx

    ts % upar(irp_Told,1)                           = state % T_old
    ts % upar(irp_dcvdt,1)                          = state % dcvdt
    ts % upar(irp_dcpdt,1)                          = state % dcpdt

    ts % yd(:,1) = state % ydot

    ts % J(:,:,1) = state % jac

    if (state % self_heat) then
       ts % upar(irp_self_heat,1) = ONE
    else
       ts % upar(irp_self_heat,1) = -ONE
    endif

  end subroutine burn_to_vbdf


  ! Given a bdf_ts, set up a burn state.

  subroutine vbdf_to_burn(ts, state)

    !$acc routine seq

    use network, only: nspec, nspec_evolve, aion, aion_inv, nrates
    use rpar_indices, only: irp_dens, irp_nspec, irp_cp, irp_cv, irp_abar, irp_zbar, &
                            irp_ye, irp_eta, irp_cs, irp_dx, &
                            irp_Told, irp_dcvdt, irp_dcpdt, irp_self_heat, &
                            n_not_evolved
    use burn_type_module, only: burn_t, net_itemp, net_ienuc
    use bl_constants_module, only: ZERO
    use extern_probin_module, only: integrate_molar_fraction

    implicit none

    type (burn_t) :: state
    type (bdf_ts) :: ts

    integer :: n

    state % rho      = ts % upar(irp_dens,1)
    state % T        = ts % y(net_itemp,1)
    state % e        = ts % y(net_ienuc,1)

    if (integrate_molar_fraction) then
       state % xn(1:nspec_evolve) = ts % y(1:nspec_evolve,1) * aion(1:nspec_evolve)
       state % xn(nspec_evolve+1:nspec) = &
            ts % upar(irp_nspec:irp_nspec+n_not_evolved-1,1) * aion(nspec_evolve+1:nspec)
    else
       state % xn(1:nspec_evolve) = ts % y(1:nspec_evolve,1)
       state % xn(nspec_evolve+1:nspec) = &
            ts % upar(irp_nspec:irp_nspec+n_not_evolved-1,1)
    endif

    state % cp       = ts % upar(irp_cp,1)
    state % cv       = ts % upar(irp_cv,1)
    state % abar     = ts % upar(irp_abar,1)
    state % zbar     = ts % upar(irp_zbar,1)
    state % y_e      = ts % upar(irp_ye,1)
    state % eta      = ts % upar(irp_eta,1)
    state % cs       = ts % upar(irp_cs,1)
    state % dx       = ts % upar(irp_dx,1)

    state % T_old    = ts % upar(irp_Told,1)
    state % dcvdt    = ts % upar(irp_dcvdt,1)
    state % dcpdt    = ts % upar(irp_dcpdt,1)

    state % ydot = ts % yd(:,1)

    state % jac = ts % J(:,:,1)

    if (ts % upar(irp_self_heat,1) > ZERO) then
       state % self_heat = .true.
    else
       state % self_heat = .false.
    endif

    state % time = ts % t

  end subroutine vbdf_to_burn

end module bdf_type_module
