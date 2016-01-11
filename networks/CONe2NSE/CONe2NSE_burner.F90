module actual_burner_module

  implicit none

contains

  subroutine actual_burner_init()

    use detonation_module, only: detonation_init

    implicit none

    call detonation_init()

  end subroutine actual_burner_init



  subroutine actual_burner(state_in, state_out, dt, time)

    use eos_module
!    use Hydro_interface, ONLY : Hydro_detectShock
    use network
    use detonation_module
    use NSE_data
    use fundamental_constants_module
    use meth_params_module, only: small_temp
    use prob_params_module, only: dx_level
    use amrinfo_module, only: amr_level
    use castro_util_module, only: position, volume

    implicit none

    type (eos_t) :: state_in, state_out
    double precision :: dt, time
    
    integer                    :: i, j, k, iref, l, n
    integer                    :: istat
    double precision           :: dx

    logical :: shock

    logical :: ignite_detonation

    double precision :: loc(3), dvol, dist, radius

    double precision :: dens, temp, eint, pres
    double precision :: xc12init, xne22init
    double precision :: flame, flamedot
    double precision :: phi_fa, phi_aq, phi_qn
    double precision :: ye, dyi_qn, dqbar_qn
    double precision :: qdot, edotnu
    double precision :: phi_fa_det

    double precision, parameter :: inflame_threshold = 1.e-6
    double precision :: cgsMeVperGram = 9.6485e17

    double precision :: dti

    ! for saving some initial states
    double precision :: phi_fa_i, phi_aq_i, phi_qn_i
    double precision :: eint_i, pres_i
    double precision :: qbar_i, ye_i, yi_i, dye_n

    double precision :: ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a

    ! for calculating time evolution of properties
    double precision :: phi_aq_dot, phi_qn_dot
    double precision :: r_c
    double precision :: yi_q_est
    double precision :: tau_nsqe_inv, tau_nse_inv, maxExp, testExp  
    double precision :: expq

    double precision :: qbar

    ! for predictions of nse final state
    double precision :: s
    double precision :: qbar_finnse_d, sumyi_finnse_d, temp_finnse_d, edot_d, yedot_d
    double precision :: qbar_finnse_p, sumyi_finnse_p, temp_finnse_p, edot_p, yedot_p
    double precision :: qbar_finnse, sumyi_finnse, temp_finnse, edot, yedot

    ! for calculation of reaction rate
    double precision :: nasigmav, t9, t9a

    ! for estimating reaction temperatures
    double precision :: qbar_unburned, yi_unburned, qbar_burned, yi_burned, x_bnf, qbar_noflame
    double precision :: yi_noflame, ye_noflame, eint_noflame

    type (eos_t) :: eos_state

    logical, parameter :: useShockBurn = .false.

100 format("*** igniting at ",f7.5,4(1x,es10.3)," ***")

    dx = dx_level(1,amr_level) ! also assume square grid

    ! shock detect if burning is turned off in shocks
    if (thermalReact .and. (.NOT. useShockBurn)) then
       ! call Hydro_detectShock(solnData, shock, blkLimits, blkLimitsGC, (/0,0,0/), &
       !      xCoord,yCoord,zCoord)
    else
       shock = .false.
    endif

    ! ! --------------------------------
    ! ! check if we are near a detonation point
    ! ! --------------------------------

    ignite_detonation = .false.
    if (autoDDT) radius = IgnRad

    loc = position(state_in % idx(1), state_in % idx(2), state_in % idx(3))

    do l = 1, IgnNum

       dist = sqrt( sum( (IgnLoc(l,:) - loc)**2 ) )

       if (.not. autoDDT) radius = IgnR(l)
       if ( IgnTime(l) < time + dt .and. IgnTime(l) >= time &
            .and. dist <= radius ) then
          write (6,100) IgnTime(l), IgnLoc(l,1), IgnLoc(l,2), IgnLoc(l,3), radius
          ignite_detonation = .true.
          phi_fa_det = ONE
       endif
    enddo

    ! ! --------------------------------
    ! ! evolve progress variables and update NSE grid quantities
    ! ! --------------------------------

    if (shock) then
       qdot = ZERO
       edotnu = ZERO
       return
    endif

    !--------------------------------------------------------
    ! initialize some information
    ! including getting unburned properties and saving some initial data
    !--------------------------------------------------------

    ! inverse timestep
    dti    = ONE / dt

    dens = state_in % rho
    temp = state_in % T
    eint = state_in % e
    pres = state_in % p

    flame     = state_in % aux(UFLAM)
    flamedot  = state_in % aux(UFLDT)
    
    xc12init  = state_in % xn(iC12)
    xne22init = state_in % xn(iNe22)

    ye        = state_in % aux(UYE)

    phi_fa    = state_in % aux(UPHFA)
    phi_aq    = state_in % aux(UPHAQ)
    phi_qn    = state_in % aux(UPHQN)
    dyi_qn    = state_in % aux(UDYQN)
    dqbar_qn  = state_in % aux(UDQQN)

    call paraFuelAshProperties(xc12init, xne22init, ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a)

    ! truncate advection errors which might have pushed out of limits
    phi_fa = max(ZERO,min(ONE,phi_fa))
    phi_aq = max(ZERO,min(phi_fa,phi_aq))
    phi_qn = max(ZERO,min(phi_aq,phi_qn))

    dyi_qn = max(ZERO,min(ONE,dyi_qn))

    ! initial plasma properties
    qbar_i  = (ONE-phi_fa)*qbar_f + (phi_fa-phi_aq)*qbar_a + dqbar_qn
    ye_i  = ye
    yi_i  = (ONE-phi_fa)*yi_f + (phi_fa-phi_aq)*yi_a + dyi_qn

    dye_n = max( ye - (ONE-phi_qn)*ye_a, ZERO )

    pres_i   = pres
    eint_i   = eint
    phi_fa_i = phi_fa
    phi_aq_i = phi_aq
    phi_qn_i = phi_qn

    !--------------------------------------------------------------------------
    ! Calculate properties (qbar, yi, temp) of NSE final state
    !--------------------------------------------------------------------------
    ! We predict the properties of the final NSE state from the current 
    ! (partially burned) state, assuming the burn will occur at either
    ! constant density (isochoric) or constant pressure (isobaric).
    ! The properties have been pre-calculated and tabulated for speed.  Inside a
    ! flame we use an interpolation between the isobaric and isochoric predictions
    ! In addtion to qbar, yi, and temp of the final state, we also obtain
    ! edot (neutrino losses) and yedot (neutronization rate) for that state.
    ! These are used in the evolution of partially burned zones.
    !--------------------------------------------------------------------------
    if ( flame >= 0.9999 .or. flame < inflame_threshold ) then
       ! isochoric prediction
       call NSE_finalAtDens(qbar_finnse, sumyi_finnse, temp_finnse, edot, yedot, &
            ye_i, dens, eint_i - qbar_i*cgsMeVperGram)

    else if (flame >= 0.99) then
       ! crossover region, linearly average
       call NSE_finalAtDens(qbar_finnse_d, sumyi_finnse_d, temp_finnse_d, edot_d, yedot_d, &
            ye_i, dens, eint_i - qbar_i*cgsMeVperGram)
       call NSE_finalAtPres(qbar_finnse_p, sumyi_finnse_p, temp_finnse_p, edot_p, yedot_p, &
            ye_i, pres_i, eint_i + pres_i/dens - qbar_i*cgsMeVperGram)
       s = (flame - 0.99)/(0.9999-0.99)
       qbar_finnse = s*qbar_finnse_d + (ONE-s)*qbar_finnse_p
       sumyi_finnse = ONE/(s/sumyi_finnse_d + (ONE-s)/sumyi_finnse_p)
       temp_finnse = s*temp_finnse_d + (ONE-s)*temp_finnse_p
       edot = s*edot_d + (ONE-s)*edot_p
       yedot = s*yedot_d + (ONE-s)*yedot_p

    else
       ! isobaric prediction
       call NSE_finalAtPres(qbar_finnse, sumyi_finnse, temp_finnse, edot, yedot, &
            ye_i, pres_i, eint_i + pres_i/dens - qbar_i*cgsMeVperGram)
    endif

    !-------------------------------------------------------------------
    !   calculate timescale calibrated in Calder etal 2007
    !-------------------------------------------------------------------

    testExp = 182.06e9/temp_finnse - 46.054d0

    ! protect from overflow for portability
    maxExp = log(HUGE(ONE))
    if (testExp .GE. maxExp) then
       tau_nsqe_inv = ZERO
    else
       tau_nsqe_inv = exp(-testExp)
    endif

    if ( flame >= 0.9999 .or. &
         ( thermalReact .and. &
         ( flame < inflame_threshold .or. (phi_fa > flame + 0.1) ) ) ) then
       ! use local temperature to get NSE timescale outside flame
       ! (particularly in detonations)

       testExp = 201.0d9/temp - 46.77d0
    else
       ! in flame with no apparent thermal burning, use estimated final temperature
       ! this is just a rate, so a discontinuous crossover shouldn't be a big deal

       testExp = 201.0d9/temp_finnse - 46.77d0
    endif
    if (testExp .GE. maxExp) then
       tau_nse_inv = ZERO
    else
       tau_nse_inv = exp(-testExp)
    endif


    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------

    ! if we are in the flame, only include thermal reactions if there is nearby
    ! thermal burning.  This prevents flame burning from causing thermal reactions
    ! without outside influence (such as an arriving detonation).
    if (thermalReact .and. &
         .not. ( flame >= inflame_threshold .and. state_in % react_proximity > ONE ) ) then
       ! away from flame or in flame near other reactions, use thermal reactions

       ! first need a temperature
       if (flame < inflame_threshold .or. phi_fa > (ONE-inflame_threshold) ) then
          ! away from flame or with nearly fully burned material ...
          ! burn based on nuclear reaction rate determined by t9
          t9 = temp/1.e9
       else
          ! in flame but not yet fully burned...
          ! estimate temperature of material if flame had not burned

          ! catch up corner case in which the flame just started burning
          !   (prevents a divide by zero risk below)
          phi_fa = max(phi_fa,inflame_threshold)
          ! x_bnf = fraction of material which was burned, but not by the flame
          x_bnf = max(ZERO, phi_fa-flame)
          qbar_unburned =  qbar_f
          yi_unburned =  yi_f
          qbar_burned = ( qbar_i - qbar_unburned*(ONE-phi_fa) ) / phi_fa
          yi_burned = ( yi_i - yi_unburned*(ONE-phi_fa) ) / phi_fa
          qbar_noflame = qbar_unburned*(ONE-x_bnf) + qbar_burned*x_bnf
          yi_noflame = yi_unburned*(ONE-x_bnf) + yi_burned*x_bnf
          !ye_noflame   = 0.5*(1-x_bnf) + ye_oldash_i*x_bnf
          ! just use current ye instead
          !ye_noflame = ye_i
          ye_noflame = ye_f
          eint_noflame = eint_i - cgsMeVperGram*(qbar_i-qbar_noflame)

          if (yi_noflame < ZERO) then
             write (6,*) ' neg sumy', phi_fa, flame, dqbar_qn
          endif
          ! find a floor energy for local material
          eos_state % rho  = dens
          eos_state % T    = small_temp
          eos_state % abar = ONE / yi_noflame
          eos_state % zbar = ye_noflame * eos_state % abar
          eos_state % y_e  = ye_noflame
          call eos(eos_input_rt, eos_state)
          if ( eint_noflame < eos_state % e ) then
             t9 = small_temp/1.d9
          else
             eos_state % e = eint_noflame
             eos_state % T = temp ! guess value
             call eos(eos_input_re, eos_state)
             t9 = eos_state % T / 1.d9
          endif
       endif

       ! thermally burn based on temperature just estimated
       ! From Caughlin & Fowler 1988, ADNDT 40 283
       ! note this is UNSCREENED, should be screened at some point
       t9a = t9/(ONE+0.0396*t9)
       nasigmav = 4.27e26*t9a**(5.0/6)*t9**(-1.5)*exp(-84.165*t9a**(-THIRD)-2.12e-3*t9**3)
       ! integrate
       !     phi_fa_dot = dens*xc12init*(1.0-phi_fa)**2/12.0*nasigmav
       ! and add as evolution after flame
       r_c = dens*xc12init/12.0*nasigmav
       phi_fa = ONE - (ONE-phi_fa)/(ONE + r_c*dt*(ONE-phi_fa))

       ! ignition override
       if (ignite_detonation) then
          ! fully burn this zone in this timestep (will hit limiter)
          phi_fa = phi_fa_det
       endif

    endif

    ! follow flame variable in any case
    ! (second operator split step if thermal burning is activated)
    phi_fa = max(ZERO,min(ONE,phi_fa+max(ZERO,flamedot)*dt))


    !------------------------------------------------------
    ! 3 Update progress variables by direct integration, assuming
    !   phi_fa is constant (i.e. evolves much faster than others)

    expq = exp(-dt*tau_nsqe_inv)

    ! integration of     phi_aq_dot = (phi_fa-phi_aq)/tau_nsqe
    phi_aq = phi_fa - (phi_fa-phi_aq_i)*expq
    phi_aq = max(ZERO,min(phi_fa,phi_aq))
    phi_aq_dot = (phi_aq-phi_aq_i)*dti

    ! integration of     phi_qn_dot = (phi_aq-phi_qn)^2/tau_nsqe
    ! assuming phi_aq is constant (i.e. operator split)
    ! need to handle case when phi_aq-phi_qn_i is zero
    if ((phi_aq-phi_qn_i) == ZERO) then
       phi_qn = phi_aq
    else
       phi_qn = phi_aq - ONE/( ONE/(phi_aq-phi_qn_i) + tau_nse_inv*dt )
    endif
    phi_qn = max(ZERO,min(phi_aq,phi_qn))
    phi_qn_dot = (phi_qn-phi_qn_i)*dti

    ! now update partial properties
    ! operator split to mirror later treatment of dyi_qn
    dqbar_qn = dqbar_qn + phi_aq_dot*dt*236.537/28
    dqbar_qn = (phi_aq*236.537/28+phi_qn*(qbar_finnse-236.537/28))*(ONE-expq) + dqbar_qn*expq

    yi_q_est = yiion(iSi28)

    ! now update dyi_qn
    ! assuming phi_qn is constant (operator split again)
    ! two steps, first change due to change of phi_aq, then evolution
    !  better considering operator splitting of phi_qn and phi_aq evolution
    dyi_qn = dyi_qn + phi_aq_dot*dt*yi_q_est
    dyi_qn = (phi_aq*yi_q_est+phi_qn*(sumyi_finnse-yi_q_est))*(ONE-expq) + dyi_qn*expq

    ! implicit update not necessary, but use actual change in phi_qn from above
    !   ( this is like operator splitting it)
    dye_n = dye_n + ye_a*phi_qn_dot*dt + phi_qn*yedot*dt

    ! Construct final state
    qbar  = (ONE-phi_fa)*qbar_f + (phi_fa-phi_aq)*qbar_a + dqbar_qn
    ye    = (ONE-phi_qn)*ye_a + dye_n

    ! calculate energy release including neutrino losses
    ! and rest mass differences
    qdot = cgsMeVperGram*(qbar - qbar_i)*dti  &
         - phi_qn*( yedot*N_A*c_light*c_light*(m_p + m_e - m_n) + edot )

    edotnu = phi_qn*edot



    ! update energy

    state_out % e = state_in % e + qdot * dt

    dvol = volume(state_in % idx(1), state_in % idx(2), state_in % idx(3))

    neutLossThisProcStep = neutLossThisProcStep + state_out % rho * dvol * edotnu * dt

  end subroutine actual_burner  

end module actual_burner_module
