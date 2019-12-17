module reaclib_rates
  use screening_module, only: screen5, add_screening_factor, screening_init, plasma_state, fill_plasma_state
  use network

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  logical, parameter :: screen_reaclib = .true.
  
  ! Temperature coefficient arrays (numbers correspond to reaction numbers in net_info)
  real(rt)        , allocatable :: ctemp_rate(:,:)

  ! Index into ctemp_rate, dimension 2, where each rate's coefficients start
  integer, allocatable :: rate_start_idx(:)
  
  ! Reaction multiplicities-1 (how many rates contribute - 1)
  integer, allocatable :: rate_extra_mult(:)

  !$acc declare create(ctemp_rate, rate_start_idx, rate_extra_mult)
  !$acc declare copyin(screen_reaclib)

contains

  subroutine init_reaclib()
    
    allocate( ctemp_rate(7, 48) )
    ! n13__c13__weak__wc12
    ctemp_rate(:, 1) = [  &
        -6.76010000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00 ]

    ! o14__n14__weak__wc12
    ctemp_rate(:, 2) = [  &
        -4.62354000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00 ]

    ! o15__n15__weak__wc12
    ctemp_rate(:, 3) = [  &
        -5.17053000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00 ]

    ! f17__o17__weak__wc12
    ctemp_rate(:, 4) = [  &
        -4.53318000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00 ]

    ! p_c12__n13
    ctemp_rate(:, 5) = [  &
        1.75428000000000d+01, &
        -3.77849000000000d+00, &
        -5.10735000000000d+00, &
        -2.24111000000000d+00, &
        1.48883000000000d-01, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 6) = [  &
        1.71482000000000d+01, &
        0.00000000000000d+00, &
        -1.36920000000000d+01, &
        -2.30881000000000d-01, &
        4.44362000000000d+00, &
        -3.15898000000000d+00, &
        -6.66667000000000d-01 ]

    ! he4_c12__o16
    ctemp_rate(:, 7) = [  &
        6.96526000000000d+01, &
        -1.39254000000000d+00, &
        5.89128000000000d+01, &
        -1.48273000000000d+02, &
        9.08324000000000d+00, &
        -5.41041000000000d-01, &
        7.03554000000000d+01 ]

    ctemp_rate(:, 8) = [  &
        2.54634000000000d+02, &
        -1.84097000000000d+00, &
        1.03411000000000d+02, &
        -4.20567000000000d+02, &
        6.40874000000000d+01, &
        -1.24624000000000d+01, &
        1.37303000000000d+02 ]

    ! p_c13__n14
    ctemp_rate(:, 9) = [  &
        1.85155000000000d+01, &
        0.00000000000000d+00, &
        -1.37200000000000d+01, &
        -4.50018000000000d-01, &
        3.70823000000000d+00, &
        -1.70545000000000d+00, &
        -6.66667000000000d-01 ]

    ctemp_rate(:, 10) = [  &
        1.39637000000000d+01, &
        -5.78147000000000d+00, &
        0.00000000000000d+00, &
        -1.96703000000000d-01, &
        1.42126000000000d-01, &
        -2.38912000000000d-02, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 11) = [  &
        1.51825000000000d+01, &
        -1.35543000000000d+01, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 ]

    ! p_n13__o14
    ctemp_rate(:, 12) = [  &
        1.81356000000000d+01, &
        0.00000000000000d+00, &
        -1.51676000000000d+01, &
        9.55166000000000d-02, &
        3.06590000000000d+00, &
        -5.07339000000000d-01, &
        -6.66667000000000d-01 ]

    ctemp_rate(:, 13) = [  &
        1.09971000000000d+01, &
        -6.12602000000000d+00, &
        1.57122000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 ]

    ! p_n14__o15
    ctemp_rate(:, 14) = [  &
        6.73578000000000d+00, &
        -4.89100000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        6.82000000000000d-02 ]

    ctemp_rate(:, 15) = [  &
        7.65444000000000d+00, &
        -2.99800000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 16) = [  &
        2.01169000000000d+01, &
        0.00000000000000d+00, &
        -1.51930000000000d+01, &
        -4.63975000000000d+00, &
        9.73458000000000d+00, &
        -9.55051000000000d+00, &
        3.33333000000000d-01 ]

    ctemp_rate(:, 17) = [  &
        1.70100000000000d+01, &
        0.00000000000000d+00, &
        -1.51930000000000d+01, &
        -1.61954000000000d-01, &
        -7.52123000000000d+00, &
        -9.87565000000000d-01, &
        -6.66667000000000d-01 ]

    ! he4_n14__f18
    ctemp_rate(:, 18) = [  &
        1.38995000000000d+01, &
        -1.09656000000000d+01, &
        -5.62270000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 19) = [  &
        1.96838000000000d-01, &
        -5.16034000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 20) = [  &
        2.15339000000000d+01, &
        0.00000000000000d+00, &
        -3.62504000000000d+01, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -5.00000000000000d+00, &
        -6.66667000000000d-01 ]

    ! p_n15__o16
    ctemp_rate(:, 21) = [  &
        1.45444000000000d+01, &
        -1.02295000000000d+01, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        4.59037000000000d-02, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 22) = [  &
        6.59056000000000d+00, &
        -2.92315000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 23) = [  &
        2.00176000000000d+01, &
        0.00000000000000d+00, &
        -1.52400000000000d+01, &
        3.34926000000000d-01, &
        4.59088000000000d+00, &
        -4.78468000000000d+00, &
        -6.66667000000000d-01 ]

    ! p_o16__f17
    ctemp_rate(:, 24) = [  &
        1.90904000000000d+01, &
        0.00000000000000d+00, &
        -1.66960000000000d+01, &
        -1.16252000000000d+00, &
        2.67703000000000d-01, &
        -3.38411000000000d-02, &
        -6.66667000000000d-01 ]

    ! p_o17__f18
    ctemp_rate(:, 25) = [  &
        1.58929000000000d+01, &
        0.00000000000000d+00, &
        -1.64035000000000d+01, &
        4.31885000000000d+00, &
        -7.09921000000000d-01, &
        -2.00000000000000d+00, &
        -6.66667000000000d-01 ]

    ctemp_rate(:, 26) = [  &
        9.39048000000000d+00, &
        -6.22828000000000d+00, &
        0.00000000000000d+00, &
        2.31435000000000d+00, &
        -3.02835000000000d-01, &
        2.01330000000000d-02, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 27) = [  &
        -1.30770000000000d+01, &
        -7.46296000000000d-01, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 ]

    ! he4_n13__p_o16
    ctemp_rate(:, 28) = [  &
        4.04644000000000d+01, &
        0.00000000000000d+00, &
        -3.58290000000000d+01, &
        -5.30275000000000d-01, &
        -9.82462000000000d-01, &
        8.08059000000000d-02, &
        -6.66667000000000d-01 ]

    ! p_n15__he4_c12
    ctemp_rate(:, 29) = [  &
        2.08972000000000d+01, &
        -7.40600000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 30) = [  &
        -4.87347000000000d+00, &
        -2.02117000000000d+00, &
        0.00000000000000d+00, &
        3.08497000000000d+01, &
        -8.50433000000000d+00, &
        -1.54426000000000d+00, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 31) = [  &
        2.74764000000000d+01, &
        0.00000000000000d+00, &
        -1.52530000000000d+01, &
        1.59318000000000d+00, &
        2.44790000000000d+00, &
        -2.19708000000000d+00, &
        -6.66667000000000d-01 ]

    ctemp_rate(:, 32) = [  &
        -6.57522000000000d+00, &
        -1.16380000000000d+00, &
        0.00000000000000d+00, &
        2.27105000000000d+01, &
        -2.90707000000000d+00, &
        2.05754000000000d-01, &
        -1.50000000000000d+00 ]

    ! he4_o14__p_f17
    ctemp_rate(:, 33) = [  &
        4.08358000000000d+01, &
        0.00000000000000d+00, &
        -3.93880000000000d+01, &
        -1.74673000000000d+01, &
        3.53029000000000d+01, &
        -2.48162000000000d+01, &
        -6.66667000000000d-01 ]

    ctemp_rate(:, 34) = [  &
        1.63087000000000d+01, &
        -2.25100000000000d+01, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 35) = [  &
        1.11184000000000d+01, &
        -1.36000000000000d+01, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 36) = [  &
        -1.06091000000000d+02, &
        -4.53036000000000d-01, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 37) = [  &
        1.21289000000000d+01, &
        -1.20223000000000d+01, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 38) = [  &
        1.86518000000000d+01, &
        -2.60000000000000d+01, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 ]

    ! p_o17__he4_n14
    ctemp_rate(:, 39) = [  &
        1.01740000000000d+01, &
        -4.95865000000000d+00, &
        0.00000000000000d+00, &
        5.10182000000000d+00, &
        3.79373000000000d-01, &
        -6.72515000000000d-02, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 40) = [  &
        5.53360000000000d+00, &
        -2.11477000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 41) = [  &
        -7.20763000000000d+00, &
        -7.53395000000000d-01, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 42) = [  &
        1.95790000000000d+01, &
        0.00000000000000d+00, &
        -1.69078000000000d+01, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -2.00000000000000d+00, &
        -6.66667000000000d-01 ]

    ! p_f18__he4_o15
    ctemp_rate(:, 43) = [  &
        -3.17388000000000d+01, &
        -3.76432000000000d-01, &
        0.00000000000000d+00, &
        6.17380000000000d+01, &
        -1.08290000000000d+02, &
        -3.42365000000000d+01, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 44) = [  &
        6.20058000000000d+01, &
        0.00000000000000d+00, &
        -2.14023000000000d+01, &
        -8.08891000000000d+01, &
        1.34600000000000d+02, &
        -1.26504000000000d+02, &
        -6.66667000000000d-01 ]

    ctemp_rate(:, 45) = [  &
        1.75704000000000d+00, &
        -3.01675000000000d+00, &
        0.00000000000000d+00, &
        1.33223000000000d+01, &
        -1.36696000000000d+00, &
        7.57363000000000d-02, &
        -1.50000000000000d+00 ]

    ! he4_he4_he4__c12
    ctemp_rate(:, 46) = [  &
        -2.43505000000000d+01, &
        -4.12656000000000d+00, &
        -1.34900000000000d+01, &
        2.14259000000000d+01, &
        -1.34769000000000d+00, &
        8.79816000000000d-02, &
        -1.31653000000000d+01 ]

    ctemp_rate(:, 47) = [  &
        -1.17884000000000d+01, &
        -1.02446000000000d+00, &
        -2.35700000000000d+01, &
        2.04886000000000d+01, &
        -1.29882000000000d+01, &
        -2.00000000000000d+01, &
        -2.16667000000000d+00 ]

    ctemp_rate(:, 48) = [  &
        -9.71052000000000d-01, &
        0.00000000000000d+00, &
        -3.70600000000000d+01, &
        2.93493000000000d+01, &
        -1.15507000000000d+02, &
        -1.00000000000000d+01, &
        -1.33333000000000d+00 ]



    allocate( rate_start_idx(nrat_reaclib) )
    rate_start_idx(:) = [ &
      1, &
      2, &
      3, &
      4, &
      5, &
      7, &
      9, &
      12, &
      14, &
      18, &
      21, &
      24, &
      25, &
      28, &
      29, &
      33, &
      39, &
      43, &
      46 ]

    allocate( rate_extra_mult(nrat_reaclib) )
    rate_extra_mult(:) = [ &
      0, &
      0, &
      0, &
      0, &
      1, &
      1, &
      2, &
      1, &
      3, &
      2, &
      2, &
      0, &
      2, &
      0, &
      3, &
      5, &
      3, &
      2, &
      2 ]

    !$acc update device(ctemp_rate, rate_start_idx, rate_extra_mult)
    
  end subroutine init_reaclib

  subroutine term_reaclib()
    deallocate( ctemp_rate )
    deallocate( rate_start_idx )
    deallocate( rate_extra_mult )
  end subroutine term_reaclib

  subroutine net_screening_init()
    ! Adds screening factors and calls screening_init

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jc12), aion(jc12))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jc12), aion(jc12))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jc13), aion(jc13))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jn13), aion(jn13))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jn14), aion(jn14))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jn14), aion(jn14))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jn15), aion(jn15))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jo16), aion(jo16))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jo17), aion(jo17))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jn13), aion(jn13))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jo14), aion(jo14))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jf18), aion(jf18))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jhe4), aion(jhe4))


    call screening_init()    
  end subroutine net_screening_init

  subroutine reaclib_evaluate(pstate, temp, iwhich, reactvec)
    !$acc routine seq

    implicit none
    
    type(plasma_state), intent(in) :: pstate
    real(rt)        , intent(in) :: temp
    integer, intent(in) :: iwhich

    real(rt)        , intent(inout) :: reactvec(num_rate_groups+2)
    ! reactvec(1) = rate     , the reaction rate
    ! reactvec(2) = drate_dt , the Temperature derivative of rate
    ! reactvec(3) = scor     , the screening factor
    ! reactvec(4) = dscor_dt , the Temperature derivative of scor
    ! reactvec(5) = dqweak   , the weak reaction dq-value (ergs)
    !                          (This accounts for modification of the reaction Q
    !                           due to the local density and temperature of the plasma.
    !                           For Reaclib rates, this is 0.0e0_rt.)
    ! reactvec(6) = epart    , the particle energy generation rate (ergs/s)
    ! NOTE: The particle energy generation rate (returned in ergs/s)
    !       is the contribution to enuc from non-ion particles associated
    !       with the reaction.
    !       For example, this accounts for neutrino energy losses
    !       in weak reactions and/or gamma heating of the plasma
    !       from nuclear transitions in daughter nuclei.

    real(rt)          :: rate, scor ! Rate and Screening Factor
    real(rt)          :: drate_dt, dscor_dt ! Temperature derivatives
    real(rt)         :: dscor_dd
    real(rt)         :: ri, T9, T9_exp, lnirate, irate, dirate_dt, dlnirate_dt
    integer :: i, j, m, istart

    ri = 0.0e0_rt
    rate = 0.0e0_rt
    drate_dt = 0.0e0_rt
    irate = 0.0e0_rt
    dirate_dt = 0.0e0_rt
    T9 = temp/1.0e9_rt
    T9_exp = 0.0e0_rt

    ! Use reaction multiplicities to tell whether the rate is Reaclib
    m = rate_extra_mult(iwhich)

    istart = rate_start_idx(iwhich)

    do i = 0, m
       lnirate = ctemp_rate(1, istart+i) + ctemp_rate(7, istart+i) * LOG(T9)
       dlnirate_dt = ctemp_rate(7, istart+i)/T9
       do j = 2, 6
          T9_exp = (2.0e0_rt*dble(j-1)-5.0e0_rt)/3.0e0_rt 
          lnirate = lnirate + ctemp_rate(j, istart+i) * T9**T9_exp
          dlnirate_dt = dlnirate_dt + &
               T9_exp * ctemp_rate(j, istart+i) * T9**(T9_exp-1.0e0_rt)
       end do
       ! If the rate will be in the approx. interval [0.0, 1.0E-100], replace by 0.0
       ! This avoids issues with passing very large negative values to EXP
       ! and getting results between 0.0 and 1.0E-308, the limit for IEEE 754.
       ! And avoids SIGFPE in CVODE due to tiny rates.
       lnirate = max(lnirate, -230.0e0_rt)
       irate = EXP(lnirate)
       rate = rate + irate
       dirate_dt = irate * dlnirate_dt/1.0e9_rt
       drate_dt = drate_dt + dirate_dt
    end do

    reactvec(i_rate)     = rate
    reactvec(i_drate_dt) = drate_dt
    reactvec(i_scor)     = 0.0e0_rt
    reactvec(i_dscor_dt) = 0.0e0_rt
    reactvec(i_dqweak)   = 0.0e0_rt
    reactvec(i_epart)    = 0.0e0_rt

    ! write(*,*) '----------------------------------------'
    ! write(*,*) 'IWHICH: ', iwhich
    ! write(*,*) 'reactvec(i_rate)', reactvec(i_rate)
    ! write(*,*) 'reactvec(i_drate_dt)', reactvec(i_drate_dt)
    ! write(*,*) 'reactvec(i_scor)', reactvec(i_scor)    
    ! write(*,*) 'reactvec(i_dscor_dt)', reactvec(i_dscor_dt)
    ! write(*,*) 'reactvec(i_dqweak)', reactvec(i_dqweak)
    ! write(*,*) 'reactvec(i_epart)', reactvec(i_epart)
    ! write(*,*) '----------------------------------------'

  end subroutine reaclib_evaluate
  
end module reaclib_rates
