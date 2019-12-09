module reaclib_rates
  use screening_module, only: screen5, add_screening_factor, screening_init, plasma_state, fill_plasma_state
  use network
  use microphysics_type_module, only: rt

  implicit none

  logical, parameter :: screen_reaclib = .true.
  
  ! Temperature coefficient arrays (numbers correspond to reaction numbers in net_info)
  real(rt), allocatable :: ctemp_rate(:,:)

  ! Index into ctemp_rate, dimension 2, where each rate's coefficients start
  integer, allocatable :: rate_start_idx(:)
  
  ! Reaction multiplicities-1 (how many rates contribute - 1)
  integer, allocatable :: rate_extra_mult(:)

  ! Should these reactions be screened?
  logical, allocatable :: do_screening(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: ctemp_rate, rate_start_idx, rate_extra_mult, do_screening
#endif

  !$acc declare create(ctemp_rate, rate_start_idx, rate_extra_mult, do_screening)
  !$acc declare copyin(screen_reaclib)

contains

  subroutine init_reaclib()
    
    allocate( ctemp_rate(7, 18) )
    ! he4_he4_he4__c12
    ctemp_rate(:, 1) = [  &
        -2.43505000000000e+01_rt, &
        -4.12656000000000e+00_rt, &
        -1.34900000000000e+01_rt, &
        2.14259000000000e+01_rt, &
        -1.34769000000000e+00_rt, &
        8.79816000000000e-02_rt, &
        -1.31653000000000e+01_rt ]

    ctemp_rate(:, 2) = [  &
        -1.17884000000000e+01_rt, &
        -1.02446000000000e+00_rt, &
        -2.35700000000000e+01_rt, &
        2.04886000000000e+01_rt, &
        -1.29882000000000e+01_rt, &
        -2.00000000000000e+01_rt, &
        -2.16667000000000e+00_rt ]

    ctemp_rate(:, 3) = [  &
        -9.71052000000000e-01_rt, &
        0.00000000000000e+00_rt, &
        -3.70600000000000e+01_rt, &
        2.93493000000000e+01_rt, &
        -1.15507000000000e+02_rt, &
        -1.00000000000000e+01_rt, &
        -1.33333000000000e+00_rt ]

    ! he4_c12__o16
    ctemp_rate(:, 4) = [  &
        6.96526000000000e+01_rt, &
        -1.39254000000000e+00_rt, &
        5.89128000000000e+01_rt, &
        -1.48273000000000e+02_rt, &
        9.08324000000000e+00_rt, &
        -5.41041000000000e-01_rt, &
        7.03554000000000e+01_rt ]

    ctemp_rate(:, 5) = [  &
        2.54634000000000e+02_rt, &
        -1.84097000000000e+00_rt, &
        1.03411000000000e+02_rt, &
        -4.20567000000000e+02_rt, &
        6.40874000000000e+01_rt, &
        -1.24624000000000e+01_rt, &
        1.37303000000000e+02_rt ]

    ! he4_n14__f18
    ctemp_rate(:, 6) = [  &
        1.38995000000000e+01_rt, &
        -1.09656000000000e+01_rt, &
        -5.62270000000000e+00_rt, &
        0.00000000000000e+00_rt, &
        0.00000000000000e+00_rt, &
        0.00000000000000e+00_rt, &
        -1.50000000000000e+00_rt ]

    ctemp_rate(:, 7) = [  &
        1.96838000000000e-01_rt, &
        -5.16034000000000e+00_rt, &
        0.00000000000000e+00_rt, &
        0.00000000000000e+00_rt, &
        0.00000000000000e+00_rt, &
        0.00000000000000e+00_rt, &
        -1.50000000000000e+00_rt ]

    ctemp_rate(:, 8) = [  &
        2.15339000000000e+01_rt, &
        0.00000000000000e+00_rt, &
        -3.62504000000000e+01_rt, &
        0.00000000000000e+00_rt, &
        0.00000000000000e+00_rt, &
        -5.00000000000000e+00_rt, &
        -6.66667000000000e-01_rt ]

    ! he4_f18__p_ne21
    ctemp_rate(:, 9) = [  &
        4.97863000000000e+01_rt, &
        -1.84559000000000e+00_rt, &
        2.14461000000000e+01_rt, &
        -7.32520000000000e+01_rt, &
        2.42329000000000e+00_rt, &
        -7.72780000000000e-02_rt, &
        4.07604000000000e+01_rt ]

    ! p_c12__n13
    ctemp_rate(:, 10) = [  &
        1.75428000000000e+01_rt, &
        -3.77849000000000e+00_rt, &
        -5.10735000000000e+00_rt, &
        -2.24111000000000e+00_rt, &
        1.48883000000000e-01_rt, &
        0.00000000000000e+00_rt, &
        -1.50000000000000e+00_rt ]

    ctemp_rate(:, 11) = [  &
        1.71482000000000e+01_rt, &
        0.00000000000000e+00_rt, &
        -1.36920000000000e+01_rt, &
        -2.30881000000000e-01_rt, &
        4.44362000000000e+00_rt, &
        -3.15898000000000e+00_rt, &
        -6.66667000000000e-01_rt ]

    ! he4_n13__p_o16
    ctemp_rate(:, 12) = [  &
        4.04644000000000e+01_rt, &
        0.00000000000000e+00_rt, &
        -3.58290000000000e+01_rt, &
        -5.30275000000000e-01_rt, &
        -9.82462000000000e-01_rt, &
        8.08059000000000e-02_rt, &
        -6.66667000000000e-01_rt ]

    ! he4_o16__ne20
    ctemp_rate(:, 13) = [  &
        3.88571000000000e+00_rt, &
        -1.03585000000000e+01_rt, &
        0.00000000000000e+00_rt, &
        0.00000000000000e+00_rt, &
        0.00000000000000e+00_rt, &
        0.00000000000000e+00_rt, &
        -1.50000000000000e+00_rt ]

    ctemp_rate(:, 14) = [  &
        2.39030000000000e+01_rt, &
        0.00000000000000e+00_rt, &
        -3.97262000000000e+01_rt, &
        -2.10799000000000e-01_rt, &
        4.42879000000000e-01_rt, &
        -7.97753000000000e-02_rt, &
        -6.66667000000000e-01_rt ]

    ctemp_rate(:, 15) = [  &
        9.50848000000000e+00_rt, &
        -1.27643000000000e+01_rt, &
        0.00000000000000e+00_rt, &
        -3.65925000000000e+00_rt, &
        7.14224000000000e-01_rt, &
        -1.07508000000000e-03_rt, &
        -1.50000000000000e+00_rt ]

    ! he4_c14__o18
    ctemp_rate(:, 16) = [  &
        -2.38050000000000e+01_rt, &
        -2.06876000000000e+00_rt, &
        0.00000000000000e+00_rt, &
        0.00000000000000e+00_rt, &
        0.00000000000000e+00_rt, &
        0.00000000000000e+00_rt, &
        -1.50000000000000e+00_rt ]

    ctemp_rate(:, 17) = [  &
        1.84877000000000e+01_rt, &
        0.00000000000000e+00_rt, &
        -3.17222000000000e+01_rt, &
        1.13923000000000e+01_rt, &
        -9.92249000000000e+00_rt, &
        -2.00000000000000e+00_rt, &
        -6.66667000000000e-01_rt ]

    ctemp_rate(:, 18) = [  &
        1.18309000000000e+01_rt, &
        -1.03983000000000e+01_rt, &
        0.00000000000000e+00_rt, &
        -3.83188000000000e+00_rt, &
        1.64358000000000e+00_rt, &
        -1.77785000000000e-01_rt, &
        -1.50000000000000e+00_rt ]



    allocate( rate_start_idx(nrat_reaclib) )
    rate_start_idx(:) = [ &
      1, &
      4, &
      6, &
      9, &
      10, &
      12, &
      13, &
      16 ]

    allocate( rate_extra_mult(nrat_reaclib) )
    rate_extra_mult(:) = [ &
      2, &
      1, &
      2, &
      0, &
      1, &
      0, &
      2, &
      2 ]

    allocate( do_screening(nrat_reaclib) )
    do_screening(:) = [ &
      .true., &
      .true., &
      .true., &
      .true., &
      .true., &
      .true., &
      .true., &
      .true. ]

    !$acc update device(ctemp_rate, rate_start_idx, rate_extra_mult, do_screening)
    
  end subroutine init_reaclib

  subroutine term_reaclib()
    deallocate( ctemp_rate )
    deallocate( rate_start_idx )
    deallocate( rate_extra_mult )
    deallocate( do_screening )
  end subroutine term_reaclib

  subroutine net_screening_init()
    ! Adds screening factors and calls screening_init

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jhe4), aion(jhe4))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jc12), aion(jc12))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jn14), aion(jn14))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jf18), aion(jf18))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jc12), aion(jc12))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jn13), aion(jn13))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jo16), aion(jo16))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jc14), aion(jc14))


    call screening_init()    
  end subroutine net_screening_init

  subroutine reaclib_evaluate(pstate, temp, iwhich, reactvec)
    !$acc routine seq

    implicit none
    
    type(plasma_state), intent(in) :: pstate
    real(rt), intent(in) :: temp
    integer, intent(in) :: iwhich

    real(rt), intent(inout) :: reactvec(num_rate_groups+2)
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

    real(rt)  :: rate, scor ! Rate and Screening Factor
    real(rt)  :: drate_dt, dscor_dt ! Temperature derivatives
    real(rt) :: dscor_dd
    real(rt) :: ri, T9, T9_exp, lnirate, irate, dirate_dt, dlnirate_dt
    integer :: i, j, m, istart

    !$gpu
    ri = 0.0e0_rt
    rate = 0.0e0_rt
    drate_dt = 0.0e0_rt
    irate = 0.0e0_rt
    dirate_dt = 0.0e0_rt
    T9 = temp/1.0e9_rt
    T9_exp = 0.0e0_rt
    scor = 1.0e0_rt
    dscor_dt = 0.0e0_rt
    dscor_dd = 0.0e0_rt

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

    if ( screen_reaclib .and. do_screening(iwhich) ) then
       call screen5(pstate, iwhich, scor, dscor_dt, dscor_dd)
    end if

    reactvec(i_rate)     = rate
    reactvec(i_drate_dt) = drate_dt
    reactvec(i_scor)     = scor
    reactvec(i_dscor_dt) = dscor_dt
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
