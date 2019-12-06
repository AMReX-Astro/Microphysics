module reaclib_rates

  use microphysics_type_module
  use screening_module, only: add_screening_factor, &
                              screening_init, screening_finalize, &
                              plasma_state, fill_plasma_state
  use network

  implicit none

  logical, parameter :: screen_reaclib = .true.

  ! Temperature coefficient arrays (numbers correspond to reaction numbers in net_info)
  real(rt), allocatable :: ctemp_rate(:,:)

  ! Index into ctemp_rate, dimension 2, where each rate's coefficients start
  integer, allocatable :: rate_start_idx(:)

  ! Reaction multiplicities-1 (how many rates contribute - 1)
  integer, allocatable :: rate_extra_mult(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: ctemp_rate, rate_start_idx, rate_extra_mult
#endif

  !$acc declare create(ctemp_rate, rate_start_idx, rate_extra_mult)
  !$acc declare copyin(screen_reaclib)

contains

  subroutine init_reaclib()

    allocate( ctemp_rate(7, 6) )
    ! c12_c12__he4_ne20
    ctemp_rate(:, 1) = [  &
        6.12863000000000e+01_rt, &
        0.00000000000000e+00_rt, &
        -8.41650000000000e+01_rt, &
        -1.56627000000000e+00_rt, &
        -7.36084000000000e-02_rt, &
        -7.27970000000000e-02_rt, &
        -6.66667000000000e-01_rt ]

    ! c12_c12__n_mg23
    ctemp_rate(:, 2) = [  &
        -1.28056000000000e+01_rt, &
        -3.01485000000000e+01_rt, &
        0.00000000000000e+00_rt, &
        1.14826000000000e+01_rt, &
        1.82849000000000e+00_rt, &
        -3.48440000000000e-01_rt, &
        0.00000000000000e+00_rt ]

    ! c12_c12__p_na23
    ctemp_rate(:, 3) = [  &
        6.09649000000000e+01_rt, &
        0.00000000000000e+00_rt, &
        -8.41650000000000e+01_rt, &
        -1.41910000000000e+00_rt, &
        -1.14619000000000e-01_rt, &
        -7.03070000000000e-02_rt, &
        -6.66667000000000e-01_rt ]

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

    ! n__p__weak__wc12
    ctemp_rate(:, 6) = [  &
        -6.78161000000000e+00_rt, &
        0.00000000000000e+00_rt, &
        0.00000000000000e+00_rt, &
        0.00000000000000e+00_rt, &
        0.00000000000000e+00_rt, &
        0.00000000000000e+00_rt, &
        0.00000000000000e+00_rt ]



    allocate( rate_start_idx(nrat_reaclib) )
    rate_start_idx(:) = [ &
      1, &
      2, &
      3, &
      4, &
      6 ]

    allocate( rate_extra_mult(nrat_reaclib) )
    rate_extra_mult(:) = [ &
      0, &
      0, &
      0, &
      1, &
      0 ]

    !$acc update device(ctemp_rate, rate_start_idx, rate_extra_mult)

  end subroutine init_reaclib

  subroutine term_reaclib()
    deallocate( ctemp_rate )
    deallocate( rate_start_idx )
    deallocate( rate_extra_mult )
  end subroutine term_reaclib


  subroutine net_screening_init()
    ! Adds screening factors and calls screening_init

    call add_screening_factor(zion(jc12), aion(jc12), &
      zion(jc12), aion(jc12))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jc12), aion(jc12))


    call screening_init()
  end subroutine net_screening_init


  subroutine net_screening_finalize()
    ! Call screening_finalize

    call screening_finalize()

  end subroutine net_screening_finalize


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

    real(rt) :: rate, scor ! Rate and Screening Factor
    real(rt) :: drate_dt, dscor_dt ! Temperature derivatives
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
    reactvec(i_scor)     = 1.0e0_rt
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
