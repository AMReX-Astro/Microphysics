module reaclib_rates
  use screening_module, only: screen5, add_screening_factor, screening_init, plasma_state, fill_plasma_state
  use network

  implicit none

  logical, parameter :: screen_reaclib = .true.
  
  ! Temperature coefficient arrays (numbers correspond to reaction numbers in net_info)
  double precision, allocatable :: ctemp_rate(:,:)

  ! Index into ctemp_rate, dimension 2, where each rate's coefficients start
  integer, allocatable :: rate_start_idx(:)
  
  ! Reaction multiplicities-1 (how many rates contribute - 1)
  integer, allocatable :: rate_extra_mult(:)

  ! Should these reactions be screened?
  logical, allocatable :: do_screening(:)
  
  !$acc declare create(ctemp_rate, rate_start_idx, rate_extra_mult, do_screening)
  !$acc declare copyin(screen_reaclib)

contains

  subroutine init_reaclib()
    
    allocate( ctemp_rate(7, 52) )
    ! o16__he4_c12
    ctemp_rate(:, 1) = [  &
        2.79295000000000d+02, &
        -8.49515000000000d+01, &
        1.03411000000000d+02, &
        -4.20567000000000d+02, &
        6.40874000000000d+01, &
        -1.24624000000000d+01, &
        1.38803000000000d+02 ]

    ctemp_rate(:, 2) = [  &
        9.43131000000000d+01, &
        -8.45030000000000d+01, &
        5.89128000000000d+01, &
        -1.48273000000000d+02, &
        9.08324000000000d+00, &
        -5.41041000000000d-01, &
        7.18554000000000d+01 ]

    ! ne20__he4_o16
    ctemp_rate(:, 3) = [  &
        3.42658000000000d+01, &
        -6.76518000000000d+01, &
        0.00000000000000d+00, &
        -3.65925000000000d+00, &
        7.14224000000000d-01, &
        -1.07508000000000d-03, &
        0.00000000000000d+00 ]

    ctemp_rate(:, 4) = [  &
        2.86431000000000d+01, &
        -6.52460000000000d+01, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00 ]

    ctemp_rate(:, 5) = [  &
        4.86604000000000d+01, &
        -5.48875000000000d+01, &
        -3.97262000000000d+01, &
        -2.10799000000000d-01, &
        4.42879000000000d-01, &
        -7.97753000000000d-02, &
        8.33333000000000d-01 ]

    ! mg24__he4_ne20
    ctemp_rate(:, 6) = [  &
        4.93244000000000d+01, &
        -1.08114000000000d+02, &
        -4.62525000000000d+01, &
        5.58901000000000d+00, &
        7.61843000000000d+00, &
        -3.68300000000000d+00, &
        8.33333000000000d-01 ]

    ctemp_rate(:, 7) = [  &
        1.60203000000000d+01, &
        -1.20895000000000d+02, &
        0.00000000000000d+00, &
        1.69229000000000d+01, &
        -2.57325000000000d+00, &
        2.08997000000000d-01, &
        0.00000000000000d+00 ]

    ctemp_rate(:, 8) = [  &
        2.68017000000000d+01, &
        -1.17334000000000d+02, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00 ]

    ctemp_rate(:, 9) = [  &
        -1.38869000000000d+01, &
        -1.10620000000000d+02, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00 ]

    ! si28__he4_mg24
    ctemp_rate(:, 10) = [  &
        3.29006000000000d+01, &
        -1.31488000000000d+02, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00 ]

    ctemp_rate(:, 11) = [  &
        -2.56886000000000d+01, &
        -1.28693000000000d+02, &
        2.13721000000000d+01, &
        3.77649000000000d+01, &
        -4.10635000000000d+00, &
        2.49618000000000d-01, &
        0.00000000000000d+00 ]

    ! s32__he4_si28
    ctemp_rate(:, 12) = [  &
        7.28130000000000d+01, &
        -8.06260000000000d+01, &
        -5.94896000000000d+01, &
        4.47205000000000d+00, &
        -4.78989000000000d+00, &
        5.57201000000000d-01, &
        8.33333000000000d-01 ]

    ! ar36__he4_s32
    ctemp_rate(:, 13) = [  &
        7.38164000000000d+01, &
        -7.70627000000000d+01, &
        -6.53709000000000d+01, &
        5.68294000000000d+00, &
        -5.00388000000000d+00, &
        5.71407000000000d-01, &
        8.33333000000000d-01 ]

    ! ca40__he4_ar36
    ctemp_rate(:, 14) = [  &
        7.72826000000000d+01, &
        -8.16916000000000d+01, &
        -7.10046000000000d+01, &
        4.06560000000000d+00, &
        -5.26509000000000d+00, &
        6.83546000000000d-01, &
        8.33333000000000d-01 ]

    ! ti44__he4_ca40
    ctemp_rate(:, 15) = [  &
        7.86991000000000d+01, &
        -5.94974000000000d+01, &
        -7.64273000000000d+01, &
        3.87451000000000d+00, &
        -3.61477000000000d+00, &
        3.67451000000000d-01, &
        8.33333000000000d-01 ]

    ! cr48__he4_ti44
    ctemp_rate(:, 16) = [  &
        8.97573000000000d+01, &
        -8.93041000000000d+01, &
        -8.16670000000000d+01, &
        -1.06333000000000d+01, &
        -6.72613000000000d-01, &
        1.61209000000000d-01, &
        8.33333000000000d-01 ]

    ! fe52__he4_cr48
    ctemp_rate(:, 17) = [  &
        9.01474000000000d+01, &
        -9.21090000000000d+01, &
        -8.67459000000000d+01, &
        -9.79373000000000d+00, &
        -7.72169000000000d-01, &
        1.55883000000000d-01, &
        8.33333000000000d-01 ]

    ! ni56__he4_fe52
    ctemp_rate(:, 18) = [  &
        9.16226000000000d+01, &
        -9.28010000000000d+01, &
        -9.16819000000000d+01, &
        -9.51885000000000d+00, &
        -5.33014000000000d-01, &
        8.92607000000000d-02, &
        8.33333000000000d-01 ]

    ! zn60__he4_ni56
    ctemp_rate(:, 19) = [  &
        8.60619000000000d+01, &
        -3.14367000000000d+01, &
        -9.64898000000000d+01, &
        6.47209000000000d+00, &
        -5.20290000000000d+00, &
        5.33391000000000d-01, &
        8.33333000000000d-01 ]

    ! c12__he4_he4_he4
    ctemp_rate(:, 20) = [  &
        2.23940000000000d+01, &
        -8.85493000000000d+01, &
        -1.34900000000000d+01, &
        2.14259000000000d+01, &
        -1.34769000000000d+00, &
        8.79816000000000d-02, &
        -1.01653000000000d+01 ]

    ctemp_rate(:, 21) = [  &
        3.49561000000000d+01, &
        -8.54472000000000d+01, &
        -2.35700000000000d+01, &
        2.04886000000000d+01, &
        -1.29882000000000d+01, &
        -2.00000000000000d+01, &
        8.33330000000000d-01 ]

    ctemp_rate(:, 22) = [  &
        4.57734000000000d+01, &
        -8.44227000000000d+01, &
        -3.70600000000000d+01, &
        2.93493000000000d+01, &
        -1.15507000000000d+02, &
        -1.00000000000000d+01, &
        1.66667000000000d+00 ]

    ! he4_c12__o16
    ctemp_rate(:, 23) = [  &
        6.96526000000000d+01, &
        -1.39254000000000d+00, &
        5.89128000000000d+01, &
        -1.48273000000000d+02, &
        9.08324000000000d+00, &
        -5.41041000000000d-01, &
        7.03554000000000d+01 ]

    ctemp_rate(:, 24) = [  &
        2.54634000000000d+02, &
        -1.84097000000000d+00, &
        1.03411000000000d+02, &
        -4.20567000000000d+02, &
        6.40874000000000d+01, &
        -1.24624000000000d+01, &
        1.37303000000000d+02 ]

    ! he4_o16__ne20
    ctemp_rate(:, 25) = [  &
        3.88571000000000d+00, &
        -1.03585000000000d+01, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 26) = [  &
        2.39030000000000d+01, &
        0.00000000000000d+00, &
        -3.97262000000000d+01, &
        -2.10799000000000d-01, &
        4.42879000000000d-01, &
        -7.97753000000000d-02, &
        -6.66667000000000d-01 ]

    ctemp_rate(:, 27) = [  &
        9.50848000000000d+00, &
        -1.27643000000000d+01, &
        0.00000000000000d+00, &
        -3.65925000000000d+00, &
        7.14224000000000d-01, &
        -1.07508000000000d-03, &
        -1.50000000000000d+00 ]

    ! he4_ne20__mg24
    ctemp_rate(:, 28) = [  &
        -8.79827000000000d+00, &
        -1.27809000000000d+01, &
        0.00000000000000d+00, &
        1.69229000000000d+01, &
        -2.57325000000000d+00, &
        2.08997000000000d-01, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 29) = [  &
        1.98307000000000d+00, &
        -9.22026000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 30) = [  &
        -3.87055000000000d+01, &
        -2.50605000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 31) = [  &
        2.45058000000000d+01, &
        0.00000000000000d+00, &
        -4.62525000000000d+01, &
        5.58901000000000d+00, &
        7.61843000000000d+00, &
        -3.68300000000000d+00, &
        -6.66667000000000d-01 ]

    ! he4_mg24__si28
    ctemp_rate(:, 32) = [  &
        -5.05494000000000d+01, &
        -1.28332000000000d+01, &
        2.13721000000000d+01, &
        3.77649000000000d+01, &
        -4.10635000000000d+00, &
        2.49618000000000d-01, &
        -1.50000000000000d+00 ]

    ctemp_rate(:, 33) = [  &
        8.03977000000000d+00, &
        -1.56290000000000d+01, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        0.00000000000000d+00, &
        -1.50000000000000d+00 ]

    ! he4_si28__s32
    ctemp_rate(:, 34) = [  &
        4.79212000000000d+01, &
        0.00000000000000d+00, &
        -5.94896000000000d+01, &
        4.47205000000000d+00, &
        -4.78989000000000d+00, &
        5.57201000000000d-01, &
        -6.66667000000000d-01 ]

    ! he4_s32__ar36
    ctemp_rate(:, 35) = [  &
        4.89010000000000d+01, &
        0.00000000000000d+00, &
        -6.53709000000000d+01, &
        5.68294000000000d+00, &
        -5.00388000000000d+00, &
        5.71407000000000d-01, &
        -6.66667000000000d-01 ]

    ! he4_ar36__ca40
    ctemp_rate(:, 36) = [  &
        5.23486000000000d+01, &
        0.00000000000000d+00, &
        -7.10046000000000d+01, &
        4.06560000000000d+00, &
        -5.26509000000000d+00, &
        6.83546000000000d-01, &
        -6.66667000000000d-01 ]

    ! he4_ca40__ti44
    ctemp_rate(:, 37) = [  &
        5.37500000000000d+01, &
        0.00000000000000d+00, &
        -7.64273000000000d+01, &
        3.87451000000000d+00, &
        -3.61477000000000d+00, &
        3.67451000000000d-01, &
        -6.66667000000000d-01 ]

    ! he4_ti44__cr48
    ctemp_rate(:, 38) = [  &
        6.47958000000000d+01, &
        0.00000000000000d+00, &
        -8.16670000000000d+01, &
        -1.06333000000000d+01, &
        -6.72613000000000d-01, &
        1.61209000000000d-01, &
        -6.66667000000000d-01 ]

    ! he4_cr48__fe52
    ctemp_rate(:, 39) = [  &
        6.51754000000000d+01, &
        0.00000000000000d+00, &
        -8.67459000000000d+01, &
        -9.79373000000000d+00, &
        -7.72169000000000d-01, &
        1.55883000000000d-01, &
        -6.66667000000000d-01 ]

    ! he4_fe52__ni56
    ctemp_rate(:, 40) = [  &
        6.66417000000000d+01, &
        0.00000000000000d+00, &
        -9.16819000000000d+01, &
        -9.51885000000000d+00, &
        -5.33014000000000d-01, &
        8.92607000000000d-02, &
        -6.66667000000000d-01 ]

    ! he4_ni56__zn60
    ctemp_rate(:, 41) = [  &
        6.10733000000000d+01, &
        0.00000000000000d+00, &
        -9.64898000000000d+01, &
        6.47209000000000d+00, &
        -5.20290000000000d+00, &
        5.33391000000000d-01, &
        -6.66667000000000d-01 ]

    ! c12_c12__he4_ne20
    ctemp_rate(:, 42) = [  &
        6.12863000000000d+01, &
        0.00000000000000d+00, &
        -8.41650000000000d+01, &
        -1.56627000000000d+00, &
        -7.36084000000000d-02, &
        -7.27970000000000d-02, &
        -6.66667000000000d-01 ]

    ! c12_o16__he4_mg24
    ctemp_rate(:, 43) = [  &
        4.85341000000000d+01, &
        3.72040000000000d-01, &
        -1.33413000000000d+02, &
        5.01572000000000d+01, &
        -3.15987000000000d+00, &
        1.78251000000000d-02, &
        -2.37027000000000d+01 ]

    ! o16_o16__he4_si28
    ctemp_rate(:, 44) = [  &
        9.72435000000000d+01, &
        -2.68514000000000d-01, &
        -1.19324000000000d+02, &
        -3.22497000000000d+01, &
        1.46214000000000d+00, &
        -2.00893000000000d-01, &
        1.32148000000000d+01 ]

    ! he4_ne20__c12_c12
    ctemp_rate(:, 45) = [  &
        6.14748000000000d+01, &
        -5.36267000000000d+01, &
        -8.41650000000000d+01, &
        -1.56627000000000d+00, &
        -7.36084000000000d-02, &
        -7.27970000000000d-02, &
        -6.66667000000000d-01 ]

    ! c12_ne20__he4_si28
    ctemp_rate(:, 46) = [  &
        -3.08905000000000d+02, &
        -4.72175000000000d+01, &
        5.14197000000000d+02, &
        -2.00896000000000d+02, &
        -6.42713000000000d+00, &
        7.58256000000000d-01, &
        2.36359000000000d+02 ]

    ! he4_mg24__c12_o16
    ctemp_rate(:, 47) = [  &
        4.95738000000000d+01, &
        -7.82020000000000d+01, &
        -1.33413000000000d+02, &
        5.01572000000000d+01, &
        -3.15987000000000d+00, &
        1.78251000000000d-02, &
        -2.37027000000000d+01 ]

    ! he4_si28__c12_ne20
    ctemp_rate(:, 48) = [  &
        -3.07762000000000d+02, &
        -1.86722000000000d+02, &
        5.14197000000000d+02, &
        -2.00896000000000d+02, &
        -6.42713000000000d+00, &
        7.58256000000000d-01, &
        2.36359000000000d+02 ]

    ! he4_si28__o16_o16
    ctemp_rate(:, 49) = [  &
        9.77904000000000d+01, &
        -1.11595000000000d+02, &
        -1.19324000000000d+02, &
        -3.22497000000000d+01, &
        1.46214000000000d+00, &
        -2.00893000000000d-01, &
        1.32148000000000d+01 ]

    ! he4_he4_he4__c12
    ctemp_rate(:, 50) = [  &
        -2.43505000000000d+01, &
        -4.12656000000000d+00, &
        -1.34900000000000d+01, &
        2.14259000000000d+01, &
        -1.34769000000000d+00, &
        8.79816000000000d-02, &
        -1.31653000000000d+01 ]

    ctemp_rate(:, 51) = [  &
        -1.17884000000000d+01, &
        -1.02446000000000d+00, &
        -2.35700000000000d+01, &
        2.04886000000000d+01, &
        -1.29882000000000d+01, &
        -2.00000000000000d+01, &
        -2.16667000000000d+00 ]

    ctemp_rate(:, 52) = [  &
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
      3, &
      6, &
      10, &
      12, &
      13, &
      14, &
      15, &
      16, &
      17, &
      18, &
      19, &
      20, &
      23, &
      25, &
      28, &
      32, &
      34, &
      35, &
      36, &
      37, &
      38, &
      39, &
      40, &
      41, &
      42, &
      43, &
      44, &
      45, &
      46, &
      47, &
      48, &
      49, &
      50 ]

    allocate( rate_extra_mult(nrat_reaclib) )
    rate_extra_mult(:) = [ &
      1, &
      2, &
      3, &
      1, &
      0, &
      0, &
      0, &
      0, &
      0, &
      0, &
      0, &
      0, &
      2, &
      1, &
      2, &
      3, &
      1, &
      0, &
      0, &
      0, &
      0, &
      0, &
      0, &
      0, &
      0, &
      0, &
      0, &
      0, &
      0, &
      0, &
      0, &
      0, &
      0, &
      2 ]

    allocate( do_screening(nrat_reaclib) )
    do_screening(:) = [ &
      .false., &
      .false., &
      .false., &
      .false., &
      .false., &
      .false., &
      .false., &
      .false., &
      .false., &
      .false., &
      .false., &
      .false., &
      .false., &
      .true., &
      .true., &
      .true., &
      .true., &
      .true., &
      .true., &
      .true., &
      .true., &
      .true., &
      .true., &
      .true., &
      .true., &
      .true., &
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
      zion(jc12), aion(jc12))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jo16), aion(jo16))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jne20), aion(jne20))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmg24), aion(jmg24))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsi28), aion(jsi28))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(js32), aion(js32))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jar36), aion(jar36))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jca40), aion(jca40))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jti44), aion(jti44))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcr48), aion(jcr48))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jfe52), aion(jfe52))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jni56), aion(jni56))

    call add_screening_factor(zion(jc12), aion(jc12), &
      zion(jc12), aion(jc12))

    call add_screening_factor(zion(jc12), aion(jc12), &
      zion(jo16), aion(jo16))

    call add_screening_factor(zion(jo16), aion(jo16), &
      zion(jo16), aion(jo16))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jne20), aion(jne20))

    call add_screening_factor(zion(jc12), aion(jc12), &
      zion(jne20), aion(jne20))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmg24), aion(jmg24))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsi28), aion(jsi28))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsi28), aion(jsi28))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jhe4), aion(jhe4))


    call screening_init()    
  end subroutine net_screening_init

  subroutine reaclib_evaluate(pstate, temp, iwhich, reactvec)
    !$acc routine seq

    implicit none
    
    type(plasma_state), intent(in) :: pstate
    double precision, intent(in) :: temp
    integer, intent(in) :: iwhich

    double precision, intent(inout) :: reactvec(num_rate_groups+2)
    ! reactvec(1) = rate     , the reaction rate
    ! reactvec(2) = drate_dt , the Temperature derivative of rate
    ! reactvec(3) = scor     , the screening factor
    ! reactvec(4) = dscor_dt , the Temperature derivative of scor
    ! reactvec(5) = dqweak   , the weak reaction dq-value (ergs)
    !                          (This accounts for modification of the reaction Q
    !                           due to the local density and temperature of the plasma.
    !                           For Reaclib rates, this is 0.0d0.)
    ! reactvec(6) = epart    , the particle energy generation rate (ergs/s)
    ! NOTE: The particle energy generation rate (returned in ergs/s)
    !       is the contribution to enuc from non-ion particles associated
    !       with the reaction.
    !       For example, this accounts for neutrino energy losses
    !       in weak reactions and/or gamma heating of the plasma
    !       from nuclear transitions in daughter nuclei.

    double precision  :: rate, scor ! Rate and Screening Factor
    double precision  :: drate_dt, dscor_dt ! Temperature derivatives
    double precision :: dscor_dd
    double precision :: ri, T9, T9_exp, lnirate, irate, dirate_dt, dlnirate_dt
    integer :: i, j, m, istart

    ri = 0.0d0
    rate = 0.0d0
    drate_dt = 0.0d0
    irate = 0.0d0
    dirate_dt = 0.0d0
    T9 = temp/1.0d9
    T9_exp = 0.0d0
    scor = 1.0d0
    dscor_dt = 0.0d0
    dscor_dd = 0.0d0

    ! Use reaction multiplicities to tell whether the rate is Reaclib
    m = rate_extra_mult(iwhich)

    istart = rate_start_idx(iwhich)

    do i = 0, m
       lnirate = ctemp_rate(1, istart+i) + ctemp_rate(7, istart+i) * LOG(T9)
       dlnirate_dt = ctemp_rate(7, istart+i)/T9
       do j = 2, 6
          T9_exp = (2.0d0*dble(j-1)-5.0d0)/3.0d0 
          lnirate = lnirate + ctemp_rate(j, istart+i) * T9**T9_exp
          dlnirate_dt = dlnirate_dt + &
               T9_exp * ctemp_rate(j, istart+i) * T9**(T9_exp-1.0d0)
       end do
       ! If the rate will be in the approx. interval [0.0, 1.0E-100], replace by 0.0
       ! This avoids issues with passing very large negative values to EXP
       ! and getting results between 0.0 and 1.0E-308, the limit for IEEE 754.
       ! And avoids SIGFPE in CVODE due to tiny rates.
       lnirate = max(lnirate, -230.0d0)
       irate = EXP(lnirate)
       rate = rate + irate
       dirate_dt = irate * dlnirate_dt/1.0d9
       drate_dt = drate_dt + dirate_dt
    end do

    if ( screen_reaclib .and. do_screening(iwhich) ) then
       call screen5(pstate, iwhich, scor, dscor_dt, dscor_dd)
    end if

    reactvec(i_rate)     = rate
    reactvec(i_drate_dt) = drate_dt
    reactvec(i_scor)     = scor
    reactvec(i_dscor_dt) = dscor_dt
    reactvec(i_dqweak)   = 0.0d0
    reactvec(i_epart)    = 0.0d0

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
