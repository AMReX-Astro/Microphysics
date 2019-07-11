module reaclib_rates

  use amrex_fort_module, only: rt => amrex_real
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

    implicit none

    integer :: unit, ireaclib, icoeff

    allocate( ctemp_rate(7, number_reaclib_sets) )
    allocate( rate_start_idx(nrat_reaclib) )
    allocate( rate_extra_mult(nrat_reaclib) )

    open(newunit=unit, file='reaclib_rate_metadata.dat')

    do ireaclib = 1, number_reaclib_sets
       do icoeff = 1, 7
          read(unit, *) ctemp_rate(icoeff, ireaclib)
       enddo
    enddo

    do ireaclib = 1, nrat_reaclib
       read(unit, *) rate_start_idx(ireaclib)
    enddo

    do ireaclib = 1, nrat_reaclib
       read(unit, *) rate_extra_mult(ireaclib)
    enddo

    close(unit)

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
      zion(jp), aion(jp))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jd), aion(jd))

    call add_screening_factor(zion(jd), aion(jd), &
      zion(jd), aion(jd))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jhe3), aion(jhe3))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jhe3), aion(jhe3))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jbe7), aion(jbe7))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jc12), aion(jc12))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jc12), aion(jc12))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jn13), aion(jn13))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jn14), aion(jn14))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jn14), aion(jn14))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jn15), aion(jn15))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jn15), aion(jn15))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jo14), aion(jo14))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jo15), aion(jo15))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jo16), aion(jo16))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jo16), aion(jo16))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jo17), aion(jo17))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jo17), aion(jo17))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jo18), aion(jo18))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jf17), aion(jf17))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jf18), aion(jf18))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jf19), aion(jf19))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jf20), aion(jf20))

    call add_screening_factor(zion(jd), aion(jd), &
      zion(jhe3), aion(jhe3))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jhe4), aion(jhe4))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jhe4), aion(jhe4))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jli7), aion(jli7))

    call add_screening_factor(zion(jc12), aion(jc12), &
      zion(jc12), aion(jc12))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jn13), aion(jn13))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jf17), aion(jf17))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jf18), aion(jf18))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jne20), aion(jne20))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jne20), aion(jne20))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jne21), aion(jne21))

    call add_screening_factor(zion(jhe3), aion(jhe3), &
      zion(jhe3), aion(jhe3))

    call add_screening_factor(zion(jd), aion(jd), &
      zion(jbe7), aion(jbe7))

    call add_screening_factor(zion(jhe3), aion(jhe3), &
      zion(jbe7), aion(jbe7))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jo15), aion(jo15))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jo18), aion(jo18))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jf19), aion(jf19))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jf20), aion(jf20))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jne18), aion(jne18))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jne18), aion(jne18))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jne19), aion(jne19))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jne19), aion(jne19))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jne21), aion(jne21))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jc13), aion(jc13))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jf16), aion(jf16))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jne22), aion(jne22))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jne22), aion(jne22))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jna21), aion(jna21))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jna21), aion(jna21))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jna22), aion(jna22))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jna22), aion(jna22))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jna23), aion(jna23))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jna23), aion(jna23))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jna24), aion(jna24))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jna24), aion(jna24))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jmg22), aion(jmg22))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmg22), aion(jmg22))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jna19), aion(jna19))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jna20), aion(jna20))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jna20), aion(jna20))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jmg23), aion(jmg23))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmg23), aion(jmg23))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jmg24), aion(jmg24))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmg24), aion(jmg24))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jmg25), aion(jmg25))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmg25), aion(jmg25))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jmg26), aion(jmg26))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmg26), aion(jmg26))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jal25), aion(jal25))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jal25), aion(jal25))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jal26), aion(jal26))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jal26), aion(jal26))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jal27), aion(jal27))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jal27), aion(jal27))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jal28), aion(jal28))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jal28), aion(jal28))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsi26), aion(jsi26))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsi26), aion(jsi26))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jal23), aion(jal23))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jal23), aion(jal23))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jal24), aion(jal24))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jal24), aion(jal24))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jmg21), aion(jmg21))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmg21), aion(jmg21))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsi27), aion(jsi27))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsi27), aion(jsi27))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsi28), aion(jsi28))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsi28), aion(jsi28))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsi29), aion(jsi29))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsi29), aion(jsi29))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsi30), aion(jsi30))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsi30), aion(jsi30))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jp29), aion(jp29))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jp29), aion(jp29))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jp30), aion(jp30))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jp30), aion(jp30))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jp31), aion(jp31))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jp31), aion(jp31))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jp32), aion(jp32))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jp32), aion(jp32))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(js30), aion(js30))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(js30), aion(js30))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jp27), aion(jp27))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jp27), aion(jp27))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsi24), aion(jsi24))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jp28), aion(jp28))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jp28), aion(jp28))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsi25), aion(jsi25))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsi25), aion(jsi25))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jal22), aion(jal22))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(js31), aion(js31))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(js31), aion(js31))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(js32), aion(js32))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(js32), aion(js32))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(js33), aion(js33))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(js33), aion(js33))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(js34), aion(js34))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(js34), aion(js34))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcl33), aion(jcl33))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcl33), aion(jcl33))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcl34), aion(jcl34))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcl34), aion(jcl34))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcl35), aion(jcl35))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcl35), aion(jcl35))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcl36), aion(jcl36))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcl36), aion(jcl36))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcl31), aion(jcl31))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcl31), aion(jcl31))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jar34), aion(jar34))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jar34), aion(jar34))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(js28), aion(js28))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(js28), aion(js28))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcl32), aion(jcl32))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcl32), aion(jcl32))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(js29), aion(js29))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(js29), aion(js29))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jp26), aion(jp26))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jar35), aion(jar35))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jar35), aion(jar35))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jar36), aion(jar36))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jar36), aion(jar36))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jar37), aion(jar37))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jar37), aion(jar37))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jar38), aion(jar38))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jar38), aion(jar38))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jk37), aion(jk37))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jk37), aion(jk37))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jk38), aion(jk38))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jk38), aion(jk38))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jk39), aion(jk39))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jk39), aion(jk39))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jk40), aion(jk40))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jk40), aion(jk40))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jar32), aion(jar32))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jar32), aion(jar32))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jk35), aion(jk35))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jk35), aion(jk35))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jca38), aion(jca38))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jca38), aion(jca38))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcl29), aion(jcl29))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jar33), aion(jar33))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jar33), aion(jar33))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jk36), aion(jk36))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jk36), aion(jk36))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcl30), aion(jcl30))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcl30), aion(jcl30))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jca39), aion(jca39))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jca39), aion(jca39))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jca40), aion(jca40))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jca40), aion(jca40))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcl37), aion(jcl37))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcl37), aion(jcl37))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jca41), aion(jca41))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jca41), aion(jca41))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jca42), aion(jca42))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jca42), aion(jca42))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsc41), aion(jsc41))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsc41), aion(jsc41))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsc42), aion(jsc42))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsc42), aion(jsc42))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsc43), aion(jsc43))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsc43), aion(jsc43))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsc44), aion(jsc44))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsc44), aion(jsc44))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jca36), aion(jca36))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jca36), aion(jca36))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jk33), aion(jk33))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jk33), aion(jk33))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsc39), aion(jsc39))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsc39), aion(jsc39))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jti42), aion(jti42))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jti42), aion(jti42))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jk34), aion(jk34))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jk34), aion(jk34))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jca37), aion(jca37))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jca37), aion(jca37))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsc40), aion(jsc40))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsc40), aion(jsc40))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jar31), aion(jar31))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jti43), aion(jti43))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jti43), aion(jti43))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jti44), aion(jti44))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jti44), aion(jti44))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jk41), aion(jk41))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jk41), aion(jk41))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jti45), aion(jti45))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jti45), aion(jti45))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jti46), aion(jti46))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jti46), aion(jti46))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jv45), aion(jv45))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jv45), aion(jv45))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jv46), aion(jv46))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jv46), aion(jv46))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jca43), aion(jca43))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jca43), aion(jca43))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jv47), aion(jv47))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jv47), aion(jv47))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jv48), aion(jv48))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jv48), aion(jv48))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jca44), aion(jca44))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jca44), aion(jca44))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jti40), aion(jti40))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jti40), aion(jti40))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsc37), aion(jsc37))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsc37), aion(jsc37))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jca34), aion(jca34))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jv43), aion(jv43))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jv43), aion(jv43))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcr46), aion(jcr46))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcr46), aion(jcr46))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jca35), aion(jca35))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jca35), aion(jca35))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsc38), aion(jsc38))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsc38), aion(jsc38))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jti41), aion(jti41))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jti41), aion(jti41))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jv44), aion(jv44))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jv44), aion(jv44))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcr47), aion(jcr47))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcr47), aion(jcr47))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcr48), aion(jcr48))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcr48), aion(jcr48))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsc45), aion(jsc45))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsc45), aion(jsc45))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcr49), aion(jcr49))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcr49), aion(jcr49))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcr50), aion(jcr50))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcr50), aion(jcr50))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jmn49), aion(jmn49))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmn49), aion(jmn49))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jmn50), aion(jmn50))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmn50), aion(jmn50))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jti47), aion(jti47))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jti47), aion(jti47))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jmn51), aion(jmn51))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmn51), aion(jmn51))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jti48), aion(jti48))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jti48), aion(jti48))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jmn52), aion(jmn52))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmn52), aion(jmn52))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcr44), aion(jcr44))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcr44), aion(jcr44))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jv41), aion(jv41))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jv41), aion(jv41))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jti38), aion(jti38))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jmn47), aion(jmn47))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmn47), aion(jmn47))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jfe50), aion(jfe50))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jfe50), aion(jfe50))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsc36), aion(jsc36))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jti39), aion(jti39))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jti39), aion(jti39))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jv42), aion(jv42))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jv42), aion(jv42))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcr45), aion(jcr45))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcr45), aion(jcr45))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jmn48), aion(jmn48))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmn48), aion(jmn48))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jfe51), aion(jfe51))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jfe51), aion(jfe51))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jfe52), aion(jfe52))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jfe52), aion(jfe52))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jv49), aion(jv49))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jv49), aion(jv49))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jfe53), aion(jfe53))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jfe54), aion(jfe54))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jco53), aion(jco53))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jco54), aion(jco54))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcr51), aion(jcr51))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcr51), aion(jcr51))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jco55), aion(jco55))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcr52), aion(jcr52))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcr52), aion(jcr52))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jmn45), aion(jmn45))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmn45), aion(jmn45))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jfe48), aion(jfe48))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jfe48), aion(jfe48))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcr42), aion(jcr42))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jco51), aion(jco51))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jv40), aion(jv40))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcr43), aion(jcr43))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcr43), aion(jcr43))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jmn46), aion(jmn46))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmn46), aion(jmn46))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jfe49), aion(jfe49))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jfe49), aion(jfe49))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jco52), aion(jco52))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jmn53), aion(jmn53))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jfe55), aion(jfe55))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jfe46), aion(jfe46))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jfe46), aion(jfe46))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jco49), aion(jco49))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jmn44), aion(jmn44))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmn44), aion(jmn44))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jfe47), aion(jfe47))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jfe47), aion(jfe47))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jco50), aion(jco50))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jmn55), aion(jmn55))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jco47), aion(jco47))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jfe45), aion(jfe45))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jco48), aion(jco48))


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
    !                           For Reaclib rates, this is 0.0d0.)
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

    ri = 0.0d0
    rate = 0.0d0
    drate_dt = 0.0d0
    irate = 0.0d0
    dirate_dt = 0.0d0
    T9 = temp/1.0d9
    T9_exp = 0.0d0

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

    reactvec(i_rate)     = rate
    reactvec(i_drate_dt) = drate_dt
    reactvec(i_scor)     = 1.0d0
    reactvec(i_dscor_dt) = 0.0d0
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
