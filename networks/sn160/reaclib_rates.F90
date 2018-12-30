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

    call add_screening_factor(zion(jd), aion(jd), &
      zion(jhe4), aion(jhe4))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jhe3), aion(jhe3))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jhe3), aion(jhe3))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jli6), aion(jli6))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jli6), aion(jli6))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jli7), aion(jli7))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jbe7), aion(jbe7))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jbe9), aion(jbe9))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jb11), aion(jb11))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jc12), aion(jc12))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jc12), aion(jc12))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jc13), aion(jc13))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jc14), aion(jc14))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jc14), aion(jc14))

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

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jo18), aion(jo18))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jf17), aion(jf17))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jf17), aion(jf17))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jf18), aion(jf18))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jf18), aion(jf18))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jf19), aion(jf19))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jf19), aion(jf19))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jne19), aion(jne19))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jne20), aion(jne20))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jne20), aion(jne20))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jne21), aion(jne21))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jne21), aion(jne21))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jne22), aion(jne22))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jne22), aion(jne22))

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

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jal25), aion(jal25))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jal26), aion(jal26))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jal27), aion(jal27))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jal27), aion(jal27))

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
      zion(jsi31), aion(jsi31))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsi31), aion(jsi31))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsi32), aion(jsi32))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsi32), aion(jsi32))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jp29), aion(jp29))

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
      zion(jp33), aion(jp33))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jp33), aion(jp33))

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
      zion(js35), aion(js35))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(js35), aion(js35))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(js36), aion(js36))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(js36), aion(js36))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcl33), aion(jcl33))

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
      zion(jcl37), aion(jcl37))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcl37), aion(jcl37))

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
      zion(jar39), aion(jar39))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jar39), aion(jar39))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jar40), aion(jar40))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jar40), aion(jar40))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jk39), aion(jk39))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jk39), aion(jk39))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jk40), aion(jk40))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jk40), aion(jk40))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jk41), aion(jk41))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jk41), aion(jk41))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jca40), aion(jca40))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jca41), aion(jca41))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jca42), aion(jca42))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jca42), aion(jca42))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jca43), aion(jca43))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jca43), aion(jca43))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jca44), aion(jca44))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jca44), aion(jca44))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jca45), aion(jca45))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jca45), aion(jca45))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jca46), aion(jca46))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jca46), aion(jca46))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jca47), aion(jca47))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jca47), aion(jca47))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jca48), aion(jca48))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsc43), aion(jsc43))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsc43), aion(jsc43))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsc44), aion(jsc44))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsc44), aion(jsc44))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsc45), aion(jsc45))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsc45), aion(jsc45))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsc46), aion(jsc46))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsc46), aion(jsc46))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsc47), aion(jsc47))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsc47), aion(jsc47))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsc48), aion(jsc48))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsc48), aion(jsc48))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jsc49), aion(jsc49))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jti44), aion(jti44))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jti45), aion(jti45))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jti45), aion(jti45))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jti46), aion(jti46))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jti46), aion(jti46))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jti47), aion(jti47))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jti47), aion(jti47))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jti48), aion(jti48))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jti48), aion(jti48))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jti49), aion(jti49))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jti49), aion(jti49))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jti50), aion(jti50))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jti50), aion(jti50))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jti51), aion(jti51))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jv46), aion(jv46))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jv47), aion(jv47))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jv47), aion(jv47))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jv48), aion(jv48))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jv48), aion(jv48))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jv49), aion(jv49))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jv49), aion(jv49))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jv50), aion(jv50))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jv50), aion(jv50))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jv51), aion(jv51))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jv51), aion(jv51))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jv52), aion(jv52))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcr48), aion(jcr48))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcr49), aion(jcr49))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcr49), aion(jcr49))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcr50), aion(jcr50))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcr50), aion(jcr50))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcr51), aion(jcr51))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcr51), aion(jcr51))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcr52), aion(jcr52))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcr52), aion(jcr52))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcr53), aion(jcr53))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcr53), aion(jcr53))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcr54), aion(jcr54))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcr54), aion(jcr54))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmn50), aion(jmn50))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jmn51), aion(jmn51))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmn51), aion(jmn51))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jmn52), aion(jmn52))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmn52), aion(jmn52))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jmn53), aion(jmn53))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmn53), aion(jmn53))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jmn54), aion(jmn54))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmn54), aion(jmn54))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jmn55), aion(jmn55))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmn55), aion(jmn55))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jfe52), aion(jfe52))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jfe52), aion(jfe52))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jfe53), aion(jfe53))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jfe53), aion(jfe53))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jfe54), aion(jfe54))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jfe54), aion(jfe54))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jfe55), aion(jfe55))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jfe55), aion(jfe55))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jfe56), aion(jfe56))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jfe56), aion(jfe56))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jfe57), aion(jfe57))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jfe57), aion(jfe57))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jfe58), aion(jfe58))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jfe58), aion(jfe58))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jco53), aion(jco53))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jco54), aion(jco54))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jco55), aion(jco55))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jco55), aion(jco55))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jco56), aion(jco56))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jco56), aion(jco56))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jco57), aion(jco57))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jco57), aion(jco57))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jco58), aion(jco58))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jco58), aion(jco58))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jco59), aion(jco59))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jco59), aion(jco59))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jni56), aion(jni56))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jni56), aion(jni56))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jni57), aion(jni57))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jni57), aion(jni57))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jni58), aion(jni58))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jni58), aion(jni58))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jni59), aion(jni59))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jni59), aion(jni59))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jni60), aion(jni60))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jni60), aion(jni60))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jni61), aion(jni61))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jni61), aion(jni61))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jni62), aion(jni62))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jni62), aion(jni62))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jni63), aion(jni63))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jni64), aion(jni64))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcu58), aion(jcu58))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcu58), aion(jcu58))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcu59), aion(jcu59))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcu59), aion(jcu59))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcu60), aion(jcu60))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcu60), aion(jcu60))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcu61), aion(jcu61))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcu62), aion(jcu62))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcu63), aion(jcu63))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcu64), aion(jcu64))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jcu65), aion(jcu65))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jzn59), aion(jzn59))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jzn60), aion(jzn60))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jzn61), aion(jzn61))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jzn62), aion(jzn62))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jzn63), aion(jzn63))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jga62), aion(jga62))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jga63), aion(jga63))

    call add_screening_factor(zion(jd), aion(jd), &
      zion(jhe3), aion(jhe3))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jhe4), aion(jhe4))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jhe4), aion(jhe4))

    call add_screening_factor(zion(jd), aion(jd), &
      zion(jli6), aion(jli6))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jli7), aion(jli7))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jbe7), aion(jbe7))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jbe9), aion(jbe9))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jb10), aion(jb10))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jb10), aion(jb10))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jb11), aion(jb11))

    call add_screening_factor(zion(jc12), aion(jc12), &
      zion(jc12), aion(jc12))

    call add_screening_factor(zion(jd), aion(jd), &
      zion(jc13), aion(jc13))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jc13), aion(jc13))

    call add_screening_factor(zion(jd), aion(jd), &
      zion(jc14), aion(jc14))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jn13), aion(jn13))

    call add_screening_factor(zion(jc12), aion(jc12), &
      zion(jo16), aion(jo16))

    call add_screening_factor(zion(jo16), aion(jo16), &
      zion(jo16), aion(jo16))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jne18), aion(jne18))

    call add_screening_factor(zion(jc12), aion(jc12), &
      zion(jne20), aion(jne20))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jna21), aion(jna21))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jmg23), aion(jmg23))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jal26), aion(jal26))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jk37), aion(jk37))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jk38), aion(jk38))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jca40), aion(jca40))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jca41), aion(jca41))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jca48), aion(jca48))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jsc49), aion(jsc49))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jti51), aion(jti51))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jv52), aion(jv52))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jni63), aion(jni63))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcu57), aion(jcu57))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcu61), aion(jcu61))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcu62), aion(jcu62))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jcu63), aion(jcu63))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jzn60), aion(jzn60))

    call add_screening_factor(zion(jhe4), aion(jhe4), &
      zion(jzn61), aion(jzn61))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jzn64), aion(jzn64))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jzn65), aion(jzn65))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jzn66), aion(jzn66))

    call add_screening_factor(zion(jp), aion(jp), &
      zion(jga64), aion(jga64))

    call add_screening_factor(zion(jhe3), aion(jhe3), &
      zion(jhe3), aion(jhe3))

    call add_screening_factor(zion(jd), aion(jd), &
      zion(jli7), aion(jli7))

    call add_screening_factor(zion(jd), aion(jd), &
      zion(jbe7), aion(jbe7))

    call add_screening_factor(zion(jhe3), aion(jhe3), &
      zion(jli7), aion(jli7))

    call add_screening_factor(zion(jhe3), aion(jhe3), &
      zion(jbe7), aion(jbe7))


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
