module actual_eos_module

    use eos_type_module

    character (len=64), public :: eos_name = "helmholtz"

    ! Runtime parameters
    logical, allocatable :: do_coulomb
    logical, allocatable :: input_is_constant

    !..for the tables, in general
    integer, parameter, private :: imax = 541, jmax = 201
    integer, allocatable :: itmax, jtmax
    double precision, allocatable :: d(:), t(:)

    double precision, allocatable :: tlo, thi, tstp, tstpi
    double precision, allocatable :: dlo, dhi, dstp, dstpi

    double precision, allocatable :: ttol, dtol

    !..for the helmholtz free energy tables
    double precision, allocatable :: f(:,:,:)

    !..for the pressure derivative with density tables
    double precision, allocatable :: dpdf(:,:,:)

    !..for chemical potential tables
    double precision, allocatable :: ef(:,:,:)

    !..for the number density tables
    double precision, allocatable :: xf(:,:,:)

    !..for storing the differences
    double precision, allocatable :: dt_sav(:), dt2_sav(:),          &
                                     dti_sav(:), dt2i_sav(:),        &
                                     dd_sav(:), dd2_sav(:),          &
                                     ddi_sav(:), dd2i_sav(:)

#ifdef CUDA
    attributes(managed) :: do_coulomb, input_is_constant
    attributes(managed) :: itmax, jtmax
    attributes(managed) :: d, t
    attributes(managed) :: tlo, thi, tstp, tstpi
    attributes(managed) :: dlo, dhi, dstp, dstpi
    attributes(managed) :: ttol, dtol
    attributes(managed) :: f
    attributes(managed) :: dpdf
    attributes(managed) :: ef
    attributes(managed) :: xf
    attributes(managed) :: dt_sav, dt2_sav, dti_sav, dt2i_sav
    attributes(managed) :: dd_sav, dd2_sav, ddi_sav, dd2i_sav
#endif

    integer, parameter          :: max_newton = 100

    ! 2006 CODATA physical constants
private
    ! Math constants
    double precision, parameter :: pi       = 3.1415926535897932384d0

    ! Physical constants
    double precision, parameter :: h       = 6.6260689633d-27
    double precision, parameter :: hbar    = 0.5d0 * h/pi
    double precision, parameter :: qe      = 4.8032042712d-10
    double precision, parameter :: avo_eos = 6.0221417930d23
    double precision, parameter :: avo_eosi = 1.0d0 / avo_eos
    double precision, parameter :: clight  = 2.99792458d10
    double precision, parameter :: kerg    = 1.380650424d-16
    double precision, parameter :: kergi   = 1.0d0 / kerg
    double precision, parameter :: ev2erg_eos  = 1.60217648740d-12
    double precision, parameter :: kev     = kerg/ev2erg_eos
    double precision, parameter :: amu     = 1.66053878283d-24
    double precision, parameter :: me_eos  = 9.1093821545d-28
    double precision, parameter :: rbohr   = hbar*hbar/(me_eos * qe * qe)
    double precision, parameter :: fine    = qe*qe/(hbar*clight)

#ifdef RADIATION
    double precision, parameter :: ssol    = 0.0d0
#else
    double precision, parameter :: ssol    = 5.67051d-5
#endif
    double precision, parameter :: asol    = 4.0d0 * ssol / clight
    double precision, parameter :: weinlam = h*clight/(kerg * 4.965114232d0)
    double precision, parameter :: weinfre = 2.821439372d0*kerg/h

    ! Astronomical constants
    double precision, parameter :: ly      = 9.460528d17
    double precision, parameter :: pc      = 3.261633d0 * ly

    ! Some other useful combinations of the constants
    double precision, parameter :: sioncon = (2.0d0 * pi * amu * kerg)/(h*h)
    double precision, parameter :: forth   = 4.0d0/3.0d0
    double precision, parameter :: forpi   = 4.0d0 * pi
    double precision, parameter :: kergavo = kerg * avo_eos
    double precision, parameter :: ikavo   = 1.0d0/kergavo
    double precision, parameter :: asoli3  = asol/3.0d0
    double precision, parameter :: light2  = clight * clight

    ! Constants used for the Coulomb corrections
    double precision, parameter :: a1    = -0.898004d0
    double precision, parameter :: b1    =  0.96786d0
    double precision, parameter :: c1    =  0.220703d0
    double precision, parameter :: d1    = -0.86097d0
    double precision, parameter :: e1    =  2.5269d0
    double precision, parameter :: a2    =  0.29561d0
    double precision, parameter :: b2    =  1.9885d0
    double precision, parameter :: c2    =  0.288675d0
    double precision, parameter :: onethird = 1.0d0/3.0d0
    double precision, parameter :: esqu = qe * qe

    !$acc declare &
    !$acc create(tlo, thi, dlo, dhi) &
    !$acc create(tstp, tstpi, dstp, dstpi) &
    !$acc create(ttol, dtol) &
    !$acc create(itmax, jtmax, d, t) &
    !$acc create(f) &
    !$acc create(dpdf) &
    !$acc create(ef, xf)  &
    !$acc create(dt_sav, dt2_sav, dti_sav, dt2i_sav) &
    !$acc create(dd_sav, dd2_sav, ddi_sav, dd2i_sav) &
    !$acc create(do_coulomb, input_is_constant)

public actual_eos, actual_eos_init, actual_eos_finalize

contains

    !  Frank Timmes Helmholtz based Equation of State
    !  http://cococubed.asu.edu/

    !..given a temperature temp [K], density den [g/cm**3], and a composition
    !..characterized by abar and zbar, this routine returns most of the other
    !..thermodynamic quantities. of prime interest is the pressure [erg/cm**3],
    !..specific thermal energy [erg/gr], the entropy [erg/g/K], along with
    !..their derivatives with respect to temperature, density, abar, and zbar.
    !..other quantites such the normalized chemical potential eta (plus its
    !..derivatives), number density of electrons and positron pair (along
    !..with their derivatives), adiabatic indices, specific heats, and
    !..relativistically correct sound speed are also returned.
    !..
    !..this routine assumes planckian photons, an ideal gas of ions,
    !..and an electron-positron gas with an arbitrary degree of relativity
    !..and degeneracy. interpolation in a table of the helmholtz free energy
    !..is used to return the electron-positron thermodynamic quantities.
    !..all other derivatives are analytic.
    !..
    !..references: cox & giuli chapter 24 ; timmes & swesty apj 1999

    AMREX_DEVICE subroutine actual_eos(input, state)

        !$acc routine seq

        use amrex_constants_module, only: ZERO, HALF, TWO

        implicit none

        !..input arguments
        integer,      intent(in   ) :: input
        type (eos_t), intent(inout) :: state

        !..rows to store EOS data
        double precision :: temp_row, den_row

        !..declare local variables
        logical :: single_iter, double_iter, converged
        integer :: var, dvar, var1, var2, iter
        double precision :: v_want
        double precision :: v1_want, v2_want
        double precision :: xnew, xtol, dvdx, smallx, error, v
        double precision :: v1, v2, dv1dt, dv1dr, dv2dt,dv2dr, delr, error1, error2, told, rold, tnew, rnew, v1i, v2i

        double precision :: x,xi,y,zz,zzi,deni,tempi,presi,xni,dxnidd,dxnida, &
                            dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt, &
                            dpraddt,deraddd,deraddt,dpiondd,dpiondt, &
                            deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt, &
                            kt,ktinv,prad,erad,srad,pion,eion, &
                            sion,pele,eele,sele,pres,ener,entr,dpresdd, &
                            dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp,cvi, &
                            gam1,chit,chid,nabad,etaele, &
                            xnefer,s,si, &
                            temp,den,abar,zbar,ytot1,ye, &
                            enth,denthdd,denthdt,dpe,dpdr_e

#ifdef EXTRA_THERMO
        !..for the abar derivatives
        double precision :: dpionda,deionda, &
                            dpepda,deepda,dsepda,    &
                            dpresda,denerda


        !..for the zbar derivatives
        double precision :: dpepdz,deepdz,dsepdz,    &
                            dpresdz,denerdz
#endif


        !..for the interpolations
        integer          :: iat,jat
        double precision :: free,df_d,df_t,df_tt,df_dt
        double precision :: xt,xd,mxt,mxd,z,din,fi(36),fwtr(6),wdt(16)
        double precision :: sit(6),sid(6),dsit(6),dsid(6),ddsit(6)

        !..for the coulomb corrections
        double precision :: dsdd,dsda,lami,inv_lami,lamida,lamidd,     &
                            plasg,plasgi,plasgdd,plasgdt,plasgda,plasgdz,     &
                            ecoul,decouldd,decouldt,decoulda,decouldz, &
                            pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz, &
                            scoul,dscouldd,dscouldt,dscoulda,dscouldz

        double precision :: p_temp, e_temp

        double precision :: smallt, smalld

        !$gpu

        call eos_get_small_temp(smallt)
        call eos_get_small_dens(smalld)

        temp_row = state % T
        den_row  = state % rho
        abar     = state % abar
        zbar     = state % zbar
        ytot1    = 1.0d0 / abar
        ye       = state % y_e
        dxnidd   = avo_eos * ytot1

        ! Initial setup for iterations

        single_iter = .false.
        double_iter = .false.

        if (input .eq. eos_input_rt) then

          ! Nothing to do here.

        elseif (input .eq. eos_input_rh) then

          single_iter = .true.
          v_want = state % h
          var  = ienth
          dvar = itemp

        elseif (input .eq. eos_input_tp) then

          single_iter = .true.
          v_want = state % p
          var  = ipres
          dvar = idens

        elseif (input .eq. eos_input_rp) then

          single_iter = .true.
          v_want = state % p
          var  = ipres
          dvar = itemp

        elseif (input .eq. eos_input_re) then

          single_iter = .true.
          v_want = state % e
          var  = iener
          dvar = itemp

        elseif (input .eq. eos_input_ps) then

          double_iter = .true.
          v1_want = state % p
          v2_want = state % s
          var1 = ipres
          var2 = ientr

        elseif (input .eq. eos_input_ph) then

          double_iter = .true.
          v1_want = state % p
          v2_want = state % h
          var1 = ipres
          var2 = ienth

        elseif (input .eq. eos_input_th) then

          single_iter = .true.
          v_want = state % h
          var  = ienth
          dvar = idens

        endif

        converged = .false.

        if (input .eq. eos_input_rt) converged = .true.

        do iter = 1, max_newton

           temp  = temp_row
           den   =  den_row

           din   = ye * den

           !..initialize
           deni    = 1.0d0/den
           tempi   = 1.0d0/temp
           kt      = kerg * temp
           ktinv   = kergi * tempi

           !..radiation section:
           prad    = asoli3 * temp * temp * temp * temp
           dpraddt = 4.0d0 * prad*tempi

           erad    = 3.0d0 * prad*deni
           deraddd = -erad*deni
           deraddt = 3.0d0 * dpraddt*deni

           srad    = (prad*deni + erad)*tempi
           dsraddd = (-prad*deni*deni + deraddd)*tempi ! + tempi * dpraddd * deni = 0
           dsraddt = (dpraddt*deni + deraddt - srad)*tempi

           !..ion section:
           xni     = avo_eos * ytot1 * den
           dxnida  = -xni * ytot1

           pion    = xni * kt
           dpiondd = dxnidd * kt
           dpiondt = xni * kerg
#ifdef EXTRA_THERMO
           dpionda = dxnida * kt
#endif

           eion    = 1.5d0 * pion*deni
           deiondd = (1.5d0 * dpiondd - eion)*deni
           deiondt = 1.5d0 * dpiondt*deni
#ifdef EXTRA_THERMO
           deionda = 1.5d0 * dpionda*deni
#endif

           x       = abar*abar*sqrt(abar) * deni * avo_eosi
           s       = sioncon * temp
           z       = x * s * sqrt(s)
           y       = log(z)
           sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
           dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi &
                - kergavo * deni * ytot1
           dsiondt = (dpiondt*deni + deiondt)*tempi -  &
                (pion*deni + eion) * tempi*tempi  &
                + 1.5d0 * kergavo * tempi*ytot1

           !..electron-positron section:

           !..hash locate this temperature and density
           jat = int((log10(temp) - tlo)*tstpi) + 1
           jat = max(1,min(jat,jtmax-1))
           iat = int((log10(din) - dlo)*dstpi) + 1
           iat = max(1,min(iat,itmax-1))

           !..access the table locations only once
           fi(1:9)   = f(1:9, iat  ,jat  ) ! f, ft, ftt, fd, fdd, fdt, fddt, fdtt, fddtt
           fi(10:18) = f(1:9, iat+1,jat  )
           fi(19:27) = f(1:9, iat  ,jat+1)
           fi(28:36) = f(1:9, iat+1,jat+1)

           !..various differences
           xt  = max( (temp - t(jat))*dti_sav(jat), 0.0d0)
           xd  = max( (din - d(iat))*ddi_sav(iat), 0.0d0)
           mxt = 1.0d0 - xt
           mxd = 1.0d0 - xd

           !..the six density and six temperature basis functions
           sit(1) = psi0(xt)
           sit(2) = psi1(xt)*dt_sav(jat)
           sit(3) = psi2(xt)*dt2_sav(jat)

           sit(4) =  psi0(mxt)
           sit(5) = -psi1(mxt)*dt_sav(jat)
           sit(6) =  psi2(mxt)*dt2_sav(jat)

           sid(1) =   psi0(xd)
           sid(2) =   psi1(xd)*dd_sav(iat)
           sid(3) =   psi2(xd)*dd2_sav(iat)

           sid(4) =  psi0(mxd)
           sid(5) = -psi1(mxd)*dd_sav(iat)
           sid(6) =  psi2(mxd)*dd2_sav(iat)

           !..derivatives of the weight functions
           dsit(1) =   dpsi0(xt)*dti_sav(jat)
           dsit(2) =   dpsi1(xt)
           dsit(3) =   dpsi2(xt)*dt_sav(jat)

           dsit(4) = -dpsi0(mxt)*dti_sav(jat)
           dsit(5) =  dpsi1(mxt)
           dsit(6) = -dpsi2(mxt)*dt_sav(jat)

           dsid(1) =   dpsi0(xd)*ddi_sav(iat)
           dsid(2) =   dpsi1(xd)
           dsid(3) =   dpsi2(xd)*dd_sav(iat)

           dsid(4) = -dpsi0(mxd)*ddi_sav(iat)
           dsid(5) =  dpsi1(mxd)
           dsid(6) = -dpsi2(mxd)*dd_sav(iat)

           !..second derivatives of the weight functions
           ddsit(1) =   ddpsi0(xt)*dt2i_sav(jat)
           ddsit(2) =   ddpsi1(xt)*dti_sav(jat)
           ddsit(3) =   ddpsi2(xt)

           ddsit(4) =  ddpsi0(mxt)*dt2i_sav(jat)
           ddsit(5) = -ddpsi1(mxt)*dti_sav(jat)
           ddsit(6) =  ddpsi2(mxt)

           ! This array saves some subexpressions that go into
           ! computing the biquintic polynomial. Instead of explicitly
           ! constructing it in full, we'll use these subexpressions
           ! and then compute the result as
           ! (table data * temperature terms) * density terms.

           fwtr(1:6) = fwt(fi, sit)

           !..the free energy
           free  = sum(fwtr * sid)

           !..derivative with respect to density
           df_d  = sum(fwtr * dsid)

           fwtr(1:6) = fwt(fi, dsit)

           !..derivative with respect to temperature
           df_t  = sum(fwtr * sid)

           !..derivative with respect to temperature and density
           df_dt = sum(fwtr * dsid)

           fwtr(1:6) = fwt(fi, ddsit)

           !..derivative with respect to temperature**2
           df_tt = sum(fwtr * sid)

           !..now get the pressure derivative with density, chemical potential, and
           !..electron positron number densities
           !..get the interpolation weight functions
           sit(1) = xpsi0(xt)
           sit(2) = xpsi1(xt)*dt_sav(jat)

           sit(3) = xpsi0(mxt)
           sit(4) = -xpsi1(mxt)*dt_sav(jat)

           sid(1) = xpsi0(xd)
           sid(2) = xpsi1(xd)*dd_sav(iat)

           sid(3) = xpsi0(mxd)
           sid(4) = -xpsi1(mxd)*dd_sav(iat)

           !..derivatives of weight functions
           dsit(1) = xdpsi0(xt)*dti_sav(jat)
           dsit(2) = xdpsi1(xt)

           dsit(3) = -xdpsi0(mxt)*dti_sav(jat)
           dsit(4) = xdpsi1(mxt)

           dsid(1) = xdpsi0(xd)*ddi_sav(iat)
           dsid(2) = xdpsi1(xd)

           dsid(3) = -xdpsi0(mxd)*ddi_sav(iat)
           dsid(4) = xdpsi1(mxd)

           ! Reuse subexpressions that would go into computing h3.
           wdt(1:4)   = sid(1) * sit(1:4)
           wdt(5:8)   = sid(2) * sit(1:4)
           wdt(9:12)  = sid(3) * sit(1:4)
           wdt(13:16) = sid(4) * sit(1:4)

           !..look in the pressure derivative only once
           fi(1:4)   = dpdf(1:4, iat  ,jat  )
           fi(5:8)   = dpdf(1:4, iat+1,jat  )
           fi(9:12)  = dpdf(1:4, iat  ,jat+1)
           fi(13:16) = dpdf(1:4, iat+1,jat+1)

           !..pressure derivative with density
           dpepdd  = h3(fi, wdt)
           dpepdd  = max(ye * dpepdd,0.0d0)

           !..look in the electron chemical potential table only once
           fi(1:4)   = ef(1:4,iat  ,jat  )
           fi(5:8)   = ef(1:4,iat+1,jat  )
           fi(9:12)  = ef(1:4,iat  ,jat+1)
           fi(13:16) = ef(1:4,iat+1,jat+1)

           !..electron chemical potential etaele
           etaele  = h3(fi, wdt)

           !..look in the number density table only once
           fi(1:4)   = xf(1:4,iat  ,jat  )
           fi(5:8)   = xf(1:4,iat+1,jat  )
           fi(9:12)  = xf(1:4,iat  ,jat+1)
           fi(13:16) = xf(1:4,iat+1,jat+1)

           !..electron + positron number densities
           xnefer   = h3(fi, wdt)

           wdt(1:4)   = dsid(1) * sit(1:4)
           wdt(5:8)   = dsid(2) * sit(1:4)
           wdt(9:12)  = dsid(3) * sit(1:4)
           wdt(13:16) = dsid(4) * sit(1:4)

           !..derivative with respect to density
           x = h3(fi, wdt)
           x = max(x,0.0d0)

           !..the desired electron-positron thermodynamic quantities
           x       = din * din
           pele    = x * df_d
           dpepdt  = x * df_dt
           s       = dpepdd/ye - 2.0d0 * din * df_d
#ifdef EXTRA_THERMO
           dpepda  = -ytot1 * (2.0d0 * pele + s * din)
           dpepdz  = den*ytot1*(2.0d0 * din * df_d  +  s)
#endif

           x       = ye * ye
           sele    = -df_t * ye
           dsepdt  = -df_tt * ye
           dsepdd  = -df_dt * x
#ifdef EXTRA_THERMO
           dsepda  = ytot1 * (ye * df_dt * din - sele)
           dsepdz  = -ytot1 * (ye * df_dt * den  + df_t)
#endif

           eele    = ye*free + temp * sele
           deepdt  = temp * dsepdt
           deepdd  = x * df_d + temp * dsepdd
#ifdef EXTRA_THERMO
           deepda  = -ye * ytot1 * (free +  df_d * din) + temp * dsepda
           deepdz  = ytot1* (free + ye * df_d * den) + temp * dsepdz
#endif

           !..coulomb section:
           !..initialize
           pcoul    = 0.0d0
           dpcouldd = 0.0d0
           dpcouldt = 0.0d0
           dpcoulda = 0.0d0
           dpcouldz = 0.0d0
           ecoul    = 0.0d0
           decouldd = 0.0d0
           decouldt = 0.0d0
           decoulda = 0.0d0
           decouldz = 0.0d0
           scoul    = 0.0d0
           dscouldd = 0.0d0
           dscouldt = 0.0d0
           dscoulda = 0.0d0
           dscouldz = 0.0d0

           if ( do_coulomb ) then

              !..uniform background corrections only
              !..from yakovlev & shalybkov 1989
              !..lami is the average ion seperation
              !..plasg is the plasma coupling parameter
              z        = forth * pi
              s        = z * xni
              si       = 1.0d0 / s
              dsdd     = z * dxnidd
              dsda     = z * dxnida

              lami     = 1.0d0 * si**onethird
              inv_lami = 1.0d0/lami
              z        = -onethird * lami
              lamidd   = z * dsdd * si
              lamida   = z * dsda * si

              plasg    = zbar*zbar*esqu*ktinv*inv_lami
              plasgi   = 1.0d0 / plasg
              z        = -plasg * inv_lami
              plasgdd  = z * lamidd
              plasgda  = z * lamida
              plasgdt  = -plasg * tempi
              plasgdz  = 2.0d0 * plasg/zbar

              !...yakovlev & shalybkov 1989 equations 82, 85, 86, 87
              if (plasg .ge. 1.0D0) then

                 x        = plasg**(0.25d0)
                 xi       = 1.0d0 / x
                 y        = avo_eos * ytot1 * kerg
                 ecoul    = y * temp * (a1 * plasg + b1 * x + c1 * xi + d1)
                 pcoul    = onethird * den * ecoul
                 scoul    = -y * (3.0d0 * b1 * x - 5.0d0 * c1 * xi &
                      + d1 * (log(plasg) - 1.0d0) - e1)

                 y        = avo_eos*ytot1*kt*(a1 + 0.25d0 * plasgi * (b1*x - c1 * xi))
                 decouldd = y * plasgdd
                 decouldt = y * plasgdt + ecoul * tempi
                 decoulda = y * plasgda - ecoul * ytot1
                 decouldz = y * plasgdz

                 y        = onethird * den
                 dpcouldd = onethird * ecoul + y*decouldd
                 dpcouldt = y * decouldt
                 dpcoulda = y * decoulda
                 dpcouldz = y * decouldz

                 y        = -avo_eos * kerg * ytot1 * plasgi * &
                      (0.75d0 * b1 * x + 1.25d0 *c1 * xi + d1)
                 dscouldd = y * plasgdd
                 dscouldt = y * plasgdt
                 dscoulda = y * plasgda - scoul * ytot1
                 dscouldz = y * plasgdz

              !...yakovlev & shalybkov 1989 equations 102, 103, 104
              else if (plasg .lt. 1.0D0) then

                 x        = plasg*sqrt(plasg)
                 y        = plasg**b2
                 z        = c2 * x - onethird * a2 * y
                 pcoul    = -pion * z
                 ecoul    = 3.0d0 * pcoul * deni
                 scoul    = -avo_eos * ytot1 * kerg * (c2 * x -a2 * (b2-1.0d0) / b2 * y)

                 s        = 1.5d0 * c2 * x * plasgi - onethird * a2 * b2 * y * plasgi
                 dpcouldd = -dpiondd*z - pion*s*plasgdd
                 dpcouldt = -dpiondt*z - pion*s*plasgdt
#ifdef EXTRA_THERMO
                 dpcoulda = -dpionda*z - pion*s*plasgda
                 dpcouldz = - pion*s*plasgdz ! -dpiondz*z = 0
#endif

                 s        = 3.0d0 * deni
                 decouldd = s * dpcouldd - ecoul * deni
                 decouldt = s * dpcouldt
                 decoulda = s * dpcoulda
                 decouldz = s * dpcouldz

                 s        = -avo_eos * kerg * ytot1 * plasgi * &
                      (1.5d0*c2*x-a2*(b2-1.0d0)*y)
                 dscouldd = s * plasgdd
                 dscouldt = s * plasgdt
                 dscoulda = s * plasgda - scoul * ytot1
                 dscouldz = s * plasgdz

              end if

              ! Disable Coulomb corrections if they cause
              ! the energy or pressure to go negative.

              p_temp = prad + pion + pele + pcoul
              e_temp = erad + eion + eele + ecoul

              if (p_temp .le. ZERO .or. e_temp .le. ZERO) then

                 pcoul    = 0.0d0
                 dpcouldd = 0.0d0
                 dpcouldt = 0.0d0
                 dpcoulda = 0.0d0
                 dpcouldz = 0.0d0
                 ecoul    = 0.0d0
                 decouldd = 0.0d0
                 decouldt = 0.0d0
                 decoulda = 0.0d0
                 decouldz = 0.0d0
                 scoul    = 0.0d0
                 dscouldd = 0.0d0
                 dscouldt = 0.0d0
                 dscoulda = 0.0d0
                 dscouldz = 0.0d0

              end if
           end if

           !..sum all the components
           pres    = prad + pion + pele + pcoul
           ener    = erad + eion + eele + ecoul
           entr    = srad + sion + sele + scoul

           dpresdd = dpiondd + dpepdd + dpcouldd ! + dpraddd = 0
           dpresdt = dpraddt + dpiondt + dpepdt + dpcouldt
#ifdef EXTRA_THERMO
           dpresda = dpionda + dpepda + dpcoulda ! + dpradda = 0
           dpresdz = dpepdz + dpcouldz ! + dpraddz + dpiondz = 0
#endif
           denerdd = deraddd + deiondd + deepdd + decouldd
           denerdt = deraddt + deiondt + deepdt + decouldt
#ifdef EXTRA_THERMO
           denerda = deionda + deepda + decoulda ! + deradda = 0
           denerdz = deepdz + decouldz ! + deraddz + deiondz = 0
#endif

           dentrdd = dsraddd + dsiondd + dsepdd + dscouldd
           dentrdt = dsraddt + dsiondt + dsepdt + dscouldt

           presi = 1.0d0 / pres

           !..the temperature and density exponents (c&g 9.81 9.82)
           !..the specific heat at constant volume (c&g 9.92)
           !..the third adiabatic exponent (c&g 9.93)
           !..the specific heat at constant pressure (c&g 9.98)
           zz    = pres*deni
           zzi   = den * presi
           chit  = temp * presi * dpresdt
           chid  = dpresdd*zzi
           cv    = denerdt
           cvi   = 1.0d0 / cv
           x     = zz * chit * tempi * cvi
           gam1  = chit*x + chid
           nabad = x/gam1
           cp    = cv * gam1/chid
           z     = 1.0d0 + (ener + light2)*zzi

           dpe = dpresdt * cvi
           dpdr_e = dpresdd - dpresdt * denerdd * cvi

           enth = ener + pres * deni
           denthdd = denerdd + dpresdd * deni - pres * deni**2
           denthdt = denerdt + dpresdt * deni

           if (converged) then

              exit

           elseif (single_iter) then

              if (dvar .eq. itemp) then

                 x = temp_row
                 smallx = smallt
                 xtol = ttol

                 if (var .eq. ipres) then
                    v    = pres
                    dvdx = dpresdt
                 elseif (var .eq. iener) then
                    v    = ener
                    dvdx = denerdt
                 elseif (var .eq. ientr) then
                    v    = entr
                    dvdx = dentrdt
                 elseif (var .eq. ienth) then
                    v    = enth
                    dvdx = denthdt
                 else
                    exit
                 endif

              else ! dvar == density

                 x = den_row
                 smallx = smalld
                 xtol = dtol

                 if (var .eq. ipres) then
                    v    = pres
                    dvdx = dpresdd
                 elseif (var .eq. iener) then
                    v    = ener
                    dvdx = denerdd
                 elseif (var .eq. ientr) then
                    v    = entr
                    dvdx = dentrdd
                 elseif (var .eq. ienth) then
                    v    = enth
                    dvdx = denthdd
                 else
                    exit
                 endif

              endif

              ! Now do the calculation for the next guess for T/rho

              xnew = x - (v - v_want) / dvdx

              ! Don't let the temperature/density change by more than a factor of two
              xnew = max(0.5 * x, min(xnew, 2.0 * x))

              ! Don't let us freeze/evacuate
              xnew = max(smallx, xnew)

              ! Store the new temperature/density

              if (dvar .eq. itemp) then
                 temp_row = xnew
              else
                 den_row  = xnew
              endif

              ! Compute the error from the last iteration

              error = abs( (xnew - x) / x )

              if (error .lt. xtol) converged = .true.

           elseif (double_iter) then

              ! Figure out which variables we're using

              told = temp_row
              rold = den_row

              if (var1 .eq. ipres) then
                 v1    = pres
                 dv1dt = dpresdt
                 dv1dr = dpresdd
              elseif (var1 .eq. iener) then
                 v1    = ener
                 dv1dt = denerdt
                 dv1dr = denerdd
              elseif (var1 .eq. ientr) then
                 v1    = entr
                 dv1dt = dentrdt
                 dv1dr = dentrdd
              elseif (var1 .eq. ienth) then
                 v1    = enth
                 dv1dt = denthdt
                 dv1dr = denthdd
              else
                 exit
              endif

              if (var2 .eq. ipres) then
                 v2    = pres
                 dv2dt = dpresdt
                 dv2dr = dpresdd
              elseif (var2 .eq. iener) then
                 v2    = ener
                 dv2dt = denerdt
                 dv2dr = denerdd
              elseif (var2 .eq. ientr) then
                 v2    = entr
                 dv2dt = dentrdt
                 dv2dr = dentrdd
              elseif (var2 .eq. ienth) then
                 v2    = enth
                 dv2dt = denthdt
                 dv2dr = denthdd
              else
                 exit
              endif

              ! Two functions, f and g, to iterate over
              v1i = v1_want - v1
              v2i = v2_want - v2

              !
              ! 0 = f + dfdr * delr + dfdt * delt
              ! 0 = g + dgdr * delr + dgdt * delt
              !

              ! note that dfi/dT = - df/dT
              delr = (-v1i*dv2dt + v2i*dv1dt) / (dv2dr*dv1dt - dv2dt*dv1dr)

              rnew = rold + delr

              tnew = told + (v1i - dv1dr*delr) / dv1dt

              ! Don't let the temperature or density change by more
              ! than a factor of two
              tnew = max(HALF * told, min(tnew, TWO * told))
              rnew = max(HALF * rold, min(rnew, TWO * rold))

              ! Don't let us freeze or evacuate
              tnew = max(smallt, tnew)
              rnew = max(smalld, rnew)

              ! Store the new temperature and density
              den_row  = rnew
              temp_row = tnew

              ! Compute the errors
              error1 = abs( (rnew - rold) / rold )
              error2 = abs( (tnew - told) / told )

              if (error1 .LT. dtol .and. error2 .LT. ttol) converged = .true.

           endif

        enddo

        state % T    = temp_row
        state % rho  = den_row

        state % p    = pres
        state % dpdT = dpresdt
        state % dpdr = dpresdd

#ifdef EXTRA_THERMO
        state % dpdA = dpresda
        state % dpdZ = dpresdz
#endif

        state % dpde = dpe
        state % dpdr_e = dpdr_e

        state % e    = ener
        state % dedT = denerdt
        state % dedr = denerdd

#ifdef EXTRA_THERMO
        state % dedA = denerda
        state % dedZ = denerdz
#endif

        state % s    = entr
        state % dsdT = dentrdt
        state % dsdr = dentrdd

        state % h    = enth
        state % dhdR = denthdd
        state % dhdT = denthdt

        state % pele = pele
        state % ppos = 0.0d0

        state % xne = xnefer
        state % xnp = 0.0d0

        state % eta = etaele

        state % cv   = cv
        state % cp   = cp
        state % gam1 = gam1

        ! Take care of final housekeeping.

        ! Count the positron contribution in the electron quantities.

        state % xne  = state % xne  + state % xnp
        state % pele = state % pele + state % ppos

        ! Use the non-relativistic version of the sound speed, cs = sqrt(gam_1 * P / rho).
        ! This replaces the relativistic version that comes out of helmeos.

        state % cs = sqrt(state % gam1 * state % p / state % rho)

        if (input_is_constant) then

          if (input .eq. eos_input_rh) then

            state % h = v_want

          elseif (input .eq. eos_input_tp) then

            state % p = v_want

          elseif (input .eq. eos_input_rp) then

            state % p = v_want

          elseif (input .eq. eos_input_re) then

            state % e = v_want

          elseif (input .eq. eos_input_ps) then

            state % p = v1_want
            state % s = v2_want

          elseif (input .eq. eos_input_ph) then

            state % p = v1_want
            state % h = v2_want

          elseif (input .eq. eos_input_th) then

            state % h = v_want

          endif

        endif

    end subroutine actual_eos



    subroutine actual_eos_init

        use amrex_error_module
        use extern_probin_module, only: eos_input_is_constant, use_eos_coulomb, eos_ttol, eos_dtol
        use amrex_paralleldescriptor_module, only: parallel_bcast => amrex_pd_bcast, amrex_pd_ioprocessor

        implicit none

        double precision :: dth, dt2, dti, dt2i
        double precision :: dd, dd2, ddi, dd2i
        double precision :: tsav, dsav
        integer :: i, j
        integer :: status

        ! Allocate managed module variables

        allocate(do_coulomb)
        allocate(input_is_constant)
        allocate(itmax)
        allocate(jtmax)
        allocate(d(imax))
        allocate(t(jmax))
        allocate(tlo)
        allocate(thi)
        allocate(tstp)
        allocate(tstpi)
        allocate(dlo)
        allocate(dhi)
        allocate(dstp)
        allocate(dstpi)
        allocate(ttol)
        allocate(dtol)
        allocate(f(9,imax,jmax))
        allocate(dpdf(4,imax,jmax))
        allocate(ef(4,imax,jmax))
        allocate(xf(4,imax,jmax))
        allocate(dt_sav(jmax))
        allocate(dt2_sav(jmax))
        allocate(dti_sav(jmax))
        allocate(dt2i_sav(jmax))
        allocate(dd_sav(imax))
        allocate(dd2_sav(imax))
        allocate(ddi_sav(imax))
        allocate(dd2i_sav(imax))

        ! Read in the runtime parameters

        input_is_constant = eos_input_is_constant
        do_coulomb = use_eos_coulomb
        ttol = eos_ttol
        dtol = eos_dtol

        if (amrex_pd_ioprocessor()) then
           print *, ''
           if (do_coulomb) then
              print *, "Initializing Helmholtz EOS and using Coulomb corrections."
           else
              print *, "Initializing Helmholtz EOS without using Coulomb corrections."
           endif
           print *, ''
        endif

        !..   read the helmholtz free energy table
        itmax = imax
        jtmax = jmax
        tlo   = 3.0d0
        thi   = 13.0d0
        tstp  = (thi - tlo)/float(jmax-1)
        tstpi = 1.0d0/tstp
        dlo   = -12.0d0
        dhi   = 15.0d0
        dstp  = (dhi - dlo)/float(imax-1)
        dstpi = 1.0d0/dstp

        do j=1,jmax
           tsav = tlo + (j-1)*tstp
           t(j) = 10.0d0**(tsav)
           do i=1,imax
              dsav = dlo + (i-1)*dstp
              d(i) = 10.0d0**(dsav)
           end do
        end do

        if (amrex_pd_ioprocessor()) then

           !..   open the table
           open(unit=2,file='helm_table.dat',status='old',iostat=status,action='read')
           if (status > 0) then

              call amrex_error('actual_eos_init: Failed to open helm_table.dat')

           endif

           ! Note that in the below, the indices are read in slightly out of numerical
           ! order. This is so that they match up with how they are actually used in
           ! the calculation of h5 and h3.

           !...  read in the free energy table
           do j=1,jmax
              do i=1,imax
                 read(2,*) f(1,i,j),f(4,i,j),f(2,i,j),f(5,i,j),f(3,i,j),f(6,i,j), &
                      f(7,i,j),f(8,i,j),f(9,i,j)
              end do
           end do

           !..   read the pressure derivative with density table
           do j = 1, jmax
              do i = 1, imax
                 read(2,*) dpdf(1,i,j),dpdf(3,i,j),dpdf(2,i,j),dpdf(4,i,j)
              end do
           end do

           !..   read the electron chemical potential table
           do j = 1, jmax
              do i = 1, imax
                 read(2,*) ef(1,i,j),ef(3,i,j),ef(2,i,j),ef(4,i,j)
              end do
           end do

           !..   read the number density table
           do j = 1, jmax
              do i = 1, imax
                 read(2,*) xf(1,i,j),xf(3,i,j),xf(2,i,j),xf(4,i,j)
              end do
           end do

        end if

        call parallel_bcast(f)
        call parallel_bcast(dpdf)
        call parallel_bcast(ef)
        call parallel_bcast(xf)

        !..   construct the temperature and density deltas and their inverses
        do j = 1, jmax-1
           dth         = t(j+1) - t(j)
           dt2         = dth * dth
           dti         = 1.0d0/dth
           dt2i        = 1.0d0/dt2
           dt_sav(j)   = dth
           dt2_sav(j)  = dt2
           dti_sav(j)  = dti
           dt2i_sav(j) = dt2i
        end do
        do i = 1, imax-1
           dd          = d(i+1) - d(i)
           dd2         = dd * dd
           ddi         = 1.0d0/dd
           dd2i        = 1.0d0/dd2
           dd_sav(i)   = dd
           dd2_sav(i)  = dd2
           ddi_sav(i)  = ddi
           dd2i_sav(i) = dd2i
        end do

        if (amrex_pd_ioprocessor()) then
           close(unit=2)
        endif

        ! Set up the minimum and maximum possible densities.

        mintemp = 10.d0**tlo
        maxtemp = 10.d0**thi
        mindens = 10.d0**dlo
        maxdens = 10.d0**dhi

        !$acc update device(mintemp, maxtemp, mindens, maxdens)

        !$acc update &
        !$acc device(tlo, thi, dlo, dhi) &
        !$acc device(tstp, tstpi, dstp, dstpi) &
        !$acc device(itmax, jtmax, d, t) &
        !$acc device(f) &
        !$acc device(dpdf) &
        !$acc device(ef, xf)  &
        !$acc device(dt_sav, dt2_sav, dti_sav, dt2i_sav) &
        !$acc device(dd_sav, dd2_sav, ddi_sav, dd2i_sav) &
        !$acc device(do_coulomb, input_is_constant)

    end subroutine actual_eos_init



    ! quintic hermite polynomial functions
    ! psi0 and its derivatives
    AMREX_DEVICE pure function psi0(z) result(psi0r)
    !$acc routine seq
      double precision, intent(in) :: z
      double precision :: psi0r
      !$gpu
      psi0r = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
    end function psi0

    AMREX_DEVICE pure function dpsi0(z) result(dpsi0r) 
    !$acc routine seq
      double precision, intent(in) :: z
      double precision :: dpsi0r
      !$gpu
      dpsi0r = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
    end function dpsi0

    AMREX_DEVICE pure function ddpsi0(z) result(ddpsi0r)
    !$acc routine seq
      double precision, intent(in) :: z
      double precision :: ddpsi0r
      !$gpu
      ddpsi0r = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)
    end function ddpsi0

    ! psi1 and its derivatives
    AMREX_DEVICE pure function psi1(z) result(psi1r)
    !$acc routine seq
      double precision, intent(in) :: z
      double precision :: psi1r
      !$gpu
      psi1r = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
    end function psi1

    AMREX_DEVICE pure function dpsi1(z) result(dpsi1r)
    !$acc routine seq
      double precision, intent(in) :: z
      double precision :: dpsi1r
      !$gpu
      dpsi1r = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
    end function dpsi1

    AMREX_DEVICE pure function ddpsi1(z) result(ddpsi1r)
    !$acc routine seq
      double precision, intent(in) :: z
      double precision :: ddpsi1r
      !$gpu
      ddpsi1r = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)
    end function ddpsi1

    ! psi2  and its derivatives
    AMREX_DEVICE pure function psi2(z) result(psi2r)
      !$acc routine seq
      double precision, intent(in) :: z
      double precision :: psi2r
      !$gpu
      psi2r = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
    end function psi2

    AMREX_DEVICE pure function dpsi2(z) result(dpsi2r)
      !$acc routine seq
      double precision, intent(in) :: z
      double precision :: dpsi2r
      !$gpu
      dpsi2r = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
    end function dpsi2

    AMREX_DEVICE pure function ddpsi2(z) result(ddpsi2r)
      !$acc routine seq
      double precision, intent(in) :: z
      double precision :: ddpsi2r
      !$gpu
      ddpsi2r = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)
    end function ddpsi2

    AMREX_DEVICE pure function fwt(fi, wt) result(fwtr)
      !$acc routine seq
      double precision, intent(in) :: fi(36), wt(6)
      double precision :: fwtr(6)

      !$gpu

      fwtr(1) = fi( 1)*wt(1) + fi( 2)*wt(2) + fi( 3)*wt(3) + fi(19)*wt(4) + fi(20)*wt(5) + fi(21)*wt(6)
      fwtr(2) = fi( 4)*wt(1) + fi( 6)*wt(2) + fi( 8)*wt(3) + fi(22)*wt(4) + fi(24)*wt(5) + fi(26)*wt(6)
      fwtr(3) = fi( 5)*wt(1) + fi( 7)*wt(2) + fi( 9)*wt(3) + fi(23)*wt(4) + fi(25)*wt(5) + fi(27)*wt(6)
      fwtr(4) = fi(10)*wt(1) + fi(11)*wt(2) + fi(12)*wt(3) + fi(28)*wt(4) + fi(29)*wt(5) + fi(30)*wt(6)
      fwtr(5) = fi(13)*wt(1) + fi(15)*wt(2) + fi(17)*wt(3) + fi(31)*wt(4) + fi(33)*wt(5) + fi(35)*wt(6)
      fwtr(6) = fi(14)*wt(1) + fi(16)*wt(2) + fi(18)*wt(3) + fi(32)*wt(4) + fi(34)*wt(5) + fi(36)*wt(6)

    end function fwt


    ! cubic hermite polynomial functions
    ! psi0 & derivatives
    AMREX_DEVICE pure function xpsi0(z) result(xpsi0r)
      !$acc routine seq
      double precision, intent(in) :: z
      double precision :: xpsi0r
      !$gpu
      xpsi0r = z * z * (2.0d0*z - 3.0d0) + 1.0
    end function xpsi0

    AMREX_DEVICE pure function xdpsi0(z) result(xdpsi0r)
      !$acc routine seq
      double precision, intent(in) :: z
      double precision :: xdpsi0r
      !$gpu
      xdpsi0r = z * (6.0d0*z - 6.0d0)
    end function xdpsi0


    ! psi1 & derivatives
    AMREX_DEVICE pure function xpsi1(z) result(xpsi1r)
      !$acc routine seq
      double precision, intent(in) :: z
      double precision :: xpsi1r
      !$gpu
      xpsi1r = z * ( z * (z - 2.0d0) + 1.0d0)
    end function xpsi1

    AMREX_DEVICE pure function xdpsi1(z) result(xdpsi1r)
      !$acc routine seq
      double precision, intent(in) :: z
      double precision :: xdpsi1r
      !$gpu
      xdpsi1r = z * (3.0d0*z - 4.0d0) + 1.0d0
    end function xdpsi1

    ! bicubic hermite polynomial function
    AMREX_DEVICE pure function h3(fi, wdt) result(h3r)
      !$acc routine seq
      double precision, intent(in) :: fi(16), wdt(16)
      double precision :: h3r

      !$gpu

      h3r = fi( 1)*wdt( 1) + fi( 2)*wdt( 2) + fi( 9)*wdt( 3) + fi(10)*wdt( 4) + &
            fi( 3)*wdt( 5) + fi( 4)*wdt( 6) + fi(11)*wdt( 7) + fi(12)*wdt( 8) + &
            fi( 5)*wdt( 9) + fi( 6)*wdt(10) + fi(13)*wdt(11) + fi(14)*wdt(12) + &
            fi( 7)*wdt(13) + fi( 8)*wdt(14) + fi(15)*wdt(15) + fi(16)*wdt(16)

    end function h3

    subroutine actual_eos_finalize

      implicit none

      ! Deallocate managed module variables

      deallocate(do_coulomb)
      deallocate(input_is_constant)
      deallocate(itmax)
      deallocate(jtmax)
      deallocate(d)
      deallocate(t)
      deallocate(tlo)
      deallocate(thi)
      deallocate(tstp)
      deallocate(tstpi)
      deallocate(dlo)
      deallocate(dhi)
      deallocate(dstp)
      deallocate(dstpi)
      deallocate(ttol)
      deallocate(dtol)
      deallocate(f)
      deallocate(dpdf)
      deallocate(ef)
      deallocate(xf)
      deallocate(dt_sav)
      deallocate(dt2_sav)
      deallocate(dti_sav)
      deallocate(dt2i_sav)
      deallocate(dd_sav)
      deallocate(dd2_sav)
      deallocate(ddi_sav)
      deallocate(dd2i_sav)

    end subroutine actual_eos_finalize

end module actual_eos_module
