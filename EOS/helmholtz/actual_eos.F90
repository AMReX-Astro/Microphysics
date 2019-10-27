module actual_eos_module

    use helmholtz_constants_module
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
    double precision, allocatable :: f(:,:), fd(:,:),                &
                                     ft(:,:), fdd(:,:), ftt(:,:),    &
                                     fdt(:,:), fddt(:,:), fdtt(:,:), &
                                     fddtt(:,:)

    !..for the pressure derivative with density ables
    double precision, allocatable :: dpdf(:,:), dpdfd(:,:),          &
                                     dpdft(:,:), dpdfdt(:,:)

    !..for chemical potential tables
    double precision, allocatable :: ef(:,:), efd(:,:),              &
                                     eft(:,:), efdt(:,:)

    !..for the number density tables
    double precision, allocatable :: xf(:,:), xfd(:,:),              &
                                     xft(:,:), xfdt(:,:)

    !..for storing the differences
    double precision, allocatable :: dt_sav(:), dt2_sav(:),          &
                                     dti_sav(:), dt2i_sav(:),        &
                                     dd_sav(:), dd2_sav(:),          &
                                     ddi_sav(:), dd2i_sav(:)

#ifdef AMREX_USE_CUDA
    attributes(managed) :: do_coulomb, input_is_constant
    attributes(managed) :: itmax, jtmax
    attributes(managed) :: d, t
    attributes(managed) :: tlo, thi, tstp, tstpi
    attributes(managed) :: dlo, dhi, dstp, dstpi
    attributes(managed) :: ttol, dtol
    attributes(managed) :: f, fd, ft, fdd, ftt, fdt, fddt, fdtt, fddtt
    attributes(managed) :: dpdf, dpdfd, dpdft, dpdfdt
    attributes(managed) :: ef, efd, eft, efdt
    attributes(managed) :: xf, xfd, xft, xfdt
    attributes(managed) :: dt_sav, dt2_sav, dti_sav, dt2i_sav
    attributes(managed) :: dd_sav, dd2_sav, ddi_sav, dd2i_sav
#endif

    integer, parameter          :: max_newton = 100

    !$acc declare &
    !$acc create(tlo, thi, dlo, dhi) &
    !$acc create(tstp, tstpi, dstp, dstpi) &
    !$acc create(ttol, dtol) &
    !$acc create(itmax, jtmax, d, t) &
    !$acc create(f, fd, ft, fdd, ftt, fdt, fddt, fdtt, fddtt) &
    !$acc create(dpdf, dpdfd, dpdft, dpdfdt) &
    !$acc create(ef, efd, eft, efdt, xf, xfd, xft, xfdt)  &
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

    subroutine actual_eos(input, state)

        !$acc routine seq
      
        use amrex_constants_module, only: ZERO, HALF, TWO
        use coulomb_module, only: apply_coulomb_corrections
        use helmholtz_radiation_module, only: apply_radiation
        use helmholtz_ions_module, only: apply_ions

        implicit none

        !..input arguments
        integer,      intent(in   ) :: input
        type (eos_t), intent(inout) :: state

        !..rows to store EOS data
        double precision :: temp_row, &
                            den_row, &
                            abar_row, &
                            zbar_row, &
                            ye_row, &
                            etot_row, &
                            ptot_row, &
                            cv_row, &
                            cp_row,  &
                            xne_row, &
                            xnp_row, &
                            etaele_row, &
                            pele_row, &
                            ppos_row, &
                            dpd_row,  &
                            dpt_row, &
                            dpa_row, &
                            dpz_row,  &
                            ded_row, &
                            det_row, &
                            dea_row,  &
                            dez_row,  &
                            stot_row, &
                            dsd_row, &
                            dst_row, &
                            htot_row, &
                            dhd_row, &
                            dht_row, &
                            dpe_row, &
                            dpdr_e_row, &
                            gam1_row, &
                            cs_row

        !..declare local variables

        logical :: single_iter, double_iter, converged
        integer :: var, dvar, var1, var2, iter
        double precision :: v_want
        double precision :: v1_want, v2_want
        double precision :: xnew, xtol, dvdx, smallx, error, v
        double precision :: v1, v2, dv1dt, dv1dr, dv2dt,dv2dr, delr, error1, error2, told, rold, tnew, rnew, v1i, v2i

        double precision :: x,y,zz,zzi,deni,tempi,xni,dxnidd,dxnida, &
                            dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt, &
                            dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt, &
                            deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt, &
                            dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion, &
                            sion,xnem,pele,eele,sele,pres,ener,entr,dpresdd, &
                            dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp, &
                            gam1,gam2,gam3,chit,chid,nabad,sound,etaele, &
                            detadt,detadd,xnefer,dxnedt,dxnedd,s, &
                            temp,den,abar,zbar,ytot1,ye

        !..for the abar derivatives
        double precision :: dpradda,deradda,dsradda, &
                            dpionda,deionda,dsionda, &
                            dpepda,deepda,dsepda,    &
                            dpresda,denerda,dentrda, &
                            detada,dxneda


        !..for the zbar derivatives
        double precision :: dpraddz,deraddz,dsraddz, &
                            dpiondz,deiondz,dsiondz, &
                            dpepdz,deepdz,dsepdz,    &
                            dpresdz,denerdz,dentrdz ,&
                            detadz,dxnedz

        !..for the interpolations
        integer          :: iat,jat
        double precision :: free,df_d,df_t,df_tt,df_dt
        double precision :: xt,xd,mxt,mxd, &
                            si0t,si1t,si2t,si0mt,si1mt,si2mt, &
                            si0d,si1d,si2d,si0md,si1md,si2md, &
                            dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, &
                            dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, &
                            ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt, &
                            z,din,fi(36)

        !..for the coulomb corrections
        double precision :: ecoul,decouldd,decouldt,decoulda,decouldz, &
                            pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz, &
                            scoul,dscouldd,dscouldt,dscoulda,dscouldz

        double precision :: smallt, smalld

        !$gpu

        call eos_get_small_temp(smallt)
        call eos_get_small_dens(smalld)

        temp_row = state % T
        den_row  = state % rho
        abar_row = state % abar
        zbar_row = state % zbar
        ye_row   = state % y_e

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

        ptot_row = 0.0d0
        dpt_row = 0.0d0
        dpd_row = 0.0d0
        dpa_row = 0.0d0
        dpz_row = 0.0d0
        dpe_row = 0.0d0
        dpdr_e_row = 0.0d0

        etot_row = 0.0d0
        det_row = 0.0d0
        ded_row = 0.0d0
        dea_row = 0.0d0
        dez_row = 0.0d0

        stot_row = 0.0d0
        dst_row = 0.0d0
        dsd_row = 0.0d0

        htot_row = 0.0d0
        dhd_row = 0.0d0
        dht_row = 0.0d0

        pele_row = 0.0d0
        ppos_row = 0.0d0

        xne_row = 0.0d0
        xnp_row = 0.0d0

        etaele_row = 0.0d0

        cv_row = 0.0d0
        cp_row = 0.0d0
        cs_row = 0.0d0
        gam1_row = 0.0d0

        converged = .false.

        if (input .eq. eos_input_rt) converged = .true.

        do iter = 1, max_newton

           temp  = temp_row
           den   =  den_row
           abar  = abar_row
           zbar  = zbar_row

           ytot1 = 1.0d0 / abar
           ye    = ye_row
           din   = ye * den

           !..initialize
           deni    = 1.0d0/den
           tempi   = 1.0d0/temp
           kt      = kerg * temp
           ktinv   = 1.0d0/kt

           call apply_radiation(deni, temp, tempi, &
                                prad, dpraddd, dpraddt, dpradda, dpraddz, &
                                erad, deraddd, deraddt, deradda, deraddz, &
                                srad, dsraddd, dsraddt, dsradda, dsraddz)

           call apply_ions(den, deni, temp, tempi, kt, abar, ytot1, &
                           xni, dxnidd, dxnida, &
                           pion, dpiondd, dpiondt, dpionda, dpiondz, &
                           eion, deiondd, deiondt, deionda, deiondz, &
                           sion, dsiondd, dsiondt, dsionda, dsiondz)

           call apply_electrons(den, temp, ye, ytot1, xni, zbar, &
                                pele, dpepdt, dpepdd, dpepda, dpepdz, &
                                sele, dsepdt, dsepdd, dsepda, dsepdz, &
                                eele, deepdt, deepdd, deepda, deepdz, &
                                etaele, detadt, detadd, xnefer)

           if (do_coulomb) then

              call apply_coulomb_corrections(den, temp, kt, ktinv, abar, zbar, ytot1, xni, &
                                             dpiondd, dpiondt, dxnidd, dxnida, &
                                             prad, pion, pele, erad, eion, eele, &
                                             ecoul, decouldd, decouldt, decoulda, decouldz, &
                                             pcoul, dpcouldd, dpcouldt, dpcoulda, dpcouldz, &
                                             scoul, dscouldd, dscouldt, dscoulda, dscouldz)

           end if

           !..sum all the components
           pres    = prad + pion + pele + pcoul
           ener    = erad + eion + eele + ecoul
           entr    = srad + sion + sele + scoul

           dpresdd = dpraddd + dpiondd + dpepdd + dpcouldd
           dpresdt = dpraddt + dpiondt + dpepdt + dpcouldt
#ifdef EXTRA_THERMO
           dpresda = dpradda + dpionda + dpepda + dpcoulda
           dpresdz = dpraddz + dpiondz + dpepdz + dpcouldz
#endif
           denerdd = deraddd + deiondd + deepdd + decouldd
           denerdt = deraddt + deiondt + deepdt + decouldt
#ifdef EXTRA_THERMO
           denerda = deradda + deionda + deepda + decoulda
           denerdz = deraddz + deiondz + deepdz + decouldz
#endif

           dentrdd = dsraddd + dsiondd + dsepdd + dscouldd
           dentrdt = dsraddt + dsiondt + dsepdt + dscouldt
#ifdef EXTRA_THERMO
           dentrda = dsradda + dsionda + dsepda + dscoulda
           dentrdz = dsraddz + dsiondz + dsepdz + dscouldz
#endif

           !..the temperature and density exponents (c&g 9.81 9.82)
           !..the specific heat at constant volume (c&g 9.92)
           !..the third adiabatic exponent (c&g 9.93)
           !..the first adiabatic exponent (c&g 9.97)
           !..the second adiabatic exponent (c&g 9.105)
           !..the specific heat at constant pressure (c&g 9.98)
           !..and relativistic formula for the sound speed (c&g 14.29)
           zz    = pres*deni
           zzi   = den/pres
           chit  = temp/pres * dpresdt
           chid  = dpresdd*zzi
           cv    = denerdt
           x     = zz * chit/(temp * cv)
           gam3  = x + 1.0d0
           gam1  = chit*x + chid
           nabad = x/gam1
           gam2  = 1.0d0/(1.0d0 - nabad)
           cp    = cv * gam1/chid
           z     = 1.0d0 + (ener + light2)*zzi
           sound = clight * sqrt(gam1/z)

           !..maxwell relations; each is zero if the consistency is perfect
           x   = den * den
           dse = temp*dentrdt/denerdt - 1.0d0
           dpe = (denerdd*x + temp*dpresdt)/pres - 1.0d0
           dsp = -dentrdd*x/dpresdt - 1.0d0

           ptot_row = pres
           dpt_row = dpresdt
           dpd_row = dpresdd
#ifdef EXTRA_THERMO
           dpa_row = dpresda
           dpz_row = dpresdz
#endif
           dpe_row = dpresdt / denerdt
           dpdr_e_row = dpresdd - dpresdt * denerdd / denerdt

           etot_row = ener
           det_row = denerdt
           ded_row = denerdd
#ifdef EXTRA_THERMO
           dea_row = denerda
           dez_row = denerdz
#endif

           stot_row = entr
           dst_row = dentrdt
           dsd_row = dentrdd

           htot_row = ener + pres / den
           dhd_row = denerdd + dpresdd / den - pres / den**2
           dht_row = denerdt + dpresdt / den

           pele_row = pele
           ppos_row = 0.0d0

           xne_row = xnefer
           xnp_row = 0.0d0

           etaele_row = etaele

           cv_row = cv
           cp_row = cp
           cs_row = sound
           gam1_row = gam1

           if (converged) then

              exit

           elseif (single_iter) then

              if (dvar .eq. itemp) then

                 x = temp_row
                 smallx = smallt
                 xtol = ttol

                 if (var .eq. ipres) then
                    v    = ptot_row
                    dvdx = dpt_row
                 elseif (var .eq. iener) then
                    v    = etot_row
                    dvdx = det_row
                 elseif (var .eq. ientr) then
                    v    = stot_row
                    dvdx = dst_row
                 elseif (var .eq. ienth) then
                    v    = htot_row
                    dvdx = dht_row
                 else
                    exit
                 endif

              else ! dvar == density

                 x = den_row
                 smallx = smalld
                 xtol = dtol

                 if (var .eq. ipres) then
                    v    = ptot_row
                    dvdx = dpd_row
                 elseif (var .eq. iener) then
                    v    = etot_row
                    dvdx = ded_row
                 elseif (var .eq. ientr) then
                    v    = stot_row
                    dvdx = dsd_row
                 elseif (var .eq. ienth) then
                    v    = htot_row
                    dvdx = dhd_row
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
                 v1    = ptot_row
                 dv1dt = dpt_row
                 dv1dr = dpd_row
              elseif (var1 .eq. iener) then
                 v1    = etot_row
                 dv1dt = det_row
                 dv1dr = ded_row
              elseif (var1 .eq. ientr) then
                 v1    = stot_row
                 dv1dt = dst_row
                 dv1dr = dsd_row
              elseif (var1 .eq. ienth) then
                 v1    = htot_row
                 dv1dt = dht_row
                 dv1dr = dhd_row
              else
                 exit
              endif

              if (var2 .eq. ipres) then
                 v2    = ptot_row
                 dv2dt = dpt_row
                 dv2dr = dpd_row
              elseif (var2 .eq. iener) then
                 v2    = etot_row
                 dv2dt = det_row
                 dv2dr = ded_row
              elseif (var2 .eq. ientr) then
                 v2    = stot_row
                 dv2dt = dst_row
                 dv2dr = dsd_row
              elseif (var2 .eq. ienth) then
                 v2    = htot_row
                 dv2dt = dht_row
                 dv2dr = dhd_row
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

        state % p    = ptot_row
        state % dpdT = dpt_row
        state % dpdr = dpd_row

#ifdef EXTRA_THERMO
        state % dpdA = dpa_row
        state % dpdZ = dpz_row
#endif

        state % dpde = dpe_row
        state % dpdr_e = dpdr_e_row

        state % e    = etot_row
        state % dedT = det_row
        state % dedr = ded_row

#ifdef EXTRA_THERMO
        state % dedA = dea_row
        state % dedZ = dez_row
#endif

        state % s    = stot_row
        state % dsdT = dst_row
        state % dsdr = dsd_row

        state % h    = htot_row
        state % dhdR = dhd_row
        state % dhdT = dht_row

        state % pele = pele_row
        state % ppos = ppos_row

        state % xne = xne_row
        state % xnp = xnp_row

        state % eta = etaele_row

        state % cv   = cv_row
        state % cp   = cp_row
        state % gam1 = gam1_row
        ! state % cs   = cs_row

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




    subroutine apply_electrons(den, temp, ye, ytot1, xni, zbar, &
                               pele, dpepdt, dpepdd, dpepda, dpepdz, &
                               sele, dsepdt, dsepdd, dsepda, dsepdz, &
                               eele, deepdt, deepdd, deepda, deepdz, &
                               etaele, detadt, detadd, xnefer)

      implicit none

      double precision, intent(in   ) :: den, temp, ye, ytot1, xni, zbar
      double precision, intent(inout) :: pele, dpepdt, dpepdd, dpepda, dpepdz
      double precision, intent(inout) :: sele, dsepdt, dsepdd, dsepda, dsepdz
      double precision, intent(inout) :: eele, deepdt, deepdd, deepda, deepdz
      double precision, intent(inout) :: etaele, detadt, detadd, xnefer

      double precision :: dxnedt, dxnedd, xnem, din, x, s

      !..for the interpolations
      integer          :: iat,jat
      double precision :: free,df_d,df_t,df_tt,df_dt
      double precision :: xt,xd,mxt,mxd, &
                          si0t,si1t,si2t,si0mt,si1mt,si2mt, &
                          si0d,si1d,si2d,si0md,si1md,si2md, &
                          dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, &
                          dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, &
                          ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt, &
                          fi(36)

      !..assume complete ionization
      xnem = xni * zbar

      !..enter the table with ye*den
      din = ye*den

      !..hash locate this temperature and density
      jat = int((log10(temp) - tlo)*tstpi) + 1
      jat = max(1,min(jat,jtmax-1))
      iat = int((log10(din) - dlo)*dstpi) + 1
      iat = max(1,min(iat,itmax-1))

      !..access the table locations only once
      fi(1)  = f(iat,jat)
      fi(2)  = f(iat+1,jat)
      fi(3)  = f(iat,jat+1)
      fi(4)  = f(iat+1,jat+1)
      fi(5)  = ft(iat,jat)
      fi(6)  = ft(iat+1,jat)
      fi(7)  = ft(iat,jat+1)
      fi(8)  = ft(iat+1,jat+1)
      fi(9)  = ftt(iat,jat)
      fi(10) = ftt(iat+1,jat)
      fi(11) = ftt(iat,jat+1)
      fi(12) = ftt(iat+1,jat+1)
      fi(13) = fd(iat,jat)
      fi(14) = fd(iat+1,jat)
      fi(15) = fd(iat,jat+1)
      fi(16) = fd(iat+1,jat+1)
      fi(17) = fdd(iat,jat)
      fi(18) = fdd(iat+1,jat)
      fi(19) = fdd(iat,jat+1)
      fi(20) = fdd(iat+1,jat+1)
      fi(21) = fdt(iat,jat)
      fi(22) = fdt(iat+1,jat)
      fi(23) = fdt(iat,jat+1)
      fi(24) = fdt(iat+1,jat+1)
      fi(25) = fddt(iat,jat)
      fi(26) = fddt(iat+1,jat)
      fi(27) = fddt(iat,jat+1)
      fi(28) = fddt(iat+1,jat+1)
      fi(29) = fdtt(iat,jat)
      fi(30) = fdtt(iat+1,jat)
      fi(31) = fdtt(iat,jat+1)
      fi(32) = fdtt(iat+1,jat+1)
      fi(33) = fddtt(iat,jat)
      fi(34) = fddtt(iat+1,jat)
      fi(35) = fddtt(iat,jat+1)
      fi(36) = fddtt(iat+1,jat+1)

      !..various differences
      xt  = max( (temp - t(jat))*dti_sav(jat), 0.0d0)
      xd  = max( (din - d(iat))*ddi_sav(iat), 0.0d0)
      mxt = 1.0d0 - xt
      mxd = 1.0d0 - xd

      !..the six density and six temperature basis functions
      si0t =   psi0(xt)
      si1t =   psi1(xt)*dt_sav(jat)
      si2t =   psi2(xt)*dt2_sav(jat)

      si0mt =  psi0(mxt)
      si1mt = -psi1(mxt)*dt_sav(jat)
      si2mt =  psi2(mxt)*dt2_sav(jat)

      si0d =   psi0(xd)
      si1d =   psi1(xd)*dd_sav(iat)
      si2d =   psi2(xd)*dd2_sav(iat)

      si0md =  psi0(mxd)
      si1md = -psi1(mxd)*dd_sav(iat)
      si2md =  psi2(mxd)*dd2_sav(iat)

      !..derivatives of the weight functions
      dsi0t =   dpsi0(xt)*dti_sav(jat)
      dsi1t =   dpsi1(xt)
      dsi2t =   dpsi2(xt)*dt_sav(jat)

      dsi0mt = -dpsi0(mxt)*dti_sav(jat)
      dsi1mt =  dpsi1(mxt)
      dsi2mt = -dpsi2(mxt)*dt_sav(jat)

      dsi0d =   dpsi0(xd)*ddi_sav(iat)
      dsi1d =   dpsi1(xd)
      dsi2d =   dpsi2(xd)*dd_sav(iat)

      dsi0md = -dpsi0(mxd)*ddi_sav(iat)
      dsi1md =  dpsi1(mxd)
      dsi2md = -dpsi2(mxd)*dd_sav(iat)

      !..second derivatives of the weight functions
      ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
      ddsi1t =   ddpsi1(xt)*dti_sav(jat)
      ddsi2t =   ddpsi2(xt)

      ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
      ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
      ddsi2mt =  ddpsi2(mxt)

      !..the free energy
      free  = h5(fi, &
                 si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
                 si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

      !..derivative with respect to density
      df_d  = h5(fi, &
                 si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
                 dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)

      !..derivative with respect to temperature
      df_t = h5(fi, &
                dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
                si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

      !..derivative with respect to temperature**2
      df_tt = h5(fi, &
                 ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
                 si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

      !..derivative with respect to temperature and density
      df_dt = h5(fi, &
                 dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
                 dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)

      !..now get the pressure derivative with density, chemical potential, and
      !..electron positron number densities
      !..get the interpolation weight functions
      si0t   =  xpsi0(xt)
      si1t   =  xpsi1(xt)*dt_sav(jat)

      si0mt  =  xpsi0(mxt)
      si1mt  =  -xpsi1(mxt)*dt_sav(jat)

      si0d   =  xpsi0(xd)
      si1d   =  xpsi1(xd)*dd_sav(iat)

      si0md  =  xpsi0(mxd)
      si1md  =  -xpsi1(mxd)*dd_sav(iat)

      !..derivatives of weight functions
      dsi0t  = xdpsi0(xt)*dti_sav(jat)
      dsi1t  = xdpsi1(xt)

      dsi0mt = -xdpsi0(mxt)*dti_sav(jat)
      dsi1mt = xdpsi1(mxt)

      dsi0d  = xdpsi0(xd)*ddi_sav(iat)
      dsi1d  = xdpsi1(xd)

      dsi0md = -xdpsi0(mxd)*ddi_sav(iat)
      dsi1md = xdpsi1(mxd)

      !..look in the pressure derivative only once
      fi(1)  = dpdf(iat,jat)
      fi(2)  = dpdf(iat+1,jat)
      fi(3)  = dpdf(iat,jat+1)
      fi(4)  = dpdf(iat+1,jat+1)
      fi(5)  = dpdft(iat,jat)
      fi(6)  = dpdft(iat+1,jat)
      fi(7)  = dpdft(iat,jat+1)
      fi(8)  = dpdft(iat+1,jat+1)
      fi(9)  = dpdfd(iat,jat)
      fi(10) = dpdfd(iat+1,jat)
      fi(11) = dpdfd(iat,jat+1)
      fi(12) = dpdfd(iat+1,jat+1)
      fi(13) = dpdfdt(iat,jat)
      fi(14) = dpdfdt(iat+1,jat)
      fi(15) = dpdfdt(iat,jat+1)
      fi(16) = dpdfdt(iat+1,jat+1)

      !..pressure derivative with density
      dpepdd  = h3(fi, &
                   si0t,   si1t,   si0mt,   si1mt, &
                   si0d,   si1d,   si0md,   si1md)
      dpepdd  = max(ye * dpepdd,0.0d0)

      !..look in the electron chemical potential table only once
      fi(1)  = ef(iat,jat)
      fi(2)  = ef(iat+1,jat)
      fi(3)  = ef(iat,jat+1)
      fi(4)  = ef(iat+1,jat+1)
      fi(5)  = eft(iat,jat)
      fi(6)  = eft(iat+1,jat)
      fi(7)  = eft(iat,jat+1)
      fi(8)  = eft(iat+1,jat+1)
      fi(9)  = efd(iat,jat)
      fi(10) = efd(iat+1,jat)
      fi(11) = efd(iat,jat+1)
      fi(12) = efd(iat+1,jat+1)
      fi(13) = efdt(iat,jat)
      fi(14) = efdt(iat+1,jat)
      fi(15) = efdt(iat,jat+1)
      fi(16) = efdt(iat+1,jat+1)

      !..electron chemical potential etaele
      etaele  = h3(fi, &
                   si0t,   si1t,   si0mt,   si1mt, &
                   si0d,   si1d,   si0md,   si1md)

      !..derivative with respect to density
      x       = h3(fi, &
                   si0t,   si1t,   si0mt,   si1mt, &
                   dsi0d,  dsi1d,  dsi0md,  dsi1md)
      detadd  = ye * x

      !..derivative with respect to temperature
      detadt  = h3(fi, &
                   dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
                   si0d,   si1d,   si0md,   si1md)

#ifdef EXTRA_THERMO
      !..derivative with respect to abar and zbar
      detada = -x * din * ytot1
      detadz =  x * den * ytot1
#endif

      !..look in the number density table only once
      fi(1)  = xf(iat,jat)
      fi(2)  = xf(iat+1,jat)
      fi(3)  = xf(iat,jat+1)
      fi(4)  = xf(iat+1,jat+1)
      fi(5)  = xft(iat,jat)
      fi(6)  = xft(iat+1,jat)
      fi(7)  = xft(iat,jat+1)
      fi(8)  = xft(iat+1,jat+1)
      fi(9)  = xfd(iat,jat)
      fi(10) = xfd(iat+1,jat)
      fi(11) = xfd(iat,jat+1)
      fi(12) = xfd(iat+1,jat+1)
      fi(13) = xfdt(iat,jat)
      fi(14) = xfdt(iat+1,jat)
      fi(15) = xfdt(iat,jat+1)
      fi(16) = xfdt(iat+1,jat+1)

      !..electron + positron number densities
      xnefer   = h3(fi, &
                    si0t,   si1t,   si0mt,   si1mt, &
                    si0d,   si1d,   si0md,   si1md)

      !..derivative with respect to density
      x        = h3(fi, &
                    si0t,   si1t,   si0mt,   si1mt, &
                    dsi0d,  dsi1d,  dsi0md,  dsi1md)
      x = max(x,0.0d0)
      dxnedd   = ye * x

      !..derivative with respect to temperature
      dxnedt   = h3(fi, &
                    dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
                    si0d,   si1d,   si0md,   si1md)

#ifdef EXTRA_THERMO
      !..derivative with respect to abar and zbar
      dxneda = -x * din * ytot1
      dxnedz =  x  * den * ytot1
#endif

      !..the desired electron-positron thermodynamic quantities

      !..dpepdd at high temperatures and low densities is below the
      !..floating point limit of the subtraction of two large terms.
      !..since dpresdd doesn't enter the maxwell relations at all, use the
      !..bicubic interpolation done above instead of this one
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
      
    end subroutine apply_electrons



    subroutine actual_eos_init

        use amrex_error_module
        use extern_probin_module, only: eos_input_is_constant, use_eos_coulomb, eos_ttol, eos_dtol
#ifndef COMPILE_WITH_F2PY
        use amrex_paralleldescriptor_module, only: parallel_bcast => amrex_pd_bcast, amrex_pd_ioprocessor
#endif

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
        allocate(f(imax,jmax))
        allocate(fd(imax,jmax))
        allocate(ft(imax,jmax))
        allocate(fdd(imax,jmax))
        allocate(ftt(imax,jmax))
        allocate(fdt(imax,jmax))
        allocate(fddt(imax,jmax))
        allocate(fdtt(imax,jmax))
        allocate(fddtt(imax,jmax))
        allocate(dpdf(imax,jmax))
        allocate(dpdfd(imax,jmax))
        allocate(dpdft(imax,jmax))
        allocate(dpdfdt(imax,jmax))
        allocate(ef(imax,jmax))
        allocate(efd(imax,jmax))
        allocate(eft(imax,jmax))
        allocate(efdt(imax,jmax))
        allocate(xf(imax,jmax))
        allocate(xfd(imax,jmax))
        allocate(xft(imax,jmax))
        allocate(xfdt(imax,jmax))
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

#ifndef COMPILE_WITH_F2PY
        if (amrex_pd_ioprocessor()) then
#endif
           print *, ''
           if (do_coulomb) then
              print *, "Initializing Helmholtz EOS and using Coulomb corrections."
           else
              print *, "Initializing Helmholtz EOS without using Coulomb corrections."
           endif
           print *, ''
#ifndef COMPILE_WITH_F2PY
        endif
#endif

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

#ifndef COMPILE_WITH_F2PY
        if (amrex_pd_ioprocessor()) then
#endif

           !..   open the table
           open(unit=2,file='helm_table.dat',status='old',iostat=status,action='read')
           if (status > 0) then

              call amrex_error('actual_eos_init: Failed to open helm_table.dat')

           endif

           !...  read in the free energy table
           do j=1,jmax
              do i=1,imax
                 read(2,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j), &
                      fddt(i,j),fdtt(i,j),fddtt(i,j)
              end do
           end do

           !..   read the pressure derivative with density table
           do j = 1, jmax
              do i = 1, imax
                 read(2,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
              end do
           end do

           !..   read the electron chemical potential table
           do j = 1, jmax
              do i = 1, imax
                 read(2,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
              end do
           end do

           !..   read the number density table
           do j = 1, jmax
              do i = 1, imax
                 read(2,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
              end do
           end do
#ifndef COMPILE_WITH_F2PY
        end if
#endif

#ifndef COMPILE_WITH_F2PY
        call parallel_bcast(f)
        call parallel_bcast(fd)
        call parallel_bcast(ft)
        call parallel_bcast(fdd)
        call parallel_bcast(ftt)
        call parallel_bcast(fdt)
        call parallel_bcast(fddt)
        call parallel_bcast(fdtt)
        call parallel_bcast(fddtt)
        call parallel_bcast(dpdf)
        call parallel_bcast(dpdfd)
        call parallel_bcast(dpdft)
        call parallel_bcast(dpdfdt)
        call parallel_bcast(ef)
        call parallel_bcast(efd)
        call parallel_bcast(eft)
        call parallel_bcast(efdt)
        call parallel_bcast(xf)
        call parallel_bcast(xfd)
        call parallel_bcast(xft)
        call parallel_bcast(xfdt)
#endif

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

#ifndef COMPILE_WITH_F2PY
        if (amrex_pd_ioprocessor()) then
#endif
           close(unit=2)
#ifndef COMPILE_WITH_F2PY
        endif
#endif

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
        !$acc device(f, fd, ft, fdd, ftt, fdt, fddt, fdtt, fddtt) &
        !$acc device(dpdf, dpdfd, dpdft, dpdfdt) &
        !$acc device(ef, efd, eft, efdt, xf, xfd, xft, xfdt)  &
        !$acc device(dt_sav, dt2_sav, dti_sav, dt2i_sav) &
        !$acc device(dd_sav, dd2_sav, ddi_sav, dd2i_sav) &
        !$acc device(do_coulomb, input_is_constant)

    end subroutine actual_eos_init



    ! quintic hermite polynomial functions
    ! psi0 and its derivatives
    pure function psi0(z) result(psi0r)
      !$acc routine seq
      double precision, intent(in) :: z
      double precision :: psi0r
      !$gpu
      psi0r = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
    end function psi0

    pure function dpsi0(z) result(dpsi0r)
      !$acc routine seq
      double precision, intent(in) :: z
      double precision :: dpsi0r
      !$gpu
      dpsi0r = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
    end function dpsi0

    pure function ddpsi0(z) result(ddpsi0r)
      !$acc routine seq
      double precision, intent(in) :: z
      double precision :: ddpsi0r
      !$gpu
      ddpsi0r = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)
    end function ddpsi0

    ! psi1 and its derivatives
    pure function psi1(z) result(psi1r)
      !$acc routine seq
      double precision, intent(in) :: z
      double precision :: psi1r
      !$gpu
      psi1r = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
    end function psi1

    pure function dpsi1(z) result(dpsi1r)
      !$acc routine seq
      double precision, intent(in) :: z
      double precision :: dpsi1r
      !$gpu
      dpsi1r = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
    end function dpsi1

    pure function ddpsi1(z) result(ddpsi1r)
      !$acc routine seq
      double precision, intent(in) :: z
      double precision :: ddpsi1r
      !$gpu
      ddpsi1r = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)
    end function ddpsi1

    ! psi2  and its derivatives
    pure function psi2(z) result(psi2r)
      !$acc routine seq
      double precision, intent(in) :: z
      double precision :: psi2r
      !$gpu
      psi2r = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
    end function psi2

    pure function dpsi2(z) result(dpsi2r)
      !$acc routine seq
      double precision, intent(in) :: z
      double precision :: dpsi2r
      !$gpu
      dpsi2r = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
    end function dpsi2

    pure function ddpsi2(z) result(ddpsi2r)
      !$acc routine seq
      double precision, intent(in) :: z
      double precision :: ddpsi2r
      !$gpu
      ddpsi2r = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)
    end function ddpsi2


    ! biquintic hermite polynomial function
    pure function h5(fi,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md) result(h5r)
      !$acc routine seq
      double precision, intent(in) :: fi(36)
      double precision, intent(in) :: w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md
      double precision :: h5r

      !$gpu

      h5r =  fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t &
           + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt &
           + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t &
           + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt &
           + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t &
           + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt &
           + fi(13) *w1d*w0t   + fi(14) *w1md*w0t &
           + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt &
           + fi(17) *w2d*w0t   + fi(18) *w2md*w0t &
           + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt &
           + fi(21) *w1d*w1t   + fi(22) *w1md*w1t &
           + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt &
           + fi(25) *w2d*w1t   + fi(26) *w2md*w1t &
           + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt &
           + fi(29) *w1d*w2t   + fi(30) *w1md*w2t &
           + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt &
           + fi(33) *w2d*w2t   + fi(34) *w2md*w2t &
           + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt
    end function h5


    ! cubic hermite polynomial functions
    ! psi0 & derivatives
    pure function xpsi0(z) result(xpsi0r)
      !$acc routine seq
      double precision, intent(in) :: z
      double precision :: xpsi0r
      !$gpu
      xpsi0r = z * z * (2.0d0*z - 3.0d0) + 1.0
    end function xpsi0

    pure function xdpsi0(z) result(xdpsi0r)
      !$acc routine seq
      double precision, intent(in) :: z
      double precision :: xdpsi0r
      !$gpu
      xdpsi0r = z * (6.0d0*z - 6.0d0)
    end function xdpsi0


    ! psi1 & derivatives
    pure function xpsi1(z) result(xpsi1r)
      !$acc routine seq
      double precision, intent(in) :: z
      double precision :: xpsi1r
      !$gpu
      xpsi1r = z * ( z * (z - 2.0d0) + 1.0d0)
    end function xpsi1

    pure function xdpsi1(z) result(xdpsi1r)
      !$acc routine seq
      double precision, intent(in) :: z
      double precision :: xdpsi1r
      !$gpu
      xdpsi1r = z * (3.0d0*z - 4.0d0) + 1.0d0
    end function xdpsi1

    ! bicubic hermite polynomial function
    pure function h3(fi,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) result(h3r)
      !$acc routine seq
      double precision, intent(in) :: fi(36)
      double precision, intent(in) :: w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md
      double precision :: h3r
      !$gpu
      h3r =   fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t &
           + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt &
           + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t &
           + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt &
           + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t &
           + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt &
           + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t &
           + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt
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
      deallocate(fd)
      deallocate(ft)
      deallocate(fdd)
      deallocate(ftt)
      deallocate(fdt)
      deallocate(fddt)
      deallocate(fdtt)
      deallocate(fddtt)
      deallocate(dpdf)
      deallocate(dpdfd)
      deallocate(dpdft)
      deallocate(dpdfdt)
      deallocate(ef)
      deallocate(efd)
      deallocate(eft)
      deallocate(efdt)
      deallocate(xf)
      deallocate(xfd)
      deallocate(xft)
      deallocate(xfdt)
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
