module actual_eos_module

  use eos_type_module

  use amrex_fort_module, only : rt => amrex_real
  character (len=64), public :: eos_name = "helmholtz"

  ! Runtime parameters
  logical, allocatable :: do_coulomb
  logical, allocatable :: input_is_constant

  !..for the tables, in general
  integer, parameter, private :: imax = 541, jmax = 201
  integer, allocatable :: itmax, jtmax
  real(rt)        , allocatable :: d(:), t(:)

  real(rt)        , allocatable :: tlo, thi, tstp, tstpi
  real(rt)        , allocatable :: dlo, dhi, dstp, dstpi

  real(rt)        , allocatable :: ttol, dtol

  !..for the helmholtz free energy tables
  real(rt)        , allocatable :: f(:,:,:)

  !..for the pressure derivative with density tables
  real(rt)        , allocatable :: dpdf(:,:,:)

  !..for chemical potential tables
  real(rt)        , allocatable :: ef(:,:,:)

  !..for the number density tables
  real(rt)        , allocatable :: xf(:,:,:)

  !..for storing the differences
  real(rt)        , allocatable :: dt_sav(:), dt2_sav(:),          &
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
  attributes(managed) :: f, dpdf, ef, xf
  attributes(managed) :: dt_sav, dt2_sav, dti_sav, dt2i_sav
  attributes(managed) :: dd_sav, dd2_sav, ddi_sav, dd2i_sav
#endif

  ! 2006 CODATA physical constants
  real(rt)        , parameter :: h       = 6.6260689633e-27_rt
  real(rt)        , parameter :: avo_eos = 6.0221417930e23_rt
  real(rt)        , parameter :: kerg    = 1.380650424e-16_rt
  real(rt)        , parameter :: amu     = 1.66053878283e-24_rt

  !$acc declare &
  !$acc create(tlo, thi, dlo, dhi) &
  !$acc create(tstp, tstpi, dstp, dstpi) &
  !$acc create(ttol, dtol) &
  !$acc create(itmax, jtmax, d, t) &
  !$acc create(f, dpdf, ef, xf) &
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

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    !..input arguments
    integer,      intent(in   ) :: input
    type (eos_t), intent(inout) :: state

    integer, parameter :: max_newton = 100

    logical :: single_iter, converged
    integer :: var, dvar, var1, var2, iter
    real(rt)         :: v_want, v1_want, v2_want

    !$gpu

    ! Initial setup for iterations.

    call prepare_for_iterations(input, state, single_iter, v_want, v1_want, v2_want, var, dvar, var1, var2)

    converged = .false.

    ! Only take a single step if we're coming in with both rho and T;
    ! in this call we just want the EOS to fill the other thermodynamic quantities.

    if (input .eq. eos_input_rt) converged = .true.

    ! Iterate until converged.

    do iter = 1, max_newton

       ! Radiation must come first since it initializes the
       ! state instead of adding to it.

       call apply_radiation(state)

       call apply_ions(state)

       call apply_electrons(state)

       if (do_coulomb) then
          call apply_coulomb_corrections(state)
       end if

       ! Calculate enthalpy the usual way, h = e + p / rho.

       state % h = state % e + state % p / state % rho
       state % dhdr = state % dedr + state % dpdr / state % rho - state % p / state % rho**2
       state % dhdT = state % dedT + state % dpdT / state % rho

       if (converged) then
          exit
       elseif (single_iter) then
          call single_iter_update(state, var, dvar, v_want, converged)
       else
          call double_iter_update(state, var1, var2, v1_want, v2_want, converged) 
       endif

    enddo

    call finalize_state(input, state, v_want, v1_want, v2_want)

  end subroutine actual_eos



  subroutine apply_electrons(state)

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type(eos_t), intent(inout) :: state

    real(rt)         :: pele, dpepdt, dpepdd, dpepda, dpepdz
    real(rt)         :: sele, dsepdt, dsepdd, dsepda, dsepdz
    real(rt)         :: eele, deepdt, deepdd, deepda, deepdz
    real(rt)         :: xni, xnem, din, x, s, etaele, xnefer, ytot1

    !..for the interpolations
    integer          :: iat, jat
    real(rt)         :: free, df_d, df_t, df_tt, df_dt
    real(rt)         :: xt, xd, mxt, mxd
    real(rt)         :: fi(36)
    real(rt)         :: wdt(16), sit(6), sid(6), dsit(6), dsid(6), fwtr(6), ddsit(6)

    !$gpu

    !..assume complete ionization
    ytot1 = 1.0e0_rt / state % abar
    xni  = avo_eos * ytot1 * state % rho
    xnem = xni * state % zbar

    ! define mu -- the total mean molecular weight including both electrons and ions
    state % mu = 1.0e0_rt / (1.0e0_rt / state % abar + 1.0e0_rt / state % mu_e)

    !..enter the table with ye*den
    din = state % y_e * state % rho

    !..hash locate this temperature and density
    jat = int((log10(state % T) - tlo) * tstpi) + 1
    jat = max(1, min(jat, jtmax-1))
    iat = int((log10(din) - dlo) * dstpi) + 1
    iat = max(1, min(iat, itmax-1))

    !..access the table locations only once
    fi(1:9)   = f(1:9, iat  ,jat  ) ! f, ft, ftt, fd, fdd, fdt, fddt, fdtt, fddtt
    fi(10:18) = f(1:9, iat+1,jat  )
    fi(19:27) = f(1:9, iat  ,jat+1)
    fi(28:36) = f(1:9, iat+1,jat+1)

    !..various differences
    xt  = max( (state % T - t(jat)) * dti_sav(jat), 0.0e0_rt)
    xd  = max( (din - d(iat)) * ddi_sav(iat), 0.0e0_rt)
    mxt = 1.0e0_rt - xt
    mxd = 1.0e0_rt - xd

    !..the six density and six temperature basis functions
    sit(1) = psi0(xt)
    sit(2) = psi1(xt) * dt_sav(jat)
    sit(3) = psi2(xt) * dt2_sav(jat)

    sit(4) =  psi0(mxt)
    sit(5) = -psi1(mxt) * dt_sav(jat)
    sit(6) =  psi2(mxt) * dt2_sav(jat)

    sid(1) =   psi0(xd)
    sid(2) =   psi1(xd) * dd_sav(iat)
    sid(3) =   psi2(xd) * dd2_sav(iat)

    sid(4) =  psi0(mxd)
    sid(5) = -psi1(mxd) * dd_sav(iat)
    sid(6) =  psi2(mxd) * dd2_sav(iat)

    !..derivatives of the weight functions
    dsit(1) =   dpsi0(xt) * dti_sav(jat)
    dsit(2) =   dpsi1(xt)
    dsit(3) =   dpsi2(xt) * dt_sav(jat)

    dsit(4) = -dpsi0(mxt) * dti_sav(jat)
    dsit(5) =  dpsi1(mxt)
    dsit(6) = -dpsi2(mxt) * dt_sav(jat)

    dsid(1) =   dpsi0(xd) * ddi_sav(iat)
    dsid(2) =   dpsi1(xd)
    dsid(3) =   dpsi2(xd) * dd_sav(iat)

    dsid(4) = -dpsi0(mxd) * ddi_sav(iat)
    dsid(5) =  dpsi1(mxd)
    dsid(6) = -dpsi2(mxd) * dd_sav(iat)

    !..second derivatives of the weight functions
    ddsit(1) =   ddpsi0(xt) * dt2i_sav(jat)
    ddsit(2) =   ddpsi1(xt) * dti_sav(jat)
    ddsit(3) =   ddpsi2(xt)

    ddsit(4) =  ddpsi0(mxt) * dt2i_sav(jat)
    ddsit(5) = -ddpsi1(mxt) * dti_sav(jat)
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
    sit(2) = xpsi1(xt) * dt_sav(jat)

    sit(3) = xpsi0(mxt)
    sit(4) = -xpsi1(mxt) * dt_sav(jat)

    sid(1) = xpsi0(xd)
    sid(2) = xpsi1(xd) * dd_sav(iat)

    sid(3) = xpsi0(mxd)
    sid(4) = -xpsi1(mxd) * dd_sav(iat)

    !..derivatives of weight functions
    dsit(1) = xdpsi0(xt) * dti_sav(jat)
    dsit(2) = xdpsi1(xt)

    dsit(3) = -xdpsi0(mxt) * dti_sav(jat)
    dsit(4) = xdpsi1(mxt)

    dsid(1) = xdpsi0(xd) * ddi_sav(iat)
    dsid(2) = xdpsi1(xd)

    dsid(3) = -xdpsi0(mxd) * ddi_sav(iat)
    dsid(4) = xdpsi1(mxd)

    ! Reuse subexpressions that would go into computing the
    ! cubic interpolation.
    wdt(1:4)   = sid(1) * sit(1:4)
    wdt(5:8)   = sid(2) * sit(1:4)
    wdt(9:12)  = sid(3) * sit(1:4)
    wdt(13:16) = sid(4) * sit(1:4)

    ! Read in the tabular data for the pressure derivatives.
    ! We have some freedom in how we store it in the local
    ! array. We choose here to index it such that we can
    ! immediately evaluate the cubic interpolant below as
    ! fi * wdt, which ensures that we have the right combination
    ! of grid points and derivatives at grid points to evaluate
    ! the interpolation correctly. Alternate indexing schemes are
    ! possible if we were to reorder wdt.
    fi([ 1,  2,  5,  6]) = dpdf(1:4, iat,   jat  )
    fi([ 9, 10, 13, 14]) = dpdf(1:4, iat+1, jat  )
    fi([ 3,  4,  7,  8]) = dpdf(1:4, iat,   jat+1)
    fi([11, 12, 15, 16]) = dpdf(1:4, iat+1, jat+1)

    !..pressure derivative with density
    dpepdd = sum(fi(1:16) * wdt)
    dpepdd = max(state % y_e * dpepdd, 0.0e0_rt)

    ! Read in the tabular data for the electron chemical potential.
    fi([ 1,  2,  5,  6]) = ef(1:4,iat  ,jat  )
    fi([ 9, 10, 13, 14]) = ef(1:4,iat+1,jat  )
    fi([ 3,  4,  7,  8]) = ef(1:4,iat  ,jat+1)
    fi([11, 12, 15, 16]) = ef(1:4,iat+1,jat+1)

    !..electron chemical potential etaele
    etaele = sum(fi(1:16) * wdt)

    ! Read in the tabular data for the number density.
    fi([ 1,  2,  5,  6]) = xf(1:4,iat  ,jat  )
    fi([ 9, 10, 13, 14]) = xf(1:4,iat+1,jat  )
    fi([ 3,  4,  7,  8]) = xf(1:4,iat  ,jat+1)
    fi([11, 12, 15, 16]) = xf(1:4,iat+1,jat+1)

    !..electron + positron number densities
    xnefer = sum(fi(1:16) * wdt)

    wdt(1:4)   = dsid(1) * sit(1:4)
    wdt(5:8)   = dsid(2) * sit(1:4)
    wdt(9:12)  = dsid(3) * sit(1:4)
    wdt(13:16) = dsid(4) * sit(1:4)

    !..derivative with respect to density
    x = sum(fi(1:16) * wdt)
    x = max(x, 0.0e0_rt)

    !..the desired electron-positron thermodynamic quantities

    !..dpepdd at high temperatures and low densities is below the
    !..floating point limit of the subtraction of two large terms.
    !..since state % dpdr doesn't enter the maxwell relations at all, use the
    !..bicubic interpolation done above instead of this one
    x       = din * din
    pele    = x * df_d
    dpepdt  = x * df_dt
    s       = dpepdd/state % y_e - 2.0e0_rt * din * df_d
#ifdef EXTRA_THERMO
    dpepda  = -ytot1 * (2.0e0_rt * pele + s * din)
    dpepdz  = state % rho*ytot1*(2.0e0_rt * din * df_d  +  s)
#endif

    x       = state % y_e * state % y_e
    sele    = -df_t * state % y_e
    dsepdt  = -df_tt * state % y_e
    dsepdd  = -df_dt * x
#ifdef EXTRA_THERMO
    dsepda  = ytot1 * (state % y_e * df_dt * din - sele)
    dsepdz  = -ytot1 * (state % y_e * df_dt * state % rho  + df_t)
#endif

    eele    = state % y_e*free + state % T * sele
    deepdt  = state % T * dsepdt
    deepdd  = x * df_d + state % T * dsepdd
#ifdef EXTRA_THERMO
    deepda  = -state % y_e * ytot1 * (free +  df_d * din) + state % T * dsepda
    deepdz  = ytot1* (free + state % y_e * df_d * state % rho) + state % T * dsepdz
#endif

    state % p    = state % p + pele
    state % dpdT = state % dpdT + dpepdt
    state % dpdr = state % dpdr + dpepdd
#ifdef EXTRA_THERMO
    state % dpdA = state % dpdA + dpepda
    state % dpdZ = state % dpdZ + dpepdz
#endif

    state % s    = state % s + sele
    state % dsdT = state % dsdT + dsepdt
    state % dsdr = state % dsdr + dsepdd

    state % e    = state % e + eele
    state % dedT = state % dedT + deepdt
    state % dedr = state % dedr + deepdd
#ifdef EXTRA_THERMO
    state % dedA = state % dedA + deepda
    state % dedZ = state % dedZ + deepdz
#endif

    state % eta = etaele
    state % xne = xnefer
    state % xnp = 0.0e0_rt

    state % pele = pele
    state % ppos = 0.0e0_rt

  end subroutine apply_electrons



  subroutine apply_ions(state)

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type(eos_t), intent(inout) :: state

    real(rt)         :: pion, dpiondd, dpiondt, dpionda, dpiondz
    real(rt)         :: eion, deiondd, deiondt, deionda, deiondz
    real(rt)         :: sion, dsiondd, dsiondt, dsionda, dsiondz

    real(rt)         :: xni, dxnidd, dxnida, kt, ytot1, deni, tempi

    real(rt)         :: s, x, y, z

    real(rt)        , parameter :: pi      = 3.1415926535897932384e0_rt
    real(rt)        , parameter :: sioncon = (2.0e0_rt * pi * amu * kerg)/(h*h)
    real(rt)        , parameter :: kergavo = kerg * avo_eos

    !$gpu

    deni = 1.0e0_rt / state % rho
    tempi = 1.0e0_rt / state % T

    ytot1   = 1.0e0_rt / state % abar
    xni     = avo_eos * ytot1 * state % rho
    dxnidd  = avo_eos * ytot1
    dxnida  = -xni * ytot1

    kt = kerg * state % T

    pion    = xni * kt
    dpiondd = dxnidd * kt
    dpiondt = xni * kerg
#ifdef EXTRA_THERMO
    dpionda = dxnida * kt
    dpiondz = 0.0e0_rt
#endif

    eion    = 1.5e0_rt * pion*deni
    deiondd = (1.5e0_rt * dpiondd - eion)*deni
    deiondt = 1.5e0_rt * dpiondt*deni
#ifdef EXTRA_THERMO
    deionda = 1.5e0_rt * dpionda*deni
    deiondz = 0.0e0_rt
#endif

    x       = state % abar*state % abar*sqrt(state % abar) * deni/avo_eos
    s       = sioncon * state % T
    z       = x * s * sqrt(s)
    y       = log(z)
    sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
    dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi &
              - kergavo * deni * ytot1
    dsiondt = (dpiondt*deni + deiondt)*tempi -  &
              (pion*deni + eion) * tempi*tempi  &
              + 1.5e0_rt * kergavo * tempi*ytot1
    x       = avo_eos*kerg/state % abar
#ifdef EXTRA_THERMO
    dsionda = (dpionda*deni + deionda)*tempi  &
              + kergavo*ytot1*ytot1* (2.5e0_rt - y)
    dsiondz = 0.0e0_rt
#endif

    state % p    = state % p + pion
    state % dpdT = state % dpdT + dpiondt
    state % dpdr = state % dpdr + dpiondd
#ifdef EXTRA_THERMO
    state % dpdA = state % dpdA + dpionda
    state % dpdZ = state % dpdZ + dpiondz
#endif

    state % e    = state % e + eion
    state % dedT = state % dedT + deiondt
    state % dedr = state % dedr + deiondd
#ifdef EXTRA_THERMO
    state % dedA = state % dedA + deionda
    state % dedZ = state % dedZ + deiondz
#endif

    state % s    = state % s + sion
    state % dsdT = state % dsdT + dsiondt
    state % dsdr = state % dsdr + dsiondd

  end subroutine apply_ions



  
  subroutine apply_radiation(state)

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type(eos_t), intent(inout) :: state

    real(rt)         :: prad, dpraddd, dpraddt, dpradda, dpraddz
    real(rt)         :: erad, deraddd, deraddt, deradda, deraddz
    real(rt)         :: srad, dsraddd, dsraddt, dsradda, dsraddz

    real(rt)         :: deni, tempi

    real(rt)        , parameter :: clight  = 2.99792458e10_rt
#ifdef RADIATION
    real(rt)        , parameter :: ssol    = 0.0e0_rt
#else
    real(rt)        , parameter :: ssol    = 5.67051e-5_rt
#endif
    real(rt)        , parameter :: asol    = 4.0e0_rt * ssol / clight
    real(rt)        , parameter :: asoli3  = asol/3.0e0_rt

    !$gpu
    
    deni = 1.0e0_rt / state % rho
    tempi = 1.0e0_rt / state % T

    prad    = asoli3 * state % T * state % T * state % T * state % T
    dpraddd = 0.0e0_rt
    dpraddt = 4.0e0_rt * prad*tempi
#ifdef EXTRA_THERMO
    dpradda = 0.0e0_rt
    dpraddz = 0.0e0_rt
#endif

    erad    = 3.0e0_rt * prad*deni
    deraddd = -erad*deni
    deraddt = 3.0e0_rt * dpraddt*deni
#ifdef EXTRA_THERMO
    deradda = 0.0e0_rt
    deraddz = 0.0e0_rt
#endif

    srad    = (prad*deni + erad)*tempi
    dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
    dsraddt = (dpraddt*deni + deraddt - srad)*tempi
#ifdef EXTRA_THERMO
    dsradda = 0.0e0_rt
    dsraddz = 0.0e0_rt
#endif

    ! Note that unlike the other terms, radiation
    ! sets these terms instead of adding to them,
    ! since it comes first.

    state % p    = prad
    state % dpdr = dpraddd
    state % dpdT = dpraddt
#ifdef EXTRA_THERMO
    state % dpdA = dpradda
    state % dpdZ = dpraddz
#endif

    state % e    = erad
    state % dedr = deraddd
    state % dedT = deraddt
#ifdef EXTRA_THERMO
    state % dedA = deradda
    state % dedZ = deraddz
#endif

    state % s    = srad
    state % dsdr = dsraddd
    state % dsdT = dsraddt

  end subroutine apply_radiation



  subroutine apply_coulomb_corrections(state)

    use amrex_constants_module, only: ZERO

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type(eos_t), intent(inout) :: state

    real(rt)         :: ecoul, decouldd, decouldt, decoulda, decouldz
    real(rt)         :: pcoul, dpcouldd, dpcouldt, dpcoulda, dpcouldz
    real(rt)         :: scoul, dscouldd, dscouldt, dscoulda, dscouldz

    real(rt)         :: dsdd, dsda, lami, inv_lami, lamida, lamidd
    real(rt)         :: plasg, plasgdd, plasgdt, plasgda, plasgdz
    real(rt)         :: pion, dpiondt, dpiondd, dpionda, dpiondz, xni, dxnidd, dxnida

    real(rt)         :: kt, ktinv, ytot1
    real(rt)         :: s, x, y, z
    real(rt)         :: p_temp, e_temp

    ! Constants used for the Coulomb corrections
    real(rt)        , parameter :: a1 = -0.898004e0_rt
    real(rt)        , parameter :: b1 =  0.96786e0_rt
    real(rt)        , parameter :: c1 =  0.220703e0_rt
    real(rt)        , parameter :: d1 = -0.86097e0_rt
    real(rt)        , parameter :: e1 =  2.5269e0_rt
    real(rt)        , parameter :: a2 =  0.29561e0_rt
    real(rt)        , parameter :: b2 =  1.9885e0_rt
    real(rt)        , parameter :: c2 =  0.288675e0_rt
    real(rt)        , parameter :: qe   = 4.8032042712e-10_rt
    real(rt)        , parameter :: esqu = qe * qe
    real(rt)        , parameter :: onethird = 1.0e0_rt/3.0e0_rt
    real(rt)        , parameter :: forth = 4.0e0_rt/3.0e0_rt
    real(rt)        , parameter :: pi    = 3.1415926535897932384e0_rt

    !$gpu

    pcoul    = ZERO
    dpcouldd = ZERO
    dpcouldt = ZERO
    dpcoulda = ZERO
    dpcouldz = ZERO
    ecoul    = ZERO
    decouldd = ZERO
    decouldt = ZERO
    decoulda = ZERO
    decouldz = ZERO
    scoul    = ZERO
    dscouldd = ZERO
    dscouldt = ZERO
    dscoulda = ZERO
    dscouldz = ZERO

    !..uniform background corrections only
    !..from yakovlev & shalybkov 1989
    !..lami is the average ion seperation
    !..plasg is the plasma coupling parameter

    ytot1 = 1.0e0_rt / state % abar
    xni     = avo_eos * ytot1 * state % rho
    dxnidd  = avo_eos * ytot1
    dxnida  = -xni * ytot1

    kt      = kerg * state % T
    ktinv   = 1.0e0_rt/kt

    z        = forth * pi
    s        = z * xni
    dsdd     = z * dxnidd
    dsda     = z * dxnida

    lami     = 1.0e0_rt/s**onethird
    inv_lami = 1.0e0_rt/lami
    z        = -onethird * lami
    lamidd   = z * dsdd/s
    lamida   = z * dsda/s

    plasg    = state % zbar*state % zbar*esqu*ktinv*inv_lami
    z        = -plasg * inv_lami
    plasgdd  = z * lamidd
    plasgda  = z * lamida
    plasgdt  = -plasg*ktinv * kerg
    plasgdz  = 2.0e0_rt * plasg/state % zbar

    !...yakovlev & shalybkov 1989 equations 82, 85, 86, 87
    if (plasg .ge. 1.0e0_rt) then
       x        = plasg**(0.25e0_rt)
       y        = avo_eos * ytot1 * kerg
       ecoul    = y * state % T * (a1*plasg + b1*x + c1/x + d1)
       pcoul    = onethird * state % rho * ecoul
       scoul    = -y * (3.0e0_rt*b1*x - 5.0e0_rt*c1/x &
                  + d1 * (log(plasg) - 1.0e0_rt) - e1)

       y        = avo_eos*ytot1*kt*(a1 + 0.25e0_rt/plasg*(b1*x - c1/x))
       decouldd = y * plasgdd
       decouldt = y * plasgdt + ecoul/state % T
       decoulda = y * plasgda - ecoul/state % abar
       decouldz = y * plasgdz

       y        = onethird * state % rho
       dpcouldd = onethird * ecoul + y*decouldd
       dpcouldt = y * decouldt
       dpcoulda = y * decoulda
       dpcouldz = y * decouldz

       y        = -avo_eos*kerg/(state % abar*plasg)* &
                  (0.75e0_rt*b1*x+1.25e0_rt*c1/x+d1)
       dscouldd = y * plasgdd
       dscouldt = y * plasgdt
       dscoulda = y * plasgda - scoul/state % abar
       dscouldz = y * plasgdz

       !...yakovlev & shalybkov 1989 equations 102, 103, 104
    else if (plasg .lt. 1.0e0_rt) then

       pion    = xni * kt
       dpiondd = dxnidd * kt
       dpiondt = xni * kerg
#ifdef EXTRA_THERMO
       dpionda = dxnida * kt
       dpiondz = 0.0e0_rt
#endif

       x        = plasg*sqrt(plasg)
       y        = plasg**b2
       z        = c2 * x - onethird * a2 * y
       pcoul    = -pion * z
       ecoul    = 3.0e0_rt * pcoul/state % rho
       scoul    = -avo_eos/state % abar*kerg*(c2*x -a2*(b2-1.0e0_rt)/b2*y)

       s        = 1.5e0_rt*c2*x/plasg - onethird*a2*b2*y/plasg
       dpcouldd = -dpiondd*z - pion*s*plasgdd
       dpcouldt = -dpiondt*z - pion*s*plasgdt
#ifdef EXTRA_THERMO
       dpcoulda = -dpionda*z - pion*s*plasgda
       dpcouldz = -dpiondz*z - pion*s*plasgdz
#endif

       s        = 3.0e0_rt/state % rho
       decouldd = s * dpcouldd - ecoul/state % rho
       decouldt = s * dpcouldt
       decoulda = s * dpcoulda
       decouldz = s * dpcouldz

       s        = -avo_eos*kerg/(state % abar*plasg)* &
                  (1.5e0_rt*c2*x-a2*(b2-1.0e0_rt)*y)
       dscouldd = s * plasgdd
       dscouldt = s * plasgdt
       dscoulda = s * plasgda - scoul/state % abar
       dscouldz = s * plasgdz
    end if

    ! Disable Coulomb corrections if they cause
    ! the energy or pressure to go negative.

    p_temp = state % p + pcoul
    e_temp = state % e + ecoul

    if (p_temp .le. ZERO .or. e_temp .le. ZERO) then

       pcoul    = 0.0e0_rt
       dpcouldd = 0.0e0_rt
       dpcouldt = 0.0e0_rt
       dpcoulda = 0.0e0_rt
       dpcouldz = 0.0e0_rt
       ecoul    = 0.0e0_rt
       decouldd = 0.0e0_rt
       decouldt = 0.0e0_rt
       decoulda = 0.0e0_rt
       decouldz = 0.0e0_rt
       scoul    = 0.0e0_rt
       dscouldd = 0.0e0_rt
       dscouldt = 0.0e0_rt
       dscoulda = 0.0e0_rt
       dscouldz = 0.0e0_rt

    end if

    state % p    = state % p + pcoul
    state % dpdr = state % dpdr + dpcouldd
    state % dpdT = state % dpdT + dpcouldt
#ifdef EXTRA_THERMO
    state % dpdA = state % dpdA + dpcoulda
    state % dpdZ = state % dpdZ + dpcouldz
#endif

    state % e    = state % e + ecoul
    state % dedr = state % dedr + decouldd
    state % dedT = state % dedT + decouldt
#ifdef EXTRA_THERMO
    state % dedA = state % dedA + decoulda
    state % dedZ = state % dedZ + decouldz
#endif

    state % s    = state % s + scoul
    state % dsdr = state % dsdr + dscouldd
    state % dsdT = state % dsdT + dscouldt

  end subroutine apply_coulomb_corrections



  subroutine prepare_for_iterations(input, state, single_iter, v_want, v1_want, v2_want, var, dvar, var1, var2)

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer,          intent(in   ) :: input
    type(eos_t),      intent(in   ) :: state
    logical,          intent(inout) :: single_iter
    real(rt)        , intent(inout) :: v_want, v1_want, v2_want
    integer,          intent(inout) :: var, dvar, var1, var2

    !$gpu

    single_iter = .true.

    if (input .eq. eos_input_rt) then

       ! Nothing to do here.

    elseif (input .eq. eos_input_rh) then

       v_want = state % h
       var  = ienth
       dvar = itemp

    elseif (input .eq. eos_input_tp) then

       v_want = state % p
       var  = ipres
       dvar = idens

    elseif (input .eq. eos_input_rp) then

       v_want = state % p
       var  = ipres
       dvar = itemp

    elseif (input .eq. eos_input_re) then

       v_want = state % e
       var  = iener
       dvar = itemp

    elseif (input .eq. eos_input_ps) then

       single_iter = .false.
       v1_want = state % p
       v2_want = state % s
       var1 = ipres
       var2 = ientr

    elseif (input .eq. eos_input_ph) then

       single_iter = .false.
       v1_want = state % p
       v2_want = state % h
       var1 = ipres
       var2 = ienth

    elseif (input .eq. eos_input_th) then

       v_want = state % h
       var  = ienth
       dvar = idens

    endif

  end subroutine prepare_for_iterations



  subroutine single_iter_update(state, var, dvar, v_want, converged)

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type(eos_t),      intent(inout) :: state
    integer,          intent(in   ) :: var, dvar
    real(rt)        , intent(in   ) :: v_want
    logical,          intent(inout) :: converged

    real(rt)         :: x, xnew, v, dvdx, xtol, smallx, error

    !$gpu

    if (dvar .eq. itemp) then

       x = state % T
       smallx = mintemp
       xtol = ttol

       if (var .eq. ipres) then
          v    = state % p
          dvdx = state % dpdT
       elseif (var .eq. iener) then
          v    = state % e
          dvdx = state % dedT
       elseif (var .eq. ientr) then
          v    = state % s
          dvdx = state % dsdT
       elseif (var .eq. ienth) then
          v    = state % h
          dvdx = state % dhdT
       endif

    else ! dvar == density

       x = state % rho
       smallx = mindens
       xtol = dtol

       if (var .eq. ipres) then
          v    = state % p
          dvdx = state % dpdr
       elseif (var .eq. iener) then
          v    = state % e
          dvdx = state % dedr
       elseif (var .eq. ientr) then
          v    = state % s
          dvdx = state % dsdr
       elseif (var .eq. ienth) then
          v    = state % h
          dvdx = state % dhdr
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
       state % T = xnew
    else
       state % rho  = xnew
    endif

    ! Compute the error from the last iteration

    error = abs( (xnew - x) / x )

    if (error .lt. xtol) converged = .true.

  end subroutine single_iter_update



  subroutine double_iter_update(state, var1, var2, v1_want, v2_want, converged)

    use amrex_constants_module, only: HALF, TWO

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type(eos_t),      intent(inout) :: state
    integer,          intent(in   ) :: var1, var2
    real(rt)        , intent(in   ) :: v1_want, v2_want
    logical,          intent(inout) :: converged

    real(rt)         :: told, rold, delr, rnew, tnew
    real(rt)         :: v1, dv1dt, dv1dr, v2, dv2dt, dv2dr, v1i, v2i
    real(rt)         :: error1, error2

    !$gpu

    ! Figure out which variables we're using

    told = state % T
    rold = state % rho

    if (var1 .eq. ipres) then
       v1    = state % p
       dv1dt = state % dpdT
       dv1dr = state % dpdr
    elseif (var1 .eq. iener) then
       v1    = state % e
       dv1dt = state % dedT
       dv1dr = state % dedr
    elseif (var1 .eq. ientr) then
       v1    = state % s
       dv1dt = state % dsdT
       dv1dr = state % dsdr
    elseif (var1 .eq. ienth) then
       v1    = state % h
       dv1dt = state % dhdT
       dv1dr = state % dhdr
    endif

    if (var2 .eq. ipres) then
       v2    = state % p
       dv2dt = state % dpdT
       dv2dr = state % dpdr
    elseif (var2 .eq. iener) then
       v2    = state % e
       dv2dt = state % dedT
       dv2dr = state % dedr
    elseif (var2 .eq. ientr) then
       v2    = state % s
       dv2dt = state % dsdT
       dv2dr = state % dsdr
    elseif (var2 .eq. ienth) then
       v2    = state % h
       dv2dt = state % dhdT
       dv2dr = state % dhdr
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
    tnew = max(mintemp, tnew)
    rnew = max(mindens, rnew)

    ! Store the new temperature and density
    state % rho = rnew
    state % T = tnew

    ! Compute the errors
    error1 = abs( (rnew - rold) / rold )
    error2 = abs( (tnew - told) / told )

    if (error1 .LT. dtol .and. error2 .LT. ttol) converged = .true.

  end subroutine double_iter_update



  subroutine finalize_state(input, state, v_want, v1_want, v2_want)

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer,          intent(in   ) :: input
    type(eos_t),      intent(inout) :: state
    real(rt)        , intent(in   ) :: v_want, v1_want, v2_want

    real(rt)         :: chit, chid

    !$gpu

    ! Calculate some remaining derivatives
    state % dpde = state % dpdT / state % dedT
    state % dpdr_e = state % dpdr - state % dpdT * state % dedr / state % dedT

    ! Specific heats and Gamma_1
    chit = state % T / state % p * state % dpdT
    chid = state % dpdr * state % rho / state % p

    state % cv = state % dedT
    state % gam1 = (chit * (state % p / state % rho)) * (chit / (state % T * state % cv)) + chid
    state % cp = state % cv * state % gam1 / chid

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

  end subroutine finalize_state



  subroutine actual_eos_init

    use amrex_error_module
    use extern_probin_module, only: eos_input_is_constant, use_eos_coulomb, eos_ttol, eos_dtol
#ifndef COMPILE_WITH_F2PY
    use amrex_paralleldescriptor_module, only: parallel_bcast => amrex_pd_bcast, amrex_pd_ioprocessor
#endif

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    real(rt)         :: dth, dt2, dti, dt2i
    real(rt)         :: dd, dd2, ddi, dd2i
    real(rt)         :: tsav, dsav
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
    tlo   = 3.0e0_rt
    thi   = 13.0e0_rt
    tstp  = (thi - tlo)/float(jmax-1)
    tstpi = 1.0e0_rt/tstp
    dlo   = -12.0e0_rt
    dhi   = 15.0e0_rt
    dstp  = (dhi - dlo)/float(imax-1)
    dstpi = 1.0e0_rt/dstp

    do j=1,jmax
       tsav = tlo + (j-1)*tstp
       t(j) = 10.0e0_rt**(tsav)
       do i=1,imax
          dsav = dlo + (i-1)*dstp
          d(i) = 10.0e0_rt**(dsav)
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
             read(2,*) f(1,i,j),f(4,i,j),f(2,i,j),f(5,i,j),f(3,i,j),f(6,i,j), &
                       f(7,i,j),f(8,i,j),f(9,i,j)
          end do
       end do

       !..   read the pressure derivative with density table
       do j = 1, jmax
          do i = 1, imax
             read(2,*) dpdf(1,i,j), dpdf(3,i,j), dpdf(2,i,j), dpdf(4,i,j)
          end do
       end do

       !..   read the electron chemical potential table
       do j = 1, jmax
          do i = 1, imax
             read(2,*) ef(1,i,j), ef(3,i,j), ef(2,i,j), ef(4,i,j)
          end do
       end do

       !..   read the number density table
       do j = 1, jmax
          do i = 1, imax
             read(2,*) xf(1,i,j), xf(3,i,j), xf(2,i,j), xf(4,i,j)
          end do
       end do
#ifndef COMPILE_WITH_F2PY
    end if
#endif

#ifndef COMPILE_WITH_F2PY
    call parallel_bcast(f)
    call parallel_bcast(dpdf)
    call parallel_bcast(ef)
    call parallel_bcast(xf)
#endif

    !..   construct the temperature and density deltas and their inverses
    do j = 1, jmax-1
       dth         = t(j+1) - t(j)
       dt2         = dth * dth
       dti         = 1.0e0_rt/dth
       dt2i        = 1.0e0_rt/dt2
       dt_sav(j)   = dth
       dt2_sav(j)  = dt2
       dti_sav(j)  = dti
       dt2i_sav(j) = dt2i
    end do
    do i = 1, imax-1
       dd          = d(i+1) - d(i)
       dd2         = dd * dd
       ddi         = 1.0e0_rt/dd
       dd2i        = 1.0e0_rt/dd2
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

    mintemp = 10.e0_rt**tlo
    maxtemp = 10.e0_rt**thi
    mindens = 10.e0_rt**dlo
    maxdens = 10.e0_rt**dhi

    !$acc update device(mintemp, maxtemp, mindens, maxdens)

    !$acc update &
    !$acc device(tlo, thi, dlo, dhi) &
    !$acc device(tstp, tstpi, dstp, dstpi) &
    !$acc device(itmax, jtmax, d, t) &
    !$acc device(f, dpdf, ef, xf) &
    !$acc device(dt_sav, dt2_sav, dti_sav, dt2i_sav) &
    !$acc device(dd_sav, dd2_sav, ddi_sav, dd2i_sav) &
    !$acc device(do_coulomb, input_is_constant)

  end subroutine actual_eos_init



  ! quintic hermite polynomial functions
  ! psi0 and its derivatives
  pure function psi0(z) result(psi0r)
    !$acc routine seq
    real(rt)        , intent(in) :: z
    real(rt)         :: psi0r
    !$gpu
    psi0r = z**3 * ( z * (-6.0e0_rt*z + 15.0e0_rt) -10.0e0_rt) + 1.0e0_rt
  end function psi0

  pure function dpsi0(z) result(dpsi0r)
    !$acc routine seq
    real(rt)        , intent(in) :: z
    real(rt)         :: dpsi0r
    !$gpu
    dpsi0r = z**2 * ( z * (-30.0e0_rt*z + 60.0e0_rt) - 30.0e0_rt)
  end function dpsi0

  pure function ddpsi0(z) result(ddpsi0r)
    !$acc routine seq
    real(rt)        , intent(in) :: z
    real(rt)         :: ddpsi0r
    !$gpu
    ddpsi0r = z* ( z*( -120.0e0_rt*z + 180.0e0_rt) -60.0e0_rt)
  end function ddpsi0

  ! psi1 and its derivatives
  pure function psi1(z) result(psi1r)
    !$acc routine seq
    real(rt)        , intent(in) :: z
    real(rt)         :: psi1r
    !$gpu
    psi1r = z* ( z**2 * ( z * (-3.0e0_rt*z + 8.0e0_rt) - 6.0e0_rt) + 1.0e0_rt)
  end function psi1

  pure function dpsi1(z) result(dpsi1r)
    !$acc routine seq
    real(rt)        , intent(in) :: z
    real(rt)         :: dpsi1r
    !$gpu
    dpsi1r = z*z * ( z * (-15.0e0_rt*z + 32.0e0_rt) - 18.0e0_rt) +1.0e0_rt
  end function dpsi1

  pure function ddpsi1(z) result(ddpsi1r)
    !$acc routine seq
    real(rt)        , intent(in) :: z
    real(rt)         :: ddpsi1r
    !$gpu
    ddpsi1r = z * (z * (-60.0e0_rt*z + 96.0e0_rt) -36.0e0_rt)
  end function ddpsi1

  ! psi2  and its derivatives
  pure function psi2(z) result(psi2r)
    !$acc routine seq
    real(rt)        , intent(in) :: z
    real(rt)         :: psi2r
    !$gpu
    psi2r = 0.5e0_rt*z*z*( z* ( z * (-z + 3.0e0_rt) - 3.0e0_rt) + 1.0e0_rt)
  end function psi2

  pure function dpsi2(z) result(dpsi2r)
    !$acc routine seq
    real(rt)        , intent(in) :: z
    real(rt)         :: dpsi2r
    !$gpu
    dpsi2r = 0.5e0_rt*z*( z*(z*(-5.0e0_rt*z + 12.0e0_rt) - 9.0e0_rt) + 2.0e0_rt)
  end function dpsi2

  pure function ddpsi2(z) result(ddpsi2r)
    !$acc routine seq
    real(rt)        , intent(in) :: z
    real(rt)         :: ddpsi2r
    !$gpu
    ddpsi2r = 0.5e0_rt*(z*( z * (-20.0e0_rt*z + 36.0e0_rt) - 18.0e0_rt) + 2.0e0_rt)
  end function ddpsi2

  pure function fwt(fi, wt) result(fwtr)
    !$acc routine seq
    real(rt)        , intent(in) :: fi(36), wt(6)
    real(rt)         :: fwtr(6)

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
  pure function xpsi0(z) result(xpsi0r)
    !$acc routine seq
    real(rt)        , intent(in) :: z
    real(rt)         :: xpsi0r
    !$gpu
    xpsi0r = z * z * (2.0e0_rt*z - 3.0e0_rt) + 1.0
  end function xpsi0

  pure function xdpsi0(z) result(xdpsi0r)
    !$acc routine seq
    real(rt)        , intent(in) :: z
    real(rt)         :: xdpsi0r
    !$gpu
    xdpsi0r = z * (6.0e0_rt*z - 6.0e0_rt)
  end function xdpsi0


  ! psi1 & derivatives
  pure function xpsi1(z) result(xpsi1r)
    !$acc routine seq
    real(rt)        , intent(in) :: z
    real(rt)         :: xpsi1r
    !$gpu
    xpsi1r = z * ( z * (z - 2.0e0_rt) + 1.0e0_rt)
  end function xpsi1

  pure function xdpsi1(z) result(xdpsi1r)
    !$acc routine seq
    real(rt)        , intent(in) :: z
    real(rt)         :: xdpsi1r
    !$gpu
    xdpsi1r = z * (3.0e0_rt*z - 4.0e0_rt) + 1.0e0_rt
  end function xdpsi1



  subroutine actual_eos_finalize

    use amrex_fort_module, only : rt => amrex_real
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
