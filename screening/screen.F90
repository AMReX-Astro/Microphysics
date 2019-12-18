module screening_module

  use amrex_constants_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer :: nscreen = 0

  real(rt)        , parameter :: fact       = 1.25992104989487e0_rt
  real(rt)        , parameter :: co2        = THIRD * 4.248719e3_rt
  real(rt)        , parameter :: gamefx     = 0.3e0_rt
  real(rt)        , parameter :: gamefs     = 0.8e0_rt
  real(rt)        , parameter :: h12_max    = 300.e0_rt

  real(rt)        , allocatable, save :: z1scr(:)
  real(rt)        , allocatable, save :: z2scr(:)
  real(rt)        , allocatable, save :: a1scr(:)
  real(rt)        , allocatable, save :: a2scr(:)

  ! zs13    = (z1+z2)**(1./3.)
  ! zhat    = combination of z1 and z2 raised to the 5/3 power
  ! zhat2   = combination of z1 and z2 raised to the 5/12 power
  ! lzav    = log of effective charge
  ! aznut   = combination of a1,z1,a2,z2 raised to 1/3 power

  real(rt)        , allocatable, save :: zs13(:)
  real(rt)        , allocatable, save :: zs13inv(:)
  real(rt)        , allocatable, save :: zhat(:)
  real(rt)        , allocatable, save :: zhat2(:)
  real(rt)        , allocatable, save :: lzav(:)
  real(rt)        , allocatable, save :: aznut(:)

  private :: fact, co2, gamefx, gamefs
  private :: z1scr, z2scr, a1scr, a2scr
  private :: zs13, zs13inv, zhat, zhat2, lzav, aznut

  type :: plasma_state

     real(rt)         :: qlam0z
     real(rt)         :: qlam0zdt
     real(rt)         :: qlam0zdd

     real(rt)         :: taufac
     real(rt)         :: taufacdt

     real(rt)         :: aa
     real(rt)         :: daadt
     real(rt)         :: daadd

  end type plasma_state

#ifdef AMREX_USE_CUDA
  attributes(managed) :: z1scr, z2scr, a1scr, a2scr
  attributes(managed) :: zs13, zs13inv, zhat, zhat2, lzav, aznut
#endif

  !$acc declare &
  !$acc create(nscreen) &
  !$acc create(fact, co2, gamefx, gamefs, h12_max) &
  !$acc create(z1scr, z2scr, a1scr, a2scr) &
  !$acc create(zs13, zs13inv, zhat, zhat2, lzav, aznut)

contains

  subroutine screening_init()

    implicit none

    integer :: i

    ! This routine assumes that we have already filled z1scr, z2scr,
    ! a1scr, and a2scr.

    allocate(zs13(nscreen))
    allocate(zs13inv(nscreen))
    allocate(zhat(nscreen))
    allocate(zhat2(nscreen))
    allocate(lzav(nscreen))
    allocate(aznut(nscreen))

    do i = 1, nscreen

       zs13(i)    = (z1scr(i) + z2scr(i))**THIRD
       zs13inv(i) = ONE/zs13(i)
       zhat(i)    = (z1scr(i) + z2scr(i))**FIVE3RD  - z1scr(i)**FIVE3RD - z2scr(i)**FIVE3RD
       zhat2(i)   = (z1scr(i) + z2scr(i))**FIVE12TH - z1scr(i)**FIVE12TH -z2scr(i)**FIVE12TH
       lzav(i)    = FIVE3RD * log(z1scr(i)*z2scr(i)/(z1scr(i) + z2scr(i)))
       aznut(i)   = (z1scr(i)**2 * z2scr(i)**2 * a1scr(i)*a2scr(i) / (a1scr(i) + a2scr(i)))**THIRD

    enddo

    !$acc update device(zs13, zs13inv, zhat, zhat2, lzav, aznut)

  end subroutine screening_init

  subroutine screening_finalize()

    implicit none

    ! Deallocate the screening buffers.

    if (allocated(z1scr)) then
       deallocate(z1scr)
    end if

    if (allocated(z2scr)) then
       deallocate(z2scr)
    end if

    if (allocated(a1scr)) then
       deallocate(a1scr)
    end if

    if (allocated(a2scr)) then
       deallocate(a2scr)
    end if

    if (allocated(zs13)) then
       deallocate(zs13)
    end if

    if (allocated(zs13inv)) then
       deallocate(zs13inv)
    end if

    if (allocated(zhat)) then
       deallocate(zhat)
    end if

    if (allocated(zhat2)) then
       deallocate(zhat2)
    end if

    if (allocated(lzav)) then
       deallocate(lzav)
    end if

    if (allocated(aznut)) then
       deallocate(aznut)
    end if

  end subroutine screening_finalize

  subroutine add_screening_factor(z1, a1, z2, a2)

    ! this is only called at initialization

    real(rt)         :: z1, a1, z2, a2

    real(rt)        , allocatable :: z1scr_temp(:), a1scr_temp(:)
    real(rt)        , allocatable :: z2scr_temp(:), a2scr_temp(:)

    ! Deallocate the buffers, then reallocate
    ! them with a size one larger. This is
    ! admittedly wasteful, but since this routine
    ! only happens at initialization the cost
    ! does not matter, and it's important for keeping
    ! memory usage down that the screening size is
    ! only as large as it needs to be.

    if (nscreen > 0) then

       allocate(z1scr_temp(nscreen))
       z1scr_temp(:) = z1scr(:)
       deallocate(z1scr)

       allocate(z2scr_temp(nscreen))
       z2scr_temp(:) = z2scr(:)
       deallocate(z2scr)

       allocate(a1scr_temp(nscreen))
       a1scr_temp(:) = a1scr(:)
       deallocate(a1scr)

       allocate(a2scr_temp(nscreen))
       a2scr_temp(:) = a2scr(:)
       deallocate(a2scr)

    else

       allocate(z1scr(1))
       allocate(z2scr(1))
       allocate(a1scr(1))
       allocate(a2scr(1))

    end if

    nscreen = nscreen + 1

    if (nscreen > 1) then

       allocate(z1scr(nscreen))
       z1scr(1:nscreen-1) = z1scr_temp(:)
       deallocate(z1scr_temp)

       allocate(z2scr(nscreen))
       z2scr(1:nscreen-1) = z2scr_temp(:)
       deallocate(z2scr_temp)

       allocate(a1scr(nscreen))
       a1scr(1:nscreen-1) = a1scr_temp(:)
       deallocate(a1scr_temp)

       allocate(a2scr(nscreen))
       a2scr(1:nscreen-1) = a2scr_temp(:)
       deallocate(a2scr_temp)

    end if

    z1scr(nscreen) = z1
    a1scr(nscreen) = a1
    z2scr(nscreen) = z2
    a2scr(nscreen) = a2

    !$acc update device(nscreen, z1scr, a1scr, z2scr, a2scr)

  end subroutine add_screening_factor


  subroutine fill_plasma_state(state, temp, dens, y)

    !$acc routine seq

    use network, only: nspec, zion

    ! Input variables

    type (plasma_state) :: state
    real(rt)         :: temp, dens, y(nspec)

    ! Local variables

    real(rt)         :: abar, zbar, z2bar
    real(rt)         :: ytot, rr, tempi, dtempi, deni
    real(rt)         :: pp, qq, dppdt, xni
!    real(rt)         :: dppdd

    !$gpu

    abar   = ONE / sum(y)
    zbar   = sum(zion * y) * abar
    z2bar  = sum(zion**2 * y) * abar

    ytot             = ONE / abar
    rr               = dens * ytot
    tempi            = ONE / temp
    dtempi           = -tempi * tempi
    deni             = ONE / dens

    pp               = sqrt(rr*tempi*(z2bar + zbar))
    qq               = HALF/pp *(z2bar + zbar)
    dppdt            = qq*rr*dtempi
    !dppdd            = qq * ytot * tempi

    state % qlam0z   = 1.88e8_rt * tempi * pp
    state % qlam0zdt = 1.88e8_rt * (dtempi*pp + tempi*dppdt)
    !state % qlam0zdd = 1.88e8_rt * tempi * dppdd

    state % taufac   = co2 * tempi**THIRD
    state % taufacdt = -THIRD * state % taufac * tempi

    qq               = rr * zbar
    xni              = qq**THIRD
    !dxnidd           = THIRD * xni * deni

    state % aa       = 2.27493e5_rt * tempi * xni
    state % daadt    = 2.27493e5_rt * dtempi * xni
    !state % daadd    = 2.27493e5_rt * tempi * dxnidd

  end subroutine fill_plasma_state


  subroutine screen5(state,jscreen,scor,scordt,scordd)

    !$acc routine seq

    use amrex_constants_module, only: M_PI

    implicit none

    ! this subroutine calculates screening factors and their derivatives
    ! for nuclear reaction rates in the weak, intermediate and strong regimes.
    ! based on graboske, dewit, grossman and cooper apj 181 457 1973 for
    ! weak screening. based on alastuey and jancovici apj 226 1034 1978,
    ! with plasma parameters from itoh et al apj 234 1079 1979, for strong
    ! screening.

    ! input:
    ! state   = plasma state (T, rho, abar, zbar, etc.)
    ! jscreen = counter of which reaction is being calculated

    ! output:
    ! scor    = screening correction
    ! scordt  = derivative of screening correction with temperature
    ! scordd  = derivative of screening correction with density


    ! declare the pass
    integer             :: jscreen
    type (plasma_state) :: state
    real(rt)            :: scor, scordt, scordd


    ! local variables
    real(rt)         :: z1, a1, z2, a2

    real(rt)         :: bb,cc,dccdt, &
                        qq,dqqdt,rr,drrdt, &
                        ss,dssdt,tt,dttdt,uu,duudt, &
                        vv,dvvdt,a3,da3, &
                        h12w,dh12wdt,h12,dh12dt, &
                        gamp,gampdt, &
                        gamef,gamefdt, &
                        tau12,tau12dt,alph12,alph12dt, &
                        xlgfac,dxlgfacdt, &
                        gamp14,gamp14dt,dgamma

!    real(rt)         :: dccdd,dqqdd,dvvdd,drrdd,dssdd,dttdd,duudd
!    real(rt)         :: dh12dd,dh12wdd,dh12xdd,alph12dd
!    real(rt)         :: gampdd,gamefdd,dxlgcfacdd,gamp14dd

    !$gpu

    ! Get the ion data based on the input index

    z1 = z1scr(jscreen)
    a1 = a1scr(jscreen)
    z2 = z2scr(jscreen)
    a2 = a2scr(jscreen)

    ! calculate individual screening factors
    bb       = z1 * z2
    gamp     = state % aa
    gampdt   = state % daadt
    !gampdd   = state % daadd

    qq       = fact * bb * zs13inv(jscreen)
    gamef    = qq * gamp
    gamefdt  = qq * gampdt
    !gamefdd  = qq * gampdd

    tau12    = state % taufac * aznut(jscreen)
    tau12dt  = state % taufacdt * aznut(jscreen)

    qq       = ONE/tau12
    alph12   = gamef * qq
    alph12dt = (gamefdt - alph12*tau12dt) * qq
    !alph12dd = gamefdd * qq



    ! limit alph12 to 1.6 to prevent unphysical behavior.
    ! this should really be replaced by a pycnonuclear reaction rate formula
    if (alph12 .gt. 1.6_rt) then
       alph12   = 1.6e0_rt
       alph12dt = ZERO
       !alph12dd = ZERO

       gamef    = 1.6e0_rt * tau12
       gamefdt  = 1.6e0_rt * tau12dt
       !gamefdd  = ZERO

       qq       = zs13(jscreen)/(fact * bb)
       gamp     = gamef * qq
       gampdt   = gamefdt * qq
       !gampdd   = ZERO
    end if



    ! weak screening regime
    h12w    = bb * state % qlam0z
    dh12wdt = bb * state % qlam0zdt
    !dh12wdd = bb * qlam0zdd

    h12     = h12w
    dh12dt  = dh12wdt
    !dh12dd  = dh12wdd



    ! intermediate and strong sceening regime
    if (gamef .gt. gamefx) then

       gamp14   = gamp**FOURTH
       rr       = ONE/gamp
       qq       = 0.25e0_rt*gamp14*rr
       gamp14dt = qq * gampdt
       !gamp14dd = qq * gampdd

       cc       =   0.896434e0_rt * gamp * zhat(jscreen) &
            - 3.44740e0_rt  * gamp14 * zhat2(jscreen) &
            - 0.5551e0_rt   * (log(gamp) + lzav(jscreen)) &
            - 2.996e0_rt

       dccdt    =   0.896434e0_rt * gampdt * zhat(jscreen) &
            - 3.44740e0_rt  * gamp14dt * zhat2(jscreen) &
            - 0.5551e0_rt*rr*gampdt

       !dccdd    =   0.896434e0_rt * gampdd * zhat(jscreen) &
       !     - 3.44740e0_rt  * gamp14dd * zhat2(jscreen) &
       !     - 0.5551e0_rt*rr*gampdd

       a3     = alph12 * alph12 * alph12
       da3    = 3.0e0_rt * alph12 * alph12

       qq     = 0.014e0_rt + 0.0128e0_rt*alph12
       dqqdt  = 0.0128e0_rt*alph12dt
       !dqqdd  = 0.0128e0_rt*alph12dd

       rr     = FIVE32ND - alph12*qq
       drrdt  = -(alph12dt*qq + alph12*dqqdt)
       !drrdd  = -(alph12dd*qq + alph12*dqqdd)

       ss     = tau12*rr
       dssdt  = tau12dt*rr + tau12*drrdt
       !dssdd  = tau12*drrdd

       tt     =  -0.0098e0_rt + 0.0048e0_rt*alph12
       dttdt  = 0.0048e0_rt*alph12dt
       !dttdd  = 0.0048e0_rt*alph12dd

       uu     =  0.0055e0_rt + alph12*tt
       duudt  = alph12dt*tt + alph12*dttdt
       !duudd  = alph12dd*tt + alph12*dttdd

       vv   = gamef * alph12 * uu
       dvvdt= gamefdt*alph12*uu + gamef*alph12dt*uu + gamef*alph12*duudt
       !dvvdd= gamefdd*alph12*uu + gamef*alph12dd*uu + gamef*alph12*duudd

       h12     = cc - a3 * (ss + vv)
       rr      = da3 * (ss + vv)
       dh12dt  = dccdt - rr*alph12dt - a3*(dssdt + dvvdt)
       !dh12dd  = dccdd - rr*alph12dd - a3*(dssdd + dvvdd)

       rr     =  ONE - 0.0562e0_rt*a3
       ss     =  -0.0562e0_rt*da3
       drrdt  = ss*alph12dt
       !drrdd  = ss*alph12dd

       if (rr .ge. 0.77e0_rt) then
          xlgfac    = rr
          dxlgfacdt = drrdt
          !dxlgfacdd = drrdd
       else
          xlgfac    = 0.77e0_rt
          dxlgfacdt = ZERO
          !dxlgfacdd = ZERO
       end if


       h12    = log(xlgfac) + h12
       rr     = ONE/xlgfac
       dh12dt = rr*dxlgfacdt + dh12dt
       !dh12dd = rr*dxlgfacdd + dh12dd

       if (gamef .le. gamefs) then
          dgamma  = 1.0e0_rt/(gamefs - gamefx)

          rr     =  dgamma*(gamefs - gamef)
          drrdt  = -dgamma*gamefdt
          !drrdd  = -dgamma*gamefdd

          ss     = dgamma*(gamef - gamefx)
          dssdt  = dgamma*gamefdt
          !dssdd  = dgamma*gamefdd

          vv     = h12

          h12    = h12w*rr + vv*ss
          dh12dt = dh12wdt*rr + h12w*drrdt + dh12dt*ss + vv*dssdt
          !dh12dd = dh12wdd*rr + h12w*drrdd + dh12dd*ss + vv*dssdd
       end if


       ! end of intermediate and strong screening if
    end if


    ! machine limit the output
    ! further limit to avoid the pycnonuclear regime
    h12    = max(min(h12, h12_max), ZERO)
    scor   = exp(h12)
    if (h12 .eq. h12_max) then
       scordt = ZERO
       !scordd = ZERO
    else
       scordt = scor * dh12dt
       !scordd = scor * dh12dd
    end if

  end subroutine screen5


  subroutine screenz (t,d,z1,z2,a1,a2,ymass,scfac,dscfacdt)

    !$acc routine seq

    use network, only: aion, zion, nspec

    implicit none

    real(rt)         :: t, d, z1, z2, a1, a2
    real(rt)         :: ymass(nspec)
    real(rt)         :: scfac
    real(rt)         :: dscfacdt

    ! this subroutine calculates screening factors for nuclear reaction
    ! rates in the weak, intermediate , and strong regimes given the
    ! temperature (t--degk), the density (d--g/cc), the atomic numbers
    ! and weights of the elements in the reaction channel with the
    ! largest coulomb barrier (z1,z2,a1,a2), and the mean plasma
    ! parameters calculated in main and passed over in common aver:
    ! (mean atomic number--zbar, mean square of the atomic
    ! number--z2bar, mean atomic weight--abar, and total number of moles
    ! of nuclei per gram--ytot1).  the unscreened rate is to be
    ! multiplied by the dimensionless the treatment is based on
    ! graboske, dewit, grossman, and cooper ap j. 181,457 (1973) for
    ! weak screening and on alastuey and jancovici, ap.j. 226, 1034,
    ! 1978, with plasma parameters from itoh, totsuji, setsuo, and
    ! dewitt, ap.j. 234, 1079,1979, for strong screening (rkw
    ! modification).

    !.... last revision 15 nov 1982

    real(rt)         abar, zbar, ytot1, z2bar, theta

    integer iy

    real(rt)         qlam0, ztilda, qlam0z, gamp, taufac
    real(rt)         dqlam0dt, dqlam0zdt, dgampdt, dtaufacdt
    real(rt)         zhat, zhat2, gamef, tau12, alph12
    real(rt)         dgamefdt, dtau12dt, dalph12dt
    real(rt)         h12w, h12, c, h12fac
    real(rt)         dh12wdt, dh12dt, dcdt

    !$gpu

    ! calculate averages for screening routine
    ! nb  y = x/a with x the mass fraction
    ! zi and ai are the nuclear chage and atomic mass number
    ! respectively

    ! this part came in through a common block in Kepler -- do it
    ! directly here
    abar=0.e0_rt
    zbar=0.e0_rt
    ytot1=0.e0_rt
    z2bar=0.e0_rt

    do iy = 1, nspec
       ytot1 = ytot1 + ymass(iy)
       z2bar = z2bar + zion(iy)**2*ymass(iy)
       abar = abar + aion(iy)*ymass(iy)
       zbar = zbar + zion(iy)*ymass(iy)
    enddo

    z2bar=z2bar/ytot1
    abar=abar/ytot1
    zbar=zbar/ytot1

    ! resume original Kepler screen...
    theta=1.e0_rt
    ytot1=1.e0_rt/abar

    !.... calculate average plasma parameters
    !....
    if ((z1*z2) > 0.e0_rt) then

       qlam0=1.88e+8_rt*sqrt(d/(abar*t**3))
       dqlam0dt=-1.5e+0_rt*qlam0 / t

       ztilda=sqrt(z2bar+zbar*theta)

       qlam0z=qlam0*ztilda
       dqlam0zdt=ztilda*dqlam0dt

       gamp=2.27493e+5_rt*(d*zbar*ytot1)**THIRD/t
       dgampdt=-gamp/t

       taufac=4.248719e+3_rt/t**THIRD
       dtaufacdt=-THIRD*taufac/t

       !.... calculate screening factor
       !.... approx. for strong screening only good for alpha .lt. 1.6

       zhat=(z1+z2)**FIVE3RD-z1**FIVE3RD-z2**FIVE3RD
       zhat2=(z1+z2)**FIVE12TH-z1**FIVE12TH-z2**FIVE12TH

       gamef=2.e0_rt**THIRD*gamp*z1*z2/(z1+z2)**THIRD
       dgamefdt=gamef*dgampdt/gamp

       tau12=taufac*(z1**2*z2**2*a1*a2/(a1+a2))**THIRD
       dtau12dt=tau12*dtaufacdt/taufac

       alph12=3.e0_rt*gamef/tau12
       dalph12dt=alph12*(dgamefdt/gamef - dtau12dt/tau12)

       !....
       !.... limit alph12 to 1.6 to prevent unphysical behavior
       !.... (h dec. as rho inc.) at high rho.  this should really
       !.... be replaced by a pycnonuclear reaction rate formula.
       !....
       if (alph12 > 1.6e0_rt) then

          alph12=1.6e0_rt
          dalph12dt=0.0e0_rt

          gamef=1.6e0_rt*tau12/3.e0_rt
          dgamefdt=gamef*dtau12dt/tau12

          gamp=gamef*(z1+z2)**THIRD/(2.e0_rt**THIRD*z1*z2)
          dgampdt=gamp*dgamefdt/gamef

       endif

       h12w=z1*z2*qlam0z
       dh12wdt=h12w*dqlam0zdt/qlam0z

       h12=h12w
       dh12dt=dh12wdt

       if (gamef > 0.3e0_rt) then

          c=0.896434e0_rt*gamp*zhat-3.44740e0_rt*gamp**FOURTH*zhat2- &
               0.5551e0_rt*(log(gamp)+FIVE3RD*log(z1*z2/(z1+z2)))-2.996e0_rt

          dcdt=0.896434e0_rt*dgampdt*zhat- &
               3.44740e0_rt*FOURTH*gamp**(FOURTH-1.0e0_rt)*zhat2*dgampdt- &
               0.5551e0_rt*dgampdt/gamp

          h12=c-(tau12/3.e0_rt)*(5.e0_rt*alph12**3/32.e0_rt-0.014e0_rt*alph12**4 &
               -0.0128e0_rt*alph12**5)-gamef*(0.0055e0_rt*alph12**4 &
               -0.0098e0_rt*alph12**5+0.0048e0_rt*alph12**6)

          dh12dt=dcdt - ((dtau12dt*alph12**3 + 3.0e0_rt*tau12*alph12**2* &
               dalph12dt)*(5.e0_rt/32.e0_rt - 0.014e0_rt*alph12 - &
               0.0128e0_rt*alph12**2) + tau12*alph12**3*dalph12dt*(-0.014e0_rt &
               - 2.e0_rt*0.0128e0_rt*alph12))/3.e0_rt -(dgamefdt*alph12**4 + 4.e0_rt &
               *gamef*alph12**3*dalph12dt)*(0.0055e0_rt - 0.0098e0_rt*alph12 - &
               0.0048e0_rt*alph12**2) - gamef*alph12**4*dalph12dt*(-0.0098e0_rt &
               + 2.e0_rt*0.0048e0_rt*alph12)

          h12fac=0.77e0_rt

          h12=log(max(1.e+0_rt-0.0562e+0_rt*alph12**3,h12fac))+h12
          if (1.e+0_rt-0.0562e+0_rt*alph12**3 .gt. h12fac) then
             dh12dt=(-3.e0_rt*0.0562e0_rt*alph12**2*dalph12dt)/ &
                  (1.e0_rt-0.0562e0_rt*alph12**3) + dh12dt
          endif

          if(gamef <= 0.8e0_rt) then

             h12=h12w*((0.8e0_rt-gamef)/0.5e0_rt)+h12*((gamef-0.3e0_rt)/0.5e0_rt)
             dh12dt=((dh12wdt*(0.8e0_rt-gamef) - h12w*dgamefdt + dh12dt* &
                  (gamef-0.3e0_rt) + h12*dgamefdt)/0.5e0_rt)
          endif
       endif

       if (h12 .gt. h12_max) then
          h12 = h12_max
          dh12dt=0.e0_rt
       endif

       if (h12.lt.0.e0_rt) then
          h12=0.e0_rt
          dh12dt=0.e0_rt
       endif

       scfac=exp(h12)
       dscfacdt=scfac*dh12dt

    else
       scfac=1.e0_rt
       dscfacdt=0.e0_rt
    endif

    return

  end subroutine screenz


end module screening_module
