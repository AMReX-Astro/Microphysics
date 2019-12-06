module actual_conductivity_module

  use microphysics_type_module, only: rt

  implicit none

  character (len=64), public :: cond_name = "stellar"


contains

  subroutine actual_conductivity_init()

    implicit none

  end subroutine actual_conductivity_init


  subroutine actual_conductivity(state)

    use eos_type_module, only: eos_t
    use network, only : zion, aion, nspec

    implicit none

    type(eos_t), intent(inout) :: state

    !..this routine is sig99, it approximates the thermal transport coefficients.
    !..
    !..input:
    !..temp   = temperature temp (in K) = state % T
    !..den    = density den (in g/cm**3) = state % rho
    !..nspec  = number of isotopes in the composition
    !..xmass  = mass fractions of the composition = state % xn
    !..zion   = number of protons in each isotope (charge of each isotope)
    !..aion   = number of protons + neutrons in each isotope (atomic weight)
    !..pep    = electron-positron pressure (in erg/cm**3) = state % pele
    !..xne    = electron-positron number density (in 1/cm**3) = state % xne
    !..eta    = electron degeneracy parameter (chemical potential / k T) = state % eta

    !..output:
    !..orad   = radiation contribution to the opacity (in cm**2/g)
    !..ocond  = conductive contribution to the opacity (in cm**2/g)
    !..opac   = the total opacity (in cm**2/g)
    !..conductivity = thermal conductivity (in erg/cm/K/sec)
    !..
    !..declare the pass

    real(rt) :: orad, ocond, opac

    !..declare the internal variables
    integer  :: iz,i
    real(rt) :: xmu,t6,ytot1,ymass,abar,zbar,w(6),xh, &
                xhe,xz,xkc,xkap,xkb,xka,dbar,oiben1,d0log,xka1, &
                xkw,xkaz,dbar1log,dbar2log,oiben2,t4,t4r,t44,t45, &
                t46,ck1,ck3,ck2,ck4,ck5,ck6,xkcx,xkcy,xkcz,ochrs, &
                th,fact,facetax,faceta,ocompt,tcut,cutfac,xkf, &
                dlog10,zdel,zdell10,eta0,eta02,thpl,thpla,cfac1, &
                cfac2,oh,pefac,pefacl,pefacal,dnefac,wpar2,walf, &
                walf10,thx,thy,thc,xmas,ymas,wfac,cint, &
                vie,cie,tpe,yg,xrel,beta2,jy,vee,cee,ov1,ov,drel, &
                drel10,drelim,x,x1,x2, &
                alfa,beta

    !..various physical and derived constants
    !..con2 = con1*sqrt(4*pi*e*e/me)
    !..meff = hbar/(me*c)*(3*pi**2)**(1/3)
    !..weid = (pi*kerg)**2/(3*me)
    !..iec  = 4*e**4*me/(3*pi*hbar**3)
    !..xec  = hbar/kerg*e*sqrt(4*pi/me)

    real(rt), parameter :: third  = 1.0e0_rt/3.0e0_rt
    real(rt), parameter :: twoth  = 2.0e0_rt * third
    real(rt), parameter :: pi     = 3.1415926535897932384e0_rt
    real(rt), parameter :: avo    = 6.0221367e23_rt
    real(rt), parameter :: c      = 2.99792458e10_rt
    real(rt), parameter :: ssol   = 5.67050407222e-5_rt
    real(rt), parameter :: asol   = 4.0e0_rt*ssol/c
    real(rt), parameter :: zbound = 0.1e0_rt
    real(rt), parameter :: t7peek = 1.0e20_rt
    real(rt), parameter :: k2c    = 4.0e0_rt/3.0e0_rt*asol*c
    real(rt), parameter :: meff   = 1.194648642401440e-10_rt
    real(rt), parameter :: weid   = 6.884326138694269e-5_rt
    real(rt), parameter :: iec    = 1.754582332329132e16_rt
    real(rt), parameter :: xec    = 4.309054377592449e-7_rt
    real(rt), parameter :: rt3    = 1.7320508075688772e0_rt
    real(rt), parameter :: con2   = 1.07726359439811217e-7_rt


    !..switches for the iben & christy regimes
    real(rt), parameter :: t6_switch1 = 0.5e0_rt
    real(rt), parameter :: t6_switch2 = 0.9e0_rt
    ! parameter        (t6_switch1 = 1.0e0_rt, t6_switch2 = 1.5e0_rt)

    !$gpu

    !..initialize
    opac      = 0.0e0_rt
    orad      = 0.0e0_rt
    ocond     = 0.0e0_rt
    oiben1    = 0.0e0_rt
    oiben2    = 0.0e0_rt
    ochrs     = 0.0e0_rt
    oh        = 0.0e0_rt
    ov        = 0.0e0_rt
    zbar      = 0.0e0_rt
    ytot1     = 0.0e0_rt


    !..set the composition variables
    do i=1,6
       w(i) = 0.0e0_rt
    enddo
    do i = 1,nspec
       iz      = min(3,max(1,int(zion(i))))
       ymass   = state % xn(i)/aion(i)
       w(iz)   = w(iz) + state % xn(i)
       w(iz+3) = w(iz+3) + zion(i) * zion(i) * ymass
       zbar    = zbar + zion(i) * ymass
       ytot1   = ytot1 + ymass
    enddo
    abar = 1.0e0_rt/ytot1
    zbar = zbar * abar
    t6   = state % T * 1.0e-6_rt
    xh   = w(1)
    xhe  = w(2)
    xz   = w(3)


    !..radiative section:
    !..from iben apj 196 525 1975
    if (xh .lt. 1.0e-5_rt) then
       xmu      = max(1.0e-99_rt, w(4)+w(5)+w(6)-1.0e0_rt)
       xkc      = (2.019e-4_rt*state % rho/t6**1.7e0_rt)**(2.425e0_rt)
       xkap     = 1.0e0_rt + xkc * (1.0e0_rt + xkc/24.55e0_rt)
       xkb      = 3.86e0_rt + 0.252e0_rt*sqrt(xmu) + 0.018e0_rt*xmu
       xka      = 3.437e0_rt * (1.25e0_rt + 0.488e0_rt*sqrt(xmu) + 0.092e0_rt*xmu)
       dbar     = exp(-xka + xkb*log(t6))
       oiben1   = xkap * (state % rho/dbar)**(0.67e0_rt)
    end if

    if ( .not.((xh.ge.1.0e-5_rt) .and. (t6.lt.t6_switch1)) .and. &
         .not.((xh.lt.1.0e-5_rt) .and. (xz.gt.zbound)) ) then
       if (t6 .gt. t6_switch1) then
          d0log = -(3.868e0_rt + 0.806e0_rt*xh) + 1.8e0_rt*log(t6)
       else
          d0log = -(3.868e0_rt + 0.806e0_rt*xh) + (3.42e0_rt - 0.52e0_rt*xh)*log(t6)
       endif
       xka1 = 2.809e0_rt * exp(-(1.74e0_rt  - 0.755e0_rt*xh) &
            * (log10(t6) - 0.22e0_rt + 0.1375e0_rt*xh)**2)
       xkw  = 4.05e0_rt * exp(-(0.306e0_rt  - 0.04125e0_rt*xh) &
            * (log10(t6) - 0.18e0_rt + 0.1625e0_rt*xh)**2)
       xkaz = 50.0e0_rt*xz*xka1 * exp(-0.5206e0_rt*((log(state % rho)-d0log)/xkw)**2)
       dbar2log = -(4.283e0_rt + 0.7196e0_rt*xh) + 3.86e0_rt*log(t6)
       dbar1log = -5.296e0_rt + 4.833e0_rt*log(t6)
       if (dbar2log .lt. dbar1log) dbar1log = dbar2log
       oiben2   = (state % rho/exp(dbar1log))**(0.67e0_rt) * exp(xkaz)
    end if

    !..from christy apj 144 108 1966
    if ((t6.lt.t6_switch2) .and. (xh .ge. 1._rt)) then
       t4    = state % T * 1.0e-4_rt
       t4r   = sqrt(t4)
       t44   = t4**4
       t45   = t44 * t4
       t46   = t45 * t4
       ck1   = 2.0e6_rt/t44 + 2.1e0_rt*t46
       ck3   = 4.0e-3_rt/t44 + 2.0e-4_rt/(state % rho)**(0.25e0_rt)
       ck2   = 4.5e0_rt*t46 + 1.0e0_rt/(t4*ck3)
       ck4   = 1.4e3_rt*t4 + t46
       ck5   = 1.0e6_rt + 0.1e0_rt*t46
       ck6   = 20.0e0_rt*t4 + 5.0e0_rt*t44 + t45
       xkcx  = xh*(t4r/ck1 + 1.0e0_rt/ck2)
       xkcy  = xhe*(1.0e0_rt/ck4 + 1.5e0_rt/ck5)
       xkcz  = xz*(t4r/ck6)
       ochrs = state % pele * (xkcx + xkcy + xkcz)
    end if

    !..opacity in presence of hydrogen
    if (xh .ge. 1.0e-5_rt) then
       if (t6 .lt. t6_switch1) then
          orad   = ochrs
       else if (t6 .le. t6_switch2) then
          orad   = 2.0e0_rt*(ochrs*(1.5e0_rt - t6) + oiben2*(t6 - 1.0e0_rt))
       else
          orad   = oiben2
       end if

       !..opacity in absence of hydrogen
    else
       if (xz .gt. zbound) then
          orad   = oiben1
       else
          orad   = oiben1*(xz/zbound) + oiben2*((zbound-xz)/zbound)
       end if
    end if

    !..add in the compton scattering opacity, weaver et al. apj 1978 225 1021
    th      = min(511.0e0_rt, state % T * 8.617e-8_rt)
    fact    = 1.0e0_rt + 2.75e-2_rt*th - 4.88e-5_rt*th*th
    facetax = 1.0e100_rt
    if (state % eta .le. 500.0_rt) facetax = exp(0.522e0_rt*state % eta - 1.563e0_rt)
    faceta  = 1.0e0_rt + facetax
    ocompt  = 6.65205e-25_rt/(fact * faceta) * state % xne/state % rho
    orad    = orad   + ocompt

    !..cutoff radiative opacity when 4kt/hbar is less than the plasma
    !..frequency
    tcut = con2 * sqrt(state % xne)
    if (state % T .lt. tcut) then
       if (tcut .gt. 200.0_rt*state % T) then
          orad   = orad * 2.658e86_rt
       else
          cutfac   = exp(tcut/state % T - 1.0e0_rt)
          orad     = orad * cutfac
       end if
    end if

    !..fudge molecular opacity for low temps
    xkf    = t7peek * state % rho * (state % T * 1.0e-7_rt)**4
    orad   = xkf * orad/(xkf + orad)


    !..conductivity section:
    !..drel is the dividing line between nondegenerate and degenerate regions,
    !..taken from clayton eq. 2-34. if the density is larger than drel, then
    !..use the degenerate expressions. if the density is smaller than
    !..drelim, use the non-degenerate formulas. in between drel and drelim,
    !..apply a smooth blending of the two.

    dlog10   = log10(state % rho)

    drel     =  2.4e-7_rt * zbar/abar * state % T * sqrt(state % T)
    if (state % T .le. 1.0e5_rt) drel = drel * 15.0e0_rt
    drel10   =  log10(drel)
    drelim   =  drel10 + 1.0e0_rt


    !..from iben apj 196 525 1975 for non-degenerate regimes
    if (dlog10 .lt. drelim) then
       zdel    = state % xne/(avo*t6*sqrt(t6))
       zdell10 = log10(zdel)
       eta0    = exp(-1.20322e0_rt + twoth * log(zdel))
       eta02   = eta0*eta0

       !..thpl factor
       if (zdell10 .lt. 0.645_rt) then
          thpl    = -7.5668e0_rt + log(zdel * (1.0e0_rt + 0.024417e0_rt*zdel))
       else
          if (zdell10 .lt. 2.5_rt) then
             thpl   = -7.58110e0_rt + log(zdel*(1.0e0_rt + 0.02804e0_rt*zdel))
             if (zdell10 .ge. 2.0_rt) then
                thpla = thpl
                thpl  = -11.0742e0_rt + log(zdel**2 * (1.0e0_rt + 9.376e0_rt/eta02))
                thpl  = 2.0e0_rt*((2.5e0_rt-zdell10)*thpla + (zdell10-2.0e0_rt)*thpl)
             end if
          else
             thpl   = -11.0742e0_rt + log(zdel**2 * (1.0e0_rt + 9.376e0_rt/eta02))
          end if
       end if

       !..pefac and walf factors
       if (zdell10 .lt. 2.0_rt) then
          pefac   = 1.0e0_rt + 0.021876e0_rt*zdel
          if (zdell10 .gt. 1.5_rt) then
             pefacal   = log(pefac)
             pefacl    = log(0.4e0_rt * eta0 + 1.64496e0_rt/eta0)
             cfac1     = 2.0e0_rt - zdell10
             cfac2     = zdell10 - 1.5e0_rt
             pefac     = exp(2.0e0_rt * (cfac1*pefacal + cfac2*pefacl))
          end if
       else
          pefac   = 0.4e0_rt * eta0 + 1.64496e0_rt/eta0
       end if

       if (zdel.lt.40.0_rt) then
          dnefac = 1.0e0_rt + zdel * (3.4838e-4_rt * zdel - 2.8966e-2_rt)
       else
          dnefac = 1.5e0_rt/eta0 * (1.0e0_rt - 0.8225e0_rt/eta02)
       endif
       wpar2  = 9.24735e-3_rt * zdel * &
            (state % rho*avo*(w(4)+w(5)+w(6))/state % xne + dnefac)/(sqrt(t6)*pefac)
       walf   = 0.5e0_rt * log(wpar2)
       walf10 = 0.5e0_rt * log10(wpar2)

       !..thx, thy and thc factors
       if (walf10 .le. -3.0_rt) then
          thx   = exp(2.413e0_rt - 0.124e0_rt*walf)
       else if (walf10 .le. -1.0_rt) then
          thx   = exp(0.299e0_rt - walf*(0.745e0_rt + 0.0456e0_rt*walf))
       else
          thx   = exp(0.426e0_rt - 0.558e0_rt*walf)
       end if

       if (walf10 .le. -3.0_rt) then
          thy   = exp(2.158e0_rt - 0.111e0_rt*walf)
       else if (walf10 .le. 0.0_rt) then
          thy   = exp(0.553e0_rt - walf*(0.55e0_rt + 0.0299e0_rt*walf))
       else
          thy   = exp(0.553e0_rt - 0.6e0_rt*walf)
       end if

       if (walf10 .le. -2.5_rt) then
          thc   = exp(2.924e0_rt - 0.1e0_rt*walf)
       else if (walf10 .le. 0.5_rt) then
          thc   = exp(1.6740e0_rt - walf*(0.511e0_rt + 0.0338e0_rt*walf))
       else
          thc   = exp(1.941e0_rt - 0.785e0_rt*walf)
       end if

       oh   = (xh*thx + xhe*thy + w(6)*third*thc) / (t6*exp(thpl))
    end if

    !..from yakovlev & urpin soviet astro 1980 24 303 and
    !..potekhin et al. 1997 aa 323 415 for degenerate regimes
    if (dlog10 .gt. drel10) then
       xmas   = meff * (state % xne)**third
       ymas   = sqrt(1.0e0_rt + xmas*xmas)
       wfac   = weid * state % T/ymas * state % xne
       cint   = 1.0e0_rt

       !..ion-electron collision frequency and the thermal conductivity
       vie   = iec * zbar * ymas * cint
       cie   = wfac/vie

       !..electron-electron collision frequency and thermal conductivity
       tpe  = xec * sqrt(state % xne/ymas)
       yg   = rt3 * tpe/state % T
       xrel = 1.009e0_rt * (zbar/abar * state % rho * 1.0e-6_rt)**third
       beta2 = xrel**2/(1.0e0_rt + xrel**2)
       jy   = (1.0e0_rt + 6.0e0_rt/(5.0e0_rt*xrel**2) + 2.0e0_rt/(5.0e0_rt*xrel**4)) &
            * ( yg**3 / (3.0e0_rt * (1.0e0_rt + 0.07414_rt * yg)**3) &
            * log((2.81e0_rt - 0.810_rt*beta2 + yg)/yg)  &
            + pi**5/6.0e0_rt * (yg/(13.91e0_rt + yg))**4 )
       vee = 0.511e0_rt * (state % T)**2 * xmas/ymas**2 * sqrt(xmas/ymas) * jy
       cee = wfac/vee

       !..total electron thermal conductivity and conversion to an opacity
       ov1   = cie * cee/(cee + cie)
       ov    = k2c/(ov1*state % rho) * (state % T)**3
    end if

    !..blend the opacities in the intermediate region
    if (dlog10 .le. drel10) then
       ocond   = oh
    else if (dlog10 .gt. drel10  .and. dlog10 .lt. drelim) then
       x        = state % rho
       x1       = 10.0e0_rt**drel10
       x2       = 10.0e0_rt**drelim
       alfa     = (x-x2)/(x1-x2)
       beta     = (x-x1)/(x2-x1)
       ocond    = alfa*oh + beta*ov
    else if (dlog10 .ge. drelim) then
       ocond   = ov
    end if

    !..total opacity
    opac    = orad * ocond / (ocond + orad)

    state % conductivity = k2c * (state % T)**3 / (opac * state % rho)

    return
  end subroutine actual_conductivity

end module actual_conductivity_module
