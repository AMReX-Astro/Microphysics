module actual_conductivity_module

  use bl_types

  implicit none

contains

  subroutine actual_conductivity_init()

    implicit none

  end subroutine actual_conductivity_init


  subroutine actual_conductivity(eos_state, conductivity)

    use eos_type_module, only: eos_t
    use network, only : zion, aion, nspec

    implicit none
    
    type(eos_t), intent(in   ) :: eos_state
    real (dp_t)  , intent(inout) :: conductivity
    
    real(dp_t) ::xmass_temp(nspec)
    real(dp_t) :: orad, ocond, opac

    ! this is done to be compatible with interface to sig99
    xmass_temp(:) = eos_state%xn(:)

    call sig99(eos_state%T, eos_state%rho, xmass_temp, &
               zion, aion, nspec, &
               eos_state%pele, eos_state%xne, eos_state%eta, &
               orad, ocond, opac, conductivity)

  end subroutine actual_conductivity


  subroutine sig99(temp,den,xmass,zion,aion,ionmax,pep,xne,eta, &
                   orad,ocond,opac,conductivity)
    implicit none
    
    !..this routine approximates the thermal transport coefficients.
    !..
    !..input:
    !..temp   = temperature temp (in K)
    !..den    = density den (in g/cm**3)
    !..ionmax = number of isotopes in the composition
    !..xmass  = mass fractions of the composition
    !..zion   = number of protons in each isotope (charge of each isotope)
    !..aion   = number of protons + neutrons in each isotope (atomic weight)
    !..pep    = electron-positron pressure (in erg/cm**3)
    !..xne    = electron-positron number density (in 1/cm**3)
    !..eta    = electron degeneracy parameter (chemical potential / k T)
    
    !..output:
    !..orad   = radiation contribution to the opacity (in cm**2/g)
    !..ocond  = conductive contribution to the opacity (in cm**2/g)
    !..opac   = the total opacity (in cm**2/g)
    !..conductivity = thermal conductivity (in erg/cm/K/sec)
    !..
    !..declare the pass
    integer          ionmax
    double precision temp,den,xmass(ionmax),zion(ionmax),aion(ionmax), &
         pep,xne,eta,opac,conductivity

    !..declare the internal variables
    integer          iz,i
    double precision xmu,t6,orad,ocond,ytot1,ymass,abar,zbar,w(6),xh, &
                     xhe,xz,xkc,xkap,xkb,xka,dbar,oiben1,d0log,xka1, &
                     xkw,xkaz,dbar1log,dbar2log,oiben2,t4,t4r,t44,t45, &
                     t46,ck1,ck3,ck2,ck4,ck5,ck6,xkcx,xkcy,xkcz,ochrs, &
                     th,fact,facetax,faceta,ocompt,tcut,cutfac,xkf, &
                     dlog10,zdel,zdell10,eta0,eta02,thpl,thpla,cfac1, &
                     cfac2,oh,pefac,pefacl,pefacal,dnefac,wpar2,walf, &
                     walf10,thx,thy,thc,xmas,ymas,wfac,cint, &
                     vie,cie,tpe,yg,xrel,beta2,jy,vee,cee,ov1,ov,drel, &
                     drel10,drelim,t6_switch1,t6_switch2,x,x1,x2, &
                     alfa,beta


    !..various physical and derived constants
    !..con2 = con1*sqrt(4*pi*e*e/me)  
    !..meff = hbar/(me*c)*(3*pi**2)**(1/3) 
    !..weid = (pi*kerg)**2/(3*me)     
    !..iec  = 4*e**4*me/(3*pi*hbar**3) 
    !..xec  = hbar/kerg*e*sqrt(4*pi/me) 

    double precision third,twoth,pi,avo,c,ssol,asol,zbound,t7peek, &
                     con2,k2c,meff,weid,iec,xec,rt3
    parameter (third = 1.0d0/3.0d0, &
              twoth  = 2.0d0 * third, & 
              pi     = 3.1415926535897932384d0, &
              avo    = 6.0221367d23, &
              c      = 2.99792458d10, &
              ssol   = 5.67050407222d-5, &
              asol   = 4.0d0*ssol/c, &
              zbound = 0.1d0, &
              t7peek = 1.0d20, &
              k2c    = 4.0d0/3.0d0*asol*c, &
              meff   = 1.194648642401440d-10, &
              weid   = 6.884326138694269d-5, &
              iec    = 1.754582332329132d16, & 
              xec    = 4.309054377592449d-7, &
              rt3    = 1.7320508075688772d0, &
              con2   = 1.07726359439811217d-7)


    !..switches for the iben & christy regimes
    parameter        (t6_switch1 = 0.5d0, t6_switch2 = 0.9d0)
    ! parameter        (t6_switch1 = 1.0d0, t6_switch2 = 1.5d0)


    !..initialize
    opac      = 0.0d0
    orad      = 0.0d0
    ocond     = 0.0d0
    oiben1    = 0.0d0
    oiben2    = 0.0d0
    ochrs     = 0.0d0
    oh        = 0.0d0
    ov        = 0.0d0
    zbar      = 0.0d0
    ytot1     = 0.0d0


    !..set the composition variables 
    do i=1,6
       w(i) = 0.0d0
    enddo
    do i = 1,ionmax 
       iz      = min(3,max(1,int(zion(i))))
       ymass   = xmass(i)/aion(i)
       w(iz)   = w(iz) + xmass(i) 
       w(iz+3) = w(iz+3) + zion(i) * zion(i) * ymass
       zbar    = zbar + zion(i) * ymass
       ytot1   = ytot1 + ymass
    enddo
    abar = 1.0d0/ytot1
    zbar = zbar * abar
    t6   = temp * 1.0d-6 
    xh   = w(1)
    xhe  = w(2)
    xz   = w(3)


    !..radiative section:
    !..from iben apj 196 525 1975
    if (xh .lt. 1.0e-5) then 
       xmu      = max(1.0d-99, w(4)+w(5)+w(6)-1.0d0)    
       xkc      = (2.019e-4*den/t6**1.7d0)**(2.425d0) 
       xkap     = 1.0d0 + xkc * (1.0d0 + xkc/24.55d0)
       xkb      = 3.86d0 + 0.252d0*sqrt(xmu) + 0.018d0*xmu 
       xka      = 3.437d0 * (1.25d0 + 0.488d0*sqrt(xmu) + 0.092d0*xmu)
       dbar     = exp(-xka + xkb*log(t6))  
       oiben1   = xkap * (den/dbar)**(0.67d0)
    end if

    if ( .not.((xh.ge.1.0e-5) .and. (t6.lt.t6_switch1)) .and. &
         .not.((xh.lt.1.0e-5) .and. (xz.gt.zbound)) ) then 
       if (t6 .gt. t6_switch1) then
          d0log = -(3.868d0 + 0.806d0*xh) + 1.8d0*log(t6)
       else 
          d0log = -(3.868d0 + 0.806d0*xh) + (3.42d0 - 0.52d0*xh)*log(t6)
       endif
       xka1 = 2.809d0 * exp(-(1.74d0  - 0.755d0*xh) &
            * (log10(t6) - 0.22d0 + 0.1375d0*xh)**2)  
       xkw  = 4.05d0 * exp(-(0.306d0  - 0.04125d0*xh) &
            * (log10(t6) - 0.18d0 + 0.1625d0*xh)**2) 
       xkaz = 50.0d0*xz*xka1 * exp(-0.5206d0*((log(den)-d0log)/xkw)**2)
       dbar2log = -(4.283d0 + 0.7196d0*xh) + 3.86d0*log(t6)
       dbar1log = -5.296d0 + 4.833d0*log(t6)  
       if (dbar2log .lt. dbar1log) dbar1log = dbar2log
       oiben2   = (den/exp(dbar1log))**(0.67d0) * exp(xkaz) 
    end if

    !..from christy apj 144 108 1966
    if ((t6.lt.t6_switch2) .and. (xh .ge. 1.0e-5)) then 
       t4    = temp * 1.0d-4 
       t4r   = sqrt(t4)  
       t44   = t4**4
       t45   = t44 * t4 
       t46   = t45 * t4 
       ck1   = 2.0d6/t44 + 2.1d0*t46  
       ck3   = 4.0d-3/t44 + 2.0d-4/den**(0.25d0) 
       ck2   = 4.5d0*t46 + 1.0d0/(t4*ck3)    
       ck4   = 1.4d3*t4 + t46  
       ck5   = 1.0d6 + 0.1d0*t46  
       ck6   = 20.0d0*t4 + 5.0d0*t44 + t45  
       xkcx  = xh*(t4r/ck1 + 1.0d0/ck2)    
       xkcy  = xhe*(1.0d0/ck4 + 1.5d0/ck5)   
       xkcz  = xz*(t4r/ck6)
       ochrs = pep * (xkcx + xkcy + xkcz)
    end if

    !..opacity in presence of hydrogen  
    if (xh .ge. 1.0e-5) then 
       if (t6 .lt. t6_switch1) then 
          orad   = ochrs  
       else if (t6 .le. t6_switch2) then 
          orad   = 2.0d0*(ochrs*(1.5d0 - t6) + oiben2*(t6 - 1.0d0)) 
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
    th      = min(511.0d0, temp * 8.617d-8)
    fact    = 1.0d0 + 2.75d-2*th - 4.88d-5*th*th   
    facetax = 1.0d100    
    if (eta .le. 500.0) facetax = exp(0.522d0*eta - 1.563d0)  
    faceta  = 1.0d0 + facetax  
    ocompt  = 6.65205d-25/(fact * faceta) * xne/den
    orad    = orad   + ocompt 
    
    !..cutoff radiative opacity when 4kt/hbar is less than the plasma
    !..frequency
    tcut = con2 * sqrt(xne) 
    if (temp .lt. tcut) then 
       if (tcut .gt. 200.0*temp) then 
          orad   = orad * 2.658d86  
       else 
          cutfac   = exp(tcut/temp - 1.0d0)  
          orad     = orad * cutfac
       end if
    end if

    !..fudge molecular opacity for low temps
    xkf    = t7peek * den * (temp * 1.0d-7)**4
    orad   = xkf * orad/(xkf + orad)


    !..conductivity section:
    !..drel is the dividing line between nondegenerate and degenerate regions,
    !..taken from clayton eq. 2-34. if the density is larger than drel, then
    !..use the degenerate expressions. if the density is smaller than
    !..drelim, use the non-degenerate formulas. in between drel and drelim,
    !..apply a smooth blending of the two.
    
    dlog10   = log10(den)    

    drel     =  2.4d-7 * zbar/abar * temp * sqrt(temp)
    if (temp .le. 1.0d5) drel = drel * 15.0d0
    drel10   =  log10(drel)
    drelim   =  drel10 + 1.0d0


    !..from iben apj 196 525 1975 for non-degenerate regimes
    if (dlog10 .lt. drelim) then 
       zdel    = xne/(avo*t6*sqrt(t6))
       zdell10 = log10(zdel)    
       eta0    = exp(-1.20322d0 + twoth * log(zdel))  
       eta02   = eta0*eta0
       
       !..thpl factor 
       if (zdell10 .lt. 0.645) then 
          thpl    = -7.5668d0 + log(zdel * (1.0d0 + 0.024417d0*zdel))  
       else 
          if (zdell10 .lt. 2.5) then 
             thpl   = -7.58110d0 + log(zdel*(1.0d0 + 0.02804d0*zdel))  
             if (zdell10 .ge. 2.0) then 
                thpla = thpl 
                thpl  = -11.0742d0 + log(zdel**2 * (1.0d0 + 9.376d0/eta02))
                thpl  = 2.0d0*((2.5d0-zdell10)*thpla + (zdell10-2.0d0)*thpl)   
             end if
          else 
             thpl   = -11.0742d0 + log(zdel**2 * (1.0d0 + 9.376d0/eta02))
          end if
       end if

       !..pefac and walf factors
       if (zdell10 .lt. 2.0) then 
          pefac   = 1.0d0 + 0.021876d0*zdel 
          if (zdell10 .gt. 1.5) then 
             pefacal   = log(pefac) 
             pefacl    = log(0.4d0 * eta0 + 1.64496d0/eta0)  
             cfac1     = 2.0d0 - zdell10  
             cfac2     = zdell10 - 1.5d0  
             pefac     = exp(2.0d0 * (cfac1*pefacal + cfac2*pefacl))    
          end if
       else 
          pefac   = 0.4d0 * eta0 + 1.64496d0/eta0
       end if
     
       if (zdel.lt.40.0) then
          dnefac = 1.0d0 + zdel * (3.4838d-4 * zdel - 2.8966d-2)
       else
          dnefac = 1.5d0/eta0 * (1.0d0 - 0.8225d0/eta02)
       endif
       wpar2  = 9.24735d-3 * zdel * &
            (den*avo*(w(4)+w(5)+w(6))/xne + dnefac)/(sqrt(t6)*pefac)
       walf   = 0.5d0 * log(wpar2)    
       walf10 = 0.5d0 * log10(wpar2)    
       
       !..thx, thy and thc factors
       if (walf10 .le. -3.0) then 
          thx   = exp(2.413d0 - 0.124d0*walf) 
       else if (walf10 .le. -1.0) then 
          thx   = exp(0.299d0 - walf*(0.745d0 + 0.0456d0*walf))   
       else 
          thx   = exp(0.426d0 - 0.558d0*walf) 
       end if
       
       if (walf10 .le. -3.0) then 
          thy   = exp(2.158d0 - 0.111d0*walf)  
       else if (walf10 .le. 0.0) then 
          thy   = exp(0.553d0 - walf*(0.55d0 + 0.0299d0*walf))  
       else 
          thy   = exp(0.553d0 - 0.6d0*walf)   
       end if

       if (walf10 .le. -2.5) then 
          thc   = exp(2.924d0 - 0.1d0*walf)   
       else if (walf10 .le. 0.5) then 
          thc   = exp(1.6740d0 - walf*(0.511d0 + 0.0338d0*walf))  
       else 
          thc   = exp(1.941d0 - 0.785d0*walf) 
       end if
       
       oh   = (xh*thx + xhe*thy + w(6)*third*thc) / (t6*exp(thpl)) 
    end if

    !..from yakovlev & urpin soviet astro 1980 24 303 and 
    !..potekhin et al. 1997 aa 323 415 for degenerate regimes
    if (dlog10 .gt. drel10) then 
       xmas   = meff * xne**third 
       ymas   = sqrt(1.0d0 + xmas*xmas) 
       wfac   = weid * temp/ymas * xne 
       cint   = 1.0d0 
       
       !..ion-electron collision frequency and the thermal conductivity 
       vie   = iec * zbar * ymas * cint
       cie   = wfac/vie 
       
       !..electron-electron collision frequency and thermal conductivity 
       tpe  = xec * sqrt(xne/ymas)
       yg   = rt3 * tpe/temp 
       xrel = 1.009d0 * (zbar/abar * den * 1.0d-6)**third
       beta2 = xrel**2/(1.0d0 + xrel**2)
       jy   = (1.0d0 + 6.0d0/(5.0d0*xrel**2) + 2.0d0/(5.0d0*xrel**4)) &
            * ( yg**3 / (3.0d0 * (1.0d0 + 0.07414 * yg)**3) &
            * log((2.81d0 - 0.810*beta2 + yg)/yg)  &  
            + pi**5/6.0d0 * (yg/(13.91d0 + yg))**4 )
       vee = 0.511d0 * temp**2 * xmas/ymas**2 * sqrt(xmas/ymas) * jy 
       cee = wfac/vee

       !..total electron thermal conductivity and conversion to an opacity 
       ov1   = cie * cee/(cee + cie)
       ov    = k2c/(ov1*den) * temp**3 
    end if

    !..blend the opacities in the intermediate region
    if (dlog10 .le. drel10) then 
       ocond   = oh 
    else if (dlog10 .gt. drel10  .and. dlog10 .lt. drelim) then 
       x        = den
       x1       = 10.0d0**drel10
       x2       = 10.0d0**drelim
       alfa     = (x-x2)/(x1-x2)
       beta     = (x-x1)/(x2-x1)
       ocond    = alfa*oh + beta*ov
    else if (dlog10 .ge. drelim) then 
       ocond   = ov 
    end if

    !..total opacity
    opac    = orad * ocond / (ocond + orad)
    
    conductivity = k2c * temp**3 / (opac * den)

    return 
  end subroutine sig99

end module actual_conductivity_module
