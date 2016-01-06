      subroutine screen5(temp,den,zbar,abar,z2bar, &
                         z1,a1,z2,a2,jscreen,init, &
                         scor,scordt,scordd)

        use bl_constants_module, only: M_PI
        use actual_network, only: nrates
        implicit none

! this subroutine calculates screening factors and their derivatives
! for nuclear reaction rates in the weak, intermediate and strong regimes.
! based on graboske, dewit, grossman and cooper apj 181 457 1973 for
! weak screening. based on alastuey and jancovici apj 226 1034 1978,
! with plasma parameters from itoh et al apj 234 1079 1979, for strong
! screening.

! input:
! temp    = temperature
! den     = density
! zbar    = mean charge per nucleus
! abar    = mean number of nucleons per nucleus
! z2bar   = mean square charge per nucleus
! z1 a1   = charge and number in the entrance channel
! z2 a2   = charge and number in the exit channel
! jscreen = counter of which reaction is being calculated
! init    = flag to compute the more expensive functions just once

! output:
! scor    = screening correction
! scordt  = derivative of screening correction with temperature
! scordd  = derivative of screening correction with density


! declare the pass        
      integer          jscreen,init
      double precision temp,den,zbar,abar,z2bar,z1,a1,z2,a2, &
                       scor,scordt,scordd


! local variables
      double precision aa,daadt,daadd,bb,cc,dccdt,dccdd, &
                       pp,dppdt,dppdd,qq,dqqdt,dqqdd,rr,drrdt,drrdd, &
                       ss,dssdt,dssdd,tt,dttdt,dttdd,uu,duudt,duudd, &
                       vv,dvvdt,dvvdd,a3,da3,tempi,dtempi,deni, &
                       qlam0z,qlam0zdt,qlam0zdd, &
                       h12w,dh12wdt,dh12wdd,h12,dh12dt,dh12dd, &
                       h12x,dh12xdt,dh12xdd,alfa,beta, &
                       taufac,taufacdt,gamp,gampdt,gampdd, &
                       gamef,gamefdt,gamefdd, &
                       tau12,tau12dt,alph12,alph12dt,alph12dd, &
                       xlgfac,dxlgfacdt,dxlgfacdd, &
                       gamp14,gamp14dt,gamp14dd, &
                       xni,dxnidd,ytot, &
                       temp_old,den_old,zbar_old,abar_old


! screening variables
! zs13    = (z1+z2)**(1./3.)
! zhat    = combination of z1 and z2 raised to the 5/3 power
! zhat2   = combination of z1 and z2 raised to the 5/12 power
! lzav    = log of effective charge
! aznut   = combination of a1,z1,a2,z2 raised to 1/3 power


      integer          nscreen_max
      parameter        (nscreen_max = 2*nrates + 40)

      double precision zs13(nscreen_max),zhat(nscreen_max), &
                       zhat2(nscreen_max),lzav(nscreen_max), &
                       aznut(nscreen_max),zs13inv(nscreen_max)


! parameter fact is the cube root of 2
      double precision  x13,x14,x53,x532,x512,fact,co2,gamefx,gamefs, &
                        blend_frac
      parameter        (x13    = 1.0d0/3.0d0, &
                        x14    = 1.0d0/4.0d0, &
                        x53    = 5.0d0/3.0d0, &
                        x532   = 5.0d0/32.0d0, &
                        x512   = 5.0d0/12.0d0, &
                        fact   = 1.25992104989487d0, &
                        co2    = x13 * 4.248719d3, &
                        gamefx = 0.3d0, &
                        gamefs = 0.8d0, &
                        blend_frac = 0.05d0)


!      data     temp_old/-1.0d0/, den_old/-1.0d0/, &
!               zbar_old/-1.0d0/, abar_old/-1.0d0/




! compute and store the more expensive screening factors

! MZ: caching this stuff is not threadsafe -- we should precompute this at
! initialization and store it in a module variable

!      if (init .eq. 1) then
       if (jscreen .gt. nscreen_max) &
       stop 'jscreen > nscreen_max in screen5'
       zs13(jscreen)    = (z1 + z2)**x13
       zs13inv(jscreen) = 1.0d0/zs13(jscreen)
       zhat(jscreen)    = (z1 + z2)**x53  - z1**x53 - z2**x53
       zhat2(jscreen)   = (z1 + z2)**x512 - z1**x512 -z2**x512
       lzav(jscreen)    = x53 * log(z1*z2/(z1 + z2))
       aznut(jscreen)   = (z1**2 * z2**2 * a1*a2 / (a1 + a2))**x13
!      endif


! MZ: caching this stuff is not threadsafe.  We should create a
! derived type to hold the plasma parameters and fill it once in the
! routine that calls the screening for each rate (since they all have
! the same T and rho)

! calculate average plasma, if need be
!      if (temp_old .ne. temp .or. &
!          den_old  .ne. den  .or. &
!          zbar_old  .ne. zbar  .or. &
!          abar_old  .ne. abar ) then

!       temp_old = temp
!       den_old  = den
!       zbar_old  = zbar
!       abar_old  = abar

       ytot     = 1.0d0/abar
       rr       = den * ytot
       tempi   = 1.0d0/temp
       dtempi  = -tempi*tempi
       deni    = 1.0d0/den

       pp       = sqrt(rr*tempi*(z2bar + zbar))
       qq       = 0.5d0/pp *(z2bar + zbar)
       dppdt    = qq*rr*dtempi
       !dppdd    = qq*ytot*tempi

       qlam0z   = 1.88d8 * tempi * pp
       qlam0zdt = 1.88d8 * (dtempi*pp + tempi*dppdt)
       !qlam0zdd = 1.88d8 * tempi * dppdd

       taufac   = co2 * tempi**x13
       taufacdt = -x13*taufac*tempi

       qq      = rr*zbar
       xni     = qq**x13
       !dxnidd  = x13 * xni * deni

       aa     = 2.27493d5 * tempi * xni
       daadt  = 2.27493d5 * dtempi * xni
       !daadd  = 2.27493d5 * tempi * dxnidd
!      end if


! calculate individual screening factors
      bb       = z1 * z2
      gamp     = aa
      gampdt   = daadt
      !gampdd   = daadd

      qq       = fact * bb * zs13inv(jscreen)
      gamef    = qq * gamp
      gamefdt  = qq * gampdt
      !gamefdd  = qq * gampdd

      tau12    = taufac * aznut(jscreen)
      tau12dt  = taufacdt * aznut(jscreen)

      qq       = 1.0d0/tau12
      alph12   = gamef * qq
      alph12dt = (gamefdt - alph12*tau12dt) * qq
      !alph12dd = gamefdd * qq



! limit alph12 to 1.6 to prevent unphysical behavior.
! this should really be replaced by a pycnonuclear reaction rate formula
      if (alph12 .gt. 1.6) then
       alph12   = 1.6d0
       alph12dt = 0.0d0
       !alph12dd = 0.0d0

       gamef    = 1.6d0 * tau12
       gamefdt  = 1.6d0 * tau12dt
       !gamefdd  = 0.0d0

       qq       = zs13(jscreen)/(fact * bb)
       gamp     = gamef * qq
       gampdt   = gamefdt * qq
       !gampdd   = 0.0d0
      end if



! weak screening regime
      h12w    = bb * qlam0z
      dh12wdt = bb * qlam0zdt
      !dh12wdd = bb * qlam0zdd

      h12     = h12w
      dh12dt  = dh12wdt
      !dh12dd  = dh12wdd



! intermediate and strong sceening regime
      if (gamef .gt. gamefx) then

       gamp14   = gamp**x14
       rr       = 1.0d0/gamp
       qq       = 0.25d0*gamp14*rr
       gamp14dt = qq * gampdt
       !gamp14dd = qq * gampdd

       cc       =   0.896434d0 * gamp * zhat(jscreen) &
                  - 3.44740d0  * gamp14 * zhat2(jscreen) &
                  - 0.5551d0   * (log(gamp) + lzav(jscreen)) &
                  - 2.996d0

       dccdt    =   0.896434d0 * gampdt * zhat(jscreen) &
                  - 3.44740d0  * gamp14dt * zhat2(jscreen) &
                  - 0.5551d0*rr*gampdt

       dccdd    =   0.896434d0 * gampdd * zhat(jscreen) &
                  - 3.44740d0  * gamp14dd * zhat2(jscreen) &
                  - 0.5551d0*rr*gampdd

       a3     = alph12 * alph12 * alph12
       da3    = 3.0d0 * alph12 * alph12

       qq     = 0.014d0 + 0.0128d0*alph12
       dqqdt  = 0.0128d0*alph12dt
       !dqqdd  = 0.0128d0*alph12dd

       rr     = x532 - alph12*qq
       drrdt  = -(alph12dt*qq + alph12*dqqdt)
       !drrdd  = -(alph12dd*qq + alph12*dqqdd)

       ss     = tau12*rr
       dssdt  = tau12dt*rr + tau12*drrdt
       !dssdd  = tau12*drrdd

       tt     =  -0.0098d0 + 0.0048d0*alph12
       dttdt  = 0.0048d0*alph12dt
       !dttdd  = 0.0048d0*alph12dd

       uu     =  0.0055d0 + alph12*tt
       duudt  = alph12dt*tt + alph12*dttdt
       !duudd  = alph12dd*tt + alph12*dttdd

       vv   = gamef * alph12 * uu
       dvvdt= gamefdt*alph12*uu + gamef*alph12dt*uu + gamef*alph12*duudt
       !dvvdd= gamefdd*alph12*uu + gamef*alph12dd*uu + gamef*alph12*duudd

       h12     = cc - a3 * (ss + vv)
       rr      = da3 * (ss + vv)
       dh12dt  = dccdt - rr*alph12dt - a3*(dssdt + dvvdt)
       !dh12dd  = dccdd - rr*alph12dd - a3*(dssdd + dvvdd)

       rr     =  1.0d0 - 0.0562d0*a3
       ss     =  -0.0562d0*da3
       drrdt  = ss*alph12dt
       !drrdd  = ss*alph12dd

       if (rr .ge. 0.77d0) then
        xlgfac    = rr
        dxlgfacdt = drrdt
        !dxlgfacdd = drrdd
       else
        xlgfac    = 0.77d0
        dxlgfacdt = 0.0d0
        !dxlgfacdd = 0.0d0
       end if


       h12    = log(xlgfac) + h12
       rr     = 1.0d0/xlgfac
       dh12dt = rr*dxlgfacdt + dh12dt
       !dh12dd = rr*dxlgfacdd + dh12dd


       if (gamef .le. gamefs) then
        rr     =  2.0d0*(gamefs - gamef)
        drrdt  = -2.0d0*gamefdt
        !drrdd  = -2.0d0*gamefdd

        ss     = 2.0d0*(gamef - gamefx)
        dssdt  = 2.0d0*gamefdt
        !dssdd  = 2.0d0*gamefdd


! store current values for possible blending
        h12x    = h12
        dh12xdt = dh12dt
        !dh12xdd = dh12dd

        vv     = h12
        h12    = h12w*rr + vv*ss
        dh12dt = dh12wdt*rr + h12w*drrdt + dh12dt*ss + vv*dssdt
        !dh12dd = dh12wdd*rr + h12w*drrdd + dh12dd*ss + vv*dssdd

! blend the transition region - from bill paxton
       if (gamefs - gamef .lt. blend_frac*(gamefs - gamefx)) then
         alfa   = (gamefs - gamef) / (blend_frac*(gamefs - gamefx))
         alfa   = 0.5d0 * (1d0 - cos(M_PI*alfa))
         beta   = 1.0d0 - alfa
         h12    = alfa * h12 + beta * h12x
         dh12dt = alfa * dh12dt + beta * dh12xdt
         !dh12dd = alfa * dh12dd + beta * dh12xdd
        end if
       end if


! end of intermediate and strong screening if
      end if


! machine limit the output
      h12    = max(min(h12,300.0d0),0.0d0)
      scor   = exp(h12)
      if (h12 .eq. 300.0d0) then
         scordt = 0.0d0
       !scordd = 0.0d0
      else
         scordt = scor * dh12dt
       !scordd = scor * dh12dd
      end if

!      write(6,111) 'weak =',h12w,' total =',h12,
!     1             ' 1-ratio =',1.0d0-h12w/h12,' correction',scor
! 111  format(1x,4(a,1pe13.6))
!      read(5,*)

      return
      end









