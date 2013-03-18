      program drive_aprox13
      include 'implno.dek'
      include 'const.dek'
      include 'timers.dek'
      include 'burn_common.dek'
      include 'network.dek'

! this program exercises the aprox13 network

! declare
      character*40     summary
      integer          i,j,nok,nbad

      double precision tstart,tstep,conserv,tin,din,ein,vin,zin,xin(18), &
                       tout,dout,eout,xout(18),edum



! initialize the network and eos
      call init_aprox13
      call read_helm_table


! keep coming back to here
20    continue

      call net_input(tstart,tstep,tin,din,vin,zin,ein,xin)


! start the clock
      call zsecond(timzer)


! a message
        write(6,*)
        write(6,*) 'starting integration'
        write(6,*)

! burn it

        call burner(tstart,tstep, &
                    tin,din,vin,zin,ein,xin, &
                    tout,dout,eout,xout, &
                    conserv,nok,nbad)


!      edum = ev2erg*1.0d6*avo * sum((xout(1:ionmax) - xin(1:ionmax))/aion(1:ionmax)*bion(1:ionmax))
!      write(6,123) edum
!      write(6,123) edum - sneut*tstep
!123   format(1x,1pe14.6)



! output a summary of the integration

       call net_summary(tstep,tin,din,ein, &
                        tout,dout,eout,conserv, &
                        nbad,nok,xout)


! back for another input point
      goto 20
      end








      subroutine burner(beg,tstep, &
                        tin,din,vin,zin,ein,xin, &
                        tout,dout,eout,xout, &
                        conserv,nok,nbad)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'



! input:
! beg       starting time
! tstep     time over which to integrate
! tin       initial temperature
! din       initial density
! vin       initial velocity                                                                                           
! zin       initial position
! ein       initial internal energy
! xin(1:19) initial composition vector

! output:
! tout       final temperature
! dout       final density
! eout       final internal energy
! xout(1:19) final composition vector
! conserv    1 - sum of mass fractions
! nok        number of good time steps taken
! nbad       number of bad timesteps attempted


! declare the pass
      integer          nok,nbad
      double precision beg,tstep,tin,din,vin,zin,ein,xin(*), &
                       tout,dout,eout,xout(*),conserv


! local variables
      integer          i
      double precision abar,zbar,wbar,ye,xcess

! for the integration driver
      integer          kount
      double precision stptry,stpmin,tend,ys2(abignet*nzmax), &
                       odescal,tol

! usually adequate
      parameter        (tol     = 1.0d-6, &
                        odescal = 1.0d-8)

! for very accurate integrations
!      parameter        (tol     = 1.0d-10, &
!                        odescal = 1.0d-12)


      external         aprox13,saprox13,baprox13,daprox13


!      external         forder_ma28
!      external         forder_umf
!      external         forder_y12m
!      external         forder_ludcmp
!      external         forder_leqs
!      external         forder_lapack
!      external         forder_gift
!      external         forder_biconj

!      external         rosen_ma28
!      external         rosen_umf
!      external         rosen_y12m
!      external         rosen_ludcmp
!      external         rosen_leqs
!      external         rosen_lapack
!      external         rosen_gift
!      external         rosen_biconj

      external         stifbs_ma28
!      external         stifbs_umf
!      external         stifbs_y12m
!      external         stifbs_ludcmp
!      external         stifbs_leqs
!      external         stifbs_lapack
!      external         stifbs_gift
!      external         stifbs_biconj





! set the initial condition


! load the mass fractions
       xmass(ionbeg:ionend) = xin(ionbeg:ionend)


! get abar, zbar and a few other composition variables
       call azbar(xmass(ionbeg),aion(ionbeg),zion(ionbeg),wion(ionbeg),ionmax, &
                  ymass(ionbeg),abar,zbar,wbar,ye,xcess)


! stuff the initial conditions into ys2
       ys2(ionbeg:ionend) = ymass(ionbeg:ionend)
       ys2(iener) = ein
       ys2(itemp) = tin
       ys2(iden)  = din
       ys2(ivelx) = vin
       ys2(iposx) = zin


! single step (tend=tstep), hydrostatic, or expansion ending times.
! the variable tstep has two meanings here. tstep in single step mode
! is the size of the time step to try. tstep in hydrostatic or expansion
! mode is the ending integration time. the integration driver really
! gets some exercise if tstep is large in single step mode.

      tend = tstep
      if (one_step) then
       stptry = tstep
       stpmin = tstep * 1.0d-20
      else
       stptry = max(beg * 1.0d-10,1.0d-16)
       stpmin = stptry * 1.0d-12
      end if



! integrate the aprox13 network
       call netint(beg,stptry,stpmin,tend,ys2, &
                   tol,neqs,nok,nbad,kount,odescal, &
!                  aprox13,saprox13,baprox13,forder_ma28)
!                  aprox13,saprox13,baprox13,forder_umf)
!                  aprox13,saprox13,baprox13,forder_y12m)
!                  aprox13,daprox13,baprox13,forder_ludcmp)
!                  aprox13,daprox13,baprox13,forder_leqs)
!                  aprox13,daprox13,baprox13,forder_lapack)
!                  aprox13,daprox13,baprox13,forder_gift)
!                  aprox13,saprox13,baprox13,forder_biconj)
!                  aprox13,saprox13,baprox13,rosen_ma28)
!                  aprox13,saprox13,baprox13,rosen_umf)
!                  aprox13,saprox13,baprox13,rosen_y12m)
!                  aprox13,daprox13,baprox13,rosen_ludcmp)
!                  aprox13,daprox13,baprox13,rosen_leqs)
!                  aprox13,daprox13,baprox13,rosen_lapack)
!                  aprox13,daprox13,baprox13,rosen_gift)
!                  aprox13,saprox13,baprox13,rosen_biconj)
                  aprox13,saprox13,baprox13,stifbs_ma28)
!                  aprox13,saprox13,baprox13,stifbs_umf)
!                  aprox13,saprox13,baprox13,stifbs_y12m)
!                  aprox13,daprox13,baprox13,stifbs_ludcmp)
!                  aprox13,daprox13,baprox13,stifbs_leqs)
!                  aprox13,daprox13,baprox13,stifbs_lapack)
!                  aprox13,daprox13,baprox13,stifbs_gift)
!                  aprox13,saprox13,baprox13,stifbs_biconj)




! set the output
      do i=ionbeg,ionend
       xout(i) = ys2(i) * aion(i)
      enddo
      tout = ys2(itemp)
      dout = ys2(iden)
      eout = ys2(iener)

      conserv = 0.0d0
      do i=ionbeg,ionend
       conserv = conserv + xout(i)
      enddo
      conserv = 1.0d0 - conserv

      return
      end



!---------------------------------------------------------------------







!---------------------------------------------------------------------
! this file contains aprox13 network

! routine aprox13 sets up the odes
! routine rhs evaluates the right hand sides
! routine daprox13 sets up the dense aprox13 jacobian
! routine baprox13 builds the nonzero locations for saprox13
! routine saprox13 sets up the sparse aprox13 jacobian
! routine aprox13rat generates the reaction rates for routine aprox13
! routine aprox13tab generates the raw rates using table interpolation
! routine screen_aprox13 applies screening corrections to the raw rates
! routine init_aprox13 initializes the aprox13 network




      subroutine aprox13(tt,y,dydt)
      include 'implno.dek'
      include 'const.dek'
      include 'burn_common.dek'
      include 'network.dek'
      include 'vector_eos.dek'
      include 'cjdet.dek'

! this routine sets up the system of ode's for the aprox13 nuclear reactions.
! this is an alpha chain + heavy ion network with (a,p)(p,g) links
!
! isotopes: he4,  c12,  o16,  ne20, mg24, si28, s32,
!           ar36, ca40, ti44, cr48, fe52, ni56



! declare the pass
      double precision tt,y(1),dydt(1)


! local variables
      logical          deriva
      parameter        (deriva = .false.)
      integer          i
      double precision enuc,taud,taut,z,denom,suma,sumz,ww, &
                       zbarxx,z2barxx,ytot1,abar,zbar,ye,z2bar, &
                       velx,posx,cs,dpde, &
                       snuda,snudz,combo,phi,dtdp



! positive definite mass fractions
      do i=ionbeg,ionend
       y(i) = min(1.0d0,max(y(i),1.0d-30))
      enddo


! generate abar and zbar for this composition
      abar = 1.0d0/sum(y(ionbeg:ionend))
      zbar = sum(zion(ionbeg:ionend)*y(ionbeg:ionend)) * abar
      ye    = zbar/abar


! positive definite temperatures and densities
       y(itemp) = min(1.0d11,max(y(itemp),1.0d4))
       y(iden)  = min(1.0d11,max(y(iden),1.0d-10))


! for evolution equations with the network
      if (pure_network .eq. 0) then

! get the new temperature and density if need be
       if (trho_hist) call update2(tt,y(itemp),y(iden))

       if (self_heat_const_pres) then
        jlo_eos = 1
        jhi_eos = 1
        ptot_row(1) = bpres
        temp_row(1) = y(itemp)
        abar_row(1) = abar
        zbar_row(1) = zbar

        den_row(1)  = y(iden)
        call invert_helm_pt_quiet
        y(iden) = den_row(1)
       end if


! positive definite temperatures and densities
       y(itemp) = min(1.0d11,max(y(itemp),1.0d4))
       y(iden)  = min(1.0d11,max(y(iden),1.0d-10))


! set the common block temperature and density
       btemp = y(itemp)
       bden  = y(iden)


! for pure network
      else if (pure_network .eq. 1) then
       if (trho_hist) call update2(tt,btemp,bden)
      end if


! get the reaction rates
       if (use_tables .eq. 1) then
        call aprox13tab
       else
        call aprox13rat
       end if


! do the screening here because the corrections depend on the composition

       call screen_aprox13(y)



! get the right hand side of the odes

       call rhs(y,ratdum,dydt,deriva)



! if we are doing a pure network, we are done

      if (pure_network .eq. 1) return


! instantaneous energy generation rate

       call ener_gener_rate(dydt,enuc)
       sdot = enuc


! get the neutrino losses
       call sneut5(btemp,bden,abar,zbar, &
                   sneut,dsneutdt,dsneutdd,snuda,snudz)


! append an energy equation
       dydt(iener) = enuc - sneut



! the type of temperature and density ode's depend
! on the burn mode:


! hydrostatic or single step cases
       if (hydrostatic  .or.  one_step  .or.  trho_hist) then
        dydt(itemp) = 0.0d0
        dydt(iden)  = 0.0d0
        dydt(ivelx) = 0.0d0
        dydt(iposx) = 0.0d0



! adiabatic expansion or contraction
       else if (expansion) then

        taud = 446.0d0/sqrt(den0)
        taut = 3.0d0 * taud

        dydt(itemp) = -psi * y(itemp)/taut
        dydt(iden)  = -psi * y(iden)/taud
        dydt(ivelx) = 0.0d0
        dydt(iposx) = 0.0d0



! self heating
       else if (self_heat_const_den) then


! call an eos
       temp_row(1) = y(itemp)
       den_row(1)  = y(iden)
       abar_row(1) = abar
       zbar_row(1) = zbar
       jlo_eos = 1
       jhi_eos = 1

       call helmeos


! temperature equation that is self-consistent with an eos
       z           = 1.0d0/cv_row(1)
       dydt(itemp) = z*dydt(iener)


! density, velocity, and position equations
       dydt(iden)  = 0.0d0
       dydt(ivelx) = 0.0d0
       dydt(iposx) = 0.0d0



! detonation
       else if (detonation) then


! map the rest of the input vector
        velx = y(ivelx)
        posx = y(iposx)


! call an eos
        temp_row(1) = btemp
        den_row(1)  = bden
        abar_row(1) = abar
        zbar_row(1) = zbar
        jlo_eos = 1
        jhi_eos = 1

        call helmeos


! for de/dy and dp/dy
        suma = 0.0d0
        do i=ionbeg,ionend
         suma = suma - dydt(i)
        enddo

        sumz = 0.0d0
        do i=ionbeg,ionend
         sumz = sumz + (zion(i) - zbar)*dydt(i)
        enddo


! the possibly singular denominator
        cs      = cs_row(1)
        denom   = velx*velx - cs*cs

! the function phi
        dpde = dpt_row(1)/det_row(1)
        z    = suma*dpa_row(1)*abar*abar + sumz*dpz_row(1)*abar
        ww   = suma*dea_row(1)*abar*abar + sumz*dez_row(1)*abar
        phi  = dpde*(dydt(iener) - ww) - z


! a common combination
        if (denom .ne. 0.0) then
         combo = phi/denom
        else
         combo = 0.0d0
        end if


! position equation
        dydt(iposx) = velx


! density equation
        dydt(iden) = combo


! velocity equations
        dydt(ivelx) = -velx/bden*dydt(iden)


! temperature equation
        dtdp          = 1.0d0/dpt_row(1)
        ww            = suma*dpa_row(1)*abar*abar + sumz*dpz_row(1)*abar
        dydt(itemp) = dtdp*((velx*velx &
                         - dpd_row(1))*dydt(iden) - ww)

       end if


      return
      end





      subroutine rhs(y,rate,dydt,deriva)
      include 'implno.dek'
      include 'network.dek'

! evaluates the right hand side of the aprox13 odes

! declare the pass
      logical          deriva
      double precision y(1),rate(1),dydt(1)


! local variables
      integer          i
      double precision sixth
      parameter        (sixth = 1.0d0/6.0d0)


! zero the abundance odes
      dydt(ionbeg:ionend) = 0.0d0


! set up the system of odes:
! he4 reactions
! heavy ion reactions
      dydt(ihe4) = &
             0.5d0 * y(ic12) * y(ic12) * rate(ir1212) &
           + 0.5d0 * y(ic12) * y(io16) * rate(ir1216) &
           + 0.56d0* 0.5d0 * y(io16)*y(io16) * rate(ir1616)

! (a,g) and (g,a) reactions
      dydt(ihe4) =  dydt(ihe4) &
           - 0.5d0 * y(ihe4)*y(ihe4)*y(ihe4)*rate(ir3a) &
           + 3.0d0 * y(ic12) * rate(irg3a) &
           - y(ihe4)  * y(ic12) * rate(ircag) &
           + y(io16)  * rate(iroga) &
           - y(ihe4)  * y(io16) * rate(iroag) &
           + y(ine20) * rate(irnega) &
           - y(ihe4)  * y(ine20) * rate(irneag) &
           + y(img24) * rate(irmgga) &
           - y(ihe4)  * y(img24)* rate(irmgag) &
           + y(isi28) * rate(irsiga) &
           - y(ihe4)  * y(isi28)*rate(irsiag) &
           + y(is32)  * rate(irsga)

      dydt(ihe4) =  dydt(ihe4) &
           - y(ihe4)  * y(is32) * rate(irsag) &
           + y(iar36) * rate(irarga) &
           - y(ihe4)  * y(iar36)*rate(irarag) &
           + y(ica40) * rate(ircaga) &
           - y(ihe4)  * y(ica40)*rate(ircaag) &
           + y(iti44) * rate(irtiga) &
           - y(ihe4)  * y(iti44)*rate(irtiag) &
           + y(icr48) * rate(ircrga) &
           - y(ihe4)  * y(icr48)*rate(ircrag) &
           + y(ife52) * rate(irfega) &
           - y(ihe4)  * y(ife52) * rate(irfeag) &
           + y(ini56) * rate(irniga)


! (a,p)(p,g) and (g,p)(p,a) reactions

      if (.not.deriva) then
       dydt(ihe4) =  dydt(ihe4) &
           + 0.34d0*0.5d0*y(io16)*y(io16)*rate(irs1)*rate(ir1616) &
           - y(ihe4)  * y(img24) * rate(irmgap)*(1.0d0-rate(irr1)) &
           + y(isi28) * rate(irsigp) * rate(irr1) &
           - y(ihe4)  * y(isi28) * rate(irsiap)*(1.0d0-rate(irs1)) &
           + y(is32)  * rate(irsgp) * rate(irs1) &
           - y(ihe4)  * y(is32) * rate(irsap)*(1.0d0-rate(irt1)) &
           + y(iar36) * rate(irargp) * rate(irt1) &
           - y(ihe4)  * y(iar36) * rate(irarap)*(1.0d0-rate(iru1)) &
           + y(ica40) * rate(ircagp) * rate(iru1) &
           - y(ihe4)  * y(ica40) * rate(ircaap)*(1.0d0-rate(irv1)) &
           + y(iti44) * rate(irtigp) * rate(irv1)

       dydt(ihe4) =  dydt(ihe4) &
           - y(ihe4)  * y(iti44) * rate(irtiap)*(1.0d0-rate(irw1)) &
           + y(icr48) * rate(ircrgp) * rate(irw1) &
           - y(ihe4)  * y(icr48) * rate(ircrap)*(1.0d0-rate(irx1)) &
           + y(ife52) * rate(irfegp) * rate(irx1) &
           - y(ihe4)  * y(ife52) * rate(irfeap)*(1.0d0-rate(iry1)) &
           + y(ini56) * rate(irnigp) * rate(iry1)

      else
       dydt(ihe4) =  dydt(ihe4) &
           + 0.34d0*0.5d0*y(io16)*y(io16)* &
             (ratdum(irs1)*rate(ir1616) + rate(irs1)*ratdum(ir1616)) &
           - y(ihe4)*y(img24)*(rate(irmgap)*(1.0d0 - ratdum(irr1)) &
                                     - ratdum(irmgap)*rate(irr1)) &
           + y(isi28) * (ratdum(irsigp) * rate(irr1) + &
                            rate(irsigp) * ratdum(irr1)) &
           - y(ihe4)*y(isi28)*(rate(irsiap)*(1.0d0 - ratdum(irs1)) &
                                     - ratdum(irsiap)*rate(irs1)) &
           + y(is32)  * (ratdum(irsgp) * rate(irs1) + &
                            rate(irsgp) * ratdum(irs1))

       dydt(ihe4) =  dydt(ihe4) &
           - y(ihe4)*y(is32)*(rate(irsap)*(1.0d0 - ratdum(irt1)) &
                                    - ratdum(irsap)*rate(irt1)) &
           + y(iar36) * (ratdum(irargp) * rate(irt1) + &
                            rate(irargp) * ratdum(irt1)) &
           - y(ihe4)*y(iar36)*(rate(irarap)*(1.0d0 - ratdum(iru1)) &
                                     - ratdum(irarap)*rate(iru1)) &
           + y(ica40) * (ratdum(ircagp) * rate(iru1) + &
                            rate(ircagp) * ratdum(iru1)) &
           - y(ihe4)*y(ica40)*(rate(ircaap)*(1.0d0-ratdum (irv1)) &
                                     - ratdum(ircaap)*rate(irv1)) &
           + y(iti44) * (ratdum(irtigp) * rate(irv1) + &
                            rate(irtigp) * ratdum(irv1))

       dydt(ihe4) =  dydt(ihe4) &
           - y(ihe4)*y(iti44)*(rate(irtiap)*(1.0d0 - ratdum(irw1)) &
                                     - ratdum(irtiap)*rate(irw1)) &
           + y(icr48) * (ratdum(ircrgp) * rate(irw1) + &
                            rate(ircrgp) * ratdum(irw1)) &
           - y(ihe4)*y(icr48)*(rate(ircrap)*(1.0d0 - ratdum(irx1)) &
                                     - ratdum(ircrap)*rate(irx1)) &
           + y(ife52) * (ratdum(irfegp) * rate(irx1) + &
                            rate(irfegp) * ratdum(irx1)) &
           - y(ihe4)*y(ife52)*(rate(irfeap)*(1.0d0 - ratdum(iry1)) &
                                    - ratdum(irfeap)*rate(iry1)) &
           + y(ini56) * (ratdum(irnigp) * rate(iry1) + &
                            rate(irnigp) * ratdum(iry1))
      end if




! c12 reactions
      dydt(ic12) = -y(ic12) * y(ic12) * rate(ir1212) &
           - y(ic12) * y(io16) * rate(ir1216) &
           + sixth * y(ihe4)*y(ihe4)*y(ihe4)*rate(ir3a) &
           - y(ic12) * rate(irg3a) &
           - y(ic12) * y(ihe4) * rate(ircag) &
           + y(io16) * rate(iroga)

! o16 reactions
      dydt(io16) = -y(ic12) * y(io16) * rate(ir1216) &
           - y(io16) * y(io16) * rate(ir1616) &
           + y(ic12) * y(ihe4) * rate(ircag) &
           - y(io16) * y(ihe4) * rate(iroag) &
           - y(io16) * rate(iroga) &
           + y(ine20) * rate(irnega)

! ne20 reactions
      dydt(ine20) =  0.5d0 * y(ic12) * y(ic12) * rate(ir1212) &
           + y(io16) * y(ihe4) * rate(iroag) &
           - y(ine20) * y(ihe4) * rate(irneag) &
           - y(ine20) * rate(irnega) &
           + y(img24) * rate(irmgga)


! mg24 reactions
      dydt(img24)  = 0.5d0 * y(ic12) * y(io16) * rate(ir1216) &
           + y(ine20) * y(ihe4) * rate(irneag) &
           - y(img24) * y(ihe4) * rate(irmgag) &
           - y(img24) * rate(irmgga) &
           + y(isi28) * rate(irsiga)

      if (.not.deriva) then
       dydt(img24)  = dydt(img24) &
           - y(img24) * y(ihe4) * rate(irmgap)*(1.0d0-rate(irr1)) &
           + y(isi28) * rate(irr1) * rate(irsigp)

      else
       dydt(img24)  = dydt(img24) &
           - y(img24)*y(ihe4)*(rate(irmgap)*(1.0d0 - ratdum(irr1)) &
                                     - ratdum(irmgap)*rate(irr1)) &
           + y(isi28) * (ratdum(irr1) * rate(irsigp) + &
                            rate(irr1) * ratdum(irsigp))
      end if



! si28 reactions
      dydt(isi28)  =  0.5d0 * y(ic12) * y(io16) * rate(ir1216) &
           + 0.56d0*0.5d0*y(io16)*y(io16)*rate(ir1616) &
           + y(img24) * y(ihe4) * rate(irmgag) &
           - y(isi28) * y(ihe4) * rate(irsiag) &
           - y(isi28) * rate(irsiga) &
           + y(is32)  * rate(irsga)

      if (.not.deriva) then
       dydt(isi28)  = dydt(isi28) &
           + 0.34d0*0.5d0*y(io16)*y(io16)*rate(irs1)*rate(ir1616) &
           + y(img24) * y(ihe4) * rate(irmgap)*(1.0d0-rate(irr1)) &
           - y(isi28) * rate(irr1) * rate(irsigp) &
           - y(isi28) * y(ihe4) * rate(irsiap)*(1.0d0-rate(irs1)) &
           + y(is32)  * rate(irs1) * rate(irsgp)

      else
       dydt(isi28)  = dydt(isi28) &
           + 0.34d0*0.5d0*y(io16)*y(io16)* &
              (ratdum(irs1)*rate(ir1616) + rate(irs1)*ratdum(ir1616)) &
           + y(img24)*y(ihe4)*(rate(irmgap)*(1.0d0 - ratdum(irr1)) &
                                     - ratdum(irmgap)*rate(irr1)) &
           - y(isi28)*(ratdum(irr1) * rate(irsigp) + &
                          rate(irr1) * ratdum(irsigp)) &
           - y(isi28)*y(ihe4)*(rate(irsiap)*(1.0d0 - ratdum(irs1)) &
                                     - ratdum(irsiap)*rate(irs1)) &
           + y(is32)*(ratdum(irs1) * rate(irsgp) + &
                          rate(irs1) * ratdum(irsgp))
      end if



! s32 reactions
      dydt(is32) = 0.1d0 * 0.5d0 *y(io16)*y(io16)*rate(ir1616) &
            +  y(isi28) * y(ihe4) * rate(irsiag) &
            - y(is32) * y(ihe4) * rate(irsag) &
            - y(is32) * rate(irsga) &
            + y(iar36) * rate(irarga)

      if (.not.deriva) then
       dydt(is32)  = dydt(is32) &
            + 0.34d0*0.5d0*y(io16)**2*rate(ir1616)*(1.0d0-rate(irs1)) &
            + y(isi28) * y(ihe4) * rate(irsiap)*(1.0d0-rate(irs1)) &
            - y(is32) * rate(irs1) * rate(irsgp) &
            - y(is32) * y(ihe4) * rate(irsap)*(1.0d0-rate(irt1)) &
            + y(iar36) * rate(irt1) * rate(irargp)
      else
       dydt(is32)  = dydt(is32) &
            + 0.34d0*0.5d0*y(io16)**2 * &
           (rate(ir1616)*(1.0d0-ratdum(irs1))-ratdum(ir1616)*rate(irs1)) &
            + y(isi28)*y(ihe4)*(rate(irsiap)*(1.0d0-ratdum(irs1)) &
                                     - ratdum(irsiap)*rate(irs1)) &
            - y(is32)*(ratdum(irs1) * rate(irsgp) + &
                          rate(irs1) * ratdum(irsgp)) &
            - y(is32)*y(ihe4)*(rate(irsap)*(1.0d0-ratdum(irt1)) &
                                     - ratdum(irsap)*rate(irt1)) &
            + y(iar36)*(ratdum(irt1) * rate(irargp) + &
                           rate(irt1) * ratdum(irargp))
      end if


! ar36 reactions
      dydt(iar36) =  y(is32)  * y(ihe4) * rate(irsag) &
            - y(iar36) * y(ihe4) * rate(irarag) &
            - y(iar36) * rate(irarga) &
            + y(ica40) * rate(ircaga)

      if (.not.deriva) then
       dydt(iar36)  = dydt(iar36) &
            + y(is32)  * y(ihe4) * rate(irsap)*(1.0d0-rate(irt1)) &
            - y(iar36) * rate(irt1) * rate(irargp) &
            - y(iar36) * y(ihe4) * rate(irarap)*(1.0d0-rate(iru1)) &
            + y(ica40) * rate(ircagp) * rate(iru1)

      else
       dydt(iar36)  = dydt(iar36) &
            + y(is32)*y(ihe4)*(rate(irsap)*(1.0d0 - ratdum(irt1)) &
                                     - ratdum(irsap)*rate(irt1)) &
            - y(iar36)*(ratdum(irt1) * rate(irargp) + &
                           rate(irt1) * ratdum(irargp)) &
            - y(iar36)*y(ihe4)*(rate(irarap)*(1.0d0-ratdum(iru1)) &
                                      - ratdum(irarap)*rate(iru1)) &
            + y(ica40)*(ratdum(ircagp) * rate(iru1) + &
                           rate(ircagp) * ratdum(iru1))
      end if


! ca40 reactions
      dydt(ica40) =  y(iar36) * y(ihe4) * rate(irarag) &
            - y(ica40) * y(ihe4) * rate(ircaag) &
            - y(ica40) * rate(ircaga) &
            + y(iti44) * rate(irtiga)


      if (.not.deriva) then
       dydt(ica40)  = dydt(ica40) &
            + y(iar36) * y(ihe4) * rate(irarap)*(1.0d0-rate(iru1)) &
            - y(ica40) * rate(ircagp) * rate(iru1) &
            - y(ica40) * y(ihe4) * rate(ircaap)*(1.0d0-rate(irv1)) &
            + y(iti44) * rate(irtigp) * rate(irv1)

      else
       dydt(ica40)  = dydt(ica40) &
            + y(iar36)*y(ihe4)*(rate(irarap)*(1.0d0-ratdum(iru1)) &
                                      - ratdum(irarap)*rate(iru1)) &
            - y(ica40)*(ratdum(ircagp) * rate(iru1) + &
                           rate(ircagp) * ratdum(iru1)) &
            - y(ica40)*y(ihe4)*(rate(ircaap)*(1.0d0-ratdum(irv1)) &
                                      - ratdum(ircaap)*rate(irv1)) &
            + y(iti44)*(ratdum(irtigp) * rate(irv1) + &
                           rate(irtigp) * ratdum(irv1))
      end if



! ti44 reactions
      dydt(iti44) =  y(ica40) * y(ihe4) * rate(ircaag) &
            - y(iti44) * y(ihe4) * rate(irtiag) &
            - y(iti44) * rate(irtiga) &
            + y(icr48) * rate(ircrga)


      if (.not.deriva) then
       dydt(iti44)  = dydt(iti44) &
            + y(ica40) * y(ihe4) * rate(ircaap)*(1.0d0-rate(irv1)) &
            - y(iti44) * rate(irv1) * rate(irtigp) &
            - y(iti44) * y(ihe4) * rate(irtiap)*(1.0d0-rate(irw1)) &
            + y(icr48) * rate(irw1) * rate(ircrgp)
      else
       dydt(iti44)  = dydt(iti44) &
            + y(ica40)*y(ihe4)*(rate(ircaap)*(1.0d0-ratdum(irv1)) &
                                      - ratdum(ircaap)*rate(irv1)) &
            - y(iti44)*(ratdum(irv1) * rate(irtigp) + &
                           rate(irv1) * ratdum(irtigp)) &
            - y(iti44)*y(ihe4)*(rate(irtiap)*(1.0d0-ratdum(irw1)) &
                                      - ratdum(irtiap)*rate(irw1)) &
            + y(icr48)*(ratdum(irw1) * rate(ircrgp) + &
                           rate(irw1) * ratdum(ircrgp))
      end if


! cr48 reactions
      dydt(icr48) =  y(iti44) * y(ihe4) * rate(irtiag) &
            - y(icr48) * y(ihe4) * rate(ircrag) &
            - y(icr48) * rate(ircrga) &
            + y(ife52) * rate(irfega)

      if (.not.deriva) then
       dydt(icr48)  = dydt(icr48) &
            + y(iti44) * y(ihe4) * rate(irtiap)*(1.0d0-rate(irw1)) &
            - y(icr48) * rate(irw1) * rate(ircrgp) &
            - y(icr48) * y(ihe4) * rate(ircrap)*(1.0d0-rate(irx1)) &
            + y(ife52) * rate(irx1) * rate(irfegp)

      else
       dydt(icr48)  = dydt(icr48) &
            + y(iti44)*y(ihe4)*(rate(irtiap)*(1.0d0-ratdum(irw1)) &
                                      - ratdum(irtiap)*rate(irw1)) &
            - y(icr48)*(ratdum(irw1) * rate(ircrgp) + &
                           rate(irw1) * ratdum(ircrgp)) &
            - y(icr48)*y(ihe4)*(rate(ircrap)*(1.0d0-ratdum(irx1)) &
                                      - ratdum(ircrap)*rate(irx1)) &
            + y(ife52)*(ratdum(irx1) * rate(irfegp) + &
                           rate(irx1) * ratdum(irfegp))
      end if


! fe52 reactions
      dydt(ife52) =  y(icr48) * y(ihe4) * rate(ircrag) &
            - y(ife52) * y(ihe4) * rate(irfeag) &
            - y(ife52) * rate(irfega) &
            + y(ini56) * rate(irniga)

      if (.not.deriva) then
       dydt(ife52)  = dydt(ife52) &
            + y(icr48) * y(ihe4) * rate(ircrap)*(1.0d0-rate(irx1)) &
            - y(ife52) * rate(irx1) * rate(irfegp) &
            - y(ife52) * y(ihe4) * rate(irfeap)*(1.0d0-rate(iry1)) &
            + y(ini56) * rate(iry1) * rate(irnigp)

      else
       dydt(ife52)  = dydt(ife52) &
            + y(icr48)*y(ihe4)*(rate(ircrap)*(1.0d0-ratdum(irx1)) &
                                      - ratdum(ircrap)*rate(irx1)) &
            - y(ife52)*(ratdum(irx1) * rate(irfegp) + &
                           rate(irx1) * ratdum(irfegp)) &
            - y(ife52)*y(ihe4)*(rate(irfeap)*(1.0d0-ratdum(iry1)) &
                                      - ratdum(irfeap)*rate(iry1)) &
            + y(ini56)*(ratdum(iry1) * rate(irnigp) + &
                           rate(iry1) * ratdum(irnigp))
      end if


! ni56 reactions
      dydt(ini56) =  y(ife52) * y(ihe4) * rate(irfeag) &
            - y(ini56) * rate(irniga)

      if (.not.deriva) then
       dydt(ini56)  = dydt(ini56) &
            + y(ife52) * y(ihe4) * rate(irfeap)*(1.0d0-rate(iry1)) &
            - y(ini56) * rate(iry1) * rate(irnigp)

      else
       dydt(ini56)  = dydt(ini56) &
            + y(ife52)*y(ihe4)*(rate(irfeap)*(1.0d0-ratdum(iry1)) &
                                      - ratdum(irfeap)*rate(iry1)) &
            - y(ini56)*(ratdum(iry1) * rate(irnigp) + &
                           rate(iry1) * ratdum(irnigp))
      end if


      return
      end







      subroutine aprox13rat
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'

! this routine generates unscreened
! nuclear reaction rates for the aprox13 network.

! declare
      integer          i
      double precision rrate,drratedt,drratedd


! zero the rates
      do i=1,nrat
       ratraw(i) = 0.0d0
      enddo
      do i=1,nrat
       dratrawdt(i) = 0.0d0
      enddo
      do i=1,nrat
       dratrawdd(i) = 0.0d0
      enddo

      if (btemp .lt. 1.0e6) return


! get the temperature factors
      call tfactors(btemp)



! c12(a,g)o16
      call rate_c12ag(btemp,bden, &
           ratraw(ircag),dratrawdt(ircag),dratrawdd(ircag), &
           ratraw(iroga),dratrawdt(iroga),dratrawdd(iroga))

! triple alpha to c12
      call rate_tripalf(btemp,bden, &
           ratraw(ir3a),dratrawdt(ir3a),dratrawdd(ir3a), &
           ratraw(irg3a),dratrawdt(irg3a),dratrawdd(irg3a))

! c12 + c12
      call rate_c12c12(btemp,bden, &
           ratraw(ir1212),dratrawdt(ir1212),dratrawdd(ir1212), &
           rrate,drratedt,drratedd)

! c12 + o16
      call rate_c12o16(btemp,bden, &
           ratraw(ir1216),dratrawdt(ir1216),dratrawdd(ir1216), &
           rrate,drratedt,drratedd)

! o16 + o16
      call rate_o16o16(btemp,bden, &
           ratraw(ir1616),dratrawdt(ir1616),dratrawdd(ir1616), &
           rrate,drratedt,drratedd)

! o16(a,g)ne20
      call rate_o16ag(btemp,bden, &
           ratraw(iroag),dratrawdt(iroag),dratrawdd(iroag), &
           ratraw(irnega),dratrawdt(irnega),dratrawdd(irnega))

! ne20(a,g)mg24
      call rate_ne20ag(btemp,bden, &
           ratraw(irneag),dratrawdt(irneag),dratrawdd(irneag), &
           ratraw(irmgga),dratrawdt(irmgga),dratrawdd(irmgga))

! mg24(a,g)si28
      call rate_mg24ag(btemp,bden, &
           ratraw(irmgag),dratrawdt(irmgag),dratrawdd(irmgag), &
           ratraw(irsiga),dratrawdt(irsiga),dratrawdd(irsiga))

! mg24(a,p)al27
      call rate_mg24ap(btemp,bden, &
           ratraw(irmgap),dratrawdt(irmgap),dratrawdd(irmgap), &
           ratraw(iralpa),dratrawdt(iralpa),dratrawdd(iralpa))

! al27(p,g)si28
      call rate_al27pg(btemp,bden, &
           ratraw(iralpg),dratrawdt(iralpg),dratrawdd(iralpg), &
           ratraw(irsigp),dratrawdt(irsigp),dratrawdd(irsigp))

! si28(a,g)s32
      call rate_si28ag(btemp,bden, &
           ratraw(irsiag),dratrawdt(irsiag),dratrawdd(irsiag), &
           ratraw(irsga),dratrawdt(irsga),dratrawdd(irsga))

! si28(a,p)p31
      call rate_si28ap(btemp,bden, &
           ratraw(irsiap),dratrawdt(irsiap),dratrawdd(irsiap), &
           ratraw(irppa),dratrawdt(irppa),dratrawdd(irppa))

! p31(p,g)s32
      call rate_p31pg(btemp,bden, &
           ratraw(irppg),dratrawdt(irppg),dratrawdd(irppg), &
           ratraw(irsgp),dratrawdt(irsgp),dratrawdd(irsgp))

! s32(a,g)ar36
      call rate_s32ag(btemp,bden, &
           ratraw(irsag),dratrawdt(irsag),dratrawdd(irsag), &
           ratraw(irarga),dratrawdt(irarga),dratrawdd(irarga))

! s32(a,p)cl35
      call rate_s32ap(btemp,bden, &
           ratraw(irsap),dratrawdt(irsap),dratrawdd(irsap), &
           ratraw(irclpa),dratrawdt(irclpa),dratrawdd(irclpa))

! cl35(p,g)ar36
      call rate_cl35pg(btemp,bden, &
           ratraw(irclpg),dratrawdt(irclpg),dratrawdd(irclpg), &
           ratraw(irargp),dratrawdt(irargp),dratrawdd(irargp))

! ar36(a,g)ca40
      call rate_ar36ag(btemp,bden, &
           ratraw(irarag),dratrawdt(irarag),dratrawdd(irarag), &
           ratraw(ircaga),dratrawdt(ircaga),dratrawdd(ircaga))

! ar36(a,p)k39
      call rate_ar36ap(btemp,bden, &
           ratraw(irarap),dratrawdt(irarap),dratrawdd(irarap), &
           ratraw(irkpa),dratrawdt(irkpa),dratrawdd(irkpa))

! k39(p,g)ca40
      call rate_k39pg(btemp,bden, &
           ratraw(irkpg),dratrawdt(irkpg),dratrawdd(irkpg), &
           ratraw(ircagp),dratrawdt(ircagp),dratrawdd(ircagp))

! ca40(a,g)ti44
      call rate_ca40ag(btemp,bden, &
           ratraw(ircaag),dratrawdt(ircaag),dratrawdd(ircaag), &
           ratraw(irtiga),dratrawdt(irtiga),dratrawdd(irtiga))

! ca40(a,p)sc43
      call rate_ca40ap(btemp,bden, &
           ratraw(ircaap),dratrawdt(ircaap),dratrawdd(ircaap), &
           ratraw(irscpa),dratrawdt(irscpa),dratrawdd(irscpa))

! sc43(p,g)ti44
      call rate_sc43pg(btemp,bden, &
           ratraw(irscpg),dratrawdt(irscpg),dratrawdd(irscpg), &
           ratraw(irtigp),dratrawdt(irtigp),dratrawdd(irtigp))

! ti44(a,g)cr48
      call rate_ti44ag(btemp,bden, &
           ratraw(irtiag),dratrawdt(irtiag),dratrawdd(irtiag), &
           ratraw(ircrga),dratrawdt(ircrga),dratrawdd(ircrga))

! ti44(a,p)v47
      call rate_ti44ap(btemp,bden, &
           ratraw(irtiap),dratrawdt(irtiap),dratrawdd(irtiap), &
           ratraw(irvpa),dratrawdt(irvpa),dratrawdd(irvpa))

! v47(p,g)cr48
      call rate_v47pg(btemp,bden, &
           ratraw(irvpg),dratrawdt(irvpg),dratrawdd(irvpg), &
           ratraw(ircrgp),dratrawdt(ircrgp),dratrawdd(ircrgp))

! cr48(a,g)fe52
      call rate_cr48ag(btemp,bden, &
           ratraw(ircrag),dratrawdt(ircrag),dratrawdd(ircrag), &
           ratraw(irfega),dratrawdt(irfega),dratrawdd(irfega))

! cr48(a,p)mn51
      call rate_cr48ap(btemp,bden, &
           ratraw(ircrap),dratrawdt(ircrap),dratrawdd(ircrap), &
           ratraw(irmnpa),dratrawdt(irmnpa),dratrawdd(irmnpa))

! mn51(p,g)fe52
      call rate_mn51pg(btemp,bden, &
           ratraw(irmnpg),dratrawdt(irmnpg),dratrawdd(irmnpg), &
           ratraw(irfegp),dratrawdt(irfegp),dratrawdd(irfegp))

! fe52(a,g)ni56
      call rate_fe52ag(btemp,bden, &
           ratraw(irfeag),dratrawdt(irfeag),dratrawdd(irfeag), &
           ratraw(irniga),dratrawdt(irniga),dratrawdd(irniga))

! fe52(a,p)co55
      call rate_fe52ap(btemp,bden, &
           ratraw(irfeap),dratrawdt(irfeap),dratrawdd(irfeap), &
           ratraw(ircopa),dratrawdt(ircopa),dratrawdd(ircopa))

! co55(p,g)ni56
      call rate_co55pg(btemp,bden, &
           ratraw(ircopg),dratrawdt(ircopg),dratrawdd(ircopg), &
           ratraw(irnigp),dratrawdt(irnigp),dratrawdd(irnigp))


! write out the rates
!      write(6,133) btemp,bden
! 133  format(1x,1pe12.4)
!      do i=1,nrat
!       write(6,134) ratnam(i),ratraw(i)
! 134   format(1x,a,'  ',1pe14.6,1pe11.3,1pe14.6)
!      enddo
!      read(5,*)



! for a strict alpha chain with only (a,g) and (g,a) reactions
! shut down the (a,p) (p,g) and (g,p) (p,a) rates
!      ratraw(irmgap)    = 1.0d-40
!      dratrawdt(irmgap) = 0.0d0
!      dratrawdd(irmgap) = 0.0d0
!      ratraw(iralpa)    = 0.0d0
!      dratrawdt(iralpa) = 0.0d0
!      dratrawdd(iralpa) = 0.0d0
!      ratraw(iralpg) = 1.0d-40
!      dratrawdt(iralpg) = 0.0d0
!      dratrawdd(iralpg) = 0.0d0
!      ratraw(irsigp) = 1.0d-40
!      dratrawdt(irsigp) = 0.0d0
!      dratrawdd(irsigp) = 0.0d0

!      ratraw(irsiap) = 1.0d-40
!      dratrawdt(irsiap) = 0.0d0
!      dratrawdd(irsiap) = 0.0d0
!      ratraw(irppa)  = 0.0d0
!      dratrawdt(irppa)  = 0.0d0
!      dratrawdd(irppa)  = 0.0d0
!      ratraw(irppg)  = 1.0d-40
!      dratrawdt(irppg)  = 0.0d0
!      dratrawdd(irppg)  = 0.0d0
!      ratraw(irsgp)  = 1.0d-40
!      dratrawdt(irsgp)  = 0.0d0
!      dratrawdd(irsgp)  = 0.0d0c

!      ratraw(irsap)  = 1.0d-40
!      dratrawdt(irsap)  = 0.0d0
!      dratrawdd(irsap)  = 0.0d0
!      ratraw(irclpa) = 0.0d0
!      dratrawdt(irclpa) = 0.0d0
!      dratrawdd(irclpa) = 0.0d0
!      ratraw(irclpg) = 1.0d-40
!      dratrawdt(irclpg) = 0.0d0
!      dratrawdd(irclpg) = 0.0d0
!      ratraw(irargp) = 1.0d-40
!      dratrawdt(irargp) = 0.0d0
!      dratrawdd(irargp) = 0.0d0

!      ratraw(irarap) = 1.0d-40
!      dratrawdt(irarap) = 0.0d0
!      dratrawdd(irarap) = 0.0d0
!      ratraw(irkpa)  = 0.0d0
!      dratrawdt(irkpa)  = 0.0d0
!      dratrawdd(irkpa)  = 0.0d0
!      ratraw(irkpg)  = 1.0d-40
!      dratrawdt(irkpg)  = 0.0d0
!      dratrawdd(irkpg)  = 0.0d0
!      ratraw(ircagp) = 1.0d-40
!      dratrawdt(ircagp) = 0.0d0
!      dratrawdd(ircagp) = 0.0d0

!      ratraw(ircaap) = 1.0d-40
!      dratrawdt(ircaap) = 0.0d0
!      dratrawdd(ircaap) = 0.0d0
!      ratraw(irscpa) = 0.0d0
!      dratrawdt(irscpa) = 0.0d0
!      dratrawdd(irscpa) = 0.0d0
!      ratraw(irscpg) = 1.0d-40
!      dratrawdt(irscpg) = 0.0d0
!      dratrawdd(irscpg) = 0.0d0
!      ratraw(irtigp) = 1.0d-40
!      dratrawdt(irtigp) = 0.0d0
!      dratrawdd(irtigp) = 0.0d0

!      ratraw(irtiap) = 1.0d-40
!      dratrawdt(irtiap) = 0.0d0
!      dratrawdd(irtiap) = 0.0d0
!      ratraw(irvpa)  = 0.0d0
!      dratrawdt(irvpa)  = 0.0d0
!      dratrawdd(irvpa)  = 0.0d0
!      ratraw(irvpg)  = 1.0d-40
!      dratrawdt(irvpg)  = 0.0d0
!      dratrawdd(irvpg)  = 0.0d0
!      ratraw(ircrgp) = 1.0d-40
!      dratrawdt(ircrgp) = 0.0d0
!      dratrawdd(ircrgp) = 0.0d0

!      ratraw(ircrap) = 1.0d-40
!      dratrawdt(ircrap) = 0.0d0
!      dratrawdd(ircrap) = 0.0d0
!      ratraw(irmnpa) = 0.0d0
!      dratrawdt(irmnpa) = 0.0d0
!      dratrawdd(irmnpa) = 0.0d0
!      ratraw(irmnpg) = 1.0d-40
!      dratrawdt(irmnpg) = 0.0d0
!      dratrawdd(irmnpg) = 0.0d0
!      ratraw(irfegp) = 1.0d-40
!      dratrawdt(irfegp) = 0.0d0
!      dratrawdd(irfegp) = 0.0d0

!      ratraw(irfeap) = 1.0d-40
!      dratrawdt(irfeap) = 0.0d0
!      dratrawdd(irfeap) = 0.0d0
!      ratraw(ircopa) = 0.0d0
!      dratrawdt(ircopa) = 0.0d0
!      dratrawdd(ircopa) = 0.0d0
!      ratraw(ircopg) = 1.0d-40
!      dratrawdt(ircopg) = 0.0d0
!      dratrawdd(ircopg) = 0.0d0
!      ratraw(irnigp) = 1.0d-40
!      dratrawdt(irnigp) = 0.0d0
!      dratrawdd(irnigp) = 0.0d0

      return
      end






      subroutine screen_aprox13(y)
      include 'implno.dek'
      include 'burn_common.dek'
      include 'network.dek'

! this routine computes the screening factors
! and applies them to the raw reaction rates,
! producing the final reaction rates used by the
! right hand sides and jacobian matrix elements

! this routine assumes screen_on = 1 or = 0 has been set at a higher
! level presumably in the top level driver


! declare the pass
      double precision y(*)

! local variables
      integer          i,jscr,init
      double precision sc1a,sc1adt,sc1add,sc2a,sc2adt,sc2add, &
                       sc3a,sc3adt,sc3add

      double precision abar,zbar,z2bar,ytot1,zbarxx,z2barxx, &
                       denom,denomdt,denomdd, &
                       r1,r1dt,r1dd,s1,s1dt,s1dd,t1,t1dt,t1dd, &
                       u1,u1dt,u1dd,v1,v1dt,v1dd,w1,w1dt,w1dd, &
                       x1,x1dt,x1dd,y1,y1dt,y1dd,zz

      data             init/1/


! initialize
      do i=1,nrat
       ratdum(i)    = ratraw(i)
       dratdumdt(i) = dratrawdt(i)
       dratdumdd(i) = dratrawdd(i)
       scfac(i)     = 1.0d0
       dscfacdt(i)  = 0.0d0
       dscfacdt(i)  = 0.0d0
      end do



! if screening is on
      if (screen_on .ne. 0) then



! with the passed composition, compute abar,zbar and other variables
       zbarxx  = 0.0d0
       z2barxx = 0.0d0
       ytot1   = 0.0d0
       do i=ionbeg,ionend
        ytot1    = ytot1 + y(i)
        zbarxx   = zbarxx + zion(i) * y(i)
        z2barxx  = z2barxx + zion(i) * zion(i) * y(i)
       enddo
       abar   = 1.0d0/ytot1
       zbar   = zbarxx * abar
       z2bar  = z2barxx * abar



! first the always fun triple alpha and its inverse
      jscr = 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ihe4),aion(ihe4),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ihe4),aion(ihe4),4.0d0,8.0d0, &
                   jscr,init,sc2a,sc2adt,sc2add)

      sc3a   = sc1a * sc2a
      sc3adt = sc1adt*sc2a + sc1a*sc2adt
      sc3add = sc1add*sc2a + sc1a*sc2add

      ratdum(ir3a)    = ratraw(ir3a) * sc3a
      dratdumdt(ir3a) = dratrawdt(ir3a)*sc3a + ratraw(ir3a)*sc3adt
      dratdumdd(ir3a) = dratrawdd(ir3a)*sc3a + ratraw(ir3a)*sc3add

      scfac(ir3a)     = sc3a
      dscfacdt(ir3a)  = sc3adt
      dscfacdd(ir3a)  = sc3add


! c12 to o16
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ic12),aion(ic12),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ircag)     = ratraw(ircag) * sc1a
      dratdumdt(ircag)  = dratrawdt(ircag)*sc1a + ratraw(ircag)*sc1adt
      dratdumdd(ircag)  = dratrawdd(ircag)*sc1a + ratraw(ircag)*sc1add

      scfac(ircag)      = sc1a
      dscfacdt(ircag)   = sc1adt
      dscfacdt(ircag)   = sc1add


! c12 + c12
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ic12),aion(ic12),zion(ic12),aion(ic12), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ir1212)    = ratraw(ir1212) * sc1a
      dratdumdt(ir1212) = dratrawdt(ir1212)*sc1a + ratraw(ir1212)*sc1adt
      dratdumdd(ir1212) = dratrawdd(ir1212)*sc1a + ratraw(ir1212)*sc1add

      scfac(ir1212)     = sc1a
      dscfacdt(ir1212)  = sc1adt
      dscfacdd(ir1212)  = sc1add



! c12 + o16
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ic12),aion(ic12),zion(io16),aion(io16), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ir1216)    = ratraw(ir1216) * sc1a
      dratdumdt(ir1216) = dratrawdt(ir1216)*sc1a + ratraw(ir1216)*sc1adt
      dratdumdd(ir1216) = dratrawdd(ir1216)*sc1a + ratraw(ir1216)*sc1add

      scfac(ir1216)     = sc1a
      dscfacdt(ir1216)  = sc1adt
      dscfacdd(ir1216)  = sc1add



! o16 + o16
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(io16),aion(io16),zion(io16),aion(io16), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ir1616)    = ratraw(ir1616) * sc1a
      dratdumdt(ir1616) = dratrawdt(ir1616)*sc1a + ratraw(ir1616)*sc1adt
      dratdumdd(ir1616) = dratrawdd(ir1616)*sc1a + ratraw(ir1616)*sc1add

      scfac(ir1616)     = sc1a
      dscfacdt(ir1616)  = sc1adt
      dscfacdd(ir1616)  = sc1add



! o16 to ne20
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(io16),aion(io16),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(iroag)    = ratraw(iroag) * sc1a
      dratdumdt(iroag) = dratrawdt(iroag)*sc1a + ratraw(iroag)*sc1adt
      dratdumdd(iroag) = dratrawdd(iroag)*sc1a + ratraw(iroag)*sc1add

      scfac(iroag)     = sc1a
      dscfacdt(iroag)  = sc1adt
      dscfacdd(iroag)  = sc1add



! ne20 to mg24
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ine20),aion(ine20),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irneag)    = ratraw(irneag) * sc1a
      dratdumdt(irneag) = dratrawdt(irneag)*sc1a + ratraw(irneag)*sc1adt
      dratdumdd(irneag) = dratrawdd(irneag)*sc1a + ratraw(irneag)*sc1add

      scfac(irneag)     = sc1a
      dscfacdt(irneag)  = sc1adt
      dscfacdd(irneag)  = sc1add


! mg24 to si28
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(img24),aion(img24),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irmgag)    = ratraw(irmgag) * sc1a
      dratdumdt(irmgag) = dratrawdt(irmgag)*sc1a + ratraw(irmgag)*sc1adt
      dratdumdd(irmgag) = dratrawdd(irmgag)*sc1a + ratraw(irmgag)*sc1add

      scfac(irmgag)     = sc1a
      dscfacdt(irmgag)  = sc1adt
      dscfacdd(irmgag)  = sc1add

      ratdum(irmgap)    = ratraw(irmgap) * sc1a
      dratdumdt(irmgap) = dratrawdt(irmgap)*sc1a + ratraw(irmgap)*sc1adt
      dratdumdd(irmgap) = dratrawdd(irmgap)*sc1a + ratraw(irmgap)*sc1add

      scfac(irmgap)     = sc1a
      dscfacdt(irmgap)  = sc1adt
      dscfacdd(irmgap)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   13.0d0,27.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(iralpa)    = ratraw(iralpa) * sc1a
      dratdumdt(iralpa) = dratrawdt(iralpa)*sc1a + ratraw(iralpa)*sc1adt
      dratdumdd(iralpa) = dratrawdd(iralpa)*sc1a + ratraw(iralpa)*sc1add

      scfac(iralpa)     = sc1a
      dscfacdt(iralpa)  = sc1adt
      dscfacdd(iralpa)  = sc1add

      ratdum(iralpg)    = ratraw(iralpg) * sc1a
      dratdumdt(iralpg) = dratrawdt(iralpg)*sc1a + ratraw(iralpg)*sc1adt
      dratdumdd(iralpg) = dratrawdd(iralpg)*sc1a + ratraw(iralpg)*sc1add

      scfac(iralpg)     = sc1a
      dscfacdt(iralpg)  = sc1adt
      dscfacdd(iralpg)  = sc1add



! si28 to s32
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(isi28),aion(isi28),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irsiag)    = ratraw(irsiag) * sc1a
      dratdumdt(irsiag) = dratrawdt(irsiag)*sc1a + ratraw(irsiag)*sc1adt
      dratdumdd(irsiag) = dratrawdd(irsiag)*sc1a + ratraw(irsiag)*sc1add

      scfac(irsiag)     = sc1a
      dscfacdt(irsiag)  = sc1adt
      dscfacdd(irsiag)  = sc1add


      ratdum(irsiap)    = ratraw(irsiap) * sc1a
      dratdumdt(irsiap) = dratrawdt(irsiap)*sc1a + ratraw(irsiap)*sc1adt
      dratdumdd(irsiap) = dratrawdd(irsiap)*sc1a + ratraw(irsiap)*sc1add

      scfac(irsiap)     = sc1a
      dscfacdt(irsiap)  = sc1adt
      dscfacdd(irsiap)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   15.0d0,31.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irppa)     = ratraw(irppa) * sc1a
      dratdumdt(irppa)  = dratrawdt(irppa)*sc1a  + ratraw(irppa)*sc1adt
      dratdumdd(irppa)  = dratrawdd(irppa)*sc1a  + ratraw(irppa)*sc1add

      scfac(irppa)      = sc1a
      dscfacdt(irppa)   = sc1adt
      dscfacdd(irppa)   = sc1add

      ratdum(irppg)     = ratraw(irppg) * sc1a
      dratdumdt(irppg)  = dratrawdt(irppg)*sc1a + ratraw(irppg)*sc1adt
      dratdumdd(irppg)  = dratrawdd(irppg)*sc1a + ratraw(irppg)*sc1add

      scfac(irppg)      = sc1a
      dscfacdt(irppg)   = sc1adt
      dscfacdd(irppg)   = sc1add



! s32 to ar36
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(is32),aion(is32),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irsag)     = ratraw(irsag) * sc1a
      dratdumdt(irsag)  = dratrawdt(irsag)*sc1a + ratraw(irsag)*sc1adt
      dratdumdd(irsag)  = dratrawdd(irsag)*sc1a + ratraw(irsag)*sc1add

      scfac(irsag)      = sc1a
      dscfacdt(irsag)   = sc1adt
      dscfacdd(irsag)   = sc1add

      ratdum(irsap)     = ratraw(irsap) * sc1a
      dratdumdt(irsap)  = dratrawdt(irsap)*sc1a + ratraw(irsap)*sc1adt
      dratdumdd(irsap)  = dratrawdd(irsap)*sc1a + ratraw(irsap)*sc1add

      scfac(irsap)      = sc1a
      dscfacdt(irsap)   = sc1adt
      dscfacdd(irsap)   = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   17.0d0,35.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irclpa)    = ratraw(irclpa) * sc1a
      dratdumdt(irclpa) = dratrawdt(irclpa)*sc1a + ratraw(irclpa)*sc1adt
      dratdumdd(irclpa) = dratrawdd(irclpa)*sc1a + ratraw(irclpa)*sc1add

      scfac(irclpa)     = sc1a
      dscfacdt(irclpa)  = sc1adt
      dscfacdt(irclpa)  = sc1add

      ratdum(irclpg)    = ratraw(irclpg) * sc1a
      dratdumdt(irclpg) = dratrawdt(irclpg)*sc1a + ratraw(irclpg)*sc1adt
      dratdumdd(irclpg) = dratrawdd(irclpg)*sc1a + ratraw(irclpg)*sc1add

      scfac(irclpg)     = sc1a
      dscfacdt(irclpg)  = sc1adt
      dscfacdd(irclpg)  = sc1add



! ar36 to ca40
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(iar36),aion(iar36),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irarag)    = ratraw(irarag) * sc1a
      dratdumdt(irarag) = dratrawdt(irarag)*sc1a + ratraw(irarag)*sc1adt
      dratdumdd(irarag) = dratrawdd(irarag)*sc1a + ratraw(irarag)*sc1add

      scfac(irarag)     = sc1a
      dscfacdt(irarag)  = sc1adt
      dscfacdd(irarag)  = sc1add

      ratdum(irarap)    = ratraw(irarap) * sc1a
      dratdumdt(irarap) = dratrawdt(irarap)*sc1a + ratraw(irarap)*sc1adt
      dratdumdd(irarap) = dratrawdd(irarap)*sc1a + ratraw(irarap)*sc1add

      scfac(irarap)     = sc1a
      dscfacdt(irarap)  = sc1adt
      dscfacdd(irarap)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   19.0d0,39.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irkpa)     = ratraw(irkpa) * sc1a
      dratdumdt(irkpa)  = dratrawdt(irkpa)*sc1a  + ratraw(irkpa)*sc1adt
      dratdumdd(irkpa)  = dratrawdd(irkpa)*sc1a  + ratraw(irkpa)*sc1add

      scfac(irkpa)      = sc1a
      dscfacdt(irkpa)   = sc1adt
      dscfacdd(irkpa)   = sc1add

      ratdum(irkpg)     = ratraw(irkpg) * sc1a
      dratdumdt(irkpg)  = dratrawdt(irkpg)*sc1a  + ratraw(irkpg)*sc1adt
      dratdumdd(irkpg)  = dratrawdd(irkpg)*sc1a  + ratraw(irkpg)*sc1add

      scfac(irkpg)      = sc1a
      dscfacdt(irkpg)   = sc1adt
      dscfacdd(irkpg)   = sc1add



! ca40 to ti44
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ica40),aion(ica40),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ircaag)    = ratraw(ircaag) * sc1a
      dratdumdt(ircaag) = dratrawdt(ircaag)*sc1a + ratraw(ircaag)*sc1adt
      dratdumdd(ircaag) = dratrawdd(ircaag)*sc1a + ratraw(ircaag)*sc1add

      scfac(ircaag)     = sc1a
      dscfacdt(ircaag)  = sc1adt
      dscfacdd(ircaag)  = sc1add

      ratdum(ircaap)    = ratraw(ircaap) * sc1a
      dratdumdt(ircaap) = dratrawdt(ircaap)*sc1a + ratraw(ircaap)*sc1adt
      dratdumdd(ircaap) = dratrawdd(ircaap)*sc1a + ratraw(ircaap)*sc1add

      scfac(ircaap)     = sc1a
      dscfacdt(ircaap)  = sc1adt
      dscfacdd(ircaap)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   21.0d0,43.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irscpa)    = ratraw(irscpa) * sc1a
      dratdumdt(irscpa) = dratrawdt(irscpa)*sc1a + ratraw(irscpa)*sc1adt
      dratdumdd(irscpa) = dratrawdd(irscpa)*sc1a + ratraw(irscpa)*sc1add

      scfac(irscpa)     = sc1a
      dscfacdt(irscpa)  = sc1adt
      dscfacdd(irscpa)  = sc1add

      ratdum(irscpg)    = ratraw(irscpg) * sc1a
      dratdumdt(irscpg) = dratrawdt(irscpg)*sc1a + ratraw(irscpg)*sc1adt
      dratdumdd(irscpg) = dratrawdd(irscpg)*sc1a + ratraw(irscpg)*sc1add

      scfac(irscpg)     = sc1a
      dscfacdt(irscpg)  = sc1adt
      dscfacdd(irscpg)  = sc1add



! ti44 to cr48
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(iti44),aion(iti44),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irtiag)    = ratraw(irtiag) * sc1a
      dratdumdt(irtiag) = dratrawdt(irtiag)*sc1a + ratraw(irtiag)*sc1adt
      dratdumdd(irtiag) = dratrawdd(irtiag)*sc1a + ratraw(irtiag)*sc1add

      scfac(irtiag)     = sc1a
      dscfacdt(irtiag)  = sc1adt
      dscfacdd(irtiag)  = sc1add

      ratdum(irtiap)    = ratraw(irtiap) * sc1a
      dratdumdt(irtiap) = dratrawdt(irtiap)*sc1a + ratraw(irtiap)*sc1adt
      dratdumdd(irtiap) = dratrawdd(irtiap)*sc1a + ratraw(irtiap)*sc1add

      scfac(irtiap)  = sc1a
      dscfacdt(irtiap)  = sc1adt
      dscfacdd(irtiap)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   23.0d0,47.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irvpa)     = ratraw(irvpa) * sc1a
      dratdumdt(irvpa)  = dratrawdt(irvpa)*sc1a  + ratraw(irvpa)*sc1adt
      dratdumdd(irvpa)  = dratrawdd(irvpa)*sc1a  + ratraw(irvpa)*sc1add

      scfac(irvpa)      = sc1a
      dscfacdt(irvpa)   = sc1adt
      dscfacdd(irvpa)   = sc1add

      ratdum(irvpg)     = ratraw(irvpg) * sc1a
      dratdumdt(irvpg)  = dratrawdt(irvpg)*sc1a  + ratraw(irvpg)*sc1adt
      dratdumdd(irvpg)  = dratrawdd(irvpg)*sc1a  + ratraw(irvpg)*sc1add

      scfac(irvpg)      = sc1a
      dscfacdt(irvpg)   = sc1adt
      dscfacdd(irvpg)   = sc1add



! cr48 to fe52
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(icr48),aion(icr48),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ircrag)    = ratraw(ircrag) * sc1a
      dratdumdt(ircrag) = dratrawdt(ircrag)*sc1a + ratraw(ircrag)*sc1adt
      dratdumdd(ircrag) = dratrawdd(ircrag)*sc1a + ratraw(ircrag)*sc1add

      scfac(ircrag)     = sc1a
      dscfacdt(ircrag)  = sc1adt
      dscfacdd(ircrag)  = sc1add

      ratdum(ircrap)    = ratraw(ircrap) * sc1a
      dratdumdt(ircrap) = dratrawdt(ircrap)*sc1a + ratraw(ircrap)*sc1adt
      dratdumdd(ircrap) = dratrawdd(ircrap)*sc1a + ratraw(ircrap)*sc1add

      scfac(ircrap)     = sc1a
      dscfacdt(ircrap)  = sc1adt
      dscfacdd(ircrap)  = sc1add


      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   25.0d0,51.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irmnpa)    = ratraw(irmnpa) * sc1a
      dratdumdt(irmnpa) = dratrawdt(irmnpa)*sc1a + ratraw(irmnpa)*sc1adt
      dratdumdd(irmnpa) = dratrawdd(irmnpa)*sc1a + ratraw(irmnpa)*sc1add

      scfac(irmnpa)     = sc1a
      dscfacdt(irmnpa)  = sc1adt
      dscfacdd(irmnpa)  = sc1add

      ratdum(irmnpg)    = ratraw(irmnpg) * sc1a
      dratdumdt(irmnpg) = dratrawdt(irmnpg)*sc1a + ratraw(irmnpg)*sc1adt
      dratdumdd(irmnpg) = dratrawdd(irmnpg)*sc1a + ratraw(irmnpg)*sc1add

      scfac(irmnpg)     = sc1a
      dscfacdt(irmnpg)  = sc1adt
      dscfacdd(irmnpg)  = sc1add


! fe52 to ni56
      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   zion(ife52),aion(ife52),zion(ihe4),aion(ihe4), &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(irfeag)    = ratraw(irfeag) * sc1a
      dratdumdt(irfeag) = dratrawdt(irfeag)*sc1a + ratraw(irfeag)*sc1adt
      dratdumdd(irfeag) = dratrawdd(irfeag)*sc1a + ratraw(irfeag)*sc1add

      scfac(irfeag)     = sc1a
      dscfacdt(irfeag)  = sc1adt
      dscfacdd(irfeag)  = sc1add

      ratdum(irfeap) = ratraw(irfeap) * sc1a
      dratdumdt(irfeap) = dratrawdt(irfeap)*sc1a + ratraw(irfeap)*sc1adt
      dratdumdd(irfeap) = dratrawdd(irfeap)*sc1a + ratraw(irfeap)*sc1add

      scfac(irfeap)     = sc1a
      dscfacdt(irfeap)  = sc1adt
      dscfacdd(irfeap)  = sc1add

      jscr = jscr + 1
      call screen5(btemp,bden,zbar,abar,z2bar, &
                   27.0d0,55.0d0,1.0d0,1.0d0, &
                   jscr,init,sc1a,sc1adt,sc1add)

      ratdum(ircopa)    = ratraw(ircopa) * sc1a
      dratdumdt(ircopa) = dratrawdt(ircopa)*sc1a + ratraw(ircopa)*sc1adt
      dratdumdd(ircopa) = dratrawdd(ircopa)*sc1a + ratraw(ircopa)*sc1add

      scfac(ircopa)     = sc1a
      dscfacdt(ircopa)  = sc1adt
      dscfacdd(ircopa)  = sc1add

      ratdum(ircopg)    = ratraw(ircopg) * sc1a
      dratdumdt(ircopg) = dratrawdt(ircopg)*sc1a + ratraw(ircopg)*sc1adt
      dratdumdd(ircopg) = dratrawdd(ircopg)*sc1a + ratraw(ircopg)*sc1add

      scfac(ircopg)     = sc1a
      dscfacdt(ircopg)  = sc1adt
      dscfacdd(ircopg)  = sc1add


! reset the screen initialization flag
      init = 0

! end of screening block if
      end if


! now form those lovely dummy proton link rates

      ratdum(irr1)     = 0.0d0
      dratdumdt(irr1)  = 0.0d0
      dratdumdd(irr1)  = 0.0d0
      denom    = ratdum(iralpa) + ratdum(iralpg)
      denomdt  = dratdumdt(iralpa) + dratdumdt(iralpg)
      denomdd  = dratdumdd(iralpa) + dratdumdd(iralpg)
      if (denom .ne. 0.0) then
       zz = 1.0d0/denom
       ratdum(irr1)    = ratdum(iralpa)*zz
       dratdumdt(irr1) = (dratdumdt(iralpa) - ratdum(irr1)*denomdt)*zz
       dratdumdd(irr1) = (dratdumdd(iralpa) - ratdum(irr1)*denomdd)*zz
      end if

      ratdum(irs1)     = 0.0d0
      dratdumdt(irs1)  = 0.0d0
      dratdumdd(irs1)  = 0.0d0
      denom    = ratdum(irppa) + ratdum(irppg)
      denomdt  = dratdumdt(irppa) + dratdumdt(irppg)
      denomdd  = dratdumdd(irppa) + dratdumdd(irppg)
      if (denom .ne. 0.0) then
       zz = 1.0d0/denom
       ratdum(irs1)    = ratdum(irppa)*zz
       dratdumdt(irs1) = (dratdumdt(irppa) - ratdum(irs1)*denomdt)*zz
       dratdumdd(irs1) = (dratdumdd(irppa) - ratdum(irs1)*denomdd)*zz
      end if

      ratdum(irt1)     = 0.0d0
      dratdumdt(irt1)  = 0.0d0
      dratdumdd(irt1)  = 0.0d0
      denom    = ratdum(irclpa) + ratdum(irclpg)
      denomdt  = dratdumdt(irclpa) + dratdumdt(irclpg)
      denomdd  = dratdumdd(irclpa) + dratdumdd(irclpg)
      if (denom .ne. 0.0) then
       zz = 1.0d0/denom
       ratdum(irt1)    = ratdum(irclpa)*zz
       dratdumdt(irt1) = (dratdumdt(irclpa) - ratdum(irt1)*denomdt)*zz
       dratdumdd(irt1) = (dratdumdd(irclpa) - ratdum(irt1)*denomdd)*zz
      end if

      ratdum(iru1)     = 0.0d0
      dratdumdt(iru1)  = 0.0d0
      dratdumdd(iru1)  = 0.0d0
      denom    = ratdum(irkpa) + ratdum(irkpg)
      denomdt  = dratdumdt(irkpa) + dratdumdt(irkpg)
      denomdd  = dratdumdd(irkpa) + dratdumdd(irkpg)
      if (denom .ne. 0.0) then
       zz   = 1.0d0/denom
       ratdum(iru1)   = ratdum(irkpa)*zz
       dratdumdt(iru1) = (dratdumdt(irkpa) - ratdum(iru1)*denomdt)*zz
       dratdumdd(iru1) = (dratdumdd(irkpa) - ratdum(iru1)*denomdd)*zz
      end if

      ratdum(irv1)     = 0.0d0
      dratdumdt(irv1)  = 0.0d0
      dratdumdd(irv1)  = 0.0d0
      denom    = ratdum(irscpa) + ratdum(irscpg)
      denomdt  = dratdumdt(irscpa) + dratdumdt(irscpg)
      denomdd  = dratdumdd(irscpa) + dratdumdd(irscpg)
      if (denom .ne. 0.0) then
       zz  = 1.0d0/denom
       ratdum(irv1)    = ratdum(irscpa)*zz
       dratdumdt(irv1) = (dratdumdt(irscpa) - ratdum(irv1)*denomdt)*zz
       dratdumdd(irv1) = (dratdumdd(irscpa) - ratdum(irv1)*denomdd)*zz
      end if

      ratdum(irw1)    = 0.0d0
      dratdumdt(irw1) = 0.0d0
      dratdumdd(irw1) = 0.0d0
      denom    = ratdum(irvpa) + ratdum(irvpg)
      denomdt  = dratdumdt(irvpa) + dratdumdt(irvpg)
      denomdd  = dratdumdd(irvpa) + dratdumdd(irvpg)
      if (denom .ne. 0.0) then
       zz = 1.0d0/denom
       ratdum(irw1)    = ratdum(irvpa)*zz
       dratdumdt(irw1) = (dratdumdt(irvpa) - ratdum(irw1)*denomdt)*zz
       dratdumdd(irw1) = (dratdumdd(irvpa) - ratdum(irw1)*denomdd)*zz
      end if

      ratdum(irx1)    = 0.0d0
      dratdumdt(irx1) = 0.0d0
      dratdumdd(irx1) = 0.0d0
      denom    = ratdum(irmnpa) + ratdum(irmnpg)
      denomdt  = dratdumdt(irmnpa) + dratdumdt(irmnpg)
      denomdd  = dratdumdd(irmnpa) + dratdumdd(irmnpg)
      if (denom .ne. 0.0) then
       zz = 1.0d0/denom
       ratdum(irx1)    = ratdum(irmnpa)*zz
       dratdumdt(irx1) = (dratdumdt(irmnpa) - ratdum(irx1)*denomdt)*zz
       dratdumdd(irx1) = (dratdumdd(irmnpa) - ratdum(irx1)*denomdd)*zz
      endif

      ratdum(iry1)    = 0.0d0
      dratdumdt(iry1) = 0.0d0
      dratdumdd(iry1) = 0.0d0
      denom    = ratdum(ircopa) + ratdum(ircopg)
      denomdt  = dratdumdt(ircopa) + dratdumdt(ircopg)
      denomdd  = dratdumdd(ircopa) + dratdumdd(ircopg)
      if (denom .ne. 0.0) then
       zz = 1.0d0/denom
       ratdum(iry1)    = ratdum(ircopa)*zz
       dratdumdt(iry1) = (dratdumdt(ircopa) - ratdum(iry1)*denomdt)*zz
       dratdumdd(iry1) = (dratdumdd(ircopa) - ratdum(iry1)*denomdd)*zz
      end if





! debugs
!      do i=1,nrat
!       if (ratdum(i) .lt. 0.0) then
!        write(6,110) i,ratnam(i),ratraw(i),scfac(i),ratdum(i)
! 110    format(1x,i4,' ',a,' ',1p3e12.4)
!        stop 'negative rate'
!       end if
!      enddo

!      do i=1,4
!       write(6,111) i,ratnam(i),ratraw(i),scfac(i),ratdum(i)
!       write(6,111) i,ratnam(i),ratdum(i),dratdumdt(i),dratdumdd(i)
! 111   format(1x,i4,' ',a,' ',1p3e14.6)
!      enddo
!      read(5,*)

      return
      end







      subroutine init_aprox13
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'
!
! this routine initializes stuff for the aprox13 network
!
! declare
      integer          i
      double precision  mev2erg,mev2gr
      parameter        (mev2erg = ev2erg*1.0d6, &
                        mev2gr  = mev2erg/clight**2)


! for easy zeroing of the isotope pointers
      integer          isotp(nisotp)
      equivalence      (isotp(1),ih1)

! for easy zeroing of the rate pointers
      integer          rts(numrates)
      equivalence      (rts(1),ir3a)


! zero all the isotope pointers
      do i=1,nisotp
       isotp(i)   = 0
      enddo

! zero all the rate pointers
      do i=1,numrates
       rts(i) = 0
      enddo


! set the size of the network and the number of rates
      idnet   = idaprox13

      ionmax  = 13

!      iener   = 1
!      itemp   = 2
!      iden    = ionmax + 3
!      ivelx   = ionmax + 4
!      iposx   = ionmax + 5
!      neqs    = iposx


      ionmax  = 13
      iener   = ionmax + 1
      itemp   = ionmax + 2
      iden    = ionmax + 3
      ivelx   = ionmax + 4
      iposx   = ionmax + 5
      neqs    = iposx


      nrat    = 67
      netname = 'aprox13'


! set the id numbers of the elements

!      ionbeg = 3
!      ionend = 15

      ionbeg = 1
      ionend = 13

      ihe4  = ionbeg
      ic12  = ionbeg + 1
      io16  = ionbeg + 2
      ine20 = ionbeg + 3
      img24 = ionbeg + 4
      isi28 = ionbeg + 5
      is32  = ionbeg + 6
      iar36 = ionbeg + 7
      ica40 = ionbeg + 8
      iti44 = ionbeg + 9
      icr48 = ionbeg + 10
      ife52 = ionbeg + 11
      ini56 = ionbeg + 12


! set the names of the elements
      ionam(ihe4)  = 'he4  '
      ionam(ic12)  = 'c12  '
      ionam(io16)  = 'o16  '
      ionam(ine20) = 'ne20 '
      ionam(img24) = 'mg24 '
      ionam(isi28) = 'si28 '
      ionam(is32)  = 's32  '
      ionam(iar36) = 'ar36 '
      ionam(ica40) = 'ca40 '
      ionam(iti44) = 'ti44 '
      ionam(icr48) = 'cr48 '
      ionam(ife52) = 'fe52 '
      ionam(ini56) = 'ni56 '

      ionam(iener) = 'ener '
      ionam(itemp) = 'temp '
      ionam(iden)  = 'den  '
      ionam(ivelx) = 'velx '
      ionam(iposx) = 'posx '


! set the number of nucleons in the element
      aion(ihe4)  = 4.0d0
      aion(ic12)  = 12.0d0
      aion(io16)  = 16.0d0
      aion(ine20) = 20.0d0
      aion(img24) = 24.0d0
      aion(isi28) = 28.0d0
      aion(is32)  = 32.0d0
      aion(iar36) = 36.0d0
      aion(ica40) = 40.0d0
      aion(iti44) = 44.0d0
      aion(icr48) = 48.0d0
      aion(ife52) = 52.0d0
      aion(ini56) = 56.0d0

! set the number of protons in the element
      zion(ihe4)  = 2.0d0
      zion(ic12)  = 6.0d0
      zion(io16)  = 8.0d0
      zion(ine20) = 10.0d0
      zion(img24) = 12.0d0
      zion(isi28) = 14.0d0
      zion(is32)  = 16.0d0
      zion(iar36) = 18.0d0
      zion(ica40) = 20.0d0
      zion(iti44) = 22.0d0
      zion(icr48) = 24.0d0
      zion(ife52) = 26.0d0
      zion(ini56) = 28.0d0


! set the binding energy of the element
      bion(ihe4)  =  28.29603d0
      bion(ic12)  =  92.16294d0
      bion(io16)  = 127.62093d0
      bion(ine20) = 160.64788d0
      bion(img24) = 198.25790d0
      bion(isi28) = 236.53790d0
      bion(is32)  = 271.78250d0
      bion(iar36) = 306.72020d0
      bion(ica40) = 342.05680d0
      bion(iti44) = 375.47720d0
      bion(icr48) = 411.46900d0
      bion(ife52) = 447.70800d0
      bion(ini56) = 484.00300d0


! set the number of neutrons and mass
      do i=ionbeg,ionend
       nion(i) = aion(i) - zion(i)
      enddo

! mass of each isotope
      do i = ionbeg,ionend
       mion(i) = nion(i)*mn + zion(i)*(mp+me) - bion(i)*mev2gr
      enddo

! molar mass
      do i = ionbeg,ionend
       wion(i) = avo * mion(i)
      enddo

! a common approximation
      do i = ionbeg,ionend
       wion(i) = aion(i)
      enddo


! set the partition functions - statistical weights, ground-state only here
      do i=ionbeg,ionend
       wpart(i) = 1.0d0
      enddo

! set the id numbers of the reaction rates
      ir3a   = 1
      irg3a  = 2
      ircag  = 3
      iroga  = 4
      ir1212 = 5
      ir1216 = 6
      ir1616 = 7
      iroag  = 8
      irnega = 9
      irneag = 10
      irmgga = 11
      irmgag = 12
      irsiga = 13
      irmgap = 14
      iralpa = 15
      iralpg = 16
      irsigp = 17
      irsiag = 18
      irsga  = 19
      irsiap = 20
      irppa  = 21
      irppg  = 22
      irsgp  = 23
      irsag  = 24
      irarga = 25
      irsap  = 26
      irclpa = 27
      irclpg = 28
      irargp = 29
      irarag = 30
      ircaga = 31
      irarap = 32
      irkpa  = 33
      irkpg  = 34
      ircagp = 35
      ircaag = 36
      irtiga = 37
      ircaap = 38
      irscpa = 39
      irscpg = 40
      irtigp = 41
      irtiag = 42
      ircrga = 43
      irtiap = 44
      irvpa  = 45
      irvpg  = 46
      ircrgp = 47
      ircrag = 48
      irfega = 49
      ircrap = 50
      irmnpa = 51
      irmnpg = 52
      irfegp = 53
      irfeag = 54
      irniga = 55
      irfeap = 56
      ircopa = 57
      ircopg = 58
      irnigp = 59

      irr1   = 60
      irs1   = 61
      irt1   = 62
      iru1   = 63
      irv1   = 64
      irw1   = 65
      irx1   = 66
      iry1   = 67


! set the names of the reaction rates
      ratnam(ir3a)   = 'r3a  '
      ratnam(irg3a)  = 'rg3a '
      ratnam(ircag)  = 'rcag '
      ratnam(ir1212) = 'r1212'
      ratnam(ir1216) = 'r1216'
      ratnam(ir1616) = 'r1616'
      ratnam(iroga)  = 'roga '
      ratnam(iroag)  = 'roag '
      ratnam(irnega) = 'rnega'
      ratnam(irneag) = 'rneag'
      ratnam(irmgga) = 'rmgga'
      ratnam(irmgag) = 'rmgag'
      ratnam(irsiga) = 'rsiga'
      ratnam(irmgap) = 'rmgap'
      ratnam(iralpa) = 'ralpa'
      ratnam(iralpg) = 'ralpg'
      ratnam(irsigp) = 'rsigp'
      ratnam(irsiag) = 'rsiag'
      ratnam(irsga)  = 'rsga '
      ratnam(irsiap) = 'rsiap'
      ratnam(irppa)  = 'rppa '
      ratnam(irppg)  = 'rppg '
      ratnam(irsgp)  = 'rsgp '
      ratnam(irsag)  = 'rsag '
      ratnam(irarga) = 'rarga'
      ratnam(irsap)  = 'rsap '
      ratnam(irclpa) = 'rclpa'
      ratnam(irclpg) = 'rclpg'
      ratnam(irargp) = 'rargp'
      ratnam(irarag) = 'rarag'
      ratnam(ircaga) = 'rcaga'
      ratnam(irarap) = 'rarap'
      ratnam(irkpa)  = 'rkpa '
      ratnam(irkpg)  = 'rkpg '
      ratnam(ircagp) = 'rcagp'
      ratnam(ircaag) = 'rcaag'
      ratnam(irtiga) = 'rtiga'
      ratnam(ircaap) = 'rcaap'
      ratnam(irscpa) = 'rscpa'
      ratnam(irscpg) = 'rscpg'
      ratnam(irtigp) = 'rtigp'
      ratnam(irtiag) = 'rtiag'
      ratnam(ircrga) = 'rcrga'
      ratnam(irtiap) = 'rtiap'
      ratnam(irvpa)  = 'rvpa '
      ratnam(irvpg)  = 'rvpg '
      ratnam(ircrgp) = 'rcrgp'
      ratnam(ircrag) = 'rcrag'
      ratnam(irfega) = 'rfega'
      ratnam(ircrap) = 'rcrap'
      ratnam(irmnpa) = 'rmnpa'
      ratnam(irmnpg) = 'rmnpg'
      ratnam(irfegp) = 'rfegp'
      ratnam(irfeag) = 'rfeag'
      ratnam(irniga) = 'rniga'
      ratnam(irfeap) = 'rfeap'
      ratnam(ircopa) = 'rcopa'
      ratnam(ircopg) = 'rcopg'
      ratnam(irnigp) = 'rnigp'

      ratnam(irr1)   = 'r1   '
      ratnam(irs1)   = 's1   '
      ratnam(irt1)   = 't1   '
      ratnam(iru1)   = 'u1   '
      ratnam(irv1)   = 'v1   '
      ratnam(irw1)   = 'w1   '
      ratnam(irx1)   = 'x1   '
      ratnam(iry1)   = 'y1   '


      return
      end
!---------------------------------------------------------------------




!---------------------------------------------------------------------
      subroutine net_input(tstart,tstep,tin,din,vin,zin,ein,xin)
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'
      include 'burn_common.dek'
      include 'network.dek'
      include 'cjdet.dek'


! declare the pass
      double precision tstart,tstep,tin,din,vin,zin,ein,xin(*)


! local variables
      character*80     string,word
      integer          i,j,k,ibtype,ictype,igues,kkase,ians,getnam
      double precision xneut,xh1,xhe4,xc12,xc13,xn14,xo16,xne20,xne22, &
                       xsi28,xfe52,xfe54,xfe56,xni56,zye,sum,abar,zbar, &
                       wbar,xcess,ye,ye_orig,xmup,xmun,qdum,a,z,xelem, &
                       andgrev,value


! bigbang specifics
      double precision fac,f1,zeta3
      parameter        (zeta3 = 1.20205690315732d0)


! popular format statements
01    format(1x,a,a,a)
02    format(1x,a,'=',1pe10.3,' ',a,'=',1pe10.3,' ', &
                a,'=',1pe10.3,' ',a,'=',1pe10.3,' ', &
                a,'=',1pe10.3)
 03   format(a)
 04   format(1x,a,'=',i2,' ',a,'=',i2,' ', &
                a,'=',i2,' ',a,'=',i2,' ', &
                a,'=',i2)



! initialize the common block variables
      call net_initialize


! inititailize local variables
      ibtype    = 0
      ictype    = 0
      tstart    = 0.0d0
      tstep     = 0.0d0
      bpres     = 0.0d0
      tin       = 0.0d0
      din       = 0.0d0
      vin       = 0.0d0
      zin       = 0.0d0
      zye       = 0.0d0
      xin(1:ionmax) = 1.0d-30


!---------------------------------------------------------------------------



! get the burn type
 10   write(6,01) 'give burning mode:'
      write(6,01) '     ibtype = 0 = stop'
      write(6,01) '              1 = onestep'
      write(6,01) '              2 = hydrostatic'
      write(6,01) '              3 = expansion'
      write(6,01) '              4 = self-heat at constant density'
      write(6,01) '              5 = self heat at constant pressure'
      write(6,01) '              6 = self-heat pressure-temp trajectory'
      write(6,01) '              7 = big bang '
      write(6,01) '              8 = detonation'
      write(6,01) '              9 = temp-den trajectory'

      read(5,*)  ibtype
      if (ibtype .lt. 0 .or. ibtype .gt. 9) goto 10

! set the burn type logical
      if (ibtype .eq. 0) then
       stop 'normal termination'
      else if (ibtype .eq. 1) then
       one_step = .true.
      else if (ibtype .eq. 2) then
       hydrostatic = .true.
      else if (ibtype .eq. 3) then
       expansion = .true.
      else if (ibtype .eq. 4) then
       self_heat_const_den = .true.
      else if (ibtype .eq. 5) then
       self_heat_const_pres = .true.
      else if (ibtype .eq. 6) then
       pt_hist = .true.
      else if (ibtype .eq. 7) then
       bbang = .true.
      else if (ibtype .eq. 8) then
       detonation = .true.
      else if (ibtype .eq. 9) then
       trho_hist = .true.
      else
       goto 10
      end if





! general options
 11   write(6,*)
      write(6,04) 'set general options:'
      write(6,04) 'screen_on',screen_on
      write(6,04) 'use_tables',use_tables
      write(6,04) 'weak_on',weak_on
      write(6,04) 'ffn_on',ffn_on
      write(6,04) 'pure_network',pure_network
      write(6,04) 'nse_analysis',nse_analysis
      write(6,04) 'allow_nse_evol',allow_nse_evol
      write(6,04) 'iprint_files',iprint_files
      write(6,04) 'iprint_screen',iprint_screen
      write(6,02) 'sthreshold',sthreshold,' set > 1 to disable'
      write(6,*)
      write(6,01) 'if these are ok, enter 1, otherwise enter 0 =>'

      read(5,*) ians
      if (ians .lt. 0 .or. ians .gt. 1) goto 11

      if (ians .eq. 0) then
 12    write(6,01) 'give the 9 integer and one real vector =>'

       read(5,*) screen_on, use_tables, weak_on, ffn_on, &
                 pure_network, nse_analysis, allow_nse_evol, &
                 iprint_files, iprint_screen, &
                 sthreshold

       if (screen_on .lt. 0 .or. screen_on .gt. 1) goto 12
       if (use_tables .lt. 0 .or. use_tables .gt. 1) goto 12
       if (weak_on .lt. 0 .or. weak_on .gt. 1) goto 12
       if (ffn_on .lt. 0 .or. ffn_on .gt. 1) goto 12
       if (pure_network .lt. 0 .or. pure_network .gt. 1) goto 12
       if (nse_analysis .lt. 0 .or. nse_analysis .gt. 1) goto 12
       if (iprint_files .lt. 0 .or. iprint_files .gt. 1) goto 12
       if (iprint_screen .lt. 0 .or. iprint_screen .gt. 1) goto 12
       goto 11
      end if



! get the bigbang parameters; set default to wmap 2008 (5 year) values
      if (bbang) then

       eta1    = 6.23e-10
       xnnu    = 3.0d0
       hubble  = 70.5d0
       cmbtemp = 2.725d0

 13    write(6,*)
       write(6,02) 'bigbang parameters:'
       write(6,02) 'eta',eta1
       write(6,02) 'number of neutrino families',xnnu
       write(6,02) 'hubble constant',hubble
       write(6,02) 'present cmb temperature',cmbtemp

       write(6,01) 'if these are ok, enter 1, otherwise enter 0 =>'
       read(5,*) ians
       if (ians .lt. 0 .or. ians .gt. 1) goto 13

       if (ians .eq. 0) then
        write(6,01) 'give eta, xnu, hubble, and cmbtemp  =>'
        read(5,*) eta1, xnnu, hubble, cmbtemp
        goto 13
       end if
      end if



! get an alternative the stopping condition; when the
! mass fraction of a given isotope falls below a given level


 14   write(6,*)
      write(6,*) 'stop when an isotope falls below a given abundance?', &
                  ' 1=yes 0=no'
      read(5,*)  ians
      if (ians .lt. 0 .or. ians .gt. 1) goto 14

      if (ians .eq. 0) then
       name_stop = 'he4 '
       xmass_stop = -1.0d30
      end if

 15   if (ians .eq. 1) then
       write(6,*) 'give the name of the isotope and the mass fraction'
       write(6,*) 'for example: c12 0.50'

       read(5,03) string
       j = 1
       i = getnam(string,word,j)
       name_stop = word(1:5)

       i = getnam(string,word,j)
       xmass_stop = value(word)

       write(6,*) name_stop,xmass_stop
      end if


! check that the name_stop isotope is in the network
      do i=1,ionmax
       if (ionam(i) .eq. name_stop) then
        id_stop = i
        goto 16
       end if
      enddo
      write(6,*)
      write(6,*) 'name_stop>',name_stop,'< not in network'
      write(6,*)
      if (ians .eq. 1) goto 15
      stop ' bad name for stopping isotope'
 16   continue




! get the initial thermodynamics
      write(6,*)
      if (self_heat_const_pres) then
       write(6,01) 'give the ending time, temperature, pressure =>'
       read(5,*)  tstep,tin,bpres

      else if (bbang) then
       write(6,01) 'give the ending time, initial temperature =>'
       read(5,*)  tstep,tin

      else if (.not. (trho_hist .or. pt_hist)) then
       write(6,01) 'give the ending time, temperature, density =>'
       read(5,*)  tstep,tin,din
      end if

! limit the temperature since the rates are invalid much above t9=100
       tin = min(1.0d11,tin)



! get the composition
      if (.not. bbang) then
 20    write(6,01) 'give initial composition:'
       write(6,01) '     ictype = 0 = leave alone; read from file'
       write(6,01) '              1 = solar abundances'
       write(6,01) '              2 = nse'
       write(6,01) '              3 = specify initial composition'

       read(5,*) ictype
       if (ictype .lt. 0 .or. ictype .gt. 3) goto 20

       if (ictype .eq. 3) then
        write(6,01) &
        'n h1 he4 c12 c13 n14 o16 ne20 ne22 si28 fe52 fe54 fe56 ni56 =>'
        read(5,*) xneut,xh1,xhe4,xc12,xc13,xn14,xo16,xne20,xne22,xsi28, &
                  xfe52,xfe54,xfe56,xni56
       end if
      end if


! get the output root file name
      write(6,*)  ' '
      write(6,01) 'give output root name, <cr> for default "foo_"=>'
      read(5,03) hfile
      if (hfile(1:2) .eq. '  ')  hfile = 'foo_'


!---------------------------------------------------------------------------



! set some more variables based on the burn type

! adiabatic expansion
! psi =  1 is an adiabatic expansion, -1 in an adiabatic implosion

      if (expansion) then
       psi       = 1.0d0
!       psi       = -1.0d0
       den0      = din
       temp0     = tin
       temp_stop = 1.0d7
!       temp_stop = 1.0d10
       if ( (psi .ge. 1.0  .and. temp_stop .ge. tin)  .or. &
            (psi .le. -1.0 .and. temp_stop .le. tin)) &
          stop 'bad adiabatic temp_stop in routine burner'



! big bang
      else if (bbang) then

! set the initial n and p abundances; equation 3 of wagoner et al 1967
       fac = exp((mn - mp)*clight**2/(kerg*tin))
       xneut = 1.0d0/(1.0d0 + fac)
       xh1   = 1.0d0 - xneut

! set the density from the temperature and eta1
       f1  = 30.0d0 * zeta3/pi**4 * asol/(kerg*avo)
       din = f1 * eta1 * tin**3


! thermodynamic profile being given
      else if (trho_hist .or. pt_hist) then
       write(6,*) 'give the trajectory file =>'
       read(5,03) trho_file
      end if


!---------------------------------------------------------------------------



! read the thermodynamic trajectory and initial abundances
! transfer the info stored in xsum and zsum from the update2 call

      if (trho_hist) then
       call update2(tstart,tin,din)
       xin(1:ionmax) = xsum(1:ionmax)
       tstart     = zwork1(1)
       tstep      = zwork1(2)
       zye        = zwork1(3)
      end if


      if (pt_hist) then
       call update3(tstart,tin,bpres)
       xin(1:ionmax) = xsum(1:ionmax)
       tstart     = zwork1(1)
       tstep      = zwork1(2)
       zye        = zwork1(3)
      end if



! massage the input composition, includes possible changes to the
! the abundances read in from the trho_hist file

! solar abundances
      if (ictype .eq. 1) then
       do i=1,ionmax
        xin(i) = andgrev(ionam(i),z,a,xelem)
       enddo
       if (iprot .ne. 0) xin(iprot) = andgrev('h1   ',z,a,xelem)


! put it in nse
      else if (ictype .eq. 2) then
       if (zye .eq. 0.0) zye   = 0.5d0
       igues = 1
       call nse(tin,din,zye,igues,1,1,xin,xmun,xmup,0)


! set the composition variables
      else if (ictype .eq. 3 .or. bbang) then
       if (ineut .ne. 0) xin(ineut) = xneut
       if (ih1   .ne. 0) xin(ih1)   = xh1
       if (iprot .ne. 0) xin(iprot) = xh1
       if (ih1 .ne. 0 .and. iprot .ne. 0) xin(iprot) = 0.0d0
       if (ihe4  .ne. 0) xin(ihe4)  = xhe4
       if (ic12  .ne. 0) xin(ic12)  = xc12
       if (ic13  .ne. 0) xin(ic13)  = xc13
       if (in14  .ne. 0) xin(in14)  = xn14
       if (io16  .ne. 0) xin(io16)  = xo16
       if (ine20 .ne. 0) xin(ine20) = xne20
       if (ine22 .ne. 0) xin(ine22) = xne22
       if (isi28 .ne. 0) xin(isi28) = xsi28
       if (ife52 .ne. 0) xin(ife52) = xfe52
       if (ife54 .ne. 0) xin(ife54) = xfe54
       if (ife56 .ne. 0) xin(ife56) = xfe56
       if (ini56 .ne. 0) xin(ini56) = xni56


! hardcode something here

!if (ih1 .ne. 0)   xin(ih1)=       7.0572558936810803E-01
!if (ih2 .ne. 0)   xin(ih2)=       4.8010000000000003E-05
!if (ihe3 .ne. 0)  xin(ihe3)=      2.9291000000000001E-05
!if (ihe4 .ne. 0)  xin(ihe4)=      2.7521000000000001E-01
!if (ili7 .ne. 0)  xin(ili7)=      9.3489999999999999E-09
!if (ic12 .ne. 0)  xin(ic12)=      3.0324000000000002E-03
!if (ic13 .ne. 0)  xin(ic13)=      3.6501000000000002E-05
!if (in14 .ne. 0)  xin(in14)=      1.1049000000000000E-03
!if (in15 .ne. 0)  xin(in15)=      4.3633999999999996E-06
!if (io16 .ne. 0)  xin(io16)=      9.5917999999999993E-03
!if (io17 .ne. 0)  xin(io17)=      3.8873000000000000E-06
!if (io18 .ne. 0)  xin(io18)=      2.1673000000000001E-05
!if (if19 .ne. 0)  xin(if19)=      4.0515000000000000E-07
!if (ine20 .ne. 0) xin(ine20)=     1.6188999999999999E-03
!if (img24 .ne. 0) xin(img24)=     3.5722704328918775E-03


!if (ihe4 .ne. 0)  xin(ihe4)=     5.45516e-07 
!if (ic12 .ne. 0)  xin(ic12)=     0.491254d0
!if (io16 .ne. 0)  xin(io16)=      0.494312d0
!if (ine20 .ne. 0) xin(ine20)=     0.0142452d0
!if (img24 .ne. 0) xin(img24)=     0.000187270d0
!if (isi28 .ne. 0) xin(isi28)=     9.08493e-07
!if (is32 .ne. 0) xin(is32)=      4.30614e-11
!if (iar36 .ne. 0) xin(iar36)=     5.21367e-16
!if (ica40 .ne. 0) xin(ica40)=     1.06910e-20
!if (iti44 .ne. 0) xin(iti44)=     1.00000e-20
!if (icr48 .ne. 0) xin(icr48)=     1.00000e-20
!if (ife52 .ne. 0) xin(ife52)=     1.00000e-20
!if (ini56 .ne. 0) xin(ini56)=     1.00000e-20


      end if


! write out the input composition so far
!      write(6,02) (ionam(i),xin(i), i=1,ionmax)
!      read(5,*)


! normalize the composition
      do i=1,ionmax
       xin(i) = min(1.0d0,max(xin(i),1.0d-30))
      end do
      sum = 0.0d0
       do i=1,ionmax
        sum = sum + xin(i)
       enddo
      sum = 1.0d0/sum
      do i=1,ionmax
       xin(i) = min(1.0d0,max(xin(i) * sum,1.0d-30))
      enddo

!      write(6,*) 'post norm', ionmax,xin(ih1)
!      write(6,02) (ionam(i),xin(i), i=1,ionmax)
!      read(5,*)

!---------------------------------------------------------------------------


! get the ye of the initial compositon
        call azbar(xin,aion,zion,wion,ionmax, &
                   zwork1,abar,zbar,wbar,ye_orig,xcess)

!       write(6,123) abar,zbar
!123    format(1x,1p2e12.3)
!       read(5,*)



! modify the composition if ye_orig is less than 0.55
!        if (ye_orig .le. 0.55) then
!
! set the mass fraction of fe58 to set the desired ye
!         ye_want = 0.495d0
!         ye_want = 0.50d0
!         if (ye_want .eq. 0.5) then
!          xin(ife58) = 0.0d0
!         else
!          xin(ife58) = (ye_orig - ye_want) /
!     1                  (ye_orig - zion(ife58)/aion(ife58))
!         end if
!
! reset the mass fractions of everything else
!         sum = 1.0d0 - xin(ife58)
!         do i=1,ionmax
!          if (i .ne. ife58) xin(i) = xin(i) * sum
!         enddo
!        end if


!---------------------------------------------------------------------------


! modify for a detonation
! get the chapman-jouget solution
       if (detonation) then
        kkase = 1
         mach  = 0.0d0
        do i=1,ionmax
         xmass_up(i) = xin(i)
        enddo
        temp_up = tin
        den_up  = din
        call cjsolve(kkase,xmass_up,temp_up,den_up,mach, &
                    qburn_cj,xmass_cj,ener_up,pres_up,cs_up, &
                    vel_det,vel_cj,temp_cj,den_cj,ener_cj,pres_cj,cs_cj)


        write(6,*)  ' '
        write(6,63) 'cj state (should be sonic with vel_mat = cs_cj):'
        write(6,61) 'temp_cj',temp_cj,'den_cj ',den_cj, &
                    'pres_cj',pres_cj
        write(6,61) 'cs_cj  ',cs_cj, &
                    'vel_mat',vel_cj,'vel_det',vel_det
        write(6,61) 'mach_cj',vel_cj/cs_cj,'qburn_cj',qburn_cj

 63     format(1x,a)
 61     format(1x,a7,'=',1pe10.3,' ',a7,'=',1pe10.3,' ', &
                 a7,'=',1pe10.3,' ',a4,'=',1pe10.3)


        write(6,*) ' '
        write(6,*) 'top 10 cj nse mass fractions:'
        call indexx(ionmax,xmass_cj,izwork1)
        write(6,02) (ionam(izwork1(i)), &
                   xmass_cj(izwork1(i)), i=ionmax,ionmax-9,-1)


! get shock solution
        kkase = 4
        mach_sh = vel_det/cs_up
        call cjsolve(kkase,xmass_up,temp_up,den_up,mach_sh, &
                    qdum,xmass_up,ener_up,pres_up,cs_up, &
                    vel_det,vel_sh,temp_sh,den_sh,ener_sh,pres_sh,cs_sh)


! reset the initial conditions for znd detonations
        tin      = temp_sh
        din      = den_sh
        vin      = vel_sh
        zin      = 1.0e-16*vel_sh
        den_stop = 1.00d0 * den_cj

        write(6,*)
        write(6,*) 'resetting initial conditions for a detonation to:'
        write(6,64) 'tin=',tin,' din=',din,' vin=',vin,' zin=',zin
 64     format(1x,4(a,1pe12.4) )
       end if


!---------------------------------------------------------------------------


! get the abundance variables for the final mixture
        call azbar(xin,aion,zion,wion,ionmax, &
                   zwork1,abar,zbar,wbar,ye_orig,xcess)


! get the thermodynamic state
      temp_row(1) = tin
      den_row(1)  = din
      ptot_row(1) = bpres
      abar_row(1) = abar
      zbar_row(1) = zbar
      jlo_eos = 1
      jhi_eos = 1

!      write(6,*) tin,abar,zbar

      if (self_heat_const_pres .or. pt_hist) then
       den_row(1)  = bpres * abar/(avo * kerg * tin)
       call invert_helm_pt
       din = den_row(1)

!       write(6,778) bpres,din
!       read(5,*)

      else
       call helmeos
       bpres = ptot_row(1)
      endif

      ein   = etot_row(1)


!---------------------------------------------------------------------------


! write out the final input
        write(6,*)
        write(6,02) 'tstart',tstart,'tstep',tstep
        write(6,02) 'tin',tin,'din',din,'bpres',bpres,'ein',ein

! largest mass fractions
        call indexx(ionmax,xin,izwork1)
        j = min(20,ionmax)
        k = max(ionmax-19,1)
        write(6,*) j,' largest mass fractions'
        do i=ionmax,k,-1
         if (xin(izwork1(i)) .gt. 1.0e-12) &
            write(6,02) ionam(izwork1(i)),xin(izwork1(i))
        end do

! nonconservation, abar, zbar of the mixture
        sum = 0.0d0
         do i=1,ionmax
          sum = sum + xin(i)
         enddo
        write(6,02) '1-sum',1.0d0 - sum
        write(6,02) 'abar',abar,'zbar',zbar,'ye',zbar/abar
        write(6,*)

!        read(5,*)



! there is probably a better place for this
! if requested, adjust the number of equations being solved
      if (pure_network .eq. 1) then
       neqs  = ionmax
       btemp = tin
       bden  = din
      end if

      return
      end
!---------------------------------------------------------------------







!---------------------------------------------------------------------
!
! this routine contains auxillary network routine

! routines for a tree construction to mark nonzero matrix locations
! routine screen6 computes screening factors
! routine screen5 computes screening factors
! routine snupp computes neutrino loss rates for the pp chain
! routine snucno computes neutrino loss rates for the cno cycles
! routine sneut5 computes neutrino loss rates
! routine ifermi12 does an inverse fermi integral of order 1/2
! routine zfermim12 does an inverse fermi integral of order -1/2

! routine ecapnuc02 computes electron capture rates
! routine ecapnuc computes electron capture rates
! routine mazurek computes ni56 electron capture rates
! routine time_scales computes various timescales
! routine ener_gener_rate computes the instantaneous energy generation rate





      subroutine sneut5(temp,den,abar,zbar, &
                        snu,dsnudt,dsnudd,dsnuda,dsnudz)
      include 'implno.dek'
      include 'const.dek'

! this routine computes neutrino losses from the analytic fits of
! itoh et al. apjs 102, 411, 1996, and also returns their derivatives.

! input:
! temp = temperature
! den  = density
! abar = mean atomic weight
! zbar = mean charge

! output:
! snu    = total neutrino loss rate in erg/g/sec
! dsnudt = derivative of snu with temperature
! dsnudd = derivative of snu with density
! dsnuda = derivative of snu with abar
! dsnudz = derivative of snu with zbar


! declare the pass
      double precision temp,den,abar,zbar, &
                       snu,dsnudt,dsnudd,dsnuda,dsnudz

! local variables
      double precision spair,spairdt,spairdd,spairda,spairdz, &
                       splas,splasdt,splasdd,splasda,splasdz, &
                       sphot,sphotdt,sphotdd,sphotda,sphotdz, &
                       sbrem,sbremdt,sbremdd,sbremda,sbremdz, &
                       sreco,srecodt,srecodd,srecoda,srecodz

      double precision t9,xl,xldt,xlp5,xl2,xl3,xl4,xl5,xl6,xl7,xl8,xl9, &
                       xlmp5,xlm1,xlm2,xlm3,xlm4,xlnt,cc,den6,tfermi, &
                       a0,a1,a2,a3,b1,b2,c00,c01,c02,c03,c04,c05,c06, &
                       c10,c11,c12,c13,c14,c15,c16,c20,c21,c22,c23,c24, &
                       c25,c26,dd00,dd01,dd02,dd03,dd04,dd05,dd11,dd12, &
                       dd13,dd14,dd15,dd21,dd22,dd23,dd24,dd25,b,c,d,f0, &
                       f1,deni,tempi,abari,zbari,f2,f3,z,xmue,ye, &
                       dum,dumdt,dumdd,dumda,dumdz, &
                       gum,gumdt,gumdd,gumda,gumdz


! pair production
      double precision rm,rmdd,rmda,rmdz,rmi,gl,gldt, &
                       zeta,zetadt,zetadd,zetada,zetadz,zeta2,zeta3, &
                       xnum,xnumdt,xnumdd,xnumda,xnumdz, &
                       xden,xdendt,xdendd,xdenda,xdendz, &
                       fpair,fpairdt,fpairdd,fpairda,fpairdz, &
                       qpair,qpairdt,qpairdd,qpairda,qpairdz

! plasma
      double precision gl2,gl2dt,gl2dd,gl2da,gl2dz,gl12,gl32,gl72,gl6, &
                       ft,ftdt,ftdd,ftda,ftdz,fl,fldt,fldd,flda,fldz, &
                       fxy,fxydt,fxydd,fxyda,fxydz

! photo
      double precision tau,taudt,cos1,cos2,cos3,cos4,cos5,sin1,sin2, &
                       sin3,sin4,sin5,last,xast, &
                       fphot,fphotdt,fphotdd,fphotda,fphotdz, &
                       qphot,qphotdt,qphotdd,qphotda,qphotdz

! brem
      double precision t8,t812,t832,t82,t83,t85,t86,t8m1,t8m2,t8m3,t8m5, &
                       t8m6, &
                       eta,etadt,etadd,etada,etadz,etam1,etam2,etam3, &
                       fbrem,fbremdt,fbremdd,fbremda,fbremdz, &
                       gbrem,gbremdt,gbremdd,gbremda,gbremdz, &
                       u,gm1,gm2,gm13,gm23,gm43,gm53,v,w,fb,gt,gb, &
                       fliq,fliqdt,fliqdd,fliqda,fliqdz, &
                       gliq,gliqdt,gliqdd,gliqda,gliqdz

! recomb
      double precision ifermi12,zfermim12,nu,nudt,nudd,nuda,nudz, &
                       nu2,nu3,bigj,bigjdt,bigjdd,bigjda,bigjdz



! numerical constants
      double precision fac1,fac2,fac3,oneth,twoth,con1,sixth,iln10
      parameter        (fac1   = 5.0d0 * pi / 3.0d0, &
                        fac2   = 10.0d0 * pi, &
                        fac3   = pi / 5.0d0, &
                        oneth  = 1.0d0/3.0d0, &
                        twoth  = 2.0d0/3.0d0, &
                        con1   = 1.0d0/5.9302d0, &
                        sixth  = 1.0d0/6.0d0, &
                        iln10  = 4.342944819032518d-1)


! theta is sin**2(theta_weinberg) = 0.2319 plus/minus 0.00005 (1996)
! xnufam is the number of neutrino flavors = 3.02 plus/minus 0.005 (1998)
! change theta and xnufam if need be, and the changes will automatically
! propagate through the routine. cv and ca are the vector and axial currents.

      double precision theta,xnufam,cv,ca,cvp,cap,tfac1,tfac2,tfac3, &
                       tfac4,tfac5,tfac6
      parameter        (theta  = 0.2319d0, &
                        xnufam = 3.0d0, &
                        cv     = 0.5d0 + 2.0d0 * theta, &
                        cvp    = 1.0d0 - cv, &
                        ca     = 0.5d0, &
                        cap    = 1.0d0 - ca, &
                        tfac1  = cv*cv + ca*ca + &
                                 (xnufam-1.0d0) * (cvp*cvp+cap*cap), &
                        tfac2  = cv*cv - ca*ca + &
                                 (xnufam-1.0d0) * (cvp*cvp - cap*cap), &
                        tfac3  = tfac2/tfac1, &
                        tfac4  = 0.5d0 * tfac1, &
                        tfac5  = 0.5d0 * tfac2, &
                        tfac6  = cv*cv + 1.5d0*ca*ca + (xnufam - 1.0d0)* &
                                 (cvp*cvp + 1.5d0*cap*cap))



! initialize
      spair   = 0.0d0
      spairdt = 0.0d0
      spairdd = 0.0d0
      spairda = 0.0d0
      spairdz = 0.0d0

      splas   = 0.0d0
      splasdt = 0.0d0
      splasdd = 0.0d0
      splasda = 0.0d0
      splasdz = 0.0d0

      sphot   = 0.0d0
      sphotdt = 0.0d0
      sphotdd = 0.0d0
      sphotda = 0.0d0
      sphotdz = 0.0d0

      sbrem   = 0.0d0
      sbremdt = 0.0d0
      sbremdd = 0.0d0
      sbremda = 0.0d0
      sbremdz = 0.0d0

      sreco   = 0.0d0
      srecodt = 0.0d0
      srecodd = 0.0d0
      srecoda = 0.0d0
      srecodz = 0.0d0

      snu     = 0.0d0
      dsnudt  = 0.0d0
      dsnudd  = 0.0d0
      dsnuda  = 0.0d0
      dsnudz  = 0.0d0

      if (temp .lt. 1.0e7) return


! to avoid lots of divisions
      deni  = 1.0d0/den
      tempi = 1.0d0/temp
      abari = 1.0d0/abar
      zbari = 1.0d0/zbar


! some composition variables
      ye    = zbar*abari
      xmue  = abar*zbari




! some frequent factors
      t9     = temp * 1.0d-9
      xl     = t9 * con1
      xldt   = 1.0d-9 * con1
      xlp5   = sqrt(xl)
      xl2    = xl*xl
      xl3    = xl2*xl
      xl4    = xl3*xl
      xl5    = xl4*xl
      xl6    = xl5*xl
      xl7    = xl6*xl
      xl8    = xl7*xl
      xl9    = xl8*xl
      xlmp5  = 1.0d0/xlp5
      xlm1   = 1.0d0/xl
      xlm2   = xlm1*xlm1
      xlm3   = xlm1*xlm2
      xlm4   = xlm1*xlm3

      rm     = den*ye
      rmdd   = ye
      rmda   = -rm*abari
      rmdz   = den*abari
      rmi    = 1.0d0/rm

      a0     = rm * 1.0d-9
      a1     = a0**oneth
      zeta   = a1 * xlm1
      zetadt = -a1 * xlm2 * xldt
      a2     = oneth * a1*rmi * xlm1
      zetadd = a2 * rmdd
      zetada = a2 * rmda
      zetadz = a2 * rmdz

      zeta2 = zeta * zeta
      zeta3 = zeta2 * zeta




! pair neutrino section
! for reactions like e+ + e- => nu_e + nubar_e

! equation 2.8
      gl   = 1.0d0 - 13.04d0*xl2 +133.5d0*xl4 +1534.0d0*xl6 +918.6d0*xl8
      gldt = xldt*(-26.08d0*xl +534.0d0*xl3 +9204.0d0*xl5 +7348.8d0*xl7)

! equation 2.7

      a1     = 6.002d19 + 2.084d20*zeta + 1.872d21*zeta2
      a2     = 2.084d20 + 2.0d0*1.872d21*zeta

      if (t9 .lt. 10.0) then
       b1     = exp(-5.5924d0*zeta)
       b2     = -b1*5.5924d0
      else
       b1     = exp(-4.9924d0*zeta)
       b2     = -b1*4.9924d0
      end if

      xnum   = a1 * b1
      c      = a2*b1 + a1*b2
      xnumdt = c*zetadt
      xnumdd = c*zetadd
      xnumda = c*zetada
      xnumdz = c*zetadz

      if (t9 .lt. 10.0) then
       a1   = 9.383d-1*xlm1 - 4.141d-1*xlm2 + 5.829d-2*xlm3
       a2   = -9.383d-1*xlm2 + 2.0d0*4.141d-1*xlm3 - 3.0d0*5.829d-2*xlm4
      else
       a1   = 1.2383d0*xlm1 - 8.141d-1*xlm2
       a2   = -1.2383d0*xlm2 + 2.0d0*8.141d-1*xlm3
      end if

      b1   = 3.0d0*zeta2

      xden   = zeta3 + a1
      xdendt = b1*zetadt + a2*xldt
      xdendd = b1*zetadd
      xdenda = b1*zetada
      xdendz = b1*zetadz

      a1      = 1.0d0/xden
      fpair   = xnum*a1
      fpairdt = (xnumdt - fpair*xdendt)*a1
      fpairdd = (xnumdd - fpair*xdendd)*a1
      fpairda = (xnumda - fpair*xdenda)*a1
      fpairdz = (xnumdz - fpair*xdendz)*a1


! equation 2.6
      a1     = 10.7480d0*xl2 + 0.3967d0*xlp5 + 1.005d0
      a2     = xldt*(2.0d0*10.7480d0*xl + 0.5d0*0.3967d0*xlmp5)
      xnum   = 1.0d0/a1
      xnumdt = -xnum*xnum*a2

      a1     = 7.692d7*xl3 + 9.715d6*xlp5
      a2     = xldt*(3.0d0*7.692d7*xl2 + 0.5d0*9.715d6*xlmp5)

      c      = 1.0d0/a1
      b1     = 1.0d0 + rm*c

      xden   = b1**(-0.3d0)

      d      = -0.3d0*xden/b1
      xdendt = -d*rm*c*c*a2
      xdendd = d*rmdd*c
      xdenda = d*rmda*c
      xdendz = d*rmdz*c

      qpair   = xnum*xden
      qpairdt = xnumdt*xden + xnum*xdendt
      qpairdd = xnum*xdendd
      qpairda = xnum*xdenda
      qpairdz = xnum*xdendz



! equation 2.5
      a1    = exp(-2.0d0*xlm1)
      a2    = a1*2.0d0*xlm2*xldt

      spair   = a1*fpair
      spairdt = a2*fpair + a1*fpairdt
      spairdd = a1*fpairdd
      spairda = a1*fpairda
      spairdz = a1*fpairdz

      a1      = spair
      spair   = gl*a1
      spairdt = gl*spairdt + gldt*a1
      spairdd = gl*spairdd
      spairda = gl*spairda
      spairdz = gl*spairdz

      a1      = tfac4*(1.0d0 + tfac3 * qpair)
      a2      = tfac4*tfac3

      a3      = spair
      spair   = a1*a3
      spairdt = a1*spairdt + a2*qpairdt*a3
      spairdd = a1*spairdd + a2*qpairdd*a3
      spairda = a1*spairda + a2*qpairda*a3
      spairdz = a1*spairdz + a2*qpairdz*a3




! plasma neutrino section
! for collective reactions like gamma_plasmon => nu_e + nubar_e
! equation 4.6

      a1   = 1.019d-6*rm
      a2   = a1**twoth
      a3   = twoth*a2/a1

      b1   =  sqrt(1.0d0 + a2)
      b2   = 1.0d0/b1

      c00  = 1.0d0/(temp*temp*b1)

      gl2   = 1.1095d11 * rm * c00

      gl2dt = -2.0d0*gl2*tempi
      d     = rm*c00*b2*0.5d0*b2*a3*1.019d-6
      gl2dd = 1.1095d11 * (rmdd*c00  - d*rmdd)
      gl2da = 1.1095d11 * (rmda*c00  - d*rmda)
      gl2dz = 1.1095d11 * (rmdz*c00  - d*rmdz)


      gl    = sqrt(gl2)
      gl12  = sqrt(gl)
      gl32  = gl * gl12
      gl72  = gl2 * gl32
      gl6   = gl2 * gl2 * gl2


! equation 4.7
      ft   = 2.4d0 + 0.6d0*gl12 + 0.51d0*gl + 1.25d0*gl32
      gum  = 1.0d0/gl2
      a1   =(0.25d0*0.6d0*gl12 +0.5d0*0.51d0*gl +0.75d0*1.25d0*gl32)*gum
      ftdt = a1*gl2dt
      ftdd = a1*gl2dd
      ftda = a1*gl2da
      ftdz = a1*gl2dz


! equation 4.8
      a1   = 8.6d0*gl2 + 1.35d0*gl72
      a2   = 8.6d0 + 1.75d0*1.35d0*gl72*gum

      b1   = 225.0d0 - 17.0d0*gl + gl2
      b2   = -0.5d0*17.0d0*gl*gum + 1.0d0

      c    = 1.0d0/b1
      fl   = a1*c

      d    = (a2 - fl*b2)*c
      fldt = d*gl2dt
      fldd = d*gl2dd
      flda = d*gl2da
      fldz = d*gl2dz


! equation 4.9 and 4.10
      cc   = log10(2.0d0*rm)
      xlnt = log10(temp)

      xnum   = sixth * (17.5d0 + cc - 3.0d0*xlnt)
      xnumdt = -iln10*0.5d0*tempi
      a2     = iln10*sixth*rmi
      xnumdd = a2*rmdd
      xnumda = a2*rmda
      xnumdz = a2*rmdz

      xden   = sixth * (-24.5d0 + cc + 3.0d0*xlnt)
      xdendt = iln10*0.5d0*tempi
      xdendd = a2*rmdd
      xdenda = a2*rmda
      xdendz = a2*rmdz


! equation 4.11
      if (abs(xnum) .gt. 0.7d0  .or.  xden .lt. 0.0d0) then
       fxy   = 1.0d0
       fxydt = 0.0d0
       fxydd = 0.0d0
       fxydz = 0.0d0
       fxyda = 0.0d0

      else

       a1  = 0.39d0 - 1.25d0*xnum - 0.35d0*sin(4.5d0*xnum)
       a2  = -1.25d0 - 4.5d0*0.35d0*cos(4.5d0*xnum)

       b1  = 0.3d0 * exp(-1.0d0*(4.5d0*xnum + 0.9d0)**2)
       b2  = -b1*2.0d0*(4.5d0*xnum + 0.9d0)*4.5d0

       c   = min(0.0d0, xden - 1.6d0 + 1.25d0*xnum)
       if (c .eq. 0.0) then
        dumdt = 0.0d0
        dumdd = 0.0d0
        dumda = 0.0d0
        dumdz = 0.0d0
       else
        dumdt = xdendt + 1.25d0*xnumdt
        dumdd = xdendd + 1.25d0*xnumdd
        dumda = xdenda + 1.25d0*xnumda
        dumdz = xdendz + 1.25d0*xnumdz
       end if

       d   = 0.57d0 - 0.25d0*xnum
       a3  = c/d
       c00 = exp(-1.0d0*a3**2)

       f1  = -c00*2.0d0*a3/d
       c01 = f1*(dumdt + a3*0.25d0*xnumdt)
       c02 = f1*(dumdd + a3*0.25d0*xnumdd)
       c03 = f1*(dumda + a3*0.25d0*xnumda)
       c04 = f1*(dumdz + a3*0.25d0*xnumdz)

       fxy   = 1.05d0 + (a1 - b1)*c00
       fxydt = (a2*xnumdt -  b2*xnumdt)*c00 + (a1-b1)*c01
       fxydd = (a2*xnumdd -  b2*xnumdd)*c00 + (a1-b1)*c02
       fxyda = (a2*xnumda -  b2*xnumda)*c00 + (a1-b1)*c03
       fxydz = (a2*xnumdz -  b2*xnumdz)*c00 + (a1-b1)*c04

      end if



! equation 4.1 and 4.5
      splas   = (ft + fl) * fxy
      splasdt = (ftdt + fldt)*fxy + (ft+fl)*fxydt
      splasdd = (ftdd + fldd)*fxy + (ft+fl)*fxydd
      splasda = (ftda + flda)*fxy + (ft+fl)*fxyda
      splasdz = (ftdz + fldz)*fxy + (ft+fl)*fxydz

      a2      = exp(-gl)
      a3      = -0.5d0*a2*gl*gum

      a1      = splas
      splas   = a2*a1
      splasdt = a2*splasdt + a3*gl2dt*a1
      splasdd = a2*splasdd + a3*gl2dd*a1
      splasda = a2*splasda + a3*gl2da*a1
      splasdz = a2*splasdz + a3*gl2dz*a1

      a2      = gl6
      a3      = 3.0d0*gl6*gum

      a1      = splas
      splas   = a2*a1
      splasdt = a2*splasdt + a3*gl2dt*a1
      splasdd = a2*splasdd + a3*gl2dd*a1
      splasda = a2*splasda + a3*gl2da*a1
      splasdz = a2*splasdz + a3*gl2dz*a1


      a2      = 0.93153d0 * 3.0d21 * xl9
      a3      = 0.93153d0 * 3.0d21 * 9.0d0*xl8*xldt

      a1      = splas
      splas   = a2*a1
      splasdt = a2*splasdt + a3*a1
      splasdd = a2*splasdd
      splasda = a2*splasda
      splasdz = a2*splasdz




! photoneutrino process section
! for reactions like e- + gamma => e- + nu_e + nubar_e
!                    e+ + gamma => e+ + nu_e + nubar_e
! equation 3.8 for tau, equation 3.6 for cc,
! and table 2 written out for speed
      if (temp .ge. 1.0d7  .and. temp .lt. 1.0d8) then
       tau  =  log10(temp * 1.0d-7)
       cc   =  0.5654d0 + tau
       c00  =  1.008d11
       c01  =  0.0d0
       c02  =  0.0d0
       c03  =  0.0d0
       c04  =  0.0d0
       c05  =  0.0d0
       c06  =  0.0d0
       c10  =  8.156d10
       c11  =  9.728d8
       c12  = -3.806d9
       c13  = -4.384d9
       c14  = -5.774d9
       c15  = -5.249d9
       c16  = -5.153d9
       c20  =  1.067d11
       c21  = -9.782d9
       c22  = -7.193d9
       c23  = -6.936d9
       c24  = -6.893d9
       c25  = -7.041d9
       c26  = -7.193d9
       dd01 =  0.0d0
       dd02 =  0.0d0
       dd03 =  0.0d0
       dd04 =  0.0d0
       dd05 =  0.0d0
       dd11 = -1.879d10
       dd12 = -9.667d9
       dd13 = -5.602d9
       dd14 = -3.370d9
       dd15 = -1.825d9
       dd21 = -2.919d10
       dd22 = -1.185d10
       dd23 = -7.270d9
       dd24 = -4.222d9
       dd25 = -1.560d9

      else if (temp .ge. 1.0d8  .and. temp .lt. 1.0d9) then
       tau   =  log10(temp * 1.0d-8)
       cc   =  1.5654d0
       c00  =  9.889d10
       c01  = -4.524d8
       c02  = -6.088d6
       c03  =  4.269d7
       c04  =  5.172d7
       c05  =  4.910d7
       c06  =  4.388d7
       c10  =  1.813d11
       c11  = -7.556d9
       c12  = -3.304d9
       c13  = -1.031d9
       c14  = -1.764d9
       c15  = -1.851d9
       c16  = -1.928d9
       c20  =  9.750d10
       c21  =  3.484d10
       c22  =  5.199d9
       c23  = -1.695d9
       c24  = -2.865d9
       c25  = -3.395d9
       c26  = -3.418d9
       dd01 = -1.135d8
       dd02 =  1.256d8
       dd03 =  5.149d7
       dd04 =  3.436d7
       dd05 =  1.005d7
       dd11 =  1.652d9
       dd12 = -3.119d9
       dd13 = -1.839d9
       dd14 = -1.458d9
       dd15 = -8.956d8
       dd21 = -1.549d10
       dd22 = -9.338d9
       dd23 = -5.899d9
       dd24 = -3.035d9
       dd25 = -1.598d9

      else if (temp .ge. 1.0d9) then
       tau  =  log10(t9)
       cc   =  1.5654d0
       c00  =  9.581d10
       c01  =  4.107d8
       c02  =  2.305d8
       c03  =  2.236d8
       c04  =  1.580d8
       c05  =  2.165d8
       c06  =  1.721d8
       c10  =  1.459d12
       c11  =  1.314d11
       c12  = -1.169d11
       c13  = -1.765d11
       c14  = -1.867d11
       c15  = -1.983d11
       c16  = -1.896d11
       c20  =  2.424d11
       c21  = -3.669d9
       c22  = -8.691d9
       c23  = -7.967d9
       c24  = -7.932d9
       c25  = -7.987d9
       c26  = -8.333d9
       dd01 =  4.724d8
       dd02 =  2.976d8
       dd03 =  2.242d8
       dd04 =  7.937d7
       dd05 =  4.859d7
       dd11 = -7.094d11
       dd12 = -3.697d11
       dd13 = -2.189d11
       dd14 = -1.273d11
       dd15 = -5.705d10
       dd21 = -2.254d10
       dd22 = -1.551d10
       dd23 = -7.793d9
       dd24 = -4.489d9
       dd25 = -2.185d9
      end if

      taudt = iln10*tempi


! equation 3.7, compute the expensive trig functions only one time
      cos1 = cos(fac1*tau)
      cos2 = cos(fac1*2.0d0*tau)
      cos3 = cos(fac1*3.0d0*tau)
      cos4 = cos(fac1*4.0d0*tau)
      cos5 = cos(fac1*5.0d0*tau)
      last = cos(fac2*tau)

      sin1 = sin(fac1*tau)
      sin2 = sin(fac1*2.0d0*tau)
      sin3 = sin(fac1*3.0d0*tau)
      sin4 = sin(fac1*4.0d0*tau)
      sin5 = sin(fac1*5.0d0*tau)
      xast = sin(fac2*tau)

      a0 = 0.5d0*c00 &
           + c01*cos1 + dd01*sin1 + c02*cos2 + dd02*sin2 &
           + c03*cos3 + dd03*sin3 + c04*cos4 + dd04*sin4 &
           + c05*cos5 + dd05*sin5 + 0.5d0*c06*last

      f0 =  taudt*fac1*(-c01*sin1 + dd01*cos1 - c02*sin2*2.0d0 &
           + dd02*cos2*2.0d0 - c03*sin3*3.0d0 + dd03*cos3*3.0d0 &
           - c04*sin4*4.0d0 + dd04*cos4*4.0d0 &
           - c05*sin5*5.0d0 + dd05*cos5*5.0d0) &
           - 0.5d0*c06*xast*fac2*taudt

      a1 = 0.5d0*c10 &
           + c11*cos1 + dd11*sin1 + c12*cos2 + dd12*sin2 &
           + c13*cos3 + dd13*sin3 + c14*cos4 + dd14*sin4 &
           + c15*cos5 + dd15*sin5 + 0.5d0*c16*last

      f1 = taudt*fac1*(-c11*sin1 + dd11*cos1 - c12*sin2*2.0d0 &
           + dd12*cos2*2.0d0 - c13*sin3*3.0d0 + dd13*cos3*3.0d0 &
           - c14*sin4*4.0d0 + dd14*cos4*4.0d0 - c15*sin5*5.0d0 &
           + dd15*cos5*5.0d0) - 0.5d0*c16*xast*fac2*taudt

      a2 = 0.5d0*c20 &
           + c21*cos1 + dd21*sin1 + c22*cos2 + dd22*sin2 &
           + c23*cos3 + dd23*sin3 + c24*cos4 + dd24*sin4 &
           + c25*cos5 + dd25*sin5 + 0.5d0*c26*last

      f2 = taudt*fac1*(-c21*sin1 + dd21*cos1 - c22*sin2*2.0d0 &
           + dd22*cos2*2.0d0 - c23*sin3*3.0d0 + dd23*cos3*3.0d0 &
           - c24*sin4*4.0d0 + dd24*cos4*4.0d0 - c25*sin5*5.0d0 &
           + dd25*cos5*5.0d0) - 0.5d0*c26*xast*fac2*taudt

! equation 3.4
      dum   = a0 + a1*zeta + a2*zeta2
      dumdt = f0 + f1*zeta + a1*zetadt + f2*zeta2 + 2.0d0*a2*zeta*zetadt
      dumdd = a1*zetadd + 2.0d0*a2*zeta*zetadd
      dumda = a1*zetada + 2.0d0*a2*zeta*zetada
      dumdz = a1*zetadz + 2.0d0*a2*zeta*zetadz

      z      = exp(-cc*zeta)

      xnum   = dum*z
      xnumdt = dumdt*z - dum*z*cc*zetadt
      xnumdd = dumdd*z - dum*z*cc*zetadd
      xnumda = dumda*z - dum*z*cc*zetada
      xnumdz = dumdz*z - dum*z*cc*zetadz

      xden   = zeta3 + 6.290d-3*xlm1 + 7.483d-3*xlm2 + 3.061d-4*xlm3

      dum    = 3.0d0*zeta2
      xdendt = dum*zetadt - xldt*(6.290d-3*xlm2 &
               + 2.0d0*7.483d-3*xlm3 + 3.0d0*3.061d-4*xlm4)
      xdendd = dum*zetadd
      xdenda = dum*zetada
      xdendz = dum*zetadz

      dum      = 1.0d0/xden
      fphot   = xnum*dum
      fphotdt = (xnumdt - fphot*xdendt)*dum
      fphotdd = (xnumdd - fphot*xdendd)*dum
      fphotda = (xnumda - fphot*xdenda)*dum
      fphotdz = (xnumdz - fphot*xdendz)*dum


! equation 3.3
      a0     = 1.0d0 + 2.045d0 * xl
      xnum   = 0.666d0*a0**(-2.066d0)
      xnumdt = -2.066d0*xnum/a0 * 2.045d0*xldt

      dum    = 1.875d8*xl + 1.653d8*xl2 + 8.449d8*xl3 - 1.604d8*xl4
      dumdt  = xldt*(1.875d8 + 2.0d0*1.653d8*xl + 3.0d0*8.449d8*xl2 &
               - 4.0d0*1.604d8*xl3)

      z      = 1.0d0/dum
      xden   = 1.0d0 + rm*z
      xdendt =  -rm*z*z*dumdt
      xdendd =  rmdd*z
      xdenda =  rmda*z
      xdendz =  rmdz*z

      z      = 1.0d0/xden
      qphot = xnum*z
      qphotdt = (xnumdt - qphot*xdendt)*z
      dum      = -qphot*z
      qphotdd = dum*xdendd
      qphotda = dum*xdenda
      qphotdz = dum*xdendz

! equation 3.2
      sphot   = xl5 * fphot
      sphotdt = 5.0d0*xl4*xldt*fphot + xl5*fphotdt
      sphotdd = xl5*fphotdd
      sphotda = xl5*fphotda
      sphotdz = xl5*fphotdz

      a1      = sphot
      sphot   = rm*a1
      sphotdt = rm*sphotdt
      sphotdd = rm*sphotdd + rmdd*a1
      sphotda = rm*sphotda + rmda*a1
      sphotdz = rm*sphotdz + rmdz*a1

      a1      = tfac4*(1.0d0 - tfac3 * qphot)
      a2      = -tfac4*tfac3

      a3      = sphot
      sphot   = a1*a3
      sphotdt = a1*sphotdt + a2*qphotdt*a3
      sphotdd = a1*sphotdd + a2*qphotdd*a3
      sphotda = a1*sphotda + a2*qphotda*a3
      sphotdz = a1*sphotdz + a2*qphotdz*a3

      if (sphot .le. 0.0) then
       sphot   = 0.0d0
       sphotdt = 0.0d0
       sphotdd = 0.0d0
       sphotda = 0.0d0
       sphotdz = 0.0d0
      end if





! bremsstrahlung neutrino section
! for reactions like e- + (z,a) => e- + (z,a) + nu + nubar
!                    n  + n     => n + n + nu + nubar
!                    n  + p     => n + p + nu + nubar
! equation 4.3

      den6   = den * 1.0d-6
      t8     = temp * 1.0d-8
      t812   = sqrt(t8)
      t832   = t8 * t812
      t82    = t8*t8
      t83    = t82*t8
      t85    = t82*t83
      t86    = t85*t8
      t8m1   = 1.0d0/t8
      t8m2   = t8m1*t8m1
      t8m3   = t8m2*t8m1
      t8m5   = t8m3*t8m2
      t8m6   = t8m5*t8m1


      tfermi = 5.9302d9*(sqrt(1.0d0+1.018d0*(den6*ye)**twoth)-1.0d0)

! "weak" degenerate electrons only
      if (temp .gt. 0.3d0 * tfermi) then

! equation 5.3
       dum   = 7.05d6 * t832 + 5.12d4 * t83
       dumdt = (1.5d0*7.05d6*t812 + 3.0d0*5.12d4*t82)*1.0d-8

       z     = 1.0d0/dum
       eta   = rm*z
       etadt = -rm*z*z*dumdt
       etadd = rmdd*z
       etada = rmda*z
       etadz = rmdz*z

       etam1 = 1.0d0/eta
       etam2 = etam1 * etam1
       etam3 = etam2 * etam1


! equation 5.2
       a0    = 23.5d0 + 6.83d4*t8m2 + 7.81d8*t8m5
       f0    = (-2.0d0*6.83d4*t8m3 - 5.0d0*7.81d8*t8m6)*1.0d-8
       xnum  = 1.0d0/a0

       dum   = 1.0d0 + 1.47d0*etam1 + 3.29d-2*etam2
       z     = -1.47d0*etam2 - 2.0d0*3.29d-2*etam3
       dumdt = z*etadt
       dumdd = z*etadd
       dumda = z*etada
       dumdz = z*etadz

       c00   = 1.26d0 * (1.0d0+etam1)
       z     = -1.26d0*etam2
       c01   = z*etadt
       c02   = z*etadd
       c03   = z*etada
       c04   = z*etadz

       z      = 1.0d0/dum
       xden   = c00*z
       xdendt = (c01 - xden*dumdt)*z
       xdendd = (c02 - xden*dumdd)*z
       xdenda = (c03 - xden*dumda)*z
       xdendz = (c04 - xden*dumdz)*z

       fbrem   = xnum + xden
       fbremdt = -xnum*xnum*f0 + xdendt
       fbremdd = xdendd
       fbremda = xdenda
       fbremdz = xdendz


! equation 5.9
       a0    = 230.0d0 + 6.7d5*t8m2 + 7.66d9*t8m5
       f0    = (-2.0d0*6.7d5*t8m3 - 5.0d0*7.66d9*t8m6)*1.0d-8

       z     = 1.0d0 + rm*1.0d-9
       dum   = a0*z
       dumdt = f0*z
       z     = a0*1.0d-9
       dumdd = z*rmdd
       dumda = z*rmda
       dumdz = z*rmdz

       xnum   = 1.0d0/dum
       z      = -xnum*xnum
       xnumdt = z*dumdt
       xnumdd = z*dumdd
       xnumda = z*dumda
       xnumdz = z*dumdz

       c00   = 7.75d5*t832 + 247.0d0*t8**(3.85d0)
       dd00  = (1.5d0*7.75d5*t812 + 3.85d0*247.0d0*t8**(2.85d0))*1.0d-8

       c01   = 4.07d0 + 0.0240d0 * t8**(1.4d0)
       dd01  = 1.4d0*0.0240d0*t8**(0.4d0)*1.0d-8

       c02   = 4.59d-5 * t8**(-0.110d0)
       dd02  = -0.11d0*4.59d-5 * t8**(-1.11d0)*1.0d-8

       z     = den**(0.656d0)
       dum   = c00*rmi  + c01  + c02*z
       dumdt = dd00*rmi + dd01 + dd02*z
       z     = -c00*rmi*rmi
       dumdd = z*rmdd + 0.656d0*c02*den**(-0.454d0)
       dumda = z*rmda
       dumdz = z*rmdz

       xden  = 1.0d0/dum
       z      = -xden*xden
       xdendt = z*dumdt
       xdendd = z*dumdd
       xdenda = z*dumda
       xdendz = z*dumdz

       gbrem   = xnum + xden
       gbremdt = xnumdt + xdendt
       gbremdd = xnumdd + xdendd
       gbremda = xnumda + xdenda
       gbremdz = xnumdz + xdendz


! equation 5.1
       dum    = 0.5738d0*zbar*ye*t86*den
       dumdt  = 0.5738d0*zbar*ye*6.0d0*t85*den*1.0d-8
       dumdd  = 0.5738d0*zbar*ye*t86
       dumda  = -dum*abari
       dumdz  = 0.5738d0*2.0d0*ye*t86*den

       z       = tfac4*fbrem - tfac5*gbrem
       sbrem   = dum * z
       sbremdt = dumdt*z + dum*(tfac4*fbremdt - tfac5*gbremdt)
       sbremdd = dumdd*z + dum*(tfac4*fbremdd - tfac5*gbremdd)
       sbremda = dumda*z + dum*(tfac4*fbremda - tfac5*gbremda)
       sbremdz = dumdz*z + dum*(tfac4*fbremdz - tfac5*gbremdz)




! liquid metal with c12 parameters (not too different for other elements)
! equation 5.18 and 5.16

      else
       u     = fac3 * (log10(den) - 3.0d0)
       a0    = iln10*fac3*deni

! compute the expensive trig functions of equation 5.21 only once
       cos1 = cos(u)
       cos2 = cos(2.0d0*u)
       cos3 = cos(3.0d0*u)
       cos4 = cos(4.0d0*u)
       cos5 = cos(5.0d0*u)

       sin1 = sin(u)
       sin2 = sin(2.0d0*u)
       sin3 = sin(3.0d0*u)
       sin4 = sin(4.0d0*u)
       sin5 = sin(5.0d0*u)

! equation 5.21
       fb =  0.5d0 * 0.17946d0  + 0.00945d0*u + 0.34529d0 &
             - 0.05821d0*cos1 - 0.04969d0*sin1 &
             - 0.01089d0*cos2 - 0.01584d0*sin2 &
             - 0.01147d0*cos3 - 0.00504d0*sin3 &
             - 0.00656d0*cos4 - 0.00281d0*sin4 &
             - 0.00519d0*cos5

       c00 =  a0*(0.00945d0 &
             + 0.05821d0*sin1       - 0.04969d0*cos1 &
             + 0.01089d0*sin2*2.0d0 - 0.01584d0*cos2*2.0d0 &
             + 0.01147d0*sin3*3.0d0 - 0.00504d0*cos3*3.0d0 &
             + 0.00656d0*sin4*4.0d0 - 0.00281d0*cos4*4.0d0 &
             + 0.00519d0*sin5*5.0d0)


! equation 5.22
       ft =  0.5d0 * 0.06781d0 - 0.02342d0*u + 0.24819d0 &
             - 0.00944d0*cos1 - 0.02213d0*sin1 &
             - 0.01289d0*cos2 - 0.01136d0*sin2 &
             - 0.00589d0*cos3 - 0.00467d0*sin3 &
             - 0.00404d0*cos4 - 0.00131d0*sin4 &
             - 0.00330d0*cos5

       c01 = a0*(-0.02342d0 &
             + 0.00944d0*sin1       - 0.02213d0*cos1 &
             + 0.01289d0*sin2*2.0d0 - 0.01136d0*cos2*2.0d0 &
             + 0.00589d0*sin3*3.0d0 - 0.00467d0*cos3*3.0d0 &
             + 0.00404d0*sin4*4.0d0 - 0.00131d0*cos4*4.0d0 &
             + 0.00330d0*sin5*5.0d0)


! equation 5.23
       gb =  0.5d0 * 0.00766d0 - 0.01259d0*u + 0.07917d0 &
             - 0.00710d0*cos1 + 0.02300d0*sin1 &
             - 0.00028d0*cos2 - 0.01078d0*sin2 &
             + 0.00232d0*cos3 + 0.00118d0*sin3 &
             + 0.00044d0*cos4 - 0.00089d0*sin4 &
             + 0.00158d0*cos5

       c02 = a0*(-0.01259d0 &
             + 0.00710d0*sin1       + 0.02300d0*cos1 &
             + 0.00028d0*sin2*2.0d0 - 0.01078d0*cos2*2.0d0 &
             - 0.00232d0*sin3*3.0d0 + 0.00118d0*cos3*3.0d0 &
             - 0.00044d0*sin4*4.0d0 - 0.00089d0*cos4*4.0d0 &
             - 0.00158d0*sin5*5.0d0)


! equation 5.24
       gt =  -0.5d0 * 0.00769d0  - 0.00829d0*u + 0.05211d0 &
             + 0.00356d0*cos1 + 0.01052d0*sin1 &
             - 0.00184d0*cos2 - 0.00354d0*sin2 &
             + 0.00146d0*cos3 - 0.00014d0*sin3 &
             + 0.00031d0*cos4 - 0.00018d0*sin4 &
             + 0.00069d0*cos5

       c03 = a0*(-0.00829d0 &
             - 0.00356d0*sin1       + 0.01052d0*cos1 &
             + 0.00184d0*sin2*2.0d0 - 0.00354d0*cos2*2.0d0 &
             - 0.00146d0*sin3*3.0d0 - 0.00014d0*cos3*3.0d0 &
             - 0.00031d0*sin4*4.0d0 - 0.00018d0*cos4*4.0d0 &
             - 0.00069d0*sin5*5.0d0)


       dum   = 2.275d-1 * zbar * zbar*t8m1 * (den6*abari)**oneth
       dumdt = -dum*tempi
       dumdd = oneth*dum*deni
       dumda = -oneth*dum*abari
       dumdz = 2.0d0*dum*zbari

       gm1   = 1.0d0/dum
       gm2   = gm1*gm1
       gm13  = gm1**oneth
       gm23  = gm13 * gm13
       gm43  = gm13*gm1
       gm53  = gm23*gm1


! equation 5.25 and 5.26
       v  = -0.05483d0 - 0.01946d0*gm13 + 1.86310d0*gm23 - 0.78873d0*gm1
       a0 = oneth*0.01946d0*gm43 - twoth*1.86310d0*gm53 + 0.78873d0*gm2

       w  = -0.06711d0 + 0.06859d0*gm13 + 1.74360d0*gm23 - 0.74498d0*gm1
       a1 = -oneth*0.06859d0*gm43 - twoth*1.74360d0*gm53 + 0.74498d0*gm2


! equation 5.19 and 5.20
       fliq   = v*fb + (1.0d0 - v)*ft
       fliqdt = a0*dumdt*(fb - ft)
       fliqdd = a0*dumdd*(fb - ft) + v*c00 + (1.0d0 - v)*c01
       fliqda = a0*dumda*(fb - ft)
       fliqdz = a0*dumdz*(fb - ft)

       gliq   = w*gb + (1.0d0 - w)*gt
       gliqdt = a1*dumdt*(gb - gt)
       gliqdd = a1*dumdd*(gb - gt) + w*c02 + (1.0d0 - w)*c03
       gliqda = a1*dumda*(gb - gt)
       gliqdz = a1*dumdz*(gb - gt)


! equation 5.17
       dum    = 0.5738d0*zbar*ye*t86*den
       dumdt  = 0.5738d0*zbar*ye*6.0d0*t85*den*1.0d-8
       dumdd  = 0.5738d0*zbar*ye*t86
       dumda  = -dum*abari
       dumdz  = 0.5738d0*2.0d0*ye*t86*den

       z       = tfac4*fliq - tfac5*gliq
       sbrem   = dum * z
       sbremdt = dumdt*z + dum*(tfac4*fliqdt - tfac5*gliqdt)
       sbremdd = dumdd*z + dum*(tfac4*fliqdd - tfac5*gliqdd)
       sbremda = dumda*z + dum*(tfac4*fliqda - tfac5*gliqda)
       sbremdz = dumdz*z + dum*(tfac4*fliqdz - tfac5*gliqdz)

      end if




! recombination neutrino section
! for reactions like e- (continuum) => e- (bound) + nu_e + nubar_e
! equation 6.11 solved for nu
      xnum   = 1.10520d8 * den * ye /(temp*sqrt(temp))
      xnumdt = -1.50d0*xnum*tempi
      xnumdd = xnum*deni
      xnumda = -xnum*abari
      xnumdz = xnum*zbari

! the chemical potential
      nu   = ifermi12(xnum)

! a0 is d(nu)/d(xnum)
      a0 = 1.0d0/(0.5d0*zfermim12(nu))
      nudt = a0*xnumdt
      nudd = a0*xnumdd
      nuda = a0*xnumda
      nudz = a0*xnumdz

      nu2  = nu * nu
      nu3  = nu2 * nu

! table 12
      if (nu .ge. -20.0  .and. nu .lt. 0.0) then
       a1 = 1.51d-2
       a2 = 2.42d-1
       a3 = 1.21d0
       b  = 3.71d-2
       c  = 9.06e-1
       d  = 9.28d-1
       f1 = 0.0d0
       f2 = 0.0d0
       f3 = 0.0d0
      else if (nu .ge. 0.0  .and. nu .le. 10.0) then
       a1 = 1.23d-2
       a2 = 2.66d-1
       a3 = 1.30d0
       b  = 1.17d-1
       c  = 8.97e-1
       d  = 1.77d-1
       f1 = -1.20d-2
       f2 = 2.29d-2
       f3 = -1.04d-3
      end if


! equation 6.7, 6.13 and 6.14
      if (nu .ge. -20.0  .and.  nu .le. 10.0) then

       zeta   = 1.579d5*zbar*zbar*tempi
       zetadt = -zeta*tempi
       zetadd = 0.0d0
       zetada = 0.0d0
       zetadz = 2.0d0*zeta*zbari

       c00    = 1.0d0/(1.0d0 + f1*nu + f2*nu2 + f3*nu3)
       c01    = f1 + f2*2.0d0*nu + f3*3.0d0*nu2
       dum    = zeta*c00
       dumdt  = zetadt*c00 + zeta*c01*nudt
       dumdd  = zeta*c01*nudd
       dumda  = zeta*c01*nuda
       dumdz  = zetadz*c00 + zeta*c01*nudz


       z      = 1.0d0/dum
       dd00   = dum**(-2.25)
       dd01   = dum**(-4.55)
       c00    = a1*z + a2*dd00 + a3*dd01
       c01    = -(a1*z + 2.25*a2*dd00 + 4.55*a3*dd01)*z


       z      = exp(c*nu)
       dd00   = b*z*(1.0d0 + d*dum)
       gum    = 1.0d0 + dd00
       gumdt  = dd00*c*nudt + b*z*d*dumdt
       gumdd  = dd00*c*nudd + b*z*d*dumdd
       gumda  = dd00*c*nuda + b*z*d*dumda
       gumdz  = dd00*c*nudz + b*z*d*dumdz


       z   = exp(nu)
       a1  = 1.0d0/gum

       bigj   = c00 * z * a1
       bigjdt = c01*dumdt*z*a1 + c00*z*nudt*a1 - c00*z*a1*a1 * gumdt
       bigjdd = c01*dumdd*z*a1 + c00*z*nudd*a1 - c00*z*a1*a1 * gumdd
       bigjda = c01*dumda*z*a1 + c00*z*nuda*a1 - c00*z*a1*a1 * gumda
       bigjdz = c01*dumdz*z*a1 + c00*z*nudz*a1 - c00*z*a1*a1 * gumdz


! equation 6.5
       z     = exp(zeta + nu)
       dum   = 1.0d0 + z
       a1    = 1.0d0/dum
       a2    = 1.0d0/bigj

       sreco   = tfac6 * 2.649d-18 * ye * zbar**13 * den * bigj*a1
       srecodt = sreco*(bigjdt*a2 - z*(zetadt + nudt)*a1)
       srecodd = sreco*(1.0d0*deni + bigjdd*a2 - z*(zetadd + nudd)*a1)
       srecoda = sreco*(-1.0d0*abari + bigjda*a2 - z*(zetada+nuda)*a1)
       srecodz = sreco*(14.0d0*zbari + bigjdz*a2 - z*(zetadz+nudz)*a1)

      end if


! convert from erg/cm^3/s to erg/g/s
! comment these out to duplicate the itoh et al plots

      spair   = spair*deni
      spairdt = spairdt*deni
      spairdd = spairdd*deni - spair*deni
      spairda = spairda*deni
      spairdz = spairdz*deni

      splas   = splas*deni
      splasdt = splasdt*deni
      splasdd = splasdd*deni - splas*deni
      splasda = splasda*deni
      splasdz = splasdz*deni

      sphot   = sphot*deni
      sphotdt = sphotdt*deni
      sphotdd = sphotdd*deni - sphot*deni
      sphotda = sphotda*deni
      sphotdz = sphotdz*deni

      sbrem   = sbrem*deni
      sbremdt = sbremdt*deni
      sbremdd = sbremdd*deni - sbrem*deni
      sbremda = sbremda*deni
      sbremdz = sbremdz*deni

      sreco   = sreco*deni
      srecodt = srecodt*deni
      srecodd = srecodd*deni - sreco*deni
      srecoda = srecoda*deni
      srecodz = srecodz*deni


! the total neutrino loss rate
      snu    =  splas + spair + sphot + sbrem + sreco
      dsnudt =  splasdt + spairdt + sphotdt + sbremdt + srecodt
      dsnudd =  splasdd + spairdd + sphotdd + sbremdd + srecodd
      dsnuda =  splasda + spairda + sphotda + sbremda + srecoda
      dsnudz =  splasdz + spairdz + sphotdz + sbremdz + srecodz

      return
      end






      double precision function ifermi12(f)
      include 'implno.dek'

! this routine applies a rational function expansion to get the inverse
! fermi-dirac integral of order 1/2 when it is equal to f.
! maximum error is 4.19d-9.   reference: antia apjs 84,101 1993

! declare
      integer          i,m1,k1,m2,k2
      double precision f,an,a1(12),b1(12),a2(12),b2(12),rn,den,ff


! load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /0.5d0, 4, 3, 6, 5/
      data  (a1(i),i=1,5)/ 1.999266880833d4,   5.702479099336d3, &
           6.610132843877d2,   3.818838129486d1, &
           1.0d0/
      data  (b1(i),i=1,4)/ 1.771804140488d4,  -2.014785161019d3, &
           9.130355392717d1,  -1.670718177489d0/
      data  (a2(i),i=1,7)/-1.277060388085d-2,  7.187946804945d-2, &
                          -4.262314235106d-1,  4.997559426872d-1, &
                          -1.285579118012d0,  -3.930805454272d-1, &
           1.0d0/
      data  (b2(i),i=1,6)/-9.745794806288d-3,  5.485432756838d-2, &
                          -3.299466243260d-1,  4.077841975923d-1, &
                          -1.145531476975d0,  -6.067091689181d-2/


      if (f .lt. 4.0d0) then
       rn  = f + a1(m1)
       do i=m1-1,1,-1
        rn  = rn*f + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*f + b1(i)
       enddo
       ifermi12 = log(f * rn/den)

      else
       ff = 1.0d0/f**(1.0d0/(1.0d0 + an))
       rn = ff + a2(m2)
       do i=m2-1,1,-1
        rn = rn*ff + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*ff + b2(i)
       enddo
       ifermi12 = rn/(den*ff)
      end if
      return
      end






      double precision function zfermim12(x)
      include 'implno.dek'

! this routine applies a rational function expansion to get the fermi-dirac
! integral of order -1/2 evaluated at x. maximum error is 1.23d-12.
! reference: antia apjs 84,101 1993

! declare
      integer          i,m1,k1,m2,k2
      double precision x,an,a1(12),b1(12),a2(12),b2(12),rn,den,xx

! load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /-0.5d0, 7, 7, 11, 11/
      data  (a1(i),i=1,8)/ 1.71446374704454d7,    3.88148302324068d7, &
                           3.16743385304962d7,    1.14587609192151d7, &
                           1.83696370756153d6,    1.14980998186874d5, &
                           1.98276889924768d3,    1.0d0/
      data  (b1(i),i=1,8)/ 9.67282587452899d6,    2.87386436731785d7, &
                           3.26070130734158d7,    1.77657027846367d7, &
                           4.81648022267831d6,    6.13709569333207d5, &
                           3.13595854332114d4,    4.35061725080755d2/
      data (a2(i),i=1,12)/-4.46620341924942d-15, -1.58654991146236d-12, &
                          -4.44467627042232d-10, -6.84738791621745d-8, &
                          -6.64932238528105d-6,  -3.69976170193942d-4, &
                          -1.12295393687006d-2,  -1.60926102124442d-1, &
                          -8.52408612877447d-1,  -7.45519953763928d-1, &
                           2.98435207466372d0,    1.0d0/
      data (b2(i),i=1,12)/-2.23310170962369d-15, -7.94193282071464d-13, &
                          -2.22564376956228d-10, -3.43299431079845d-8, &
                          -3.33919612678907d-6,  -1.86432212187088d-4, &
                          -5.69764436880529d-3,  -8.34904593067194d-2, &
                          -4.78770844009440d-1,  -4.99759250374148d-1, &
                           1.86795964993052d0,    4.16485970495288d-1/


      if (x .lt. 2.0d0) then
       xx = exp(x)
       rn = xx + a1(m1)
       do i=m1-1,1,-1
        rn = rn*xx + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*xx + b1(i)
       enddo
       zfermim12 = xx * rn/den
!
      else
       xx = 1.0d0/(x*x)
       rn = xx + a2(m2)
       do i=m2-1,1,-1
        rn = rn*xx + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*xx + b2(i)
       enddo
       zfermim12 = sqrt(x)*rn/den
      end if
      return
      end




!---------------------------------------------------------------------
      subroutine net_output(kount,x,y,derivs)
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'
      include 'burn_common.dek'
      include 'network.dek'
      include 'cjdet.dek'

! writes the output

! declare the pass
      external         derivs
      integer          kount
      double precision x,y(*)


! local variables
      character*8      atim
      character*9      adat
      character*80     string
      integer          k,kk,j,lop,ilop,jrem,kb,ke,nn,lenstr
      double precision sum,xcons,ycons,yex,ydum(abignet), &
                       dydt_dum(nzmax*abignet),xdum(abignet), &
                       abar,zbar,wbar,ye,xcess,zero,tdum,ddum,pdum, &
                       ener,denerdt,zc12,xc12,ff, &
                       chem_pot(nzmax*abignet),chem_sum, &
                       ydum_sav(nzmax*abignet),posx,velx
      parameter        (zero = 0.0d0)


! for nse
      integer          igues
      double precision xmun,xmup,t9,tau_nse,tau_qse,taud


! popular format statements
01    format(1x,'*',t13,a,t33,a,t47,a,t61,a,t75,a,t89,a, &
                    t103,a,t117,a,t131,a,t145,a,t159,a)
03    format(a30,i4.4,a2,i8,a)
04    format(1x,i6,1pe20.12,1p15e14.6)
05    format(1x,i6,1pe20.12,1p12e14.6)
07    format(1x,'* ',a,5(a,1pe11.3))



!      write(6,*) kount,neqs,nzone

! initialize the files with their headers
      if (kount .eq. 1) then



! for every spatial zone
       do k=1,max(1,nzone)
        kk = neqs*(k-1)

! logical unit 22 records the energetics
        write(string,03) hfile,0,'_z',k,'.dat'
        call sqeeze(string)
        call today(adat,atim)
        open (unit=22, file=string, status='unknown')


! logical unit 23 records the thermodynamics
        write(string,03) hfile,1,'_z',k,'.dat'
        call sqeeze(string)
        open (unit=23, file=string, status='unknown')


         write(22,01) adat,atim
         write(23,01) adat,atim

        if (one_step) then
         write(22,07) 'one_step:','  btemp=',btemp,' bden=',bden
         write(23,07) 'one_step:','  btemp=',btemp,' bden=',bden

        else if (hydrostatic) then
         write(22,07) 'hydrostatic:','  btemp=',btemp,' bden=',bden
         write(23,07) 'hydrostatic:','  btemp=',btemp,' bden=',bden

        else if (expansion) then
         write(22,07) 'expansion:','  temp0=',temp0,' den0=',den0, &
                     ' temp_stop=',temp_stop
         write(23,07) 'expansion:','  temp0=',temp0,' den0=',den0, &
                     ' temp_stop=',temp_stop

        else if (self_heat_const_den) then
         write(22,07) 'self_heat:','   temp0=',y(itemp+kk), &
                                  '   den0=',y(iden+kk)
         write(23,07) 'self_heat:','   temp0=',y(itemp+kk), &
                                  '   den0=',y(iden+kk)

        else if (self_heat_const_pres) then
         write(22,07) 'self_heat:','   temp0=',y(itemp+kk), &
                                  '   den0=',y(iden+kk)
         write(23,07) 'self_heat:','   temp0=',y(itemp+kk), &
                                  '   den0=',y(iden+kk)

        else if (bbang) then
         write(22,07) 'big bang:','   eta=',eta1, &
                                  '   Nnu=',xnnu, &
                                  '   H0=',hubble, &
                                  '   Tcmb=',cmbtemp
         write(23,07) 'big_bang:','   eta=',eta1, &
                                  '   Nnu=',xnnu, &
                                  '   H0=',hubble, &
                                  '   Tcmb=',cmbtemp

        else if (detonation) then
         write(22,07) 'detonation:',' temp0=',temp_up, &
                                  '   den0=',den_up, &
                                  '   pres0=',pres_up, &
                                  '   mach=',mach_sh
         write(22,07) '           ',' temp_sh=',temp_sh, &
                                  '   den_sh=',den_sh, &
                                  '   pres_sh=',pres_sh, &
                                  '   vel_sh=',vel_sh, &
                                  '   cs_sh=',cs_sh
         write(22,07) '           ',' temp_cj=',temp_cj, &
                                  '   den_cj=',den_cj, &
                                  '   pres_cj=',pres_cj, &
                                  '   vel_cj=',vel_cj, &
                                  '   cs_cj=',cs_cj

         write(23,07) 'detonation:','  temp0=',temp_up, &
                                  '   den0=',den_up, &
                                  '   pres0=',pres_up, &
                                  '   mach=',mach_sh
         write(23,07) '           ',' temp_sh=',temp_sh, &
                                  '   den_sh=',den_sh, &
                                  '   pres_sh=',pres_sh, &
                                  '   vel_sh=',vel_sh, &
                                  '   cs_sh=',cs_sh
         write(23,07) '           ',' temp_cj=',temp_cj, &
                                  '   den_cj=',den_cj, &
                                  '   pres_cj=',pres_cj, &
                                  '   vel_cj=',vel_cj, &
                                  '   cs_cj=',cs_cj

        else if (trho_hist) then
         call update2(zero,tdum,ddum)
         write(22,07) 'trho_hist:','  mass interior =',mint, &
                                   '  shell mass =',mshell
         write(23,07) 'trho_hist:',' mass interior =',mint, &
                                   '  shell mass =',mshell

        else if (pt_hist) then
         call update3(zero,tdum,pdum)
         write(22,07) 'pt_hist:  ','  mass interior =',mint, &
                                   '  shell mass =',mshell
         write(23,07) 'pt_hist:  ',' mass interior =',mint, &
                                   '  shell mass =',mshell
        end if


        write(22,01) 'time','temp','den','ener','sdot','sneut', &
                     's-snu','ye','1-sum'

        write(23,01) 'time','pos','vel','temp','den','pres','ener', &
                     'entr','cs'


! close up the files
        close(unit=22)
        close(unit=23)


! end of spatial loop
       enddo



! if we are doing an nse analysis, we'll write out another file

       if (nse_analysis .eq. 1) then

! for every spatial zone
       do k=1,max(1,nzone)
        kk = neqs*(k-1)

! logical unit 25 records the nse analysis
        write(string,03) hfile,0,'_z',k,'_nse.dat'
        call sqeeze(string)
        call today(adat,atim)
        open (unit=25, file=string, status='unknown')

        write(25,01) adat,atim

        if (one_step) then
         write(25,07) 'one_step:','  btemp=',btemp,' bden=',bden

        else if (hydrostatic) then
         write(25,07) 'hydrostatic:','  btemp=',btemp,' bden=',bden

        else if (expansion) then
         write(25,07) 'expansion:','  temp0=',temp0,' den0=',den0, &
                     ' temp_stop=',temp_stop

        else if (self_heat_const_den) then
         write(25,07) 'self_heat:','   temp0=',y(itemp+kk), &
                                  '   den0=',y(iden+kk)

        else if (self_heat_const_pres) then
         write(25,07) 'self_heat:','   temp0=',y(itemp+kk), &
                                  '   den0=',y(iden+kk)

        else if (bbang) then
         write(25,07) 'big bang:','   eta=',eta1, &
                                  '   Nnu=',xnnu, &
                                  '   H0=',hubble, &
                                  '   Tcmb=',cmbtemp

        else if (detonation) then
         write(25,07) 'detonation:','  temp0=',temp_up, &
                                  '   den0=',den_up, &
                                  '   pres0=',pres_up, &
                                  '   mach=',mach_sh
         write(25,07) '           ',' temp_sh=',temp_sh, &
                                  '   den_sh=',den_sh, &
                                  '   pres_sh=',pres_sh, &
                                  '   vel_sh=',vel_sh, &
                                  '   cs_sh=',cs_sh
         write(25,07) '           ',' temp_cj=',temp_cj, &
                                  '   den_cj=',den_cj, &
                                  '   pres_cj=',pres_cj, &
                                  '   vel_cj=',vel_cj, &
                                  '   cs_cj=',cs_cj

        else if (trho_hist) then
         call update2(zero,tdum,ddum)
         write(25,07) 'trho_hist:','  mass interior =',mint, &
                                   '  shell mass =',mshell

        else if (pt_hist) then
         call update3(zero,tdum,pdum)
         write(25,07) 'pt_hist:  ','  mass interior =',mint, &
                                   '  shell mass =',mshell
        end if


        write(25,01) 'time','temp','den','ye','tqse','tnse','delta', &
                     '1-sum'


! close up the files
        close(unit=25)

! end of spatial loop
       enddo

       end if



! done writing thermodynamic headers



! for every spatial zone
       do k=1,max(1,nzone)
        kk = neqs*(k-1)


! write out the isotopic mass fractions in blocks of 8
! lop is how many groups of 8 exist; jrem is the remainder
        lop  = ionmax/8
        jrem  = ionmax - 8*lop
        do ilop = 1,lop+1
         kb = 1 + 8*(ilop-1)
         ke = 8 + 8*(ilop-1)
         if (ilop .eq. lop+1  .and. jrem .eq. 0) goto 50
         if (ilop .eq. lop+1) ke = ionmax


! logical unit 34 records the abundance evolution
! open the output file
         write(string,03) hfile,ilop+1,'_z',k,'.dat'
         call sqeeze(string)
         open (unit=34, file=string, status='unknown')

         write(34,01) adat,atim


        if (one_step) then
         write(34,07) 'one_step:','  btemp=',btemp,' bden=',bden

        else if (hydrostatic) then
         write(34,07) 'hydrostatic:','  btemp=',btemp,' bden=',bden

        else if (expansion) then
         write(34,07) 'expansion:','  temp0=',temp0,' den0=',den0, &
                                  ' temp_stop=',temp_stop

        else if (self_heat_const_den) then
         write(34,07) 'self_heat:','   temp0=',y(itemp+kk), &
                                  '   den0=',y(iden+kk)

        else if (bbang) then
         write(34,07) 'big bang:','   eta=',eta1, &
                                  '   Nnu=',xnnu, &
                                  '   H0=',hubble, &
                                  '   Tcmb=',cmbtemp

        else if (detonation) then
         write(34,07) 'detonation:','  temp0=',temp_up, &
                                  '   den0=',den_up, &
                                  '   pres0=',pres_up, &
                                  '   mach=',mach_sh
         write(34,07) '           ',' temp_sh=',temp_sh, &
                                  '   den_sh=',den_sh, &
                                  '   pres_sh=',pres_sh, &
                                  '   vel_sh=',vel_sh, &
                                  '   cs_sh=',cs_sh
         write(34,07) '           ',' temp_cj=',temp_cj, &
                                  '   den_cj=',den_cj, &
                                  '   pres_cj=',pres_cj, &
                                  '   vel_cj=',vel_cj, &
                                  '   cs_cj=',cs_cj


        else if (trho_hist) then
         call update2(zero,tdum,ddum)
         write(34,07) 'trho_hist:','  mass interior =',mint, &
                                   '  shell mass =',mshell

        else if (pt_hist) then
         call update3(zero,tdum,ddum)
         write(34,07) 'pt_hist:  ','  mass interior =',mint, &
                                   '  shell mass =',mshell

        end if

        write(34,01) 'time',(ionam(nn), nn=kb,ke)

        close(unit=34)
 50     continue
       enddo

! end of the spatial loop
      enddo



! if we are doing an nse analysis, we'll write out another
! set of abundance file

       if (nse_analysis .eq. 1) then


! for every spatial zone
       do k=1,max(1,nzone)
        kk = neqs*(k-1)

! write out the isotopic mass fractions in blocks of 8
! lop is how many groups of 8 exist; jrem is the remainder
        lop  = ionmax/8
        jrem  = ionmax - 8*lop
        do ilop = 1,lop+1
         kb = 1 + 8*(ilop-1)
         ke = 8 + 8*(ilop-1)
         if (ilop .eq. lop+1  .and. jrem .eq. 0) goto 60
         if (ilop .eq. lop+1) ke = ionmax


! logical unit 35 records the abundance evolution
! open the output file
         write(string,03) hfile,ilop+1,'_z',k,'_nse.dat'
         call sqeeze(string)
         open (unit=35, file=string, status='unknown')

         write(35,01) adat,atim


        if (one_step) then
         write(35,07) 'one_step:','  btemp=',btemp,' bden=',bden

        else if (hydrostatic) then
         write(35,07) 'hydrostatic:','  btemp=',btemp,' bden=',bden

        else if (expansion) then
         write(35,07) 'expansion:','  temp0=',temp0,' den0=',den0, &
                                  ' temp_stop=',temp_stop

        else if (self_heat_const_den) then
         write(35,07) 'self_heat:','   temp0=',y(itemp+kk), &
                                  '   den0=',y(iden+kk)

        else if (bbang) then
         write(35,07) 'big bang:','   eta=',eta1, &
                                  '   Nnu=',xnnu, &
                                  '   H0=',hubble, &
                                  '   Tcmb=',cmbtemp

        else if (detonation) then
         write(35,07) 'detonation:','  temp0=',temp_up, &
                                  '   den0=',den_up, &
                                  '   pres0=',pres_up, &
                                  '   mach=',mach_sh
         write(35,07) '           ',' temp_sh=',temp_sh, &
                                  '   den_sh=',den_sh, &
                                  '   pres_sh=',pres_sh, &
                                  '   vel_sh=',vel_sh, &
                                  '   cs_sh=',cs_sh
         write(35,07) '           ',' temp_cj=',temp_cj, &
                                  '   den_cj=',den_cj, &
                                  '   pres_cj=',pres_cj, &
                                  '   vel_cj=',vel_cj, &
                                  '   cs_cj=',cs_cj


        else if (trho_hist) then
         call update2(zero,tdum,ddum)
         write(35,07) 'trho_hist:','  mass interior =',mint, &
                                   '  shell mass =',mshell

        else if (pt_hist) then
         call update3(zero,tdum,ddum)
         write(35,07) 'pt_hist:  ','  mass interior =',mint, &
                                   '  shell mass =',mshell

        end if

        write(35,01) 'time',(ionam(nn), nn=kb,ke)

        close(unit=35)
 60     continue
       enddo

! end of the spatial loop and nse analyis test if
      enddo
      end if

!       write(6,*) 'wrote mass fraction headers'

! end of the file initialization
      end if

!      write(6,*) 'done with initialization'







! normal execution starts here

! for any time point

! for every spatial zone
      do k=1,max(1,nzone)
       kk = neqs*(k-1)

! open the files in append mode (f77) or position mode (f90)

! energetics file
       write(string,03) hfile,0,'_z',k,'.dat'
       call sqeeze(string)
!       open (unit=22, file=string, status='old', access='append')
       open (unit=22, file=string, status='old', position='append')


! thermodynamics file
       write(string,03) hfile,1,'_z',k,'.dat'
       call sqeeze(string)
!       open (unit=23, file=string, status='old', access='append')
       open (unit=23, file=string, status='old', position='append')


! form the mass fractions
       do j=1,ionmax
        xdum(j) = min(1.0d0,max(y(j+kk)*aion(j),1.0d-30))
       enddo


! mass conservation
       sum = 0.0d0
       do j=1,ionmax
        sum = sum + xdum(j)
       enddo
       sum = 1.0d0 - sum
       xcons = sum


! y sum
!       sum = 0.0d0
!       do j=1,ionmax
!        if (zion(j) .gt. 2.0) then
!         sum = sum + max(y(j+kk),1.0d-30)
!        endif
!       enddo
!       ycons = sum


! get ye using normalized mass fractions
       sum = 0.0d0
       do j=1,ionmax
        sum = sum + xdum(j)
       enddo
       sum = 1.0d0/sum
       do j=1,ionmax
        xdum(j) = min(1.0d0,max(sum*xdum(j),1.0d-30))
       enddo


! get abar, zbar and a few other composition variables
       call azbar(xdum(ionbeg),aion(ionbeg),zion(ionbeg),wion(ionbeg),ionmax, &
                  ydum(ionbeg),abar,zbar,wbar,yex,xcess)




! get the right hand sides, exact energy generation rate and so on
       if (nse_on .eq. 0) then
        call derivs(x,y,dydt_dum)
        if (pure_network .eq. 0) then
         ener = y(iener + kk)
         denerdt = dydt_dum(iener + kk)
        else
         ener = 0.0d0
         denerdt = 0.0d0
        end if
       else
        sdot    = 0.0d0
        sneut   = 0.0d0
        ener    = 0.0d0
        denerdt = 0.0d0
       end if


! call an eos
       if (pure_network .eq. 0) then
        temp_row(1) = y(itemp+kk)
        den_row(1)  = y(iden+kk)
       else
        temp_row(1) = btemp
        den_row(1)  = bden
       end if
       if (trho_hist) call update2(x,temp_row(1),den_row(1))
       abar_row(1) = abar
       zbar_row(1) = zbar
       jlo_eos = 1
       jhi_eos = 1

      if (pt_hist) then
       call update3(x,temp_row(1),bpres)
       den_row(1)  = bpres * abar/(avo * kerg * temp_row(1))
       call invert_helm_pt
      else
       call helmeos
!       call eosfxt
      end if


! figure some time scales
       call time_scales(temp_row(1),den_row(1),taud,tau_nse,tau_qse)


! compute the chemical potentials
!       do j=1,ionmax
!        chem_pot(j) = abar*((zion(j) - zbar)*deionz_row(1) &
!                          + (aion(j) - abar)*deiona_row(1))
!       end do
!       sum = 0.0d0
!       do j=1,ionmax
!        sum = sum + chem_pot(j) * dydt_dum(j)
!       end do
!       chem_sum = sum



! and write what we found


! total c12+c12 rate, mass fraction of c12, function
!       zc12 = ratdum(ir1212n) + ratdum(ir1212p) + ratdum(ir1212a)
!       xc12 = y(ic12)*aion(ic12)
!       ff   = sdot/(y(ic12)**2 * zc12) * 2.0d0/3.0d0

       write(22,05) kount,x,temp_row(1),den_row(1), &
                    ener,sdot,sneut,denerdt,yex,xcons


!                    chem_sum,chem_sum/denerdt
!     2              xc12,zc12/den_row(1),ff


       if (iposx .eq. 0) then 
        posx = 0.0d0
       else
        posx = y(iposx+kk)
       end if
       if (ivelx .eq. 0) then 
        velx = 0.0d0
       else
        velx = y(ivelx+kk)
       end if

       write(23,05) kount,x,posx,velx, &
                   temp_row(1),den_row(1),ptot_row(1), &
                   ener,stot_row(1),cs_row(1)



! close up the files
       close(unit=22)
       close(unit=23)

! end of spatial loop
      end do


!      write(6,*) 'done with thermo file'




! for every spatial zone
      do k=1,max(1,nzone)
       kk = neqs*(k-1)

! write out the isotopic mass fractions in blocks of 8
! lop is how many groups of 8 exist; jrem is the remainder
       lop  = ionmax/8
       jrem  = ionmax - 8*lop
       do ilop = 1,lop+1
        kb = 1 + 8*(ilop-1)
        ke = 8 + 8*(ilop-1)
        if (ilop .eq. lop+1  .and. jrem .eq. 0) goto 70
        if (ilop .eq. lop+1) ke = ionmax

! open the output file in append mode (f77) or position mode (f90)
! abundance evolution file
        write(string,03) hfile,ilop+1,'_z',k,'.dat'
        call sqeeze(string)
!        open (unit=34, file=string, status='old', access='append')
        open (unit=34, file=string, status='old', position='append')

        write(34,04) kount,x,(y(nn+kk)*aion(nn), nn=kb,ke)
!        write(34,04) kount,x,(y(nn+kk), nn=kb,ke)

        close(unit=34)
70      continue
       enddo


! end of spatial zone loop
      enddo

!      write(6,*) 'done with mass fractions file'




! start of nse analysis

      if (nse_analysis .eq. 1) then

! for every spatial zone
       do k=1,max(1,nzone)
        kk = neqs*(k-1)


! open the files in append mode (f77) or position mode (f90)

! nse analysis file
       write(string,03) hfile,0,'_z',k,'_nse.dat'
       call sqeeze(string)
!       open (unit=25, file=string, status='old', access='append')
       open (unit=25, file=string, status='old', position='append')


! form the mass fractions
       do j=1,ionmax
        xdum(j) = min(1.0d0,max(y(j+kk)*aion(j),1.0d-30))
       enddo



! normalized mass fractions
       sum = 0.0d0
       do j=1,ionmax
        sum = sum + xdum(j)
       enddo
       xcons = 1.0d0 - sum
       sum = 1.0d0/sum
       do j=1,ionmax
        xdum(j) = min(1.0d0,max(sum*xdum(j),1.0d-30))
       enddo


! get abar, zbar and a few other composition variables
       call azbar(xdum(ionbeg),aion(ionbeg),zion(ionbeg),wion(ionbeg),ionmax, &
                  ydum(ionbeg),abar,zbar,wbar,yex,xcess)



! set the temperature and density
       if (pure_network .eq. 0) then
        temp_row(1) = y(itemp+kk)
        den_row(1)  = y(iden+kk)
       else
        temp_row(1) = btemp
        den_row(1)  = bden
       end if
       if (trho_hist) call update2(x,temp_row(1),den_row(1))
       if (pt_hist) then
        call update3(x,temp_row(1),bpres)
        den_row(1)  = bpres * abar/(avo * kerg * temp_row(1))
        call invert_helm_pt
       end if


! with the temperature, density, and ye
! compute the nse state if the temperature is high enough

       if (temp_row(1) .gt. 2.0e9) then
        igues = 1
        call nse(temp_row(1),den_row(1),yex,igues,1,1,xsum,xmun,xmup,0)
       else
        do j=1,ionmax
         xsum(j) = 1.0e20
        enddo
       end if

! figure delta on the top 20 nse mass fractions
       call indexx(ionmax,xsum(ionbeg),izwork1(ionbeg))
       sum = 0.0d0
       kb  = 0
       do j = ionmax, max(1,ionmax-19), -1
        if (xsum(izwork1(j)) .ge. 1.0e-6) then
         kb = kb + 1
         tdum = (xsum(izwork1(j)) - xdum(izwork1(j)))/xsum(izwork1(j))
!         tdum = (xsum(izwork1(j)) - xdum(izwork1(j)))**2
         sum  = sum + tdum
        end if
       enddo
       sum = sum/float(kb)
!       sum = sqrt(sum/kb)


! figure the time scales
       call time_scales(temp_row(1),den_row(1),taud,tau_nse,tau_qse)



! write out what we got
       write(25,05) kount,x,temp_row(1),den_row(1),yex, &
                   tau_qse,tau_nse,sum,xcons

! close up the files
       close(unit=25)


! write out the isotopic mass fractions in blocks of 8
! lop is how many groups of 8 exist; jrem is the remainder
       lop  = ionmax/8
       jrem  = ionmax - 8*lop
       do ilop = 1,lop+1
        kb = 1 + 8*(ilop-1)
        ke = 8 + 8*(ilop-1)
        if (ilop .eq. lop+1  .and. jrem .eq. 0) goto 80
        if (ilop .eq. lop+1) ke = ionmax


! open the output file in append mode (f77) or position mode (f90)
! abundance evolution file

        write(string,03) hfile,ilop+1,'_z',k,'_nse.dat'
        call sqeeze(string)
!        open (unit=35, file=string, status='old', access='append')
        open (unit=35, file=string, status='old', position='append')

        write(35,04) kount,x,(xsum(nn+kk), nn=kb,ke)

        close(unit=35)
80      continue
       enddo


! end of spatial zone loop
      enddo

! end of the nse analysis if
      end if



      return
      end
!---------------------------------------------------------------------











!---------------------------------------------------------------------
! reaction rate library

! torch rates
! li7(t,n)   a(an,g)    be9(p,d)    be9(p,n)    b10(a,n)   b11(a,n)
! n14(p,a)   c11(p,g)   c12(a,n)    c13(a,n)    c13(p,n)   c14(a,g)
! c14(p,n)   c14(p,g)   o16(p,a)    n14(p,n)    n14(a,n)   n15(p,n)
! n15(a,n)   n15(a,g)   o14(a,g)    o17(a,g)    o17(a,n)   o18(a,g)
! o18(a,n)   ne20(p,a)  f18(p,g)    f19(p,g)    f19(p,n)   f19(a,p)
! na22(n,a)  ne20(p,g)  na23(p,a)   ne20(n,g)   ne21(p,g)  ne21(a,g)
! ne22(p,g)  ne22(a,g)  na22(n,p)   ne22(a,n)   na21(p,g)  mg24(p,a)
! ne21(a,n)  na22(p,g)  na23(p,g)   na23(p,n)   mg24(p,g)  al27(p,a)
! mg25(p,g)  mg25(a,p)  mg25(a,g)   mg25(a,n)   mg26(p,g)  mg26(a,g)
! mg26(a,n)  al25(p,g)  al26(p,g)   al27(a,n)   si27(p,g)  si28(p,g)
! si29(p,g)  si30(p,g)

! bigbang rates:
! n(e-nu)p   p(e-,nu)n  d(p,n)      d(n,g)      d(d,p)     d(d,n)
! t(p,n)     d(d,g)     t(p,g)      t(d,n)      t(t,2n)    he3(d,p)
! he3(t,d)   he3(t,np)  he4(np,g)   he4(d,g)    he4(t,n)   li6(p,he3)
! li6(n,g)   li7(d,n)   lit(t,2n)   li7(he3,np) li6(p,g)   li7(p,n)
! be7(d,p)   be7(t,np)  be7(3he,2p) li6(a,g)    li7(a,n)   be9(p,g)
! b10(p,a)   li7(a,g)   b11(p,a)    be7(a,g)    b11(p,n)   b8(a,p)
! b10(p,g)   c11(n,a)   be9(a,n)    b11(p,g)    b11(a,p)

! pp123 rates:
! p(p,e+nu)  p(n,g)     d(p,g)      he3(n,g)    he3+he3    he3(a,g)
! be7(e-,nu) be7(p,g)   li7(p,g)    li7(p,a)    b8(e+,nu)

! cno rates:
! c12(p,g)   n13(e-nu)  c13(p,g)    n14(p,g)    o15(e-nu)  n14(a,g)
! n15(p,g)   n15(p,a)   o16(p,g)    o17(p,a)    o17(p,g)   o18(p,a)
! o18(p,g)   f17(e-nu)  f18(e-nu)   f19(p,a)

! hot cno rates
! n13(p,g)   o14(e-nu)  o14(a,p)    o15(a,g)    f17(p,g)   ne18(e-nu)
! f18(p,a)   ne18(a,p)  ne19(p,g)   ne19(e-nu)  si26(a,p)

! alfa chain rates:
! a(aa,g)    c12(a,g)   c12+c12     c12+o16     o16+o16    o16(a,g)
! ne20(a,g)  ne20(a,g)  mg24(a,g)   mg24(a,p)   al27(p,g)  si28(a,g)
! si28(a,p)  p31(p,g)   s32(a,g)    s32(a,p)    cl35(p,g)  ar36(a,g)
! ar36(a,p)  k39(p,g)   ca40(a,g)   ca40(a,p)   sc43(p,g)  ti44(a,g)
! ti44(a,p)  v47(p,g)   cr48(a,g)   cr(a,p)     mn51(p,g)  fe52(a,g)
! fe52(a,p)  co55(p,g)

! photodisintegration rates:
! fe52(n,g) fe53(n,g)  fe54(p,g)






      subroutine tfactors(temp)
      include 'implno.dek'
      include 'tfactors.dek'

! sets various popular temperature factors into common block
! this routine must be called before any of the rates are called

! declare the pass
      double precision temp

! all these are in common block

      t9    = temp * 1.0d-9
      t92   = t9*t9
      t93   = t9*t92
      t94   = t9*t93
      t95   = t9*t94
      t96   = t9*t95

      t912  = sqrt(t9)
      t932  = t9*t912
      t952  = t9*t932
      t972  = t9*t952

      t913  = t9**oneth
      t923  = t913*t913
      t943  = t9*t913
      t953  = t9*t923
      t973  = t953*t923
      t9113 = t973*t943

      t914  = t9**(0.25d0)
      t934  = t914*t914*t914
      t954  = t9*t914
      t974  = t9*t934

      t915  = t9**onefif
      t935  = t915*t915*t915
      t945  = t915 * t935
      t965  = t9 * t915

      t916  = t9**onesix
      t976  = t9 * t916
      t9i76 = 1.0d0/t976

      t917  = t9**onesev
      t927  = t917*t917
      t947  = t927*t927

      t918  = sqrt(t914)
      t938  = t918*t918*t918
      t958  = t938*t918*t918

      t9i   = 1.0d0/t9
      t9i2  = t9i*t9i
      t9i3  = t9i2*t9i

      t9i12 = 1.0d0/t912
      t9i32 = t9i*t9i12
      t9i52 = t9i*t9i32
      t9i72 = t9i*t9i52

      t9i13 = 1.0d0/t913
      t9i23 = t9i13*t9i13
      t9i43 = t9i*t9i13
      t9i53 = t9i*t9i23

      t9i14 = 1.0d0/t914
      t9i34 = t9i14*t9i14*t9i14
      t9i54 = t9i*t9i14

      t9i15 = 1.0d0/t915
      t9i35 = t9i15*t9i15*t9i15
      t9i45 = t9i15 * t9i35
      t9i65 = t9i*t9i15

      t9i17 = 1.0d0/t917
      t9i27 = t9i17*t9i17
      t9i47 = t9i27*t9i27

      t9i18 = 1.0d0/t918
      t9i38 = t9i18*t9i18*t9i18
      t9i58 = t9i38*t9i18*t9i18

      return
      end







!---------------------------------------------------------------------







!---------------------------------------------------------------------
      subroutine net_initialize
      include 'implno.dek'
      include 'network.dek'

! initializes quantities

! local variables
      integer   i


! general options
      screen_on      = 1
      use_tables     = 1
      weak_on        = 1
      ffn_on         = 0
      pure_network   = 0
      nse_analysis   = 0
      allow_nse_evol = 0


! printing information
      iprint_files  = 1
      iprint_screen = 1


! inititailize the burn type logicals
      one_step             = .false.
      hydrostatic          = .false.
      expansion            = .false.
      self_heat_const_den  = .false.
      self_heat_const_pres = .false.
      pt_hist              = .false.
      bbang                = .false.
      detonation           = .false.
      trho_hist            = .false.


! adiabatic expanion off
      psi       = 0.0d0
      temp_stop = 1.0d30


! mass fractions above sthreshold are written to the summary file
      sthreshold = 1.0d30

      return
      end
!---------------------------------------------------------------------







!---------------------------------------------------------------------
      subroutine net_summary(tstep,tin,din,ein,tout,dout,eout,conserv, &
                             nbad,nok,xout)
      include 'implno.dek'
      include 'timers.dek'
      include 'vector_eos.dek'
      include 'burn_common.dek'
      include 'network.dek'


! writes out a summary of the network run

! declare the pass
      integer          nbad,nok
      double precision tstep,tin,din,ein,tout,dout,eout,conserv, &
                       xout(*)

! local variables
      character*80     summary
      integer          i,j,k,lenstr,ioff
      double precision abar,zbar,wbar,ye,xcess


! popular format statements
 01   format(a,'summary.dat')
 02   format(1x,a,'=',1pe10.3,' ',a,'=',1pe10.3,' ', &
                a,'=',1pe10.3,' ',a,'=',1pe10.3,' ', &
                a,'=',1pe10.3)
 03   format(1x,a,1pe20.12)
 04   format(1x,a,':',/, &
             1x,3(a,1pe20.12),/, &
             1x,3(a,1pe20.12),/, &
             1x,2(a,1pe11.3),2(a,i5))
 08   format(1x,a,1pe10.3,a)
 09   format(1x,a,i2,a)



! construct the file name and open it
       write(summary,01) hfile(1:lenstr(hfile,80))
       call sqeeze(summary)
       open(unit=41,file=summary,status='unknown')


       write(6,*) ' '
       write(6,04) netname, &
                   ' tin =',tin,' din =',din,' ein =',ein, &
                   ' tout=',tout,' dout=',dout,' eout=',eout, &
                   ' dener=',(eout - ein),' sum =',conserv, &
                   ' nbad=',nbad,' nok=',nok
       write(6,*) ' '

       write(41,*) ' '
       write(41,04) netname, &
                   ' tin =',tin,' din =',din,' ein =',ein, &
                   ' tout=',tout,' dout=',dout,' eout=',eout, &
                   ' enuc=',(eout - ein)/tstep,' sum =',conserv, &
                   ' nbad=',nbad,' nok=',nok
       write(41,*) ' '



! write out the biggest mass fractions
       call indexx(ionmax,xout(ionbeg),izwork1(ionbeg))
       ioff = ionbeg - 1

       if (sthreshold .le. 1  .and. sthreshold .gt. 0.0) then
        do i=ionmax,1,-1
         if (xout(izwork1(i)+ioff) .lt. sthreshold) then
          k = i + 1
          write(6,08)  'mass fractions larger than ',sthreshold
          write(41,08) 'mass fractions larger than ',sthreshold
          goto 20
         end if
        end do
       else
        j = min(20,ionmax)
        k = max(ionmax-19,ionbeg)
        write(6,09)  'top ',j,' mass fractions:'
        write(41,09) 'top ',j,' mass fractions:'
       end if

 20   continue


       write(6,02) (ionam(izwork1(i)+ioff),xout(izwork1(i)+ioff), i=ionmax,k,-1)
       if (iprot .ne. 0 .and. ineut .ne. 0) then
        write(6,02) ionam(iprot),xout(iprot), &
                    ionam(ineut),xout(ineut), &
                    ionam(ihe4),xout(ihe4)
       end if
       write(6,*) ' '

       write(41,02) (ionam(izwork1(i)+ioff),xout(izwork1(i)+ioff), i=ionmax,k,-1)
       if (iprot .ne. 0 .and. ineut .ne. 0) then
        write(41,02) ionam(iprot),xout(iprot), &
                     ionam(ineut),xout(ineut), &
                     ionam(ihe4),xout(ihe4)
       end if
       write(41,*) ' '



! end the clock
      call zsecond(timtot)
      timtot = timtot - timzer
      call timlap(timtot,hours,minuts,secs,msecs)
      write(6,100) hours,minuts,secs,msecs
      write(41,100) hours,minuts,secs,msecs
 100  format(1x,'cpu time : ',i2.2,' hrs  ',i2.2,' min  ', &
                              i2.2,' sec  ',i6,' usec',/,/)


! close up shop
      close(unit=41)
      return
      end
!---------------------------------------------------------------------







!---------------------------------------------------------------------
! this file contains routines that sort, search and select parts of arrays:
!
! index and rank makers:
! routine indexx constructs a sort index for a real array



      subroutine indexx(n,arr,indx)
      include 'implno.dek'
!
! indexes an array arr(1:n). that is it outputs the array indx(1:n) such
! that arr(indx(j)) is in ascending order for j=1...n. the input quantities
! are not changed.
!
! declare
      integer          n,indx(n),m,nstack
      parameter        (m=7, nstack = 50)
      integer          i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
      double precision arr(n),a
!
! initialize
      do 11 j=1,n
       indx(j) = j
11    continue
      jstack = 0
      l      = 1
      ir     = n
!
! insertion sort when subbarray small enough
1     if (ir - l .lt. m) then
       do 13 j=l+1,ir
        indxt = indx(j)
        a     = arr(indxt)
        do 12 i=j-1,l,-1
         if (arr(indx(i)) .le. a) go to 2
         indx(i+1) = indx(i)
12      continue
        i = l - 1
2       indx(i+1) = indxt
13     continue
!
! pop stack and begin a new round of partitioning
       if (jstack .eq. 0) return
       ir     = istack(jstack)
       l      = istack(jstack-1)
       jstack = jstack - 2
!
! choose median of left, center and right elements as partitioning element
! also rearrange so that a(l+1) < a(l) < a(ir)
      else
       k         = (l + ir)/2
       itemp     = indx(k)
       indx(k)   = indx(l+1)
       indx(l+1) = itemp

       if (arr(indx(l)) .gt. arr(indx(ir))) then
        itemp    = indx(l)
        indx(l)  = indx(ir)
        indx(ir) = itemp
       end if


       if(arr(indx(l+1)).gt.arr(indx(ir)))then
        itemp=indx(l+1)
        indx(l+1)=indx(ir)
        indx(ir)=itemp
       endif
       if(arr(indx(l)).gt.arr(indx(l+1)))then
        itemp=indx(l)
        indx(l)=indx(l+1)
        indx(l+1)=itemp
       endif

!
! initialize pointers for partitioning
       i     = l + 1
       j     = ir
       indxt = indx(l+1)
       a     = arr(indxt)
3      continue
       i = i + 1
       if (arr(indx(i)) .lt. a) go to 3
4      continue
       j = j - 1
       if (arr(indx(j)) .gt. a) go to 4
       if (j .lt. i) go to 5
       itemp   = indx(i)
       indx(i) = indx(j)
       indx(j) = itemp
       go to 3
!
5      indx(l+1) = indx(j)
       indx(j)   = indxt
       jstack    = jstack + 2
!
! push pointers to larger subarray on stack
       if (jstack .gt. nstack) stop 'jstack > nstack in routine indexx'
       if (ir - i + 1  .ge.  j - l) then
        istack(jstack)   = ir
        istack(jstack-1) = i
        ir               = j - 1
       else
        istack(jstack)   = j-1
        istack(jstack-1) = l
        l                = i
       end if
      end if
      go to 1
      end
!---------------------------------------------------------------------








!---------------------------------------------------------------------

      subroutine azbar(xmass,aion,zion,wion,ionmax, &
                       ymass,abar,zbar,wbar,ye,nxcess)
      include 'implno.dek'

! this routine calculates composition variables

! input:
! mass fractions               = xmass(1:ionmax)  dimensionless
! number of nucleons           = aion(1:ionmax)   dimensionless
! charge of nucleus            = zion(1:ionmax)   dimensionless
! atomic weight or molar mass  = wion(1:ionmax)    g/mole
! number of isotopes           = ionmax
!
! output:
! molar abundances        = ymass(1:ionmax)   mole/g
! mean number of nucleons = abar              dimensionless
! mean nucleon charge     = zbar              dimensionless
! mean weight             = wbar              g/mole
! electron fraction       = ye                mole/g
! neutron excess          = xcess


! declare the pass
      integer          ionmax
      double precision xmass(ionmax),aion(ionmax),zion(ionmax), &
                       wion(ionmax),ymass(ionmax),abar,zbar,wbar, &
                       ye,nxcess


! local variables
      double precision asum,sum1


! molar abundances
      ymass(1:ionmax) = xmass(1:ionmax)/wion(1:ionmax)

! mean molar mass
      wbar  = 1.0d0/sum(ymass(1:ionmax))

! mean number of nucleons
      sum1  = sum(aion(1:ionmax)*ymass(1:ionmax))
      abar  = wbar * sum1

! mean charge
      ye  = sum(zion(1:ionmax)*ymass(1:ionmax))
      zbar  = wbar * ye

! neutron excess
      nxcess = sum1 - 2.0d0 * ye

      return
      end
!---------------------------------------------------------------------




