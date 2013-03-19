module rhs_module

  use network_indices
  implicit none

contains

  subroutine aprox13(tt,temp,dens,y,dydt,enucdot)

    ! this routine sets up the system of ode's for the aprox13
    ! nuclear reactions.  this is an alpha chain + heavy ion network
    ! with (a,p)(p,g) links
    !     
    ! isotopes: he4,  c12,  o16,  ne20, mg24, si28, s32,
    !           ar36, ca40, ti44, cr48, fe52, ni56
    
    

    ! declare the pass
    double precision :: tt,temp,dens,y(nspec)
    double precision :: dydt(nspec)
    

    ! local variables
    logical          deriva
    parameter        (deriva = .false.)
    
    double precision :: abar, zbar, ye
    
    
    ! positive definite mass fractions
    do i=ionbeg,ionend
       y(i) = min(1.0d0,max(y(i),smallx))
    enddo


    ! generate abar and zbar for this composition
    abar = 1.0d0/sum(y(ionbeg:ionend))
    zbar = sum(zion(ionbeg:ionend)*y(ionbeg:ionend)) * abar
    ye    = zbar/abar
    

    ! get the raw reaction rates
    call aprox13rat(temp, dens, ratraw, dratrawdt, dratrawdd)
    

    ! do the screening here because the corrections depend on the composition
    call screen_aprox13(temp, dens, y, &
                        ratraw, dratrawdt, dratrawdd, &
                        ratdum, dratdumdt, dratdumdd, &
                        scfac, dscfacdt, dscfacdd)
    

    ! get the right hand side of the odes
    call rhs(y,ratdum,dydt,deriva)
    

    ! instantaneous energy generation rate
    call ener_gener_rate(dydt,enuc)
    
    ! get the neutrino losses
    call sneut5(temp,dens,abar,zbar, &
                sneut,dsneutdt,dsneutdd,snuda,snudz)


    ! append an energy equation
    enucdot = enuc - sneut

    return
  end subroutine aprox13


  subroutine rhs(y,rate,dydt,deriva)

    ! evaluates the right hand side of the aprox13 odes

    ! deriva is used in forming the analytic Jacobian to get
    ! the derivative wrt A

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
  end subroutine rhs


  subroutine aprox13rat(btemp, bden, ratraw, dratrawdt, dratrawdd)

    ! this routine generates unscreened
    ! nuclear reaction rates for the aprox13 network.

    ! declare
    integer          i
    double precision rrate,drratedt,drratedd


    ! zero the rates
    do i=1,nrat
       ratraw(i) = 0.0d0
       dratrawdt(i) = 0.0d0
       dratrawdd(i) = 0.0d0
    enddo
  
    if (btemp .lt. 1.0d6) return


    ! get the temperature factors
    tf = get_tfactors(btemp)


    ! c12(a,g)o16
    call rate_c12ag(tf,bden, &
         ratraw(ircag),dratrawdt(ircag),dratrawdd(ircag), &
         ratraw(iroga),dratrawdt(iroga),dratrawdd(iroga))
    
    ! triple alpha to c12
    call rate_tripalf(tf,bden, &
         ratraw(ir3a),dratrawdt(ir3a),dratrawdd(ir3a), &
         ratraw(irg3a),dratrawdt(irg3a),dratrawdd(irg3a))
    
    ! c12 + c12
    call rate_c12c12(tf,bden, &
         ratraw(ir1212),dratrawdt(ir1212),dratrawdd(ir1212), &
         rrate,drratedt,drratedd)
    
    ! c12 + o16
    call rate_c12o16(tf,bden, &
         ratraw(ir1216),dratrawdt(ir1216),dratrawdd(ir1216), &
         rrate,drratedt,drratedd)
    
    ! o16 + o16
    call rate_o16o16(tf,bden, &
         ratraw(ir1616),dratrawdt(ir1616),dratrawdd(ir1616), &
         rrate,drratedt,drratedd)
    
    ! o16(a,g)ne20
    call rate_o16ag(tf,bden, &
         ratraw(iroag),dratrawdt(iroag),dratrawdd(iroag), &
         ratraw(irnega),dratrawdt(irnega),dratrawdd(irnega))
    
    ! ne20(a,g)mg24
    call rate_ne20ag(tf,bden, &
         ratraw(irneag),dratrawdt(irneag),dratrawdd(irneag), &
         ratraw(irmgga),dratrawdt(irmgga),dratrawdd(irmgga))
    
    ! mg24(a,g)si28
    call rate_mg24ag(tf,bden, &
         ratraw(irmgag),dratrawdt(irmgag),dratrawdd(irmgag), &
         ratraw(irsiga),dratrawdt(irsiga),dratrawdd(irsiga))
    
    ! mg24(a,p)al27
    call rate_mg24ap(tf,bden, &
         ratraw(irmgap),dratrawdt(irmgap),dratrawdd(irmgap), &
         ratraw(iralpa),dratrawdt(iralpa),dratrawdd(iralpa))
    
    ! al27(p,g)si28
    call rate_al27pg(tf,bden, &
         ratraw(iralpg),dratrawdt(iralpg),dratrawdd(iralpg), &
         ratraw(irsigp),dratrawdt(irsigp),dratrawdd(irsigp))
    
    ! si28(a,g)s32
    call rate_si28ag(tf,bden, &
         ratraw(irsiag),dratrawdt(irsiag),dratrawdd(irsiag), &
         ratraw(irsga),dratrawdt(irsga),dratrawdd(irsga))
    
    ! si28(a,p)p31
    call rate_si28ap(tf,bden, &
         ratraw(irsiap),dratrawdt(irsiap),dratrawdd(irsiap), &
         ratraw(irppa),dratrawdt(irppa),dratrawdd(irppa))
    
    ! p31(p,g)s32
    call rate_p31pg(tf,bden, &
         ratraw(irppg),dratrawdt(irppg),dratrawdd(irppg), &
         ratraw(irsgp),dratrawdt(irsgp),dratrawdd(irsgp))
    
    ! s32(a,g)ar36
    call rate_s32ag(tf,bden, &
         ratraw(irsag),dratrawdt(irsag),dratrawdd(irsag), &
         ratraw(irarga),dratrawdt(irarga),dratrawdd(irarga))
    
    ! s32(a,p)cl35
    call rate_s32ap(tf,bden, &
         ratraw(irsap),dratrawdt(irsap),dratrawdd(irsap), &
         ratraw(irclpa),dratrawdt(irclpa),dratrawdd(irclpa))
    
    ! cl35(p,g)ar36
    call rate_cl35pg(tf,bden, &
         ratraw(irclpg),dratrawdt(irclpg),dratrawdd(irclpg), &
         ratraw(irargp),dratrawdt(irargp),dratrawdd(irargp))
    
    ! ar36(a,g)ca40
    call rate_ar36ag(tf,bden, &
         ratraw(irarag),dratrawdt(irarag),dratrawdd(irarag), &
         ratraw(ircaga),dratrawdt(ircaga),dratrawdd(ircaga))
    
    ! ar36(a,p)k39
    call rate_ar36ap(tf,bden, &
         ratraw(irarap),dratrawdt(irarap),dratrawdd(irarap), &
         ratraw(irkpa),dratrawdt(irkpa),dratrawdd(irkpa))
    
    ! k39(p,g)ca40
    call rate_k39pg(tf,bden, &
         ratraw(irkpg),dratrawdt(irkpg),dratrawdd(irkpg), &
         ratraw(ircagp),dratrawdt(ircagp),dratrawdd(ircagp))
    
    ! ca40(a,g)ti44
    call rate_ca40ag(tf,bden, &
         ratraw(ircaag),dratrawdt(ircaag),dratrawdd(ircaag), &
         ratraw(irtiga),dratrawdt(irtiga),dratrawdd(irtiga))
    
    ! ca40(a,p)sc43
    call rate_ca40ap(tf,bden, &
         ratraw(ircaap),dratrawdt(ircaap),dratrawdd(ircaap), &
         ratraw(irscpa),dratrawdt(irscpa),dratrawdd(irscpa))
    
    ! sc43(p,g)ti44
    call rate_sc43pg(tf,bden, &
         ratraw(irscpg),dratrawdt(irscpg),dratrawdd(irscpg), &
         ratraw(irtigp),dratrawdt(irtigp),dratrawdd(irtigp))
    
    ! ti44(a,g)cr48
    call rate_ti44ag(tf,bden, &
         ratraw(irtiag),dratrawdt(irtiag),dratrawdd(irtiag), &
         ratraw(ircrga),dratrawdt(ircrga),dratrawdd(ircrga))

    ! ti44(a,p)v47
    call rate_ti44ap(tf,bden, &
         ratraw(irtiap),dratrawdt(irtiap),dratrawdd(irtiap), &
         ratraw(irvpa),dratrawdt(irvpa),dratrawdd(irvpa))
    
    ! v47(p,g)cr48
    call rate_v47pg(tf,bden, &
         ratraw(irvpg),dratrawdt(irvpg),dratrawdd(irvpg), &
         ratraw(ircrgp),dratrawdt(ircrgp),dratrawdd(ircrgp))
    
    ! cr48(a,g)fe52
    call rate_cr48ag(tf,bden, &
         ratraw(ircrag),dratrawdt(ircrag),dratrawdd(ircrag), &
         ratraw(irfega),dratrawdt(irfega),dratrawdd(irfega))
    
    ! cr48(a,p)mn51
    call rate_cr48ap(tf,bden, &
         ratraw(ircrap),dratrawdt(ircrap),dratrawdd(ircrap), &
         ratraw(irmnpa),dratrawdt(irmnpa),dratrawdd(irmnpa))
    
    ! mn51(p,g)fe52
    call rate_mn51pg(tf,bden, &
         ratraw(irmnpg),dratrawdt(irmnpg),dratrawdd(irmnpg), &
         ratraw(irfegp),dratrawdt(irfegp),dratrawdd(irfegp))
    
    ! fe52(a,g)ni56
    call rate_fe52ag(tf,bden, &
         ratraw(irfeag),dratrawdt(irfeag),dratrawdd(irfeag), &
         ratraw(irniga),dratrawdt(irniga),dratrawdd(irniga))
    
    ! fe52(a,p)co55
    call rate_fe52ap(tf,bden, &
         ratraw(irfeap),dratrawdt(irfeap),dratrawdd(irfeap), &
         ratraw(ircopa),dratrawdt(ircopa),dratrawdd(ircopa))
    
    ! co55(p,g)ni56
    call rate_co55pg(tf,bden, &
         ratraw(ircopg),dratrawdt(ircopg),dratrawdd(ircopg), &
         ratraw(irnigp),dratrawdt(irnigp),dratrawdd(irnigp))
    
    
    return
  end subroutine aprox13rat


  subroutine screen_aprox13(btemp, bden, y, &
                            ratraw, dratrawdt, dratrawdd, &
                            ratdum, dratdumdt, dratdumdd, &
                            scfac, dscfacdt, dscfacdd)
    
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
       dscfacdd(i)  = 0.0d0
    end do


    ! always screen
    
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
    
    return
  end subroutine screen_aprox13



  subroutine ener_gener_rate(dydt,enuc)
    include 'implno.dek'
    include 'const.dek'
    include 'network.dek'

    ! computes the instantaneous energy generation rate

    ! declare the pass
    double precision dydt(1),enuc
    
    ! local variables
    integer          i
    
    ! conversion factors for the nuclear energy generation rate detlap
    ! is the mass excess of the proton in amu detlan is the mass excess
    ! of the neutron in amu
    
    double precision enuc_conv,enuc_conv2,deltap,deltan
    parameter        (enuc_conv  = ev2erg*1.0d6*avo, &
                      enuc_conv2 = -avo*clight*clight, &
                      deltap     = 7.288969d0, &
                      deltan     = 8.071323d0)
    ! instantaneous energy generation rate

    ! this form misses n <-> p differences 
    
    ! enuc = 0.0d0 
    ! do i=1,ionmax 
    !   enuc = enuc + dydt(i) * bion(i) 
    ! enddo 
    ! enuc = enuc * enuc_conv


    ! this form gets the n <-> p differences 
    
    ! enuc = 0.0d0 
    ! do i=1,ionmax
    !      enuc = enuc + dydt(i) * (bion(i) - zion(i)*deltap - nion(i)*deltan)
    ! enddo 
    ! enuc = enuc * enuc_conv
    
    ! this form is closest to e = m c**2 and gives the same results as
    ! the form above
    
    enuc = 0.0d0
    do i=1,ionmax
       enuc = enuc + dydt(i) * mion(i)
    enddo
    enuc = enuc * enuc_conv2
    
    return
  end subroutine ener_gener_rate

end module rhs_module


