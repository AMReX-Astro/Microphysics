module actual_rhs_module

  use network
  use eos_type_module
  use burn_type_module
  use temperature_integration_module, only: temperature_rhs, temperature_jac
  use sneut_module, only: sneut5
  use rate_type_module

  implicit none

  double precision, parameter :: c54 = 56.0d0/54.0d0

contains

  subroutine actual_rhs_init()

    use aprox_rates_module, only: rates_init
    use screening_module, only: screening_init

    implicit none

    call rates_init()

    call set_up_screening_factors()

    call screening_init()

  end subroutine actual_rhs_init



  subroutine actual_rhs(state)

    implicit none

    ! This routine sets up the system of ode's for the aprox19
    ! nuclear reactions.  This is an alpha chain + heavy ion network
    ! with (a,p)(p,g) links, as well as Fe54 for photodisintegration
    ! and hydrogen and nitrogen for PP and CNO burning.
    !
    ! Isotopes: h1,   he3,  he4,  c12,  n14,  o16,  ne20, mg24, si28, s32,
    !           ar36, ca40, ti44, cr48, fe52, fe54, ni56, neut, prot

    type (burn_t)    :: state
    type (rate_t)    :: rr

    logical          :: deriva = .false.

    double precision :: sneut, dsneutdt, dsneutdd, snuda, snudz
    double precision :: enuc

    double precision :: rho, temp, abar, zbar

    double precision :: y(nspec)

    call evaluate_rates(state, rr)

    ! Get the data from the state

    rho  = state % rho
    temp = state % T
    abar = state % abar
    zbar = state % zbar
    y    = state % xn * aion_inv

    ! Call the RHS to actually get dydt.

    call rhs(y, rr % rates(1,:), rr % rates(1,:), state % ydot(1:nspec), deriva)

    ! Instantaneous energy generation rate -- this needs molar fractions

    call ener_gener_rate(state % ydot(1:nspec), enuc)

    ! Get the neutrino losses

    call sneut5(temp, rho, abar, zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)

    ! Append the energy equation (this is erg/g/s)

    state % ydot(net_ienuc) = enuc - sneut

    ! Append the temperature equation

    call temperature_rhs(state)

  end subroutine actual_rhs



  ! Analytical Jacobian

  subroutine actual_jac(state)

    use bl_types
    use bl_constants_module, only: ZERO
    use eos_module

    implicit none

    type (burn_t)    :: state
    type (rate_t)    :: rr

    logical          :: deriva = .true.

    double precision :: b1, sneut, dsneutdt, dsneutdd, snuda, snudz

    integer          :: j

    double precision :: rho, temp, abar, zbar
    double precision :: y(nspec)

    state % jac(:,:) = ZERO

    call evaluate_rates(state, rr)

    ! Get the data from the state

    rho  = state % rho
    temp = state % T
    abar = state % abar
    zbar = state % zbar
    y    = state % xn * aion_inv

    ! Species Jacobian elements with respect to other species

    call dfdy_isotopes_aprox19(y, state % jac(1:nspec,1:nspec), rr % rates(1,:), &
                               rr % rates(3,:), rr % rates(4,:))

    ! Energy generation rate Jacobian elements with respect to species

    do j = 1, nspec
       call ener_gener_rate(state % jac(1:nspec,j), state % jac(net_ienuc,j))
    enddo

    ! Account for the thermal neutrino losses

    call sneut5(temp, rho, abar, zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)

    do j = 1, nspec
       b1 = ((aion(j) - abar) * abar * snuda + (zion(j) - zbar) * abar * snudz)
       state % jac(net_ienuc,j) = state % jac(net_ienuc,j) - b1
    enddo

    ! Evaluate the Jacobian elements with respect to temperature by
    ! calling the RHS using d(ratdum) / dT

    call rhs(y, rr % rates(2,:), rr % rates(1,:), state % jac(1:nspec,net_itemp), deriva)

    call ener_gener_rate(state % jac(1:nspec,net_itemp), state % jac(net_ienuc,net_itemp))
    state % jac(net_ienuc,net_itemp) = state % jac(net_ienuc,net_itemp) - dsneutdt

    ! Temperature Jacobian elements

    call temperature_jac(state)

  end subroutine actual_jac



  subroutine evaluate_rates(state, rr)

    implicit none

    type (burn_t), intent(in)    :: state
    type (rate_t), intent(out)   :: rr

    double precision :: ratraw(nrates), dratrawdt(nrates), dratrawdd(nrates)
    double precision :: ratdum(nrates), dratdumdt(nrates), dratdumdd(nrates)
    double precision :: dratdumdy1(nrates), dratdumdy2(nrates)
    double precision :: scfac(nrates),  dscfacdt(nrates),  dscfacdd(nrates)

    double precision :: rho, temp, abar, zbar

    double precision :: y(nspec)

    ! Get the data from the state

    rho  = state % rho
    temp = state % T
    abar = state % abar
    zbar = state % zbar
    y    = state % xn * aion_inv

    ! Get the raw reaction rates
    call aprox19rat(temp, rho, ratraw, dratrawdt, dratrawdd)

    ! Weak screening rates
    call weak_aprox19(y, state, ratraw, dratrawdt, dratrawdd)

    ! Do the screening here because the corrections depend on the composition
    call screen_aprox19(temp, rho, y,                 &
                        ratraw, dratrawdt, dratrawdd, &
                        ratdum, dratdumdt, dratdumdd, &
                        dratdumdy1, dratdumdy2,       &
                        scfac,  dscfacdt,  dscfacdd)

    ! Save the rate data, for the Jacobian later if we need it.

    rr % rates(1,:) = ratdum
    rr % rates(2,:) = dratdumdt
    rr % rates(3,:) = dratdumdy1
    rr % rates(4,:) = dratdumdy2

    rr % T_eval = temp

  end subroutine evaluate_rates


  ! Evaluates the right hand side of the aprox19 ODEs

  subroutine rhs(y, rate, ratdum, dydt, deriva)

    use bl_constants_module, only: ZERO, SIXTH
    use microphysics_math_module, only: esum

    implicit none

    ! deriva is used in forming the analytic Jacobian to get
    ! the derivative wrt A

    logical          :: deriva
    double precision :: y(nspec),rate(nrates),ratdum(nrates),dydt(nspec)

    ! local variables

    double precision :: a(20)

    dydt(1:nspec) = ZERO

    ! hydrogen reactions
    a(1) = -1.5d0 * y(ih1) * y(ih1) * rate(irpp)
    a(2) =  y(ihe3) * y(ihe3) * rate(ir33)
    a(3) = -y(ihe3) * y(ihe4) * rate(irhe3ag)
    a(4) = -2.0d0 * y(ic12) * y(ih1) * rate(ircpg)
    a(5) = -2.0d0 * y(in14) * y(ih1) * rate(irnpg)
    a(6) = -2.0d0 * y(io16) * y(ih1) * rate(iropg)
    a(7) = -3.0d0 * y(ih1) * rate(irpen)

    dydt(ih1) = dydt(ih1) + esum(a,7)



    ! he3 reactions
    a(1)  =  0.5d0 * y(ih1) * y(ih1) * rate(irpp)
    a(2)  = -y(ihe3) * y(ihe3) * rate(ir33)
    a(3)  = -y(ihe3) * y(ihe4) * rate(irhe3ag)
    a(4)  =  y(ih1) * rate(irpen)

    dydt(ihe3) = dydt(ihe3) + esum(a,4)


    ! he4 reactions
    ! heavy ion reactions
    a(1)  = 0.5d0 * y(ic12) * y(ic12) * rate(ir1212)
    a(2)  = 0.5d0 * y(ic12) * y(io16) * rate(ir1216)
    a(3)  = 0.56d0 * 0.5d0 * y(io16) * y(io16) * rate(ir1616)

    dydt(ihe4) =  dydt(ihe4) + esum(a,3)


    ! (a,g) and (g,a) reactions
    a(1)  = -0.5d0 * y(ihe4) * y(ihe4) * y(ihe4) * rate(ir3a)
    a(2)  =  3.0d0 * y(ic12) * rate(irg3a)
    a(3)  = -y(ihe4)  * y(ic12) * rate(ircag)
    a(4)  =  y(io16)  * rate(iroga)
    a(5)  = -y(ihe4)  * y(io16) * rate(iroag)
    a(6)  =  y(ine20) * rate(irnega)
    a(7)  = -y(ihe4)  * y(ine20) * rate(irneag)
    a(8)  =  y(img24) * rate(irmgga)
    a(9)  = -y(ihe4)  * y(img24)* rate(irmgag)
    a(10) =  y(isi28) * rate(irsiga)
    a(11) = -y(ihe4)  * y(isi28)*rate(irsiag)
    a(12) =  y(is32)  * rate(irsga)

    dydt(ihe4) =  dydt(ihe4) + esum(a,12)

    a(1)  = -y(ihe4)  * y(is32) * rate(irsag)
    a(2)  =  y(iar36) * rate(irarga)
    a(3)  = -y(ihe4)  * y(iar36)*rate(irarag)
    a(4)  =  y(ica40) * rate(ircaga)
    a(5)  = -y(ihe4)  * y(ica40)*rate(ircaag)
    a(6)  =  y(iti44) * rate(irtiga)
    a(7)  = -y(ihe4)  * y(iti44)*rate(irtiag)
    a(8)  =  y(icr48) * rate(ircrga)
    a(9)  = -y(ihe4)  * y(icr48)*rate(ircrag)
    a(10) =  y(ife52) * rate(irfega)
    a(11) = -y(ihe4)  * y(ife52) * rate(irfeag)
    a(12) =  y(ini56) * rate(irniga)

    dydt(ihe4) =  dydt(ihe4) + esum(a,12)


    ! (a,p)(p,g) and (g,p)(p,a) reactions

    if (.not.deriva) then
       a(1)  =  0.34d0*0.5d0*y(io16)*y(io16)*rate(irs1)*rate(ir1616)
       a(2)  = -y(ihe4)  * y(img24) * rate(irmgap)*(1.0d0-rate(irr1))
       a(3)  =  y(isi28) * rate(irsigp) * rate(irr1)
       a(4)  = -y(ihe4)  * y(isi28) * rate(irsiap)*(1.0d0-rate(irs1))
       a(5)  =  y(is32)  * rate(irsgp) * rate(irs1)
       a(6)  = -y(ihe4)  * y(is32) * rate(irsap)*(1.0d0-rate(irt1))
       a(7)  =  y(iar36) * rate(irargp) * rate(irt1)
       a(8)  = -y(ihe4)  * y(iar36) * rate(irarap)*(1.0d0-rate(iru1))
       a(9)  =  y(ica40) * rate(ircagp) * rate(iru1)
       a(10) = -y(ihe4)  * y(ica40) * rate(ircaap)*(1.0d0-rate(irv1))
       a(11) =  y(iti44) * rate(irtigp) * rate(irv1)
       a(12) = -y(ihe4)  * y(iti44) * rate(irtiap)*(1.0d0-rate(irw1))
       a(13) =  y(icr48) * rate(ircrgp) * rate(irw1)
       a(14) = -y(ihe4)  * y(icr48) * rate(ircrap)*(1.0d0-rate(irx1))
       a(15) =  y(ife52) * rate(irfegp) * rate(irx1)

       dydt(ihe4) =  dydt(ihe4) + esum(a,15)

    else
       a(1)  =  0.34d0*0.5d0*y(io16)*y(io16) * ratdum(irs1) * rate(ir1616)
       a(2)  =  0.34d0*0.5d0*y(io16)*y(io16) * rate(irs1) * ratdum(ir1616)
       a(3)  = -y(ihe4)*y(img24) * rate(irmgap)*(1.0d0 - ratdum(irr1))
       a(4)  =  y(ihe4)*y(img24) * ratdum(irmgap)*rate(irr1)
       a(5)  =  y(isi28) * ratdum(irsigp) * rate(irr1)
       a(6)  =  y(isi28) * rate(irsigp) * ratdum(irr1)
       a(7)  = -y(ihe4)*y(isi28) * rate(irsiap)*(1.0d0 - ratdum(irs1))
       a(8)  =  y(ihe4)*y(isi28) * ratdum(irsiap) * rate(irs1)
       a(9)  =  y(is32)  * ratdum(irsgp) * rate(irs1)
       a(10) =  y(is32)  * rate(irsgp) * ratdum(irs1)

       dydt(ihe4) =  dydt(ihe4) + esum(a,10)

       a(1)  = -y(ihe4)*y(is32) * rate(irsap)*(1.0d0 - ratdum(irt1))
       a(2)  =  y(ihe4)*y(is32) * ratdum(irsap)*rate(irt1)
       a(3)  =  y(iar36) * ratdum(irargp) * rate(irt1)
       a(4)  =  y(iar36) * rate(irargp) * ratdum(irt1)
       a(5)  = -y(ihe4)*y(iar36) * rate(irarap)*(1.0d0 - ratdum(iru1))
       a(6)  =  y(ihe4)*y(iar36) * ratdum(irarap)*rate(iru1)
       a(7)  =  y(ica40) * ratdum(ircagp) * rate(iru1)
       a(8)  =  y(ica40) * rate(ircagp) * ratdum(iru1)
       a(9)  = -y(ihe4)*y(ica40) * rate(ircaap)*(1.0d0-ratdum (irv1))
       a(10) =  y(ihe4)*y(ica40) * ratdum(ircaap)*rate(irv1)
       a(11) =  y(iti44) * ratdum(irtigp) * rate(irv1)
       a(12) =  y(iti44) * rate(irtigp) * ratdum(irv1)

       dydt(ihe4) =  dydt(ihe4) + esum(a,12)

       a(1)  = -y(ihe4)*y(iti44) * rate(irtiap)*(1.0d0 - ratdum(irw1))
       a(2)  =  y(ihe4)*y(iti44) * ratdum(irtiap)*rate(irw1)
       a(3)  =  y(icr48) * ratdum(ircrgp) * rate(irw1)
       a(4)  =  y(icr48) * rate(ircrgp) * ratdum(irw1)
       a(5)  = -y(ihe4)*y(icr48) * rate(ircrap)*(1.0d0 - ratdum(irx1))
       a(6)  =  y(ihe4)*y(icr48) * ratdum(ircrap)*rate(irx1)
       a(7)  =  y(ife52) * ratdum(irfegp) * rate(irx1)
       a(8)  =  y(ife52) * rate(irfegp) * ratdum(irx1)

       dydt(ihe4) =  dydt(ihe4) + esum(a,8)
    end if


    ! photodisintegration reactions
    a(1) =  y(ife54) * y(iprot) * y(iprot) * rate(ir5f54)
    a(2) = -y(ife52) * y(ihe4) * rate(ir6f54)
    a(3) = -y(ife52) * y(ihe4) * y(iprot) * rate(ir7f54)
    a(4) =  y(ini56) * y(iprot) * rate(ir8f54)
    a(5) = -y(ihe4) * rate(iralf1)
    a(6) =  y(ineut)*y(ineut) * y(iprot)*y(iprot) * rate(iralf2)

    dydt(ihe4) =  dydt(ihe4) + esum(a,6)


    ! ppchain
    a(1) = 0.5d0 * y(ihe3) * y(ihe3) * rate(ir33)
    a(2) = y(ihe3) * y(ihe4) * rate(irhe3ag)

    dydt(ihe4) =  dydt(ihe4) + sum(a(1:2))


    ! cno cycles
    a(1) = y(io16) * y(ih1) * rate(iropg)
    a(2) = -y(ihe4) * y(in14) * rate(irnag) * 1.5d0

    dydt(ihe4) =  dydt(ihe4) + sum(a(1:2))

    if (.not. deriva) then
       a(1) = y(in14) * y(ih1) * rate(ifa) * rate(irnpg)
       dydt(ihe4) =  dydt(ihe4) + a(1)
    else
       a(1) = y(in14) * y(ih1) * rate(ifa) * ratdum(irnpg)
       a(2) = y(in14) * y(ih1) * ratdum(ifa) * rate(irnpg)
       dydt(ihe4) =  dydt(ihe4) + sum(a(1:2))
    end if


    ! c12 reactions
    a(1) = -y(ic12) * y(ic12) * rate(ir1212)
    a(2) = -y(ic12) * y(io16) * rate(ir1216)
    a(3) =  SIXTH * y(ihe4) * y(ihe4) * y(ihe4) * rate(ir3a)
    a(4) = -y(ic12) * rate(irg3a)
    a(5) = -y(ic12) * y(ihe4) * rate(ircag)
    a(6) =  y(io16) * rate(iroga)
    a(7) = -y(ic12) * y(ih1) * rate(ircpg)

    dydt(ic12) =  dydt(ic12) + esum(a,7)

    if (.not. deriva) then
       a(1) =  y(in14) * y(ih1) * rate(ifa) * rate(irnpg)
       dydt(ic12) =  dydt(ic12) + a(1)
    else
       a(1) =  y(in14) * y(ih1) * rate(ifa) * ratdum(irnpg)
       a(2) =  y(in14) * y(ih1) * ratdum(ifa) * rate(irnpg)
       dydt(ic12) =  dydt(ic12) + sum(a(1:2))
    end if


    ! n14 reactions
    a(1) =  y(ic12) * y(ih1) * rate(ircpg)
    a(2) = -y(in14) * y(ih1) * rate(irnpg)
    a(3) =  y(io16) * y(ih1) * rate(iropg)
    a(4) = -y(ihe4) * y(in14) * rate(irnag)

    dydt(in14) =  dydt(in14) + esum(a,4)


    ! o16 reactions
    a(1) = -y(ic12) * y(io16) * rate(ir1216)
    a(2) = -y(io16) * y(io16) * rate(ir1616)
    a(3) =  y(ic12) * y(ihe4) * rate(ircag)
    a(4) = -y(io16) * y(ihe4) * rate(iroag)
    a(5) = -y(io16) * rate(iroga)
    a(6) =  y(ine20) * rate(irnega)
    a(7) = -y(io16) * y(ih1) * rate(iropg)

    dydt(io16) =  dydt(io16) + esum(a,7)

    if (.not. deriva) then
       a(1) =  y(in14) * y(ih1) * rate(ifg) * rate(irnpg)
       dydt(io16) =  dydt(io16) + a(1)
    else
       a(1) =  y(in14) * y(ih1) * rate(ifg) * ratdum(irnpg)
       a(2) =  y(in14) * y(ih1) * ratdum(ifg) * rate(irnpg)
       dydt(io16) =  dydt(io16) + sum(a(1:2))
    end if


    ! ne20 reactions
    a(1) =  0.5d0 * y(ic12) * y(ic12) * rate(ir1212)
    a(2) =  y(io16) * y(ihe4) * rate(iroag)
    a(3) = -y(ine20) * y(ihe4) * rate(irneag)
    a(4) = -y(ine20) * rate(irnega)
    a(5) =  y(img24) * rate(irmgga)
    a(6) =  y(in14) * y(ihe4) * rate(irnag)

    dydt(ine20) =  dydt(ine20) + esum(a,6)


    ! mg24 reactions
    a(1) =  0.5d0 * y(ic12) * y(io16) * rate(ir1216)
    a(2) =  y(ine20) * y(ihe4) * rate(irneag)
    a(3) = -y(img24) * y(ihe4) * rate(irmgag)
    a(4) = -y(img24) * rate(irmgga)
    a(5) =  y(isi28) * rate(irsiga)

    dydt(img24) =  dydt(img24) + esum(a,5)

    if (.not.deriva) then
       a(1) = -y(img24) * y(ihe4) * rate(irmgap)*(1.0d0-rate(irr1))
       a(2) =  y(isi28) * rate(irr1) * rate(irsigp)

       dydt(img24) =  dydt(img24) + sum(a(1:2))

    else
       a(1) = -y(img24)*y(ihe4) * rate(irmgap)*(1.0d0 - ratdum(irr1))
       a(2) =  y(img24)*y(ihe4) * ratdum(irmgap)*rate(irr1)
       a(3) =  y(isi28) * ratdum(irr1) * rate(irsigp)
       a(4) =  y(isi28) * rate(irr1) * ratdum(irsigp)

       dydt(img24) =  dydt(img24) + esum(a,4)
    end if


    ! si28 reactions
    a(1) =  0.5d0 * y(ic12) * y(io16) * rate(ir1216)
    a(2) =  0.56d0 * 0.5d0*y(io16) * y(io16) * rate(ir1616)
    a(3) =  y(img24) * y(ihe4) * rate(irmgag)
    a(4) = -y(isi28) * y(ihe4) * rate(irsiag)
    a(5) = -y(isi28) * rate(irsiga)
    a(6) =  y(is32)  * rate(irsga)

    dydt(isi28) =  dydt(isi28) + esum(a,6)

    if (.not.deriva) then
       a(1) =  0.34d0*0.5d0*y(io16)*y(io16)*rate(irs1)*rate(ir1616)
       a(2) =  y(img24) * y(ihe4) * rate(irmgap)*(1.0d0-rate(irr1))
       a(3) = -y(isi28) * rate(irr1) * rate(irsigp)
       a(4) = -y(isi28) * y(ihe4) * rate(irsiap)*(1.0d0-rate(irs1))
       a(5) =  y(is32)  * rate(irs1) * rate(irsgp)

       dydt(isi28) =  dydt(isi28) + esum(a,5)

    else
       a(1)  =  0.34d0*0.5d0*y(io16)*y(io16) * ratdum(irs1)*rate(ir1616)
       a(2)  =  0.34d0*0.5d0*y(io16)*y(io16) * rate(irs1)*ratdum(ir1616)
       a(3)  =  y(img24)*y(ihe4) * rate(irmgap)*(1.0d0 - ratdum(irr1))
       a(4)  = -y(img24)*y(ihe4) * ratdum(irmgap)*rate(irr1)
       a(5)  = -y(isi28) * ratdum(irr1) * rate(irsigp)
       a(6)  = -y(isi28) * rate(irr1) * ratdum(irsigp)
       a(7)  = -y(isi28)*y(ihe4) * rate(irsiap)*(1.0d0 - ratdum(irs1))
       a(8)  =  y(isi28)*y(ihe4) * ratdum(irsiap)*rate(irs1)
       a(9)  = y(is32) * ratdum(irs1) * rate(irsgp)
       a(10) = y(is32) * rate(irs1) * ratdum(irsgp)

       dydt(isi28) =  dydt(isi28) + esum(a,10)
    end if


    ! s32 reactions
    a(1) =  0.1d0 * 0.5d0*y(io16) * y(io16) * rate(ir1616)
    a(2) =  y(isi28) * y(ihe4) * rate(irsiag)
    a(3) = -y(is32) * y(ihe4) * rate(irsag)
    a(4) = -y(is32) * rate(irsga)
    a(5) =  y(iar36) * rate(irarga)

    dydt(is32) =  dydt(is32) + esum(a,5)

    if (.not.deriva) then
       a(1) =  0.34d0*0.5d0*y(io16)*y(io16)* rate(ir1616)*(1.0d0-rate(irs1))
       a(2) =  y(isi28) * y(ihe4) * rate(irsiap)*(1.0d0-rate(irs1))
       a(3) = -y(is32) * rate(irs1) * rate(irsgp)
       a(4) = -y(is32) * y(ihe4) * rate(irsap)*(1.0d0-rate(irt1))
       a(5) =  y(iar36) * rate(irt1) * rate(irargp)

       dydt(is32) =  dydt(is32) + esum(a,5)

    else
       a(1)  =  0.34d0*0.5d0*y(io16)*y(io16) * rate(ir1616)*(1.0d0-ratdum(irs1))
       a(2)  = -0.34d0*0.5d0*y(io16)*y(io16) * ratdum(ir1616)*rate(irs1)
       a(3)  =  y(isi28)*y(ihe4) * rate(irsiap)*(1.0d0-ratdum(irs1))
       a(4)  = -y(isi28)*y(ihe4) * ratdum(irsiap)*rate(irs1)
       a(5)  = -y(is32) * ratdum(irs1) * rate(irsgp)
       a(6)  = -y(is32) * rate(irs1) * ratdum(irsgp)
       a(7)  = -y(is32)*y(ihe4) * rate(irsap)*(1.0d0-ratdum(irt1))
       a(8)  =  y(is32)*y(ihe4) * ratdum(irsap)*rate(irt1)
       a(9)  =  y(iar36) * ratdum(irt1) * rate(irargp)
       a(10) =  y(iar36) * rate(irt1) * ratdum(irargp)

       dydt(is32) =  dydt(is32) + esum(a,10)
    end if


    ! ar36 reactions
    a(1) =  y(is32)  * y(ihe4) * rate(irsag)
    a(2) = -y(iar36) * y(ihe4) * rate(irarag)
    a(3) = -y(iar36) * rate(irarga)
    a(4) =  y(ica40) * rate(ircaga)

    dydt(iar36) =  dydt(iar36) + esum(a,4)

    if (.not.deriva) then
       a(1) = y(is32)  * y(ihe4) * rate(irsap)*(1.0d0-rate(irt1))
       a(2) = -y(iar36) * rate(irt1) * rate(irargp)
       a(3) = -y(iar36) * y(ihe4) * rate(irarap)*(1.0d0-rate(iru1))
       a(4) =  y(ica40) * rate(ircagp) * rate(iru1)

       dydt(iar36) =  dydt(iar36) + esum(a,4)

    else
       a(1) =  y(is32)*y(ihe4) * rate(irsap)*(1.0d0 - ratdum(irt1))
       a(2) = -y(is32)*y(ihe4) * ratdum(irsap)*rate(irt1)
       a(3) = -y(iar36) * ratdum(irt1) * rate(irargp)
       a(4) = -y(iar36) * rate(irt1) * ratdum(irargp)
       a(5) = -y(iar36)*y(ihe4) * rate(irarap)*(1.0d0-ratdum(iru1))
       a(6) =  y(iar36)*y(ihe4) * ratdum(irarap)*rate(iru1)
       a(7) =  y(ica40) * ratdum(ircagp) * rate(iru1)
       a(8) =  y(ica40) * rate(ircagp) * ratdum(iru1)

       dydt(iar36) =  dydt(iar36) + esum(a,8)
    end if


    ! ca40 reactions
    a(1) =  y(iar36) * y(ihe4) * rate(irarag)
    a(2) = -y(ica40) * y(ihe4) * rate(ircaag)
    a(3) = -y(ica40) * rate(ircaga)
    a(4) =  y(iti44) * rate(irtiga)

    dydt(ica40) =  dydt(ica40) + esum(a,4)

    if (.not.deriva) then
       a(1) =  y(iar36) * y(ihe4) * rate(irarap)*(1.0d0-rate(iru1))
       a(2) = -y(ica40) * rate(ircagp) * rate(iru1)
       a(3) = -y(ica40) * y(ihe4) * rate(ircaap)*(1.0d0-rate(irv1))
       a(4) =  y(iti44) * rate(irtigp) * rate(irv1)

       dydt(ica40) =  dydt(ica40) + esum(a,4)

    else
       a(1) =  y(iar36)*y(ihe4) * rate(irarap)*(1.0d0-ratdum(iru1))
       a(2) = -y(iar36)*y(ihe4) * ratdum(irarap)*rate(iru1)
       a(3) = -y(ica40) * ratdum(ircagp) * rate(iru1)
       a(4) = -y(ica40) * rate(ircagp) * ratdum(iru1)
       a(5) = -y(ica40)*y(ihe4) * rate(ircaap)*(1.0d0-ratdum(irv1))
       a(6) =  y(ica40)*y(ihe4) * ratdum(ircaap)*rate(irv1)
       a(7) =  y(iti44) * ratdum(irtigp) * rate(irv1)
       a(8) =  y(iti44) * rate(irtigp) * ratdum(irv1)

       dydt(ica40) =  dydt(ica40) + esum(a,8)
    end if


    ! ti44 reactions
    a(1) =  y(ica40) * y(ihe4) * rate(ircaag)
    a(2) = -y(iti44) * y(ihe4) * rate(irtiag)
    a(3) = -y(iti44) * rate(irtiga)
    a(4) =  y(icr48) * rate(ircrga)

    dydt(iti44) =  dydt(iti44) + esum(a,4)

    if (.not.deriva) then
       a(1) =  y(ica40) * y(ihe4) * rate(ircaap)*(1.0d0-rate(irv1))
       a(2) = -y(iti44) * rate(irv1) * rate(irtigp)
       a(3) = -y(iti44) * y(ihe4) * rate(irtiap)*(1.0d0-rate(irw1))
       a(4) =  y(icr48) * rate(irw1) * rate(ircrgp)

       dydt(iti44) =  dydt(iti44) + esum(a,4)

    else
       a(1) =  y(ica40)*y(ihe4) * rate(ircaap)*(1.0d0-ratdum(irv1))
       a(2) = -y(ica40)*y(ihe4) * ratdum(ircaap)*rate(irv1)
       a(3) = -y(iti44) * ratdum(irv1) * rate(irtigp)
       a(4) = -y(iti44) * rate(irv1) * ratdum(irtigp)
       a(5) = -y(iti44)*y(ihe4) * rate(irtiap)*(1.0d0-ratdum(irw1))
       a(6) =  y(iti44)*y(ihe4) * ratdum(irtiap)*rate(irw1)
       a(7) =  y(icr48) * ratdum(irw1) * rate(ircrgp)
       a(8) =  y(icr48) * rate(irw1) * ratdum(ircrgp)

       dydt(iti44) =  dydt(iti44) + esum(a,8)
    end if


    ! cr48 reactions
    a(1) =  y(iti44) * y(ihe4) * rate(irtiag)
    a(2) = -y(icr48) * y(ihe4) * rate(ircrag)
    a(3) = -y(icr48) * rate(ircrga)
    a(4) =  y(ife52) * rate(irfega)

    dydt(icr48) =  dydt(icr48) + esum(a,4)

    if (.not.deriva) then
       a(1) =  y(iti44) * y(ihe4) * rate(irtiap)*(1.0d0-rate(irw1))
       a(2) = -y(icr48) * rate(irw1) * rate(ircrgp)
       a(3) = -y(icr48) * y(ihe4) * rate(ircrap)*(1.0d0-rate(irx1))
       a(4) =  y(ife52) * rate(irx1) * rate(irfegp)

       dydt(icr48) =  dydt(icr48) + esum(a,4)

    else
       a(1) =  y(iti44)*y(ihe4) * rate(irtiap)*(1.0d0-ratdum(irw1))
       a(2) = -y(iti44)*y(ihe4) * ratdum(irtiap)*rate(irw1)
       a(3) = -y(icr48) * ratdum(irw1) * rate(ircrgp)
       a(4) = -y(icr48) * rate(irw1) * ratdum(ircrgp)
       a(5) = -y(icr48)*y(ihe4) * rate(ircrap)*(1.0d0-ratdum(irx1))
       a(6) =  y(icr48)*y(ihe4) * ratdum(ircrap)*rate(irx1)
       a(7) =  y(ife52) * ratdum(irx1) * rate(irfegp)
       a(8) =  y(ife52) * rate(irx1) * ratdum(irfegp)

       dydt(icr48) =  dydt(icr48) + esum(a,8)
    end if



    ! fe52 reactions
    a(1) =  y(icr48) * y(ihe4) * rate(ircrag)
    a(2) = -y(ife52) * y(ihe4) * rate(irfeag)
    a(3) = -y(ife52) * rate(irfega)
    a(4) =  y(ini56) * rate(irniga)

    dydt(ife52) =  dydt(ife52) + esum(a,4)

    if (.not.deriva) then
       a(1) =  y(icr48) * y(ihe4) * rate(ircrap)*(1.0d0-rate(irx1))
       a(2) = -y(ife52) * rate(irx1) * rate(irfegp)
       dydt(ife52) =  dydt(ife52) + sum(a(1:2))

    else
       a(1) =  y(icr48)*y(ihe4) * rate(ircrap)*(1.0d0-ratdum(irx1))
       a(2) = -y(icr48)*y(ihe4) * ratdum(ircrap)*rate(irx1)
       a(3) = -y(ife52) * ratdum(irx1) * rate(irfegp)
       a(4) = -y(ife52) * rate(irx1) * ratdum(irfegp)
       dydt(ife52) =  dydt(ife52) + esum(a,4)
    end if


    ! photodisintegration reactions
    a(1) =  y(ife54) * rate(ir1f54)
    a(2) = -y(ife52) * y(ineut) * y(ineut) * rate(ir2f54)
    a(3) =  y(ife54) * y(iprot) * y(iprot) * rate(ir5f54)
    a(4) = -y(ife52) * y(ihe4) * rate(ir6f54)
    a(5) = -y(ife52) * y(ihe4) * y(iprot) * rate(ir7f54)
    a(6) =  y(ini56) * y(iprot) * rate(ir8f54)

    dydt(ife52) =  dydt(ife52) + esum(a,6)


    ! fe54 reactions
    a(1)  = -y(ife54) * rate(ir1f54)
    a(2)  =  y(ife52) * y(ineut) * y(ineut) * rate(ir2f54)
    a(3)  = -y(ife54) * y(iprot) * y(iprot) * rate(ir3f54)
    a(4)  =  y(ini56) * rate(ir4f54)
    a(5)  = -y(ife54) * y(iprot) * y(iprot) * rate(ir5f54)
    a(6)  =  y(ife52) * y(ihe4) * rate(ir6f54)
    a(7)  =  c54 * y(ini56) * rate(irn56ec)

    dydt(ife54) =  dydt(ife54) + esum(a,7)


    ! ni56 reactions
    a(1) =  y(ife52) * y(ihe4) * rate(irfeag)
    a(2) = -y(ini56) * rate(irniga)
    a(3) = -y(ini56) * rate(irn56ec)

    dydt(ini56) =  dydt(ini56) + esum(a,3)


    ! photodisintegration reactions
    a(1) =  y(ife54) * y(iprot) * y(iprot) * rate(ir3f54)
    a(2) = -y(ini56) * rate(ir4f54)
    a(3) =  y(ife52) * y(ihe4)* y(iprot) * rate(ir7f54)
    a(4) = -y(ini56) * y(iprot) * rate(ir8f54)

    dydt(ini56) =  dydt(ini56) + esum(a,4)


    ! photodisintegration neutrons
    a(1) =  2.0d0 * y(ife54) * rate(ir1f54)
    a(2) = -2.0d0 * y(ife52) * y(ineut) * y(ineut) * rate(ir2f54)
    a(3) =  2.0d0 * y(ihe4) * rate(iralf1)
    a(4) = -2.0d0 * y(ineut)*y(ineut) * y(iprot)*y(iprot) * rate(iralf2)
    a(5) =  y(iprot) * rate(irpen)
    a(6) = -y(ineut) * rate(irnep)

    dydt(ineut) =  dydt(ineut) + esum(a,6)


    ! photodisintegration protons
    a(1)  = -2.0d0 * y(ife54) * y(iprot) * y(iprot) * rate(ir3f54)
    a(2)  =  2.0d0 * y(ini56) * rate(ir4f54)
    a(3)  = -2.0d0 * y(ife54) * y(iprot) * y(iprot) * rate(ir5f54)
    a(4)  =  2.0d0 * y(ife52) * y(ihe4) * rate(ir6f54)
    a(5)  =  2.0d0 * y(ihe4) * rate(iralf1)
    a(6)  = -2.0d0 * y(ineut)*y(ineut) * y(iprot)*y(iprot) * rate(iralf2)
    a(7)  = -y(iprot) * rate(irpen)
    a(8)  =  y(ineut) * rate(irnep)

    dydt(iprot) =  dydt(iprot) + esum(a,8)

  end subroutine rhs



  subroutine aprox19rat(btemp, bden, ratraw, dratrawdt, dratrawdd)

    ! this routine generates unscreened
    ! nuclear reaction rates for the aprox19 network.

    use tfactors_module
    use aprox_rates_module
    use bl_constants_module, only: ZERO
    use extern_probin_module, only: use_c12ag_deboer17

    double precision :: btemp, bden
    double precision :: ratraw(nrates), dratrawdt(nrates), dratrawdd(nrates)

    integer          :: i
    double precision :: rrate,drratedt,drratedd
    double precision :: ff1,dff1dt,dff1dd,ff2,dff2dt,dff2dd,tot,dtotdt,invtot
    !double precision :: dtotdd
    type (tf_t)      :: tf

    do i=1,nrates
       ratraw(i)    = ZERO
       dratrawdt(i) = ZERO
       dratrawdd(i) = ZERO
    enddo

    if (btemp .lt. 1.0d6) return


    ! get the temperature factors
    call get_tfactors(btemp, tf)


    ! p(p,e+nu)d
    call rate_pp(tf,bden, &
                 ratraw(irpp),dratrawdt(irpp),dratrawdd(irpp), &
                 rrate,drratedt,drratedd)

    ! p(n,g)2d
    call rate_png(tf,bden, &
                  ratraw(irhng),dratrawdt(irhng),dratrawdd(irhng), &
                  ratraw(irdgn),dratrawdt(irdgn),dratrawdd(irdgn))

    ! d(p,g)he3
    call rate_dpg(tf,bden, &
                  ratraw(irdpg),dratrawdt(irdpg),dratrawdd(irdpg), &
                  ratraw(irhegp),dratrawdt(irhegp),dratrawdd(irhegp))

    ! he3(n,g)he4
    call rate_he3ng(tf,bden, &
                    ratraw(irheng),dratrawdt(irheng),dratrawdd(irheng), &
                    ratraw(irhegn),dratrawdt(irhegn),dratrawdd(irhegn))

    ! he3(he3,2p)he4
    call rate_he3he3(tf,bden, &
                     ratraw(ir33),dratrawdt(ir33),dratrawdd(ir33), &
                     rrate,drratedt,drratedd)

    ! he3(he4,g)be7
    call rate_he3he4(tf,bden, &
                     ratraw(irhe3ag),dratrawdt(irhe3ag),dratrawdd(irhe3ag), &
                     rrate,drratedt,drratedd)

    ! triple alpha to c12
    call rate_tripalf(tf,bden, &
                      ratraw(ir3a),dratrawdt(ir3a),dratrawdd(ir3a), &
                      ratraw(irg3a),dratrawdt(irg3a),dratrawdd(irg3a))

    ! c12(p,g)13n
    call rate_c12pg(tf,bden, &
                    ratraw(ircpg),dratrawdt(ircpg),dratrawdd(ircpg), &
                    rrate,drratedt,drratedd)

    ! n14(p,g)o15
    call rate_n14pg(tf,bden, &
                    ratraw(irnpg),dratrawdt(irnpg),dratrawdd(irnpg), &
                    rrate,drratedt,drratedd)


    ! fraction fg of n15 goes (pg) to o16, fraction fa of n15 goes (pa) to c12
    ! result is  n14+2p  goes to o16    at a rate = rnpg*fg
    !                    goes to c12+a  at a rate = rnpg*fa

    call rate_n15pg(tf,bden, &
                    ff1,dff1dt,dff1dd, &
                    rrate,drratedt,drratedd)

    call rate_n15pa(tf,bden, &
                    ff2,dff2dt,dff2dd, &
                    rrate,drratedt,drratedd)

    tot            = ff1 + ff2
    dtotdt         = dff1dt + dff2dt
!    dtotdd         = dff1dd + dff2dd
    invtot         = 1.0d0/tot

    ratraw(ifa)    = ff2 * invtot
    dratrawdt(ifa) = dff2dt * invtot - ff2 * invtot*invtot * dtotdt
!    dratrawdd(ifa) = dff2dd * invtot - ff2 * invtot*invtot * dtotdd

    ratraw(ifg)    = 1.0d0 - ratraw(ifa)
    dratrawdt(ifg) = -dratrawdt(ifa)
!    dratrawdd(ifg) = -dratrawdd(ifa)


    ! o16(p,g)f17
    call rate_o16pg(tf,bden, &
                    ratraw(iropg),dratrawdt(iropg),dratrawdd(iropg), &
                    rrate,drratedt,drratedd)

    ! n14(a,g)f18
    call rate_n14ag(tf,bden, &
                    ratraw(irnag),dratrawdt(irnag),dratrawdd(irnag), &
                    rrate,drratedt,drratedd)

    ! Determine which c12(a,g)o16 rate to use
    if (use_c12ag_deboer17) then
    ! deboer + 2017 c12(a,g)o16 rate
       call rate_c12ag_deboer17(tf,bden, &
                    ratraw(ircag),dratrawdt(ircag),dratrawdd(ircag), &
                    ratraw(iroga),dratrawdt(iroga),dratrawdd(iroga))
    else
    ! 1.7 times cf88 c12(a,g)o16 rate
       call rate_c12ag(tf,bden, &
                    ratraw(ircag),dratrawdt(ircag),dratrawdd(ircag), &
                    ratraw(iroga),dratrawdt(iroga),dratrawdd(iroga))
    endif

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

    ! fe52(n,g)fe53
    call rate_fe52ng(tf,bden, &
                     ratraw(ir52ng),dratrawdt(ir52ng),dratrawdd(ir52ng), &
                     ratraw(ir53gn),dratrawdt(ir53gn),dratrawdd(ir53gn))

    ! fe53(n,g)fe54
    call rate_fe53ng(tf,bden, &
                     ratraw(ir53ng),dratrawdt(ir53ng),dratrawdd(ir53ng), &
                     ratraw(ir54gn),dratrawdt(ir54gn),dratrawdd(ir54gn))

    ! fe54(p,g)co55
    call rate_fe54pg(tf,bden, &
                     ratraw(irfepg),dratrawdt(irfepg),dratrawdd(irfepg), &
                     ratraw(ircogp),dratrawdt(ircogp),dratrawdd(ircopg))

  end subroutine aprox19rat



  ! electron capture rates on nucleons for aprox19
  ! note they are composition dependent

  subroutine weak_aprox19(y, state, ratraw, dratrawdt, dratrawdd)

    use aprox_rates_module, only: ecapnuc, langanke

    implicit none

    double precision :: y(nspec)
    type (burn_t)    :: state
    double precision :: ratraw(nrates), dratrawdt(nrates), dratrawdd(nrates)

    double precision :: xx, spen, snep

    ! initialize
    ratraw(irpen)      = 0.0d0
    dratrawdt(irpen)   = 0.0d0
    !dratrawdd(irpen)   = 0.0d0
    ratraw(irnep)      = 0.0d0
    dratrawdt(irnep)   = 0.0d0
    !dratrawdd(irnep)   = 0.0d0
    ratraw(irn56ec)    = 0.0d0
    dratrawdt(irn56ec) = 0.0d0
    !dratrawdd(irn56ec) = 0.0d0

    if (state % T .lt. 1.0e6  .or. state % rho .lt. 1.0e-9) return

    ! get the p <-> n electron capture rates
    call ecapnuc(state % eta, state % T, ratraw(irpen), ratraw(irnep), spen, snep)

    ! ni56 electron capture rate
    call langanke(state % T, state % rho, y(ini56), state % y_e, ratraw(irn56ec), xx)

  end subroutine weak_aprox19



  subroutine screen_aprox19(btemp, bden, y, &
                            ratraw, dratrawdt, dratrawdd, &
                            ratdum, dratdumdt, dratdumdd, &
                            dratdumdy1, dratdumdy2, &
                            scfac, dscfacdt, dscfacdd)

    use bl_constants_module, only: ZERO, ONE
    use screening_module, only: screen5, plasma_state, fill_plasma_state

    ! this routine computes the screening factors
    ! and applies them to the raw reaction rates,
    ! producing the final reaction rates used by the
    ! right hand sides and jacobian matrix elements

    double precision :: btemp, bden
    double precision :: y(nspec)
    double precision :: ratraw(nrates), dratrawdt(nrates), dratrawdd(nrates)
    double precision :: ratdum(nrates), dratdumdt(nrates), dratdumdd(nrates)
    double precision :: dratdumdy1(nrates), dratdumdy2(nrates)
    double precision :: scfac(nrates),  dscfacdt(nrates),  dscfacdd(nrates)

    integer          :: i, jscr
    double precision :: sc1a,sc1adt,sc1add,sc2a,sc2adt,sc2add, &
                        sc3a,sc3adt,denom,denomdt,xx,zz
    !double precision :: sc3add, denomdd

    type (plasma_state) :: state

    ! initialize
    do i = 1, nrates
       ratdum(i)     = ratraw(i)
       dratdumdt(i)  = dratrawdt(i)
       dratdumdd(i)  = dratrawdd(i)
       dratdumdy1(i) = ZERO
       dratdumdy2(i) = ZERO
       scfac(i)      = ONE
       dscfacdt(i)   = ZERO
       dscfacdd(i)   = ZERO
    enddo



    ! Set up the state data, which is the same for all screening factors.

    call fill_plasma_state(state, btemp, bden, y(1:nspec))



    ! first the always fun triple alpha and its inverse
    jscr = 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    jscr = jscr + 1
    call screen5(state,jscr,sc2a,sc2adt,sc2add)

    sc3a   = sc1a * sc2a
    sc3adt = sc1adt*sc2a + sc1a*sc2adt
    !sc3add = sc1add*sc2a + sc1a*sc2add

    ratdum(ir3a)    = ratraw(ir3a) * sc3a
    dratdumdt(ir3a) = dratrawdt(ir3a)*sc3a + ratraw(ir3a)*sc3adt
    !dratdumdd(ir3a) = dratrawdd(ir3a)*sc3a + ratraw(ir3a)*sc3add

    scfac(ir3a) = sc3a
    dscfacdt(ir3a)  = sc3adt
    !dscfacdd(ir3a)  = sc3add

    ratdum(irg3a)    = ratraw(irg3a) * sc3a
    dratdumdt(irg3a) = dratrawdt(irg3a)*sc3a + ratraw(irg3a)*sc3adt
    !dratdumdd(irg3a) = dratrawdd(irg3a)*sc3a + ratraw(irg3a)*sc3add

    scfac(irg3a)  = sc3a
    dscfacdt(irg3a)  = sc3adt
    !dscfacdd(irg3a)  = sc3add


    ! c12 to o16
    ! c12(a,g)o16
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ratdum(ircag)     = ratraw(ircag) * sc1a
    dratdumdt(ircag)  = dratrawdt(ircag)*sc1a + ratraw(ircag)*sc1adt
    !dratdumdd(ircag)  = dratrawdd(ircag)*sc1a + ratraw(ircag)*sc1add

    scfac(ircag)  = sc1a
    dscfacdt(ircag)   = sc1adt
    !dscfacdd(ircag)   = sc1add

    ratdum(iroga)     = ratraw(iroga) * sc1a
    dratdumdt(iroga)  = dratrawdt(iroga)*sc1a + ratraw(iroga)*sc1adt
    !dratdumdd(iroga)  = dratrawdd(iroga)*sc1a + ratraw(iroga)*sc1add

    scfac(iroga)  = sc1a
    dscfacdt(iroga)   = sc1adt
    !dscfacdd(iroga)   = sc1add


    ! c12 + c12
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ratdum(ir1212)    = ratraw(ir1212) * sc1a
    dratdumdt(ir1212) = dratrawdt(ir1212)*sc1a + ratraw(ir1212)*sc1adt
    !dratdumdd(ir1212) = dratrawdd(ir1212)*sc1a + ratraw(ir1212)*sc1add

    scfac(ir1212)     = sc1a
    dscfacdt(ir1212)  = sc1adt
    !dscfacdd(ir1212)  = sc1add



    ! c12 + o16
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ratdum(ir1216)    = ratraw(ir1216) * sc1a
    dratdumdt(ir1216) = dratrawdt(ir1216)*sc1a + ratraw(ir1216)*sc1adt
    !dratdumdd(ir1216) = dratrawdd(ir1216)*sc1a + ratraw(ir1216)*sc1add

    scfac(ir1216)     = sc1a
    dscfacdt(ir1216)  = sc1adt
    !dscfacdd(ir1216)  = sc1add


    ! o16 + o16
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ratdum(ir1616)    = ratraw(ir1616) * sc1a
    dratdumdt(ir1616) = dratrawdt(ir1616)*sc1a + ratraw(ir1616)*sc1adt
    !dratdumdd(ir1616) = dratrawdd(ir1616)*sc1a + ratraw(ir1616)*sc1add

    scfac(ir1616)     = sc1a
    dscfacdt(ir1616)  = sc1adt
    !dscfacdd(ir1616)  = sc1add


    ! o16 to ne20
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! o16(a,g)ne20
    ratdum(iroag)    = ratraw(iroag) * sc1a
    dratdumdt(iroag) = dratrawdt(iroag)*sc1a + ratraw(iroag)*sc1adt
    !dratdumdd(iroag) = dratrawdd(iroag)*sc1a + ratraw(iroag)*sc1add

    scfac(iroag)  = sc1a
    dscfacdt(iroag)  = sc1adt
    !dscfacdd(iroag)  = sc1add

    ratdum(irnega)    = ratraw(irnega) * sc1a
    dratdumdt(irnega) = dratrawdt(irnega)*sc1a + ratraw(irnega)*sc1adt
    !dratdumdd(irnega) = dratrawdd(irnega)*sc1a + ratraw(irnega)*sc1add

    scfac(irnega)  = sc1a
    dscfacdt(irnega)  = sc1adt
    !dscfacdd(irnega)  = sc1add


    ! ne20 to mg24
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! ne20(a,g)mg24
    ratdum(irneag)    = ratraw(irneag) * sc1a
    dratdumdt(irneag) = dratrawdt(irneag)*sc1a + ratraw(irneag)*sc1adt
    !dratdumdd(irneag) = dratrawdd(irneag)*sc1a + ratraw(irneag)*sc1add

    scfac(irneag) = sc1a
    dscfacdt(irneag)  = sc1adt
    !dscfacdd(irneag)  = sc1add

    ratdum(irmgga)    = ratraw(irmgga) * sc1a
    dratdumdt(irmgga) = dratrawdt(irmgga)*sc1a + ratraw(irmgga)*sc1adt
    !dratdumdd(irmgga) = dratrawdd(irmgga)*sc1a + ratraw(irmgga)*sc1add

    scfac(irmgga) = sc1a
    dscfacdt(irmgga)  = sc1adt
    !dscfacdd(irmgga)  = sc1add


    ! mg24 to si28
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! mg24(a,g)si28
    ratdum(irmgag)    = ratraw(irmgag) * sc1a
    dratdumdt(irmgag) = dratrawdt(irmgag)*sc1a + ratraw(irmgag)*sc1adt
    !dratdumdd(irmgag) = dratrawdd(irmgag)*sc1a + ratraw(irmgag)*sc1add

    scfac(irmgag) = sc1a
    dscfacdt(irmgag)  = sc1adt
    !dscfacdd(irmgag)  = sc1add

    ratdum(irsiga)    = ratraw(irsiga) * sc1a
    dratdumdt(irsiga) = dratrawdt(irsiga)*sc1a + ratraw(irsiga)*sc1adt
    !dratdumdd(irsiga) = dratrawdd(irsiga)*sc1a + ratraw(irsiga)*sc1add

    scfac(irsiga) = sc1a
    dscfacdt(irsiga)  = sc1adt
    !dscfacdd(irsiga)  = sc1add


    ! mg24(a,p)al27
    ratdum(irmgap)    = ratraw(irmgap) * sc1a
    dratdumdt(irmgap) = dratrawdt(irmgap)*sc1a + ratraw(irmgap)*sc1adt
    !dratdumdd(irmgap) = dratrawdd(irmgap)*sc1a + ratraw(irmgap)*sc1add

    scfac(irmgap)     = sc1a
    dscfacdt(irmgap)  = sc1adt
    !dscfacdd(irmgap)  = sc1add

    ratdum(iralpa)    = ratraw(iralpa) * sc1a
    dratdumdt(iralpa) = dratrawdt(iralpa)*sc1a + ratraw(iralpa)*sc1adt
    !dratdumdd(iralpa) = dratrawdd(iralpa)*sc1a + ratraw(iralpa)*sc1add

    scfac(iralpa)     = sc1a
    dscfacdt(iralpa)  = sc1adt
    !dscfacdd(iralpa)  = sc1add


    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ! al27(p,g)si28
    ratdum(iralpg)    = ratraw(iralpg) * sc1a
    dratdumdt(iralpg) = dratrawdt(iralpg)*sc1a + ratraw(iralpg)*sc1adt
    !dratdumdd(iralpg) = dratrawdd(iralpg)*sc1a + ratraw(iralpg)*sc1add

    scfac(iralpg)     = sc1a
    dscfacdt(iralpg)  = sc1adt
    !dscfacdd(iralpg)  = sc1add

    ratdum(irsigp)    = ratraw(irsigp) * sc1a
    dratdumdt(irsigp) = dratrawdt(irsigp)*sc1a + ratraw(irsigp)*sc1adt
    !dratdumdd(irsigp) = dratrawdd(irsigp)*sc1a + ratraw(irsigp)*sc1add

    scfac(irsigp)     = sc1a
    dscfacdt(irsigp)  = sc1adt
    !dscfacdd(irsigp)  = sc1add



    ! si28 to s32
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! si28(a,g)s32
    ratdum(irsiag)    = ratraw(irsiag) * sc1a
    dratdumdt(irsiag) = dratrawdt(irsiag)*sc1a + ratraw(irsiag)*sc1adt
    !dratdumdd(irsiag) = dratrawdd(irsiag)*sc1a + ratraw(irsiag)*sc1add

    scfac(irsiag)     = sc1a
    dscfacdt(irsiag)  = sc1adt
    !dscfacdd(irsiag)  = sc1add

    ratdum(irsga)    = ratraw(irsga) * sc1a
    dratdumdt(irsga) = dratrawdt(irsga)*sc1a + ratraw(irsga)*sc1adt
    !dratdumdd(irsga) = dratrawdd(irsga)*sc1a + ratraw(irsga)*sc1add

    scfac(irsga)     = sc1a
    dscfacdt(irsga)  = sc1adt
    !dscfacdd(irsga)  = sc1add


    ! si28(a,p)p31
    ratdum(irsiap)    = ratraw(irsiap) * sc1a
    dratdumdt(irsiap) = dratrawdt(irsiap)*sc1a + ratraw(irsiap)*sc1adt
    !dratdumdd(irsiap) = dratrawdd(irsiap)*sc1a + ratraw(irsiap)*sc1add

    scfac(irsiap)     = sc1a
    dscfacdt(irsiap)  = sc1adt
    !dscfacdd(irsiap)  = sc1add

    ratdum(irppa)     = ratraw(irppa) * sc1a
    dratdumdt(irppa)  = dratrawdt(irppa)*sc1a  + ratraw(irppa)*sc1adt
    !dratdumdd(irppa)  = dratrawdd(irppa)*sc1a  + ratraw(irppa)*sc1add

    scfac(irppa)      = sc1a
    dscfacdt(irppa)   = sc1adt
    !dscfacdd(irppa)   = sc1add


    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ! p31(p,g)s32
    ratdum(irppg)     = ratraw(irppg) * sc1a
    dratdumdt(irppg)  = dratrawdt(irppg)*sc1a + ratraw(irppg)*sc1adt
    !dratdumdd(irppg)  = dratrawdd(irppg)*sc1a + ratraw(irppg)*sc1add

    scfac(irppg)      = sc1a
    dscfacdt(irppg)   = sc1adt
    !dscfacdd(irppg)   = sc1add

    ratdum(irsgp)     = ratraw(irsgp) * sc1a
    dratdumdt(irsgp)  = dratrawdt(irsgp)*sc1a + ratraw(irsgp)*sc1adt
    !dratdumdd(irsgp)  = dratrawdd(irsgp)*sc1a + ratraw(irsgp)*sc1add

    scfac(irsgp)      = sc1a
    dscfacdt(irsgp)   = sc1adt
    !dscfacdd(irsgp)   = sc1add



    ! s32 to ar36
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! s32(a,g)ar36
    ratdum(irsag)     = ratraw(irsag) * sc1a
    dratdumdt(irsag)  = dratrawdt(irsag)*sc1a + ratraw(irsag)*sc1adt
    !dratdumdd(irsag)  = dratrawdd(irsag)*sc1a + ratraw(irsag)*sc1add

    scfac(irsag)      = sc1a
    dscfacdt(irsag)   = sc1adt
    !dscfacdd(irsag)   = sc1add

    ratdum(irarga)     = ratraw(irarga) * sc1a
    dratdumdt(irarga)  = dratrawdt(irarga)*sc1a + ratraw(irarga)*sc1adt
    !dratdumdd(irarga)  = dratrawdd(irarga)*sc1a + ratraw(irarga)*sc1add

    scfac(irarga)      = sc1a
    dscfacdt(irarga)   = sc1adt
    !dscfacdd(irarga)   = sc1add

    ! s32(a,p)cl35
    ratdum(irsap)     = ratraw(irsap) * sc1a
    dratdumdt(irsap)  = dratrawdt(irsap)*sc1a + ratraw(irsap)*sc1adt
    !dratdumdd(irsap)  = dratrawdd(irsap)*sc1a + ratraw(irsap)*sc1add

    scfac(irsap)      = sc1a
    dscfacdt(irsap)   = sc1adt
    !dscfacdd(irsap)   = sc1add

    ratdum(irclpa)    = ratraw(irclpa) * sc1a
    dratdumdt(irclpa) = dratrawdt(irclpa)*sc1a + ratraw(irclpa)*sc1adt
    !dratdumdd(irclpa) = dratrawdd(irclpa)*sc1a + ratraw(irclpa)*sc1add

    scfac(irclpa)     = sc1a
    dscfacdt(irclpa)  = sc1adt
    !dscfacdd(irclpa)  = sc1add


    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ! cl35(p,g)ar36
    ratdum(irclpg)    = ratraw(irclpg) * sc1a
    dratdumdt(irclpg) = dratrawdt(irclpg)*sc1a + ratraw(irclpg)*sc1adt
    !dratdumdd(irclpg) = dratrawdd(irclpg)*sc1a + ratraw(irclpg)*sc1add

    scfac(irclpg)     = sc1a
    dscfacdt(irclpg)  = sc1adt
    !dscfacdd(irclpg)  = sc1add

    ratdum(irargp)    = ratraw(irargp) * sc1a
    dratdumdt(irargp) = dratrawdt(irargp)*sc1a + ratraw(irargp)*sc1adt
    !dratdumdd(irargp) = dratrawdd(irargp)*sc1a + ratraw(irargp)*sc1add

    scfac(irargp)     = sc1a
    dscfacdt(irargp)  = sc1adt
    !dscfacdd(irargp)  = sc1add



    ! ar36 to ca40
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! ar36(a,g)ca40
    ratdum(irarag)    = ratraw(irarag) * sc1a
    dratdumdt(irarag) = dratrawdt(irarag)*sc1a + ratraw(irarag)*sc1adt
    !dratdumdd(irarag) = dratrawdd(irarag)*sc1a + ratraw(irarag)*sc1add

    scfac(irarag)     = sc1a
    dscfacdt(irarag)  = sc1adt
    !dscfacdd(irarag)  = sc1add

    ratdum(ircaga)    = ratraw(ircaga) * sc1a
    dratdumdt(ircaga) = dratrawdt(ircaga)*sc1a + ratraw(ircaga)*sc1adt
    !dratdumdd(ircaga) = dratrawdd(ircaga)*sc1a + ratraw(ircaga)*sc1add

    scfac(ircaga)     = sc1a
    dscfacdt(ircaga)  = sc1adt
    !dscfacdd(ircaga)  = sc1add


    ! ar36(a,p)k39
    ratdum(irarap)    = ratraw(irarap) * sc1a
    dratdumdt(irarap) = dratrawdt(irarap)*sc1a + ratraw(irarap)*sc1adt
    !dratdumdd(irarap) = dratrawdd(irarap)*sc1a + ratraw(irarap)*sc1add

    scfac(irarap)     = sc1a
    dscfacdt(irarap)  = sc1adt
    !dscfacdd(irarap)  = sc1add

    ratdum(irkpa)     = ratraw(irkpa) * sc1a
    dratdumdt(irkpa)  = dratrawdt(irkpa)*sc1a  + ratraw(irkpa)*sc1adt
    !dratdumdd(irkpa)  = dratrawdd(irkpa)*sc1a  + ratraw(irkpa)*sc1add

    scfac(irkpa)      = sc1a
    dscfacdt(irkpa)   = sc1adt
    !dscfacdd(irkpa)   = sc1add


    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ! k39(p,g)ca40
    ratdum(irkpg)     = ratraw(irkpg) * sc1a
    dratdumdt(irkpg)  = dratrawdt(irkpg)*sc1a  + ratraw(irkpg)*sc1adt
    !dratdumdd(irkpg)  = dratrawdd(irkpg)*sc1a  + ratraw(irkpg)*sc1add

    scfac(irkpg)      = sc1a
    dscfacdt(irkpg)   = sc1adt
    !dscfacdd(irkpg)   = sc1add

    ratdum(ircagp)     = ratraw(ircagp) * sc1a
    dratdumdt(ircagp)  = dratrawdt(ircagp)*sc1a  + ratraw(ircagp)*sc1adt
    !dratdumdd(ircagp)  = dratrawdd(ircagp)*sc1a  + ratraw(ircagp)*sc1add

    scfac(ircagp)      = sc1a
    dscfacdt(ircagp)   = sc1adt
    !dscfacdd(ircagp)   = sc1add



    ! ca40 to ti44
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! ca40(a,g)ti44
    ratdum(ircaag)    = ratraw(ircaag) * sc1a
    dratdumdt(ircaag) = dratrawdt(ircaag)*sc1a + ratraw(ircaag)*sc1adt
    !dratdumdd(ircaag) = dratrawdd(ircaag)*sc1a + ratraw(ircaag)*sc1add

    scfac(ircaag)     = sc1a
    dscfacdt(ircaag)  = sc1adt
    !dscfacdd(ircaag)  = sc1add

    ratdum(irtiga)    = ratraw(irtiga) * sc1a
    dratdumdt(irtiga) = dratrawdt(irtiga)*sc1a + ratraw(irtiga)*sc1adt
    !dratdumdd(irtiga) = dratrawdd(irtiga)*sc1a + ratraw(irtiga)*sc1add

    scfac(irtiga)     = sc1a
    dscfacdt(irtiga)  = sc1adt
    !dscfacdd(irtiga)  = sc1add


    ! ca40(a,p)sc43
    ratdum(ircaap)    = ratraw(ircaap) * sc1a
    dratdumdt(ircaap) = dratrawdt(ircaap)*sc1a + ratraw(ircaap)*sc1adt
    !dratdumdd(ircaap) = dratrawdd(ircaap)*sc1a + ratraw(ircaap)*sc1add

    scfac(ircaap)     = sc1a
    dscfacdt(ircaap)  = sc1adt
    !dscfacdd(ircaap)  = sc1add

    ratdum(irscpa)    = ratraw(irscpa) * sc1a
    dratdumdt(irscpa) = dratrawdt(irscpa)*sc1a + ratraw(irscpa)*sc1adt
    !dratdumdd(irscpa) = dratrawdd(irscpa)*sc1a + ratraw(irscpa)*sc1add

    scfac(irscpa)     = sc1a
    dscfacdt(irscpa)  = sc1adt
    !dscfacdd(irscpa)  = sc1add


    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ! sc43(p,g)ti44
    ratdum(irscpg)    = ratraw(irscpg) * sc1a
    dratdumdt(irscpg) = dratrawdt(irscpg)*sc1a + ratraw(irscpg)*sc1adt
    !dratdumdd(irscpg) = dratrawdd(irscpg)*sc1a + ratraw(irscpg)*sc1add

    scfac(irscpg)     = sc1a
    dscfacdt(irscpg)  = sc1adt
    !dscfacdd(irscpg)  = sc1add

    ratdum(irtigp)    = ratraw(irtigp) * sc1a
    dratdumdt(irtigp) = dratrawdt(irtigp)*sc1a + ratraw(irtigp)*sc1adt
    !dratdumdd(irtigp) = dratrawdd(irtigp)*sc1a + ratraw(irtigp)*sc1add

    scfac(irtigp)     = sc1a
    dscfacdt(irtigp)  = sc1adt
    !dscfacdd(irtigp)  = sc1add



    ! ti44 to cr48
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! ti44(a,g)cr48
    ratdum(irtiag)    = ratraw(irtiag) * sc1a
    dratdumdt(irtiag) = dratrawdt(irtiag)*sc1a + ratraw(irtiag)*sc1adt
    !dratdumdd(irtiag) = dratrawdd(irtiag)*sc1a + ratraw(irtiag)*sc1add

    scfac(irtiag)     = sc1a
    dscfacdt(irtiag)  = sc1adt
    !dscfacdd(irtiag)  = sc1add

    ratdum(ircrga)    = ratraw(ircrga) * sc1a
    dratdumdt(ircrga) = dratrawdt(ircrga)*sc1a + ratraw(ircrga)*sc1adt
    !dratdumdd(ircrga) = dratrawdd(ircrga)*sc1a + ratraw(ircrga)*sc1add

    scfac(ircrga)     = sc1a
    dscfacdt(ircrga)  = sc1adt
    !dscfacdd(ircrga)  = sc1add

    ! ti44(a,p)v47
    ratdum(irtiap)    = ratraw(irtiap) * sc1a
    dratdumdt(irtiap) = dratrawdt(irtiap)*sc1a + ratraw(irtiap)*sc1adt
    !dratdumdd(irtiap) = dratrawdd(irtiap)*sc1a + ratraw(irtiap)*sc1add

    scfac(irtiap)  = sc1a
    dscfacdt(irtiap)  = sc1adt
    !dscfacdd(irtiap)  = sc1add

    ratdum(irvpa)     = ratraw(irvpa) * sc1a
    dratdumdt(irvpa)  = dratrawdt(irvpa)*sc1a  + ratraw(irvpa)*sc1adt
    !dratdumdd(irvpa)  = dratrawdd(irvpa)*sc1a  + ratraw(irvpa)*sc1add

    scfac(irvpa)      = sc1a
    dscfacdt(irvpa)   = sc1adt
    !dscfacdd(irvpa)   = sc1add


    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ! v47(p,g)cr48
    ratdum(irvpg)     = ratraw(irvpg) * sc1a
    dratdumdt(irvpg)  = dratrawdt(irvpg)*sc1a  + ratraw(irvpg)*sc1adt
    !dratdumdd(irvpg)  = dratrawdd(irvpg)*sc1a  + ratraw(irvpg)*sc1add

    scfac(irvpg)      = sc1a
    dscfacdt(irvpg)   = sc1adt
    !dscfacdd(irvpg)   = sc1add

    ratdum(ircrgp)     = ratraw(ircrgp) * sc1a
    dratdumdt(ircrgp)  = dratrawdt(ircrgp)*sc1a  + ratraw(ircrgp)*sc1adt
    !dratdumdd(ircrgp)  = dratrawdd(ircrgp)*sc1a  + ratraw(ircrgp)*sc1add

    scfac(ircrgp)      = sc1a
    dscfacdt(ircrgp)   = sc1adt
    !dscfacdd(ircrgp)   = sc1add



    ! cr48 to fe52
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! cr48(a,g)fe52
    ratdum(ircrag)    = ratraw(ircrag) * sc1a
    dratdumdt(ircrag) = dratrawdt(ircrag)*sc1a + ratraw(ircrag)*sc1adt
    !dratdumdd(ircrag) = dratrawdd(ircrag)*sc1a + ratraw(ircrag)*sc1add

    scfac(ircrag)     = sc1a
    dscfacdt(ircrag)  = sc1adt
    !dscfacdd(ircrag)  = sc1add

    ratdum(irfega)    = ratraw(irfega) * sc1a
    dratdumdt(irfega) = dratrawdt(irfega)*sc1a + ratraw(irfega)*sc1adt
    !dratdumdd(irfega) = dratrawdd(irfega)*sc1a + ratraw(irfega)*sc1add

    scfac(irfega)     = sc1a
    dscfacdt(irfega)  = sc1adt
    !dscfacdd(irfega)  = sc1add


    ! cr48(a,p)mn51
    ratdum(ircrap)    = ratraw(ircrap) * sc1a
    dratdumdt(ircrap) = dratrawdt(ircrap)*sc1a + ratraw(ircrap)*sc1adt
    !dratdumdd(ircrap) = dratrawdd(ircrap)*sc1a + ratraw(ircrap)*sc1add

    scfac(ircrap)     = sc1a
    dscfacdt(ircrap)  = sc1adt
    !dscfacdd(ircrap)  = sc1add

    ratdum(irmnpa)    = ratraw(irmnpa) * sc1a
    dratdumdt(irmnpa) = dratrawdt(irmnpa)*sc1a + ratraw(irmnpa)*sc1adt
    !dratdumdd(irmnpa) = dratrawdd(irmnpa)*sc1a + ratraw(irmnpa)*sc1add

    scfac(irmnpa)     = sc1a
    dscfacdt(irmnpa)  = sc1adt
    !dscfacdd(irmnpa)  = sc1add


    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ! mn51(p,g)fe52
    ratdum(irmnpg)    = ratraw(irmnpg) * sc1a
    dratdumdt(irmnpg) = dratrawdt(irmnpg)*sc1a + ratraw(irmnpg)*sc1adt
    !dratdumdd(irmnpg) = dratrawdd(irmnpg)*sc1a + ratraw(irmnpg)*sc1add

    scfac(irmnpg)     = sc1a
    dscfacdt(irmnpg)  = sc1adt
    !dscfacdd(irmnpg)  = sc1add

    ratdum(irfegp)    = ratraw(irfegp) * sc1a
    dratdumdt(irfegp) = dratrawdt(irfegp)*sc1a + ratraw(irfegp)*sc1adt
    !dratdumdd(irfegp) = dratrawdd(irfegp)*sc1a + ratraw(irfegp)*sc1add

    scfac(irfegp)     = sc1a
    dscfacdt(irfegp)  = sc1adt
    !dscfacdd(irfegp)  = sc1add



    ! fe52 to ni56
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! fe52(a,g)ni56
    ratdum(irfeag)    = ratraw(irfeag) * sc1a
    dratdumdt(irfeag) = dratrawdt(irfeag)*sc1a + ratraw(irfeag)*sc1adt
    !dratdumdd(irfeag) = dratrawdd(irfeag)*sc1a + ratraw(irfeag)*sc1add

    scfac(irfeag)     = sc1a
    dscfacdt(irfeag)  = sc1adt
    !dscfacdd(irfeag)  = sc1add

    ratdum(irniga)    = ratraw(irniga) * sc1a
    dratdumdt(irniga) = dratrawdt(irniga)*sc1a + ratraw(irniga)*sc1adt
    !dratdumdd(irniga) = dratrawdd(irniga)*sc1a + ratraw(irniga)*sc1add

    scfac(irniga)     = sc1a
    dscfacdt(irniga)  = sc1adt
    !dscfacdd(irniga)  = sc1add


    ! fe52(a,p)co55
    ratdum(irfeap) = ratraw(irfeap) * sc1a
    dratdumdt(irfeap) = dratrawdt(irfeap)*sc1a + ratraw(irfeap)*sc1adt
    !dratdumdd(irfeap) = dratrawdd(irfeap)*sc1a + ratraw(irfeap)*sc1add

    scfac(irfeap)     = sc1a
    dscfacdt(irfeap)  = sc1adt
    !dscfacdd(irfeap)  = sc1add

    ratdum(ircopa)    = ratraw(ircopa) * sc1a
    dratdumdt(ircopa) = dratrawdt(ircopa)*sc1a + ratraw(ircopa)*sc1adt
    !dratdumdd(ircopa) = dratrawdd(ircopa)*sc1a + ratraw(ircopa)*sc1add

    scfac(ircopa)     = sc1a
    dscfacdt(ircopa)  = sc1adt
    !dscfacdd(ircopa)  = sc1add


    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! co55(p,g)ni56
    ratdum(ircopg)    = ratraw(ircopg) * sc1a
    dratdumdt(ircopg) = dratrawdt(ircopg)*sc1a + ratraw(ircopg)*sc1adt
    !dratdumdd(ircopg) = dratrawdd(ircopg)*sc1a + ratraw(ircopg)*sc1add

    scfac(ircopg)     = sc1a
    dscfacdt(ircopg)  = sc1adt
    !dscfacdd(ircopg)  = sc1add

    ratdum(irnigp)    = ratraw(irnigp) * sc1a
    dratdumdt(irnigp) = dratrawdt(irnigp)*sc1a + ratraw(irnigp)*sc1adt
    !dratdumdd(irnigp) = dratrawdd(irnigp)*sc1a + ratraw(irnigp)*sc1add

    scfac(irnigp)     = sc1a
    dscfacdt(irnigp)  = sc1adt
    !dscfacdd(irniga)  = sc1add


    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ! fe54(p,g)co55
    ratdum(irfepg)    = ratraw(irfepg) * sc1a
    dratdumdt(irfepg) = dratrawdt(irfepg)*sc1a + ratraw(irfepg)*sc1adt
    !dratdumdd(irfepg) = dratrawdd(irfepg)*sc1a + ratraw(irfepg)*sc1add

    scfac(irfepg)     = sc1a
    dscfacdt(irfepg)  = sc1adt
    !dscfacdd(irfepg)  = sc1add

    ratdum(ircogp)    = ratraw(ircogp) * sc1a
    dratdumdt(ircogp) = dratrawdt(ircogp)*sc1a + ratraw(ircogp)*sc1adt
    !dratdumdd(ircogp) = dratrawdd(ircogp)*sc1a + ratraw(ircogp)*sc1add

    scfac(ircogp)     = sc1a
    dscfacdt(ircogp)  = sc1adt
    !dscfacdd(ircogp)  = sc1add



    ! d(p,g)he4
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! d(p,g)he3
    ratdum(irdpg)    = ratraw(irdpg) * sc1a
    dratdumdt(irdpg) = dratrawdt(irdpg)*sc1a + ratraw(irdpg)*sc1adt
    !dratdumdd(irdpg) = dratrawdd(irdpg)*sc1a + ratraw(irdpg)*sc1add

    scfac(irdpg)     = sc1a
    dscfacdt(irdpg)  = sc1adt
    !dscfacdd(irdpg) = sc1add

    ratdum(irhegp)    = ratraw(irhegp) * sc1a
    dratdumdt(irhegp) = dratrawdt(irhegp)*sc1a + ratraw(irhegp)*sc1adt
    !dratdumdd(irhegp) = dratrawdd(irhegp)*sc1a + ratraw(irhegp)*sc1add

    scfac(irhegp)     = sc1a
    dscfacdt(irhegp)  = sc1adt
    !dscfacdd(irhegp) = sc1add


    ! pp
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ratdum(irpp)    = ratraw(irpp) * sc1a
    dratdumdt(irpp) = dratrawdt(irpp)*sc1a + ratraw(irpp)*sc1adt
    !dratdumdd(irpp) = dratrawdd(irpp)*sc1a + ratraw(irpp)*sc1add

    scfac(irpp)     = sc1a
    dscfacdt(irpp)  = sc1adt
    !dscfacdd(irpp)  = sc1add


    ! he3 + he3
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ! he3(he3,2p)he4
    ratdum(ir33)    = ratraw(ir33) * sc1a
    dratdumdt(ir33) = dratrawdt(ir33)*sc1a + ratraw(ir33)*sc1adt
    !dratdumdd(ir33) = dratrawdd(ir33)*sc1a + ratraw(ir33)*sc1add

    scfac(ir33)     = sc1a
    dscfacdt(ir33)  = sc1adt
!    dscfacdd(ir33)  = sc1add


    ! he3 + he4
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ! he3(he4,g)be7
    ratdum(irhe3ag)    = ratraw(irhe3ag) * sc1a
    dratdumdt(irhe3ag) = dratrawdt(irhe3ag)*sc1a &
         + ratraw(irhe3ag)*sc1adt
!    dratdumdd(irhe3ag) = dratrawdd(irhe3ag)*sc1a &
!         + ratraw(irhe3ag)*sc1add

    scfac(irhe3ag)     = sc1a
    dscfacdt(irhe3ag)  = sc1adt
!    dscfacdd(irhe3ag)  = sc1add


    ! cno cycles

    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ! c12(p,g)13n
    ratdum(ircpg)    = ratraw(ircpg) * sc1a
    dratdumdt(ircpg) = dratrawdt(ircpg)*sc1a + ratraw(ircpg)*sc1adt
    !dratdumdd(ircpg) = dratrawdd(ircpg)*sc1a + ratraw(ircpg)*sc1add

    scfac(ircpg)     = sc1a
    dscfacdt(ircpg)  = sc1adt
    !dscfacdd(ircpg)  = sc1add


    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ! n14(p,g)o15
    ratdum(irnpg)    = ratraw(irnpg) * sc1a
    dratdumdt(irnpg) = dratrawdt(irnpg)*sc1a + ratraw(irnpg)*sc1adt
    !dratdumdd(irnpg) = dratrawdd(irnpg)*sc1a + ratraw(irnpg)*sc1add

    scfac(irnpg)     = sc1a
    dscfacdt(irnpg)  = sc1adt
    !dscfacdd(irnpg)  = sc1add


    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! o16(p,g)f17
    ratdum(iropg)    = ratraw(iropg) * sc1a
    dratdumdt(iropg) = dratrawdt(iropg)*sc1a+ratraw(iropg)*sc1adt
    !dratdumdd(iropg) = dratrawdd(iropg)*sc1a+ratraw(iropg)*sc1add

    scfac(iropg)     = sc1a
    dscfacdt(iropg)  = sc1adt
    !dscfacdd(iropg)  = sc1add

    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ! n14(a,g)f18
    ratdum(irnag)    = ratraw(irnag) * sc1a
    dratdumdt(irnag) = dratrawdt(irnag)*sc1a + ratraw(irnag)*sc1adt
    !dratdumdd(irnag) = dratrawdd(irnag)*sc1a + ratraw(irnag)*sc1add

    scfac(irnag)     = sc1a
    dscfacdt(irnag)  = sc1adt
    !dscfacdd(irnag)  = sc1add



    ! now form those lovely equilibrium rates

    ! mg24(a,p)27al(p,g)28si
    ratdum(irr1)     = ZERO
    dratdumdt(irr1)  = ZERO
    !dratdumdd(irr1)  = ZERO
    denom    = ratdum(iralpa) + ratdum(iralpg)
    denomdt  = dratdumdt(iralpa) + dratdumdt(iralpg)
    !denomdd  = dratdumdd(iralpa) + dratdumdd(iralpg)
    if (denom .gt. 1.0d-50) then
       zz = 1.0d0/denom
       ratdum(irr1)    = ratdum(iralpa)*zz
       dratdumdt(irr1) = (dratdumdt(iralpa) - ratdum(irr1)*denomdt)*zz
       !dratdumdd(irr1) = (dratdumdd(iralpa) - ratdum(irr1)*denomdd)*zz
    end if

    ! si28(a,p)p31(p,g)s32
    ratdum(irs1)     = ZERO
    dratdumdt(irs1)  = ZERO
    !dratdumdd(irs1)  = ZERO
    denom    = ratdum(irppa) + ratdum(irppg)
    denomdt  = dratdumdt(irppa) + dratdumdt(irppg)
    !denomdd  = dratdumdd(irppa) + dratdumdd(irppg)
    if (denom .gt. 1.0d-50) then
       zz = 1.0d0/denom
       ratdum(irs1)    = ratdum(irppa)*zz
       dratdumdt(irs1) = (dratdumdt(irppa) - ratdum(irs1)*denomdt)*zz
       !dratdumdd(irs1) = (dratdumdd(irppa) - ratdum(irs1)*denomdd)*zz
    end if

    ! s32(a,p)cl35(p,g)ar36
    ratdum(irt1)     = ZERO
    dratdumdt(irt1)  = ZERO
    !dratdumdd(irt1)  = ZERO
    denom    = ratdum(irclpa) + ratdum(irclpg)
    denomdt  = dratdumdt(irclpa) + dratdumdt(irclpg)
    !denomdd  = dratdumdd(irclpa) + dratdumdd(irclpg)
    if (denom .gt. 1.0d-50) then
       zz = 1.0d0/denom
       ratdum(irt1)    = ratdum(irclpa)*zz
       dratdumdt(irt1) = (dratdumdt(irclpa) - ratdum(irt1)*denomdt)*zz
       !dratdumdd(irt1) = (dratdumdd(irclpa) - ratdum(irt1)*denomdd)*zz
    end if

    ! ar36(a,p)k39(p,g)ca40
    ratdum(iru1)     = ZERO
    dratdumdt(iru1)  = ZERO
    !dratdumdd(iru1)  = ZERO
    denom    = ratdum(irkpa) + ratdum(irkpg)
    denomdt  = dratdumdt(irkpa) + dratdumdt(irkpg)
    !denomdd  = dratdumdd(irkpa) + dratdumdd(irkpg)
    if (denom .gt. 1.0d-50) then
       zz   = 1.0d0/denom
       ratdum(iru1)   = ratdum(irkpa)*zz
       dratdumdt(iru1) = (dratdumdt(irkpa) - ratdum(iru1)*denomdt)*zz
       !dratdumdd(iru1) = (dratdumdd(irkpa) - ratdum(iru1)*denomdd)*zz
    end if

    ! ca40(a,p)sc43(p,g)ti44
    ratdum(irv1)     = ZERO
    dratdumdt(irv1)  = ZERO
    !dratdumdd(irv1)  = ZERO
    denom    = ratdum(irscpa) + ratdum(irscpg)
    denomdt  = dratdumdt(irscpa) + dratdumdt(irscpg)
    !denomdd  = dratdumdd(irscpa) + dratdumdd(irscpg)
    if (denom .gt. 1.0d-50) then
       zz  = 1.0d0/denom
       ratdum(irv1)    = ratdum(irscpa)*zz
       dratdumdt(irv1) = (dratdumdt(irscpa) - ratdum(irv1)*denomdt)*zz
       !dratdumdd(irv1) = (dratdumdd(irscpa) - ratdum(irv1)*denomdd)*zz
    end if

    ! ti44(a,p)v47(p,g)cr48
    ratdum(irw1)    = ZERO
    dratdumdt(irw1) = ZERO
    !dratdumdd(irw1) = ZERO
    denom    = ratdum(irvpa) + ratdum(irvpg)
    denomdt  = dratdumdt(irvpa) + dratdumdt(irvpg)
    !denomdd  = dratdumdd(irvpa) + dratdumdd(irvpg)
    if (denom .gt. 1.0d-50) then
       zz = 1.0d0/denom
       ratdum(irw1)    = ratdum(irvpa)*zz
       dratdumdt(irw1) = (dratdumdt(irvpa) - ratdum(irw1)*denomdt)*zz
       !dratdumdd(irw1) = (dratdumdd(irvpa) - ratdum(irw1)*denomdd)*zz
    end if

    ! cr48(a,p)mn51(p,g)fe52
    ratdum(irx1)    = ZERO
    dratdumdt(irx1) = ZERO
    !dratdumdd(irx1) = ZERO
    denom    = ratdum(irmnpa) + ratdum(irmnpg)
    denomdt  = dratdumdt(irmnpa) + dratdumdt(irmnpg)
    !denomdd  = dratdumdd(irmnpa) + dratdumdd(irmnpg)
    if (denom .gt. 1.0d-50) then
       zz = 1.0d0/denom
       ratdum(irx1)    = ratdum(irmnpa)*zz
       dratdumdt(irx1) = (dratdumdt(irmnpa) - ratdum(irx1)*denomdt)*zz
       !dratdumdd(irx1) = (dratdumdd(irmnpa) - ratdum(irx1)*denomdd)*zz
    endif


    ! fe52(n,g)fe53(n,g)fe54 equilibrium links

    ratdum(ir1f54) = ZERO
    dratdumdy1(ir1f54) = ZERO
    dratdumdt(ir1f54)  = ZERO
    !dratdumdd(ir1f54)  = ZERO
    ratdum(ir2f54) = ZERO
    dratdumdy1(ir2f54) = ZERO
    dratdumdt(ir2f54)  = ZERO
    !dratdumdd(ir2f54)  = ZERO

    denom   = ratdum(ir53gn) + y(ineut)*ratdum(ir53ng)
    denomdt = dratdumdt(ir53gn) + y(ineut)*dratdumdt(ir53ng)
    !denomdd = dratdumdd(ir53gn) + y(ineut)*dratdumdd(ir53ng)

    if (denom .gt. 1.0d-50 .and. btemp .gt. 1.5e9) then
       zz      = 1.0d0/denom

       ratdum(ir1f54)     = ratdum(ir54gn)*ratdum(ir53gn)*zz
       dratdumdy1(ir1f54) = -ratdum(ir1f54)*zz * ratdum(ir53ng)
       dratdumdt(ir1f54)  = dratdumdt(ir54gn)*ratdum(ir53gn)*zz + ratdum(ir54gn)*dratdumdt(ir53gn)*zz - ratdum(ir1f54)*zz*denomdt
       !dratdumdd(ir1f54)  = dratdumdd(ir54gn)*ratdum(ir53gn)*zz + ratdum(ir54gn)*dratdumdd(ir53gn)*zz - ratdum(ir1f54)*zz*denomdd

       ratdum(ir2f54)     = ratdum(ir52ng)*ratdum(ir53ng)*zz
       dratdumdy1(ir2f54) = -ratdum(ir2f54)*zz * ratdum(ir53ng)
       dratdumdt(ir2f54)  = dratdumdt(ir52ng)*ratdum(ir53ng)*zz + ratdum(ir52ng)*dratdumdt(ir53ng)*zz - ratdum(ir2f54)*zz*denomdt
       !dratdumdd(ir2f54)  = dratdumdd(ir52ng)*ratdum(ir53ng)*zz + ratdum(ir52ng)*dratdumdd(ir53ng)*zz - ratdum(ir2f54)*zz*denomdd
    end if


    ! fe54(p,g)co55(p,g)ni56 equilibrium links r3f54 r4f54
    ! fe52(a,p)co55(g,p)fe54 equilibrium links r5f54 r6f54
    ! fe52(a,p)co55(p,g)ni56 equilibrium links r7f54 r8f54

    ratdum(ir3f54) = ZERO
    dratdumdy1(ir3f54) = ZERO
    dratdumdt(ir3f54) = ZERO
    !dratdumdd(ir3f54) = ZERO
    ratdum(ir4f54) = ZERO
    dratdumdy1(ir4f54) = ZERO
    dratdumdt(ir4f54) = ZERO
    !dratdumdd(ir4f54) = ZERO
    ratdum(ir5f54) = ZERO
    dratdumdy1(ir5f54) = ZERO
    dratdumdt(ir5f54) = ZERO
    !dratdumdd(ir5f54) = ZERO
    ratdum(ir6f54) = ZERO
    dratdumdy1(ir6f54) = ZERO
    dratdumdt(ir6f54) = ZERO
    !dratdumdd(ir6f54) = ZERO
    ratdum(ir7f54) = ZERO
    dratdumdy1(ir7f54) = ZERO
    dratdumdt(ir7f54) = ZERO
    !dratdumdd(ir7f54) = ZERO
    ratdum(ir8f54) = ZERO
    dratdumdy1(ir8f54) = ZERO
    dratdumdt(ir8f54) = ZERO
    !dratdumdd(ir8f54) = ZERO

    denom   = ratdum(ircogp) + y(iprot)*(ratdum(ircopg) + ratdum(ircopa))

    if (denom .gt. 1.0d-50 .and. btemp .gt. 1.5e9) then

       denomdt = dratdumdt(ircogp) + y(iprot)*(dratdumdt(ircopg) + dratdumdt(ircopa))
       !denomdd = dratdumdd(ircogp) + y(iprot)*(dratdumdd(ircopg) + dratdumdd(ircopa))

       zz      = 1.0d0/denom

       ratdum(ir3f54)     = ratdum(irfepg) * ratdum(ircopg) * zz
       dratdumdy1(ir3f54) = -ratdum(ir3f54) * zz * (ratdum(ircopg) + ratdum(ircopa))
       dratdumdt(ir3f54)  = dratdumdt(irfepg) * ratdum(ircopg) * zz &
                          + ratdum(irfepg) * dratdumdt(ircopg) * zz - ratdum(ir3f54)*zz*denomdt
       !dratdumdd(ir3f54)  = dratdumdd(irfepg) * ratdum(ircopg) * zz &
       !                    + ratdum(irfepg) * dratdumdd(ircopg) * zz - ratdum(ir3f54)*zz*denomdd

       ratdum(ir4f54)     = ratdum(irnigp) * ratdum(ircogp) * zz
       dratdumdy1(ir4f54) = -ratdum(ir4f54) * zz * (ratdum(ircopg)+ratdum(ircopa))
       dratdumdt(ir4f54)  = dratdumdt(irnigp) * ratdum(ircogp) * zz &
                          + ratdum(irnigp) * dratdumdt(ircogp) * zz - ratdum(ir4f54)*zz*denomdt
       !dratdumdd(ir4f54)  = dratdumdd(irnigp) * ratdum(ircogp) * zz &
       !                    + ratdum(irnigp) * dratdumdd(ircogp) * zz  - ratdum(ir4f54)*zz*denomdd

       ratdum(ir5f54)     = ratdum(irfepg) * ratdum(ircopa) * zz
       dratdumdy1(ir5f54) = -ratdum(ir5f54) * zz * (ratdum(ircopg)+ratdum(ircopa))
       dratdumdt(ir5f54)  = dratdumdt(irfepg) * ratdum(ircopa) * zz &
                          + ratdum(irfepg) * dratdumdt(ircopa) * zz - ratdum(ir5f54) * zz * denomdt
       !dratdumdd(ir5f54)  = dratdumdd(irfepg) * ratdum(ircopa) * zz &
       !                    + ratdum(irfepg) * dratdumdd(ircopa) * zz - ratdum(ir5f54) * zz * denomdd

       ratdum(ir6f54)     = ratdum(irfeap) * ratdum(ircogp) * zz
       dratdumdy1(ir6f54) = -ratdum(ir6f54) * zz * (ratdum(ircopg)+ratdum(ircopa))
       dratdumdt(ir6f54)  = dratdumdt(irfeap) * ratdum(ircogp) * zz &
                          + ratdum(irfeap) * dratdumdt(ircogp) * zz - ratdum(ir6f54) * zz * denomdt
       !dratdumdd(ir6f54)  = dratdumdd(irfeap) * ratdum(ircogp) * zz &
       !                    + ratdum(irfeap) * dratdumdd(ircogp) * zz - ratdum(ir6f54) * zz * denomdd

       ratdum(ir7f54)     = ratdum(irfeap) * ratdum(ircopg) * zz
       dratdumdy1(ir7f54) = -ratdum(ir7f54) * zz * (ratdum(ircopg)+ratdum(ircopa))
       dratdumdt(ir7f54)  = dratdumdt(irfeap) * ratdum(ircopg) * zz &
                          + ratdum(irfeap) * dratdumdt(ircopg) * zz - ratdum(ir7f54) * zz * denomdt
       !dratdumdd(ir7f54)  = dratdumdd(irfeap) * ratdum(ircopg) * zz &
       !                   + ratdum(irfeap) * dratdumdd(ircopg) * zz - ratdum(ir7f54) * zz * denomdd

       ratdum(ir8f54)     = ratdum(irnigp) * ratdum(ircopa) * zz
       dratdumdy1(ir8f54) = -ratdum(ir8f54) * zz * (ratdum(ircopg)+ratdum(ircopa))
       dratdumdt(ir8f54)  = dratdumdt(irnigp) * ratdum(ircopa) * zz &
                          + ratdum(irnigp) * dratdumdt(ircopa) * zz - ratdum(ir8f54) * zz * denomdt
       !dratdumdd(ir8f54)  = dratdumdd(irnigp) * ratdum(ircopa) * zz &
       !                   + ratdum(irnigp) * dratdumdd(ircopa) * zz - ratdum(ir8f54) * zz * denomdd
    end if



    ! p(n,g)h2(n,g)3h(p,g)he4   photodisintegrated n and p back to he4 equilibrium links
    ! p(n,g)h2(p,g)he3(n,g)he4

    ratdum(iralf1)     = ZERO
    dratdumdy1(iralf1) = ZERO
    dratdumdy2(iralf1) = ZERO
    dratdumdt(iralf1)  = ZERO
    !dratdumdd(iralf1)  = ZERO
    ratdum(iralf2)     = ZERO
    dratdumdy1(iralf2) = ZERO
    dratdumdy2(iralf2) = ZERO
    dratdumdt(iralf2)  = ZERO
    !dratdumdd(iralf2)  = ZERO

    denom  = ratdum(irhegp)*ratdum(irdgn) + y(ineut)*ratdum(irheng)*ratdum(irdgn) + y(ineut)*y(iprot)*ratdum(irheng)*ratdum(irdpg)

    if (denom .gt. 1.0d-50 .and. btemp .gt. 1.5e9) then

       denomdt  = dratdumdt(irhegp)*ratdum(irdgn) + ratdum(irhegp)*dratdumdt(irdgn) &
            +  y(ineut) * (dratdumdt(irheng)*ratdum(irdgn) + ratdum(irheng)*dratdumdt(irdgn)) &
            +  y(ineut)*y(iprot) * (dratdumdt(irheng)*ratdum(irdpg) + ratdum(irheng)*dratdumdt(irdpg))

       !denomdd  = dratdumdd(irhegp)*ratdum(irdgn) + ratdum(irhegp)*dratdumdd(irdgn) &
       !     +  y(ineut) * (dratdumdd(irheng)*ratdum(irdgn) + ratdum(irheng)*dratdumdd(irdgn)) &
       !     +  y(ineut)*y(iprot) * (dratdumdd(irheng)*ratdum(irdpg) + ratdum(irheng)*dratdumdd(irdpg))

       zz = 1.0d0/denom

       ratdum(iralf1)     = ratdum(irhegn) * ratdum(irhegp) * ratdum(irdgn) * zz
       dratdumdy1(iralf1) = -ratdum(iralf1) * zz * (ratdum(irheng)*ratdum(irdgn) + y(iprot)*ratdum(irheng)*ratdum(irdpg))
       dratdumdy2(iralf1) = -ratdum(iralf1) * zz * y(ineut) * ratdum(irheng) * ratdum(irdpg)
       dratdumdt(iralf1)  = dratdumdt(irhegn)*ratdum(irhegp) * ratdum(irdgn) * zz &
            + ratdum(irhegn)*dratdumdt(irhegp)*ratdum(irdgn)*zz &
            + ratdum(irhegn)*ratdum(irhegp)*dratdumdt(irdgn)*zz &
            - ratdum(iralf1)*zz*denomdt
       !dratdumdd(iralf1)  = dratdumdd(irhegn) * ratdum(irhegp) * ratdum(irdgn) * zz &
       !     + ratdum(irhegn)*dratdumdd(irhegp)*ratdum(irdgn)*zz &
       !     + ratdum(irhegn)*ratdum(irhegp)*dratdumdd(irdgn)*zz &
       !     - ratdum(iralf1)*zz*denomdt


       ratdum(iralf2)     = ratdum(irheng)*ratdum(irdpg) * ratdum(irhng)*zz
       dratdumdy1(iralf2) = -ratdum(iralf2) * zz * (ratdum(irheng)*ratdum(irdgn) + y(iprot)*ratdum(irheng)*ratdum(irdpg))
       dratdumdy2(iralf2) = -ratdum(iralf2) * zz * y(ineut) * ratdum(irheng) * ratdum(irdpg)
       dratdumdt(iralf2)  = dratdumdt(irheng)*ratdum(irdpg) * ratdum(irhng) * zz &
            + ratdum(irheng)*dratdumdt(irdpg)*ratdum(irhng)*zz &
            + ratdum(irheng)*ratdum(irdpg)*dratdumdt(irhng)*zz &
            - ratdum(iralf2)*zz*denomdt
       !dratdumdd(iralf2)  = dratdumdd(irheng)*ratdum(irdpg) * ratdum(irhng)*zz &
       !     + ratdum(irheng)*dratdumdd(irdpg)*ratdum(irhng)*zz &
       !     + ratdum(irheng)*ratdum(irdpg)*dratdumdd(irhng)*zz &
       !     - ratdum(iralf2)*zz*denomdd
    end if



    ! he3(a,g)be7(p,g)8b(e+nu)8be(2a)
    ! beta limit he3+he4 by the 8B decay half life

    if (y(ihe4) .gt. 1.0d-30) then
       xx            = 0.896d0/y(ihe4)
       ratdum(irhe3ag)  = min(ratdum(irhe3ag),xx)
       if (ratdum(irhe3ag) .eq. xx) then
          dratdumdy1(irhe3ag) = -xx/y(ihe4)
          dratdumdt(irhe3ag)  = ZERO
          !dratdumdd(irhe3ag)  = ZERO
       else
          dratdumdy1(irhe3ag) = ZERO
       endif
    endif


    ! beta limit n14(p,g)o15(enu)o16  and o16(p,g)f17(e+nu)17o(p,a)n14

    if (y(ih1)  .gt. 1.0d-30) then
       xx = 5.68d-03/(y(ih1)*1.57d0)
       ratdum(irnpg) = min(ratdum(irnpg),xx)
       if (ratdum(irnpg) .eq. xx) then
          dratdumdy1(irnpg) = -xx/y(ih1)
          dratdumdt(irnpg)  = ZERO
          !dratdumdd(irnpg)  = ZERO
       else
          dratdumdy1(irnpg) = ZERO
       end if

       xx = 0.0105d0/y(ih1)
       ratdum(iropg) = min(ratdum(iropg),xx)
       if (ratdum(iropg) .eq. xx) then
          dratdumdy1(iropg) = -xx/y(ih1)
          dratdumdt(iropg)  = ZERO
          !dratdumdd(iropg)  = ZERO
       else
          dratdumdy1(iropg) = ZERO
       end if
    end if

  end subroutine screen_aprox19



  subroutine dfdy_isotopes_aprox19(y,dfdy,ratdum,dratdumdy1,dratdumdy2)

    use network
    use microphysics_math_module, only: esum

    implicit none

    ! this routine sets up the dense aprox19 jacobian for the isotopes

    double precision :: y(nspec), dfdy(nspec,nspec)
    double precision :: ratdum(nrates), dratdumdy1(nrates), dratdumdy2(nrates)

    double precision :: b(30)

    ! h1 jacobian elements
    ! d(h1)/d(h1)
    b(1) = -3.0d0 * y(ih1)  * ratdum(irpp)
    b(2) = -2.0d0 * y(ic12) * ratdum(ircpg)
    b(3) = -2.0d0 * y(in14) * ratdum(irnpg)
    b(4) = -2.0d0 * y(in14) * y(ih1) * dratdumdy1(irnpg)
    b(5) = -2.0d0 * y(io16) * ratdum(iropg)
    b(6) = -2.0d0 * y(io16) * y(ih1) * dratdumdy1(iropg)
    b(7) = -3.0d0 * ratdum(irpen)
    dfdy(ih1,ih1) = esum(b,7)

    ! d(h1)/d(he3)
    b(1) = 2.0d0 * y(ihe3) * ratdum(ir33)
    b(2) = -y(ihe4) * ratdum(irhe3ag)
    dfdy(ih1,ihe3) = sum(b(1:2))

    ! d(h1)/d(he4)
    b(1) = -y(ihe3) * ratdum(irhe3ag)
    b(2) = -y(ihe3) * y(ihe4) * dratdumdy1(irhe3ag)
    dfdy(ih1,ihe4) = sum(b(1:2))

    ! d(h1)/d(c12)
    dfdy(ih1,ic12) = -2.0d0 * y(ih1) * ratdum(ircpg)

    ! d(h1)/d(n14)
    dfdy(ih1,in14) = -2.0d0 * y(ih1) * ratdum(irnpg)

    ! d(h1)/d(o16)
    dfdy(ih1,io16) = -2.0d0 * y(ih1) * ratdum(iropg)


    ! he3 jacobian elements
    ! d(he3)/d(h1)
    b(1) = y(ih1) * ratdum(irpp)
    b(2) = ratdum(irpen)
    dfdy(ihe3,ih1) = sum(b(1:2))

    ! d(he3)/d(he3)
    b(1) = -2.0d0 * y(ihe3) * ratdum(ir33)
    b(2) = -y(ihe4) * ratdum(irhe3ag)
    dfdy(ihe3,ihe3) = sum(b(1:2))

    ! d(he3)/d(he4)
    b(1) = -y(ihe3) * ratdum(irhe3ag)
    b(2) = -y(ihe3) * y(ihe4) * dratdumdy1(irhe3ag)
    dfdy(ihe3,ihe4) = sum(b(1:2))


    ! he4 jacobian elements
    ! d(he4)/d(h1)
    b(1) =  y(in14) * ratdum(ifa) * ratdum(irnpg)
    b(2) =  y(in14) * y(ih1) * ratdum(ifa) * dratdumdy1(irnpg)
    b(3) =  y(io16) * ratdum(iropg)
    b(4) =  y(io16) * y(ih1) * dratdumdy1(iropg)
    dfdy(ihe4,ih1) = esum(b,4)

    ! d(he4)/d(he3)
    b(1) = y(ihe3) * ratdum(ir33)
    b(2) = y(ihe4) * ratdum(irhe3ag)
    dfdy(ihe4,ihe3) = sum(b(1:2))


    ! d(he4)/d(he4)
    b(1)  = -1.5d0 * y(ihe4) * y(ihe4) * ratdum(ir3a)
    b(2)  = -y(ic12)  * ratdum(ircag)
    b(3)  = -y(io16)  * ratdum(iroag)
    b(4)  = -y(ine20) * ratdum(irneag)
    b(5)  = -y(img24) * ratdum(irmgag)
    b(6)  = -y(isi28) * ratdum(irsiag)
    b(7)  = -y(is32)  * ratdum(irsag)
    b(8)  = -y(iar36) * ratdum(irarag)
    b(9)  = -y(ica40) * ratdum(ircaag)
    b(10) = -y(iti44) * ratdum(irtiag)
    b(11) = -y(icr48) * ratdum(ircrag)
    b(12) = -y(ife52) * ratdum(irfeag)
    b(13) = -y(img24) * ratdum(irmgap) * (1.0d0-ratdum(irr1))
    b(14) = -y(isi28) * ratdum(irsiap) * (1.0d0-ratdum(irs1))
    b(15) = -y(is32)  * ratdum(irsap)  * (1.0d0-ratdum(irt1))
    b(16) = -y(iar36) * ratdum(irarap) * (1.0d0-ratdum(iru1))
    b(17) = -y(ica40) * ratdum(ircaap) * (1.0d0-ratdum(irv1))
    b(18) = -y(iti44) * ratdum(irtiap) * (1.0d0-ratdum(irw1))
    b(19) = -y(icr48) * ratdum(ircrap) * (1.0d0-ratdum(irx1))
    b(20) = -y(ife52) * ratdum(ir6f54)
    b(21) = -y(ife52) * y(iprot) * ratdum(ir7f54)
    b(22) = -ratdum(iralf1)
    b(23) =  y(ihe3) * ratdum(irhe3ag)
    b(24) =  y(ihe3) * y(ihe4) * dratdumdy1(irhe3ag)
    b(25) = -y(in14) * ratdum(irnag) * 1.5d0
    dfdy(ihe4,ihe4) = esum(b,25)

    ! d(he4)/d(c12)
    b(1) =  y(ic12) * ratdum(ir1212)
    b(2) =  0.5d0 * y(io16) * ratdum(ir1216)
    b(3) =  3.0d0 * ratdum(irg3a)
    b(4) = -y(ihe4) * ratdum(ircag)
    dfdy(ihe4,ic12) = esum(b,4)

    ! d(he4)/d(n14)
    b(1) =  y(ih1) * ratdum(ifa) * ratdum(irnpg)
    b(2) = -y(ihe4) * ratdum(irnag) * 1.5d0
    dfdy(ihe4,in14) = sum(b(1:2))

    ! d(he4)/d(o16)
    b(1) =  0.5d0 * y(ic12) * ratdum(ir1216)
    b(2) =  1.12d0 * 0.5d0*y(io16) * ratdum(ir1616)
    b(3) =  0.68d0 * ratdum(irs1) * 0.5d0*y(io16) * ratdum(ir1616)
    b(4) =  ratdum(iroga)
    b(5) = -y(ihe4) * ratdum(iroag)
    b(6) =  y(ih1) * ratdum(iropg)
    dfdy(ihe4,io16) = esum(b,6)

    ! d(he4)/d(ne20)
    b(1) =  ratdum(irnega)
    b(2) = -y(ihe4) * ratdum(irneag)
    dfdy(ihe4,ine20) = sum(b(1:2))

    ! d(he4)/d(mg24)
    b(1) =  ratdum(irmgga)
    b(2) = -y(ihe4) * ratdum(irmgag)
    b(3) = -y(ihe4) * ratdum(irmgap) * (1.0d0-ratdum(irr1))
    dfdy(ihe4,img24) = esum(b,3)

    ! d(he4)/d(si28)
    b(1) =  ratdum(irsiga)
    b(2) = -y(ihe4) * ratdum(irsiag)
    b(3) = -y(ihe4) * ratdum(irsiap) * (1.0d0-ratdum(irs1))
    b(4) =  ratdum(irr1) * ratdum(irsigp)
    dfdy(ihe4,isi28) = esum(b,4)

    ! d(he4)/d(s32)
    b(1) =  ratdum(irsga)
    b(2) = -y(ihe4) * ratdum(irsag)
    b(3) = -y(ihe4) * ratdum(irsap) * (1.0d0-ratdum(irt1))
    b(4) =  ratdum(irs1) * ratdum(irsgp)
    dfdy(ihe4,is32) = esum(b,4)

    ! d(he4)/d(ar36)
    b(1)  =  ratdum(irarga)
    b(2)  = -y(ihe4) * ratdum(irarag)
    b(3)  = -y(ihe4) * ratdum(irarap) * (1.0d0-ratdum(iru1))
    b(4)  =  ratdum(irt1) * ratdum(irargp)
    dfdy(ihe4,iar36) = esum(b,4)

    ! d(he4)/d(ca40)
    b(1) =  ratdum(ircaga)
    b(2) = -y(ihe4) * ratdum(ircaag)
    b(3) = -y(ihe4) * ratdum(ircaap) * (1.0d0-ratdum(irv1))
    b(4) =  ratdum(iru1) * ratdum(ircagp)
    dfdy(ihe4,ica40) = esum(b,4)

    ! d(he4)/d(ti44)
    b(1) =  ratdum(irtiga)
    b(2) = -y(ihe4) * ratdum(irtiag)
    b(3) = -y(ihe4) * ratdum(irtiap) * (1.0d0-ratdum(irw1))
    b(4) =  ratdum(irv1) * ratdum(irtigp)
    dfdy(ihe4,iti44) = esum(b,4)

    ! d(he4)/d(cr48)
    b(1) =  ratdum(ircrga)
    b(2) = -y(ihe4) * ratdum(ircrag)
    b(3) = -y(ihe4) * ratdum(ircrap) * (1.0d0-ratdum(irx1))
    b(4) =  ratdum(irw1) * ratdum(ircrgp)
    dfdy(ihe4,icr48) = esum(b,4)

    ! d(he4)/d(fe52)
    b(1) =  ratdum(irfega)
    b(2) = -y(ihe4) * ratdum(irfeag)
    b(3) =  ratdum(irx1) * ratdum(irfegp)
    b(4) = -y(ihe4) * ratdum(ir6f54)
    b(5) = -y(ihe4) * y(iprot) * ratdum(ir7f54)
    dfdy(ihe4,ife52) = esum(b,5)

    ! d(he4)/d(fe54)
    dfdy(ihe4,ife54) = y(iprot) * y(iprot) * ratdum(ir5f54)

    ! d(he4)/d(ni56)
    b(1) = ratdum(irniga)
    b(2) = y(iprot) * ratdum(ir8f54)
    dfdy(ihe4,ini56) = sum(b(1:2))

    ! d(he4)/d(neut)
    b(1) = -y(ihe4) * dratdumdy1(iralf1)
    b(2) =  2.0d0 * y(ineut) * y(iprot)*y(iprot) * ratdum(iralf2)
    b(3) =  y(ineut)*y(ineut) * y(iprot)*y(iprot) * dratdumdy1(iralf2)
    dfdy(ihe4,ineut) = esum(b,3)

    ! d(he4)/d(prot)
    b(1)  =  2.0d0 * y(ife54) * y(iprot) * ratdum(ir5f54)
    b(2)  =  y(ife54) * y(iprot) * y(iprot) * dratdumdy1(ir5f54)
    b(3)  = -y(ihe4) * y(ife52) * dratdumdy1(ir6f54)
    b(4)  = -y(ife52) * y(ihe4) * ratdum(ir7f54)
    b(5)  = -y(ife52) * y(ihe4) * y(iprot) * dratdumdy1(ir7f54)
    b(6)  =  y(ini56) * ratdum(ir8f54)
    b(7)  =  y(ini56) * y(iprot) * dratdumdy1(ir8f54)
    b(8)  = -y(ihe4) * dratdumdy2(iralf1)
    b(9)  =  2.0d0 * y(ineut)*y(ineut) * y(iprot) * ratdum(iralf2)
    b(10) =  y(ineut)*y(ineut) * y(iprot)*y(iprot) * dratdumdy2(iralf2)
    dfdy(ihe4,iprot) = esum(b,10)


    ! c12 jacobian elements
    ! d(c12)/d(h1)
    b(1) = -y(ic12) * ratdum(ircpg)
    b(2) =  y(in14) * ratdum(ifa) * ratdum(irnpg)
    b(3) =  y(in14) * y(ih1) * ratdum(ifa) * dratdumdy1(irnpg)
    dfdy(ic12,ih1) = esum(b,3)

    ! d(c12)/d(he4)
    b(1) =  0.5d0 * y(ihe4) * y(ihe4) * ratdum(ir3a)
    b(2) = -y(ic12) * ratdum(ircag)
    dfdy(ic12,ihe4) = sum(b(1:2))

    ! d(c12)/d(c12)
    b(1) = -2.0d0 * y(ic12) * ratdum(ir1212)
    b(2) = -y(io16) * ratdum(ir1216)
    b(3) = -ratdum(irg3a)
    b(4) = -y(ihe4) * ratdum(ircag)
    b(5) = -y(ih1) * ratdum(ircpg)
    dfdy(ic12,ic12) = esum(b,5)

    ! d(c12)/d(n14)
    dfdy(ic12,in14) = y(ih1) * ratdum(ifa) * ratdum(irnpg)

    ! d(c12)/d(o16)
    b(1) = -y(ic12) * ratdum(ir1216)
    b(2) =  ratdum(iroga)
    dfdy(ic12,io16) = sum(b(1:2))


    ! n14 jacobian elements
    ! d(n14)/d(h1)
    b(1) =   y(ic12) * ratdum(ircpg)
    b(2) =  -y(in14) * ratdum(irnpg)
    b(3) =  -y(in14) * y(ih1) * dratdumdy1(irnpg)
    b(4) =   y(io16) * ratdum(iropg)
    b(5) =   y(io16) * y(ih1) * dratdumdy1(iropg)
    dfdy(in14,ih1) = esum(b,5)

    ! d(n14)/d(he4)
    dfdy(in14,ihe4) = -y(in14) * ratdum(irnag)

    ! d(n14)/d(c12)
    dfdy(in14,ic12) =  y(ih1) * ratdum(ircpg)

    ! d(n14)/d(n14)
    b(1) = -y(ih1) * ratdum(irnpg)
    b(2) = -y(ihe4) * ratdum(irnag)
    dfdy(in14,in14) = sum(b(1:2))

    ! d(n14)/d(o16)
    dfdy(in14,io16) = y(ih1) * ratdum(iropg)


    ! o16 jacobian elements
    ! d(o16)/d(h1)
    b(1) =  y(in14) * ratdum(ifg) * ratdum(irnpg)
    b(2) =  y(in14) * y(ih1) * ratdum(ifg) * dratdumdy1(irnpg)
    b(3) = -y(io16) * ratdum(iropg)
    b(4) = -y(io16) * y(ih1) * dratdumdy1(iropg)
    dfdy(io16,ih1) = esum(b,4)

    ! d(o16)/d(he4)
    b(1) =  y(ic12)*ratdum(ircag)
    b(2) = -y(io16)*ratdum(iroag)
    dfdy(io16,ihe4) = sum(b(1:2))

    ! d(o16)/d(c12)
    b(1) = -y(io16)*ratdum(ir1216)
    b(2) =  y(ihe4)*ratdum(ircag)
    dfdy(io16,ic12) = sum(b(1:2))

    ! d(o16)/d(n14)
    dfdy(io16,in14) = y(ih1) * ratdum(ifg) * ratdum(irnpg)

    ! d(o16)/d(o16)
    b(1) = -y(ic12) * ratdum(ir1216)
    b(2) = -2.0d0 * y(io16) * ratdum(ir1616)
    b(3) = -y(ihe4) * ratdum(iroag)
    b(4) = -ratdum(iroga)
    b(5) = -y(ih1) * ratdum(iropg)
    dfdy(io16,io16) = esum(b,5)

    ! d(o16)/d(ne20)
    dfdy(io16,ine20) = ratdum(irnega)


    ! ne20 jacobian elements
    ! d(ne20)/d(he4)
    b(1) =  y(io16) * ratdum(iroag)
    b(2) = -y(ine20) * ratdum(irneag)
    b(3) =  y(in14) * ratdum(irnag)
    dfdy(ine20,ihe4) = esum(b,3)

    ! d(ne20)/d(c12)
    dfdy(ine20,ic12) = y(ic12) * ratdum(ir1212)

    ! d(ne20)/d(n14)
    dfdy(ine20,in14) = y(ihe4) * ratdum(irnag)

    ! d(ne20)/d(o16)
    dfdy(ine20,io16) = y(ihe4) * ratdum(iroag)

    ! d(ne20)/d(ne20)
    b(1) = -y(ihe4) * ratdum(irneag)
    b(2) = -ratdum(irnega)
    dfdy(ine20,ine20) = sum(b(1:2))

    ! d(ne20)/d(mg24)
    dfdy(ine20,img24) = ratdum(irmgga)


    ! mg24 jacobian elements
    ! d(mg24)/d(he4)
    b(1) =  y(ine20) * ratdum(irneag)
    b(2) = -y(img24) * ratdum(irmgag)
    b(3) = -y(img24) * ratdum(irmgap) * (1.0d0-ratdum(irr1))
    dfdy(img24,ihe4) = esum(b,3)

    ! d(mg24)/d(c12)
    dfdy(img24,ic12) = 0.5d0 * y(io16) * ratdum(ir1216)

    ! d(mg24)/d(o16)
    dfdy(img24,io16) = 0.5d0 * y(ic12) * ratdum(ir1216)

    ! d(mg24)/d(ne20)
    dfdy(img24,ine20) = y(ihe4) * ratdum(irneag)

    ! d(mg24)/d(mg24)
    b(1) = -y(ihe4) * ratdum(irmgag)
    b(2) = -ratdum(irmgga)
    b(3) = -y(ihe4) * ratdum(irmgap) * (1.0d0-ratdum(irr1))
    dfdy(img24,img24) = esum(b,3)

    ! d(mg24)/d(si28)
    b(1) = ratdum(irsiga)
    b(2) = ratdum(irr1) * ratdum(irsigp)
    dfdy(img24,isi28) = sum(b(1:2))


    ! si28 jacobian elements
    ! d(si28)/d(he4)
    b(1) =  y(img24) * ratdum(irmgag)
    b(2) = -y(isi28) * ratdum(irsiag)
    b(3) =  y(img24) * ratdum(irmgap) * (1.0d0-ratdum(irr1))
    b(4) = -y(isi28) * ratdum(irsiap) * (1.0d0-ratdum(irs1))
    dfdy(isi28,ihe4) = esum(b,4)

    ! d(si28)/d(c12)
    dfdy(isi28,ic12) =  0.5d0 * y(io16) * ratdum(ir1216)

    ! d(si28)/d(o16)
    b(1) = 0.5d0 * y(ic12) * ratdum(ir1216)
    b(2) = 1.12d0 * 0.5d0*y(io16) * ratdum(ir1616)
    b(3) = 0.68d0 * 0.5d0*y(io16) * ratdum(irs1) * ratdum(ir1616)
    dfdy(isi28,io16) = esum(b,3)

    ! d(si28)/d(mg24)
    b(1) =  y(ihe4) * ratdum(irmgag)
    b(2) =  y(ihe4) * ratdum(irmgap) * (1.0d0-ratdum(irr1))
    dfdy(isi28,img24) = sum(b(1:2))

    ! d(si28)/d(si28)
    b(1) =  -y(ihe4) * ratdum(irsiag)
    b(2) = -ratdum(irsiga)
    b(3) = -ratdum(irr1) * ratdum(irsigp)
    b(4) = -y(ihe4) * ratdum(irsiap) * (1.0d0-ratdum(irs1))
    dfdy(isi28,isi28) = esum(b,4)

    ! d(si28)/d(s32)
    b(1) = ratdum(irsga)
    b(2) = ratdum(irs1) * ratdum(irsgp)
    dfdy(isi28,is32) = sum(b(1:2))


    ! s32 jacobian elements
    ! d(s32)/d(he4)
    b(1) =  y(isi28) * ratdum(irsiag)
    b(2) = -y(is32) * ratdum(irsag)
    b(3) =  y(isi28) * ratdum(irsiap) * (1.0d0-ratdum(irs1))
    b(4) = -y(is32) * ratdum(irsap) * (1.0d0-ratdum(irt1))
    dfdy(is32,ihe4) = esum(b,4)

    ! d(s32)/d(o16)
    b(1) = 0.68d0*0.5d0*y(io16)*ratdum(ir1616)*(1.0d0-ratdum(irs1))
    b(2) = 0.2d0 * 0.5d0*y(io16) * ratdum(ir1616)
    dfdy(is32,io16) = sum(b(1:2))

    ! d(s32)/d(si28)
    b(1)  =y(ihe4) * ratdum(irsiag)
    b(2) = y(ihe4) * ratdum(irsiap) * (1.0d0-ratdum(irs1))
    dfdy(is32,isi28) = sum(b(1:2))

    ! d(s32)/d(s32)
    b(1) = -y(ihe4) * ratdum(irsag)
    b(2) = -ratdum(irsga)
    b(3) = -ratdum(irs1) * ratdum(irsgp)
    b(4) = -y(ihe4) * ratdum(irsap) * (1.0d0-ratdum(irt1))
    dfdy(is32,is32) = esum(b,4)

    ! d(s32)/d(ar36)
    b(1) = ratdum(irarga)
    b(2) = ratdum(irt1) * ratdum(irargp)
    dfdy(is32,iar36) = sum(b(1:2))


    ! ar36 jacobian elements
    ! d(ar36)/d(he4)
    b(1) =  y(is32)  * ratdum(irsag)
    b(2) = -y(iar36) * ratdum(irarag)
    b(3) =  y(is32)  * ratdum(irsap) * (1.0d0-ratdum(irt1))
    b(4) = -y(iar36) * ratdum(irarap) * (1.0d0-ratdum(iru1))
    dfdy(iar36,ihe4) = esum(b,4)

    ! d(ar36)/d(s32)
    b(1) = y(ihe4) * ratdum(irsag)
    b(2) = y(ihe4) * ratdum(irsap) * (1.0d0-ratdum(irt1))
    dfdy(iar36,is32) = sum(b(1:2))

    ! d(ar36)/d(ar36)
    b(1) = -y(ihe4) * ratdum(irarag)
    b(2) = -ratdum(irarga)
    b(3) = -ratdum(irt1) * ratdum(irargp)
    b(4) = -y(ihe4) * ratdum(irarap) * (1.0d0-ratdum(iru1))
    dfdy(iar36,iar36) = esum(b,4)

    ! d(ar36)/d(ca40)
    b(1) = ratdum(ircaga)
    b(2) = ratdum(ircagp) * ratdum(iru1)
    dfdy(iar36,ica40) = sum(b(1:2))


    ! ca40 jacobian elements
    ! d(ca40)/d(he4)
    b(1)  =  y(iar36) * ratdum(irarag)
    b(2)  = -y(ica40) * ratdum(ircaag)
    b(3)  =  y(iar36) * ratdum(irarap)*(1.0d0-ratdum(iru1))
    b(4)  = -y(ica40) * ratdum(ircaap)*(1.0d0-ratdum(irv1))
    dfdy(ica40,ihe4) = esum(b,4)

    ! d(ca40)/d(ar36)
    b(1) =  y(ihe4) * ratdum(irarag)
    b(2) =  y(ihe4) * ratdum(irarap)*(1.0d0-ratdum(iru1))
    dfdy(ica40,iar36) = sum(b(1:2))

    ! d(ca40)/d(ca40)
    b(1) = -y(ihe4) * ratdum(ircaag)
    b(2) = -ratdum(ircaga)
    b(3) = -ratdum(ircagp) * ratdum(iru1)
    b(4) = -y(ihe4) * ratdum(ircaap)*(1.0d0-ratdum(irv1))
    dfdy(ica40,ica40) = esum(b,4)

    ! d(ca40)/d(ti44)
    b(1) = ratdum(irtiga)
    b(2) = ratdum(irtigp) * ratdum(irv1)
    dfdy(ica40,iti44) = sum(b(1:2))


    ! ti44 jacobian elements
    ! d(ti44)/d(he4)
    b(1) =  y(ica40) * ratdum(ircaag)
    b(2) = -y(iti44) * ratdum(irtiag)
    b(3) =  y(ica40) * ratdum(ircaap)*(1.0d0-ratdum(irv1))
    b(4) = -y(iti44) * ratdum(irtiap)*(1.0d0-ratdum(irw1))
    dfdy(iti44,ihe4) = esum(b,4)

    ! d(ti44)/d(ca40)
    b(1) =  y(ihe4) * ratdum(ircaag)
    b(2) =  y(ihe4) * ratdum(ircaap)*(1.0d0-ratdum(irv1))
    dfdy(iti44,ica40) = sum(b(1:2))

    ! d(ti44)/d(ti44)
    b(1) = -y(ihe4) * ratdum(irtiag)
    b(2) = -ratdum(irtiga)
    b(3) = -ratdum(irv1) * ratdum(irtigp)
    b(4) = -y(ihe4) * ratdum(irtiap)*(1.0d0-ratdum(irw1))
    dfdy(iti44,iti44) = esum(b,4)

    ! d(ti44)/d(cr48)
    b(1) = ratdum(ircrga)
    b(2) = ratdum(irw1) * ratdum(ircrgp)
    dfdy(iti44,icr48) = sum(b(1:2))


    ! cr48 jacobian elements
    ! d(cr48)/d(he4)
    b(1) =  y(iti44) * ratdum(irtiag)
    b(2) = -y(icr48) * ratdum(ircrag)
    b(3) =  y(iti44) * ratdum(irtiap)*(1.0d0-ratdum(irw1))
    b(4) = -y(icr48) * ratdum(ircrap)*(1.0d0-ratdum(irx1))
    dfdy(icr48,ihe4) = esum(b,4)

    ! d(cr48)/d(ti44)
    b(1) =  y(ihe4) * ratdum(irtiag)
    b(2) =  y(ihe4) * ratdum(irtiap)*(1.0d0-ratdum(irw1))
    dfdy(icr48,iti44) = sum(b(1:2))

    ! d(cr48)/d(cr48)
    b(1) = -y(ihe4) * ratdum(ircrag)
    b(2) = -ratdum(ircrga)
    b(3) = -ratdum(irw1) * ratdum(ircrgp)
    b(4) = -y(ihe4) * ratdum(ircrap)*(1.0d0-ratdum(irx1))
    dfdy(icr48,icr48) = esum(b,4)

    ! d(cr48)/d(fe52)
    b(1) = ratdum(irfega)
    b(2) = ratdum(irx1) * ratdum(irfegp)
    dfdy(icr48,ife52) = sum(b(1:2))


    ! fe52 jacobian elements
    ! d(fe52)/d(he4)
    b(1) =  y(icr48) * ratdum(ircrag)
    b(2) = -y(ife52) * ratdum(irfeag)
    b(3) =  y(icr48) * ratdum(ircrap) * (1.0d0-ratdum(irx1))
    b(4) = -y(ife52) * ratdum(ir6f54)
    b(5) = -y(ife52) * y(iprot) * ratdum(ir7f54)
    dfdy(ife52,ihe4) = esum(b,5)

    ! d(fe52)/d(cr48)
    b(1) = y(ihe4) * ratdum(ircrag)
    b(2) = y(ihe4) * ratdum(ircrap) * (1.0d0-ratdum(irx1))
    dfdy(ife52,icr48) = sum(b(1:2))

    ! d(fe52)/d(fe52)
    b(1) = -y(ihe4) * ratdum(irfeag)
    b(2) = -ratdum(irfega)
    b(3) = -ratdum(irx1) * ratdum(irfegp)
    b(4) = -y(ineut) * y(ineut) * ratdum(ir2f54)
    b(5) = -y(ihe4) * ratdum(ir6f54)
    b(6) = -y(ihe4) * y(iprot) * ratdum(ir7f54)
    dfdy(ife52,ife52) = esum(b,6)

    ! d(fe52)/d(fe54)
    b(1) = ratdum(ir1f54)
    b(2) = y(iprot) * y(iprot) * ratdum(ir5f54)
    dfdy(ife52,ife54) = sum(b(1:2))

    ! d(fe52)/d(ni56)
    b(1) = ratdum(irniga)
    b(2) = y(iprot) * ratdum(ir8f54)
    dfdy(ife52,ini56) = sum(b(1:2))

    ! d(fe52)/d(neut)
    b(1) =  y(ife54) * dratdumdy1(ir1f54)
    b(2) = -2.0d0 * y(ife52) * y(ineut) * ratdum(ir2f54)
    b(3) = -y(ife52) * y(ineut) * y(ineut) * dratdumdy1(ir2f54)
    dfdy(ife52,ineut) = esum(b,3)

    ! d(fe52)/d(prot)
    b(1) =  2.0d0 * y(ife54) * y(iprot) * ratdum(ir5f54)
    b(2) =  y(ife54) * y(iprot) * y(iprot) * dratdumdy1(ir5f54)
    b(3) = -y(ihe4) * y(ife52) * dratdumdy1(ir6f54)
    b(4) = -y(ife52) * y(ihe4) * ratdum(ir7f54)
    b(5) = -y(ife52) * y(ihe4) * y(iprot) * dratdumdy1(ir7f54)
    b(6) =  y(ini56) * ratdum(ir8f54)
    b(7) =  y(ini56) * y(iprot) * dratdumdy1(ir8f54)
    dfdy(ife52,iprot) = esum(b,7)


    ! fe54 jacobian elements
    ! d(fe54)/d(he4)
    dfdy(ife52,ihe4) = y(ife52) * ratdum(ir6f54)

    ! d(fe54)/d(fe52)
    b(1) =  y(ineut) * y(ineut) * ratdum(ir2f54)
    b(2) =  y(ihe4) * ratdum(ir6f54)
    dfdy(ife54,ife52) = sum(b(1:2))

    ! d(fe54)/d(fe54)
    b(1) = -ratdum(ir1f54)
    b(2) = -y(iprot) * y(iprot) * ratdum(ir3f54)
    b(3) = -y(iprot) * y(iprot) * ratdum(ir5f54)
    dfdy(ife54,ife54) = esum(b,3)

    ! d(fe54)/d(ni56)
    b(1)  = ratdum(ir4f54)
    b(2)  = c54 * ratdum(irn56ec)
    dfdy(ife54,ini56) = sum(b(1:2))

    ! d(fe54)/d(neut)
    b(1) = -y(ife54) * dratdumdy1(ir1f54)
    b(2) =  2.0d0 * y(ife52) * y(ineut) * ratdum(ir2f54)
    b(3) =  y(ife52) * y(ineut) * y(ineut) * dratdumdy1(ir2f54)
    dfdy(ife52,ineut) = esum(b,3)

    ! d(fe54)/d(prot)
    b(1) = -2.0d0 * y(ife54) * y(iprot) * ratdum(ir3f54)
    b(2) = -y(ife54) * y(iprot) * y(iprot) * dratdumdy1(ir3f54)
    b(3) =  y(ini56) * dratdumdy1(ir4f54)
    b(4) = -2.0d0 * y(ife54) * y(iprot) * ratdum(ir5f54)
    b(5) = -y(ife54) * y(iprot) * y(iprot) * dratdumdy1(ir5f54)
    b(6) =  y(ihe4) * y(ife52) * dratdumdy1(ir6f54)
    dfdy(ife52,iprot) = esum(b,6)


    ! ni56 jacobian elements
    ! d(ni56)/d(he4)
    b(1) =  y(ife52) * ratdum(irfeag)
    b(2) =  y(ife52) * y(iprot) * ratdum(ir7f54)
    dfdy(ini56,ihe4) = sum(b(1:2))

    ! d(ni56)/d(fe52)
    b(1) = y(ihe4) * ratdum(irfeag)
    b(2) = y(ihe4)* y(iprot) * ratdum(ir7f54)
    dfdy(ini56,ife52) = sum(b(1:2))

    ! d(ni56)/d(fe54)
    dfdy(ini56,ife54) = y(iprot) * y(iprot) * ratdum(ir3f54)

    ! d(ni56)/d(ni56)
    b(1) = -ratdum(irniga)
    b(2) = -ratdum(ir4f54)
    b(3) = -y(iprot) * ratdum(ir8f54)
    b(4) = -ratdum(irn56ec)
    dfdy(ini56,ini56) = esum(b,4)

    ! d(ni56)/d(prot)
    b(1) =  2.0d0 * y(ife54) * y(iprot) * ratdum(ir3f54)
    b(2) =  y(ife54) * y(iprot) * y(iprot) * dratdumdy1(ir3f54)
    b(3) = -y(ini56) * dratdumdy1(ir4f54)
    b(4) =  y(ife52) * y(ihe4)* ratdum(ir7f54)
    b(5) =  y(ife52) * y(ihe4)* y(iprot) * dratdumdy1(ir7f54)
    b(6) = -y(ini56) * ratdum(ir8f54)
    b(7) = -y(ini56) * y(iprot) * dratdumdy1(ir8f54)
    dfdy(ini56,iprot) = esum(b,7)


    ! photodisintegration neutron jacobian elements
    ! d(neut)/d(he4)
    dfdy(ineut,ihe4) = 2.0d0 * ratdum(iralf1)

    ! d(neut)/d(fe52)
    dfdy(ineut,ife52) = -2.0d0 * y(ineut) * y(ineut) * ratdum(ir2f54)

    ! d(neut)/d(fe54)
    dfdy(ineut,ife54) = 2.0d0 * ratdum(ir1f54)

    ! d(neut)/d(neut)
    b(1)  =  2.0d0 * y(ife54) * dratdumdy1(ir1f54)
    b(2)  = -4.0d0 * y(ife52) * y(ineut) * ratdum(ir2f54)
    b(3)  = -2.0d0 * y(ife52) * y(ineut) * y(ineut) * dratdumdy1(ir2f54)
    b(4)  =  2.0d0 * y(ihe4) * dratdumdy1(iralf1)
    b(5)  = -4.0d0 * y(ineut) * y(iprot)*y(iprot) * ratdum(iralf2)
    b(6)  = -2.0d0 * y(ineut)*y(ineut) * y(iprot)*y(iprot) * dratdumdy1(iralf2)
    b(7)  = -ratdum(irnep)
    dfdy(ineut,ineut) = esum(b,7)

    ! d(neut)/d(prot)
    b(1) =  2.0d0 * y(ihe4) * dratdumdy2(iralf1)
    b(2) = -4.0d0 * y(ineut)*y(ineut) * y(iprot) * ratdum(iralf2)
    b(3) = -2.0d0 * y(ineut)*y(ineut) * y(iprot)*y(iprot) * dratdumdy2(iralf2)
    b(4) =  ratdum(irpen)
    dfdy(ineut,iprot) = esum(b,4)


    ! photodisintegration proton jacobian elements
    b(1) =  2.0d0 * y(ife52) * ratdum(ir6f54)
    b(2) =  2.0d0 * ratdum(iralf1)
    dfdy(iprot,ihe4) = sum(b(1:2))

    ! d(prot)/d(fe52)
    dfdy(iprot,ife52) = 2.0d0 * y(ihe4) * ratdum(ir6f54)

    ! d(prot)/d(fe54)
    b(1) = -2.0d0 * y(iprot) * y(iprot) * ratdum(ir3f54)
    b(2) = -2.0d0 * y(iprot) * y(iprot) * ratdum(ir5f54)
    dfdy(iprot,ife54) = sum(b(1:2))

    ! d(prot)/d(ni56)
    dfdy(iprot,ini56) = 2.0d0 * ratdum(ir4f54)

    ! d(prot)/d(neut)
    b(1) =  2.0d0 * y(ihe4) * dratdumdy1(iralf1)
    b(2) = -4.0d0 * y(ineut) * y(iprot)*y(iprot) * ratdum(iralf2)
    b(3) = -2.0d0 * y(ineut)*y(ineut) * y(iprot)*y(iprot) * dratdumdy1(iralf2)
    b(4) =  ratdum(irnep)
    dfdy(iprot,ineut) = esum(b,4)

    ! d(prot)/d(prot)
    b(1)  =  -4.0d0 * y(ife54) * y(iprot) * ratdum(ir3f54)
    b(2)  =  -2.0d0 * y(ife54) * y(iprot)*y(iprot) * dratdumdy1(ir3f54)
    b(3)  =   2.0d0 * y(ini56) * dratdumdy1(ir4f54)
    b(4)  =  -4.0d0 * y(ife54) * y(iprot) * ratdum(ir5f54)
    b(5)  =  -2.0d0 * y(ife54) * y(iprot)*y(iprot) * dratdumdy1(ir5f54)
    b(6)  =   2.0d0 * y(ihe4) * y(ife52) * dratdumdy1(ir6f54)
    b(7)  =   2.0d0 * y(ihe4) * dratdumdy2(iralf1)
    b(8)  =  -4.0d0 * y(ineut)*y(ineut) * y(iprot) * ratdum(iralf2)
    b(9)  =  -2.0d0 * y(ineut)*y(ineut) * y(iprot)*y(iprot) * dratdumdy2(iralf2)
    b(10) =  -ratdum(irpen)
    dfdy(iprot,iprot) = esum(b,10)

  end subroutine dfdy_isotopes_aprox19



  ! Computes the instantaneous energy generation rate

  subroutine ener_gener_rate(dydt, enuc)

    use actual_network, only: nspec, mion, enuc_conv2

    implicit none

    double precision :: dydt(nspec), enuc

    ! This is basically e = m c**2

    enuc = sum(dydt(:) * mion(:)) * enuc_conv2

  end subroutine ener_gener_rate



  ! Compute and store the more expensive screening factors

  subroutine set_up_screening_factors()

    use screening_module, only: add_screening_factor
    use network, only: aion, zion

    implicit none

    ! note: it is critical that these are called in the exact order
    ! that the screening calls are done in the RHS routine, since we
    ! use that order in the screening

    call add_screening_factor(zion(ihe4),aion(ihe4),zion(ihe4),aion(ihe4))

    call add_screening_factor(zion(ihe4),aion(ihe4),4.0d0,8.0d0)

    call add_screening_factor(zion(ic12),aion(ic12),zion(ihe4),aion(ihe4))

    call add_screening_factor(zion(ic12),aion(ic12),zion(ic12),aion(ic12))

    call add_screening_factor(zion(ic12),aion(ic12),zion(io16),aion(io16))

    call add_screening_factor(zion(io16),aion(io16),zion(io16),aion(io16))

    call add_screening_factor(zion(io16),aion(io16),zion(ihe4),aion(ihe4))

    call add_screening_factor(zion(ine20),aion(ine20),zion(ihe4),aion(ihe4))

    call add_screening_factor(zion(img24),aion(img24),zion(ihe4),aion(ihe4))

    call add_screening_factor(13.0d0,27.0d0,1.0d0,1.0d0)

    call add_screening_factor(zion(isi28),aion(isi28),zion(ihe4),aion(ihe4))

    call add_screening_factor(15.0d0,31.0d0,1.0d0,1.0d0)

    call add_screening_factor(zion(is32),aion(is32),zion(ihe4),aion(ihe4))

    call add_screening_factor(17.0d0,35.0d0,1.0d0,1.0d0)

    call add_screening_factor(zion(iar36),aion(iar36),zion(ihe4),aion(ihe4))

    call add_screening_factor(19.0d0,39.0d0,1.0d0,1.0d0)

    call add_screening_factor(zion(ica40),aion(ica40),zion(ihe4),aion(ihe4))

    call add_screening_factor(21.0d0,43.0d0,1.0d0,1.0d0)

    call add_screening_factor(zion(iti44),aion(iti44),zion(ihe4),aion(ihe4))

    call add_screening_factor(23.0d0,47.0d0,1.0d0,1.0d0)

    call add_screening_factor(zion(icr48),aion(icr48),zion(ihe4),aion(ihe4))

    call add_screening_factor(25.0d0,51.0d0,1.0d0,1.0d0)

    call add_screening_factor(zion(ife52),aion(ife52),zion(ihe4),aion(ihe4))

    call add_screening_factor(27.0d0,55.0d0,1.0d0,1.0d0)

    call add_screening_factor(zion(ife54),aion(ife54),1.0d0,1.0d0)

    call add_screening_factor(1.0d0,2.0d0,zion(ih1),aion(ih1))

    call add_screening_factor(zion(ih1),aion(ih1),zion(ih1),aion(ih1))

    call add_screening_factor(zion(ihe3),aion(ihe3),zion(ihe3),aion(ihe3))

    call add_screening_factor(zion(ihe3),aion(ihe3),zion(ihe4),aion(ihe4))

    call add_screening_factor(zion(ic12),aion(ic12),zion(ih1),aion(ih1))

    call add_screening_factor(zion(in14),aion(in14),zion(ih1),aion(ih1))

    call add_screening_factor(zion(io16),aion(io16),zion(ih1),aion(ih1))

    call add_screening_factor(zion(in14),aion(in14),zion(ihe4),aion(ihe4))

  end subroutine set_up_screening_factors

  subroutine update_unevolved_species(state)

    !$acc routine seq

    implicit none

    type (burn_t)    :: state

  end subroutine update_unevolved_species

end module actual_rhs_module
