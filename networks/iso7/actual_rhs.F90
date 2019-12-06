module actual_rhs_module

  use network
  use eos_type_module
  use burn_type_module
  use actual_network, only: nrates
  use rate_type_module
  use microphysics_type_module

  implicit none

  real(rt), parameter :: c54 = 56.0e0_rt/54.0e0_rt

  ! Table interpolation data

  real(rt), parameter :: tab_tlo = 6.0e0_rt, tab_thi = 10.0e0_rt
  integer, parameter :: tab_per_decade = 2000
  integer, parameter :: nrattab = int(tab_thi - tab_tlo) * tab_per_decade + 1
  integer, parameter :: tab_imax = int(tab_thi - tab_tlo) * tab_per_decade + 1
  real(rt), parameter :: tab_tstp = (tab_thi - tab_tlo) / dble(tab_imax - 1)

  real(rt), allocatable :: rattab(:,:)
  real(rt), allocatable :: drattabdt(:,:)
  real(rt), allocatable :: ttab(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: rattab, drattabdt, ttab
#endif

contains

  subroutine actual_rhs_init()

    use aprox_rates_module, only: rates_init
    use extern_probin_module, only: use_tables
    use screening_module, only: screening_init
    use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor

    implicit none

    call rates_init()

    call set_up_screening_factors()

    call screening_init()

    if (use_tables) then

       if (parallel_IOProcessor()) then
          print *, ""
          print *, "Initializing iso7 rate table"
          print *, ""
       endif

       call create_rates_table()

    endif

  end subroutine actual_rhs_init


  subroutine get_rates(state, rr)

    use extern_probin_module, only: use_tables

    implicit none

    type (burn_t), intent(in) :: state
    type (rate_t), intent(out) :: rr

    real(rt) :: rho, temp
    real(rt) :: y(nspec)

    real(rt) :: rate(nrates), dratedt(nrates)
    real(rt) :: dratedy1(irsi2ni:irni2si), dratedy2(irsi2ni:irni2si)

    !$gpu

    ! Get the data from the state

    rho  = state % rho
    temp = state % T
    y    = state % xn * aion_inv

    ! Get the raw reaction rates
    if (use_tables) then
       call iso7tab(temp, rho, rate, dratedt)
    else
       call iso7rat(temp, rho, rate, dratedt)
    endif

    ! Do the screening here because the corrections depend on the composition
    call screen_iso7(temp, rho, y,      &
                     rate, dratedt, &
                     dratedy1, dratedy2)

    ! Save the rate data, for the Jacobian later if we need it.

    rr % rates(1,:) = rate
    rr % rates(2,:) = dratedt
    rr % rates(3,irsi2ni:irni2si) = dratedy1
    rr % rates(4,irsi2ni:irni2si) = dratedy2

    rr % T_eval = temp

  end subroutine get_rates



  subroutine iso7tab(btemp, bden, rate, dratedt)

    implicit none

    real(rt) :: btemp, bden, rate(nrates), dratedt(nrates)

    integer, parameter :: mp = 4

    integer          :: j, iat
    real(rt) :: x, x1, x2, x3, x4
    real(rt) :: a, b, c, d, e, f, g, h, p, q
    real(rt) :: alfa, beta, gama, delt

    real(rt) :: dtab(nrates)

    !$gpu

    ! Set the density dependence array
    dtab(ircag)  = bden
    dtab(iroga)  = 1.0e0_rt
    dtab(ir3a)   = bden*bden
    dtab(irg3a)  = 1.0e0_rt
    dtab(ir1212) = bden
    dtab(ir1216) = bden
    dtab(ir1616) = bden
    dtab(iroag)  = bden
    dtab(irnega) = 1.0e0_rt
    dtab(irneag) = bden
    dtab(irmgga) = 1.0e0_rt
    dtab(irmgag) = bden
    dtab(irsiga) = 1.0e0_rt
    dtab(ircaag) = bden
    dtab(irtiga) = 1.0e0_rt
    dtab(irsi2ni) = 0.0e0_rt
    dtab(irni2si) = 0.0e0_rt

    ! hash locate
    iat = int((log10(btemp) - tab_tlo)/tab_tstp) + 1
    iat = max(1, min(iat - mp / 2 + 1, tab_imax - mp + 1))

    ! setup the lagrange interpolation coefficients for a cubic
    x  = btemp
    x1 = ttab(iat)
    x2 = ttab(iat+1)
    x3 = ttab(iat+2)
    x4 = ttab(iat+3)
    a  = x - x1
    b  = x - x2
    c  = x - x3
    d  = x - x4
    e  = x1 - x2
    f  = x1 - x3
    g  = x1 - x4
    h  = x2 - x3
    p  = x2 - x4
    q  = x3 - x4
    alfa =  b*c*d/(e*f*g)
    beta = -a*c*d/(e*h*p)
    gama =  a*b*d/(f*h*q)
    delt = -a*b*c/(g*p*q)

    ! crank off the raw reaction rates
    do j = 1, nrates

       rate(j) = (alfa * rattab(j,iat) &
                  + beta * rattab(j,iat+1) &
                  + gama * rattab(j,iat+2) &
                  + delt * rattab(j,iat+3) ) * dtab(j)

       dratedt(j) = (alfa * drattabdt(j,iat) &
                     + beta * drattabdt(j,iat+1) &
                     + gama * drattabdt(j,iat+2) &
                     + delt * drattabdt(j,iat+3) ) * dtab(j)

    enddo

  end subroutine iso7tab



  subroutine create_rates_table()

    implicit none

    ! Allocate memory for the tables
    allocate(rattab(nrates, nrattab))
    allocate(drattabdt(nrates, nrattab))
    allocate(ttab(nrattab))

    call set_iso7rat()

  end subroutine create_rates_table



  subroutine set_iso7rat()

    implicit none

    real(rt) :: btemp, bden, rate(nrates), dratedt(nrates)
    integer :: i, j

    bden = 1.0e0_rt

    do i = 1, tab_imax

       btemp = tab_tlo + dble(i-1) * tab_tstp
       btemp = 10.0e0_rt**(btemp)

       call iso7rat(btemp, bden, rate, dratedt)

       ttab(i) = btemp

       do j = 1, nrates

          rattab(j,i)    = rate(j)
          drattabdt(j,i) = dratedt(j)

       enddo

    enddo

  end subroutine set_iso7rat



  subroutine actual_rhs(state)

    use sneut_module, only: sneut5
    use temperature_integration_module, only: temperature_rhs

    implicit none

    ! This routine sets up the system of ODE's for the iso7
    ! nuclear reactions.  This is an alpha chain + heavy ion network
    ! with (a,p)(p,g) links up to silicon, as well as a Si <-> Ni link.
    !
    ! Isotopes: he4,  c12,  o16,  ne20, mg24, si28, ni56

    type (burn_t)    :: state
    type (rate_t)    :: rr

    logical          :: deriva

    real(rt) :: sneut, dsneutdt, dsneutdd, snuda, snudz
    real(rt) :: enuc

    real(rt) :: rho, temp, abar, zbar
    real(rt) :: y(nspec), ydot(nspec)

    !$gpu

    ! Get the data from the state

    rho  = state % rho
    temp = state % T
    abar = state % abar
    zbar = state % zbar
    y(:) = state % xn(:) * aion_inv(:)

    call get_rates(state, rr)

    ! Call the RHS to actually get dydt.

    deriva = .false.
    call rhs(y, rr, ydot, deriva, for_jacobian_tderiv = .false.)
    state % ydot(1:nspec) = ydot(1:nspec)

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

    use amrex_constants_module, only: ZERO
    use sneut_module, only: sneut5
    use temperature_integration_module, only: temperature_jac

    implicit none

    type (burn_t)    :: state
    type (rate_t)    :: rr

    logical          :: deriva

    real(rt) :: b1, sneut, dsneutdt, dsneutdd, snuda, snudz

    integer          :: j

    real(rt) :: rho, temp, abar, zbar
    real(rt) :: y(nspec), yderivs(nspec)

    !$gpu

    state % jac(:,:) = ZERO

    call get_rates(state, rr)

    ! Get the data from the state

    rho  = state % rho
    temp = state % T
    abar = state % abar
    zbar = state % zbar
    y    = state % xn * aion_inv

    ! Species Jacobian elements with respect to other species

    call dfdy_isotopes_iso7(y, state, rr)

    ! Energy generation rate Jacobian elements with respect to species

    do j = 1, nspec
       call ener_gener_rate(state % jac(1:nspec,j), state % jac(net_ienuc,j))
    enddo

    ! Account for the thermal neutrino losses

    call sneut5(temp, rho, abar, zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)

    do j = 1, nspec
       b1 = (-abar * abar * snuda + (zion(j) - zbar) * abar * snudz)
       state % jac(net_ienuc,j) = state % jac(net_ienuc,j) - b1
    enddo

    ! Evaluate the Jacobian elements with respect to temperature by
    ! calling the RHS using d(rate) / dT

    deriva = .true.
    call rhs(y, rr, yderivs, deriva, for_jacobian_tderiv = .true.)

    state % jac(1:nspec,net_itemp) = yderivs

    call ener_gener_rate(state % jac(1:nspec,net_itemp), state % jac(net_ienuc,net_itemp))
    state % jac(net_ienuc,net_itemp) = state % jac(net_ienuc,net_itemp) - dsneutdt

    ! Temperature Jacobian elements

    call temperature_jac(state)

  end subroutine actual_jac



  ! Evaluates the right hand side of the iso7 ODEs

  subroutine rhs(y, rr, dydt, deriva, for_jacobian_tderiv)

    use amrex_constants_module, only: ZERO, SIXTH
    use microphysics_math_module, only: esum5, esum6, esum15 ! function

    implicit none

    ! deriva is used in forming the analytic Jacobian to get
    ! the derivative wrt A

    logical          :: deriva, for_jacobian_tderiv
    real(rt) :: y(nspec), dydt(nspec)
    type (rate_t)    :: rr

    ! local variables
    integer          :: i, index_rate
    real(rt) :: a(15)

    if (for_jacobian_tderiv) then
       index_rate = 2
    else
       index_rate = 1
    endif

    !$gpu

    dydt(1:nspec) = ZERO

    ! set up the system of ode's :
    ! 4he reactions
    a(1)  =  3.0e0_rt * y(ic12) * rr % rates(index_rate,irg3a)
    a(2)  = -0.5e0_rt * y(ihe4) * y(ihe4) * y(ihe4) * rr % rates(index_rate,ir3a)
    a(3)  =  y(io16) * rr % rates(index_rate,iroga)
    a(4)  = -y(ic12) * y(ihe4) * rr % rates(index_rate,ircag)
    a(5)  =  0.5e0_rt * y(ic12) * y(ic12) * rr % rates(index_rate,ir1212)
    a(6)  =  0.5e0_rt * y(ic12) * y(io16) * rr % rates(index_rate,ir1216)
    a(7)  =  0.5e0_rt * y(io16) * y(io16) * rr % rates(index_rate,ir1616)
    a(8)  = -y(io16) * y(ihe4) * rr % rates(index_rate,iroag)
    a(9)  =  y(ine20) * rr % rates(index_rate,irnega)
    a(10) =  y(img24) * rr % rates(index_rate,irmgga)
    a(11) = -y(ine20) * y(ihe4) * rr % rates(index_rate,irneag)
    a(12) =  y(isi28) * rr % rates(index_rate,irsiga)
    a(13) = -y(img24) * y(ihe4) * rr % rates(index_rate,irmgag)
    a(14) = -7.0e0_rt * rr % rates(index_rate,irsi2ni) * y(ihe4)
    a(15) =  7.0e0_rt * rr % rates(index_rate,irni2si) * y(ini56)

    dydt(ihe4) = esum15(a)

    ! 12c reactions
    a(1) =  SIXTH * y(ihe4) * y(ihe4) * y(ihe4) * rr % rates(index_rate,ir3a)
    a(2) = -y(ic12) * rr % rates(index_rate,irg3a)
    a(3) =  y(io16) * rr % rates(index_rate,iroga)
    a(4) = -y(ic12) * y(ihe4) * rr % rates(index_rate,ircag)
    a(5) = -y(ic12) * y(ic12) * rr % rates(index_rate,ir1212)
    a(6) = -y(ic12) * y(io16) * rr % rates(index_rate,ir1216)

    dydt(ic12) = esum6(a)

    ! 16o reactions
    a(1) = -y(io16) * rr % rates(index_rate,iroga)
    a(2) =  y(ic12) * y(ihe4) * rr % rates(index_rate,ircag)
    a(3) = -y(ic12) * y(io16) * rr % rates(index_rate,ir1216)
    a(4) = -y(io16) * y(io16) * rr % rates(index_rate,ir1616)
    a(5) = -y(io16) * y(ihe4) * rr % rates(index_rate,iroag)
    a(6) =  y(ine20) * rr % rates(index_rate,irnega)

    dydt(io16) = esum6(a)

    ! 20ne reactions
    a(1) =  0.5e0_rt * y(ic12) * y(ic12) * rr % rates(index_rate,ir1212)
    a(2) =  y(io16) * y(ihe4) * rr % rates(index_rate,iroag)
    a(3) = -y(ine20) * rr % rates(index_rate,irnega)
    a(4) =  y(img24) * rr % rates(index_rate,irmgga)
    a(5) = -y(ine20) * y(ihe4) * rr % rates(index_rate,irneag)

    dydt(ine20) = esum5(a)

    ! 24mg reactions
    a(1) =  0.5e0_rt * y(ic12) * y(io16) * rr % rates(index_rate,ir1216)
    a(2) = -y(img24) * rr % rates(index_rate,irmgga)
    a(3) =  y(ine20) * y(ihe4) * rr % rates(index_rate,irneag)
    a(4) =  y(isi28) * rr % rates(index_rate,irsiga)
    a(5) = -y(img24) * y(ihe4) * rr % rates(index_rate,irmgag)

    dydt(img24) = esum5(a)

    ! 28si reactions
    a(1) =  0.5e0_rt * y(ic12) * y(io16) * rr % rates(index_rate,ir1216)
    a(2) =  0.5e0_rt * y(io16) * y(io16) * rr % rates(index_rate,ir1616)
    a(3) = -y(isi28) * rr % rates(index_rate,irsiga)
    a(4) =  y(img24) * y(ihe4) * rr % rates(index_rate,irmgag)
    a(5) = -rr % rates(index_rate,irsi2ni) * y(ihe4)
    a(6) =  rr % rates(index_rate,irni2si) * y(ini56)

    dydt(isi28) = esum6(a)

    ! ni56 reactions
    a(1) =  rr % rates(index_rate,irsi2ni) * y(ihe4)
    a(2) = -rr % rates(index_rate,irni2si) * y(ini56)

    dydt(ini56) = sum(a(1:2))

  end subroutine rhs



  subroutine iso7rat(btemp, bden, rate, dratedt)

    ! this routine generates unscreened
    ! nuclear reaction rates for the iso7 network.

    use tfactors_module
    use aprox_rates_module
    use amrex_constants_module, only: ZERO
    use extern_probin_module, only: use_c12ag_deboer17

    real(rt) :: btemp, bden
    real(rt) :: rate(nrates), dratedt(nrates)

    integer          :: i
    real(rt) :: rrate,drratedt,drratedd,drratede2_rt
    real(rt) :: ff1,dff1dt,dff1dd,ff2,dff2dt,dff2dd,tot,dtotdt,dtotdd,invtot
    type (tf_t)      :: tf

    !$gpu
    
    do i=1,nrates
       rate(i)    = ZERO
       dratedt(i) = ZERO
    enddo

    if (btemp .lt. 1.0e6_rt) return


    ! get the temperature factors
    call get_tfactors(btemp, tf)

    ! Determine which c12(a,g)o16 rate to use
    if (use_c12ag_deboer17) then
    ! deboer + 2017 c12(a,g)o16 rate
       call rate_c12ag_deboer17(tf,bden, &
                    rate(ircag),dratedt(ircag),drratedd, &
                    rate(iroga),dratedt(iroga),drratede2_rt)
    else
    ! 1.7 times cf88 c12(a,g)o16 rate
       call rate_c12ag(tf,bden, &
                    rate(ircag),dratedt(ircag),drratedd, &
                    rate(iroga),dratedt(iroga),drratede2_rt)
    endif

    ! triple alpha to c12
    call rate_tripalf(tf,bden, &
                      rate(ir3a),dratedt(ir3a),drratedd, &
                      rate(irg3a),dratedt(irg3a),drratede2_rt)

    ! c12 + c12
    call rate_c12c12(tf,bden, &
                     rate(ir1212),dratedt(ir1212),drratede2_rt, &
                     rrate,drratedt,drratedd)

    ! c12 + o16
    call rate_c12o16(tf,bden, &
                     rate(ir1216),dratedt(ir1216),drratede2_rt, &
                     rrate,drratedt,drratedd)

    ! 16o + 16o
    call rate_o16o16(tf,bden, &
                     rate(ir1616),dratedt(ir1616),drratede2_rt, &
                     rrate,drratedt,drratedd)

    ! o16(a,g)ne20
    call rate_o16ag(tf,bden, &
                    rate(iroag),dratedt(iroag),drratedd, &
                    rate(irnega),dratedt(irnega),drratede2_rt)

    ! ne20(a,g)mg24
    call rate_ne20ag(tf,bden, &
                     rate(irneag),dratedt(irneag),drratedd, &
                     rate(irmgga),dratedt(irmgga),drratede2_rt)

    ! mg24(a,g)si28
    call rate_mg24ag(tf,bden, &
                     rate(irmgag),dratedt(irmgag),drratedd, &
                     rate(irsiga),dratedt(irsiga),drratede2_rt)

    ! ca40(a,g)ti44
    call rate_ca40ag(tf,bden, &
                     rate(ircaag),dratedt(ircaag),drratedd, &
                     rate(irtiga),dratedt(irtiga),drratede2_rt)

  end subroutine iso7rat



  subroutine screen_iso7(btemp, bden, y, &
                         rate, dratedt, &
                         dratedy1, dratedy2)

    use amrex_constants_module, only: ZERO, ONE
    use screening_module, only: screen5, plasma_state, fill_plasma_state
    use tfactors_module

    ! this routine computes the screening factors
    ! and applies them to the raw reaction rates,
    ! producing the final reaction rates used by the
    ! right hand sides and jacobian matrix elements

    real(rt) :: btemp, bden
    real(rt) :: y(nspec)
    real(rt) :: rate(nrates), dratedt(nrates)
    real(rt) :: dratedy1(irsi2ni:irni2si), dratedy2(irsi2ni:irni2si)

    integer          :: i, jscr
    real(rt) :: sc1a,sc1adt,sc1add,sc2a,sc2adt,sc2add, &
                        sc3a,sc3adt,sc3add,abar,zbar,ye,z2bar, &
                        t992,t9i92,yeff_ca40,yeff_ca40dt,yeff_ti44,yeff_ti44dt, &
                        denom,denomdt,denomdd,xx,zz

    type (plasma_state) :: pstate
    type (tf_t)         :: tf

    !$gpu
    
    ! initialize
    dratedy1(:) = ZERO
    dratedy2(:) = ZERO

    ! get the temperature factors
    call get_tfactors(btemp, tf)

    ! Set up the state data, which is the same for all screening factors.

    call fill_plasma_state(pstate, btemp, bden, y(1:nspec))

    ! first the always fun triple alpha and its inverse
    jscr = 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)

    jscr = jscr + 1
    call screen5(pstate,jscr,sc2a,sc2adt,sc2add)

    sc3a   = sc1a * sc2a
    sc3adt = sc1adt*sc2a + sc1a*sc2adt

    dratedt(ir3a) = dratedt(ir3a) * sc3a + rate(ir3a) * sc3adt
    rate(ir3a)    = rate(ir3a) * sc3a

    ! c12 to o16
    jscr = jscr + 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)

    dratedt(ircag)  = dratedt(ircag) * sc1a + rate(ircag) * sc1adt
    rate(ircag)     = rate(ircag) * sc1a

    ! c12 + c12
    jscr = jscr + 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)

    dratedt(ir1212) = dratedt(ir1212) * sc1a + rate(ir1212) * sc1adt
    rate(ir1212)    = rate(ir1212) * sc1a

    ! c12 + o16
    jscr = jscr + 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)

    dratedt(ir1216) = dratedt(ir1216) * sc1a + rate(ir1216) * sc1adt
    rate(ir1216)    = rate(ir1216) * sc1a

    ! o16 + o16
    jscr = jscr + 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)

    dratedt(ir1616) = dratedt(ir1616) * sc1a + rate(ir1616) * sc1adt
    rate(ir1616)    = rate(ir1616) * sc1a

    ! o16 to ne20
    jscr = jscr + 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)

    dratedt(iroag) = dratedt(iroag) * sc1a + rate(iroag) * sc1adt
    rate(iroag)    = rate(iroag) * sc1a

    ! o16 to mg24
    jscr = jscr + 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)

    dratedt(irneag) = dratedt(irneag) * sc1a + rate(irneag) * sc1adt
    rate(irneag)    = rate(irneag) * sc1a

    ! mg24 to si28
    jscr = jscr + 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)

    dratedt(irmgag) = dratedt(irmgag) * sc1a + rate(irmgag) * sc1adt
    rate(irmgag)    = rate(irmgag) * sc1a

    ! ca40 to ti44
    jscr = jscr + 1
    call screen5(pstate,jscr,sc1a,sc1adt,sc1add)

    dratedt(ircaag) = dratedt(ircaag) * sc1a + rate(ircaag) * sc1adt
    rate(ircaag)    = rate(ircaag) * sc1a

    ! the publication, timmes, woosley & hoffman apjs, 129, 377
    ! has a typo on page 393, where its says "y(ic12)+y(io16) .gt. 0.004"
    ! it should be less than or equal to, since the idea is this piece
    ! gets activated during silicon buring, after all the c + o from
    ! oxygen burning is gone.

    if (tf%t9 .gt. 2.5e0_rt .and. y(ic12)+y(io16) .le. 4.0e-3_rt) then

       t992  = tf%t972 * tf%t9
       t9i92 = 1.0e0_rt/t992

       yeff_ca40   = t9i92 * exp(239.42e0_rt*tf%t9i - 74.741e0_rt)
       yeff_ca40dt = -yeff_ca40*(239.42e0_rt*tf%t9i2 + 4.5e0_rt*tf%t9i)

       yeff_ti44   = t992  * exp(-274.12e0_rt*tf%t9i + 74.914e0_rt)
       yeff_ti44dt = yeff_ti44*(274.12e0_rt*tf%t9i2 + 4.5e0_rt*tf%t9i)

       denom     = (bden * y(ihe4))**3

       rate(irsi2ni)     = yeff_ca40*denom*rate(ircaag)*y(isi28)
       dratedy1(irsi2ni) = 3.0e0_rt * rate(irsi2ni)/y(ihe4)
       dratedy2(irsi2ni) = yeff_ca40*denom*rate(ircaag)
       dratedt(irsi2ni)  = (yeff_ca40dt*rate(ircaag) &
            + yeff_ca40*dratedt(ircaag))*denom*y(isi28)*1.0e-9_rt

       if (denom .ne. 0.0e0_rt) then

          zz     = 1.0e0_rt/denom
          rate(irni2si) = min(1.0e10_rt,yeff_ti44*rate(irtiga)*zz)

          if (rate(irni2si) .eq. 1.0e10_rt) then
             dratedy1(irni2si) = 0.0e0_rt
             dratedt(irni2si)  = 0.0e0_rt
          else
             dratedy1(irni2si) = -3.0e0_rt * rate(irni2si)/y(ihe4)
             dratedt(irni2si)  = (yeff_ti44dt*rate(irtiga) &
                  + yeff_ti44*dratedt(irtiga))*zz*1.0e-9_rt
          end if
       endif
    end if

  end subroutine screen_iso7



  subroutine dfdy_isotopes_iso7(y, state, rr)

    use network
    use microphysics_math_module, only: esum3, esum4, esum8 ! function

    implicit none

    ! this routine sets up the dense iso7 jacobian for the isotopes

    type (burn_t) :: state
    real(rt) :: y(nspec)
    type (rate_t)    :: rr

    real(rt) :: b(8)

    !$gpu

    ! set up the jacobian
    ! 4he jacobian elements
    ! d(he4)/d(he4)
    b(1) = -1.5e0_rt * y(ihe4) * y(ihe4) * rr % rates(1,ir3a)
    b(2) = -y(ic12) * rr % rates(1,ircag)
    b(3) = -y(io16) * rr % rates(1,iroag)
    b(4) = -y(ine20) * rr % rates(1,irneag)
    b(5) = -y(img24) * rr % rates(1,irmgag)
    b(6) = -7.0e0_rt * rr % rates(1,irsi2ni)
    b(7) = -7.0e0_rt * rr % rates(3,irsi2ni) * y(ihe4)
    b(8) =  7.0e0_rt * rr % rates(3,irni2si) * y(ini56)

    state%jac(ihe4,ihe4) = esum8(b)

    ! d(he4)/d(c12)
    b(1) =  3.0e0_rt * rr % rates(1,irg3a)
    b(2) = -y(ihe4) * rr % rates(1,ircag)
    b(3) =  y(ic12) * rr % rates(1,ir1212)
    b(4) =  0.5e0_rt * y(io16) * rr % rates(1,ir1216)

    state%jac(ihe4,ic12) = esum4(b)

    ! d(he4)/d(o16)
    b(1) =  rr % rates(1,iroga)
    b(2) =  0.5e0_rt * y(ic12) * rr % rates(1,ir1216)
    b(3) =  y(io16) * rr % rates(1,ir1616)
    b(4) = -y(ihe4) * rr % rates(1,iroag)

    state%jac(ihe4,io16) = esum4(b)

    ! d(he4)/d(ne20)
    b(1) =  rr % rates(1,irnega)
    b(2) = -y(ihe4) * rr % rates(1,irneag)

    state%jac(ihe4,ine20) = sum(b(1:2))

    ! d(he4)/d(mg24)
    b(1) =  rr % rates(1,irmgga)
    b(2) = -y(ihe4) * rr % rates(1,irmgag)

    state%jac(ihe4,img24) = sum(b(1:2))

    ! d(he4)/d(si28)
    b(1) =  rr % rates(1,irsiga)
    b(2) = -7.0e0_rt * rr % rates(4,irsi2ni) * y(ihe4)

    state%jac(ihe4,isi28) = sum(b(1:2))

    ! d(he4)/d(ni56)
    b(1) =  7.0e0_rt * rr % rates(1,irni2si)

    state%jac(ihe4,ini56) = b(1)



    ! 12c jacobian elements
    ! d(c12)/d(he4)
    b(1) =  0.5e0_rt * y(ihe4) * y(ihe4) * rr % rates(1,ir3a)
    b(2) = -y(ic12) * rr % rates(1,ircag)

    state%jac(ic12,ihe4) = sum(b(1:2))

    ! d(c12)/d(c12)
    b(1) = -rr % rates(1,irg3a)
    b(2) = -y(ihe4) * rr % rates(1,ircag)
    b(3) = -2.0e0_rt * y(ic12) * rr % rates(1,ir1212)
    b(4) = -y(io16) * rr % rates(1,ir1216)

    state%jac(ic12,ic12) = esum4(b)

    ! d(c12)/d(o16)
    b(1) =  rr % rates(1,iroga)
    b(2) = -y(ic12) * rr % rates(1,ir1216)

    state%jac(ic12,io16) = sum(b(1:2))


    ! 16o jacobian elements
    ! d(o16)/d(he4)
    b(1) =  y(ic12) * rr % rates(1,ircag)
    b(2) = -y(io16) * rr % rates(1,iroag)

    state%jac(io16,ihe4) = sum(b(1:2))

    ! d(o16)/d(c12)
    b(1) =  y(ihe4) * rr % rates(1,ircag)
    b(2) = -y(io16) * rr % rates(1,ir1216)

    state%jac(io16,ic12) = sum(b(1:2))

    ! d(o16)/d(o16)
    b(1) = -rr % rates(1,iroga)
    b(2) = -y(ic12) * rr % rates(1,ir1216)
    b(3) = -2.0e0_rt * y(io16) * rr % rates(1,ir1616)
    b(4) = -y(ihe4) * rr % rates(1,iroag)

    state%jac(io16,io16) = esum4(b)

    ! d(o16)/d(ne20)
    b(1) =  rr % rates(1,irnega)

    state%jac(io16,ine20) = b(1)



    ! 20ne jacobian elements
    ! d(ne20)/d(he4)
    b(1) =  y(io16) * rr % rates(1,iroag) - y(ine20) * rr % rates(1,irneag)

    state%jac(ine20,ihe4) = b(1)

    ! d(ne20)/d(c12)
    b(1) =  y(ic12) * rr % rates(1,ir1212)

    state%jac(ine20,ic12) = b(1)

    ! d(ne20)/d(o16)
    b(1) =  y(ihe4) * rr % rates(1,iroag)

    state%jac(ine20,io16) = b(1)

    ! d(ne20)/d(ne20)
    b(1) = -rr % rates(1,irnega) - y(ihe4) * rr % rates(1,irneag)

    state%jac(ine20,ine20) = b(1)

    ! d(ne20)/d(mg24)
    b(1) =  rr % rates(1,irmgga)

    state%jac(ine20,img24) = b(1)



    ! 24mg jacobian elements
    ! d(mg24)/d(he4)
    b(1) =  y(ine20) * rr % rates(1,irneag)
    b(2) = -y(img24) * rr % rates(1,irmgag)

    state%jac(img24,ihe4) = sum(b(1:2))

    ! d(mg24)/d(c12)
    b(1) =  0.5e0_rt * y(io16) * rr % rates(1,ir1216)

    state%jac(img24,ic12) = b(1)

    ! d(mg24)/d(o16)
    b(1) =  0.5e0_rt * y(ic12) * rr % rates(1,ir1216)

    state%jac(img24,io16) = b(1)

    ! d(mg24)/d(ne20)
    b(1) =  y(ihe4) * rr % rates(1,irneag)

    state%jac(img24,ine20) = b(1)

    ! d(mg24)/d(mg24)
    b(1) = -rr % rates(1,irmgga)
    b(2) = -y(ihe4) * rr % rates(1,irmgag)

    state%jac(img24,img24) = sum(b(1:2))

    ! d(mg24)/d(si28)
    b(1) =  rr % rates(1,irsiga)

    state%jac(img24,isi28) = b(1)

    ! 28si jacobian elements
    ! d(si28)/d(he4)
    b(1) =  y(img24) * rr % rates(1,irmgag)
    b(2) = -rr % rates(1,irsi2ni)
    b(3) = -rr % rates(3,irsi2ni) * y(ihe4)
    b(4) =  rr % rates(3,irni2si) * y(ini56)

    state%jac(isi28,ihe4) = esum4(b)

    ! d(si28)/d(c12)
    b(1) =  0.5e0_rt * y(io16) * rr % rates(1,ir1216)

    state%jac(isi28,ic12) = b(1)

    ! d(si28)/d(o16)
    b(1) =  y(io16) * rr % rates(1,ir1616)
    b(2) =  0.5e0_rt * y(ic12) * rr % rates(1,ir1216)

    state%jac(isi28,io16) = sum(b(1:2))

    ! d(si28)/d(mg24)
    b(1) =  y(ihe4) * rr % rates(1,irmgag)

    state%jac(isi28,img24) = b(1)

    ! d(si28)/d(si28)
    b(1) = -rr % rates(1,irsiga)
    b(2) = -rr % rates(4,irsi2ni) * y(ihe4)

    state%jac(isi28,isi28) = sum(b(1:2))

    ! d(si28)/d(ni56)
    b(1) =  rr % rates(1,irni2si)

    state%jac(isi28,ini56) = b(1)

    ! ni56 jacobian elements
    ! d(ni56)/d(he4)
    b(1) =  rr % rates(1,irsi2ni)
    b(2) =  rr % rates(3,irsi2ni) * y(ihe4)
    b(3) = -rr % rates(3,irni2si) * y(ini56)

    state%jac(ini56,ihe4) = esum3(b)

    ! d(ni56)/d(si28)
    b(1) = rr % rates(4,irsi2ni) * y(ihe4)

    state%jac(ini56,isi28) = b(1)

    ! d(ni56)/d(ni56)
    b(1) = -rr % rates(1,irni2si)

    state%jac(ini56,ini56) = b(1)

  end subroutine dfdy_isotopes_iso7



  ! Computes the instantaneous energy generation rate

  subroutine ener_gener_rate(dydt, enuc)

    use actual_network, only: nspec, mion, enuc_conv2

    implicit none

    real(rt) :: dydt(nspec), enuc

    !$gpu

    ! This is basically e = m c**2

    enuc = sum(dydt(:) * mion(:)) * enuc_conv2

  end subroutine ener_gener_rate



  subroutine set_up_screening_factors()
    ! Compute and store the more expensive screening factors

    use screening_module, only: add_screening_factor
    use network, only: aion, zion

    implicit none

    ! note: we need to set these up in the same order that we evaluate the
    ! rates in actual_rhs.f90 (yes, it's ugly)
    call add_screening_factor(zion(ihe4),aion(ihe4),zion(ihe4),aion(ihe4))
    call add_screening_factor(zion(ihe4),aion(ihe4),4.0e0_rt,8.0e0_rt)
    call add_screening_factor(zion(ic12),aion(ic12),zion(ihe4),aion(ihe4))
    call add_screening_factor(zion(ic12),aion(ic12),zion(ic12),aion(ic12))
    call add_screening_factor(zion(ic12),aion(ic12),zion(io16),aion(io16))
    call add_screening_factor(zion(ic12),aion(ic12),zion(io16),aion(io16))
    call add_screening_factor(zion(io16),aion(io16),zion(ihe4),aion(ihe4))
    call add_screening_factor(zion(ine20),aion(ine20),zion(ihe4),aion(ihe4))
    call add_screening_factor(zion(img24),aion(img24),zion(ihe4),aion(ihe4))
    call add_screening_factor(20.0e0_rt,40.0e0_rt,zion(ihe4),aion(ihe4))

  end subroutine set_up_screening_factors


  subroutine update_unevolved_species(state)

    implicit none

    type (burn_t) :: state

    !$gpu

  end subroutine update_unevolved_species

end module actual_rhs_module
