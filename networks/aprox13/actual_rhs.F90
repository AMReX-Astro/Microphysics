module actual_rhs_module

  use network
  use eos_type_module
  use burn_type_module
  use actual_network, only: nrates
  use rate_type_module

  implicit none

  ! Table interpolation data

  double precision, parameter :: tab_tlo = 6.0d0, tab_thi = 10.0d0
  integer, parameter :: tab_per_decade = 500
  integer, parameter :: nrattab = int(tab_thi - tab_tlo) * tab_per_decade + 1
  integer, parameter :: tab_imax = int(tab_thi - tab_tlo) * tab_per_decade + 1
  double precision, parameter :: tab_tstp = (tab_thi - tab_tlo) / dble(tab_imax - 1)

  double precision, allocatable :: rattab(:,:)
  double precision, allocatable :: drattabdt(:,:)
  double precision, allocatable :: drattabdd(:,:)
  double precision, allocatable :: ttab(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: rattab, drattabdt, drattabdd, ttab
#endif

  !$acc declare create(rattab, drattabdt, drattabdd, ttab)

contains


  subroutine actual_rhs_init()

    use screening_module, only: screening_init
    use aprox_rates_module, only: rates_init
    use extern_probin_module, only: use_tables
    use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor

    implicit none

    call rates_init()

    call set_up_screening_factors()

    call screening_init()

    if (use_tables) then

       if (parallel_IOProcessor()) then
          print *, ""
          print *, "Initializing aprox13 rate table"
          print *, ""
       endif

       call create_rates_table()

    endif

  end subroutine actual_rhs_init



  subroutine actual_rhs(state)

    use amrex_constants_module, only: ZERO
    use sneut_module, only: sneut5
    use temperature_integration_module, only: temperature_rhs

    implicit none

    ! This routine sets up the system of ode's for the aprox13
    ! nuclear reactions.  This is an alpha chain + heavy ion network
    ! with (a,p)(p,g) links.
    !
    ! Isotopes: he4,  c12,  o16,  ne20, mg24, si28, s32,
    !           ar36, ca40, ti44, cr48, fe52, ni56

    type (burn_t)    :: state
    type (rate_t)    :: rr

    logical          :: deriva

    double precision :: sneut, dsneutdt, dsneutdd, snuda, snudz
    double precision :: enuc

    double precision :: rho, temp, abar, zbar

    double precision :: y(nspec), ydot(nspec)

    !$gpu

    call evaluate_rates(state, rr)

    rho  = state % rho
    temp = state % T

    abar = state % abar
    zbar = state % zbar

    y    = state % xn * aion_inv

    deriva = .false.

    ! Call the RHS to actually get dydt.

    ydot = ZERO
    call rhs(y, rr, ydot, deriva, for_jacobian_tderiv = .false.)
    state % ydot(1:nspec) = ydot

    ! Instantaneous energy generation rate -- this needs molar fractions

    call ener_gener_rate(ydot, enuc)

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
    use eos_module
    use sneut_module, only: sneut5
    use temperature_integration_module, only: temperature_jac
    use jacobian_sparsity_module, only: set_jac_zero, set_jac_entry, get_jac_entry

    implicit none

    type (burn_t)    :: state
    type (rate_t)    :: rr

    logical          :: deriva

    double precision :: b1, sneut, dsneutdt, dsneutdd, snuda, snudz

    integer          :: j, k

    double precision :: rho, temp, abar, zbar, scratch
    double precision :: y(nspec), yderivs(nspec)

    !$gpu

    call set_jac_zero(state)

    call evaluate_rates(state, rr)

    ! Get the data from the state

    rho  = state % rho
    temp = state % T

    abar = state % abar
    zbar = state % zbar

    y    = state % xn * aion_inv

    ! Species Jacobian elements with respect to other species
    call dfdy_isotopes_aprox13(y, state, rr)

    ! Energy generation rate Jacobian elements with respect to species

    do j = 1, nspec
       do k = 1, nspec
          call get_jac_entry(state, k, j, yderivs(k))
       enddo
       call ener_gener_rate(yderivs, scratch)
       call set_jac_entry(state, net_ienuc, j, scratch)
    enddo

    ! Account for the thermal neutrino losses

    call sneut5(temp, rho, abar, zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)

    do j = 1, nspec
       b1 = ((aion(j) - abar) * abar * snuda + (zion(j) - zbar) * abar * snudz)
       call get_jac_entry(state, net_ienuc, j, scratch)
       scratch = scratch - b1
       call set_jac_entry(state, net_ienuc, j, scratch)
    enddo

    ! Evaluate the Jacobian elements with respect to temperature by
    ! calling the RHS using d(ratdum) / dT

    deriva = .true.

    call rhs(y, rr, yderivs, deriva, for_jacobian_tderiv = .true.)

    do k = 1, nspec
       call set_jac_entry(state, k, net_itemp, yderivs(k))
    enddo

    call ener_gener_rate(yderivs, scratch)
    scratch = scratch - dsneutdt
    call set_jac_entry(state, net_ienuc, net_itemp, scratch)

    ! Temperature Jacobian elements

    call temperature_jac(state)

  end subroutine actual_jac




  subroutine evaluate_rates(state, rr)

    use extern_probin_module, only: use_tables

    implicit none

    type (burn_t), intent(in)    :: state
    type (rate_t), intent(out) :: rr

    double precision :: ratraw(nrates), dratrawdt(nrates), dratrawdd(nrates)
    double precision :: ratdum(nrates), dratdumdt(nrates), dratdumdd(nrates)
    double precision :: scfac(nrates),  dscfacdt(nrates),  dscfacdd(nrates)

    double precision :: rho, temp, abar, zbar

    double precision :: y(nspec)

    !$gpu

    rho  = state % rho
    temp = state % T
    abar = state % abar
    zbar = state % zbar
    y    = state % xn * aion_inv

    ! Get the raw reaction rates
    if (use_tables) then
       call aprox13tab(temp, rho, ratraw, dratrawdt, dratrawdd)
    else
       call aprox13rat(temp, rho, ratraw, dratrawdt, dratrawdd)
    endif

    ! Do the screening here because the corrections depend on the composition
    call screen_aprox13(temp, rho, y,                 &
                        ratraw, dratrawdt, dratrawdd, &
                        ratdum, dratdumdt, dratdumdd, &
                        scfac,  dscfacdt,  dscfacdd)

    ! Save the rate data in the state.
    rr % T_eval = temp
    rr % rates(1,:) = ratdum
    rr % rates(2,:) = dratdumdt

  end subroutine evaluate_rates



  subroutine aprox13tab(btemp, bden, ratraw, dratrawdt, dratrawdd)

    implicit none

    double precision :: btemp, bden, ratraw(nrates), dratrawdt(nrates), dratrawdd(nrates)

    integer, parameter :: mp = 4

    integer          :: j, iat
    double precision :: x, x1, x2, x3, x4
    double precision :: a, b, c, d, e, f, g, h, p, q
    double precision :: alfa, beta, gama, delt

    double precision :: dtab(nrates)

    !$gpu

    ! Set the density dependence array

    dtab(ir3a)   = bden*bden
    dtab(irg3a)  = 1.0d0
    dtab(ircag)  = bden
    dtab(iroga)  = 1.0d0
    dtab(ir1212) = bden
    dtab(ir1216) = bden
    dtab(ir1616) = bden
    dtab(iroag)  = bden
    dtab(irnega) = 1.0d0
    dtab(irneag) = bden
    dtab(irmgga) = 1.0d0
    dtab(irmgag) = bden
    dtab(irsiga) = 1.0d0
    dtab(irmgap) = bden
    dtab(iralpa) = bden
    dtab(iralpg) = bden
    dtab(irsigp) = 1.0d0
    dtab(irsiag) = bden
    dtab(irsga)  = 1.0d0
    dtab(irppa)  = bden
    dtab(irsiap) = bden
    dtab(irppg)  = bden
    dtab(irsgp)  = 1.0d0
    dtab(irsag)  = bden
    dtab(irarga) = 1.0d0
    dtab(irsap)  = bden
    dtab(irclpa) = bden
    dtab(irclpg) = bden
    dtab(irargp) = 1.0d0
    dtab(irarag) = bden
    dtab(ircaga) = 1.0d0
    dtab(irarap) = bden
    dtab(irkpa)  = bden
    dtab(irkpg)  = bden
    dtab(ircagp) = 1.0d0
    dtab(ircaag) = bden
    dtab(irtiga) = 1.0d0
    dtab(ircaap) = bden
    dtab(irscpa) = bden
    dtab(irscpg) = bden
    dtab(irtigp) = 1.0d0
    dtab(irtiag) = bden
    dtab(ircrga) = 1.0d0
    dtab(irtiap) = bden
    dtab(irvpa)  = bden
    dtab(irvpg)  = bden
    dtab(ircrgp) = 1.0d0
    dtab(ircrag) = bden
    dtab(irfega) = 1.0d0
    dtab(ircrap) = bden
    dtab(irmnpa) = bden
    dtab(irmnpg) = bden
    dtab(irfegp) = 1.0d0
    dtab(irfeag) = bden
    dtab(irniga) = 1.0d0
    dtab(irfeap) = bden
    dtab(ircopa) = bden
    dtab(ircopg) = bden
    dtab(irnigp) = 1.0d0
    dtab(irr1)   = 0.0d0
    dtab(irs1)   = 0.0d0
    dtab(irt1)   = 0.0d0
    dtab(iru1)   = 0.0d0
    dtab(irv1)   = 0.0d0
    dtab(irw1)   = 0.0d0
    dtab(irx1)   = 0.0d0
    dtab(iry1)   = 0.0d0

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

       ratraw(j) = (alfa * rattab(j,iat) &
                    + beta * rattab(j,iat+1) &
                    + gama * rattab(j,iat+2) &
                    + delt * rattab(j,iat+3) ) * dtab(j)

       dratrawdt(j) = (alfa * drattabdt(j,iat) &
                       + beta * drattabdt(j,iat+1) &
                       + gama * drattabdt(j,iat+2) &
                       + delt * drattabdt(j,iat+3) ) * dtab(j)

       dratrawdd(j) = alfa * drattabdd(j,iat) &
                    + beta * drattabdd(j,iat+1) &
                    + gama * drattabdd(j,iat+2) &
                    + delt * drattabdd(j,iat+3)

    enddo

    ! hand finish the three body reactions
    dratrawdd(ir3a) = bden * dratrawdd(ir3a)

  end subroutine aprox13tab


  ! Form the table

  subroutine create_rates_table()

    implicit none

#ifdef AMREX_USE_CUDA
    integer, parameter :: numThreads=256
    integer :: numBlocks
#endif

    ! Allocate memory for the tables
    allocate(rattab(nrates, nrattab))
    allocate(drattabdt(nrates, nrattab))
    allocate(drattabdd(nrates, nrattab))
    allocate(ttab(nrattab))

#ifdef AMREX_USE_CUDA
    numBlocks = ceiling(real(tab_imax)/numThreads)
    call set_aprox13rat<<<numBlocks, numThreads>>>()
#else
    call set_aprox13rat()
#endif

    !$acc update device(rattab, drattabdt, drattabdd, ttab)

  end subroutine create_rates_table


#ifdef AMREX_USE_CUDA
  attributes(global) &
#endif
  subroutine set_aprox13rat()
#ifdef AMREX_USE_CUDA
    use cudafor
#endif
    double precision :: btemp, bden, ratraw(nrates), dratrawdt(nrates), dratrawdd(nrates)
    integer :: i, j

    bden = 1.0d0

#ifdef AMREX_USE_CUDA
    i = blockDim%x * (blockIdx%x - 1) + threadIdx%x
    if (i .le. tab_imax) then
#else
       do i = 1, tab_imax
#endif

          btemp = tab_tlo + dble(i-1) * tab_tstp
          btemp = 10.0d0**(btemp)

#ifdef AMREX_USE_CUDA
          call aprox13rat_device(btemp, bden, ratraw, dratrawdt, dratrawdd)
#else
          call aprox13rat(btemp, bden, ratraw, dratrawdt, dratrawdd)
#endif

          ttab(i) = btemp

          do j = 1, nrates

             rattab(j,i)    = ratraw(j)
             drattabdt(j,i) = dratrawdt(j)
             drattabdd(j,i) = dratrawdd(j)

          enddo

#ifdef AMREX_USE_CUDA
    endif
#else
       enddo
#endif
  end subroutine set_aprox13rat


  ! Evaluates the right hand side of the aprox13 ODEs

  subroutine rhs(y,rr,dydt,deriva,for_jacobian_tderiv)

    use amrex_constants_module, only: ZERO, SIXTH
    use microphysics_math_module, only: esum3, esum4, esum5, esum6, esum8, esum10, esum12, esum17 ! function

    implicit none

    ! deriva is used in forming the analytic Jacobian to get
    ! the derivative wrt A

    logical          :: deriva
    double precision :: y(nspec), dydt(nspec)
    type (rate_t)    :: rr

    double precision :: a(17)

    logical          :: for_jacobian_tderiv
    integer          :: index_rate

    !$gpu

    if (for_jacobian_tderiv) then
       index_rate = 2
    else
       index_rate = 1
    endif

    dydt(1:nspec) = ZERO

    ! he4 reactions
    ! heavy ion reactions
    a(1)  = 0.5d0 * y(ic12) * y(ic12) * rr % rates(index_rate, ir1212)
    a(2)  = 0.5d0 * y(ic12) * y(io16) * rr % rates(index_rate, ir1216)
    a(3)  = 0.56d0 * 0.5d0 * y(io16) * y(io16) * rr % rates(index_rate, ir1616)

    dydt(ihe4) = dydt(ihe4) + esum3(a)

    ! (a,g) and (g,a) reactions
    a(1)  = -0.5d0 * y(ihe4) * y(ihe4) * y(ihe4) * rr % rates(index_rate, ir3a)
    a(2)  =  3.0d0 * y(ic12) * rr % rates(index_rate, irg3a)
    a(3)  = -y(ihe4)  * y(ic12) * rr % rates(index_rate, ircag)
    a(4)  =  y(io16)  * rr % rates(index_rate, iroga)
    a(5)  = -y(ihe4)  * y(io16) * rr % rates(index_rate, iroag)
    a(6)  =  y(ine20) * rr % rates(index_rate, irnega)
    a(7)  = -y(ihe4)  * y(ine20) * rr % rates(index_rate, irneag)
    a(8)  =  y(img24) * rr % rates(index_rate, irmgga)
    a(9)  = -y(ihe4)  * y(img24)* rr % rates(index_rate, irmgag)
    a(10) =  y(isi28) * rr % rates(index_rate, irsiga)
    a(11) = -y(ihe4)  * y(isi28)*rr % rates(index_rate, irsiag)
    a(12) =  y(is32)  * rr % rates(index_rate, irsga)

    dydt(ihe4) = dydt(ihe4) + esum12(a)

    a(1)  = -y(ihe4)  * y(is32) * rr % rates(index_rate, irsag)
    a(2)  =  y(iar36) * rr % rates(index_rate, irarga)
    a(3)  = -y(ihe4)  * y(iar36)*rr % rates(index_rate, irarag)
    a(4)  =  y(ica40) * rr % rates(index_rate, ircaga)
    a(5)  = -y(ihe4)  * y(ica40)*rr % rates(index_rate, ircaag)
    a(6)  =  y(iti44) * rr % rates(index_rate, irtiga)
    a(7)  = -y(ihe4)  * y(iti44)*rr % rates(index_rate, irtiag)
    a(8)  =  y(icr48) * rr % rates(index_rate, ircrga)
    a(9)  = -y(ihe4)  * y(icr48)*rr % rates(index_rate, ircrag)
    a(10) =  y(ife52) * rr % rates(index_rate, irfega)
    a(11) = -y(ihe4)  * y(ife52) * rr % rates(index_rate, irfeag)
    a(12) =  y(ini56) * rr % rates(index_rate, irniga)

    dydt(ihe4) = dydt(ihe4) + esum12(a)

    ! (a,p)(p,g) and (g,p)(p,a) reactions

    if (.not.deriva) then

       a(1)  =  0.34d0*0.5d0*y(io16)*y(io16)*rr % rates(index_rate, irs1)*rr % rates(index_rate, ir1616)
       a(2)  = -y(ihe4)  * y(img24) * rr % rates(index_rate, irmgap)*(1.0d0-rr % rates(index_rate, irr1))
       a(3)  =  y(isi28) * rr % rates(index_rate, irsigp) * rr % rates(index_rate, irr1)
       a(4)  = -y(ihe4)  * y(isi28) * rr % rates(index_rate, irsiap)*(1.0d0-rr % rates(index_rate, irs1))
       a(5)  =  y(is32)  * rr % rates(index_rate, irsgp) * rr % rates(index_rate, irs1)
       a(6)  = -y(ihe4)  * y(is32) * rr % rates(index_rate, irsap)*(1.0d0-rr % rates(index_rate, irt1))
       a(7)  =  y(iar36) * rr % rates(index_rate, irargp) * rr % rates(index_rate, irt1)
       a(8)  = -y(ihe4)  * y(iar36) * rr % rates(index_rate, irarap)*(1.0d0-rr % rates(index_rate, iru1))
       a(9)  =  y(ica40) * rr % rates(index_rate, ircagp) * rr % rates(index_rate, iru1)
       a(10) = -y(ihe4)  * y(ica40) * rr % rates(index_rate, ircaap)*(1.0d0-rr % rates(index_rate, irv1))
       a(11) =  y(iti44) * rr % rates(index_rate, irtigp) * rr % rates(index_rate, irv1)
       a(12) = -y(ihe4)  * y(iti44) * rr % rates(index_rate, irtiap)*(1.0d0-rr % rates(index_rate, irw1))
       a(13) =  y(icr48) * rr % rates(index_rate, ircrgp) * rr % rates(index_rate, irw1)
       a(14) = -y(ihe4)  * y(icr48) * rr % rates(index_rate, ircrap)*(1.0d0-rr % rates(index_rate, irx1))
       a(15) =  y(ife52) * rr % rates(index_rate, irfegp) * rr % rates(index_rate, irx1)
       a(16) = -y(ihe4)  * y(ife52) * rr % rates(index_rate, irfeap)*(1.0d0-rr % rates(index_rate, iry1))
       a(17) =  y(ini56) * rr % rates(index_rate, irnigp) * rr % rates(index_rate, iry1)

       dydt(ihe4) = dydt(ihe4) + esum17(a)

    else
       a(1)  =  0.34d0*0.5d0*y(io16)*y(io16) * rr % rates(1, irs1) * rr % rates(index_rate, ir1616)
       a(2)  =  0.34d0*0.5d0*y(io16)*y(io16) * rr % rates(index_rate, irs1) * rr % rates(1, ir1616)
       a(3)  = -y(ihe4)*y(img24) * rr % rates(index_rate, irmgap)*(1.0d0 - rr % rates(1, irr1))
       a(4)  =  y(ihe4)*y(img24) * rr % rates(1, irmgap)*rr % rates(index_rate, irr1)
       a(5)  =  y(isi28) * rr % rates(1, irsigp) * rr % rates(index_rate, irr1)
       a(6)  =  y(isi28) * rr % rates(index_rate, irsigp) * rr % rates(1, irr1)
       a(7)  = -y(ihe4)*y(isi28) * rr % rates(index_rate, irsiap)*(1.0d0 - rr % rates(1, irs1))
       a(8)  =  y(ihe4)*y(isi28) * rr % rates(1, irsiap) * rr % rates(index_rate, irs1)
       a(9)  =  y(is32)  * rr % rates(1, irsgp) * rr % rates(index_rate, irs1)
       a(10) =  y(is32)  * rr % rates(index_rate, irsgp) * rr % rates(1, irs1)

       dydt(ihe4) = dydt(ihe4) + esum10(a)

       a(1)  = -y(ihe4)*y(is32) * rr % rates(index_rate, irsap)*(1.0d0 - rr % rates(1, irt1))
       a(2)  =  y(ihe4)*y(is32) * rr % rates(1, irsap)*rr % rates(index_rate, irt1)
       a(3)  =  y(iar36) * rr % rates(1, irargp) * rr % rates(index_rate, irt1)
       a(4)  =  y(iar36) * rr % rates(index_rate, irargp) * rr % rates(1, irt1)
       a(5)  = -y(ihe4)*y(iar36) * rr % rates(index_rate, irarap)*(1.0d0 - rr % rates(1, iru1))
       a(6)  =  y(ihe4)*y(iar36) * rr % rates(1, irarap)*rr % rates(index_rate, iru1)
       a(7)  =  y(ica40) * rr % rates(1, ircagp) * rr % rates(index_rate, iru1)
       a(8)  =  y(ica40) * rr % rates(index_rate, ircagp) * rr % rates(1, iru1)
       a(9)  = -y(ihe4)*y(ica40) * rr % rates(index_rate, ircaap)*(1.0d0-rr % rates(1, irv1))
       a(10) =  y(ihe4)*y(ica40) * rr % rates(1, ircaap)*rr % rates(index_rate, irv1)
       a(11) =  y(iti44) * rr % rates(1, irtigp) * rr % rates(index_rate, irv1)
       a(12) =  y(iti44) * rr % rates(index_rate, irtigp) * rr % rates(1, irv1)

       dydt(ihe4) = dydt(ihe4) + esum12(a)

       a(1)  = -y(ihe4)*y(iti44) * rr % rates(index_rate, irtiap)*(1.0d0 - rr % rates(1, irw1))
       a(2)  =  y(ihe4)*y(iti44) * rr % rates(1, irtiap)*rr % rates(index_rate, irw1)
       a(3)  =  y(icr48) * rr % rates(1, ircrgp) * rr % rates(index_rate, irw1)
       a(4)  =  y(icr48) * rr % rates(index_rate, ircrgp) * rr % rates(1, irw1)
       a(5)  = -y(ihe4)*y(icr48) * rr % rates(index_rate, ircrap)*(1.0d0 - rr % rates(1, irx1))
       a(6)  =  y(ihe4)*y(icr48) * rr % rates(1, ircrap)*rr % rates(index_rate, irx1)
       a(7)  =  y(ife52) * rr % rates(1, irfegp) * rr % rates(index_rate, irx1)
       a(8)  =  y(ife52) * rr % rates(index_rate, irfegp) * rr % rates(1, irx1)
       a(9)  = -y(ihe4)*y(ife52) * rr % rates(index_rate, irfeap)*(1.0d0 - rr % rates(1, iry1))
       a(10) =  y(ihe4)*y(ife52) * rr % rates(1, irfeap)*rr % rates(index_rate, iry1)
       a(11) =  y(ini56) * rr % rates(1, irnigp) * rr % rates(index_rate, iry1)
       a(12) =  y(ini56) * rr % rates(index_rate, irnigp) * rr % rates(1, iry1)

       dydt(ihe4) = dydt(ihe4) + esum12(a)
    end if


    ! c12 reactions
    a(1) = -y(ic12) * y(ic12) * rr % rates(index_rate, ir1212)
    a(2) = -y(ic12) * y(io16) * rr % rates(index_rate, ir1216)
    a(3) =  SIXTH * y(ihe4) * y(ihe4) * y(ihe4) * rr % rates(index_rate, ir3a)
    a(4) = -y(ic12) * rr % rates(index_rate, irg3a)
    a(5) = -y(ic12) * y(ihe4) * rr % rates(index_rate, ircag)
    a(6) =  y(io16) * rr % rates(index_rate, iroga)

    dydt(ic12) = dydt(ic12) + esum6(a)


    ! o16 reactions
    a(1) = -y(ic12) * y(io16) * rr % rates(index_rate, ir1216)
    a(2) = -y(io16) * y(io16) * rr % rates(index_rate, ir1616)
    a(3) =  y(ic12) * y(ihe4) * rr % rates(index_rate, ircag)
    a(4) = -y(io16) * y(ihe4) * rr % rates(index_rate, iroag)
    a(5) = -y(io16) * rr % rates(index_rate, iroga)
    a(6) =  y(ine20) * rr % rates(index_rate, irnega)

    dydt(io16) = dydt(io16) + esum6(a)


    ! ne20 reactions
    a(1) =  0.5d0 * y(ic12) * y(ic12) * rr % rates(index_rate, ir1212)
    a(2) =  y(io16) * y(ihe4) * rr % rates(index_rate, iroag)
    a(3) = -y(ine20) * y(ihe4) * rr % rates(index_rate, irneag)
    a(4) = -y(ine20) * rr % rates(index_rate, irnega)
    a(5) =  y(img24) * rr % rates(index_rate, irmgga)

    dydt(ine20) = dydt(ine20) + esum5(a)


    ! mg24 reactions
    a(1) =  0.5d0 * y(ic12) * y(io16) * rr % rates(index_rate, ir1216)
    a(2) =  y(ine20) * y(ihe4) * rr % rates(index_rate, irneag)
    a(3) = -y(img24) * y(ihe4) * rr % rates(index_rate, irmgag)
    a(4) = -y(img24) * rr % rates(index_rate, irmgga)
    a(5) =  y(isi28) * rr % rates(index_rate, irsiga)

    dydt(img24) = dydt(img24) + esum5(a)

    if (.not.deriva) then
       a(1) = -y(img24) * y(ihe4) * rr % rates(index_rate, irmgap)*(1.0d0-rr % rates(index_rate, irr1))
       a(2) =  y(isi28) * rr % rates(index_rate, irr1) * rr % rates(index_rate, irsigp)

       dydt(img24) = dydt(img24) + sum(a(1:2))

    else
       a(1) = -y(img24)*y(ihe4) * rr % rates(index_rate, irmgap)*(1.0d0 - rr % rates(1, irr1))
       a(2) =  y(img24)*y(ihe4) * rr % rates(1, irmgap)*rr % rates(index_rate, irr1)
       a(3) =  y(isi28) * rr % rates(1, irr1) * rr % rates(index_rate, irsigp)
       a(4) =  y(isi28) * rr % rates(index_rate, irr1) * rr % rates(1, irsigp)

       dydt(img24) = dydt(img24) + esum4(a)
    end if



    ! si28 reactions
    a(1) =  0.5d0 * y(ic12) * y(io16) * rr % rates(index_rate, ir1216)
    a(2) =  0.56d0 * 0.5d0*y(io16) * y(io16) * rr % rates(index_rate, ir1616)
    a(3) =  y(img24) * y(ihe4) * rr % rates(index_rate, irmgag)
    a(4) = -y(isi28) * y(ihe4) * rr % rates(index_rate, irsiag)
    a(5) = -y(isi28) * rr % rates(index_rate, irsiga)
    a(6) =  y(is32)  * rr % rates(index_rate, irsga)

    dydt(isi28) = dydt(isi28) + esum6(a)

    if (.not.deriva) then

       a(1) =  0.34d0*0.5d0*y(io16)*y(io16)*rr % rates(index_rate, irs1)*rr % rates(index_rate, ir1616)
       a(2) =  y(img24) * y(ihe4) * rr % rates(index_rate, irmgap)*(1.0d0-rr % rates(index_rate, irr1))
       a(3) = -y(isi28) * rr % rates(index_rate, irr1) * rr % rates(index_rate, irsigp)
       a(4) = -y(isi28) * y(ihe4) * rr % rates(index_rate, irsiap)*(1.0d0-rr % rates(index_rate, irs1))
       a(5) =  y(is32)  * rr % rates(index_rate, irs1) * rr % rates(index_rate, irsgp)

       dydt(isi28) = dydt(isi28) + esum5(a)

    else
       a(1)  =  0.34d0*0.5d0*y(io16)*y(io16) * rr % rates(1, irs1)*rr % rates(index_rate, ir1616)
       a(2)  =  0.34d0*0.5d0*y(io16)*y(io16) * rr % rates(index_rate, irs1)*rr % rates(1, ir1616)
       a(3)  =  y(img24)*y(ihe4) * rr % rates(index_rate, irmgap)*(1.0d0 - rr % rates(1, irr1))
       a(4)  = -y(img24)*y(ihe4) * rr % rates(1, irmgap)*rr % rates(index_rate, irr1)
       a(5)  = -y(isi28) * rr % rates(1, irr1) * rr % rates(index_rate, irsigp)
       a(6)  = -y(isi28) * rr % rates(index_rate, irr1) * rr % rates(1, irsigp)
       a(7)  = -y(isi28)*y(ihe4) * rr % rates(index_rate, irsiap)*(1.0d0 - rr % rates(1, irs1))
       a(8)  =  y(isi28)*y(ihe4) * rr % rates(1, irsiap)*rr % rates(index_rate, irs1)
       a(9)  = y(is32) * rr % rates(1, irs1) * rr % rates(index_rate, irsgp)
       a(10) = y(is32) * rr % rates(index_rate, irs1) * rr % rates(1, irsgp)

       dydt(isi28) = dydt(isi28) + esum10(a)
    end if



    ! s32 reactions
    a(1) =  0.1d0 * 0.5d0*y(io16) * y(io16) * rr % rates(index_rate, ir1616)
    a(2) =  y(isi28) * y(ihe4) * rr % rates(index_rate, irsiag)
    a(3) = -y(is32) * y(ihe4) * rr % rates(index_rate, irsag)
    a(4) = -y(is32) * rr % rates(index_rate, irsga)
    a(5) =  y(iar36) * rr % rates(index_rate, irarga)

    dydt(is32) = dydt(is32) + esum5(a)


    if (.not.deriva) then
       a(1) =  0.34d0*0.5d0*y(io16)*y(io16)* rr % rates(index_rate, ir1616)*(1.0d0-rr % rates(index_rate, irs1))
       a(2) =  y(isi28) * y(ihe4) * rr % rates(index_rate, irsiap)*(1.0d0-rr % rates(index_rate, irs1))
       a(3) = -y(is32) * rr % rates(index_rate, irs1) * rr % rates(index_rate, irsgp)
       a(4) = -y(is32) * y(ihe4) * rr % rates(index_rate, irsap)*(1.0d0-rr % rates(index_rate, irt1))
       a(5) =  y(iar36) * rr % rates(index_rate, irt1) * rr % rates(index_rate, irargp)

       dydt(is32) = dydt(is32) + esum5(a)

    else
       a(1)  =  0.34d0*0.5d0*y(io16)*y(io16) * rr % rates(index_rate, ir1616)*(1.0d0-rr % rates(1, irs1))
       a(2)  = -0.34d0*0.5d0*y(io16)*y(io16) * rr % rates(1, ir1616)*rr % rates(index_rate, irs1)
       a(3)  =  y(isi28)*y(ihe4) * rr % rates(index_rate, irsiap)*(1.0d0-rr % rates(1, irs1))
       a(4)  = -y(isi28)*y(ihe4) * rr % rates(1, irsiap)*rr % rates(index_rate, irs1)
       a(5)  = -y(is32) * rr % rates(1, irs1) * rr % rates(index_rate, irsgp)
       a(6)  = -y(is32) * rr % rates(index_rate, irs1) * rr % rates(1, irsgp)
       a(7)  = -y(is32)*y(ihe4) * rr % rates(index_rate, irsap)*(1.0d0-rr % rates(1, irt1))
       a(8)  =  y(is32)*y(ihe4) * rr % rates(1, irsap)*rr % rates(index_rate, irt1)
       a(9)  =  y(iar36) * rr % rates(1, irt1) * rr % rates(index_rate, irargp)
       a(10) =  y(iar36) * rr % rates(index_rate, irt1) * rr % rates(1, irargp)

       dydt(is32) = dydt(is32) + esum10(a)
    end if


    ! ar36 reactions
    a(1) =  y(is32)  * y(ihe4) * rr % rates(index_rate, irsag)
    a(2) = -y(iar36) * y(ihe4) * rr % rates(index_rate, irarag)
    a(3) = -y(iar36) * rr % rates(index_rate, irarga)
    a(4) =  y(ica40) * rr % rates(index_rate, ircaga)

    dydt(iar36) = dydt(iar36) + esum4(a)

    if (.not.deriva) then
       a(1) = y(is32)  * y(ihe4) * rr % rates(index_rate, irsap)*(1.0d0-rr % rates(index_rate, irt1))
       a(2) = -y(iar36) * rr % rates(index_rate, irt1) * rr % rates(index_rate, irargp)
       a(3) = -y(iar36) * y(ihe4) * rr % rates(index_rate, irarap)*(1.0d0-rr % rates(index_rate, iru1))
       a(4) =  y(ica40) * rr % rates(index_rate, ircagp) * rr % rates(index_rate, iru1)

       dydt(iar36) = dydt(iar36) + esum4(a)

    else
       a(1) =  y(is32)*y(ihe4) * rr % rates(index_rate, irsap)*(1.0d0 - rr % rates(1, irt1))
       a(2) = -y(is32)*y(ihe4) * rr % rates(1, irsap)*rr % rates(index_rate, irt1)
       a(3) = -y(iar36) * rr % rates(1, irt1) * rr % rates(index_rate, irargp)
       a(4) = -y(iar36) * rr % rates(index_rate, irt1) * rr % rates(1, irargp)
       a(5) = -y(iar36)*y(ihe4) * rr % rates(index_rate, irarap)*(1.0d0-rr % rates(1, iru1))
       a(6) =  y(iar36)*y(ihe4) * rr % rates(1, irarap)*rr % rates(index_rate, iru1)
       a(7) =  y(ica40) * rr % rates(1, ircagp) * rr % rates(index_rate, iru1)
       a(8) =  y(ica40) * rr % rates(index_rate, ircagp) * rr % rates(1, iru1)

       dydt(iar36) = dydt(iar36) + esum8(a)
    end if


    ! ca40 reactions
    a(1) =  y(iar36) * y(ihe4) * rr % rates(index_rate, irarag)
    a(2) = -y(ica40) * y(ihe4) * rr % rates(index_rate, ircaag)
    a(3) = -y(ica40) * rr % rates(index_rate, ircaga)
    a(4) =  y(iti44) * rr % rates(index_rate, irtiga)

    dydt(ica40) = dydt(ica40) + esum4(a)

    if (.not.deriva) then
       a(1) =  y(iar36) * y(ihe4) * rr % rates(index_rate, irarap)*(1.0d0-rr % rates(index_rate, iru1))
       a(2) = -y(ica40) * rr % rates(index_rate, ircagp) * rr % rates(index_rate, iru1)
       a(3) = -y(ica40) * y(ihe4) * rr % rates(index_rate, ircaap)*(1.0d0-rr % rates(index_rate, irv1))
       a(4) =  y(iti44) * rr % rates(index_rate, irtigp) * rr % rates(index_rate, irv1)

       dydt(ica40) = dydt(ica40) + esum4(a)

    else
       a(1) =  y(iar36)*y(ihe4) * rr % rates(index_rate, irarap)*(1.0d0-rr % rates(1, iru1))
       a(2) = -y(iar36)*y(ihe4) * rr % rates(1, irarap)*rr % rates(index_rate, iru1)
       a(3) = -y(ica40) * rr % rates(1, ircagp) * rr % rates(index_rate, iru1)
       a(4) = -y(ica40) * rr % rates(index_rate, ircagp) * rr % rates(1, iru1)
       a(5) = -y(ica40)*y(ihe4) * rr % rates(index_rate, ircaap)*(1.0d0-rr % rates(1, irv1))
       a(6) =  y(ica40)*y(ihe4) * rr % rates(1, ircaap)*rr % rates(index_rate, irv1)
       a(7) =  y(iti44) * rr % rates(1, irtigp) * rr % rates(index_rate, irv1)
       a(8) =  y(iti44) * rr % rates(index_rate, irtigp) * rr % rates(1, irv1)

       dydt(ica40) = dydt(ica40) + esum8(a)
    end if


    ! ti44 reactions
    a(1) =  y(ica40) * y(ihe4) * rr % rates(index_rate, ircaag)
    a(2) = -y(iti44) * y(ihe4) * rr % rates(index_rate, irtiag)
    a(3) = -y(iti44) * rr % rates(index_rate, irtiga)
    a(4) =  y(icr48) * rr % rates(index_rate, ircrga)

    dydt(iti44) = dydt(iti44) + esum4(a)

    if (.not.deriva) then
       a(1) =  y(ica40) * y(ihe4) * rr % rates(index_rate, ircaap)*(1.0d0-rr % rates(index_rate, irv1))
       a(2) = -y(iti44) * rr % rates(index_rate, irv1) * rr % rates(index_rate, irtigp)
       a(3) = -y(iti44) * y(ihe4) * rr % rates(index_rate, irtiap)*(1.0d0-rr % rates(index_rate, irw1))
       a(4) =  y(icr48) * rr % rates(index_rate, irw1) * rr % rates(index_rate, ircrgp)

       dydt(iti44) = dydt(iti44) + esum4(a)

    else
       a(1) =  y(ica40)*y(ihe4) * rr % rates(index_rate, ircaap)*(1.0d0-rr % rates(1, irv1))
       a(2) = -y(ica40)*y(ihe4) * rr % rates(1, ircaap)*rr % rates(index_rate, irv1)
       a(3) = -y(iti44) * rr % rates(1, irv1) * rr % rates(index_rate, irtigp)
       a(4) = -y(iti44) * rr % rates(index_rate, irv1) * rr % rates(1, irtigp)
       a(5) = -y(iti44)*y(ihe4) * rr % rates(index_rate, irtiap)*(1.0d0-rr % rates(1, irw1))
       a(6) =  y(iti44)*y(ihe4) * rr % rates(1, irtiap)*rr % rates(index_rate, irw1)
       a(7) =  y(icr48) * rr % rates(1, irw1) * rr % rates(index_rate, ircrgp)
       a(8) =  y(icr48) * rr % rates(index_rate, irw1) * rr % rates(1, ircrgp)

       dydt(iti44) = dydt(iti44) + esum8(a)
    end if


    ! cr48 reactions
    a(1) =  y(iti44) * y(ihe4) * rr % rates(index_rate, irtiag)
    a(2) = -y(icr48) * y(ihe4) * rr % rates(index_rate, ircrag)
    a(3) = -y(icr48) * rr % rates(index_rate, ircrga)
    a(4) =  y(ife52) * rr % rates(index_rate, irfega)

    dydt(icr48) = dydt(icr48) + esum4(a)

    if (.not.deriva) then
       a(1) =  y(iti44) * y(ihe4) * rr % rates(index_rate, irtiap)*(1.0d0-rr % rates(index_rate, irw1))
       a(2) = -y(icr48) * rr % rates(index_rate, irw1) * rr % rates(index_rate, ircrgp)
       a(3) = -y(icr48) * y(ihe4) * rr % rates(index_rate, ircrap)*(1.0d0-rr % rates(index_rate, irx1))
       a(4) =  y(ife52) * rr % rates(index_rate, irx1) * rr % rates(index_rate, irfegp)

       dydt(icr48) = dydt(icr48) + esum4(a)

    else
       a(1) =  y(iti44)*y(ihe4) * rr % rates(index_rate, irtiap)*(1.0d0-rr % rates(1, irw1))
       a(2) = -y(iti44)*y(ihe4) * rr % rates(1, irtiap)*rr % rates(index_rate, irw1)
       a(3) = -y(icr48) * rr % rates(1, irw1) * rr % rates(index_rate, ircrgp)
       a(4) = -y(icr48) * rr % rates(index_rate, irw1) * rr % rates(1, ircrgp)
       a(5) = -y(icr48)*y(ihe4) * rr % rates(index_rate, ircrap)*(1.0d0-rr % rates(1, irx1))
       a(6) =  y(icr48)*y(ihe4) * rr % rates(1, ircrap)*rr % rates(index_rate, irx1)
       a(7) =  y(ife52) * rr % rates(1, irx1) * rr % rates(index_rate, irfegp)
       a(8) =  y(ife52) * rr % rates(index_rate, irx1) * rr % rates(1, irfegp)

       dydt(icr48) = dydt(icr48) + esum8(a)
    end if


    ! fe52 reactions
    a(1) =  y(icr48) * y(ihe4) * rr % rates(index_rate, ircrag)
    a(2) = -y(ife52) * y(ihe4) * rr % rates(index_rate, irfeag)
    a(3) = -y(ife52) * rr % rates(index_rate, irfega)
    a(4) =  y(ini56) * rr % rates(index_rate, irniga)

    dydt(ife52) = dydt(ife52) + esum4(a)

    if (.not.deriva) then
       a(1) =  y(icr48) * y(ihe4) * rr % rates(index_rate, ircrap)*(1.0d0-rr % rates(index_rate, irx1))
       a(2) = -y(ife52) * rr % rates(index_rate, irx1) * rr % rates(index_rate, irfegp)
       a(3) = -y(ife52) * y(ihe4) * rr % rates(index_rate, irfeap)*(1.0d0-rr % rates(index_rate, iry1))
       a(4) =  y(ini56) * rr % rates(index_rate, iry1) * rr % rates(index_rate, irnigp)

       dydt(ife52) = dydt(ife52) + esum4(a)

    else
       a(1) =  y(icr48)*y(ihe4) * rr % rates(index_rate, ircrap)*(1.0d0-rr % rates(1, irx1))
       a(2) = -y(icr48)*y(ihe4) * rr % rates(1, ircrap)*rr % rates(index_rate, irx1)
       a(3) = -y(ife52) * rr % rates(1, irx1) * rr % rates(index_rate, irfegp)
       a(4) = -y(ife52) * rr % rates(index_rate, irx1) * rr % rates(1, irfegp)
       a(5) = -y(ife52)*y(ihe4) * rr % rates(index_rate, irfeap)*(1.0d0-rr % rates(1, iry1))
       a(6) =  y(ife52)*y(ihe4) * rr % rates(1, irfeap)*rr % rates(index_rate, iry1)
       a(7) =  y(ini56) * rr % rates(1, iry1) * rr % rates(index_rate, irnigp)
       a(8) =  y(ini56) * rr % rates(index_rate, iry1) * rr % rates(1, irnigp)

       dydt(ife52) = dydt(ife52) + esum8(a)
    end if


    ! ni56 reactions
    a(1) =  y(ife52) * y(ihe4) * rr % rates(index_rate, irfeag)
    a(2) = -y(ini56) * rr % rates(index_rate, irniga)

    dydt(ini56) = dydt(ini56) + sum(a(1:2))

    if (.not.deriva) then
       a(1) =  y(ife52) * y(ihe4) * rr % rates(index_rate, irfeap)*(1.0d0-rr % rates(index_rate, iry1))
       a(2) = -y(ini56) * rr % rates(index_rate, iry1) * rr % rates(index_rate, irnigp)

       dydt(ini56) = dydt(ini56) + sum(a(1:2))

    else
       a(1) =  y(ife52)*y(ihe4) * rr % rates(index_rate, irfeap)*(1.0d0-rr % rates(1, iry1))
       a(2) = -y(ife52)*y(ihe4) * rr % rates(1, irfeap)*rr % rates(index_rate, iry1)
       a(3) = -y(ini56) * rr % rates(1, iry1) * rr % rates(index_rate, irnigp)
       a(4) = -y(ini56) * rr % rates(index_rate, iry1) * rr % rates(1, irnigp)

       dydt(ini56) = dydt(ini56) + esum4(a)
    end if

  end subroutine rhs


  subroutine aprox13rat(btemp, bden, ratraw, dratrawdt, dratrawdd)

    ! this routine generates unscreened
    ! nuclear reaction rates for the aprox13 network.

    use tfactors_module
    use aprox_rates_module
    use amrex_constants_module, only: ZERO
    use extern_probin_module, only: use_c12ag_deboer17

    implicit none

    double precision :: btemp, bden
    double precision :: ratraw(nrates), dratrawdt(nrates), dratrawdd(nrates)

    integer          :: i
    double precision :: rrate,drratedt,drratedd
    type (tf_t)      :: tf

    !$gpu

    do i=1,nrates
       ratraw(i)    = ZERO
       dratrawdt(i) = ZERO
       dratrawdd(i) = ZERO
    enddo

    if (btemp .lt. 1.0d6) return


    ! get the temperature factors
    call get_tfactors(btemp, tf)

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

  end subroutine aprox13rat



  subroutine screen_aprox13(btemp, bden, y, &
                                         ratraw, dratrawdt, dratrawdd, &
                                         ratdum, dratdumdt, dratdumdd, &
                                         scfac, dscfacdt, dscfacdd)

    use amrex_constants_module, only: ZERO, ONE
    use screening_module, only: screen5, plasma_state, fill_plasma_state

    implicit none

    ! this routine computes the screening factors
    ! and applies them to the raw reaction rates,
    ! producing the final reaction rates used by the
    ! right hand sides and jacobian matrix elements

    double precision :: btemp, bden
    double precision :: y(nspec)
    double precision :: ratraw(nrates), dratrawdt(nrates), dratrawdd(nrates)
    double precision :: ratdum(nrates), dratdumdt(nrates), dratdumdd(nrates)
    double precision :: scfac(nrates),  dscfacdt(nrates),  dscfacdd(nrates)

    integer          :: i, jscr
    double precision :: sc1a,sc1adt,sc2a,sc2adt,sc3a,sc3adt
    double precision :: sc1add,sc2add
!    double precision :: sc3add

    double precision :: denom,denomdt,zz
!    double precision :: denomdd,r1dd,s1dd,t1dd,u1dd,v1dd,w1dd,x1dd,y1dd

    type (plasma_state) :: state

    !$gpu

    ! initialize
    do i=1,nrates
       ratdum(i)    = ratraw(i)
       dratdumdt(i) = dratrawdt(i)
       dratdumdd(i) = dratrawdd(i)
       scfac(i)     = ONE
       dscfacdt(i)  = ZERO
       dscfacdd(i)  = ZERO
    end do



    ! Set up the state data, which is the same for all screening factors.

    call fill_plasma_state(state, btemp, bden, y)



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
    !dscfacdt(irclpa)  = sc1add


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



    ! now form those lovely dummy proton link rates

    ! mg24(a,p)27al(p,g)28si
    ratdum(irr1)     = 0.0d0
    dratdumdt(irr1)  = 0.0d0
    !dratdumdd(irr1)  = 0.0d0
    denom    = ratdum(iralpa) + ratdum(iralpg)
    denomdt  = dratdumdt(iralpa) + dratdumdt(iralpg)
    !denomdd  = dratdumdd(iralpa) + dratdumdd(iralpg)
    if (denom .gt. 1.0d-30) then
       zz = 1.0d0/denom
       ratdum(irr1)    = ratdum(iralpa)*zz
       dratdumdt(irr1) = (dratdumdt(iralpa) - ratdum(irr1)*denomdt)*zz
       !dratdumdd(irr1) = (dratdumdd(iralpa) - ratdum(irr1)*denomdd)*zz
    end if

    ! si28(a,p)p31(p,g)s32
    ratdum(irs1)     = 0.0d0
    dratdumdt(irs1)  = 0.0d0
    !dratdumdd(irs1)  = 0.0d0
    denom    = ratdum(irppa) + ratdum(irppg)
    denomdt  = dratdumdt(irppa) + dratdumdt(irppg)
    !denomdd  = dratdumdd(irppa) + dratdumdd(irppg)
    if (denom .gt. 1.0d-30) then
       zz = 1.0d0/denom
       ratdum(irs1)    = ratdum(irppa)*zz
       dratdumdt(irs1) = (dratdumdt(irppa) - ratdum(irs1)*denomdt)*zz
       !dratdumdd(irs1) = (dratdumdd(irppa) - ratdum(irs1)*denomdd)*zz
    end if

    ! s32(a,p)cl35(p,g)ar36
    ratdum(irt1)     = 0.0d0
    dratdumdt(irt1)  = 0.0d0
    !dratdumdd(irt1)  = 0.0d0
    denom    = ratdum(irclpa) + ratdum(irclpg)
    denomdt  = dratdumdt(irclpa) + dratdumdt(irclpg)
    !denomdd  = dratdumdd(irclpa) + dratdumdd(irclpg)
    if (denom .gt. 1.0d-30) then
       zz = 1.0d0/denom
       ratdum(irt1)    = ratdum(irclpa)*zz
       dratdumdt(irt1) = (dratdumdt(irclpa) - ratdum(irt1)*denomdt)*zz
       !dratdumdd(irt1) = (dratdumdd(irclpa) - ratdum(irt1)*denomdd)*zz
    end if

    ! ar36(a,p)k39(p,g)ca40
    ratdum(iru1)     = 0.0d0
    dratdumdt(iru1)  = 0.0d0
    !dratdumdd(iru1)  = 0.0d0
    denom    = ratdum(irkpa) + ratdum(irkpg)
    denomdt  = dratdumdt(irkpa) + dratdumdt(irkpg)
    !denomdd  = dratdumdd(irkpa) + dratdumdd(irkpg)
    if (denom .gt. 1.0d-30) then
       zz   = 1.0d0/denom
       ratdum(iru1)   = ratdum(irkpa)*zz
       dratdumdt(iru1) = (dratdumdt(irkpa) - ratdum(iru1)*denomdt)*zz
       !dratdumdd(iru1) = (dratdumdd(irkpa) - ratdum(iru1)*denomdd)*zz
    end if

    ! ca40(a,p)sc43(p,g)ti44
    ratdum(irv1)     = 0.0d0
    dratdumdt(irv1)  = 0.0d0
    !dratdumdd(irv1)  = 0.0d0
    denom    = ratdum(irscpa) + ratdum(irscpg)
    denomdt  = dratdumdt(irscpa) + dratdumdt(irscpg)
    !denomdd  = dratdumdd(irscpa) + dratdumdd(irscpg)
    if (denom .gt. 1.0d-30) then
       zz  = 1.0d0/denom
       ratdum(irv1)    = ratdum(irscpa)*zz
       dratdumdt(irv1) = (dratdumdt(irscpa) - ratdum(irv1)*denomdt)*zz
       !dratdumdd(irv1) = (dratdumdd(irscpa) - ratdum(irv1)*denomdd)*zz
    end if

    ! ti44(a,p)v47(p,g)cr48
    ratdum(irw1)    = 0.0d0
    dratdumdt(irw1) = 0.0d0
    !dratdumdd(irw1) = 0.0d0
    denom    = ratdum(irvpa) + ratdum(irvpg)
    denomdt  = dratdumdt(irvpa) + dratdumdt(irvpg)
    !denomdd  = dratdumdd(irvpa) + dratdumdd(irvpg)
    if (denom .gt. 1.0d-30) then
       zz = 1.0d0/denom
       ratdum(irw1)    = ratdum(irvpa)*zz
       dratdumdt(irw1) = (dratdumdt(irvpa) - ratdum(irw1)*denomdt)*zz
       !dratdumdd(irw1) = (dratdumdd(irvpa) - ratdum(irw1)*denomdd)*zz
    end if

    ! cr48(a,p)mn51(p,g)fe52
    ratdum(irx1)    = 0.0d0
    dratdumdt(irx1) = 0.0d0
    !dratdumdd(irx1) = 0.0d0
    denom    = ratdum(irmnpa) + ratdum(irmnpg)
    denomdt  = dratdumdt(irmnpa) + dratdumdt(irmnpg)
    !denomdd  = dratdumdd(irmnpa) + dratdumdd(irmnpg)
    if (denom .gt. 1.0d-30) then
       zz = 1.0d0/denom
       ratdum(irx1)    = ratdum(irmnpa)*zz
       dratdumdt(irx1) = (dratdumdt(irmnpa) - ratdum(irx1)*denomdt)*zz
       !dratdumdd(irx1) = (dratdumdd(irmnpa) - ratdum(irx1)*denomdd)*zz
    endif

    ! fe52(a,p)co55(p,g)ni56
    ratdum(iry1)    = 0.0d0
    dratdumdt(iry1) = 0.0d0
    !dratdumdd(iry1) = 0.0d0
    denom    = ratdum(ircopa) + ratdum(ircopg)
    denomdt  = dratdumdt(ircopa) + dratdumdt(ircopg)
    !denomdd  = dratdumdd(ircopa) + dratdumdd(ircopg)
    if (denom .gt. 1.0d-30) then
       zz = 1.0d0/denom
       ratdum(iry1)    = ratdum(ircopa)*zz
       dratdumdt(iry1) = (dratdumdt(ircopa) - ratdum(iry1)*denomdt)*zz
       !dratdumdd(iry1) = (dratdumdd(ircopa) - ratdum(iry1)*denomdd)*zz
    end if

  end subroutine screen_aprox13



  subroutine dfdy_isotopes_aprox13(y,state,rr)

    use network
    use microphysics_math_module, only: esum3, esum4, esum5, esum20 ! function
    use jacobian_sparsity_module, only: set_jac_entry

    implicit none

    ! this routine sets up the dense aprox13 jacobian for the isotopes

    type (burn_t) :: state
    double precision :: y(nspec)
    type (rate_t)    :: rr

    double precision :: b(30)

    !$gpu

    ! he4 jacobian elements
    ! d(he4)/d(he4)
    b(1)  = -1.5d0 * y(ihe4) * y(ihe4) * rr % rates(1,ir3a)
    b(2)  = -y(ic12)  * rr % rates(1,ircag)
    b(3)  = -y(io16)  * rr % rates(1,iroag)
    b(4)  = -y(ine20) * rr % rates(1,irneag)
    b(5)  = -y(img24) * rr % rates(1,irmgag)
    b(6)  = -y(isi28) * rr % rates(1,irsiag)
    b(7)  = -y(is32)  * rr % rates(1,irsag)
    b(8)  = -y(iar36) * rr % rates(1,irarag)
    b(9)  = -y(ica40) * rr % rates(1,ircaag)
    b(10) = -y(iti44) * rr % rates(1,irtiag)
    b(11) = -y(icr48) * rr % rates(1,ircrag)
    b(12) = -y(ife52) * rr % rates(1,irfeag)
    b(13) = -y(img24) * rr % rates(1,irmgap) * (1.0d0-rr % rates(1,irr1))
    b(14) = -y(isi28) * rr % rates(1,irsiap) * (1.0d0-rr % rates(1,irs1))
    b(15) = -y(is32)  * rr % rates(1,irsap)  * (1.0d0-rr % rates(1,irt1))
    b(16) = -y(iar36) * rr % rates(1,irarap) * (1.0d0-rr % rates(1,iru1))
    b(17) = -y(ica40) * rr % rates(1,ircaap) * (1.0d0-rr % rates(1,irv1))
    b(18) = -y(iti44) * rr % rates(1,irtiap) * (1.0d0-rr % rates(1,irw1))
    b(19) = -y(icr48) * rr % rates(1,ircrap) * (1.0d0-rr % rates(1,irx1))
    b(20) = -y(ife52) * rr % rates(1,irfeap) * (1.0d0-rr % rates(1,iry1))
    b(30) = esum20(b)
    call set_jac_entry(state, ihe4, ihe4, b(30))


    ! d(he4)/d(c12)
    b(1) =  y(ic12) * rr % rates(1,ir1212)
    b(2) =  0.5d0 * y(io16) * rr % rates(1,ir1216)
    b(3) =  3.0d0 * rr % rates(1,irg3a)
    b(4) = -y(ihe4) * rr % rates(1,ircag)
    b(30) = esum4(b)
    call set_jac_entry(state, ihe4, ic12, b(30))

    ! d(he4)/d(o16)
    b(1) =  0.5d0 * y(ic12) * rr % rates(1,ir1216)
    b(2) =  1.12d0 * 0.5d0*y(io16) * rr % rates(1,ir1616)
    b(3) =  0.68d0 * rr % rates(1,irs1) * 0.5d0*y(io16) * rr % rates(1,ir1616)
    b(4) =  rr % rates(1,iroga)
    b(5) = -y(ihe4) * rr % rates(1,iroag)
    b(30) = esum5(b)
    call set_jac_entry(state, ihe4, io16, b(30))

    ! d(he4)/d(ne20)
    b(1) =  rr % rates(1,irnega)
    b(2) = -y(ihe4) * rr % rates(1,irneag)
    b(30) = sum(b(1:2))
    call set_jac_entry(state, ihe4, ine20, b(30))

    ! d(he4)/d(mg24)
    b(1) =  rr % rates(1,irmgga)
    b(2) = -y(ihe4) * rr % rates(1,irmgag)
    b(3) = -y(ihe4) * rr % rates(1,irmgap) * (1.0d0-rr % rates(1,irr1))
    b(30) = esum3(b)
    call set_jac_entry(state, ihe4, img24, b(30))

    ! d(he4)/d(si28)
    b(1) =  rr % rates(1,irsiga)
    b(2) = -y(ihe4) * rr % rates(1,irsiag)
    b(3) = -y(ihe4) * rr % rates(1,irsiap) * (1.0d0-rr % rates(1,irs1))
    b(4) =  rr % rates(1,irr1) * rr % rates(1,irsigp)
    b(30) = esum4(b)
    call set_jac_entry(state, ihe4, isi28, b(30))

    ! d(he4)/d(s32)
    b(1) =  rr % rates(1,irsga)
    b(2) = -y(ihe4) * rr % rates(1,irsag)
    b(3) = -y(ihe4) * rr % rates(1,irsap) * (1.0d0-rr % rates(1,irt1))
    b(4) =  rr % rates(1,irs1) * rr % rates(1,irsgp)
    b(30) = esum4(b)
    call set_jac_entry(state, ihe4, is32, b(30))

    ! d(he4)/d(ar36)
    b(1)  =  rr % rates(1,irarga)
    b(2)  = -y(ihe4) * rr % rates(1,irarag)
    b(3)  = -y(ihe4) * rr % rates(1,irarap) * (1.0d0-rr % rates(1,iru1))
    b(4)  =  rr % rates(1,irt1) * rr % rates(1,irargp)
    b(30) = esum4(b)
    call set_jac_entry(state, ihe4, iar36, b(30))

    ! d(he4)/d(ca40)
    b(1) =  rr % rates(1,ircaga)
    b(2) = -y(ihe4) * rr % rates(1,ircaag)
    b(3) = -y(ihe4) * rr % rates(1,ircaap) * (1.0d0-rr % rates(1,irv1))
    b(4) =  rr % rates(1,iru1) * rr % rates(1,ircagp)
    b(30) = esum4(b)
    call set_jac_entry(state, ihe4, ica40, b(30))

    ! d(he4)/d(ti44)
    b(1) =  rr % rates(1,irtiga)
    b(2) = -y(ihe4) * rr % rates(1,irtiag)
    b(3) = -y(ihe4) * rr % rates(1,irtiap) * (1.0d0-rr % rates(1,irw1))
    b(4) =  rr % rates(1,irv1) * rr % rates(1,irtigp)
    b(30) = esum4(b)
    call set_jac_entry(state, ihe4, iti44, b(30))

    ! d(he4)/d(cr48)
    b(1) =  rr % rates(1,ircrga)
    b(2) = -y(ihe4) * rr % rates(1,ircrag)
    b(3) = -y(ihe4) * rr % rates(1,ircrap) * (1.0d0-rr % rates(1,irx1))
    b(4) =  rr % rates(1,irw1) * rr % rates(1,ircrgp)
    b(30) = esum4(b)
    call set_jac_entry(state, ihe4, icr48, b(30))

    ! d(he4)/d(fe52)
    b(1) =  rr % rates(1,irfega)
    b(2) = -y(ihe4) * rr % rates(1,irfeag)
    b(3) = -y(ihe4) * rr % rates(1,irfeap) * (1.0d0-rr % rates(1,iry1))
    b(4) =  rr % rates(1,irx1) * rr % rates(1,irfegp)
    b(30) = esum4(b)
    call set_jac_entry(state, ihe4, ife52, b(30))

    ! d(he4)/d(ni56)
    b(1) = rr % rates(1,irniga)
    b(2) = rr % rates(1,iry1) * rr % rates(1,irnigp)
    b(30) = sum(b(1:2))
    call set_jac_entry(state, ihe4, ini56, b(30))


    ! c12 jacobian elements
    ! d(c12)/d(he4)
    b(1) =  0.5d0 * y(ihe4) * y(ihe4) * rr % rates(1,ir3a)
    b(2) = -y(ic12) * rr % rates(1,ircag)
    b(30) = sum(b(1:2))
    call set_jac_entry(state, ic12, ihe4, b(30))

    ! d(c12)/d(c12)
    b(1) = -2.0d0 * y(ic12) * rr % rates(1,ir1212)
    b(2) = -y(io16) * rr % rates(1,ir1216)
    b(3) = -rr % rates(1,irg3a)
    b(4) = -y(ihe4) * rr % rates(1,ircag)
    b(30) = esum4(b)
    call set_jac_entry(state, ic12, ic12, b(30))

    ! d(c12)/d(o16)
    b(1) = -y(ic12) * rr % rates(1,ir1216)
    b(2) =  rr % rates(1,iroga)
    b(30) = sum(b(1:2))
    call set_jac_entry(state, ic12, io16, b(30))



    ! o16 jacobian elements
    ! d(o16)/d(he4)
    b(1) =  y(ic12)*rr % rates(1,ircag)
    b(2) = -y(io16)*rr % rates(1,iroag)
    b(30) = sum(b(1:2))
    call set_jac_entry(state, io16, ihe4, b(30))

    ! d(o16)/d(c12)
    b(1) = -y(io16)*rr % rates(1,ir1216)
    b(2) =  y(ihe4)*rr % rates(1,ircag)
    b(30) = sum(b(1:2))
    call set_jac_entry(state, io16, ic12, b(30))

    ! d(o16)/d(o16)
    b(1) = -y(ic12) * rr % rates(1,ir1216)
    b(2) = -2.0d0 * y(io16) * rr % rates(1,ir1616)
    b(3) = -y(ihe4) * rr % rates(1,iroag)
    b(4) = -rr % rates(1,iroga)
    b(30) = esum4(b)
    call set_jac_entry(state, io16, io16, b(30))

    ! d(o16)/d(ne20)
    call set_jac_entry(state, io16, ine20, rr % rates(1,irnega))



    ! ne20 jacobian elements
    ! d(ne20)/d(he4)
    b(1) =  y(io16) * rr % rates(1,iroag)
    b(2) = -y(ine20) * rr % rates(1,irneag)
    b(30) = sum(b(1:2))
    call set_jac_entry(state, ine20, ihe4, b(30))

    ! d(ne20)/d(c12)
    b(30) = y(ic12) * rr % rates(1,ir1212)
    call set_jac_entry(state, ine20, ic12, b(30))

    ! d(ne20)/d(o16)
    b(30) = y(ihe4) * rr % rates(1,iroag)
    call set_jac_entry(state, ine20, io16, b(30))

    ! d(ne20)/d(ne20)
    b(1) = -y(ihe4) * rr % rates(1,irneag)
    b(2) = -rr % rates(1,irnega)
    b(30) = sum(b(1:2))
    call set_jac_entry(state, ine20, ine20, b(30))

    ! d(ne20)/d(mg24)
    call set_jac_entry(state, ine20, img24, rr % rates(1,irmgga))


    ! mg24 jacobian elements
    ! d(mg24)/d(he4)
    b(1) =  y(ine20) * rr % rates(1,irneag)
    b(2) = -y(img24) * rr % rates(1,irmgag)
    b(3) = -y(img24) * rr % rates(1,irmgap) * (1.0d0-rr % rates(1,irr1))
    b(30) = esum3(b)
    call set_jac_entry(state, img24, ihe4, b(30))

    ! d(mg24)/d(c12)
    b(30) = 0.5d0 * y(io16) * rr % rates(1,ir1216)
    call set_jac_entry(state, img24, ic12, b(30))

    ! d(mg24)/d(o16)
    b(30) =  0.5d0 * y(ic12) * rr % rates(1,ir1216)
    call set_jac_entry(state, img24, io16, b(30))

    ! d(mg24)/d(ne20)
    b(30) = y(ihe4) * rr % rates(1,irneag)
    call set_jac_entry(state, img24, ine20, b(30))

    ! d(mg24)/d(mg24)
    b(1) = -y(ihe4) * rr % rates(1,irmgag)
    b(2) = -rr % rates(1,irmgga)
    b(3) = -y(ihe4) * rr % rates(1,irmgap) * (1.0d0-rr % rates(1,irr1))
    b(30) = esum3(b)
    call set_jac_entry(state, img24, img24, b(30))

    ! d(mg24)/d(si28)
    b(1) = rr % rates(1,irsiga)
    b(2) = rr % rates(1,irr1) * rr % rates(1,irsigp)
    b(30) = sum(b(1:2))
    call set_jac_entry(state, img24, isi28, b(30))


    ! si28 jacobian elements
    ! d(si28)/d(he4)
    b(1) =  y(img24) * rr % rates(1,irmgag)
    b(2) = -y(isi28) * rr % rates(1,irsiag)
    b(3) =  y(img24) * rr % rates(1,irmgap) * (1.0d0-rr % rates(1,irr1))
    b(4) = -y(isi28) * rr % rates(1,irsiap) * (1.0d0-rr % rates(1,irs1))
    b(30) = esum4(b)
    call set_jac_entry(state, isi28, ihe4, b(30))

    ! d(si28)/d(c12)
    b(30) = 0.5d0 * y(io16) * rr % rates(1,ir1216)
    call set_jac_entry(state, isi28, ic12, b(30))

    ! d(si28)/d(o16)
    b(1) = 0.5d0 * y(ic12) * rr % rates(1,ir1216)
    b(2) = 1.12d0 * 0.5d0*y(io16) * rr % rates(1,ir1616)
    b(3) = 0.68d0 * 0.5d0*y(io16) * rr % rates(1,irs1) * rr % rates(1,ir1616)
    b(30) = esum3(b)
    call set_jac_entry(state, isi28, io16, b(30))

    ! d(si28)/d(mg24)
    b(1) =  y(ihe4) * rr % rates(1,irmgag)
    b(2) =  y(ihe4) * rr % rates(1,irmgap) * (1.0d0-rr % rates(1,irr1))
    b(30) = sum(b(1:2))
    call set_jac_entry(state, isi28, img24, b(30))

    ! d(si28)/d(si28)
    b(1) =  -y(ihe4) * rr % rates(1,irsiag)
    b(2) = -rr % rates(1,irsiga)
    b(3) = -rr % rates(1,irr1) * rr % rates(1,irsigp)
    b(4) = -y(ihe4) * rr % rates(1,irsiap) * (1.0d0-rr % rates(1,irs1))
    b(30) = esum4(b)
    call set_jac_entry(state, isi28, isi28, b(30))

    ! d(si28)/d(s32)
    b(1) = rr % rates(1,irsga)
    b(2) = rr % rates(1,irs1) * rr % rates(1,irsgp)
    b(30) = sum(b(1:2))
    call set_jac_entry(state, isi28, is32, b(30))


    ! s32 jacobian elements
    ! d(s32)/d(he4)
    b(1) =  y(isi28) * rr % rates(1,irsiag)
    b(2) = -y(is32) * rr % rates(1,irsag)
    b(3) =  y(isi28) * rr % rates(1,irsiap) * (1.0d0-rr % rates(1,irs1))
    b(4) = -y(is32) * rr % rates(1,irsap) * (1.0d0-rr % rates(1,irt1))
    b(30) = esum4(b)
    call set_jac_entry(state, is32, ihe4, b(30))

    ! d(s32)/d(o16)
    b(1) = 0.68d0*0.5d0*y(io16)*rr % rates(1,ir1616)*(1.0d0-rr % rates(1,irs1))
    b(2) = 0.2d0 * 0.5d0*y(io16) * rr % rates(1,ir1616)
    b(30) = sum(b(1:2))
    call set_jac_entry(state, is32, io16, b(30))

    ! d(s32)/d(si28)
    b(1)  =y(ihe4) * rr % rates(1,irsiag)
    b(2) = y(ihe4) * rr % rates(1,irsiap) * (1.0d0-rr % rates(1,irs1))
    b(30) = sum(b(1:2))
    call set_jac_entry(state, is32, isi28, b(30))

    ! d(s32)/d(s32)
    b(1) = -y(ihe4) * rr % rates(1,irsag)
    b(2) = -rr % rates(1,irsga)
    b(3) = -rr % rates(1,irs1) * rr % rates(1,irsgp)
    b(4) = -y(ihe4) * rr % rates(1,irsap) * (1.0d0-rr % rates(1,irt1))
    b(30) = esum4(b)
    call set_jac_entry(state, is32, is32, b(30))

    ! d(s32)/d(ar36)
    b(1) = rr % rates(1,irarga)
    b(2) = rr % rates(1,irt1) * rr % rates(1,irargp)
    b(30) = sum(b(1:2))
    call set_jac_entry(state, is32, iar36, b(30))


    ! ar36 jacobian elements
    ! d(ar36)/d(he4)
    b(1) =  y(is32)  * rr % rates(1,irsag)
    b(2) = -y(iar36) * rr % rates(1,irarag)
    b(3) =  y(is32)  * rr % rates(1,irsap) * (1.0d0-rr % rates(1,irt1))
    b(4) = -y(iar36) * rr % rates(1,irarap) * (1.0d0-rr % rates(1,iru1))
    b(30) = esum4(b)
    call set_jac_entry(state, iar36, ihe4, b(30))

    ! d(ar36)/d(s32)
    b(1) = y(ihe4) * rr % rates(1,irsag)
    b(2) = y(ihe4) * rr % rates(1,irsap) * (1.0d0-rr % rates(1,irt1))
    b(30) = sum(b(1:2))
    call set_jac_entry(state, iar36, is32, b(30))

    ! d(ar36)/d(ar36)
    b(1) = -y(ihe4) * rr % rates(1,irarag)
    b(2) = -rr % rates(1,irarga)
    b(3) = -rr % rates(1,irt1) * rr % rates(1,irargp)
    b(4) = -y(ihe4) * rr % rates(1,irarap) * (1.0d0-rr % rates(1,iru1))
    b(30) = esum4(b)
    call set_jac_entry(state, iar36, iar36, b(30))

    ! d(ar36)/d(ca40)
    b(1) = rr % rates(1,ircaga)
    b(2) = rr % rates(1,ircagp) * rr % rates(1,iru1)
    b(30) = sum(b(1:2))
    call set_jac_entry(state, iar36, ica40, b(30))


    ! ca40 jacobian elements
    ! d(ca40)/d(he4)
    b(1)  =  y(iar36) * rr % rates(1,irarag)
    b(2)  = -y(ica40) * rr % rates(1,ircaag)
    b(3)  =  y(iar36) * rr % rates(1,irarap)*(1.0d0-rr % rates(1,iru1))
    b(4)  = -y(ica40) * rr % rates(1,ircaap)*(1.0d0-rr % rates(1,irv1))
    b(30) = esum4(b)
    call set_jac_entry(state, ica40, ihe4, b(30))

    ! d(ca40)/d(ar36)
    b(1) =  y(ihe4) * rr % rates(1,irarag)
    b(2) =  y(ihe4) * rr % rates(1,irarap)*(1.0d0-rr % rates(1,iru1))
    b(30) = sum(b(1:2))
    call set_jac_entry(state, ica40, iar36, b(30))

    ! d(ca40)/d(ca40)
    b(1) = -y(ihe4) * rr % rates(1, ircaag)
    b(2) = -rr % rates(1, ircaga)
    b(3) = -rr % rates(1, ircagp) * rr % rates(1, iru1)
    b(4) = -y(ihe4) * rr % rates(1, ircaap)*(1.0d0-rr % rates(1, irv1))
    b(30) = esum4(b)
    call set_jac_entry(state, ica40, ica40, b(30))

    ! d(ca40)/d(ti44)
    b(1) = rr % rates(1, irtiga)
    b(2) = rr % rates(1, irtigp) * rr % rates(1, irv1)
    b(30) = sum(b(1:2))
    call set_jac_entry(state, ica40, iti44, b(30))



    ! ti44 jacobian elements
    ! d(ti44)/d(he4)
    b(1) =  y(ica40) * rr % rates(1, ircaag)
    b(2) = -y(iti44) * rr % rates(1, irtiag)
    b(3) =  y(ica40) * rr % rates(1, ircaap)*(1.0d0-rr % rates(1, irv1))
    b(4) = -y(iti44) * rr % rates(1, irtiap)*(1.0d0-rr % rates(1, irw1))
    b(30) = esum4(b)
    call set_jac_entry(state, iti44, ihe4, b(30))

    ! d(ti44)/d(ca40)
    b(1) =  y(ihe4) * rr % rates(1, ircaag)
    b(2) =  y(ihe4) * rr % rates(1, ircaap)*(1.0d0-rr % rates(1, irv1))
    b(30) = sum(b(1:2))
    call set_jac_entry(state, iti44, ica40, b(30))

    ! d(ti44)/d(ti44)
    b(1) = -y(ihe4) * rr % rates(1, irtiag)
    b(2) = -rr % rates(1, irtiga)
    b(3) = -rr % rates(1, irv1) * rr % rates(1, irtigp)
    b(4) = -y(ihe4) * rr % rates(1, irtiap)*(1.0d0-rr % rates(1, irw1))
    b(30) = esum4(b)
    call set_jac_entry(state, iti44, iti44, b(30))

    ! d(ti44)/d(cr48)
    b(1) = rr % rates(1, ircrga)
    b(2) = rr % rates(1, irw1) * rr % rates(1, ircrgp)
    b(30) = sum(b(1:2))
    call set_jac_entry(state, iti44, icr48, b(30))



    ! cr48 jacobian elements
    ! d(cr48)/d(he4)
    b(1) =  y(iti44) * rr % rates(1, irtiag)
    b(2) = -y(icr48) * rr % rates(1, ircrag)
    b(3) =  y(iti44) * rr % rates(1, irtiap)*(1.0d0-rr % rates(1, irw1))
    b(4) = -y(icr48) * rr % rates(1, ircrap)*(1.0d0-rr % rates(1, irx1))
    b(30) = esum4(b)
    call set_jac_entry(state, icr48, ihe4, b(30))

    ! d(cr48)/d(ti44)
    b(1) =  y(ihe4) * rr % rates(1, irtiag)
    b(2) =  y(ihe4) * rr % rates(1, irtiap)*(1.0d0-rr % rates(1, irw1))
    b(30) = sum(b(1:2))
    call set_jac_entry(state, icr48, iti44, b(30))

    ! d(cr48)/d(cr48)
    b(1) = -y(ihe4) * rr % rates(1, ircrag)
    b(2) = -rr % rates(1, ircrga)
    b(3) = -rr % rates(1, irw1) * rr % rates(1, ircrgp)
    b(4) = -y(ihe4) * rr % rates(1, ircrap)*(1.0d0-rr % rates(1, irx1))
    b(30) = esum4(b)
    call set_jac_entry(state, icr48, icr48, b(30))

    ! d(cr48)/d(fe52)
    b(1) = rr % rates(1, irfega)
    b(2) = rr % rates(1, irx1) * rr % rates(1, irfegp)
    b(30) = sum(b(1:2))
    call set_jac_entry(state, icr48, ife52, b(30))



    ! fe52 jacobian elements
    ! d(fe52)/d(he4)
    b(1) =  y(icr48) * rr % rates(1, ircrag)
    b(2) = -y(ife52) * rr % rates(1, irfeag)
    b(3) =  y(icr48) * rr % rates(1, ircrap) * (1.0d0-rr % rates(1, irx1))
    b(4) = -y(ife52) * rr % rates(1, irfeap) * (1.0d0-rr % rates(1, iry1))
    b(30) = esum4(b)
    call set_jac_entry(state, ife52, ihe4, b(30))

    ! d(fe52)/d(cr48)
    b(1) = y(ihe4) * rr % rates(1, ircrag)
    b(2) = y(ihe4) * rr % rates(1, ircrap) * (1.0d0-rr % rates(1, irx1))
    b(30) = sum(b(1:2))
    call set_jac_entry(state, ife52, icr48, b(30))

    ! d(fe52)/d(fe52)
    b(1) = -y(ihe4) * rr % rates(1, irfeag)
    b(2) = -rr % rates(1, irfega)
    b(3) = -rr % rates(1, irx1) * rr % rates(1, irfegp)
    b(4) = -y(ihe4) * rr % rates(1, irfeap) * (1.0d0-rr % rates(1, iry1))
    b(30) = esum4(b)
    call set_jac_entry(state, ife52, ife52, b(30))

    ! d(fe52)/d(ni56)
    b(1) = rr % rates(1, irniga)
    b(2) = rr % rates(1, iry1) * rr % rates(1, irnigp)
    b(30) = sum(b(1:2))
    call set_jac_entry(state, ife52, ini56, b(30))


    ! ni56 jacobian elements
    ! d(ni56)/d(he4)
    b(1) =  y(ife52) * rr % rates(1, irfeag)
    b(2) =  y(ife52) * rr % rates(1, irfeap) * (1.0d0-rr % rates(1, iry1))
    b(30) = sum(b(1:2))
    call set_jac_entry(state, ini56, ihe4, b(30))

    ! d(ni56)/d(fe52)
    b(1) = y(ihe4) * rr % rates(1, irfeag)
    b(2) = y(ihe4) * rr % rates(1, irfeap) * (1.0d0-rr % rates(1, iry1))
    b(30) = sum(b(1:2))
    call set_jac_entry(state, ini56, ife52, b(30))

    ! d(ni56)/d(ni56)
    b(1) = -rr % rates(1, irniga)
    b(2) = -rr % rates(1, iry1) * rr % rates(1, irnigp)
    b(30) = sum(b(1:2))
    call set_jac_entry(state, ini56, ini56, b(30))

  end subroutine dfdy_isotopes_aprox13



  ! Computes the instantaneous energy generation rate

  subroutine ener_gener_rate(dydt, enuc)

    use actual_network, only: nspec, mion, enuc_conv2

    implicit none

    double precision :: dydt(nspec), enuc

    !$gpu

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

  end subroutine set_up_screening_factors

  subroutine update_unevolved_species(state)

    implicit none

    type (burn_t)    :: state

    !$gpu

  end subroutine update_unevolved_species

end module actual_rhs_module
