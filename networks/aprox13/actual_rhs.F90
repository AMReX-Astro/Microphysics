module actual_rhs_module

  use network
  use eos_type_module
  use burn_type_module
  use actual_network, only: nrates
  use rate_type_module
  use microphysics_type_module, only: rt, ZERO

  implicit none

  ! Table interpolation data

  real(rt), parameter :: tab_tlo = 6.0e0_rt, tab_thi = 10.0e0_rt
  integer, parameter :: tab_per_decade = 500
  integer, parameter :: nrattab = int(tab_thi - tab_tlo) * tab_per_decade + 1
  integer, parameter :: tab_imax = int(tab_thi - tab_tlo) * tab_per_decade + 1
  real(rt), parameter :: tab_tstp = (tab_thi - tab_tlo) / dble(tab_imax - 1)

  real(rt), allocatable :: rattab(:,:)
  real(rt), allocatable :: drattabdt(:,:)
!  real(rt), allocatable :: drattabdd(:,:)
  real(rt), allocatable :: ttab(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: rattab, drattabdt, ttab !, drattabdd
#endif

  !$acc declare create(rattab, drattabdt, ttab)
  !!$acc declare create(drattabdd)

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

    real(rt) :: sneut, dsneutdt, dsneutdd, snuda, snudz
    real(rt) :: enuc

    real(rt) :: rho, temp, abar, zbar

    real(rt) :: y(nspec), ydot(nspec)

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

    use eos_module
    use sneut_module, only: sneut5
    use temperature_integration_module, only: temperature_jac
    use jacobian_sparsity_module, only: set_jac_zero, set_jac_entry, get_jac_entry

    implicit none

    type (burn_t)    :: state
    type (rate_t)    :: rr

    logical          :: deriva

    real(rt) :: b1, sneut, dsneutdt, dsneutdd, snuda, snudz

    integer          :: j, k

    real(rt) :: rho, temp, abar, zbar, scratch
    real(rt) :: y(nspec), yderivs(nspec)

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
       b1 = (-abar * abar * snuda + (zion(j) - zbar) * abar * snudz)
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

    type (burn_t), intent(in)  :: state
    type (rate_t), intent(out) :: rr

    real(rt) :: rho, temp, abar, zbar

    real(rt) :: y(nspec)

    !$gpu

    rho  = state % rho
    temp = state % T
    abar = state % abar
    zbar = state % zbar
    y    = state % xn * aion_inv

    ! Get the raw reaction rates
    if (use_tables) then
       call aprox13tab(temp, rho, rr)
    else
       call aprox13rat(temp, rho, rr)
    endif

    ! Do the screening here because the corrections depend on the composition
    call screen_aprox13(temp, rho, y, rr)

    ! Save the rate data in the state.
    rr % T_eval = temp

  end subroutine evaluate_rates



  subroutine aprox13tab(btemp, bden, rr)

    implicit none

    real(rt), intent(in   ) :: btemp, bden
    type (rate_t),    intent(inout) :: rr

    integer, parameter :: mp = 4

    integer          :: j, iat
    real(rt) :: x, x1, x2, x3, x4
    real(rt) :: a, b, c, d, e, f, g, h, p, q
    real(rt) :: alfa, beta, gama, delt

    real(rt) :: dtab(nrates)

    !$gpu

    ! Set the density dependence array

    dtab(ir3a)   = bden*bden
    dtab(irg3a)  = 1.0e0_rt
    dtab(ircag)  = bden
    dtab(iroga)  = 1.0e0_rt
    dtab(ir1212) = bden
    dtab(ir1216) = bden
    dtab(ir1616) = bden
    dtab(iroag)  = bden
    dtab(irnega) = 1.0e0_rt
    dtab(irneag) = bden
    dtab(irmgga) = 1.0e0_rt
    dtab(irmgag) = bden
    dtab(irsiga) = 1.0e0_rt
    dtab(irmgap) = bden
    dtab(iralpa) = bden
    dtab(iralpg) = bden
    dtab(irsigp) = 1.0e0_rt
    dtab(irsiag) = bden
    dtab(irsga)  = 1.0e0_rt
    dtab(irppa)  = bden
    dtab(irsiap) = bden
    dtab(irppg)  = bden
    dtab(irsgp)  = 1.0e0_rt
    dtab(irsag)  = bden
    dtab(irarga) = 1.0e0_rt
    dtab(irsap)  = bden
    dtab(irclpa) = bden
    dtab(irclpg) = bden
    dtab(irargp) = 1.0e0_rt
    dtab(irarag) = bden
    dtab(ircaga) = 1.0e0_rt
    dtab(irarap) = bden
    dtab(irkpa)  = bden
    dtab(irkpg)  = bden
    dtab(ircagp) = 1.0e0_rt
    dtab(ircaag) = bden
    dtab(irtiga) = 1.0e0_rt
    dtab(ircaap) = bden
    dtab(irscpa) = bden
    dtab(irscpg) = bden
    dtab(irtigp) = 1.0e0_rt
    dtab(irtiag) = bden
    dtab(ircrga) = 1.0e0_rt
    dtab(irtiap) = bden
    dtab(irvpa)  = bden
    dtab(irvpg)  = bden
    dtab(ircrgp) = 1.0e0_rt
    dtab(ircrag) = bden
    dtab(irfega) = 1.0e0_rt
    dtab(ircrap) = bden
    dtab(irmnpa) = bden
    dtab(irmnpg) = bden
    dtab(irfegp) = 1.0e0_rt
    dtab(irfeag) = bden
    dtab(irniga) = 1.0e0_rt
    dtab(irfeap) = bden
    dtab(ircopa) = bden
    dtab(ircopg) = bden
    dtab(irnigp) = 1.0e0_rt
    dtab(irr1)   = 0.0e0_rt
    dtab(irs1)   = 0.0e0_rt
    dtab(irt1)   = 0.0e0_rt
    dtab(iru1)   = 0.0e0_rt
    dtab(irv1)   = 0.0e0_rt
    dtab(irw1)   = 0.0e0_rt
    dtab(irx1)   = 0.0e0_rt
    dtab(iry1)   = 0.0e0_rt

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

       rr % rates(1, j) = (alfa * rattab(j,iat  ) + &
                           beta * rattab(j,iat+1) + &
                           gama * rattab(j,iat+2) + &
                           delt * rattab(j,iat+3) ) * dtab(j)

       rr % rates(2, j) = (alfa * drattabdt(j,iat  ) + &
                           beta * drattabdt(j,iat+1) + &
                           gama * drattabdt(j,iat+2) + &
                           delt * drattabdt(j,iat+3) ) * dtab(j)

       !dratrawdd(j) = alfa * drattabdd(j,iat) &
       !             + beta * drattabdd(j,iat+1) &
       !             + gama * drattabdd(j,iat+2) &
       !             + delt * drattabdd(j,iat+3)

    enddo

    ! hand finish the three body reactions
    !dratrawdd(ir3a) = bden * dratrawdd(ir3a)

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
!    allocate(drattabdd(nrates, nrattab))
    allocate(ttab(nrattab))

#ifdef AMREX_USE_CUDA
    numBlocks = ceiling(real(tab_imax)/numThreads)
    call set_aprox13rat<<<numBlocks, numThreads>>>()
#else
    call set_aprox13rat()
#endif

    !$acc update device(rattab, drattabdt, ttab)
    !!$acc update device(drattabdd)

  end subroutine create_rates_table


#ifdef AMREX_USE_CUDA
  attributes(global) &
#endif
  subroutine set_aprox13rat()
#ifdef AMREX_USE_CUDA
    use cudafor
#endif
    real(rt) :: btemp, bden
    type (rate_t)    :: rr

    integer :: i, j

    bden = 1.0e0_rt

#ifdef AMREX_USE_CUDA
    i = blockDim%x * (blockIdx%x - 1) + threadIdx%x
    if (i .le. tab_imax) then
#else
       do i = 1, tab_imax
#endif

          btemp = tab_tlo + dble(i-1) * tab_tstp
          btemp = 10.0e0_rt**(btemp)

#ifdef AMREX_USE_CUDA
#ifdef AMREX_GPU_PRAGMA_NO_HOST
          call aprox13rat(btemp, bden, rr)
#else
          call aprox13rat_device(btemp, bden, rr)
#endif
#else
          call aprox13rat(btemp, bden, rr)
#endif

          ttab(i) = btemp

          do j = 1, nrates

             rattab(j,i)    = rr % rates(1, j)
             drattabdt(j,i) = rr % rates(2, j)
             !drattabdd(j,i) = dratrawdd(j)

          enddo

#ifdef AMREX_USE_CUDA
    endif
#else
       enddo
#endif
  end subroutine set_aprox13rat


  ! Evaluates the right hand side of the aprox13 ODEs

  subroutine rhs(y,rr,dydt,deriva,for_jacobian_tderiv)

    use microphysics_math_module, only: esum3, esum4, esum5, esum6, esum8, esum10, esum12, esum17 ! function

    implicit none

    ! deriva is used in forming the analytic Jacobian to get
    ! the derivative wrt A

    logical          :: deriva
    real(rt) :: y(nspec), dydt(nspec)
    type (rate_t)    :: rr

    real(rt) :: a(17)

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
    a(1)  = 0.5e0_rt * y(ic12) * y(ic12) * rr % rates(index_rate, ir1212)
    a(2)  = 0.5e0_rt * y(ic12) * y(io16) * rr % rates(index_rate, ir1216)
    a(3)  = 0.56e0_rt * 0.5e0_rt * y(io16) * y(io16) * rr % rates(index_rate, ir1616)

    dydt(ihe4) = dydt(ihe4) + esum3(a)

    ! (a,g) and (g,a) reactions
    a(1)  = -0.5e0_rt * y(ihe4) * y(ihe4) * y(ihe4) * rr % rates(index_rate, ir3a)
    a(2)  =  3.0e0_rt * y(ic12) * rr % rates(index_rate, irg3a)
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

       a(1)  =  0.34e0_rt*0.5e0_rt*y(io16)*y(io16)*rr % rates(index_rate, irs1)*rr % rates(index_rate, ir1616)
       a(2)  = -y(ihe4)  * y(img24) * rr % rates(index_rate, irmgap)*(1.0e0_rt-rr % rates(index_rate, irr1))
       a(3)  =  y(isi28) * rr % rates(index_rate, irsigp) * rr % rates(index_rate, irr1)
       a(4)  = -y(ihe4)  * y(isi28) * rr % rates(index_rate, irsiap)*(1.0e0_rt-rr % rates(index_rate, irs1))
       a(5)  =  y(is32)  * rr % rates(index_rate, irsgp) * rr % rates(index_rate, irs1)
       a(6)  = -y(ihe4)  * y(is32) * rr % rates(index_rate, irsap)*(1.0e0_rt-rr % rates(index_rate, irt1))
       a(7)  =  y(iar36) * rr % rates(index_rate, irargp) * rr % rates(index_rate, irt1)
       a(8)  = -y(ihe4)  * y(iar36) * rr % rates(index_rate, irarap)*(1.0e0_rt-rr % rates(index_rate, iru1))
       a(9)  =  y(ica40) * rr % rates(index_rate, ircagp) * rr % rates(index_rate, iru1)
       a(10) = -y(ihe4)  * y(ica40) * rr % rates(index_rate, ircaap)*(1.0e0_rt-rr % rates(index_rate, irv1))
       a(11) =  y(iti44) * rr % rates(index_rate, irtigp) * rr % rates(index_rate, irv1)
       a(12) = -y(ihe4)  * y(iti44) * rr % rates(index_rate, irtiap)*(1.0e0_rt-rr % rates(index_rate, irw1))
       a(13) =  y(icr48) * rr % rates(index_rate, ircrgp) * rr % rates(index_rate, irw1)
       a(14) = -y(ihe4)  * y(icr48) * rr % rates(index_rate, ircrap)*(1.0e0_rt-rr % rates(index_rate, irx1))
       a(15) =  y(ife52) * rr % rates(index_rate, irfegp) * rr % rates(index_rate, irx1)
       a(16) = -y(ihe4)  * y(ife52) * rr % rates(index_rate, irfeap)*(1.0e0_rt-rr % rates(index_rate, iry1))
       a(17) =  y(ini56) * rr % rates(index_rate, irnigp) * rr % rates(index_rate, iry1)

       dydt(ihe4) = dydt(ihe4) + esum17(a)

    else
       a(1)  =  0.34e0_rt*0.5e0_rt*y(io16)*y(io16) * rr % rates(1, irs1) * rr % rates(index_rate, ir1616)
       a(2)  =  0.34e0_rt*0.5e0_rt*y(io16)*y(io16) * rr % rates(index_rate, irs1) * rr % rates(1, ir1616)
       a(3)  = -y(ihe4)*y(img24) * rr % rates(index_rate, irmgap)*(1.0e0_rt - rr % rates(1, irr1))
       a(4)  =  y(ihe4)*y(img24) * rr % rates(1, irmgap)*rr % rates(index_rate, irr1)
       a(5)  =  y(isi28) * rr % rates(1, irsigp) * rr % rates(index_rate, irr1)
       a(6)  =  y(isi28) * rr % rates(index_rate, irsigp) * rr % rates(1, irr1)
       a(7)  = -y(ihe4)*y(isi28) * rr % rates(index_rate, irsiap)*(1.0e0_rt - rr % rates(1, irs1))
       a(8)  =  y(ihe4)*y(isi28) * rr % rates(1, irsiap) * rr % rates(index_rate, irs1)
       a(9)  =  y(is32)  * rr % rates(1, irsgp) * rr % rates(index_rate, irs1)
       a(10) =  y(is32)  * rr % rates(index_rate, irsgp) * rr % rates(1, irs1)

       dydt(ihe4) = dydt(ihe4) + esum10(a)

       a(1)  = -y(ihe4)*y(is32) * rr % rates(index_rate, irsap)*(1.0e0_rt - rr % rates(1, irt1))
       a(2)  =  y(ihe4)*y(is32) * rr % rates(1, irsap)*rr % rates(index_rate, irt1)
       a(3)  =  y(iar36) * rr % rates(1, irargp) * rr % rates(index_rate, irt1)
       a(4)  =  y(iar36) * rr % rates(index_rate, irargp) * rr % rates(1, irt1)
       a(5)  = -y(ihe4)*y(iar36) * rr % rates(index_rate, irarap)*(1.0e0_rt - rr % rates(1, iru1))
       a(6)  =  y(ihe4)*y(iar36) * rr % rates(1, irarap)*rr % rates(index_rate, iru1)
       a(7)  =  y(ica40) * rr % rates(1, ircagp) * rr % rates(index_rate, iru1)
       a(8)  =  y(ica40) * rr % rates(index_rate, ircagp) * rr % rates(1, iru1)
       a(9)  = -y(ihe4)*y(ica40) * rr % rates(index_rate, ircaap)*(1.0e0_rt-rr % rates(1, irv1))
       a(10) =  y(ihe4)*y(ica40) * rr % rates(1, ircaap)*rr % rates(index_rate, irv1)
       a(11) =  y(iti44) * rr % rates(1, irtigp) * rr % rates(index_rate, irv1)
       a(12) =  y(iti44) * rr % rates(index_rate, irtigp) * rr % rates(1, irv1)

       dydt(ihe4) = dydt(ihe4) + esum12(a)

       a(1)  = -y(ihe4)*y(iti44) * rr % rates(index_rate, irtiap)*(1.0e0_rt - rr % rates(1, irw1))
       a(2)  =  y(ihe4)*y(iti44) * rr % rates(1, irtiap)*rr % rates(index_rate, irw1)
       a(3)  =  y(icr48) * rr % rates(1, ircrgp) * rr % rates(index_rate, irw1)
       a(4)  =  y(icr48) * rr % rates(index_rate, ircrgp) * rr % rates(1, irw1)
       a(5)  = -y(ihe4)*y(icr48) * rr % rates(index_rate, ircrap)*(1.0e0_rt - rr % rates(1, irx1))
       a(6)  =  y(ihe4)*y(icr48) * rr % rates(1, ircrap)*rr % rates(index_rate, irx1)
       a(7)  =  y(ife52) * rr % rates(1, irfegp) * rr % rates(index_rate, irx1)
       a(8)  =  y(ife52) * rr % rates(index_rate, irfegp) * rr % rates(1, irx1)
       a(9)  = -y(ihe4)*y(ife52) * rr % rates(index_rate, irfeap)*(1.0e0_rt - rr % rates(1, iry1))
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
    a(1) =  0.5e0_rt * y(ic12) * y(ic12) * rr % rates(index_rate, ir1212)
    a(2) =  y(io16) * y(ihe4) * rr % rates(index_rate, iroag)
    a(3) = -y(ine20) * y(ihe4) * rr % rates(index_rate, irneag)
    a(4) = -y(ine20) * rr % rates(index_rate, irnega)
    a(5) =  y(img24) * rr % rates(index_rate, irmgga)

    dydt(ine20) = dydt(ine20) + esum5(a)


    ! mg24 reactions
    a(1) =  0.5e0_rt * y(ic12) * y(io16) * rr % rates(index_rate, ir1216)
    a(2) =  y(ine20) * y(ihe4) * rr % rates(index_rate, irneag)
    a(3) = -y(img24) * y(ihe4) * rr % rates(index_rate, irmgag)
    a(4) = -y(img24) * rr % rates(index_rate, irmgga)
    a(5) =  y(isi28) * rr % rates(index_rate, irsiga)

    dydt(img24) = dydt(img24) + esum5(a)

    if (.not.deriva) then
       a(1) = -y(img24) * y(ihe4) * rr % rates(index_rate, irmgap)*(1.0e0_rt-rr % rates(index_rate, irr1))
       a(2) =  y(isi28) * rr % rates(index_rate, irr1) * rr % rates(index_rate, irsigp)

       dydt(img24) = dydt(img24) + sum(a(1:2))

    else
       a(1) = -y(img24)*y(ihe4) * rr % rates(index_rate, irmgap)*(1.0e0_rt - rr % rates(1, irr1))
       a(2) =  y(img24)*y(ihe4) * rr % rates(1, irmgap)*rr % rates(index_rate, irr1)
       a(3) =  y(isi28) * rr % rates(1, irr1) * rr % rates(index_rate, irsigp)
       a(4) =  y(isi28) * rr % rates(index_rate, irr1) * rr % rates(1, irsigp)

       dydt(img24) = dydt(img24) + esum4(a)
    end if



    ! si28 reactions
    a(1) =  0.5e0_rt * y(ic12) * y(io16) * rr % rates(index_rate, ir1216)
    a(2) =  0.56e0_rt * 0.5e0_rt*y(io16) * y(io16) * rr % rates(index_rate, ir1616)
    a(3) =  y(img24) * y(ihe4) * rr % rates(index_rate, irmgag)
    a(4) = -y(isi28) * y(ihe4) * rr % rates(index_rate, irsiag)
    a(5) = -y(isi28) * rr % rates(index_rate, irsiga)
    a(6) =  y(is32)  * rr % rates(index_rate, irsga)

    dydt(isi28) = dydt(isi28) + esum6(a)

    if (.not.deriva) then

       a(1) =  0.34e0_rt*0.5e0_rt*y(io16)*y(io16)*rr % rates(index_rate, irs1)*rr % rates(index_rate, ir1616)
       a(2) =  y(img24) * y(ihe4) * rr % rates(index_rate, irmgap)*(1.0e0_rt-rr % rates(index_rate, irr1))
       a(3) = -y(isi28) * rr % rates(index_rate, irr1) * rr % rates(index_rate, irsigp)
       a(4) = -y(isi28) * y(ihe4) * rr % rates(index_rate, irsiap)*(1.0e0_rt-rr % rates(index_rate, irs1))
       a(5) =  y(is32)  * rr % rates(index_rate, irs1) * rr % rates(index_rate, irsgp)

       dydt(isi28) = dydt(isi28) + esum5(a)

    else
       a(1)  =  0.34e0_rt*0.5e0_rt*y(io16)*y(io16) * rr % rates(1, irs1)*rr % rates(index_rate, ir1616)
       a(2)  =  0.34e0_rt*0.5e0_rt*y(io16)*y(io16) * rr % rates(index_rate, irs1)*rr % rates(1, ir1616)
       a(3)  =  y(img24)*y(ihe4) * rr % rates(index_rate, irmgap)*(1.0e0_rt - rr % rates(1, irr1))
       a(4)  = -y(img24)*y(ihe4) * rr % rates(1, irmgap)*rr % rates(index_rate, irr1)
       a(5)  = -y(isi28) * rr % rates(1, irr1) * rr % rates(index_rate, irsigp)
       a(6)  = -y(isi28) * rr % rates(index_rate, irr1) * rr % rates(1, irsigp)
       a(7)  = -y(isi28)*y(ihe4) * rr % rates(index_rate, irsiap)*(1.0e0_rt - rr % rates(1, irs1))
       a(8)  =  y(isi28)*y(ihe4) * rr % rates(1, irsiap)*rr % rates(index_rate, irs1)
       a(9)  = y(is32) * rr % rates(1, irs1) * rr % rates(index_rate, irsgp)
       a(10) = y(is32) * rr % rates(index_rate, irs1) * rr % rates(1, irsgp)

       dydt(isi28) = dydt(isi28) + esum10(a)
    end if



    ! s32 reactions
    a(1) =  0.1e0_rt * 0.5e0_rt*y(io16) * y(io16) * rr % rates(index_rate, ir1616)
    a(2) =  y(isi28) * y(ihe4) * rr % rates(index_rate, irsiag)
    a(3) = -y(is32) * y(ihe4) * rr % rates(index_rate, irsag)
    a(4) = -y(is32) * rr % rates(index_rate, irsga)
    a(5) =  y(iar36) * rr % rates(index_rate, irarga)

    dydt(is32) = dydt(is32) + esum5(a)


    if (.not.deriva) then
       a(1) =  0.34e0_rt*0.5e0_rt*y(io16)*y(io16)* rr % rates(index_rate, ir1616)*(1.0e0_rt-rr % rates(index_rate, irs1))
       a(2) =  y(isi28) * y(ihe4) * rr % rates(index_rate, irsiap)*(1.0e0_rt-rr % rates(index_rate, irs1))
       a(3) = -y(is32) * rr % rates(index_rate, irs1) * rr % rates(index_rate, irsgp)
       a(4) = -y(is32) * y(ihe4) * rr % rates(index_rate, irsap)*(1.0e0_rt-rr % rates(index_rate, irt1))
       a(5) =  y(iar36) * rr % rates(index_rate, irt1) * rr % rates(index_rate, irargp)

       dydt(is32) = dydt(is32) + esum5(a)

    else
       a(1)  =  0.34e0_rt*0.5e0_rt*y(io16)*y(io16) * rr % rates(index_rate, ir1616)*(1.0e0_rt-rr % rates(1, irs1))
       a(2)  = -0.34e0_rt*0.5e0_rt*y(io16)*y(io16) * rr % rates(1, ir1616)*rr % rates(index_rate, irs1)
       a(3)  =  y(isi28)*y(ihe4) * rr % rates(index_rate, irsiap)*(1.0e0_rt-rr % rates(1, irs1))
       a(4)  = -y(isi28)*y(ihe4) * rr % rates(1, irsiap)*rr % rates(index_rate, irs1)
       a(5)  = -y(is32) * rr % rates(1, irs1) * rr % rates(index_rate, irsgp)
       a(6)  = -y(is32) * rr % rates(index_rate, irs1) * rr % rates(1, irsgp)
       a(7)  = -y(is32)*y(ihe4) * rr % rates(index_rate, irsap)*(1.0e0_rt-rr % rates(1, irt1))
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
       a(1) = y(is32)  * y(ihe4) * rr % rates(index_rate, irsap)*(1.0e0_rt-rr % rates(index_rate, irt1))
       a(2) = -y(iar36) * rr % rates(index_rate, irt1) * rr % rates(index_rate, irargp)
       a(3) = -y(iar36) * y(ihe4) * rr % rates(index_rate, irarap)*(1.0e0_rt-rr % rates(index_rate, iru1))
       a(4) =  y(ica40) * rr % rates(index_rate, ircagp) * rr % rates(index_rate, iru1)

       dydt(iar36) = dydt(iar36) + esum4(a)

    else
       a(1) =  y(is32)*y(ihe4) * rr % rates(index_rate, irsap)*(1.0e0_rt - rr % rates(1, irt1))
       a(2) = -y(is32)*y(ihe4) * rr % rates(1, irsap)*rr % rates(index_rate, irt1)
       a(3) = -y(iar36) * rr % rates(1, irt1) * rr % rates(index_rate, irargp)
       a(4) = -y(iar36) * rr % rates(index_rate, irt1) * rr % rates(1, irargp)
       a(5) = -y(iar36)*y(ihe4) * rr % rates(index_rate, irarap)*(1.0e0_rt-rr % rates(1, iru1))
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
       a(1) =  y(iar36) * y(ihe4) * rr % rates(index_rate, irarap)*(1.0e0_rt-rr % rates(index_rate, iru1))
       a(2) = -y(ica40) * rr % rates(index_rate, ircagp) * rr % rates(index_rate, iru1)
       a(3) = -y(ica40) * y(ihe4) * rr % rates(index_rate, ircaap)*(1.0e0_rt-rr % rates(index_rate, irv1))
       a(4) =  y(iti44) * rr % rates(index_rate, irtigp) * rr % rates(index_rate, irv1)

       dydt(ica40) = dydt(ica40) + esum4(a)

    else
       a(1) =  y(iar36)*y(ihe4) * rr % rates(index_rate, irarap)*(1.0e0_rt-rr % rates(1, iru1))
       a(2) = -y(iar36)*y(ihe4) * rr % rates(1, irarap)*rr % rates(index_rate, iru1)
       a(3) = -y(ica40) * rr % rates(1, ircagp) * rr % rates(index_rate, iru1)
       a(4) = -y(ica40) * rr % rates(index_rate, ircagp) * rr % rates(1, iru1)
       a(5) = -y(ica40)*y(ihe4) * rr % rates(index_rate, ircaap)*(1.0e0_rt-rr % rates(1, irv1))
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
       a(1) =  y(ica40) * y(ihe4) * rr % rates(index_rate, ircaap)*(1.0e0_rt-rr % rates(index_rate, irv1))
       a(2) = -y(iti44) * rr % rates(index_rate, irv1) * rr % rates(index_rate, irtigp)
       a(3) = -y(iti44) * y(ihe4) * rr % rates(index_rate, irtiap)*(1.0e0_rt-rr % rates(index_rate, irw1))
       a(4) =  y(icr48) * rr % rates(index_rate, irw1) * rr % rates(index_rate, ircrgp)

       dydt(iti44) = dydt(iti44) + esum4(a)

    else
       a(1) =  y(ica40)*y(ihe4) * rr % rates(index_rate, ircaap)*(1.0e0_rt-rr % rates(1, irv1))
       a(2) = -y(ica40)*y(ihe4) * rr % rates(1, ircaap)*rr % rates(index_rate, irv1)
       a(3) = -y(iti44) * rr % rates(1, irv1) * rr % rates(index_rate, irtigp)
       a(4) = -y(iti44) * rr % rates(index_rate, irv1) * rr % rates(1, irtigp)
       a(5) = -y(iti44)*y(ihe4) * rr % rates(index_rate, irtiap)*(1.0e0_rt-rr % rates(1, irw1))
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
       a(1) =  y(iti44) * y(ihe4) * rr % rates(index_rate, irtiap)*(1.0e0_rt-rr % rates(index_rate, irw1))
       a(2) = -y(icr48) * rr % rates(index_rate, irw1) * rr % rates(index_rate, ircrgp)
       a(3) = -y(icr48) * y(ihe4) * rr % rates(index_rate, ircrap)*(1.0e0_rt-rr % rates(index_rate, irx1))
       a(4) =  y(ife52) * rr % rates(index_rate, irx1) * rr % rates(index_rate, irfegp)

       dydt(icr48) = dydt(icr48) + esum4(a)

    else
       a(1) =  y(iti44)*y(ihe4) * rr % rates(index_rate, irtiap)*(1.0e0_rt-rr % rates(1, irw1))
       a(2) = -y(iti44)*y(ihe4) * rr % rates(1, irtiap)*rr % rates(index_rate, irw1)
       a(3) = -y(icr48) * rr % rates(1, irw1) * rr % rates(index_rate, ircrgp)
       a(4) = -y(icr48) * rr % rates(index_rate, irw1) * rr % rates(1, ircrgp)
       a(5) = -y(icr48)*y(ihe4) * rr % rates(index_rate, ircrap)*(1.0e0_rt-rr % rates(1, irx1))
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
       a(1) =  y(icr48) * y(ihe4) * rr % rates(index_rate, ircrap)*(1.0e0_rt-rr % rates(index_rate, irx1))
       a(2) = -y(ife52) * rr % rates(index_rate, irx1) * rr % rates(index_rate, irfegp)
       a(3) = -y(ife52) * y(ihe4) * rr % rates(index_rate, irfeap)*(1.0e0_rt-rr % rates(index_rate, iry1))
       a(4) =  y(ini56) * rr % rates(index_rate, iry1) * rr % rates(index_rate, irnigp)

       dydt(ife52) = dydt(ife52) + esum4(a)

    else
       a(1) =  y(icr48)*y(ihe4) * rr % rates(index_rate, ircrap)*(1.0e0_rt-rr % rates(1, irx1))
       a(2) = -y(icr48)*y(ihe4) * rr % rates(1, ircrap)*rr % rates(index_rate, irx1)
       a(3) = -y(ife52) * rr % rates(1, irx1) * rr % rates(index_rate, irfegp)
       a(4) = -y(ife52) * rr % rates(index_rate, irx1) * rr % rates(1, irfegp)
       a(5) = -y(ife52)*y(ihe4) * rr % rates(index_rate, irfeap)*(1.0e0_rt-rr % rates(1, iry1))
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
       a(1) =  y(ife52) * y(ihe4) * rr % rates(index_rate, irfeap)*(1.0e0_rt-rr % rates(index_rate, iry1))
       a(2) = -y(ini56) * rr % rates(index_rate, iry1) * rr % rates(index_rate, irnigp)

       dydt(ini56) = dydt(ini56) + sum(a(1:2))

    else
       a(1) =  y(ife52)*y(ihe4) * rr % rates(index_rate, irfeap)*(1.0e0_rt-rr % rates(1, iry1))
       a(2) = -y(ife52)*y(ihe4) * rr % rates(1, irfeap)*rr % rates(index_rate, iry1)
       a(3) = -y(ini56) * rr % rates(1, iry1) * rr % rates(index_rate, irnigp)
       a(4) = -y(ini56) * rr % rates(index_rate, iry1) * rr % rates(1, irnigp)

       dydt(ini56) = dydt(ini56) + esum4(a)
    end if

  end subroutine rhs


  subroutine aprox13rat(btemp, bden, rr)

    ! this routine generates unscreened
    ! nuclear reaction rates for the aprox13 network.

    use tfactors_module
    use aprox_rates_module
    use extern_probin_module, only: use_c12ag_deboer17

    implicit none

    real(rt), intent(in   ) :: btemp, bden
    type (rate_t),    intent(inout) :: rr

    integer          :: i
    real(rt) :: rrate,drratedt
    real(rt) :: dratedd1, dratedd2
    type (tf_t)      :: tf

    !$gpu

    do i=1,nrates
       rr % rates(1,i) = ZERO
       rr % rates(2,i) = ZERO
       !dratrawdd(i) = ZERO
    enddo

    if (btemp .lt. 1.0e6_rt) return


    ! get the temperature factors
    call get_tfactors(btemp, tf)

    ! Determine which c12(a,g)o16 rate to use
    if (use_c12ag_deboer17) then
    ! deboer + 2017 c12(a,g)o16 rate
       call rate_c12ag_deboer17(tf,bden, &
                    rr % rates(1,ircag),rr % rates(2,ircag),dratedd1, &
                    rr % rates(1,iroga),rr % rates(2,iroga),dratedd2)
    else
    ! 1.7 times cf88 c12(a,g)o16 rate
       call rate_c12ag(tf,bden, &
                    rr % rates(1,ircag),rr % rates(2,ircag),dratedd1, &
                    rr % rates(1,iroga),rr % rates(2,iroga),dratedd2)
    endif

    ! triple alpha to c12
    call rate_tripalf(tf,bden, &
                      rr % rates(1,ir3a),rr % rates(2,ir3a),dratedd1, &
                      rr % rates(1,irg3a),rr % rates(2,irg3a),dratedd2)

    ! c12 + c12
    call rate_c12c12(tf,bden, &
                     rr % rates(1,ir1212),rr % rates(2,ir1212),dratedd1, &
                     rrate,drratedt,dratedd2)

    ! c12 + o16
    call rate_c12o16(tf,bden, &
                     rr % rates(1,ir1216),rr % rates(2,ir1216),dratedd1, &
                     rrate,drratedt,dratedd2)

    ! o16 + o16
    call rate_o16o16(tf,bden, &
                     rr % rates(1,ir1616),rr % rates(2,ir1616),dratedd1, &
                     rrate,drratedt,dratedd2)

    ! o16(a,g)ne20
    call rate_o16ag(tf,bden, &
                    rr % rates(1,iroag),rr % rates(2,iroag),dratedd1, &
                    rr % rates(1,irnega),rr % rates(2,irnega),dratedd2)

    ! ne20(a,g)mg24
    call rate_ne20ag(tf,bden, &
                     rr % rates(1,irneag),rr % rates(2,irneag),dratedd1, &
                     rr % rates(1,irmgga),rr % rates(2,irmgga),dratedd2)

    ! mg24(a,g)si28
    call rate_mg24ag(tf,bden, &
                     rr % rates(1,irmgag),rr % rates(2,irmgag),dratedd1, &
                     rr % rates(1,irsiga),rr % rates(2,irsiga),dratedd2)

    ! mg24(a,p)al27
    call rate_mg24ap(tf,bden, &
                     rr % rates(1,irmgap),rr % rates(2,irmgap),dratedd1, &
                     rr % rates(1,iralpa),rr % rates(2,iralpa),dratedd2)

    ! al27(p,g)si28
    call rate_al27pg(tf,bden, &
                     rr % rates(1,iralpg),rr % rates(2,iralpg),dratedd1, &
                     rr % rates(1,irsigp),rr % rates(2,irsigp),dratedd2)

    ! si28(a,g)s32
    call rate_si28ag(tf,bden, &
                     rr % rates(1,irsiag),rr % rates(2,irsiag),dratedd1, &
                     rr % rates(1,irsga),rr % rates(2,irsga),dratedd2)

    ! si28(a,p)p31
    call rate_si28ap(tf,bden, &
                     rr % rates(1,irsiap),rr % rates(2,irsiap),dratedd1, &
                     rr % rates(1,irppa),rr % rates(2,irppa),dratedd2)

    ! p31(p,g)s32
    call rate_p31pg(tf,bden, &
                    rr % rates(1,irppg),rr % rates(2,irppg),dratedd1, &
                    rr % rates(1,irsgp),rr % rates(2,irsgp),dratedd2)

    ! s32(a,g)ar36
    call rate_s32ag(tf,bden, &
                    rr % rates(1,irsag),rr % rates(2,irsag),dratedd1, &
                    rr % rates(1,irarga),rr % rates(2,irarga),dratedd2)

    ! s32(a,p)cl35
    call rate_s32ap(tf,bden, &
                    rr % rates(1,irsap),rr % rates(2,irsap),dratedd1, &
                    rr % rates(1,irclpa),rr % rates(2,irclpa),dratedd2)

    ! cl35(p,g)ar36
    call rate_cl35pg(tf,bden, &
                     rr % rates(1,irclpg),rr % rates(2,irclpg),dratedd1, &
                     rr % rates(1,irargp),rr % rates(2,irargp),dratedd2)

    ! ar36(a,g)ca40
    call rate_ar36ag(tf,bden, &
                     rr % rates(1,irarag),rr % rates(2,irarag),dratedd1, &
                     rr % rates(1,ircaga),rr % rates(2,ircaga),dratedd2)

    ! ar36(a,p)k39
    call rate_ar36ap(tf,bden, &
                     rr % rates(1,irarap),rr % rates(2,irarap),dratedd1, &
                     rr % rates(1,irkpa),rr % rates(2,irkpa),dratedd2)

    ! k39(p,g)ca40
    call rate_k39pg(tf,bden, &
                    rr % rates(1,irkpg),rr % rates(2,irkpg),dratedd1, &
                    rr % rates(1,ircagp),rr % rates(2,ircagp),dratedd2)

    ! ca40(a,g)ti44
    call rate_ca40ag(tf,bden, &
                     rr % rates(1,ircaag),rr % rates(2,ircaag),dratedd1, &
                     rr % rates(1,irtiga),rr % rates(2,irtiga),dratedd2)

    ! ca40(a,p)sc43
    call rate_ca40ap(tf,bden, &
                     rr % rates(1,ircaap),rr % rates(2,ircaap),dratedd1, &
                     rr % rates(1,irscpa),rr % rates(2,irscpa),dratedd2)

    ! sc43(p,g)ti44
    call rate_sc43pg(tf,bden, &
                     rr % rates(1,irscpg),rr % rates(2,irscpg),dratedd1, &
                     rr % rates(1,irtigp),rr % rates(2,irtigp),dratedd2)

    ! ti44(a,g)cr48
    call rate_ti44ag(tf,bden, &
                     rr % rates(1,irtiag),rr % rates(2,irtiag),dratedd1, &
                     rr % rates(1,ircrga),rr % rates(2,ircrga),dratedd2)

    ! ti44(a,p)v47
    call rate_ti44ap(tf,bden, &
                     rr % rates(1,irtiap),rr % rates(2,irtiap),dratedd1, &
                     rr % rates(1,irvpa),rr % rates(2,irvpa),dratedd2)

    ! v47(p,g)cr48
    call rate_v47pg(tf,bden, &
                    rr % rates(1,irvpg),rr % rates(2,irvpg),dratedd1, &
                    rr % rates(1,ircrgp),rr % rates(2,ircrgp),dratedd2)

    ! cr48(a,g)fe52
    call rate_cr48ag(tf,bden, &
                     rr % rates(1,ircrag),rr % rates(2,ircrag),dratedd1, &
                     rr % rates(1,irfega),rr % rates(2,irfega),dratedd2)

    ! cr48(a,p)mn51
    call rate_cr48ap(tf,bden, &
                     rr % rates(1,ircrap),rr % rates(2,ircrap),dratedd1, &
                     rr % rates(1,irmnpa),rr % rates(2,irmnpa),dratedd2)

    ! mn51(p,g)fe52
    call rate_mn51pg(tf,bden, &
                     rr % rates(1,irmnpg),rr % rates(2,irmnpg),dratedd1, &
                     rr % rates(1,irfegp),rr % rates(2,irfegp),dratedd2)

    ! fe52(a,g)ni56
    call rate_fe52ag(tf,bden, &
                     rr % rates(1,irfeag),rr % rates(2,irfeag),dratedd1, &
                     rr % rates(1,irniga),rr % rates(2,irniga),dratedd2)

    ! fe52(a,p)co55
    call rate_fe52ap(tf,bden, &
                     rr % rates(1,irfeap),rr % rates(2,irfeap),dratedd1, &
                     rr % rates(1,ircopa),rr % rates(2,ircopa),dratedd2)

    ! co55(p,g)ni56
    call rate_co55pg(tf,bden, &
                     rr % rates(1,ircopg),rr % rates(2,ircopg),dratedd1, &
                     rr % rates(1,irnigp),rr % rates(2,irnigp),dratedd2)

  end subroutine aprox13rat



  subroutine screen_aprox13(btemp, bden, y, rr)

    use screening_module, only: screen5, plasma_state, fill_plasma_state

    implicit none

    ! this routine computes the screening factors
    ! and applies them to the raw reaction rates,
    ! producing the final reaction rates used by the
    ! right hand sides and jacobian matrix elements

    real(rt), intent(in   ) :: btemp, bden
    real(rt), intent(in   ) :: y(nspec)
    type (rate_t),    intent(inout) :: rr

    integer          :: jscr
    real(rt) :: sc1a,sc1adt,sc2a,sc2adt,sc3a,sc3adt
    real(rt) :: sc1add,sc2add
!    real(rt) :: sc3add

    real(rt) :: denom,denomdt,zz
!    real(rt) :: denomdd,r1dd,s1dd,t1dd,u1dd,v1dd,w1dd,x1dd,y1dd

    real(rt) :: ratraw

    type (plasma_state) :: state

    !$gpu

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

    ratraw = rr % rates(1,ir3a)
    rr % rates(1,ir3a) = ratraw * sc3a
    rr % rates(2,ir3a) = rr % rates(2,ir3a)*sc3a + ratraw*sc3adt
    !dratdumdd(ir3a) = dratrawdd(ir3a)*sc3a + rr % rates(1,ir3a)*sc3add

    ratraw = rr % rates(1,irg3a)
    rr % rates(1,irg3a) = ratraw * sc3a
    rr % rates(2,irg3a) = rr % rates(2,irg3a)*sc3a + ratraw*sc3adt
    !dratdumdd(irg3a) = dratrawdd(irg3a)*sc3a + rr % rates(1,irg3a)*sc3add


    ! c12 to o16
    ! c12(a,g)o16
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ratraw = rr % rates(1,ircag)
    rr % rates(1,ircag)  = ratraw * sc1a
    rr % rates(2,ircag)  = rr % rates(2,ircag)*sc1a + ratraw*sc1adt
    !dratdumdd(ircag)  = dratrawdd(ircag)*sc1a + rr % rates(1,ircag)*sc1add

    ratraw = rr % rates(1,iroga)
    rr % rates(1,iroga)  = ratraw * sc1a
    rr % rates(2,iroga)  = rr % rates(2,iroga)*sc1a + ratraw*sc1adt
    !dratdumdd(iroga)  = dratrawdd(iroga)*sc1a + rr % rates(1,iroga)*sc1add


    ! c12 + c12
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ratraw = rr % rates(1,ir1212)
    rr % rates(1,ir1212) = ratraw * sc1a
    rr % rates(2,ir1212) = rr % rates(2,ir1212)*sc1a + ratraw*sc1adt
    !dratdumdd(ir1212) = dratrawdd(ir1212)*sc1a + rr % rates(1,ir1212)*sc1add


    ! c12 + o16
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ratraw = rr % rates(1,ir1216)
    rr % rates(1,ir1216) = ratraw * sc1a
    rr % rates(2,ir1216) = rr % rates(2,ir1216)*sc1a + ratraw*sc1adt
    !dratdumdd(ir1216) = dratrawdd(ir1216)*sc1a + rr % rates(1,ir1216)*sc1add


    ! o16 + o16
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ratraw = rr % rates(1,ir1616)
    rr % rates(1,ir1616) = ratraw * sc1a
    rr % rates(2,ir1616) = rr % rates(2,ir1616)*sc1a + ratraw*sc1adt
    !dratdumdd(ir1616) = dratrawdd(ir1616)*sc1a + rr % rates(1,ir1616)*sc1add


    ! o16 to ne20
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! o16(a,g)ne20
    ratraw = rr % rates(1,iroag)
    rr % rates(1,iroag) = ratraw * sc1a
    rr % rates(2,iroag) = rr % rates(2,iroag)*sc1a + ratraw*sc1adt
    !dratdumdd(iroag) = dratrawdd(iroag)*sc1a + rr % rates(1,iroag)*sc1add

    ratraw = rr % rates(1,irnega)
    rr % rates(1,irnega) = ratraw * sc1a
    rr % rates(2,irnega) = rr % rates(2,irnega)*sc1a + ratraw*sc1adt
    !dratdumdd(irnega) = dratrawdd(irnega)*sc1a + rr % rates(1,irnega)*sc1add


    ! ne20 to mg24
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! ne20(a,g)mg24
    ratraw = rr % rates(1,irneag)
    rr % rates(1,irneag) = ratraw * sc1a
    rr % rates(2,irneag) = rr % rates(2,irneag)*sc1a + ratraw*sc1adt
    !dratdumdd(irneag) = dratrawdd(irneag)*sc1a + rr % rates(1,irneag)*sc1add

    ratraw = rr % rates(1,irmgga)
    rr % rates(1,irmgga) = ratraw * sc1a
    rr % rates(2,irmgga) = rr % rates(2,irmgga)*sc1a + ratraw*sc1adt
    !dratdumdd(irmgga) = dratrawdd(irmgga)*sc1a + rr % rates(1,irmgga)*sc1add


    ! mg24 to si28
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! mg24(a,g)si28
    ratraw = rr % rates(1,irmgag)
    rr % rates(1,irmgag) = ratraw * sc1a
    rr % rates(2,irmgag) = rr % rates(2,irmgag)*sc1a + ratraw*sc1adt
    !dratdumdd(irmgag) = dratrawdd(irmgag)*sc1a + rr % rates(1,irmgag)*sc1add

    ratraw = rr % rates(1,irsiga)
    rr % rates(1,irsiga) = ratraw * sc1a
    rr % rates(2,irsiga) = rr % rates(2,irsiga)*sc1a + ratraw*sc1adt
    !dratdumdd(irsiga) = dratrawdd(irsiga)*sc1a + rr % rates(1,irsiga)*sc1add


    ! mg24(a,p)al27
    ratraw = rr % rates(1,irmgap)
    rr % rates(1,irmgap) = ratraw * sc1a
    rr % rates(2,irmgap) = rr % rates(2,irmgap)*sc1a + ratraw*sc1adt
    !dratdumdd(irmgap) = dratrawdd(irmgap)*sc1a + rr % rates(1,irmgap)*sc1add

    ratraw = rr % rates(1,iralpa)
    rr % rates(1,iralpa) = ratraw * sc1a
    rr % rates(2,iralpa) = rr % rates(2,iralpa)*sc1a + ratraw*sc1adt
    !dratdumdd(iralpa) = dratrawdd(iralpa)*sc1a + rr % rates(1,iralpa)*sc1add


    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ! al27(p,g)si28
    ratraw = rr % rates(1,iralpg)
    rr % rates(1,iralpg) = ratraw * sc1a
    rr % rates(2,iralpg) = rr % rates(2,iralpg)*sc1a + ratraw*sc1adt
    !dratdumdd(iralpg) = dratrawdd(iralpg)*sc1a + rr % rates(1,iralpg)*sc1add

    ratraw = rr % rates(1,irsigp)
    rr % rates(1,irsigp) = ratraw * sc1a
    rr % rates(2,irsigp) = rr % rates(2,irsigp)*sc1a + ratraw*sc1adt
    !dratdumdd(irsigp) = dratrawdd(irsigp)*sc1a + rr % rates(1,irsigp)*sc1add


    ! si28 to s32
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! si28(a,g)s32
    ratraw = rr % rates(1,irsiag)
    rr % rates(1,irsiag) = ratraw * sc1a
    rr % rates(2,irsiag) = rr % rates(2,irsiag)*sc1a + ratraw*sc1adt
    !dratdumdd(irsiag) = dratrawdd(irsiag)*sc1a + rr % rates(1,irsiag)*sc1add

    ratraw = rr % rates(1,irsga)
    rr % rates(1,irsga) = ratraw * sc1a
    rr % rates(2,irsga) = rr % rates(2,irsga)*sc1a + ratraw*sc1adt
    !dratdumdd(irsga) = dratrawdd(irsga)*sc1a + rr % rates(1,irsga)*sc1add


    ! si28(a,p)p31
    ratraw = rr % rates(1,irsiap)
    rr % rates(1,irsiap) = ratraw * sc1a
    rr % rates(2,irsiap) = rr % rates(2,irsiap)*sc1a + ratraw*sc1adt
    !dratdumdd(irsiap) = dratrawdd(irsiap)*sc1a + rr % rates(1,irsiap)*sc1add

    ratraw = rr % rates(1,irppa)
    rr % rates(1,irppa)  = ratraw * sc1a
    rr % rates(2,irppa)  = rr % rates(2,irppa)*sc1a  + ratraw*sc1adt
    !dratdumdd(irppa)  = dratrawdd(irppa)*sc1a  + rr % rates(1,irppa)*sc1add


    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ! p31(p,g)s32
    ratraw = rr % rates(1,irppg)
    rr % rates(1,irppg)  = ratraw * sc1a
    rr % rates(2,irppg)  = rr % rates(2,irppg)*sc1a + ratraw*sc1adt
    !dratdumdd(irppg)  = dratrawdd(irppg)*sc1a + rr % rates(1,irppg)*sc1add

    ratraw = rr % rates(1,irsgp)
    rr % rates(1,irsgp)  = ratraw * sc1a
    rr % rates(2,irsgp)  = rr % rates(2,irsgp)*sc1a + ratraw*sc1adt
    !dratdumdd(irsgp)  = dratrawdd(irsgp)*sc1a + rr % rates(1,irsgp)*sc1add


    ! s32 to ar36
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! s32(a,g)ar36
    ratraw = rr % rates(1,irsag)
    rr % rates(1,irsag)  = ratraw * sc1a
    rr % rates(2,irsag)  = rr % rates(2,irsag)*sc1a + ratraw*sc1adt
    !dratdumdd(irsag)  = dratrawdd(irsag)*sc1a + rr % rates(1,irsag)*sc1add

    ratraw = rr % rates(1,irarga)
    rr % rates(1,irarga)  = ratraw * sc1a
    rr % rates(2,irarga)  = rr % rates(2,irarga)*sc1a + ratraw*sc1adt
    !dratdumdd(irarga)  = dratrawdd(irarga)*sc1a + rr % rates(1,irarga)*sc1add

    ! s32(a,p)cl35
    ratraw = rr % rates(1,irsap)
    rr % rates(1,irsap)  = ratraw * sc1a
    rr % rates(2,irsap)  = rr % rates(2,irsap)*sc1a + ratraw*sc1adt
    !dratdumdd(irsap)  = dratrawdd(irsap)*sc1a + rr % rates(1,irsap)*sc1add

    ratraw = rr % rates(1,irclpa)
    rr % rates(1,irclpa) = ratraw * sc1a
    rr % rates(2,irclpa) = rr % rates(2,irclpa)*sc1a + ratraw*sc1adt
    !dratdumdd(irclpa) = dratrawdd(irclpa)*sc1a + rr % rates(1,irclpa)*sc1add


    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ! cl35(p,g)ar36
    ratraw = rr % rates(1,irclpg)
    rr % rates(1,irclpg) = ratraw * sc1a
    rr % rates(2,irclpg) = rr % rates(2,irclpg)*sc1a + ratraw*sc1adt
    !dratdumdd(irclpg) = dratrawdd(irclpg)*sc1a + rr % rates(1,irclpg)*sc1add

    ratraw = rr % rates(1,irargp)
    rr % rates(1,irargp) = ratraw * sc1a
    rr % rates(2,irargp) = rr % rates(2,irargp)*sc1a + ratraw*sc1adt
    !dratdumdd(irargp) = dratrawdd(irargp)*sc1a + rr % rates(1,irargp)*sc1add


    ! ar36 to ca40
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! ar36(a,g)ca40
    ratraw = rr % rates(1,irarag)
    rr % rates(1,irarag) = ratraw * sc1a
    rr % rates(2,irarag) = rr % rates(2,irarag)*sc1a + ratraw*sc1adt
    !dratdumdd(irarag) = dratrawdd(irarag)*sc1a + rr % rates(1,irarag)*sc1add

    ratraw = rr % rates(1,ircaga)
    rr % rates(1,ircaga) = ratraw * sc1a
    rr % rates(2,ircaga) = rr % rates(2,ircaga)*sc1a + ratraw*sc1adt
    !dratdumdd(ircaga) = dratrawdd(ircaga)*sc1a + rr % rates(1,ircaga)*sc1add


    ! ar36(a,p)k39
    ratraw = rr % rates(1,irarap)
    rr % rates(1,irarap) = ratraw * sc1a
    rr % rates(2,irarap) = rr % rates(2,irarap)*sc1a + ratraw*sc1adt
    !dratdumdd(irarap) = dratrawdd(irarap)*sc1a + rr % rates(1,irarap)*sc1add

    ratraw = rr % rates(1,irkpa)
    rr % rates(1,irkpa) = ratraw * sc1a
    rr % rates(2,irkpa) = rr % rates(2,irkpa)*sc1a  + ratraw*sc1adt
    !dratdumdd(irkpa)  = dratrawdd(irkpa)*sc1a  + rr % rates(1,irkpa)*sc1add


    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ! k39(p,g)ca40
    ratraw = rr % rates(1,irkpg)
    rr % rates(1,irkpg) = ratraw * sc1a
    rr % rates(2,irkpg) = rr % rates(2,irkpg)*sc1a  + ratraw*sc1adt
    !dratdumdd(irkpg)  = dratrawdd(irkpg)*sc1a  + rr % rates(1,irkpg)*sc1add

    ratraw = rr % rates(1,ircagp)
    rr % rates(1,ircagp) = ratraw * sc1a
    rr % rates(2,ircagp) = rr % rates(2,ircagp)*sc1a  + ratraw*sc1adt
    !dratdumdd(ircagp)  = dratrawdd(ircagp)*sc1a  + rr % rates(1,ircagp)*sc1add


    ! ca40 to ti44
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! ca40(a,g)ti44
    ratraw = rr % rates(1,ircaag)
    rr % rates(1,ircaag) = ratraw * sc1a
    rr % rates(2,ircaag) = rr % rates(2,ircaag)*sc1a + ratraw*sc1adt
    !dratdumdd(ircaag) = dratrawdd(ircaag)*sc1a + rr % rates(1,ircaag)*sc1add

    ratraw = rr % rates(1,irtiga)
    rr % rates(1,irtiga) = ratraw * sc1a
    rr % rates(2,irtiga) = rr % rates(2,irtiga)*sc1a + ratraw*sc1adt
    !dratdumdd(irtiga) = dratrawdd(irtiga)*sc1a + rr % rates(1,irtiga)*sc1add


    ! ca40(a,p)sc43
    ratraw = rr % rates(1,ircaap)
    rr % rates(1,ircaap) = ratraw * sc1a
    rr % rates(2,ircaap) = rr % rates(2,ircaap)*sc1a + ratraw*sc1adt
    !dratdumdd(ircaap) = dratrawdd(ircaap)*sc1a + rr % rates(1,ircaap)*sc1add

    ratraw = rr % rates(1,irscpa)
    rr % rates(1,irscpa) = ratraw * sc1a
    rr % rates(2,irscpa) = rr % rates(2,irscpa)*sc1a + ratraw*sc1adt
    !dratdumdd(irscpa) = dratrawdd(irscpa)*sc1a + rr % rates(1,irscpa)*sc1add


    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ! sc43(p,g)ti44
    ratraw = rr % rates(1,irscpg)
    rr % rates(1,irscpg) = ratraw * sc1a
    rr % rates(2,irscpg) = rr % rates(2,irscpg)*sc1a + ratraw*sc1adt
    !dratdumdd(irscpg) = dratrawdd(irscpg)*sc1a + rr % rates(1,irscpg)*sc1add

    ratraw = rr % rates(1,irtigp)
    rr % rates(1,irtigp) = ratraw * sc1a
    rr % rates(2,irtigp) = rr % rates(2,irtigp)*sc1a + ratraw*sc1adt
    !dratdumdd(irtigp) = dratrawdd(irtigp)*sc1a + rr % rates(1,irtigp)*sc1add


    ! ti44 to cr48
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! ti44(a,g)cr48
    ratraw = rr % rates(1,irtiag)
    rr % rates(1,irtiag) = ratraw * sc1a
    rr % rates(2,irtiag) = rr % rates(2,irtiag)*sc1a + ratraw*sc1adt
    !dratdumdd(irtiag) = dratrawdd(irtiag)*sc1a + rr % rates(1,irtiag)*sc1add

    ratraw = rr % rates(1,ircrga)
    rr % rates(1,ircrga) = ratraw * sc1a
    rr % rates(2,ircrga) = rr % rates(2,ircrga)*sc1a + ratraw*sc1adt
    !dratdumdd(ircrga) = dratrawdd(ircrga)*sc1a + rr % rates(1,ircrga)*sc1add

    ! ti44(a,p)v47
    ratraw = rr % rates(1,irtiap)
    rr % rates(1,irtiap) = ratraw * sc1a
    rr % rates(2,irtiap) = rr % rates(2,irtiap)*sc1a + ratraw*sc1adt
    !dratdumdd(irtiap) = dratrawdd(irtiap)*sc1a + rr % rates(1,irtiap)*sc1add

    ratraw = rr % rates(1,irvpa)
    rr % rates(1,irvpa) = ratraw * sc1a
    rr % rates(2,irvpa) = rr % rates(2,irvpa)*sc1a  + ratraw*sc1adt
    !dratdumdd(irvpa)  = dratrawdd(irvpa)*sc1a  + rr % rates(1,irvpa)*sc1add


    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ! v47(p,g)cr48
    ratraw = rr % rates(1,irvpg)
    rr % rates(1,irvpg) = ratraw * sc1a
    rr % rates(2,irvpg) = rr % rates(2,irvpg)*sc1a  + ratraw*sc1adt
    !dratdumdd(irvpg)  = dratrawdd(irvpg)*sc1a  + rr % rates(1,irvpg)*sc1add

    ratraw = rr % rates(1,ircrgp)
    rr % rates(1,ircrgp)  = ratraw * sc1a
    rr % rates(2,ircrgp)  = rr % rates(2,ircrgp)*sc1a  + ratraw*sc1adt
    !dratdumdd(ircrgp)  = dratrawdd(ircrgp)*sc1a  + rr % rates(1,ircrgp)*sc1add


    ! cr48 to fe52
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! cr48(a,g)fe52
    ratraw = rr % rates(1,ircrag)
    rr % rates(1,ircrag) = ratraw * sc1a
    rr % rates(2,ircrag) = rr % rates(2,ircrag)*sc1a + ratraw*sc1adt
    !dratdumdd(ircrag) = dratrawdd(ircrag)*sc1a + rr % rates(1,ircrag)*sc1add

    ratraw = rr % rates(1,irfega)
    rr % rates(1,irfega) = ratraw * sc1a
    rr % rates(2,irfega) = rr % rates(2,irfega)*sc1a + ratraw*sc1adt
    !dratdumdd(irfega) = dratrawdd(irfega)*sc1a + rr % rates(1,irfega)*sc1add


    ! cr48(a,p)mn51
    ratraw = rr % rates(1,ircrap)
    rr % rates(1,ircrap) = ratraw * sc1a
    rr % rates(2,ircrap) = rr % rates(2,ircrap)*sc1a + ratraw*sc1adt
    !dratdumdd(ircrap) = dratrawdd(ircrap)*sc1a + rr % rates(1,ircrap)*sc1add

    ratraw = rr % rates(1,irmnpa)
    rr % rates(1,irmnpa) = ratraw * sc1a
    rr % rates(2,irmnpa) = rr % rates(2,irmnpa)*sc1a + ratraw*sc1adt
    !dratdumdd(irmnpa) = dratrawdd(irmnpa)*sc1a + rr % rates(1,irmnpa)*sc1add


    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)

    ! mn51(p,g)fe52
    ratraw = rr % rates(1,irmnpg)
    rr % rates(1,irmnpg) = ratraw * sc1a
    rr % rates(2,irmnpg) = rr % rates(2,irmnpg)*sc1a + ratraw*sc1adt
    !dratdumdd(irmnpg) = dratrawdd(irmnpg)*sc1a + rr % rates(1,irmnpg)*sc1add

    ratraw = rr % rates(1,irfegp)
    rr % rates(1,irfegp) = ratraw * sc1a
    rr % rates(2,irfegp) = rr % rates(2,irfegp)*sc1a + ratraw*sc1adt
    !dratdumdd(irfegp) = dratrawdd(irfegp)*sc1a + rr % rates(1,irfegp)*sc1add


    ! fe52 to ni56
    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! fe52(a,g)ni56
    ratraw = rr % rates(1,irfeag)
    rr % rates(1,irfeag) = ratraw * sc1a
    rr % rates(2,irfeag) = rr % rates(2,irfeag)*sc1a + ratraw*sc1adt
    !dratdumdd(irfeag) = dratrawdd(irfeag)*sc1a + rr % rates(1,irfeag)*sc1add

    ratraw = rr % rates(1,irniga)
    rr % rates(1,irniga) = ratraw * sc1a
    rr % rates(2,irniga) = rr % rates(2,irniga)*sc1a + ratraw*sc1adt
    !dratdumdd(irniga) = dratrawdd(irniga)*sc1a + rr % rates(1,irniga)*sc1add


    ! fe52(a,p)co55
    ratraw = rr % rates(1,irfeap)
    rr % rates(1,irfeap) = ratraw * sc1a
    rr % rates(2,irfeap) = rr % rates(2,irfeap)*sc1a + ratraw*sc1adt
    !dratdumdd(irfeap) = dratrawdd(irfeap)*sc1a + rr % rates(1,irfeap)*sc1add

    ratraw = rr % rates(1,ircopa)
    rr % rates(1,ircopa) = ratraw * sc1a
    rr % rates(2,ircopa) = rr % rates(2,ircopa)*sc1a + ratraw*sc1adt
    !dratdumdd(ircopa) = dratrawdd(ircopa)*sc1a + rr % rates(1,ircopa)*sc1add


    jscr = jscr + 1
    call screen5(state,jscr,sc1a,sc1adt,sc1add)


    ! co55(p,g)ni56
    ratraw = rr % rates(1,ircopg)
    rr % rates(1,ircopg) = ratraw * sc1a
    rr % rates(2,ircopg) = rr % rates(2,ircopg)*sc1a + ratraw*sc1adt
    !dratdumdd(ircopg) = dratrawdd(ircopg)*sc1a + rr % rates(1,ircopg)*sc1add

    ratraw = rr % rates(1,irnigp)
    rr % rates(1,irnigp) = ratraw * sc1a
    rr % rates(2,irnigp) = rr % rates(2,irnigp)*sc1a + ratraw*sc1adt
    !dratdumdd(irnigp) = dratrawdd(irnigp)*sc1a + rr % rates(1,irnigp)*sc1add


    ! now form those lovely dummy proton link rates

    ! mg24(a,p)27al(p,g)28si
    rr % rates(1,irr1)  = 0.0e0_rt
    rr % rates(2,irr1)  = 0.0e0_rt
    !dratdumdd(irr1)  = 0.0e0_rt
    denom    = rr % rates(1,iralpa) + rr % rates(1,iralpg)
    denomdt  = rr % rates(2,iralpa) + rr % rates(2,iralpg)
    !denomdd  = dratdumdd(iralpa) + dratdumdd(iralpg)
    if (denom .gt. 1.0e-30_rt) then
       zz = 1.0e0_rt/denom
       rr % rates(1,irr1) = rr % rates(1,iralpa)*zz
       rr % rates(2,irr1) = (rr % rates(2,iralpa) - rr % rates(1,irr1)*denomdt)*zz
       !dratdumdd(irr1) = (dratdumdd(iralpa) - rr % rates(1,irr1)*denomdd)*zz
    end if

    ! si28(a,p)p31(p,g)s32
    rr % rates(1,irs1)  = 0.0e0_rt
    rr % rates(2,irs1)  = 0.0e0_rt
    !dratdumdd(irs1)  = 0.0e0_rt
    denom    = rr % rates(1,irppa) + rr % rates(1,irppg)
    denomdt  = rr % rates(2,irppa) + rr % rates(2,irppg)
    !denomdd  = dratdumdd(irppa) + dratdumdd(irppg)
    if (denom .gt. 1.0e-30_rt) then
       zz = 1.0e0_rt/denom
       rr % rates(1,irs1) = rr % rates(1,irppa)*zz
       rr % rates(2,irs1) = (rr % rates(2,irppa) - rr % rates(1,irs1)*denomdt)*zz
       !dratdumdd(irs1) = (dratdumdd(irppa) - rr % rates(1,irs1)*denomdd)*zz
    end if

    ! s32(a,p)cl35(p,g)ar36
    rr % rates(1,irt1)  = 0.0e0_rt
    rr % rates(2,irt1)  = 0.0e0_rt
    !dratdumdd(irt1)  = 0.0e0_rt
    denom    = rr % rates(1,irclpa) + rr % rates(1,irclpg)
    denomdt  = rr % rates(2,irclpa) + rr % rates(2,irclpg)
    !denomdd  = dratdumdd(irclpa) + dratdumdd(irclpg)
    if (denom .gt. 1.0e-30_rt) then
       zz = 1.0e0_rt/denom
       rr % rates(1,irt1) = rr % rates(1,irclpa)*zz
       rr % rates(2,irt1) = (rr % rates(2,irclpa) - rr % rates(1,irt1)*denomdt)*zz
       !dratdumdd(irt1) = (dratdumdd(irclpa) - rr % rates(1,irt1)*denomdd)*zz
    end if

    ! ar36(a,p)k39(p,g)ca40
    rr % rates(1,iru1)  = 0.0e0_rt
    rr % rates(2,iru1)  = 0.0e0_rt
    !dratdumdd(iru1)  = 0.0e0_rt
    denom    = rr % rates(1,irkpa) + rr % rates(1,irkpg)
    denomdt  = rr % rates(2,irkpa) + rr % rates(2,irkpg)
    !denomdd  = dratdumdd(irkpa) + dratdumdd(irkpg)
    if (denom .gt. 1.0e-30_rt) then
       zz   = 1.0e0_rt/denom
       rr % rates(1,iru1)   = rr % rates(1,irkpa)*zz
       rr % rates(2,iru1) = (rr % rates(2,irkpa) - rr % rates(1,iru1)*denomdt)*zz
       !dratdumdd(iru1) = (dratdumdd(irkpa) - rr % rates(1,iru1)*denomdd)*zz
    end if

    ! ca40(a,p)sc43(p,g)ti44
    rr % rates(1,irv1)  = 0.0e0_rt
    rr % rates(2,irv1)  = 0.0e0_rt
    !dratdumdd(irv1)  = 0.0e0_rt
    denom    = rr % rates(1,irscpa) + rr % rates(1,irscpg)
    denomdt  = rr % rates(2,irscpa) + rr % rates(2,irscpg)
    !denomdd  = dratdumdd(irscpa) + dratdumdd(irscpg)
    if (denom .gt. 1.0e-30_rt) then
       zz  = 1.0e0_rt/denom
       rr % rates(1,irv1) = rr % rates(1,irscpa)*zz
       rr % rates(2,irv1) = (rr % rates(2,irscpa) - rr % rates(1,irv1)*denomdt)*zz
       !dratdumdd(irv1) = (dratdumdd(irscpa) - rr % rates(1,irv1)*denomdd)*zz
    end if

    ! ti44(a,p)v47(p,g)cr48
    rr % rates(1,irw1) = 0.0e0_rt
    rr % rates(2,irw1) = 0.0e0_rt
    !dratdumdd(irw1) = 0.0e0_rt
    denom    = rr % rates(1,irvpa) + rr % rates(1,irvpg)
    denomdt  = rr % rates(2,irvpa) + rr % rates(2,irvpg)
    !denomdd  = dratdumdd(irvpa) + dratdumdd(irvpg)
    if (denom .gt. 1.0e-30_rt) then
       zz = 1.0e0_rt/denom
       rr % rates(1,irw1) = rr % rates(1,irvpa)*zz
       rr % rates(2,irw1) = (rr % rates(2,irvpa) - rr % rates(1,irw1)*denomdt)*zz
       !dratdumdd(irw1) = (dratdumdd(irvpa) - rr % rates(1,irw1)*denomdd)*zz
    end if

    ! cr48(a,p)mn51(p,g)fe52
    rr % rates(1,irx1) = 0.0e0_rt
    rr % rates(2,irx1) = 0.0e0_rt
    !dratdumdd(irx1) = 0.0e0_rt
    denom    = rr % rates(1,irmnpa) + rr % rates(1,irmnpg)
    denomdt  = rr % rates(2,irmnpa) + rr % rates(2,irmnpg)
    !denomdd  = dratdumdd(irmnpa) + dratdumdd(irmnpg)
    if (denom .gt. 1.0e-30_rt) then
       zz = 1.0e0_rt/denom
       rr % rates(1,irx1) = rr % rates(1,irmnpa)*zz
       rr % rates(2,irx1) = (rr % rates(2,irmnpa) - rr % rates(1,irx1)*denomdt)*zz
       !dratdumdd(irx1) = (dratdumdd(irmnpa) - rr % rates(1,irx1)*denomdd)*zz
    endif

    ! fe52(a,p)co55(p,g)ni56
    rr % rates(1,iry1) = 0.0e0_rt
    rr % rates(2,iry1) = 0.0e0_rt
    !dratdumdd(iry1) = 0.0e0_rt
    denom    = rr % rates(1,ircopa) + rr % rates(1,ircopg)
    denomdt  = rr % rates(2,ircopa) + rr % rates(2,ircopg)
    !denomdd  = dratdumdd(ircopa) + dratdumdd(ircopg)
    if (denom .gt. 1.0e-30_rt) then
       zz = 1.0e0_rt/denom
       rr % rates(1,iry1) = rr % rates(1,ircopa)*zz
       rr % rates(2,iry1) = (rr % rates(2,ircopa) - rr % rates(1,iry1)*denomdt)*zz
       !dratdumdd(iry1) = (dratdumdd(ircopa) - rr % rates(1,iry1)*denomdd)*zz
    end if

  end subroutine screen_aprox13



  subroutine dfdy_isotopes_aprox13(y,state,rr)

    use network
    use microphysics_math_module, only: esum3, esum4, esum5, esum20 ! function
    use jacobian_sparsity_module, only: set_jac_entry

    implicit none

    ! this routine sets up the dense aprox13 jacobian for the isotopes

    type (burn_t) :: state
    real(rt) :: y(nspec)
    type (rate_t)    :: rr

    real(rt) :: b(30)

    !$gpu

    ! he4 jacobian elements
    ! d(he4)/d(he4)
    b(1)  = -1.5e0_rt * y(ihe4) * y(ihe4) * rr % rates(1,ir3a)
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
    b(13) = -y(img24) * rr % rates(1,irmgap) * (1.0e0_rt-rr % rates(1,irr1))
    b(14) = -y(isi28) * rr % rates(1,irsiap) * (1.0e0_rt-rr % rates(1,irs1))
    b(15) = -y(is32)  * rr % rates(1,irsap)  * (1.0e0_rt-rr % rates(1,irt1))
    b(16) = -y(iar36) * rr % rates(1,irarap) * (1.0e0_rt-rr % rates(1,iru1))
    b(17) = -y(ica40) * rr % rates(1,ircaap) * (1.0e0_rt-rr % rates(1,irv1))
    b(18) = -y(iti44) * rr % rates(1,irtiap) * (1.0e0_rt-rr % rates(1,irw1))
    b(19) = -y(icr48) * rr % rates(1,ircrap) * (1.0e0_rt-rr % rates(1,irx1))
    b(20) = -y(ife52) * rr % rates(1,irfeap) * (1.0e0_rt-rr % rates(1,iry1))
    b(30) = esum20(b)
    call set_jac_entry(state, ihe4, ihe4, b(30))


    ! d(he4)/d(c12)
    b(1) =  y(ic12) * rr % rates(1,ir1212)
    b(2) =  0.5e0_rt * y(io16) * rr % rates(1,ir1216)
    b(3) =  3.0e0_rt * rr % rates(1,irg3a)
    b(4) = -y(ihe4) * rr % rates(1,ircag)
    b(30) = esum4(b)
    call set_jac_entry(state, ihe4, ic12, b(30))

    ! d(he4)/d(o16)
    b(1) =  0.5e0_rt * y(ic12) * rr % rates(1,ir1216)
    b(2) =  1.12e0_rt * 0.5e0_rt*y(io16) * rr % rates(1,ir1616)
    b(3) =  0.68e0_rt * rr % rates(1,irs1) * 0.5e0_rt*y(io16) * rr % rates(1,ir1616)
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
    b(3) = -y(ihe4) * rr % rates(1,irmgap) * (1.0e0_rt-rr % rates(1,irr1))
    b(30) = esum3(b)
    call set_jac_entry(state, ihe4, img24, b(30))

    ! d(he4)/d(si28)
    b(1) =  rr % rates(1,irsiga)
    b(2) = -y(ihe4) * rr % rates(1,irsiag)
    b(3) = -y(ihe4) * rr % rates(1,irsiap) * (1.0e0_rt-rr % rates(1,irs1))
    b(4) =  rr % rates(1,irr1) * rr % rates(1,irsigp)
    b(30) = esum4(b)
    call set_jac_entry(state, ihe4, isi28, b(30))

    ! d(he4)/d(s32)
    b(1) =  rr % rates(1,irsga)
    b(2) = -y(ihe4) * rr % rates(1,irsag)
    b(3) = -y(ihe4) * rr % rates(1,irsap) * (1.0e0_rt-rr % rates(1,irt1))
    b(4) =  rr % rates(1,irs1) * rr % rates(1,irsgp)
    b(30) = esum4(b)
    call set_jac_entry(state, ihe4, is32, b(30))

    ! d(he4)/d(ar36)
    b(1)  =  rr % rates(1,irarga)
    b(2)  = -y(ihe4) * rr % rates(1,irarag)
    b(3)  = -y(ihe4) * rr % rates(1,irarap) * (1.0e0_rt-rr % rates(1,iru1))
    b(4)  =  rr % rates(1,irt1) * rr % rates(1,irargp)
    b(30) = esum4(b)
    call set_jac_entry(state, ihe4, iar36, b(30))

    ! d(he4)/d(ca40)
    b(1) =  rr % rates(1,ircaga)
    b(2) = -y(ihe4) * rr % rates(1,ircaag)
    b(3) = -y(ihe4) * rr % rates(1,ircaap) * (1.0e0_rt-rr % rates(1,irv1))
    b(4) =  rr % rates(1,iru1) * rr % rates(1,ircagp)
    b(30) = esum4(b)
    call set_jac_entry(state, ihe4, ica40, b(30))

    ! d(he4)/d(ti44)
    b(1) =  rr % rates(1,irtiga)
    b(2) = -y(ihe4) * rr % rates(1,irtiag)
    b(3) = -y(ihe4) * rr % rates(1,irtiap) * (1.0e0_rt-rr % rates(1,irw1))
    b(4) =  rr % rates(1,irv1) * rr % rates(1,irtigp)
    b(30) = esum4(b)
    call set_jac_entry(state, ihe4, iti44, b(30))

    ! d(he4)/d(cr48)
    b(1) =  rr % rates(1,ircrga)
    b(2) = -y(ihe4) * rr % rates(1,ircrag)
    b(3) = -y(ihe4) * rr % rates(1,ircrap) * (1.0e0_rt-rr % rates(1,irx1))
    b(4) =  rr % rates(1,irw1) * rr % rates(1,ircrgp)
    b(30) = esum4(b)
    call set_jac_entry(state, ihe4, icr48, b(30))

    ! d(he4)/d(fe52)
    b(1) =  rr % rates(1,irfega)
    b(2) = -y(ihe4) * rr % rates(1,irfeag)
    b(3) = -y(ihe4) * rr % rates(1,irfeap) * (1.0e0_rt-rr % rates(1,iry1))
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
    b(1) =  0.5e0_rt * y(ihe4) * y(ihe4) * rr % rates(1,ir3a)
    b(2) = -y(ic12) * rr % rates(1,ircag)
    b(30) = sum(b(1:2))
    call set_jac_entry(state, ic12, ihe4, b(30))

    ! d(c12)/d(c12)
    b(1) = -2.0e0_rt * y(ic12) * rr % rates(1,ir1212)
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
    b(2) = -2.0e0_rt * y(io16) * rr % rates(1,ir1616)
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
    b(3) = -y(img24) * rr % rates(1,irmgap) * (1.0e0_rt-rr % rates(1,irr1))
    b(30) = esum3(b)
    call set_jac_entry(state, img24, ihe4, b(30))

    ! d(mg24)/d(c12)
    b(30) = 0.5e0_rt * y(io16) * rr % rates(1,ir1216)
    call set_jac_entry(state, img24, ic12, b(30))

    ! d(mg24)/d(o16)
    b(30) =  0.5e0_rt * y(ic12) * rr % rates(1,ir1216)
    call set_jac_entry(state, img24, io16, b(30))

    ! d(mg24)/d(ne20)
    b(30) = y(ihe4) * rr % rates(1,irneag)
    call set_jac_entry(state, img24, ine20, b(30))

    ! d(mg24)/d(mg24)
    b(1) = -y(ihe4) * rr % rates(1,irmgag)
    b(2) = -rr % rates(1,irmgga)
    b(3) = -y(ihe4) * rr % rates(1,irmgap) * (1.0e0_rt-rr % rates(1,irr1))
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
    b(3) =  y(img24) * rr % rates(1,irmgap) * (1.0e0_rt-rr % rates(1,irr1))
    b(4) = -y(isi28) * rr % rates(1,irsiap) * (1.0e0_rt-rr % rates(1,irs1))
    b(30) = esum4(b)
    call set_jac_entry(state, isi28, ihe4, b(30))

    ! d(si28)/d(c12)
    b(30) = 0.5e0_rt * y(io16) * rr % rates(1,ir1216)
    call set_jac_entry(state, isi28, ic12, b(30))

    ! d(si28)/d(o16)
    b(1) = 0.5e0_rt * y(ic12) * rr % rates(1,ir1216)
    b(2) = 1.12e0_rt * 0.5e0_rt*y(io16) * rr % rates(1,ir1616)
    b(3) = 0.68e0_rt * 0.5e0_rt*y(io16) * rr % rates(1,irs1) * rr % rates(1,ir1616)
    b(30) = esum3(b)
    call set_jac_entry(state, isi28, io16, b(30))

    ! d(si28)/d(mg24)
    b(1) =  y(ihe4) * rr % rates(1,irmgag)
    b(2) =  y(ihe4) * rr % rates(1,irmgap) * (1.0e0_rt-rr % rates(1,irr1))
    b(30) = sum(b(1:2))
    call set_jac_entry(state, isi28, img24, b(30))

    ! d(si28)/d(si28)
    b(1) =  -y(ihe4) * rr % rates(1,irsiag)
    b(2) = -rr % rates(1,irsiga)
    b(3) = -rr % rates(1,irr1) * rr % rates(1,irsigp)
    b(4) = -y(ihe4) * rr % rates(1,irsiap) * (1.0e0_rt-rr % rates(1,irs1))
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
    b(3) =  y(isi28) * rr % rates(1,irsiap) * (1.0e0_rt-rr % rates(1,irs1))
    b(4) = -y(is32) * rr % rates(1,irsap) * (1.0e0_rt-rr % rates(1,irt1))
    b(30) = esum4(b)
    call set_jac_entry(state, is32, ihe4, b(30))

    ! d(s32)/d(o16)
    b(1) = 0.68e0_rt*0.5e0_rt*y(io16)*rr % rates(1,ir1616)*(1.0e0_rt-rr % rates(1,irs1))
    b(2) = 0.2e0_rt * 0.5e0_rt*y(io16) * rr % rates(1,ir1616)
    b(30) = sum(b(1:2))
    call set_jac_entry(state, is32, io16, b(30))

    ! d(s32)/d(si28)
    b(1)  =y(ihe4) * rr % rates(1,irsiag)
    b(2) = y(ihe4) * rr % rates(1,irsiap) * (1.0e0_rt-rr % rates(1,irs1))
    b(30) = sum(b(1:2))
    call set_jac_entry(state, is32, isi28, b(30))

    ! d(s32)/d(s32)
    b(1) = -y(ihe4) * rr % rates(1,irsag)
    b(2) = -rr % rates(1,irsga)
    b(3) = -rr % rates(1,irs1) * rr % rates(1,irsgp)
    b(4) = -y(ihe4) * rr % rates(1,irsap) * (1.0e0_rt-rr % rates(1,irt1))
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
    b(3) =  y(is32)  * rr % rates(1,irsap) * (1.0e0_rt-rr % rates(1,irt1))
    b(4) = -y(iar36) * rr % rates(1,irarap) * (1.0e0_rt-rr % rates(1,iru1))
    b(30) = esum4(b)
    call set_jac_entry(state, iar36, ihe4, b(30))

    ! d(ar36)/d(s32)
    b(1) = y(ihe4) * rr % rates(1,irsag)
    b(2) = y(ihe4) * rr % rates(1,irsap) * (1.0e0_rt-rr % rates(1,irt1))
    b(30) = sum(b(1:2))
    call set_jac_entry(state, iar36, is32, b(30))

    ! d(ar36)/d(ar36)
    b(1) = -y(ihe4) * rr % rates(1,irarag)
    b(2) = -rr % rates(1,irarga)
    b(3) = -rr % rates(1,irt1) * rr % rates(1,irargp)
    b(4) = -y(ihe4) * rr % rates(1,irarap) * (1.0e0_rt-rr % rates(1,iru1))
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
    b(3)  =  y(iar36) * rr % rates(1,irarap)*(1.0e0_rt-rr % rates(1,iru1))
    b(4)  = -y(ica40) * rr % rates(1,ircaap)*(1.0e0_rt-rr % rates(1,irv1))
    b(30) = esum4(b)
    call set_jac_entry(state, ica40, ihe4, b(30))

    ! d(ca40)/d(ar36)
    b(1) =  y(ihe4) * rr % rates(1,irarag)
    b(2) =  y(ihe4) * rr % rates(1,irarap)*(1.0e0_rt-rr % rates(1,iru1))
    b(30) = sum(b(1:2))
    call set_jac_entry(state, ica40, iar36, b(30))

    ! d(ca40)/d(ca40)
    b(1) = -y(ihe4) * rr % rates(1, ircaag)
    b(2) = -rr % rates(1, ircaga)
    b(3) = -rr % rates(1, ircagp) * rr % rates(1, iru1)
    b(4) = -y(ihe4) * rr % rates(1, ircaap)*(1.0e0_rt-rr % rates(1, irv1))
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
    b(3) =  y(ica40) * rr % rates(1, ircaap)*(1.0e0_rt-rr % rates(1, irv1))
    b(4) = -y(iti44) * rr % rates(1, irtiap)*(1.0e0_rt-rr % rates(1, irw1))
    b(30) = esum4(b)
    call set_jac_entry(state, iti44, ihe4, b(30))

    ! d(ti44)/d(ca40)
    b(1) =  y(ihe4) * rr % rates(1, ircaag)
    b(2) =  y(ihe4) * rr % rates(1, ircaap)*(1.0e0_rt-rr % rates(1, irv1))
    b(30) = sum(b(1:2))
    call set_jac_entry(state, iti44, ica40, b(30))

    ! d(ti44)/d(ti44)
    b(1) = -y(ihe4) * rr % rates(1, irtiag)
    b(2) = -rr % rates(1, irtiga)
    b(3) = -rr % rates(1, irv1) * rr % rates(1, irtigp)
    b(4) = -y(ihe4) * rr % rates(1, irtiap)*(1.0e0_rt-rr % rates(1, irw1))
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
    b(3) =  y(iti44) * rr % rates(1, irtiap)*(1.0e0_rt-rr % rates(1, irw1))
    b(4) = -y(icr48) * rr % rates(1, ircrap)*(1.0e0_rt-rr % rates(1, irx1))
    b(30) = esum4(b)
    call set_jac_entry(state, icr48, ihe4, b(30))

    ! d(cr48)/d(ti44)
    b(1) =  y(ihe4) * rr % rates(1, irtiag)
    b(2) =  y(ihe4) * rr % rates(1, irtiap)*(1.0e0_rt-rr % rates(1, irw1))
    b(30) = sum(b(1:2))
    call set_jac_entry(state, icr48, iti44, b(30))

    ! d(cr48)/d(cr48)
    b(1) = -y(ihe4) * rr % rates(1, ircrag)
    b(2) = -rr % rates(1, ircrga)
    b(3) = -rr % rates(1, irw1) * rr % rates(1, ircrgp)
    b(4) = -y(ihe4) * rr % rates(1, ircrap)*(1.0e0_rt-rr % rates(1, irx1))
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
    b(3) =  y(icr48) * rr % rates(1, ircrap) * (1.0e0_rt-rr % rates(1, irx1))
    b(4) = -y(ife52) * rr % rates(1, irfeap) * (1.0e0_rt-rr % rates(1, iry1))
    b(30) = esum4(b)
    call set_jac_entry(state, ife52, ihe4, b(30))

    ! d(fe52)/d(cr48)
    b(1) = y(ihe4) * rr % rates(1, ircrag)
    b(2) = y(ihe4) * rr % rates(1, ircrap) * (1.0e0_rt-rr % rates(1, irx1))
    b(30) = sum(b(1:2))
    call set_jac_entry(state, ife52, icr48, b(30))

    ! d(fe52)/d(fe52)
    b(1) = -y(ihe4) * rr % rates(1, irfeag)
    b(2) = -rr % rates(1, irfega)
    b(3) = -rr % rates(1, irx1) * rr % rates(1, irfegp)
    b(4) = -y(ihe4) * rr % rates(1, irfeap) * (1.0e0_rt-rr % rates(1, iry1))
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
    b(2) =  y(ife52) * rr % rates(1, irfeap) * (1.0e0_rt-rr % rates(1, iry1))
    b(30) = sum(b(1:2))
    call set_jac_entry(state, ini56, ihe4, b(30))

    ! d(ni56)/d(fe52)
    b(1) = y(ihe4) * rr % rates(1, irfeag)
    b(2) = y(ihe4) * rr % rates(1, irfeap) * (1.0e0_rt-rr % rates(1, iry1))
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

    real(rt) :: dydt(nspec), enuc

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

    call add_screening_factor(zion(ihe4),aion(ihe4),4.0e0_rt,8.0e0_rt)

    call add_screening_factor(zion(ic12),aion(ic12),zion(ihe4),aion(ihe4))

    call add_screening_factor(zion(ic12),aion(ic12),zion(ic12),aion(ic12))

    call add_screening_factor(zion(ic12),aion(ic12),zion(io16),aion(io16))

    call add_screening_factor(zion(io16),aion(io16),zion(io16),aion(io16))

    call add_screening_factor(zion(io16),aion(io16),zion(ihe4),aion(ihe4))

    call add_screening_factor(zion(ine20),aion(ine20),zion(ihe4),aion(ihe4))

    call add_screening_factor(zion(img24),aion(img24),zion(ihe4),aion(ihe4))

    call add_screening_factor(13.0e0_rt,27.0e0_rt,1.0e0_rt,1.0e0_rt)

    call add_screening_factor(zion(isi28),aion(isi28),zion(ihe4),aion(ihe4))

    call add_screening_factor(15.0e0_rt,31.0e0_rt,1.0e0_rt,1.0e0_rt)

    call add_screening_factor(zion(is32),aion(is32),zion(ihe4),aion(ihe4))

    call add_screening_factor(17.0e0_rt,35.0e0_rt,1.0e0_rt,1.0e0_rt)

    call add_screening_factor(zion(iar36),aion(iar36),zion(ihe4),aion(ihe4))

    call add_screening_factor(19.0e0_rt,39.0e0_rt,1.0e0_rt,1.0e0_rt)

    call add_screening_factor(zion(ica40),aion(ica40),zion(ihe4),aion(ihe4))

    call add_screening_factor(21.0e0_rt,43.0e0_rt,1.0e0_rt,1.0e0_rt)

    call add_screening_factor(zion(iti44),aion(iti44),zion(ihe4),aion(ihe4))

    call add_screening_factor(23.0e0_rt,47.0e0_rt,1.0e0_rt,1.0e0_rt)

    call add_screening_factor(zion(icr48),aion(icr48),zion(ihe4),aion(ihe4))

    call add_screening_factor(25.0e0_rt,51.0e0_rt,1.0e0_rt,1.0e0_rt)

    call add_screening_factor(zion(ife52),aion(ife52),zion(ihe4),aion(ihe4))

    call add_screening_factor(27.0e0_rt,55.0e0_rt,1.0e0_rt,1.0e0_rt)

  end subroutine set_up_screening_factors

  subroutine update_unevolved_species(state)

    implicit none

    type (burn_t)    :: state

    !$gpu

  end subroutine update_unevolved_species

end module actual_rhs_module
