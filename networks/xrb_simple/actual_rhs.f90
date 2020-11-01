module actual_rhs_module

  use burn_type_module, only: burn_t, net_ienuc, neqs, njrows, njcols
  use rate_type_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

  subroutine actual_rhs_init()

    implicit none

  end subroutine actual_rhs_init



  subroutine actual_rhs(state, ydot)

    use amrex_constants_module, only: ZERO
    use network, only: nspec, aion, aion_inv
    use temperature_integration_module, only: temperature_rhs
    
    implicit none

    type (burn_t), intent(in) :: state
    real(rt)        , intent(inout) :: ydot(neqs)

    type (rate_t) :: rr

    real(rt), parameter :: T2T9 = 1.0e-9_rt

    real(rt) :: dens, t9
    real(rt) :: ymol(nspec)

    ydot = ZERO

    ! several thermo vars via rpar
    dens = state % rho
    ymol = state % xn * aion_inv
    t9   = state % T * T2T9

    ! build the rates
    call make_rates(t9, dens, ymol, state, rr)

    ! set up the ODEs
    call make_ydots(ymol, t9, state, ydot(1:nspec), rr)

    call ener_gener_rate(ydot(1:nspec), ydot(net_ienuc))

    call temperature_rhs(state, ydot)

  end subroutine actual_rhs



  ! TODO - make this an analytic jacobian
  subroutine actual_jac(state, jac)

    use amrex_constants_module, only: ZERO

    implicit none

    type (burn_t) :: state
    real(rt)        , intent(inout) :: jac(njrows, njcols)

    jac = ZERO

  end subroutine actual_jac



  subroutine make_rates(t9, dens, y, state, rr)

    use amrex_constants_module, only: ZERO, THIRD, ONE, SIX
    use actual_network, only: nspec, wk14o, wk15o, &
                              ir3a, irag15, irap14, irap18, irwk14o, irwk15o

    implicit none

    real(rt), intent(in   ) :: t9, dens, y(nspec)
    type (burn_t) :: state
    type (rate_t) :: rr

    real(rt) :: t9m13, t9m23, t9m1, t9m32, t9m3
    real(rt) :: t913, t923, t943, t953, t9113, t9log

    ! common temperature factors
    t9m1 = ONE / t9
    t9m13 = t9m1**THIRD
    t9m23 = t9m13*t9m13
    t9m32 = t9m1*sqrt(t9m1)
    t9m3 = t9m1*t9m1*t9m1
    t913 = t9**THIRD
    t923 = t913*t913
    t943 = t923*t923
    t953 = t943*t913
    t9113 = t953*t953*t913
    t9log = log(t9)

    ! zero the rates
    rr % rates(:,:) = ZERO

    !**********************************************************************
    ! Start the rate calculations
    ! TODO - temperature derivatives for use in analytic jacobian
    !**********************************************************************
    ! helium burning - has been divided by 6
    rr % rates(1,ir3a) = 2.79e-8_rt*t9m3*exp(-4.403_rt*t9m1)*dens*dens/SIX

    ! 15o(ag)19ne
    rr % rates(1,irag15) = (19.0_rt * (t9**2.85_rt) * &
         exp(-7.356_rt*t9m1) + &
         0.395_rt * t9m32 * exp(-5.849_rt*t9m1))*dens

    ! 14o(ap)17f
    ! an approx good below t9=0.5  is to drop first term, then
    !  rap14 = 3.31e+04_rt*t9m32*exp(-11.733*t9m1)
    !   1  + 1.79e+07_rt*t9m32*exp(-22.609/t9)+9000.*t9113*exp(-12.517/t9)
    rr % rates(1,irap14) = 1.68e13_rt*t9m23 * &
         exp(-39.388e0_rt*t9m13-(t9/0.717e0_rt)**2) &
         *(ONE+0.011e0_rt*t913+13.117e0_rt*t923+0.971e0_rt*t9+85.295e0_rt*t943 &
         +16.06e0_rt*t943)+3.31e4_rt*t9m32*exp(-11.733e0_rt*t9m1) &
         +1.79e7_rt*t9m32*exp(-22.609e0_rt*t9m1) &
         +9.00e+03_rt*t9113*exp(-12.517e0_rt*t9m1)
    rr % rates(1,irap14) = dens* rr % rates(1,irap14)

    ! 18ne(ap)21na
    rr % rates(1,irap18) = exp(56.59577e0_rt-2.006856e0_rt*t9m1 &
         +26.05828e0_rt*t9m13 &
         -87.88732e0_rt*t913 + 3.718809e0_rt*t9 - 0.1767444e0_rt*t953 &
         + 46.971960e0_rt*t9log)
    rr % rates(1,irap18) = dens * rr % rates(1,irap18)

    ! weak rates
    rr % rates(1,irwk14o) = wk14o
    rr % rates(1,irwk15o) = wk15o

  end subroutine make_rates



  subroutine make_ydots(ymol, t9, state, dydt, rr)

    use amrex_constants_module, only: ZERO, TWO, THREE, SIX
    use actual_network, only: nspec, io14, io15, ine18, isi25, ihe4, ih1, ife56, &
                              ir3a, irag15, irap14, irap18, irwk14o, irwk15o

    implicit none

    real(rt), intent(in   ) :: ymol(nspec), t9
    type (burn_t), intent(in) :: state
    type (rate_t), intent(in) :: rr
    real(rt), intent(  out) :: dydt(nspec)

    real(rt) :: dens

    ! initialize
    dydt = ZERO
    dens = state % rho

    ! From Stan:
    !       Reaction                   Rate
    !         3a + 2p   --> 14O          3a
    !        14O +      --> 18Ne       14O(ap)17F
    !     15O + a + 6p  --> 25Si       15O(ag)19Ne
    !   18Ne + a + 3p   --> 25Si       18Ne(ap)21Na
    !         14O + p   --> 15O        14O(e+nu)14N
    !         15O + 3p  --> 14O + a    15O(e+nu)15N

    ! o14
    dydt(io14) = ymol(ihe4)**3 * rr % rates(1,ir3a) &
         + ymol(io15) * rr % rates(1,irwk15o)       &
         - ymol(io14) * ymol(ihe4) * rr % rates(1,irap14)

    ! o15
    dydt(io15) = ymol(io14) * rr % rates(1,irwk14o) &
         - ymol(io15) * ymol(ihe4) * rr % rates(1,irag15) &
         - ymol(io15) * rr % rates(1,irwk15o)

    ! ne18
    dydt(ine18) = ymol(io14) * ymol(ihe4) * rr % rates(1,irap14) &
         - ymol(ine18) * ymol(ihe4) * rr % rates(1,irap18)

    ! si25
    dydt(isi25) = ymol(io15) * ymol(ihe4) * rr % rates(1,irag15) &
         + ymol(ine18) * ymol(ihe4) * rr % rates(1,irap18)

    ! he4
    dydt(ihe4) = ymol(io15) * rr % rates(1,irwk15o) &
         - THREE * ymol(ihe4)**3 * rr % rates(1,ir3a) &
         - ymol(io14) * ymol(ihe4) * rr % rates(1,irap14) &
         - ymol(io15) * ymol(ihe4) * rr % rates(1,irag15) &
         - ymol(ine18) * ymol(ihe4) * rr % rates(1,irap18)

    ! h1
    dydt(ih1) = - TWO * ymol(ihe4)**3 * rr % rates(1,ir3a) &
         - SIX * ymol(io15) * ymol(ihe4) * rr % rates(1,irag15) &
         - THREE * ymol(ine18) * ymol(ihe4) * rr % rates(1,irap18) &
         - ymol(io14) * rr % rates(1,irwk14o) &
         - THREE * ymol(io15) * rr % rates(1,irwk15o)

    ! iron is just a tracer
    dydt(ife56) = ZERO

  end subroutine make_ydots



  subroutine ener_gener_rate(dydt, enuc)

    use actual_network, only: nspec, aion, bion

    implicit none

    real(rt)         :: dydt(nspec), enuc

    enuc = -sum(dydt(:) * aion * bion(1:nspec))

  end subroutine ener_gener_rate

end module actual_rhs_module
