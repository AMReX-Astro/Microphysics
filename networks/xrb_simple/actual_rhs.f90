module actual_rhs_module

  use burn_type_module, only: burn_t, net_ienuc
  use rate_type_module

  implicit none

contains

  subroutine actual_rhs_init()

    implicit none

  end subroutine actual_rhs_init



  subroutine actual_rhs(state)

    use bl_types, only: dp_t
    use bl_constants_module, only: ZERO
    use actual_network, only: nspec, aion
    use temperature_integration_module, only: temperature_rhs
    
    implicit none

    type (burn_t) :: state
    type (rate_t) :: rr

    real(kind=dp_t), parameter :: T2T9 = 1.0e-9_dp_t

    real(kind=dp_t) :: dens, t9
    real(kind=dp_t) :: ymol(nspec), ydot(nspec)

    ydot = ZERO

    ! several thermo vars via rpar
    dens = state % rho
    ymol = state % xn / aion
    t9   = state % T * T2T9

    ! build the rates
    call make_rates(t9, dens, ymol, state, rr)

    ! set up the ODEs
    call make_ydots(ymol, t9, state, ydot(1:nspec), rr)

    state % ydot(1:nspec) = ydot

    call ener_gener_rate(state % ydot(1:nspec), state % ydot(net_ienuc))

    call temperature_rhs(state)

  end subroutine actual_rhs



  ! TODO - make this an analytic jacobian
  subroutine actual_jac(state)

    use bl_types, only: dp_t
    use bl_constants_module, only: ZERO

    implicit none

    type (burn_t) :: state

    state % jac = ZERO

  end subroutine actual_jac



  subroutine make_rates(t9, dens, y, state, rr)

    use bl_types, only: dp_t
    use bl_constants_module, only: ZERO, THIRD, ONE, SIX
    use actual_network, only: nspec, wk14o, wk15o, &
                              ir3a, irag15, irap14, irap18, irwk14o, irwk15o

    implicit none

    real(kind=dp_t), intent(in   ) :: t9, dens, y(nspec)
    type (burn_t) :: state
    type (rate_t) :: rr

    real(kind=dp_t) :: t9m13, t9m23, t9m1, t9m32, t9m3
    real(kind=dp_t) :: t913, t923, t943, t953, t9113, t9log

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
    rr % rates(1,ir3a) = 2.79d-8*t9m3*exp(-4.403_dp_t*t9m1)*dens*dens/SIX

    ! 15o(ag)19ne
    rr % rates(1,irag15) = (19.0_dp_t * (t9**2.85_dp_t) * &
         exp(-7.356_dp_t*t9m1) + &
         0.395_dp_t * t9m32 * exp(-5.849_dp_t*t9m1))*dens

    ! 14o(ap)17f
    ! an approx good below t9=0.5  is to drop first term, then
    !  rap14 = 3.31d+04*t9m32*exp(-11.733*t9m1)
    !   1  + 1.79d+07*t9m32*exp(-22.609/t9)+9000.*t9113*exp(-12.517/t9)
    rr % rates(1,irap14) = 1.68d13*t9m23 * &
         exp(-39.388d0*t9m13-(t9/0.717d0)**2) &
         *(ONE+0.011d0*t913+13.117d0*t923+0.971d0*t9+85.295d0*t943 &
         +16.06d0*t943)+3.31d4*t9m32*exp(-11.733d0*t9m1) &
         +1.79d7*t9m32*exp(-22.609d0*t9m1) &
         +9.00d+03*t9113*exp(-12.517d0*t9m1)
    rr % rates(1,irap14) = dens* rr % rates(1,irap14)

    ! 18ne(ap)21na
    rr % rates(1,irap18) = exp(56.59577d0-2.006856d0*t9m1 &
         +26.05828d0*t9m13 &
         -87.88732d0*t913 + 3.718809d0*t9 - 0.1767444d0*t953 &
         + 46.971960d0*t9log)
    rr % rates(1,irap18) = dens * rr % rates(1,irap18)

    ! weak rates
    rr % rates(1,irwk14o) = wk14o
    rr % rates(1,irwk15o) = wk15o

  end subroutine make_rates



  subroutine make_ydots(ymol, t9, state, dydt, rr)

    use bl_types, only: dp_t
    use bl_constants_module, only: ZERO, TWO, THREE, SIX
    use actual_network, only: nspec, io14, io15, ine18, isi25, ihe4, ih1, ife56, &
                              ir3a, irag15, irap14, irap18, irwk14o, irwk15o

    implicit none

    real(kind=dp_t), intent(in   ) :: ymol(nspec), t9
    type (burn_t), intent(in) :: state
    type (rate_t), intent(in) :: rr
    real(kind=dp_t), intent(  out) :: dydt(nspec)

    real(kind=dp_t) :: dens

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

    double precision :: dydt(nspec), enuc

    enuc = -sum(dydt(:) * aion * bion(1:nspec))

  end subroutine ener_gener_rate

  subroutine update_unevolved_species(state)

    !$acc routine seq

    implicit none

    type (burn_t)    :: state

  end subroutine update_unevolved_species

end module actual_rhs_module
