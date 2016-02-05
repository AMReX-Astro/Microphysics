module actual_rhs_module

  use actual_burner_data
  use burn_type_module
  use temperature_integration_module, only: temperature_rhs, temperature_jac

  implicit none

contains

  subroutine actual_rhs(state)

    use bl_types
    use bl_constants_module
    use network
    use rates_module

    implicit none

    type (burn_t)    :: state

    double precision :: dens, t9, y(nspec)

    y = state % xn / aion

    state % ydot = ZERO

    dens = state % rho
    t9   = state % T / 10.0**9

    ! build the rates; weak rates are the wk* variables
    call make_rates(t9, dens, y(1:nspec), state)

    ! set up the ODEs for the species
    call make_ydots(y(1:nspec),t9,state,state % ydot(1:nspec),.false.)

    ! Energy release
    state % ydot(net_ienuc) = -sum( ebin(:) * state % ydot(1:nspec) )

    ! Temperature ODE
    call temperature_rhs(state)

  end subroutine actual_rhs



  subroutine make_rates(t9,dens,y,state)

    use bl_types
    use bl_constants_module
    use rates_module
    use network
    use burner_module

    implicit none

    double precision :: t9, dens, y(nspec)
    type (burn_t)    :: state

    ! locally used rates
    double precision :: rate,dratedt,wk18ne,wk19ne
    double precision :: r56pg,dr56pgdt,cutoni,dcutonidt,r57decay,dr56eff,ddr56effdt
    double precision :: t9i32

    ! some numbers from appendix C in WW81; these should probably be
    ! updated with current rates
    double precision, parameter :: Lweak = 1.05d0, & ! this is for NS
                                !  Lweak = 0.107d0, & ! this is for lower 
                                                      ! densities
                                   la2 = ONE/FIFTEEN ! mean rate from 30s to 56ni
                                                     ! from p-capture and beta 
                                                     ! decays

    type (temp_t) :: tfactors

    state % rates(:,:) = ZERO ! Zero out rates

    tfactors = calc_tfactors(t9)

    ! some common parameters
    state % rates(1,irLweak) = Lweak
    state % rates(1,irla2)   = la2

    ! weak rates first
    !
    ! 14o(beta nu)14n
    call rate_o14_to_n14(tfactors,rate,dratedt)
    state % rates(1,irwk14o) = rate
    ! 15o(beta nu)15n
    call rate_o15_to_n15(tfactors,rate,dratedt)
    state % rates(1,irwk15o) = rate
    ! 17f(beta nu)17o
    call rate_f17_to_o17(tfactors,rate,dratedt)
    state % rates(1,irwk17f) = rate
    ! these weak rates aren't needed outside of this routine
    ! 18ne(beta nu)18f
    call rate_ne18_to_f18(tfactors,wk18ne,dratedt) 
    ! 19ne(beta nu)19f
    call rate_ne19_to_f19(tfactors,wk19ne,dratedt)

    ! 12c(p,g)13n
    call rate_p_c12_to_n13(tfactors,rate,dratedt)
    state % rates(1,irpg12c) = dens*rate
    state % rates(2,irpg12c) = dens*dratedt

    ! triple alpha
    call rate_he4_he4_he4_to_c12(tfactors,rate,dratedt)
    state % rates(1,ir3a) = dens*dens*rate
    state % rates(2,ir3a) = dens*dens*dratedt

    ! 17f(p,g)18ne
    call rate_p_f17_to_ne18(tfactors,rate,dratedt)
    state % rates(1,irpg17f) = dens*rate
    state % rates(2,irpg17f) = dens*dratedt

    ! 17f(g,p)16o
    call rate_f17_to_p_o16(tfactors,rate,dratedt)
    state % rates(1,irgp17f) = rate
    state % rates(2,irgp17f) = dratedt

    ! 15o(a,g)19ne
    call rate_he4_o15_to_ne19(tfactors,rate,dratedt)
    state % rates(1,irag15o) = dens*rate
    state % rates(2,irag15o) = dens*dratedt

    ! 16o(a,g)20ne
    call rate_he4_o16_to_ne20(tfactors,rate,dratedt)
    state % rates(1,irag16o) = dens*rate
    state % rates(2,irag16o) = dens*dratedt

    ! 16o(p,g)17f
    call rate_p_o16_to_f17(tfactors,rate,dratedt)
    state % rates(1,irpg16o) = dens*rate
    state % rates(2,irpg16o) = dens*dratedt

    ! 14o(a,p)17f
    call rate_he4_o14_to_p_f17(tfactors,rate,dratedt)
    state % rates(1,irap14o) = dens*rate
    state % rates(2,irap14o) = dens*dratedt

    ! limit CNO as minimum between 14n(p,g)15o and 15o(beta nu)15n
    ! we store the limited rate in irlambCNO; this is lambda_CNO in WW81
    call rate_p_n14_to_o15(tfactors,rate,dratedt)
    state % rates(1,irlambCNO) = min(state % rates(1,irwk15o),rate*dens*y(ih1))
    if (state % rates(1,irlambCNO) < state % rates(1,irwk15o)) then
       state % rates(2,irlambCNO) = dens*y(ih1)*dratedt
       state % rates(3,dlambCNOdh1) = rate*dens
    endif

    ! 22mg(...)30s
    ! check if this proceeds via p-captures or (a,p) reactions
    ! the Lweak is from WW81, eqn C15
    ! we store the rate in irlambda1; this is the lambda1 in WW81
    call rate_he4_si26_to_p_p29(tfactors,rate,dratedt)
    state % rates(1,irlambda1) = max(state % rates(1,irLweak),dens*y(ihe4)*rate)
    if (state % rates(1,irlambda1) > state % rates(1,irLweak)) then
       state % rates(2,irlambda1) = dens*y(ihe4)*dratedt
       state % rates(3,dlambda1dhe4) = dens*rate
       ! use the sign of state % rates(1,irlambda1) to indicate the value of delta1 in WW81
       ! if delta1 = 1, then we multiply the rate by -1
       state % rates(1,irlambda1) = -ONE*state % rates(1,irlambda1)
    endif

    ! 30s(...) 56ni 
    ! check if this proceeds via p-captures or (a,p) reactions
    ! use 44ti(a,p)v47 as a typical limiting rate for the (a,p) process
    ! store this in irlambda2; this is lambda2 in WW81
    call rate_he4_ti44_to_p_v47(tfactors,rate,dratedt)
    state % rates(1,irlambda2) = max(state % rates(1,irla2),dens*y(ihe4)*rate)
    if (state % rates(1,irlambda2) > state % rates(1,irla2)) then
       state % rates(2,irlambda2) = dens*y(ihe4)*dratedt
       state % rates(3,dlambda2dhe4) = dens*rate
       ! use the sign of state % rates(1,irlambda2) to indicate th value of delta2
       ! if delta2 = 1, then we multiply the rate by -1
       state % rates(1,irlambda2) = -ONE*state % rates(1,irlambda2)
    endif

    ! form s1 from WW81; branching ratio for 18ne beta decay (wk18ne) vs (a,p)
    ! store result in irs1
    ! 18ne(a,p)21na
    call rate_he4_ne18_to_p_na21(tfactors,rate,dratedt)
    state % rates(1,irs1) = wk18ne / (wk18ne + dens*y(ihe4)*rate)
    state % rates(2,irs1) = -state % rates(1,irs1)*dens*y(ihe4)*dratedt &
         / (wk18ne + dens*y(ihe4)*rate)
    state % rates(3,drs1dhe4) = -state % rates(1,irs1)*dens*rate &
         / (wk18ne + dens*y(ihe4)*rate)

    ! form r1 from WW81; ranching ratio for 19ne beta decay (wk19ne) vs (p,g)
    ! store result in irr1
    ! 19ne(p,g)20na
    call rate_p_ne19_to_na20(tfactors,rate,dratedt)
    state % rates(1,irr1) = wk19ne / (wk19ne + dens*y(ih1)*rate)
    state % rates(2,irr1) = -state % rates(1,irr1)*dens*y(ih1)*dratedt &
         / (wk19ne + dens*y(ih1)*rate)
    state % rates(3,drr1dh1) = -state % rates(1,irr1)*dens*rate &
         / (wk19ne + dens*y(ih1)*rate)


    !....
    !....  additional coding for proton capture on 56ni to heavier elements
    !....   kludge    56ni+56p -> 2 (56ni) at a rate given by min
    !....   of 56ni(pg) and 57cu decay rate
    !....
    !....  use 56ni rate from wallace and woosley 1981
    t9i32=tfactors%t9i*sqrt(tfactors%t9i)
    r56pg=dens*(1.29e-02_dp_t*exp(-4.897_dp_t*tfactors%t9i) &
         +7.065e+03_dp_t*exp(-20.33_dp_t*tfactors%t9i))*t9i32
    dr56pgdt = -THREE*HALF*r56pg*tfactors%t9i + &
         dens*t9i32*tfactors%t9i*tfactors%t9i* &
         (4.897_dp_t*1.29e-2_dp_t*exp(-4.897_dp_t*tfactors%t9i) &
         +20.33_dp_t*7.065e3_dp_t*exp(-20.33_dp_t*tfactors%t9i))
    !....  use generic proton separation energy of 400 kev
    !....  8.02 -> 4.64
    !      cutoni=2.08d-10*dens*exp(8.02*t9m1)/t932
    cutoni=2.08e-10_dp_t*dens*exp(4.642_dp_t*tfactors%t9i)*t9i32
    dcutonidt = cutoni*tfactors%t9i*(-THREE*HALF - 4.642_dp_t*tfactors%t9i)
    r57decay=3.54_dp_t
    dr56eff=min(r56pg,cutoni*r57decay)
    !   state % rates(3,r56eff) = d56eff
    !   if (d56eff < r56pg) state % rates(3,dr56effdt) = r57decay*dcutonidt
    state % rates(3,r56eff) = ZERO
    state % rates(3,dr56effdt) = ZERO

  end subroutine make_rates



  subroutine make_ydots(ymol,t9,state,dydt,doing_dratesdt)

    use bl_types
    use bl_constants_module
    use network

    implicit none

    double precision, intent(IN   ) :: ymol(nspec), t9
    logical ,         intent(IN   ) :: doing_dratesdt
    type(burn_t),     intent(INOUT) :: state
    double precision, intent(  OUT) :: dydt(nspec)

    integer          :: rate_idx
    double precision :: ddelta1, ddelta2
    double precision :: dens

    ! initialize
    dydt = ZERO
    dens = state % rho

    ! check to see if we are doing this with the t-derivatives
    ! if so, offset our starting index in the rate groups
    rate_idx = 3
    if (doing_dratesdt) rate_idx = 2
    

    if (.not. doing_dratesdt) then
       ddelta1 = ZERO
       ddelta2 = ZERO
       ! figure out the delta's; we used negative rates to indicate delta=1
       if (state % rates(1,irlambda1) < ZERO) then
          ddelta1 = ONE
          state % rates(1,irlambda1) = -ONE*state % rates(1,irlambda1)
       endif
       if (state % rates(1,irlambda2) < ZERO) then
          ddelta2 = ONE
          state % rates(1,irlambda2) = -ONE*state % rates(1,irlambda2)
       endif
       state % rates(3,delta1) = ddelta1
       state % rates(3,delta2) = ddelta2
    endif

    ! setup ODEs
    !
    !....
    !.... 12c = 1
    !....
    dydt(ic12)=-ymol(ic12)*ymol(ih1)*state % rates(rate_idx,irpg12c) &
         +ymol(ihe4)**3*state % rates(rate_idx,ir3a)/SIX &
         +ymol(io15)*state % rates(rate_idx,irlambCNO)
    !....
    !.... 14o = 2
    !....
    dydt(io14)=-ymol(io14)*ymol(ihe4)*state % rates(rate_idx,irap14o) &
         -ymol(io14)*state % rates(rate_idx,irwk14o) &
         +ymol(ic12)*ymol(ih1)*state % rates(rate_idx,irpg12c)
    !....
    !.... 15o = 3
    !....
    dydt(io15)=ymol(io14)*state % rates(rate_idx,irwk14o) &
         -ymol(io15)*ymol(ihe4)*state % rates(rate_idx,irag15o) &
         -ymol(io15)*state % rates(rate_idx,irlambCNO) &
         +ymol(if17)*ymol(ih1)*state % rates(rate_idx,irpg17f)*state % rates(rate_idx,irs1) &
         +ymol(if17)*state % rates(rate_idx,irwk17f)
    !....
    !.... 16o = 4
    !....
    dydt(io16) = ymol(if17)*state % rates(rate_idx,irgp17f) &
         -ymol(io16)*ymol(ih1)*state % rates(rate_idx,irpg16o) &
         +ymol(io15)*ymol(ihe4)*state % rates(rate_idx,irr1)*state % rates(rate_idx,irag15o) &
         -ymol(io16)*ymol(ihe4)*state % rates(rate_idx,irag16o)
    !....
    !.... 17f = 5
    !....
    dydt(if17)=ymol(io14)*ymol(ihe4)*state % rates(rate_idx,irap14o) &
         +ymol(io16)*ymol(ih1)*state % rates(rate_idx,irpg16o) &
         -ymol(if17)*state % rates(rate_idx,irgp17f) &
         -ymol(if17)*ymol(ih1)*state % rates(rate_idx,irpg17f) &
         -ymol(if17)*state % rates(rate_idx,irwk17f)
    !....
    !.... 22mg = 6
    !....
    dydt(img22)=ymol(io16)*ymol(ihe4)*state % rates(rate_idx,irag16o) &
         +ymol(if17)*ymol(ih1)*state % rates(rate_idx,irpg17f)*(ONE-state % rates(rate_idx,irs1)) &
         +ymol(io15)*ymol(ihe4)*state % rates(rate_idx,irag15o)*(ONE-state % rates(rate_idx,irr1)) &
         -ymol(img22)*state % rates(rate_idx,irlambda1)
    !....
    !.... 30s = 7
    !....
    dydt(is30)=ymol(img22)*state % rates(rate_idx,irlambda1) &
         -ymol(is30)*state % rates(rate_idx,irlambda2)
    !....
    !.... amax (56ni) = 8  (note that WW81 have a typo -- they write lambda1 here)
    !....
    dydt(ini56)=ymol(is30)*state % rates(rate_idx,irlambda2)
    !....
    !.... 4he (alpha) = 9
    !....
    dydt(ihe4)=-ymol(ihe4)**3*HALF*state % rates(rate_idx,ir3a) &
         +ymol(io15)*state % rates(rate_idx,irlambCNO) &
         -ymol(io14)*ymol(ihe4)*state % rates(rate_idx,irap14o) &
         +ymol(if17)*ymol(ih1)*state % rates(rate_idx,irpg17f)*state % rates(rate_idx,irs1) &
         -ymol(io15)*ymol(ihe4)*state % rates(rate_idx,irag15o)*(ONE-state % rates(rate_idx,irr1)) &
         -ymol(io16)*ymol(ihe4)*state % rates(rate_idx,irag16o) &
         -ymol(if17)*ymol(ih1)*state % rates(rate_idx,irpg17f)*(ONE-state % rates(rate_idx,irs1)) &
         +ymol(if17)*state % rates(rate_idx,irwk17f) &
         -TWO*ymol(img22)*state % rates(rate_idx,irlambda1)*state % rates(3,delta1) &
         -6.5d0*ymol(is30)*state % rates(rate_idx,irlambda2)*state % rates(3,delta2)
    !....
    !.... 1h (p) = 10
    !....
    dydt(ih1)=-ymol(io14)*state % rates(rate_idx,irwk14o) &
         -ymol(io15)*state % rates(rate_idx,irlambCNO) &
         -TWO*ymol(ic12)*ymol(ih1)*state % rates(rate_idx,irpg12c) &
         +ymol(io14)*ymol(ihe4)*state % rates(rate_idx,irap14o) &
         -TWO*ymol(if17)*ymol(ih1)*state % rates(rate_idx,irpg17f)*state % rates(rate_idx,irs1) &
         +ymol(if17)*state % rates(rate_idx,irgp17f) &
         -ymol(io16)*ymol(ih1)*state % rates(rate_idx,irpg16o) &
         -ymol(io15)*ymol(ihe4)*state % rates(rate_idx,irag15o)*state % rates(rate_idx,irr1) &
         -TWO*ymol(io16)*ymol(ihe4)*state % rates(rate_idx,irag16o) &
         -THREE*ymol(io15)*ymol(ihe4)*state % rates(rate_idx,irag15o)*(ONE-state % rates(rate_idx,irr1)) &
         -ymol(if17)*ymol(ih1)*state % rates(rate_idx,irpg17f)*(ONE-state % rates(rate_idx,irs1)) &
         -TWO*ymol(if17)*state % rates(rate_idx,irwk17f) &
         -EIGHT*ymol(img22)*state % rates(rate_idx,irlambda1)*(ONE-state % rates(3,delta1)) &
         -26.e0_dp_t*ymol(is30)*state % rates(rate_idx,irlambda2)*(ONE-state % rates(3,delta2))


    if (.not. doing_dratesdt) then
       dydt(ini56) = dydt(ini56)+ymol(ini56)*ymol(ih1)*state % rates(3,r56eff)
       dydt(ih1) = dydt(ih1)-56.0d0*ymol(ini56)*ymol(ih1)*state % rates(3,r56eff)
    else
       dydt(ini56) = dydt(ini56)+ymol(ini56)*ymol(ih1)*state % rates(3,dr56effdt)
       dydt(ih1) = dydt(ih1)-56.0d0*ymol(ini56)*ymol(ih1)*state % rates(3,dr56effdt)
    endif

  end subroutine make_ydots




  subroutine actual_jac(state)

    use bl_types
    use bl_constants_module
    use network

    implicit none

    type (burn_t)    :: state

    double precision :: ymol(nspec), t9
    double precision :: psum
    integer          :: i, j

    ! initialize
    state % jac(:,:) = ZERO
    ymol = state % xn / aion
    t9 = state % T / 10.0**9

    ! carbon-12
    state % jac(ic12,ic12) = -ymol(ih1)*state % rates(1,irpg12c)
    state % jac(ic12,io15) = state % rates(1,irlambCNO)
    state % jac(ic12,ihe4) = HALF*ymol(ihe4)*ymol(ihe4)*state % rates(1,ir3a)
    state % jac(ic12,ih1)  = -ymol(ic12)*state % rates(1,irpg12c) &
         +ymol(io15)*state % rates(3,dlambCNOdh1)

    ! oxygen-14
    state % jac(io14,ic12) = ymol(ih1)*state % rates(1,irpg12c)
    state % jac(io14,io14) = -ymol(ihe4)*state % rates(1,irap14o) &
         -state % rates(1,irwk14o)
    state % jac(io14,ihe4) = -ymol(io14)*state % rates(1,irap14o)
    state % jac(io14,ih1)  = ymol(ic12)*state % rates(1,irpg12c)

    ! oxygen-15
    state % jac(io15,io14) = state % rates(1,irwk14o)
    state % jac(io15,io15) = -ymol(ihe4)*state % rates(1,irag15o) &
         -state % rates(1,irlambCNO)
    state % jac(io15,if17) = ymol(ih1)*state % rates(1,irpg17f)*state % rates(1,irs1) &
         +state % rates(1,irwk17f)
    state % jac(io15,ihe4) = -ymol(io15)*state % rates(1,irag15o) &
         +ymol(if17)*ymol(ih1)*state % rates(1,irpg17f)*state % rates(3,drs1dhe4)
    state % jac(io15,ih1)  = ymol(if17)*state % rates(1,irpg17f)*state % rates(1,irs1) &
         -ymol(io15)*state % rates(3,dlambCNOdh1)

    ! oxygen-16
    state % jac(io16,io15) = ymol(ihe4)*state % rates(1,irr1)*state % rates(1,irag15o)
    state % jac(io16,io16) = -ymol(ih1)*state % rates(1,irpg16o) &
         -ymol(ihe4)*state % rates(1,irag16o)
    state % jac(io16,if17) = state % rates(1,irgp17f)
    state % jac(io16,ihe4) = ymol(io15)*state % rates(1,irr1)*state % rates(1,irag15o) &
         -ymol(io16)*state % rates(1,irag16o)
    state % jac(io16,ih1)  = -ymol(io16)*state % rates(1,irpg16o) &
         +ymol(io15)*ymol(ihe4)*state % rates(3,drr1dh1)*state % rates(1,irag15o)

    ! flourine-17
    state % jac(if17,io14) = ymol(ihe4)*state % rates(1,irap14o)
    state % jac(if17,io16) = ymol(ih1)*state % rates(1,irpg16o)
    state % jac(if17,if17) = -state % rates(1,irgp17f) &
         -ymol(ih1)*state % rates(1,irpg17f) &
         -state % rates(1,irwk17f)
    state % jac(if17,ihe4) = ymol(io14)*state % rates(1,irap14o)
    state % jac(if17,ih1)  = ymol(io16)*state % rates(1,irpg16o) &
         -ymol(if17)*state % rates(1,irpg17f)

    ! magnesium-22
    state % jac(img22,io15) = ymol(ihe4)*state % rates(1,irag15o)*(ONE-state % rates(1,irr1))
    state % jac(img22,io16) = ymol(ihe4)*state % rates(1,irag16o)
    state % jac(img22,if17) = ymol(ih1)*state % rates(1,irpg17f)*(ONE-state % rates(1,irs1))
    state % jac(img22,img22) = -state % rates(1,irlambda1)
    state % jac(img22,ihe4) = ymol(io16)*state % rates(1,irag16o) &
         +ymol(io15)*state % rates(1,irag15o)*(ONE-state % rates(1,irr1)) &
         -ymol(if17)*ymol(ih1)*state % rates(1,irpg17f)*state % rates(3,drs1dhe4) &
         -ymol(img22)*state % rates(3,dlambda1dhe4)
    state % jac(img22,ih1)  = ymol(if17)*state % rates(1,irpg17f)*(ONE-state % rates(1,irs1)) &
         -ymol(io15)*ymol(ihe4)*state % rates(1,irag15o)*state % rates(3,drr1dh1)

    ! sulfur-30
    state % jac(is30,img22) = state % rates(1,irlambda1)
    state % jac(is30,is30)  = -state % rates(1,irlambda2)
    state % jac(is30,ihe4)  = ymol(img22)*state % rates(3,dlambda1dhe4) &
         -ymol(is30)*state % rates(3,dlambda2dhe4)

    ! nickel-56
    state % jac(ini56,is30) = state % rates(1,irlambda2)
    state % jac(ini56,ini56) = ymol(ih1)*state % rates(3,r56eff)
    state % jac(ini56,ihe4) = ymol(is30)*state % rates(3,dlambda2dhe4)
    state % jac(ini56,ih1) = ymol(ini56)*state % rates(3,r56eff)

    ! helium-4
    state % jac(ihe4,io14) = -ymol(ihe4)*state % rates(1,irap14o)
    state % jac(ihe4,io15) = state % rates(1,irlambCNO) &
         -ymol(ihe4)*state % rates(1,irag15o)*(ONE-state % rates(1,irr1))
    state % jac(ihe4,io16) = -ymol(ihe4)*state % rates(1,irag16o)
    state % jac(ihe4,if17) = ymol(ih1)*state % rates(1,irpg17f)*state % rates(1,irs1) &
         -ymol(ih1)*state % rates(1,irpg17f)*(ONE-state % rates(1,irs1)) &
         +state % rates(1,irwk17f)
    state % jac(ihe4,img22) = -TWO*state % rates(1,irlambda1)*state % rates(3,delta1)
    state % jac(ihe4,is30) = -6.5d0*state % rates(1,irlambda2)*state % rates(3,delta2)
    state % jac(ihe4,ihe4) = -THREE*ymol(ihe4)*ymol(ihe4)*HALF*state % rates(1,ir3a) &
         -ymol(io14)*state % rates(1,irap14o) &
         -ymol(io16)*state % rates(1,irag16o) &
         -ymol(io15)*state % rates(1,irag15o)*(ONE-state % rates(1,irr1)) &
         +ymol(if17)*ymol(ih1)*state % rates(1,irpg17f)*state % rates(3,drs1dhe4) &
         +ymol(if17)*ymol(ih1)*state % rates(1,irpg17f)*state % rates(3,drs1dhe4) &
         -TWO*ymol(img22)*state % rates(3,dlambda1dhe4)*state % rates(3,delta1) &
         -6.5d0*ymol(is30)*state % rates(3,dlambda2dhe4)*state % rates(3,delta2)
    state % jac(ihe4,ih1)  = ymol(if17)*state % rates(1,irpg17f)*state % rates(1,irs1) &
         -ymol(if17)*state % rates(1,irpg17f)*(ONE-state % rates(1,irs1)) &
         +ymol(io15)*state % rates(3,dlambCNOdh1) &
         +ymol(io15)*ymol(ihe4)*state % rates(1,irag15o)*state % rates(3,drr1dh1)

    ! hydrogen-1
    state % jac(ih1,ic12) = -TWO*ymol(ih1)*state % rates(1,irpg12c)
    state % jac(ih1,io14) = ymol(ihe4)*state % rates(1,irap14o) &
         -state % rates(1,irwk14o)
    state % jac(ih1,io15) = -state % rates(1,irlambCNO) &
         -ymol(ihe4)*state % rates(1,irag15o)*state % rates(1,irr1) &
         -THREE*ymol(ihe4)*state % rates(1,irag15o)*(ONE-state % rates(1,irr1))
    state % jac(ih1,io16) = -ymol(ih1)*state % rates(1,irpg16o) &
         -TWO*ymol(ihe4)*state % rates(1,irag16o)
    state % jac(ih1,if17) = -TWO*ymol(ih1)*state % rates(1,irpg17f)*state % rates(1,irs1) &
         +state % rates(1,irgp17f) &
         -ymol(ih1)*state % rates(1,irpg17f)*(ONE-state % rates(1,irs1)) &
         -TWO*state % rates(1,irwk17f)
    state % jac(ih1,img22) = -EIGHT*state % rates(1,irlambda1)*(ONE-state % rates(3,delta1))
    state % jac(ih1,is30)  = -26.d0*state % rates(1,irlambda2)*(ONE-state % rates(3,delta2))
    state % jac(ih1,ini56) = -56.0d0*ymol(ih1)*state % rates(3,r56eff)
    state % jac(ih1,ihe4) = ymol(io14)*state % rates(1,irap14o) &
         -ymol(io15)*state % rates(1,irag15o)*state % rates(1,irr1) &
         -TWO*ymol(io16)*state % rates(1,irag16o) &
         -THREE*ymol(io15)*state % rates(1,irag15o)*(ONE-state % rates(1,irr1)) &
         -ymol(if17)*ymol(ih1)*state % rates(1,irpg17f)*state % rates(3,drs1dhe4) &
         -EIGHT*ymol(img22)*state % rates(3,dlambda1dhe4)*(ONE-state % rates(3,delta1)) &
         -26.d0*ymol(is30)*state % rates(3,dlambda2dhe4)*(ONE-state % rates(3,delta2))
    state % jac(ih1,ih1)  = -TWO*ymol(ic12)*state % rates(1,irpg12c) &
         -TWO*ymol(if17)*state % rates(1,irpg17f)*state % rates(1,irs1) &
         -ymol(io16)*state % rates(1,irpg16o) &
         -ymol(if17)*state % rates(1,irpg17f)*(ONE-state % rates(1,irs1)) &
         -ymol(io15)*state % rates(3,dlambCNOdh1) &
         +TWO*ymol(io15)*ymol(ihe4)*state % rates(1,irag15o)*state % rates(3,drr1dh1) &
         -56.0d0*ymol(ini56)*state % rates(3,r56eff)

    ! temperature derivatives df(Y)/df(T)
    call make_ydots(ymol,t9,state,state % jac(1:nspec,net_itemp),.true.)

    ! Temperature Jacobian elements
    call temperature_jac(state)

  end subroutine actual_jac

end module actual_rhs_module
