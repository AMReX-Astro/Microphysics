! The f_rhs routine provides the right-hand-side for the DVODE solver.
! It deals with molar abundances throughout (we expect that the input
! vector y has molar abundances, Y = X/A) for make_rates, and
! dydt routines.  It also checks to see if the temperature has changed
! much since the last call - if so, it updates the temperature to get
! a better estimate of the reaction rates.
!
! The jac routine provides an explicit Jacobian to the DVODE solver.
!

subroutine f_rhs(n, time, y, ydot, rpar, ipar)

  use bl_types
  use bl_constants_module
  use network
  use eos_module
  use eos_type_module
  use rates_module
  use rpar_indices

  implicit none

  integer,          intent(IN   ) :: n, ipar
  double precision, intent(IN   ) :: time, y(n)
  double precision, intent(INOUT) :: rpar(n_rpar_comps)
  double precision, intent(  OUT) :: ydot(n)

  integer :: k

  double precision, parameter :: T2T9 = 1.0d-9

  double precision :: dens, cp, dhdx(nspec), T9_eos, dT_crit, t9
  double precision :: ymol(nspec)

  type (eos_t) :: eos_state
  
  ydot = ZERO

  ! several thermodynamic quantities come in via rpar
  dens = rpar(irp_dens)
  T9_eos = rpar(irp_T9_eos)
  dT_crit = rpar(irp_dTcrit)

  ymol = y(1:nspec)
  t9 = y(nspec+1)

  if (abs(t9 - T9_eos) > dT_crit*T9_eos) then
     T9_eos = t9

     eos_state%T = T9_eos/T2T9
     eos_state%rho = dens
     eos_state%xn = ymol*aion

     call eos(eos_input_rt, eos_state)
     
     rpar(irp_T9_eos) = T9_eos
     rpar(irp_cp) = eos_state%cp
     rpar(irp_dhdX:irp_dhdX+nspec-1) = eos_state%dhdX
  endif

  ! more thermodyanmics that possibly were updated
  cp   = rpar(irp_cp)
  dhdX = rpar(irp_dhdX:irp_dhdX+nspec-1)

  ! build the rates; weak rates are the wk* variables
  call make_rates(t9, dens, ymol, rpar)

  ! set up the ODEs
  call make_ydots(ymol,t9,rpar,ydot(1:nspec),.false.)

!.... t9
! shouldn't an aion be here?  putting it in breaks VODE convergence...
  ydot(n) = -sum((dhdX+ebin)*ydot(1:nspec))/cp
  ydot(n) = ydot(n) * T2T9

end subroutine f_rhs




subroutine make_rates(t9,dens,y,rpar)

  use bl_types
  use bl_constants_module
  use rates_module
  use network
  use rpar_indices
  
  implicit none

  double precision, intent(in   ) :: t9, dens, y(nspec)
  double precision, intent(inout) :: rpar(n_rpar_comps)

  ! locally used rates
  double precision :: rate,dratedt,wk18ne,wk19ne
  double precision :: r56pg,dr56pgdt,cutoni,dcutonidt,r57decay,r56eff,dr56effdt
  double precision :: t9i32

  ! some numbers from appendix C in WW81; these should probably be
  ! updated with current rates
  double precision, parameter :: Lweak = 1.05d0, & ! this is for NS
!                                Lweak = 0.107d0, & ! this is for lower 
                                                    ! densities
                                 la2 = ONE/FIFTEEN ! mean rate from 30s to 56ni
                                                   ! from p-capture and beta 
                                                   ! decays
                                
  type (temp_t) :: tfactors
  
  rpar(irp_rates:irp_drtdt-1+nrat) = ZERO ! rates and dratesdt
  rpar(irp_dlambCNOdh1:n_rpar_comps) = ZERO ! other constants

  tfactors = calc_tfactors(t9)

  ! some common parameters
  rpar(irp_rates-1+irLweak) = Lweak
  rpar(irp_rates-1+irla2)   = la2

  ! weak rates first
  !
  ! 14o(beta nu)14n
  call rate_o14_to_n14(tfactors,rate,dratedt)
  rpar(irp_rates-1+irwk14o) = rate
  ! 15o(beta nu)15n
  call rate_o15_to_n15(tfactors,rate,dratedt)
  rpar(irp_rates-1+irwk15o) = rate
  ! 17f(beta nu)17o
  call rate_f17_to_o17(tfactors,rate,dratedt)
  rpar(irp_rates-1+irwk17f) = rate
  ! these weak rates aren't needed outside of this routine
  ! 18ne(beta nu)18f
  call rate_ne18_to_f18(tfactors,wk18ne,dratedt) 
  ! 19ne(beta nu)19f
  call rate_ne19_to_f19(tfactors,wk19ne,dratedt)

  ! 12c(p,g)13n
  call rate_p_c12_to_n13(tfactors,rate,dratedt)
  rpar(irp_rates-1+irpg12c) = dens*rate
  rpar(irp_drtdt-1+irpg12c) = dens*dratedt

  ! triple alpha
  call rate_he4_he4_he4_to_c12(tfactors,rate,dratedt)
  rpar(irp_rates-1+ir3a) = dens*dens*rate
  rpar(irp_drtdt-1+ir3a) = dens*dens*dratedt

  ! 17f(p,g)18ne
  call rate_p_f17_to_ne18(tfactors,rate,dratedt)
  rpar(irp_rates-1+irpg17f) = dens*rate
  rpar(irp_drtdt-1+irpg17f) = dens*dratedt

  ! 17f(g,p)16o
  call rate_f17_to_p_o16(tfactors,rate,dratedt)
  rpar(irp_rates-1+irgp17f) = rate
  rpar(irp_drtdt-1+irgp17f) = dratedt

  ! 15o(a,g)19ne
  call rate_he4_o15_to_ne19(tfactors,rate,dratedt)
  rpar(irp_rates-1+irag15o) = dens*rate
  rpar(irp_drtdt-1+irag15o) = dens*dratedt

  ! 16o(a,g)20ne
  call rate_he4_o16_to_ne20(tfactors,rate,dratedt)
  rpar(irp_rates-1+irag16o) = dens*rate
  rpar(irp_drtdt-1+irag16o) = dens*dratedt

  ! 16o(p,g)17f
  call rate_p_o16_to_f17(tfactors,rate,dratedt)
  rpar(irp_rates-1+irpg16o) = dens*rate
  rpar(irp_drtdt-1+irpg16o) = dens*dratedt

  ! 14o(a,p)17f
  call rate_he4_o14_to_p_f17(tfactors,rate,dratedt)
  rpar(irp_rates-1+irap14o) = dens*rate
  rpar(irp_drtdt-1+irap14o) = dens*dratedt

  ! limit CNO as minimum between 14n(p,g)15o and 15o(beta nu)15n
  ! we store the limited rate in irlambCNO; this is lambda_CNO in WW81
  call rate_p_n14_to_o15(tfactors,rate,dratedt)
  rpar(irp_rates-1+irlambCNO) = min(rpar(irp_rates-1+irwk15o),rate*dens*y(ih1))
  if (rpar(irp_rates-1+irlambCNO) < rpar(irp_rates-1+irwk15o)) then
     rpar(irp_drtdt-1+irlambCNO) = dens*y(ih1)*dratedt
     rpar(irp_dlambCNOdh1) = rate*dens
  endif

  ! 22mg(...)30s
  ! check if this proceeds via p-captures or (a,p) reactions
  ! the Lweak is from WW81, eqn C15
  ! we store the rate in irlambda1; this is the lambda1 in WW81
  call rate_he4_si26_to_p_p29(tfactors,rate,dratedt)
  rpar(irp_rates-1+irlambda1) = max(rpar(irp_rates-1+irLweak),dens*y(ihe4)*rate)
  if (rpar(irp_rates-1+irlambda1) > rpar(irp_rates-1+irLweak)) then
       rpar(irp_drtdt-1+irlambda1) = dens*y(ihe4)*dratedt
       rpar(irp_dlambda1dhe4) = dens*rate
  ! use the sign of rpar(irp_rates-1+irlambda1) to indicate the value of delta1 in WW81
  ! if delta1 = 1, then we multiply the rate by -1
       rpar(irp_rates-1+irlambda1) = -ONE*rpar(irp_rates-1+irlambda1)
  endif

  ! 30s(...) 56ni 
  ! check if this proceeds via p-captures or (a,p) reactions
  ! use 44ti(a,p)v47 as a typical limiting rate for the (a,p) process
  ! store this in irlambda2; this is lambda2 in WW81
  call rate_he4_ti44_to_p_v47(tfactors,rate,dratedt)
  rpar(irp_rates-1+irlambda2) = max(rpar(irp_rates-1+irla2),dens*y(ihe4)*rate)
  if (rpar(irp_rates-1+irlambda2) > rpar(irp_rates-1+irla2)) then
       rpar(irp_drtdt-1+irlambda2) = dens*y(ihe4)*dratedt
       rpar(irp_dlambda2dhe4) = dens*rate
  ! use the sign of rpar(irp_rates-1+irlambda2) to indicate th value of delta2
  ! if delta2 = 1, then we multiply the rate by -1
       rpar(irp_rates-1+irlambda2) = -ONE*rpar(irp_rates-1+irlambda2)
  endif

  ! form s1 from WW81; branching ratio for 18ne beta decay (wk18ne) vs (a,p)
  ! store result in irs1
  ! 18ne(a,p)21na
  call rate_he4_ne18_to_p_na21(tfactors,rate,dratedt)
  rpar(irp_rates-1+irs1) = wk18ne / (wk18ne + dens*y(ihe4)*rate)
  rpar(irp_drtdt-1+irs1) = -rpar(irp_rates-1+irs1)*dens*y(ihe4)*dratedt &
       / (wk18ne + dens*y(ihe4)*rate)
  rpar(irp_drs1dhe4) = -rpar(irp_rates-1+irs1)*dens*rate &
       / (wk18ne + dens*y(ihe4)*rate)

  ! form r1 from WW81; ranching ratio for 19ne beta decay (wk19ne) vs (p,g)
  ! store result in irr1
  ! 19ne(p,g)20na
  call rate_p_ne19_to_na20(tfactors,rate,dratedt)
  rpar(irp_rates-1+irr1) = wk19ne / (wk19ne + dens*y(ih1)*rate)
  rpar(irp_drtdt-1+irr1) = -rpar(irp_rates-1+irr1)*dens*y(ih1)*dratedt &
       / (wk19ne + dens*y(ih1)*rate)
  rpar(irp_drr1dh1) = -rpar(irp_rates-1+irr1)*dens*rate &
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
  r56eff=min(r56pg,cutoni*r57decay)
!   rpar(irp_r56eff) = r56eff
!   if (r56eff < r56pg) rpar(irp_dr56effdt) = r57decay*dcutonidt
  rpar(irp_r56eff) = ZERO
  rpar(irp_dr56effdt) = ZERO

end subroutine make_rates




subroutine make_ydots(ymol,t9,rpar,dydt,doing_dratesdt)

  use bl_types
  use bl_constants_module
  use network
  use rpar_indices
  
  implicit none

  double precision, intent(IN   ) :: ymol(nspec), t9
  logical ,         intent(IN   ) :: doing_dratesdt
  double precision, intent(INOUT) :: rpar(n_rpar_comps)
  double precision, intent(  OUT) :: dydt(nspec)
  
  integer :: irp_start
  double precision :: delta1, delta2
  double precision :: dens

  ! initialize
  dydt = ZERO
  dens = rpar(irp_dens)

  ! check to see if we are doing this with the t-derivatives
  ! if so, offset our starting index in rpar
  irp_start = irp_rates
  if (doing_dratesdt) irp_start = irp_drtdt

  if (.not. doing_dratesdt) then
     delta1 = ZERO; delta2 = ZERO
     ! figure out the delta's; we used negative rates to indicate delta=1
     if (rpar(irp_rates-1+irlambda1) < ZERO) then
        delta1 = ONE
        rpar(irp_rates-1+irlambda1) = -ONE*rpar(irp_rates-1+irlambda1)
     endif
     if (rpar(irp_rates-1+irlambda2) < ZERO) then
        delta2 = ONE
        rpar(irp_rates-1+irlambda2) = -ONE*rpar(irp_rates-1+irlambda2)
     endif
     rpar(irp_delta1) = delta1
     rpar(irp_delta2) = delta2
  endif

! setup ODEs
!
!....
!.... 12c = 1
!....
      dydt(ic12)=-ymol(ic12)*ymol(ih1)*rpar(irp_start-1+irpg12c) &
           +ymol(ihe4)**3*rpar(irp_start-1+ir3a)/SIX &
           +ymol(io15)*rpar(irp_start-1+irlambCNO)
!....
!.... 14o = 2
!....
      dydt(io14)=-ymol(io14)*ymol(ihe4)*rpar(irp_start-1+irap14o) &
           -ymol(io14)*rpar(irp_start-1+irwk14o) &
           +ymol(ic12)*ymol(ih1)*rpar(irp_start-1+irpg12c)
!....
!.... 15o = 3
!....
      dydt(io15)=ymol(io14)*rpar(irp_start-1+irwk14o) &
           -ymol(io15)*ymol(ihe4)*rpar(irp_start-1+irag15o) &
           -ymol(io15)*rpar(irp_start-1+irlambCNO) &
           +ymol(if17)*ymol(ih1)*rpar(irp_start-1+irpg17f)*rpar(irp_start-1+irs1) &
           +ymol(if17)*rpar(irp_start-1+irwk17f)
!....
!.... 16o = 4
!....
      dydt(io16) = ymol(if17)*rpar(irp_start-1+irgp17f) &
           -ymol(io16)*ymol(ih1)*rpar(irp_start-1+irpg16o) &
           +ymol(io15)*ymol(ihe4)*rpar(irp_start-1+irr1)*rpar(irp_start-1+irag15o) &
           -ymol(io16)*ymol(ihe4)*rpar(irp_start-1+irag16o)
!....
!.... 17f = 5
!....
      dydt(if17)=ymol(io14)*ymol(ihe4)*rpar(irp_start-1+irap14o) &
           +ymol(io16)*ymol(ih1)*rpar(irp_start-1+irpg16o) &
           -ymol(if17)*rpar(irp_start-1+irgp17f) &
           -ymol(if17)*ymol(ih1)*rpar(irp_start-1+irpg17f) &
           -ymol(if17)*rpar(irp_start-1+irwk17f)
!....
!.... 22mg = 6
!....
      dydt(img22)=ymol(io16)*ymol(ihe4)*rpar(irp_start-1+irag16o) &
           +ymol(if17)*ymol(ih1)*rpar(irp_start-1+irpg17f)*(ONE-rpar(irp_start-1+irs1)) &
           +ymol(io15)*ymol(ihe4)*rpar(irp_start-1+irag15o)*(ONE-rpar(irp_start-1+irr1)) &
           -ymol(img22)*rpar(irp_start-1+irlambda1)
!....
!.... 30s = 7
!....
      dydt(is30)=ymol(img22)*rpar(irp_start-1+irlambda1) &
           -ymol(is30)*rpar(irp_start-1+irlambda2)
!....
!.... amax (56ni) = 8  (note that WW81 have a typo -- they write lambda1 here)
!....
      dydt(ini56)=ymol(is30)*rpar(irp_start-1+irlambda2)
!....
!.... 4he (alpha) = 9
!....
      dydt(ihe4)=-ymol(ihe4)**3*HALF*rpar(irp_start-1+ir3a) &
           +ymol(io15)*rpar(irp_start-1+irlambCNO) &
           -ymol(io14)*ymol(ihe4)*rpar(irp_start-1+irap14o) &
           +ymol(if17)*ymol(ih1)*rpar(irp_start-1+irpg17f)*rpar(irp_start-1+irs1) &
           -ymol(io15)*ymol(ihe4)*rpar(irp_start-1+irag15o)*(ONE-rpar(irp_start-1+irr1)) &
           -ymol(io16)*ymol(ihe4)*rpar(irp_start-1+irag16o) &
           -ymol(if17)*ymol(ih1)*rpar(irp_start-1+irpg17f)*(ONE-rpar(irp_start-1+irs1)) &
           +ymol(if17)*rpar(irp_start-1+irwk17f) &
           -TWO*ymol(img22)*rpar(irp_start-1+irlambda1)*rpar(irp_delta1) &
           -6.5d0*ymol(is30)*rpar(irp_start-1+irlambda2)*rpar(irp_delta2)
!....
!.... 1h (p) = 10
!....
      dydt(ih1)=-ymol(io14)*rpar(irp_start-1+irwk14o) &
           -ymol(io15)*rpar(irp_start-1+irlambCNO) &
           -TWO*ymol(ic12)*ymol(ih1)*rpar(irp_start-1+irpg12c) &
           +ymol(io14)*ymol(ihe4)*rpar(irp_start-1+irap14o) &
           -TWO*ymol(if17)*ymol(ih1)*rpar(irp_start-1+irpg17f)*rpar(irp_start-1+irs1) &
           +ymol(if17)*rpar(irp_start-1+irgp17f) &
           -ymol(io16)*ymol(ih1)*rpar(irp_start-1+irpg16o) &
           -ymol(io15)*ymol(ihe4)*rpar(irp_start-1+irag15o)*rpar(irp_start-1+irr1) &
           -TWO*ymol(io16)*ymol(ihe4)*rpar(irp_start-1+irag16o) &
           -THREE*ymol(io15)*ymol(ihe4)*rpar(irp_start-1+irag15o)*(ONE-rpar(irp_start-1+irr1)) &
           -ymol(if17)*ymol(ih1)*rpar(irp_start-1+irpg17f)*(ONE-rpar(irp_start-1+irs1)) &
           -TWO*ymol(if17)*rpar(irp_start-1+irwk17f) &
           -EIGHT*ymol(img22)*rpar(irp_start-1+irlambda1)*(ONE-rpar(irp_delta1)) &
           -26.e0_dp_t*ymol(is30)*rpar(irp_start-1+irlambda2)*(ONE-rpar(irp_delta2))


      if (.not. doing_dratesdt) then
         dydt(ini56) = dydt(ini56)+ymol(ini56)*ymol(ih1)*rpar(irp_r56eff)
         dydt(ih1) = dydt(ih1)-56.0d0*ymol(ini56)*ymol(ih1)*rpar(irp_r56eff)
      else
         dydt(ini56) = dydt(ini56)+ymol(ini56)*ymol(ih1)*rpar(irp_dr56effdt)
         dydt(ih1) = dydt(ih1)-56.0d0*ymol(ini56)*ymol(ih1)*rpar(irp_dr56effdt)
      endif

end subroutine make_ydots




subroutine jac(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)

  use bl_types
  use bl_constants_module
  use network
  use rpar_indices

  implicit none

  integer        , intent(IN   ) :: neq, ml, mu, nrpd, ipar
  double precision, intent(IN   ) :: y(neq), t
  double precision, intent(INOUT) :: rpar(*)
  double precision, intent(  OUT) :: pd(neq,neq)

  double precision, parameter :: T2T9 = 1.0e-9_dp_t
  double precision :: ymol(nspec), t9
  double precision :: cp, dhdX(nspec)
  double precision :: psum
  integer :: it9
  integer :: i, j


  ! initialize
  pd(:,:)  = ZERO
  ymol = y(1:nspec)
  it9 = neq
  t9 = y(it9)
  cp = rpar(irp_cp)
  dhdX = rpar(irp_dhdX:irp_dhdX+nspec-1)

  ! carbon-12
  pd(ic12,ic12) = -ymol(ih1)*rpar(irp_rates-1+irpg12c)
  pd(ic12,io15) = rpar(irp_rates-1+irlambCNO)
  pd(ic12,ihe4) = HALF*ymol(ihe4)*ymol(ihe4)*rpar(irp_rates-1+ir3a)
  pd(ic12,ih1)  = -ymol(ic12)*rpar(irp_rates-1+irpg12c) &
       +ymol(io15)*rpar(irp_dlambCNOdh1)

  ! oxygen-14
  pd(io14,ic12) = ymol(ih1)*rpar(irp_rates-1+irpg12c)
  pd(io14,io14) = -ymol(ihe4)*rpar(irp_rates-1+irap14o) &
       -rpar(irp_rates-1+irwk14o)
  pd(io14,ihe4) = -ymol(io14)*rpar(irp_rates-1+irap14o)
  pd(io14,ih1)  = ymol(ic12)*rpar(irp_rates-1+irpg12c)

  ! oxygen-15
  pd(io15,io14) = rpar(irp_rates-1+irwk14o)
  pd(io15,io15) = -ymol(ihe4)*rpar(irp_rates-1+irag15o) &
       -rpar(irp_rates-1+irlambCNO)
  pd(io15,if17) = ymol(ih1)*rpar(irp_rates-1+irpg17f)*rpar(irp_rates-1+irs1) &
       +rpar(irp_rates-1+irwk17f)
  pd(io15,ihe4) = -ymol(io15)*rpar(irp_rates-1+irag15o) &
       +ymol(if17)*ymol(ih1)*rpar(irp_rates-1+irpg17f)*rpar(irp_drs1dhe4)
  pd(io15,ih1)  = ymol(if17)*rpar(irp_rates-1+irpg17f)*rpar(irp_rates-1+irs1) &
       -ymol(io15)*rpar(irp_dlambCNOdh1)

  ! oxygen-16
  pd(io16,io15) = ymol(ihe4)*rpar(irp_rates-1+irr1)*rpar(irp_rates-1+irag15o)
  pd(io16,io16) = -ymol(ih1)*rpar(irp_rates-1+irpg16o) &
       -ymol(ihe4)*rpar(irp_rates-1+irag16o)
  pd(io16,if17) = rpar(irp_rates-1+irgp17f)
  pd(io16,ihe4) = ymol(io15)*rpar(irp_rates-1+irr1)*rpar(irp_rates-1+irag15o) &
       -ymol(io16)*rpar(irp_rates-1+irag16o)
  pd(io16,ih1)  = -ymol(io16)*rpar(irp_rates-1+irpg16o) &
       +ymol(io15)*ymol(ihe4)*rpar(irp_drr1dh1)*rpar(irp_rates-1+irag15o)

  ! flourine-17
  pd(if17,io14) = ymol(ihe4)*rpar(irp_rates-1+irap14o)
  pd(if17,io16) = ymol(ih1)*rpar(irp_rates-1+irpg16o)
  pd(if17,if17) = -rpar(irp_rates-1+irgp17f) &
       -ymol(ih1)*rpar(irp_rates-1+irpg17f) &
       -rpar(irp_rates-1+irwk17f)
  pd(if17,ihe4) = ymol(io14)*rpar(irp_rates-1+irap14o)
  pd(if17,ih1)  = ymol(io16)*rpar(irp_rates-1+irpg16o) &
       -ymol(if17)*rpar(irp_rates-1+irpg17f)

  ! magnesium-22
  pd(img22,io15) = ymol(ihe4)*rpar(irp_rates-1+irag15o)*(ONE-rpar(irp_rates-1+irr1))
  pd(img22,io16) = ymol(ihe4)*rpar(irp_rates-1+irag16o)
  pd(img22,if17) = ymol(ih1)*rpar(irp_rates-1+irpg17f)*(ONE-rpar(irp_rates-1+irs1))
  pd(img22,img22) = -rpar(irp_rates-1+irlambda1)
  pd(img22,ihe4) = ymol(io16)*rpar(irp_rates-1+irag16o) &
       +ymol(io15)*rpar(irp_rates-1+irag15o)*(ONE-rpar(irp_rates-1+irr1)) &
       -ymol(if17)*ymol(ih1)*rpar(irp_rates-1+irpg17f)*rpar(irp_drs1dhe4) &
       -ymol(img22)*rpar(irp_dlambda1dhe4)
  pd(img22,ih1)  = ymol(if17)*rpar(irp_rates-1+irpg17f)*(ONE-rpar(irp_rates-1+irs1)) &
       -ymol(io15)*ymol(ihe4)*rpar(irp_rates-1+irag15o)*rpar(irp_drr1dh1)

  ! sulfur-30
  pd(is30,img22) = rpar(irp_rates-1+irlambda1)
  pd(is30,is30)  = -rpar(irp_rates-1+irlambda2)
  pd(is30,ihe4)  = ymol(img22)*rpar(irp_dlambda1dhe4) &
       -ymol(is30)*rpar(irp_dlambda2dhe4)

  ! nickel-56
  pd(ini56,is30) = rpar(irp_rates-1+irlambda2)
  pd(ini56,ini56) = ymol(ih1)*rpar(irp_r56eff)
  pd(ini56,ihe4) = ymol(is30)*rpar(irp_dlambda2dhe4)
  pd(ini56,ih1) = ymol(ini56)*rpar(irp_r56eff)

  ! helium-4
  pd(ihe4,io14) = -ymol(ihe4)*rpar(irp_rates-1+irap14o)
  pd(ihe4,io15) = rpar(irp_rates-1+irlambCNO) &
       -ymol(ihe4)*rpar(irp_rates-1+irag15o)*(ONE-rpar(irp_rates-1+irr1))
  pd(ihe4,io16) = -ymol(ihe4)*rpar(irp_rates-1+irag16o)
  pd(ihe4,if17) = ymol(ih1)*rpar(irp_rates-1+irpg17f)*rpar(irp_rates-1+irs1) &
       -ymol(ih1)*rpar(irp_rates-1+irpg17f)*(ONE-rpar(irp_rates-1+irs1)) &
       +rpar(irp_rates-1+irwk17f)
  pd(ihe4,img22) = -TWO*rpar(irp_rates-1+irlambda1)*rpar(irp_delta1)
  pd(ihe4,is30) = -6.5d0*rpar(irp_rates-1+irlambda2)*rpar(irp_delta2)
  pd(ihe4,ihe4) = -THREE*ymol(ihe4)*ymol(ihe4)*HALF*rpar(irp_rates-1+ir3a) &
       -ymol(io14)*rpar(irp_rates-1+irap14o) &
       -ymol(io16)*rpar(irp_rates-1+irag16o) &
       -ymol(io15)*rpar(irp_rates-1+irag15o)*(ONE-rpar(irp_rates-1+irr1)) &
       +ymol(if17)*ymol(ih1)*rpar(irp_rates-1+irpg17f)*rpar(irp_drs1dhe4) &
       +ymol(if17)*ymol(ih1)*rpar(irp_rates-1+irpg17f)*rpar(irp_drs1dhe4) &
       -TWO*ymol(img22)*rpar(irp_dlambda1dhe4)*rpar(irp_delta1) &
       -6.5d0*ymol(is30)*rpar(irp_dlambda2dhe4)*rpar(irp_delta2)
  pd(ihe4,ih1)  = ymol(if17)*rpar(irp_rates-1+irpg17f)*rpar(irp_rates-1+irs1) &
       -ymol(if17)*rpar(irp_rates-1+irpg17f)*(ONE-rpar(irp_rates-1+irs1)) &
       +ymol(io15)*rpar(irp_dlambCNOdh1) &
       +ymol(io15)*ymol(ihe4)*rpar(irp_rates-1+irag15o)*rpar(irp_drr1dh1)

  ! hydrogen-1
  pd(ih1,ic12) = -TWO*ymol(ih1)*rpar(irp_rates-1+irpg12c)
  pd(ih1,io14) = ymol(ihe4)*rpar(irp_rates-1+irap14o) &
       -rpar(irp_rates-1+irwk14o)
  pd(ih1,io15) = -rpar(irp_rates-1+irlambCNO) &
       -ymol(ihe4)*rpar(irp_rates-1+irag15o)*rpar(irp_rates-1+irr1) &
       -THREE*ymol(ihe4)*rpar(irp_rates-1+irag15o)*(ONE-rpar(irp_rates-1+irr1))
  pd(ih1,io16) = -ymol(ih1)*rpar(irp_rates-1+irpg16o) &
       -TWO*ymol(ihe4)*rpar(irp_rates-1+irag16o)
  pd(ih1,if17) = -TWO*ymol(ih1)*rpar(irp_rates-1+irpg17f)*rpar(irp_rates-1+irs1) &
       +rpar(irp_rates-1+irgp17f) &
       -ymol(ih1)*rpar(irp_rates-1+irpg17f)*(ONE-rpar(irp_rates-1+irs1)) &
       -TWO*rpar(irp_rates-1+irwk17f)
  pd(ih1,img22) = -EIGHT*rpar(irp_rates-1+irlambda1)*(ONE-rpar(irp_delta1))
  pd(ih1,is30)  = -26.d0*rpar(irp_rates-1+irlambda2)*(ONE-rpar(irp_delta2))
  pd(ih1,ini56) = -56.0d0*ymol(ih1)*rpar(irp_r56eff)
  pd(ih1,ihe4) = ymol(io14)*rpar(irp_rates-1+irap14o) &
       -ymol(io15)*rpar(irp_rates-1+irag15o)*rpar(irp_rates-1+irr1) &
       -TWO*ymol(io16)*rpar(irp_rates-1+irag16o) &
       -THREE*ymol(io15)*rpar(irp_rates-1+irag15o)*(ONE-rpar(irp_rates-1+irr1)) &
       -ymol(if17)*ymol(ih1)*rpar(irp_rates-1+irpg17f)*rpar(irp_drs1dhe4) &
       -EIGHT*ymol(img22)*rpar(irp_dlambda1dhe4)*(ONE-rpar(irp_delta1)) &
       -26.d0*ymol(is30)*rpar(irp_dlambda2dhe4)*(ONE-rpar(irp_delta2))
  pd(ih1,ih1)  = -TWO*ymol(ic12)*rpar(irp_rates-1+irpg12c) &
       -TWO*ymol(if17)*rpar(irp_rates-1+irpg17f)*rpar(irp_rates-1+irs1) &
       -ymol(io16)*rpar(irp_rates-1+irpg16o) &
       -ymol(if17)*rpar(irp_rates-1+irpg17f)*(ONE-rpar(irp_rates-1+irs1)) &
       -ymol(io15)*rpar(irp_dlambCNOdh1) &
       +TWO*ymol(io15)*ymol(ihe4)*rpar(irp_rates-1+irag15o)*rpar(irp_drr1dh1) &
       -56.0d0*ymol(ini56)*rpar(irp_r56eff)

  ! temperature derivatives df(Y)/df(T)
  call make_ydots(ymol,t9,rpar,pd(1:nspec,it9),.true.)

  ! temperature jacobian elements
  do i = 1, neq
     psum = 0.0d0
     do j = 1, nspec
        psum = psum + (dhdX(j) + ebin(j))*pd(j,i)
     enddo
     pd(it9,i) = -psum
  enddo

  pd(it9,:) = pd(it9,:) * T2T9 / cp

end subroutine jac
