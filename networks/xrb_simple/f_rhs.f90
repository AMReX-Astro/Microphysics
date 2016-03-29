! This is the right-hand-side function for the DVODE solver.
! It deals with molar abundances throughout (we expect that the input vector
! y has molar abundances Y = X / A) for make_rates and dydt routines.  It
! also checks to see if the temperature has changed much since the last call
! - if so, a new EOS call is made to get new values of cp and dhdX.
subroutine f_rhs(n, t, y, ydot, rpar, ipar)

  use bl_types
  use bl_constants_module, only: ZERO
  use network
  use eos_module
  use eos_type_module
  use rpar_indices

  implicit none

  integer,         intent(in   ) :: n, ipar
  real(kind=dp_t), intent(in   ) :: t, y(n)
  real(kind=dp_t), intent(inout) :: rpar(n_rpar_comps)
  real(kind=dp_t), intent(  out) :: ydot(n)

  integer :: k

  real(kind=dp_t), parameter :: T2T9 = 1.0e-9_dp_t

  real(kind=dp_t) :: dens, cp, dhdX(nspec), T9_eos, dT_crit, t9
  real(kind=dp_t) :: ymol(nspec)

  type(eos_t) :: eos_state

  ydot = ZERO

  ! several thermo vars via rpar
  dens = rpar(irp_dens)
  T9_eos = rpar(irp_T9_eos)
  dT_crit = rpar(irp_dTcrit)

  ymol = y(1:nspec)
  t9   = y(nspec+1)

  if (abs(t9 - T9_eos) > dT_crit*T9_eos) then
     T9_eos = t9

     eos_state % T   = T9_eos
     eos_state % rho = dens
     eos_state % Xn  = ymol * aion

     call eos(eos_input_rt, eos_state)

     rpar(irp_T9_eos) = T9_eos
     rpar(irp_cp) = eos_state % cp
     rpar(irp_dhdX:irp_dhdX+nspec-1) = eos_state % dhdX
  endif

  cp = rpar(irp_cp)
  dhdX = rpar(irp_dhdX:irp_dhdX+nspec-1)

  ! build the rates
  call make_rates(t9, dens, ymol, rpar)

  ! set up the ODEs
  call make_ydots(ymol, t9, rpar, ydot(1:nspec))

  ! t9 ODE
  ydot(n) = -sum((dhdX+bion)*ydot(1:nspec))/cp * T2T9
  
end subroutine f_rhs


subroutine make_rates(t9, dens, y, rpar)

  use bl_types
  use bl_constants_module, only: ZERO, THIRD, ONE, SIX
  use network
  use network_indices
  use rpar_indices

  implicit none

  real(kind=dp_t), intent(in   ) :: t9, dens, y(nspec)
  real(kind=dp_t), intent(inout) :: rpar(n_rpar_comps)

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
  rpar(irp_rates:irp_drtdt+nrates-1) = zero

  !**********************************************************************
  ! Start the rate calculations
  ! TODO - temperature derivatives for use in analytic jacobian
  !**********************************************************************
  ! helium burning - has been divided by 6
  rpar(irp_rates+ir3a-1) = 2.79d-8*t9m3*exp(-4.403_dp_t*t9m1)*dens*dens/SIX

  ! 15o(ag)19ne
  rpar(irp_rates+irag15-1) = (19.0_dp_t * (t9**2.85_dp_t) * &
       exp(-7.356_dp_t*t9m1) + &
       0.395_dp_t * t9m32 * exp(-5.849_dp_t*t9m1))*dens

  ! 14o(ap)17f
  ! an approx good below t9=0.5  is to drop first term, then
  !  rap14 = 3.31d+04*t9m32*exp(-11.733*t9m1)
  !   1  + 1.79d+07*t9m32*exp(-22.609/t9)+9000.*t9113*exp(-12.517/t9)
  rpar(irp_rates+irap14-1) = 1.68d13*t9m23 * &
       exp(-39.388d0*t9m13-(t9/0.717d0)**2) &
       *(ONE+0.011d0*t913+13.117d0*t923+0.971d0*t9+85.295d0*t943 &
       +16.06d0*t943)+3.31d4*t9m32*exp(-11.733d0*t9m1) &
       +1.79d7*t9m32*exp(-22.609d0*t9m1) &
       +9.00d+03*t9113*exp(-12.517d0*t9m1)
  rpar(irp_rates+irap14-1) = dens*rpar(irp_rates+irap14-1)

  ! 18ne(ap)21na
  rpar(irp_rates+irap18-1) = exp(56.59577d0-2.006856d0*t9m1 &
       +26.05828d0*t9m13 &
       -87.88732d0*t913 + 3.718809d0*t9 - 0.1767444d0*t953 &
       + 46.971960d0*t9log)
  rpar(irp_rates+irap18-1) = dens*rpar(irp_rates+irap18-1)

  ! weak rates
  rpar(irp_rates+irwk14o-1) = wk14o
  rpar(irp_rates+irwk15o-1) = wk15o
  
end subroutine make_rates


subroutine make_ydots(ymol, t9, rpar, dydt)

  use bl_types
  use bl_constants_module, only: ZERO, TWO, THREE, SIX
  use network
  use network_indices
  use rpar_indices

  implicit none

  real(kind=dp_t), intent(in   ) :: ymol(nspec), t9
  real(kind=dp_t), intent(inout) :: rpar(n_rpar_comps)
  real(kind=dp_t), intent(  out) :: dydt(nspec)

  real(kind=dp_t) :: dens

  ! initialize
  dydt = ZERO
  dens = rpar(irp_dens)

  ! From Stan:
  !       Reaction                   Rate
  !         3a + 2p   --> 14O          3a
  !        14O +      --> 18Ne       14O(ap)17F
  !     15O + a + 6p  --> 25Si       15O(ag)19Ne
  !   18Ne + a + 3p   --> 25Si       18Ne(ap)21Na
  !         14O + p   --> 15O        14O(e+nu)14N
  !         15O + 3p  --> 14O + a    15O(e+nu)15N

  ! o14
  dydt(io14) = ymol(ihe4)**3 * rpar(irp_rates+ir3a-1) &
       + ymol(io15) * rpar(irp_rates+irwk15o-1)       &
       - ymol(io14) * ymol(ihe4) * rpar(irp_rates+irap14-1)

  ! o15
  dydt(io15) = ymol(io14) * rpar(irp_rates+irwk14o-1) &
       - ymol(io15) * ymol(ihe4) * rpar(irp_rates+irag15-1) &
       - ymol(io15) * rpar(irp_rates+irwk15o-1)

  ! ne18
  dydt(ine18) = ymol(io14) * ymol(ihe4) * rpar(irp_rates+irap14-1) &
       - ymol(ine18) * ymol(ihe4) * rpar(irp_rates+irap18-1)

  ! si25
  dydt(isi25) = ymol(io15) * ymol(ihe4) * rpar(irp_rates+irag15-1) &
       + ymol(ine18) * ymol(ihe4) * rpar(irp_rates+irap18-1)

  ! he4
  dydt(ihe4) = ymol(io15) * rpar(irp_rates+irwk15o-1) &
       - THREE * ymol(ihe4)**3 * rpar(irp_rates+ir3a-1) &
       - ymol(io14) * ymol(ihe4) * rpar(irp_rates+irap14-1) &
       - ymol(io15) * ymol(ihe4) * rpar(irp_rates+irag15-1) &
       - ymol(ine18) * ymol(ihe4) * rpar(irp_rates+irap18-1)

  ! h1
  dydt(ih1) = - TWO * ymol(ihe4)**3 * rpar(irp_rates+ir3a-1) &
       - SIX * ymol(io15) * ymol(ihe4) * rpar(irp_rates+irag15-1) &
       - THREE * ymol(ine18) * ymol(ihe4) * rpar(irp_rates+irap18-1) &
       - ymol(io14) * rpar(irp_rates+irwk14o-1) &
       - THREE * ymol(io15) * rpar(irp_rates+irwk15o-1)
  
  ! iron is just a tracer
  dydt(ife56) = ZERO
  
end subroutine make_ydots

! stub for jac
! TODO - make this an analytic jacobian
subroutine jac(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)

  use bl_types
  use bl_constants_module, only: ZERO
  
  implicit none

  integer,         intent(in   ) :: neq, ml, mu, nrpd, ipar
  real(kind=dp_t), intent(in   ) :: y(neq), t
  real(kind=dp_t), intent(inout) :: rpar(*)
  real(kind=dp_t), intent(  out) :: pd(neq,neq)

  pd(:,:) = ZERO
  
end subroutine jac
