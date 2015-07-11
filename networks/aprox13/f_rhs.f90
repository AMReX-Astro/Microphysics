! The f_rhs routine provides the right-hand-side for the DVODE solver.  It 
! converts the mass fractions into molar abundances before calling the 
! make_rates, screen, and dydt routines.  It also checks to see if the 
! temperature has changed much since the last call - if so, it updates the 
! temperature to get a better estimate of the reaction rates.
!
! The jac routine provides an explicit Jacobian to the DVODE solver.

subroutine f_rhs(n, time, y, ydot, rpar, ipar)

  use bl_types
  use bl_constants_module, only: ZERO
  use network
  use rpar_indices
  use network_indices
  use rhs_module, only: aprox13
  use eos_type_module
  use eos_module
  use extern_probin_module, only: call_eos_in_rhs
  
  implicit none

  integer,          intent(IN   ) :: n, ipar
  double precision, intent(IN   ) :: time, y(n)
  double precision, intent(INOUT) :: rpar(*)
  double precision, intent(  OUT) :: ydot(n)

  integer :: k

  type(eos_t) :: state

  ! We are integrating a system of
  !
  ! y(1:nspec)   = dX/dt  
  ! y(net_itemp) = dT/dt
  ! y(net_ienuc) = denuc/dt
  
  ydot = ZERO
  
  ! Several thermodynamic quantities come in via rpar -- note: these
  ! are evaluated at the start of the integration, so if things change
  ! dramatically, they will fall out of sync with the current
  ! thermodynamics.
  state % rho     = rpar(irp_dens)
  state % cp      = rpar(irp_cp)
  state % cv      = rpar(irp_cv)
  state % xn(:)   = y(1:nspec)
  state % dhdX(:) = rpar(irp_dhdX:irp_dhdX-1+nspec)
  state % dedX(:) = rpar(irp_dedX:irp_dedX-1+nspec)

  ! Temperature is one of the quantities that we are integrating --
  ! always use the current T.
  state % T       = y(net_itemp)

  ! Evaluate the thermodynamics -- if desired.
  ! Otherwise just do the composition calculations since
  ! that's needed to construct dY/dt. Then make sure
  ! the abundances are safe.
  
  if (call_eos_in_rhs) then
     call eos(eos_input_rt, state)     
  else
     call composition(state)
  endif

  state % xn(:) = min(1.0d0,max(state % xn(:),smallx))
  
  ! Call the aprox13 routines to get dY/dt and de/dt.

  call aprox13(time,state,ydot)

end subroutine f_rhs
  

subroutine jac(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)

  use bl_types
  use bl_constants_module, only: ZERO
  implicit none

  integer        , intent(IN   ) :: neq, ml, mu, nrpd, ipar
  double precision, intent(IN   ) :: y(neq), rpar(*), t
  double precision, intent(  OUT) :: pd(neq,neq)

  pd(:,:) = ZERO

end subroutine jac
