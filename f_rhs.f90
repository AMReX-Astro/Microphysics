! The f_rhs routine provides the right-hand-side for the DVODE solver.  It 
! converts the mass fractions into molar abundances before calling the 
! make_rates, screen, and dydt routines.  It also checks to see if the 
! temperature has changed much since the last call - if so, it updates the 
! temperature to get a better estimate of the reaction rates.
!
! The jac routine provides an explicit Jacobian to the DVODE solver.
!

subroutine f_rhs(n, t, y, ydot, rpar, ipar)

  use bl_types
  use bl_constants_module, only: ZERO
  use network
  use rpar_indices
  use network_indices
  use rhs_module, only: aprox13
  use extern_probin_module, only: do_constant_volume_burn
  implicit none

  integer,         intent(IN   ) :: n, ipar
  real(kind=dp_t), intent(IN   ) :: t, y(n)
  real(kind=dp_t), intent(INOUT) :: rpar(*)
  real(kind=dp_t), intent(  OUT) :: ydot(n)

  integer :: k

  real(kind=dp_t) :: temp

  real(kind=dp_t) :: dens, c_p, c_v, dhdX(nspec), dedX(nspec)
  real(kind=dp_t) :: ymol(nspec), dymoldt(nspec)
  real(kind=dp_t) :: denucdt
  real(kind=dp_t) :: smallx

  ! we are integrating a system of
  !
  ! y(1:nspec) = dX/dt  
  ! y(itemp)   = dT/dt
  ! y(ienuc)   = denuc/dt

  
  ydot = ZERO
  
  ! several thermodynamic quantities come in via rpar -- note: these
  ! are evaluated at the start of the integration, so if things change
  ! dramatically, they will fall out of sync with the current
  ! thermodynamcis
  dens    = rpar(irp_dens)
  c_p     = rpar(irp_cp)
  c_v     = rpar(irp_cv)
  dhdX(:) = rpar(irp_dhdX:irp_dhdX-1+nspec)
  dedX(:) = rpar(irp_dedX:irp_dedX-1+nspec)
  smallx  = rpar(irp_smallx)


  ! temperature is one of the quantities that we are integrating --
  ! always use the current T
  temp = y(itemp)


  ! calculate molar fractions from the input mass fractions
  do k = 1, nspec
     ymol(k) = y(k) / aion(k)
  enddo

  ! call the aprox13 routines to get dY/dt
  call aprox13(t,temp,dens,ymol,dymoldt,denucdt,smallx)

  ! we are evolving mass fraction equations -- setup the RHS for the
  ! species abundances
  ydot(1:nspec) = aion(1:nspec)*dymoldt(1:nspec)

  ! nuclear energy ODE -- this just gets us the total nuclear energy
  ! release when we are all done
  ydot(ienuc) = denucdt

  ! set up the temperature ODE.  For constant pressure, Dp/Dt = 0, we
  ! evolve :
  !    dT/dt = (1/c_p) [ -sum_i (xi_i omega_i) + Hnuc]
  ! 
  ! For constant volume, div{U} = 0, and we evolve:
  !    dT/dt = (1/c_v) [ -sum_i ( {e_x}_i omega_i) + Hnuc]
  !
  ! see paper III, including Eq. A3 for details.

  if (do_constant_volume_burn) then
     ydot(itemp) = 0.0d0
     do k = 1, nspec
        ydot(itemp) = ydot(itemp) - dedX(k)*ydot(k) 
     enddo
     ydot(itemp) = ydot(itemp) + denucdt

     ydot(itemp) = ydot(itemp) / c_v

  else
     ydot(itemp) = 0.0d0
     do k = 1, nspec
        ydot(itemp) = ydot(itemp) - dhdX(k)*ydot(k) 
     enddo
     ydot(itemp) = ydot(itemp) + denucdt

     ydot(itemp) = ydot(itemp) / c_p
  endif

  return
end subroutine f_rhs
  

subroutine jac(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)

  use bl_types
  use bl_constants_module, only: ZERO
  implicit none

  integer        , intent(IN   ) :: neq, ml, mu, nrpd, ipar
  real(kind=dp_t), intent(IN   ) :: y(neq), rpar(*), t
  real(kind=dp_t), intent(  OUT) :: pd(neq,neq)

  pd(:,:) = ZERO

  return
end subroutine jac
