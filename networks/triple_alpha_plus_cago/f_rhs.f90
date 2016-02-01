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
  use eos_module, only: eos_input_rt, eos
  use eos_type_module
  use rates_module
  use screen_module
  use dydt_module
  use rpar_indices

  implicit none

  integer,         intent(IN   ) :: n, ipar
  real(kind=dp_t), intent(IN   ) :: t, y(n)
  real(kind=dp_t), intent(INOUT) :: rpar(*)
  real(kind=dp_t), intent(  OUT) :: ydot(n)

  integer :: k

  real(kind=dp_t) :: temp

  real(kind=dp_t) :: dens, c_p, dhdx(nspec), T_eos, dT_crit

  real(kind=dp_t) :: rates(nrat), dratesdt(nrat)
  
  real(kind=dp_t) :: ymol(nspec)

  type (eos_t) :: eos_state

  ydot = ZERO
  
  ! several thermodynamic quantities come in via rpar
  dens    = rpar(irp_dens)
  c_p     = rpar(irp_cp)
  dhdX(:) = rpar(irp_dhdX:irp_dhdX-1+nspec)
  T_eos   = rpar(irp_Teos)
  dT_crit = rpar(irp_Tcrit)

  temp = y(n)

  if (abs(temp - T_eos) > dT_crit*T_eos) then

     T_eos = temp

     eos_state%rho  = dens
     eos_state%T = temp

     eos_state%xn(ihe4_) = y(ihe4_)
     eos_state%xn(ic12_) = y(ic12_)
     eos_state%xn(io16_) = y(io16_)
     eos_state%xn(ife56_) = rpar(irp_Y56)*aion(ife56_)

     call eos(eos_input_rt, eos_state)

     c_p = eos_state%cp
     dhdx(:) = eos_state%dhdx(:)

     rpar(irp_cp) = c_p
     rpar(irp_dhdX:irp_dhdX-1+nspec) = dhdX(:)
     rpar(irp_Teos) = temp

  endif

  ! calculate ymol(:)
  do k = 1, nevolve
     ymol(k) = y(k) / aion(k)
  enddo
  ymol(ife56_) = rpar(irp_Y56)

  ! set up the species ODEs for the reaction network
  ! species inputs are in molar fractions but come out in mass fractions
  call make_rates(temp, dens, rates, dratesdt)

  call screen(temp, dens, ymol, rates, dratesdt)

  call dydt(ymol, rates, ydot(1:nevolve))

  rpar(irp_rates:irp_rates-1+nrat) = rates(:)
  rpar(irp_drtdt:irp_drtdt-1+nrat) = dratesdt(:)


  ! set up the temperature ODE
  ! dT/dt = -(1/c_p) sum_i (xi_i + q_i) omega_i
  do k = 1, nevolve
     ydot(n) = ydot(n) - (dhdx(k) + ebin(k))*ydot(k)
  enddo

  ydot(n) = ydot(n) / c_p


  return

end subroutine f_rhs
  



subroutine jac(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)

  use bl_types
  use bl_constants_module
  use network
  use dydt_module
  use rpar_indices

  implicit none

  integer        , intent(IN   ) :: neq, ml, mu, nrpd, ipar
  real(kind=dp_t), intent(IN   ) :: y(neq), rpar(*), t
  real(kind=dp_t), intent(  OUT) :: pd(neq,neq)

  real(kind=dp_t) :: ymol(nspec)

  real(kind=dp_t) :: rates(nrat), dratesdt(nrat)

  real(kind=dp_t) :: dens, c_p, dhdx(nspec), T_eos, dT_crit

  integer :: itemp, i, j

  rates(:)    = rpar(irp_rates:irp_rates-1+nrat)
  dratesdt(:) = rpar(irp_drtdt:irp_drtdt-1+nrat) 

  ! several thermodynamic quantities come in via rpar
  dens    = rpar(irp_dens)
  c_p     = rpar(irp_cp)
  dhdX(:) = rpar(irp_dhdX:irp_dhdX-1+nspec)
  T_eos   = rpar(irp_Teos)
  dT_crit = rpar(irp_Tcrit)


  ! initialize
  pd(:,:)  = ZERO

  ymol(1:nevolve) = y(1:nevolve) / aion(1:nevolve)
  ymol(ife56_) = rpar(irp_Y56)

  itemp = neq

  ! ======================================================================
  ! THESE ARE IN TERMS OF MOLAR FRACTIONS

  ! helium jacobian elements
  pd(ihe4_,ihe4_)  = - NINE * ymol(ihe4_) * ymol(ihe4_) * rates(ir3a_) &
                     - ONE * ymol(ic12_) * rates(ircago_)
  pd(ihe4_,ic12_)  = - ONE * ymol(ihe4_) * rates(ircago_)

  ! carbon jacobian elements
  pd(ic12_,ihe4_) =   THREE * ymol(ihe4_) * ymol(ihe4_) * rates(ir3a_) &
                    - ONE * ymol(ic12_) * rates(ircago_)
  pd(ic12_,ic12_) = - ONE * ymol(ihe4_) * rates(ircago_)

  ! oxygen jacobian elements
  pd(io16_,ihe4_) = ONE * ymol(ic12_) * rates(ircago_)
  pd(io16_,ic12_) = ONE * ymol(ihe4_) * rates(ircago_)

  ! ======================================================================

  ! convert to mass fractions
  do j = 1,nevolve
     pd(j,1:nevolve) = pd(j,1:nevolve)*aion(j)
     pd(1:nevolve,j) = pd(1:nevolve,j)/aion(j)
  enddo

  ! add the temperature derivatives: df(y_i) / dT
  ! dydt automatically converts the "ydot"'s to mass fractions
  call dydt(ymol, dratesdt, pd(1:nevolve,itemp))


  ! add the temperature jacobian elements df(T) / df(y)
  do j = 1, nevolve
     do i = 1, nevolve
        pd(itemp,j) = pd(itemp,j) - (dhdX(i) + ebin(i))*pd(i,j)
     enddo
  enddo


  ! add df(T) / dT
  do i = 1, nevolve
     pd(itemp,itemp) = pd(itemp,itemp) - (dhdX(i) + ebin(i))*pd(i,itemp)
  enddo


  ! need to divide the temperature jacobian elements by cp
  pd(itemp,:) =   pd(itemp,:) / c_p


  return
end subroutine jac
