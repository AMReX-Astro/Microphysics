! The f_rhs routine provides the right-hand-side for the DVODE solver.  It 
! converts the mass fractions into molar abundances before calling the 
! make_rates, screen, and dydt routines.  It also checks to see if the 
! temperature has changed much since the last call - if so, it updates the 
! temperature to get a better estimate of the reaction rates.
!
! The jac routine provides an explicit Jacobian to the DVODE solver.

subroutine f_rhs(n, time, y, ydot, rpar, ipar)

  use bl_types
  use bl_constants_module, only: ZERO, ONE
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
  double precision, intent(INOUT) :: rpar(n_rpar_comps)
  double precision, intent(  OUT) :: ydot(n)

  double precision :: rates(nrates)

  integer :: k

  type (eos_t) :: state
 
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
  state % abar    = rpar(irp_abar)
  state % zbar    = rpar(irp_zbar)
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

  rpar(irp_dhdX:irp_dhdX-1+nspec) = state % dhdX
  rpar(irp_dedX:irp_dedX-1+nspec) = state % dedX
  rpar(irp_cp) = state % cp
  rpar(irp_cv) = state % cv
  rpar(irp_abar) = state % abar
  rpar(irp_zbar) = state % zbar
  
  call normalize_abundances(state)
  
  ! Call the aprox13 routines to get dY/dt and de/dt.

  call aprox13(time,state,ydot,rpar)

end subroutine f_rhs
  


! Analytical Jacobian

subroutine jac(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)

  use network, only: nrates, nspec
  use bl_types
  use bl_constants_module, only: ZERO
  use network_indices
  use rpar_indices
  use rhs_module
  use eos_module
  use extern_probin_module, only: do_constant_volume_burn
  
  implicit none

  integer         , intent(IN   ) :: neq, ml, mu, nrpd, ipar
  double precision, intent(IN   ) :: y(neq), rpar(n_rpar_comps), t
  double precision, intent(  OUT) :: pd(neq,neq)
   
  double precision :: rates(nrates), ydot(nspec)

  double precision :: b1,sneut,dsneutdt,dsneutdd,snuda,snudz  

  integer :: i, j
  
  double precision :: rho, temp, cv, cp, abar, zbar, dEdX(nspec), dhdX(nspec)
  
  pd(:,:) = ZERO

  rho  = rpar(irp_dens)
  temp = y(net_itemp)
  cv   = rpar(irp_cv)
  cp   = rpar(irp_cp)
  
  ! Get the data from rpar

  dhdX = rpar(irp_dhdX:irp_dhdX+nspec-1)
  dEdX = rpar(irp_dEdX:irp_dEdX+nspec-1)
  abar = rpar(irp_abar)
  zbar = rpar(irp_zbar)
  
  ! Note that this RHS has been evaluated using rates = d(ratdum) / dT
  
  ydot = rpar(irp_dydt:irp_dydt+nspec-1)
  rates = rpar(irp_dratesdt:irp_dratesdt+nrates-1)
  
  ! Species Jacobian elements with respect to other species
  
  call dfdy_isotopes_aprox13(y, pd, neq, rates)

  ! Species Jacobian elements with respect to temperature

  pd(1:nspec,net_itemp) = ydot

  if (rpar(irp_self_heat) > ZERO) then
  
     ! Energy generation rate Jacobian elements

     do j = 1, nspec
        call ener_gener_rate(pd(1:nspec,j) / aion,pd(net_ienuc,j))
     enddo
     call ener_gener_rate(pd(1:nspec,net_itemp) / aion, pd(net_ienuc,net_itemp))

     ! Account for the thermal neutrino losses
     call sneut5(T,rho,abar,zbar,sneut,dsneutdt,dsneutdd,snuda,snudz)

     do j = 1, nspec
        b1 = ((aion(j) - abar) * abar * snuda + (zion(j) - zbar) * abar * snudz)
        pd(net_ienuc,j) = pd(net_ienuc,j) - b1
     enddo
     pd(net_ienuc,net_itemp) = pd(net_ienuc,net_itemp) - dsneutdt

     ! Temperature Jacobian elements
     
     if (do_constant_volume_burn) then
        
        ! d(itemp)/d(yi)
        do j = 1, nspec
           pd(net_itemp,j) = ( pd(net_ienuc,j) - sum( dEdX(:) * pd(1:nspec,j) ) ) / cv
        enddo
           
        ! d(itemp)/d(temp)
        pd(net_itemp,net_itemp) = ( pd(net_ienuc,net_itemp) - sum( dEdX(:) * pd(1:nspec,net_itemp) ) ) / cv
     
     else

        ! d(itemp)/d(yi)
        do j = 1, nspec
           pd(net_itemp,j) = ( pd(net_ienuc,j) - sum( dhdX(:) * pd(1:nspec,j) ) ) / cp
        enddo

        ! d(itemp)/d(temp)
        pd(net_itemp,net_itemp) = ( pd(net_ienuc,net_itemp) - sum( dhdX(:) * pd(1:nspec,net_itemp) ) ) / cp

     endif
        
  endif
     
end subroutine jac
