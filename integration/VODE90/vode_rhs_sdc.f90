
! The f_rhs routine provides the right-hand-side for the DVODE solver.
! This is a generic interface that calls the specific RHS routine in the
! network you're actually using.

subroutine f_rhs(neq, time, y, ydot, rpar, ipar)

  use actual_network, only: aion, nspec_evolve
  use amrex_fort_module, only: rt => amrex_real
  use burn_type_module, only: burn_t, net_ienuc, net_itemp
  use amrex_constants_module, only: ZERO, ONE
  use actual_rhs_module, only: actual_rhs
  use vode_type_module, only: clean_state, renormalize_species, &
       rhs_to_vode, vode_to_burn
  use vode_rpar_indices

  implicit none

  integer,    intent(IN   ) :: neq, ipar
  real(rt), intent(INOUT) :: time, y(neq)
  real(rt), intent(INOUT) :: rpar(n_rpar_comps)
  real(rt), intent(  OUT) :: ydot(neq)

  type (burn_t) :: burn_state

  ydot = ZERO

  ! Fix the state as necessary.

  call clean_state(time, y, rpar)

  ! convert to the burn_t

  call vode_to_burn(time, y, rpar, burn_state)

  burn_state % time = time

  ! call the specific network to get the RHS

  call actual_rhs(burn_state)

  ! convert back to the vode type -- this will add the advective terms

  call rhs_to_vode(time, burn_state, y, ydot, rpar)

end subroutine f_rhs



! Analytical Jacobian

subroutine jac(neq, time, y, ml, mu, pd, nrpd, rpar, ipar)

  ! with VODE, we will only ever enter this routine if we choose to
  ! use an analytic Jacobian.  Otherwise, VODE will use its internal
  ! Jacobian routines.

  use network, only: aion, aion_inv, nspec_evolve
  use amrex_constants_module, only: ZERO
  use actual_rhs_module, only: actual_jac
  use burn_type_module, only: burn_t, net_ienuc
  use amrex_fort_module, only: rt => amrex_real
  use vode_rpar_indices
  use vode_type_module, only: clean_state, renormalize_species, &
       jac_to_vode, vode_to_burn


  implicit none

  integer   , intent(IN   ) :: neq, ml, mu, nrpd, ipar
  real(rt), intent(INOUT) :: y(neq), rpar(n_rpar_comps), time
  real(rt), intent(  OUT) :: pd(neq,neq)

  type (burn_t) :: state
  integer :: n

  pd(:,:) = ZERO

  ! Call the specific network routine to get the Jacobian.
  call vode_to_burn(time, y, rpar, state)
  state % time = time

  ! call the analytic Jacobian
  call actual_jac(state)

  ! convert to the system we are using
  call jac_to_vode(time, state, y, pd, rpar)

end subroutine jac
