module vode_rhs_module

  use cuvode_types_module, only : dvode_t
  contains

! The f_rhs routine provides the right-hand-side for the DVODE solver.
! This is a generic interface that calls the specific RHS routine in the
! network you're actually using.

subroutine f_rhs(time, vode_state, ydot)

  use actual_network, only: aion, nspec
  use amrex_fort_module, only: rt => amrex_real
  use burn_type_module, only: burn_t, net_ienuc, net_itemp, neqs
  use amrex_constants_module, only: ZERO, ONE
  use network_rhs_module, only: network_rhs
  use vode_type_module, only: clean_state, renormalize_species
  use vode_type_module, only: rhs_to_vode, vode_to_burn
  use vode_rpar_indices
  use cuvode_parameters_module

  implicit none

  real(rt), intent(INOUT) :: time
  type(dvode_t), intent(INOUT) :: vode_state
  real(rt), intent(  OUT) :: ydot(VODE_NEQS)

  type (burn_t) :: burn_state
  real(rt) :: ydot_react(neqs)

  !$gpu

  ydot = ZERO

  ! Fix the state as necessary.

  call clean_state(time, vode_state)

  ! convert to the burn_t

  call vode_to_burn(time, vode_state, burn_state)

  burn_state % time = time

  ! call the specific network to get the RHS

  call network_rhs(burn_state, ydot_react, vode_state % rpar(irp_t0))

  ! convert back to the vode type -- this will add the advective terms

  call rhs_to_vode(time, burn_state, ydot_react, vode_state, ydot)

end subroutine f_rhs



! Analytical Jacobian

subroutine jac(time, vode_state, ml, mu, pd, nrpd)

  ! with VODE, we will only ever enter this routine if we choose to
  ! use an analytic Jacobian.  Otherwise, VODE will use its internal
  ! Jacobian routines.

  use network, only: aion, aion_inv, nspec
  use amrex_constants_module, only: ZERO
  use network_rhs_module, only: network_jac
  use burn_type_module, only: burn_t, net_ienuc, neqs
  use amrex_fort_module, only: rt => amrex_real
  use vode_rpar_indices
  use vode_type_module, only: clean_state, renormalize_species
  use vode_type_module, only: jac_to_vode, vode_to_burn
  use cuvode_parameters_module

  implicit none

  integer   , intent(IN   ) :: ml, mu, nrpd
  real(rt), intent(IN) :: time
  real(rt), intent(  OUT) :: pd(VODE_NEQS, VODE_NEQS)
  type (dvode_t), intent(inout) :: vode_state

  type (burn_t) :: state
  real(rt) :: jac_react(neqs, neqs)
  integer :: n

  !$gpu

  pd(:,:) = ZERO

  ! Call the specific network routine to get the Jacobian.
  call vode_to_burn(time, vode_state, state)
  state % time = time

  ! call the analytic Jacobian
  call network_jac(state, jac_react, vode_state % rpar(irp_t0))

  ! convert to the system we are using
  call jac_to_vode(time, jac_react, vode_state, pd)

end subroutine jac

end module vode_rhs_module
