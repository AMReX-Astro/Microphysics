module vode_rhs_module

contains
  
  ! The f_rhs routine provides the right-hand-side for the DVODE solver.
  ! This is a generic interface that calls the specific RHS routine in the
  ! network you're actually using.

#ifdef CUDA
  attributes(device) &
#endif
  subroutine f_rhs(neq, time, y, ydot, rpar, ipar)

    !$acc routine seq
    
    use actual_network, only: aion, nspec_evolve
    use bl_types, only: dp_t
    use burn_type_module, only: burn_t, net_ienuc, net_itemp
    use bl_constants_module, only: ZERO, ONE
    use actual_rhs_module, only: actual_rhs
    use managed_probin_module, only: cu_call_eos_in_rhs, cu_dT_crit, cu_renormalize_abundances, &
         cu_burning_mode, cu_burning_mode_factor, &
         cu_integrate_temperature, cu_integrate_energy
    use vode_type_module, only: clean_state, renormalize_species, update_thermodynamics, &
                                burn_to_vode, vode_to_burn
    use rpar_indices, only: n_rpar_comps, irp_y_init, irp_t_sound

    implicit none

    integer,    intent(IN   ) :: neq, ipar(:)
    real(dp_t), intent(INOUT) :: time, y(:)
    real(dp_t), intent(INOUT) :: rpar(:)
    real(dp_t), intent(INOUT) :: ydot(:)

    type (burn_t) :: burn_state

    real(dp_t) :: limit_factor, t_sound, t_enuc

    ! We are integrating a system of
    !
    ! y(1:nspec_evolve) = dX/dt
    ! y(net_itemp) = dT/dt
    ! y(net_ienuc) = denuc/dt

    ydot = ZERO

    ! Fix the state as necessary.

    call clean_state(y, rpar)

    ! Renormalize the abundances as necessary.

    if (cu_renormalize_abundances) then
       call renormalize_species(y, rpar)
    endif

    ! Update the thermodynamics as necessary.

    call update_thermodynamics(y, rpar)

    ! Call the specific network routine to get the RHS.

    call vode_to_burn(y, rpar, burn_state)

    burn_state % time = time
    call actual_rhs(burn_state)

    ! We integrate X, not Y
    burn_state % ydot(1:nspec_evolve) = &
         burn_state % ydot(1:nspec_evolve) * aion(1:nspec_evolve)

    ! Allow temperature and energy integration to be disabled.
    if (.not. cu_integrate_temperature) then
       burn_state % ydot(net_itemp) = ZERO
    endif

    if (.not. cu_integrate_energy) then
       burn_state % ydot(net_ienuc) = ZERO
    endif

    ! For burning_mode == 3, limit the rates.
    ! Note that we are limiting with respect to the initial zone energy.

    if (cu_burning_mode == 3) then

       t_enuc = rpar(irp_y_init + net_ienuc - 1) / max(abs(burn_state % ydot(net_ienuc)), 1.d-50)
       t_sound = rpar(irp_t_sound)

       limit_factor = min(1.0d0, cu_burning_mode_factor * t_enuc / t_sound)

       burn_state % ydot = limit_factor * burn_state % ydot

    endif

    call burn_to_vode(burn_state, y, rpar, ydot = ydot)

  end subroutine f_rhs



  ! Analytical Jacobian
#ifdef CUDA
  attributes(device) &
#endif
  subroutine jac(neq, time, y, ml, mu, pd, nrpd, rpar, ipar)

    !$acc routine seq
    
    use network, only: aion, aion_inv, nspec_evolve
    use bl_constants_module, only: ZERO
    use actual_rhs_module, only: actual_jac
    use burn_type_module, only: burn_t, net_ienuc, net_itemp
    use vode_type_module, only: vode_to_burn, burn_to_vode
    use rpar_indices, only: n_rpar_comps, irp_y_init, irp_t_sound
    use bl_types, only: dp_t
    use managed_probin_module, only: cu_burning_mode, cu_burning_mode_factor, &
         cu_integrate_temperature, cu_integrate_energy

    implicit none

    integer   , intent(IN   ) :: neq, ml, mu, nrpd, ipar(:)
    real(dp_t), intent(INOUT) :: y(neq), rpar(n_rpar_comps), time
    real(dp_t), intent(  OUT) :: pd(:,:)

    type (burn_t) :: state
    real(dp_t) :: limit_factor, t_sound, t_enuc
    integer :: n

    ! Call the specific network routine to get the Jacobian.

    call vode_to_burn(y, rpar, state)
    state % time = time
    call actual_jac(state)

    ! We integrate X, not Y
    do n = 1, nspec_evolve
       state % jac(n,:) = state % jac(n,:) * aion(n)
       state % jac(:,n) = state % jac(:,n) * aion_inv(n)
    enddo

    ! Allow temperature and energy integration to be disabled.
    if (.not. cu_integrate_temperature) then
       state % jac(net_itemp,:) = ZERO
    endif

    if (.not. cu_integrate_energy) then
       state % jac(net_ienuc,:) = ZERO
    endif

    ! For burning_mode == 3, limit the rates.
    ! Note that we are limiting with respect to the initial zone energy.

    if (cu_burning_mode == 3) then

       t_enuc = rpar(irp_y_init + net_ienuc - 1) / max(abs(state % ydot(net_ienuc)), 1.d-50)
       t_sound = rpar(irp_t_sound)

       limit_factor = min(1.0d0, cu_burning_mode_factor * t_enuc / t_sound)

       state % jac = limit_factor * state % jac

    endif

    call burn_to_vode(state, y, rpar, jac = pd)

  end subroutine jac
end module vode_rhs_module
