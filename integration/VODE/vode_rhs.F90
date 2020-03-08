module vode_rhs_module

contains
  
  ! The f_rhs routine provides the right-hand-side for the DVODE solver.
  ! This is a generic interface that calls the specific RHS routine in the
  ! network you're actually using.

  subroutine f_rhs(time, y, ydot, rpar)

    !$acc routine seq
    
    use actual_network, only: aion, nspec_evolve
    use amrex_fort_module, only: rt => amrex_real
    use burn_type_module, only: burn_t, net_ienuc, net_itemp
    use amrex_constants_module, only: ZERO, ONE
    use network_rhs_module, only: network_rhs
    use extern_probin_module, only: renormalize_abundances, &
         integrate_temperature, integrate_energy, react_boost
    use vode_type_module, only: clean_state, renormalize_species, update_thermodynamics, burn_to_vode, vode_to_burn, VODE_NEQS
    use vode_rpar_indices, only: n_rpar_comps, irp_t_sound, irp_t0

    implicit none

    real(rt), intent(INOUT) :: time, y(VODE_NEQS)
    real(rt), intent(INOUT) :: rpar(n_rpar_comps)
    real(rt), intent(INOUT) :: ydot(VODE_NEQS)

    type (burn_t) :: burn_state

    !$gpu

    ! We are integrating a system of
    !
    ! y(1:nspec_evolve) = dX/dt
    ! y(net_itemp) = dT/dt
    ! y(net_ienuc) = denuc/dt

    ! Fix the state as necessary.

    call clean_state(y, rpar)

    ! Update the thermodynamics as necessary.

    call update_thermodynamics(y, rpar)

    ! Call the specific network routine to get the RHS.

    call vode_to_burn(y, rpar, burn_state)

    burn_state % time = time
    call network_rhs(burn_state, ydot, rpar(irp_t0))

    ! We integrate X, not Y
    ydot(1:nspec_evolve) = ydot(1:nspec_evolve) * aion(1:nspec_evolve)

    ! Allow temperature and energy integration to be disabled.
    if (.not. integrate_temperature) then
       ydot(net_itemp) = ZERO
    endif

    if (.not. integrate_energy) then
       ydot(net_ienuc) = ZERO
    endif

    ! apply fudge factor:
    if (react_boost > ZERO) then
       ydot(:) = react_boost * ydot(:)
    endif

    !call burn_to_vode(burn_state, y, rpar)

  end subroutine f_rhs



  ! Analytical Jacobian

  subroutine jac(time, y, ml, mu, pd, nrpd, rpar)

    !$acc routine seq
    
    use network, only: aion, aion_inv, nspec_evolve
    use amrex_constants_module, only: ZERO
    use network_rhs_module, only: network_jac
    use burn_type_module, only: burn_t, net_ienuc, net_itemp
    use vode_type_module, only: vode_to_burn, burn_to_vode, VODE_NEQS
    use vode_rpar_indices, only: n_rpar_comps, irp_t_sound, irp_t0
    use amrex_fort_module, only: rt => amrex_real
    use extern_probin_module, only: integrate_temperature, integrate_energy, react_boost

    implicit none

    integer   , intent(IN   ) :: ml, mu, nrpd
    real(rt), intent(INOUT) :: y(VODE_NEQS), rpar(n_rpar_comps), time
    real(rt), intent(  OUT) :: pd(VODE_NEQS,VODE_NEQS)

    type (burn_t) :: state
    integer :: n

    !$gpu

    ! Call the specific network routine to get the Jacobian.

    call vode_to_burn(y, rpar, state)
    state % time = time
    call network_jac(state, pd, rpar(irp_t0))

    ! We integrate X, not Y
    do n = 1, nspec_evolve
       pd(n,:) = pd(n,:) * aion(n)
       pd(:,n) = pd(:,n) * aion_inv(n)
    enddo

    ! apply fudge factor:
    if (react_boost > ZERO) then
       pd(:,:) = react_boost * pd(:,:)
    endif

    ! Allow temperature and energy integration to be disabled.
    if (.not. integrate_temperature) then
       pd(net_itemp,:) = ZERO
    endif

    if (.not. integrate_energy) then
       pd(net_ienuc,:) = ZERO
    endif

    !call burn_to_vode(state, y, rpar)

  end subroutine jac
end module vode_rhs_module
