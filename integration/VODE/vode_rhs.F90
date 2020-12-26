module vode_rhs_module

  use cuvode_types_module, only : dvode_t

  implicit none

contains
  
  ! The f_rhs routine provides the right-hand-side for the DVODE solver.
  ! This is a generic interface that calls the specific RHS routine in the
  ! network you're actually using.

  subroutine f_rhs(time, vode_state, ydot)

    !$acc routine seq
    
    use actual_network, only: aion, nspec
    use amrex_fort_module, only: rt => amrex_real
    use eos_type_module, only : eos_get_small_temp
    use burn_type_module, only: burn_t, net_ienuc, net_itemp
    use amrex_constants_module, only: ZERO, ONE
    use network_rhs_module, only: network_rhs
    use extern_probin_module, only: integrate_temperature, integrate_energy, react_boost, MAX_TEMP
    use vode_type_module, only: clean_state, update_thermodynamics, burn_to_vode, vode_to_burn, VODE_NEQS
    use vode_rpar_indices, only: n_rpar_comps, irp_t0

    implicit none

    real(rt), intent(INOUT) :: time
    type(dvode_t), intent(INOUT) :: vode_state
    real(rt), intent(INOUT) :: ydot(VODE_NEQS)

    type (burn_t) :: burn_state

    real(rt) :: small_temp

    !$gpu

    ! We are integrating a system of
    !
    ! y(1:nspec)   = dX/dt
    ! y(net_itemp) = dT/dt
    ! y(net_ienuc) = denuc/dt

    ! Only do the burn if the incoming temperature is within the temperature
    ! bounds. Otherwise set the RHS to zero and return.

    call eos_get_small_temp(small_temp)

    if (vode_state % y(net_itemp) <= small_temp .or. vode_state % y(net_itemp) >= MAX_TEMP) then

       ydot(:) = 0.0_rt

       return

    end if

    ! Fix the state as necessary.

    call clean_state(vode_state)

    ! Update the thermodynamics as necessary.

    call update_thermodynamics(vode_state)

    ! Call the specific network routine to get the RHS.

    call vode_to_burn(vode_state, burn_state)

    burn_state % time = time
    call network_rhs(burn_state, ydot, vode_state % rpar(irp_t0))

    ! We integrate X, not Y
    ydot(1:nspec) = ydot(1:nspec) * aion(1:nspec)

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

    call burn_to_vode(burn_state, vode_state)

  end subroutine f_rhs



  ! Analytical Jacobian

  subroutine jac(time, vode_state, ml, mu, pd, nrpd)

    !$acc routine seq
    
    use network, only: aion, aion_inv, nspec
    use amrex_constants_module, only: ZERO
    use network_rhs_module, only: network_jac
    use eos_type_module, only : eos_get_small_temp
    use burn_type_module, only: burn_t, net_ienuc, net_itemp
    use vode_type_module, only: vode_to_burn, burn_to_vode, VODE_NEQS
    use vode_rpar_indices, only: n_rpar_comps, irp_t0
    use amrex_fort_module, only: rt => amrex_real
    use extern_probin_module, only: integrate_temperature, integrate_energy, react_boost, MAX_TEMP

    implicit none

    integer   , intent(IN   ) :: ml, mu, nrpd
    real(rt), intent(IN) :: time
    real(rt), intent(  OUT) :: pd(VODE_NEQS,VODE_NEQS)
    type (dvode_t), intent(inout) :: vode_state

    type (burn_t) :: state
    integer :: n

    real(rt) :: small_temp

    !$gpu

    ! Only do the burn if the incoming temperature is within the temperature
    ! bounds. Otherwise set the Jacobian to zero and return.

    call eos_get_small_temp(small_temp)

    if (vode_state % y(net_itemp) <= small_temp .or. vode_state % y(net_itemp) >= MAX_TEMP) then

       pd(:,:) = 0.0_rt

       return

    end if

    ! Call the specific network routine to get the Jacobian.

    call vode_to_burn(vode_state, state)
    state % time = time
    call network_jac(state, pd, vode_state % rpar(irp_t0))

    ! We integrate X, not Y
    do n = 1, nspec
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

    call burn_to_vode(state, vode_state)

  end subroutine jac
end module vode_rhs_module
