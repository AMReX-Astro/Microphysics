module rhs_module

contains

  ! The rhs routine provides the right-hand-side for the SDC solver.
  ! This is a generic interface that calls the specific RHS routine in the
  ! network you're actually using.

  subroutine f_rhs(time, y, ydot, rpar)

    use actual_network, only: aion, nspec_evolve
    use amrex_fort_module, only : rt => amrex_real
    use burn_type_module, only: burn_t, net_ienuc, net_itemp
    use amrex_constants_module, only: ZERO, ONE
    use network_rhs_module, only: network_rhs
    use extern_probin_module, only: integrate_temperature, integrate_energy, react_boost
    use sdc_type_module, only: clean_state, renormalize_species, update_thermodynamics, &
                              burn_to_sdc, sdc_to_burn
    use sdc_rpar_indices, only: n_rpar_comps, irp_y_init, irp_t0
    use sdc_parameters_module, only : SDC_NEQS
    implicit none

    real(rt), intent(in) :: time
    real(rt), intent(INOUT) :: y(SDC_NEQS)
    real(rt), intent(INOUT) :: rpar(n_rpar_comps)
    real(rt), intent(INOUT) :: ydot(SDC_NEQS)

    type (burn_t) :: burn_state

    ! We are integrating a system of
    !
    ! y(1:nspec_evolve) = dX/dt
    ! y(net_itemp) = dT/dt
    ! y(net_ienuc) = denuc/dt

    ! Initialize the RHS to zero.

    ydot(:) = ZERO

    ! Fix the state as necessary.
    call clean_state(y, rpar)

    ! Update the thermodynamic quantities as necessary.
    call update_thermodynamics(y, rpar)

    ! Call the specific network routine to get the RHS.
    call sdc_to_burn(y, rpar, burn_state)

    burn_state % time = rpar(irp_t0) + time
    call network_rhs(burn_state)

    ! We integrate X, not Y
    burn_state % ydot(1:nspec_evolve) = &
         burn_state % ydot(1:nspec_evolve) * aion(1:nspec_evolve)

    ! Allow temperature and energy integration to be disabled.
    if (.not. integrate_temperature) then
       burn_state % ydot(net_itemp) = ZERO
    endif

    if (.not. integrate_energy) then
       burn_state % ydot(net_ienuc) = ZERO
    endif

    ! apply fudge factor:
    if (react_boost > ZERO) then
       burn_state % ydot(:) = react_boost * burn_state % ydot(:)
    endif

    call burn_to_sdc(burn_state, y, rpar, ydot = ydot)

  end subroutine f_rhs

  subroutine jac(time, y, ml, mu, pd, nrpd, rpar)

    use network, only: aion, aion_inv, nspec_evolve
    use amrex_constants_module, only: ZERO, ONE
    use amrex_fort_module, only : rt => amrex_real
    use network_rhs_module, only: network_jac
    use numerical_jac_module, only: numerical_jac
    use extern_probin_module, only: jacobian, integrate_temperature, integrate_energy, react_boost
    use burn_type_module, only: burn_t, net_ienuc, net_itemp
    use sdc_type_module, only: sdc_to_burn, burn_to_sdc
    use sdc_rpar_indices, only: n_rpar_comps, irp_y_init, irp_t0
    use sdc_parameters_module, only : SDC_NEQS

    implicit none

    integer   , intent(IN   ) :: ml, mu, nrpd
    real(rt), intent(in) :: time
    real(rt), intent(INOUT) :: y(SDC_NEQS), rpar(n_rpar_comps)
    real(rt), intent(  OUT) :: pd(SDC_NEQS,SDC_NEQS)

    type (burn_t) :: state
    integer :: n

    ! Initialize the Jacobian to zero.
    pd(:,:) = ZERO

    ! Call the specific network routine to get the Jacobian.

    call sdc_to_burn(y, rpar, state)
    state % time = rpar(irp_t0) + time

    if (jacobian == 1) then

       call network_jac(state)

       ! We integrate X, not Y
       do n = 1, nspec_evolve
          state % jac(n,:) = state % jac(n,:) * aion(n)
          state % jac(:,n) = state % jac(:,n) * aion_inv(n)
       enddo

       ! Allow temperature and energy integration to be disabled.
       if (.not. integrate_temperature) then
          state % jac(net_itemp,:) = ZERO
       endif

       if (.not. integrate_energy) then
          state % jac(net_ienuc,:) = ZERO
       endif

    else

       call numerical_jac(state)

    endif

    ! apply fudge factor:
    if (react_boost > ZERO) then
       state % jac(:,:) = react_boost * state % jac(:,:)
    endif

    call burn_to_sdc(state, y, rpar, jac = pd)

  end subroutine jac

end module rhs_module
