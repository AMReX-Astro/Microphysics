module rhs_module

    use amrex_fort_module, only : rt => amrex_real

contains

  ! The rhs routine provides the right-hand-side for the VBDF solver.
  ! This is a generic interface that calls the specific RHS routine in the
  ! network you're actually using.

  subroutine rhs(ts)

    !$acc routine seq

    use actual_network, only: aion, nspec_evolve
    use burn_type_module, only: burn_t, net_ienuc, net_itemp
    use amrex_constants_module, only: ZERO, ONE
    use network_rhs_module, only: network_rhs
    use extern_probin_module, only: renormalize_abundances, integrate_temperature, integrate_energy, react_boost
    use bdf_type_module, only: bdf_ts, clean_state, renormalize_species, update_thermodynamics, &
                               burn_to_vbdf, vbdf_to_burn
    use vbdf_rpar_indices, only: irp_y_init, irp_t_sound, irp_t0

    implicit none

    type (bdf_ts) :: ts

    type (burn_t) :: burn_state

    real(rt) :: limit_factor, t_sound, t_enuc

    ! We are integrating a system of
    !
    ! y(1:nspec_evolve) = dX/dt
    ! y(net_itemp) = dT/dt
    ! y(net_ienuc) = denuc/dt

    ! Initialize the RHS to zero.

    ts % yd(:,1) = ZERO

    ! Fix the state as necessary.

    call clean_state(ts)

    ! Renormalize the abundances as necessary.

    if (renormalize_abundances) then
       call renormalize_species(ts)
    endif

    ! Update the thermodynamics as necessary.

    call update_thermodynamics(ts)

    ! Call the specific network routine to get the RHS.

    call vbdf_to_burn(ts, burn_state)
    call network_rhs(burn_state, ts % yd(:,1), ts % upar(irp_t0,1))

    ! We integrate X not Y, so convert here
    ts % yd(1:nspec_evolve,1) = ts % yd(1:nspec_evolve,1) * aion(1:nspec_evolve)

    ! Allow temperature and energy integration to be disabled.

    if (.not. integrate_temperature) then

       ts % yd(net_itemp,1) = ZERO

    endif

    if (.not. integrate_energy) then

       ts % yd(net_ienuc,1) = ZERO

    endif

    ! apply fudge factor:
    if (react_boost > ZERO) then
       ts % yd(:,1) = react_boost * ts % yd(:,1)
    endif

    call burn_to_vbdf(burn_state, ts)

  end subroutine rhs



  ! Analytical Jacobian

  subroutine jac(ts)

    !$acc routine seq

    use network, only: aion, aion_inv, nspec_evolve
    use amrex_constants_module, only: ZERO, ONE
    use network_rhs_module, only: network_jac
    use numerical_jac_module, only: numerical_jac
    use extern_probin_module, only: jacobian, integrate_temperature, integrate_energy, react_boost
    use burn_type_module, only: burn_t, net_ienuc, net_itemp
    use bdf_type_module, only: bdf_ts, vbdf_to_burn, burn_to_vbdf
    use vbdf_rpar_indices, only: irp_y_init, irp_t_sound, irp_t0

    implicit none

    type (bdf_ts) :: ts

    type (burn_t) :: state

    real(rt) :: limit_factor, t_sound, t_enuc

    integer :: n

    ! Initialize the Jacobian to zero.

    ts % J(:,:,1) = ZERO

    ! Call the specific network routine to get the Jacobian.

    call vbdf_to_burn(ts, state)

    if (jacobian == 1) then

       call network_jac(state, ts % J(:,:,1), ts % upar(irp_t0,1))

       ! We integrate X, not Y
       do n = 1, nspec_evolve
          ts % J(n,:,1) = ts % J(n,:,1) * aion(n)
          ts % J(:,n,1) = ts % J(:,n,1) * aion_inv(n)
       enddo

       ! Allow temperature and energy integration to be disabled.
       if (.not. integrate_temperature) then
          ts % J(net_itemp,:,1) = ZERO
       endif

       if (.not. integrate_energy) then
          ts % J(net_ienuc,:,1) = ZERO
       endif

    else

       call numerical_jac(state, ts % J(:,:,1))

    endif

    ! apply fudge factor:
    if (react_boost > ZERO) then
       ts % J(:,:,1) = react_boost * ts % J(:,:,1)
    endif

    call burn_to_vbdf(state, ts)

  end subroutine jac

end module rhs_module
