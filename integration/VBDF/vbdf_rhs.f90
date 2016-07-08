module rhs_module

contains

  ! The rhs routine provides the right-hand-side for the VBDF solver.
  ! This is a generic interface that calls the specific RHS routine in the
  ! network you're actually using.

  subroutine rhs(ts)

    !$acc routine seq

    use actual_network, only: aion, nspec_evolve
    use bl_types, only: dp_t
    use burn_type_module, only: burn_t, net_ienuc, net_itemp
    use bl_constants_module, only: ZERO, ONE
    use actual_rhs_module, only: actual_rhs
    use extern_probin_module, only: renormalize_abundances, burning_mode, burning_mode_factor, &
                                    integrate_temperature, integrate_energy, integrate_molar_fraction
    use bdf_type_module, only: bdf_ts, clean_state, renormalize_species, update_thermodynamics, &
                               burn_to_vbdf, vbdf_to_burn
    use rpar_indices, only: irp_have_rates, irp_y_init, irp_t_sound

    implicit none

    type (bdf_ts) :: ts

    type (burn_t) :: burn_state

    real(dp_t) :: limit_factor, t_sound, t_enuc

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

    ! Indicate that we don't yet have valid rates.

    ts % upar(irp_have_rates,1) = -ONE

    ! Call the specific network routine to get the RHS.

    call vbdf_to_burn(ts, burn_state)
    call actual_rhs(burn_state)

    ! Allow integration of X instead of Y.

    if (.not. integrate_molar_fraction) then

       burn_state % ydot(1:nspec_evolve) = burn_state % ydot(1:nspec_evolve) * aion(1:nspec_evolve)

    endif

    ! Allow temperature and energy integration to be disabled.

    if (.not. integrate_temperature) then

       burn_state % ydot(net_itemp) = ZERO

    endif

    if (.not. integrate_energy) then

       burn_state % ydot(net_ienuc) = ZERO

    endif

    ! For burning_mode == 3, limit the rates.
    ! Note that we are limiting with respect to the initial zone energy.

    if (burning_mode == 3) then

       t_enuc = ts % upar(irp_y_init + net_ienuc - 1,1) / max(abs(burn_state % ydot(net_ienuc)), 1.d-50)
       t_sound = ts % upar(irp_t_sound,1)

       limit_factor = min(1.0d0, burning_mode_factor * t_enuc / t_sound)

       burn_state % ydot = limit_factor * burn_state % ydot

    endif

    call burn_to_vbdf(burn_state, ts)

  end subroutine rhs



  ! Analytical Jacobian

  subroutine jac(ts)

    !$acc routine seq

    use actual_network, only: aion, nspec_evolve
    use integration_data, only: aionInv
    use bl_types, only: dp_t
    use bl_constants_module, only: ZERO, ONE
    use actual_rhs_module, only: actual_jac
    use numerical_jac_module, only: numerical_jac
    use extern_probin_module, only: jacobian, burning_mode, burning_mode_factor, &
                                    integrate_temperature, integrate_energy, integrate_molar_fraction
    use burn_type_module, only: burn_t, net_ienuc, net_itemp
    use bdf_type_module, only: bdf_ts, vbdf_to_burn, burn_to_vbdf
    use rpar_indices, only: irp_have_rates, irp_y_init, irp_t_sound

    implicit none

    type (bdf_ts) :: ts

    type (burn_t) :: state

    real(dp_t) :: limit_factor, t_sound, t_enuc

    integer :: n

    ! Indicate that we don't yet have valid rates.

    ts % upar(irp_have_rates,1) = -ONE

    ! Initialize the Jacobian to zero.

    ts % J(:,:,1) = ZERO

    ! Call the specific network routine to get the Jacobian.

    call vbdf_to_burn(ts, state)

    if (jacobian == 1) then

       call actual_jac(state)

       ! Allow integration of X instead of Y.

       if (.not. integrate_molar_fraction) then

          do n = 1, nspec_evolve
             state % jac(n,:) = state % jac(n,:) * aion(n)
             state % jac(:,n) = state % jac(:,n) * aionInv(n)
          enddo

       endif

       ! Allow temperature and energy integration to be disabled.

       if (.not. integrate_temperature) then

          state % jac(net_itemp,:) = ZERO

       endif

       if (.not. integrate_energy) then

          state % jac(net_ienuc,:) = ZERO

       endif

       ! For burning_mode == 3, limit the rates.
       ! Note that we are limiting with respect to the initial zone energy.

       if (burning_mode == 3) then

          t_enuc = ts % upar(irp_y_init + net_ienuc - 1, 1) / max(abs(state % ydot(net_ienuc)), 1.d-50)
          t_sound = ts % upar(irp_t_sound,1)

          limit_factor = min(1.0d0, burning_mode_factor * t_enuc / t_sound)

          state % jac = limit_factor * state % jac

       endif

    else

       call numerical_jac(state)

    endif

    call burn_to_vbdf(state, ts)

  end subroutine jac

end module rhs_module
