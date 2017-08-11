module bs_rhs_module

contains

  ! The rhs routine provides the right-hand-side for the BS solver.
  ! This is a generic interface that calls the specific RHS routine in the
  ! network you're actually using.

  subroutine f_bs_rhs(bs)

    !$acc routine seq

    use actual_network, only: aion, nspec_evolve
    use bl_types, only: dp_t
    use burn_type_module, only: burn_t, net_ienuc, net_itemp
    use bl_constants_module, only: ZERO, ONE
    use actual_rhs_module, only: actual_rhs
    use extern_probin_module, only: burning_mode, burning_mode_factor, &
                                    integrate_temperature, integrate_energy
    use bs_type_module, only: bs_t, clean_state, renormalize_species, update_thermodynamics, &
                              burn_to_bs, bs_to_burn
    use rpar_indices, only: irp_y_init, irp_t_sound

    implicit none

    type (bs_t) :: bs

    real(dp_t) :: limit_factor, t_sound, t_enuc

    ! We are integrating a system of
    !
    ! y(1:nspec_evolve) = dX/dt
    ! y(net_itemp) = dT/dt
    ! y(net_ienuc) = denuc/dt

    ! Initialize the RHS to zero.

    bs % burn_s % ydot(:) = ZERO

    ! Fix the state as necessary.
    call clean_state(bs)

    ! Update the thermodynamic quantities as necessary.
    call update_thermodynamics(bs)

    ! Call the specific network routine to get the RHS.
    call bs_to_burn(bs)
    call actual_rhs(bs % burn_s)

    ! We integrate Y, not X
    bs % burn_s % ydot(1:nspec_evolve) = &
         bs % burn_s % ydot(1:nspec_evolve) * aion(1:nspec_evolve)

    ! Allow temperature and energy integration to be disabled.
    if (.not. integrate_temperature) then
       bs % burn_s % ydot(net_itemp) = ZERO
    endif

    if (.not. integrate_energy) then
       bs % burn_s % ydot(net_ienuc) = ZERO
    endif

    ! For burning_mode == 3, limit the rates.
    ! Note that we are limiting with respect to the initial zone energy.

    if (burning_mode == 3) then
       t_enuc = bs % upar(irp_y_init + net_ienuc - 1) / &
            max(abs(bs % burn_s % ydot(net_ienuc)), 1.d-50)
       t_sound = bs % upar(irp_t_sound)

       limit_factor = min(1.0d0, burning_mode_factor * t_enuc / t_sound)

       bs % burn_s % ydot = limit_factor * bs % burn_s % ydot

    endif

    call burn_to_bs(bs)

    ! Increment the evaluation counter.

    bs % burn_s % n_rhs = bs % burn_s % n_rhs + 1

  end subroutine f_bs_rhs

end module bs_rhs_module
