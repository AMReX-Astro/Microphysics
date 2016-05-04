module bs_convert_module

  use integration_data, only: temp_scale, dens_scale, ener_scale

  implicit none

  public

contains

  ! Given a bs_t, set up an EOS state.
  ! We could fill the energy component by storing the initial energy in
  ! ts % upar if we wanted, but we won't do that because (1) if we do call the EOS,
  ! it is always in (rho, T) mode and (2) converting back would imply subtracting
  ! off the nuclear energy from the zone's internal energy, which could lead to
  ! issues from roundoff error if the energy released from burning is small.

  subroutine bs_to_eos(state, bs)

    !$acc routine seq

    use network, only: nspec, aion
    use eos_type_module, only: eos_t
    use bs_type_module
    use burn_type_module
    use rpar_indices

    implicit none

    type (eos_t) :: state
    type (bs_t) :: bs

    state % rho     = bs % upar(irp_dens) * dens_scale
    state % T       = bs % y(net_itemp) * temp_scale
    state % xn(1:nspec_evolve) = bs % y(1:nspec_evolve) * aion(1:nspec_evolve)
    state % xn(nspec_evolve+1:nspec) = bs % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) * aion(nspec_evolve+1:nspec)
    state % cp      = bs % upar(irp_cp)
    state % cv      = bs % upar(irp_cv)
    state % abar    = bs % upar(irp_abar)
    state % zbar    = bs % upar(irp_zbar)
    state % eta     = bs % upar(irp_eta)
    state % y_e     = bs % upar(irp_ye)
    state % dhdX(:) = bs % upar(irp_dhdY:irp_dhdY-1+nspec) / aion(:)
    state % dedX(:) = bs % upar(irp_dedY:irp_dedY-1+nspec) / aion(:)

  end subroutine bs_to_eos



  ! Given an EOS state, fill the rpar and integration state data.

  subroutine eos_to_bs(state, bs)

    !$acc routine seq

    use network, only: nspec, aion
    use eos_type_module, only: eos_t
    use rpar_indices
    use integration_data, only: temp_scale, dens_scale
    use burn_type_module
    use bs_type_module

    implicit none

    type (eos_t)  :: state
    type (bs_t) :: bs

    bs % upar(irp_dens)                  = state % rho / dens_scale
    bs % y(net_itemp)                    = state % T / temp_scale
    bs % y(1:nspec_evolve)               = state % xn(1:nspec_evolve) / aion(1:nspec_evolve)
    bs % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) = state % xn(nspec_evolve+1:nspec) / aion(nspec_evolve+1:nspec)
    bs % upar(irp_cp)                    = state % cp
    bs % upar(irp_cv)                    = state % cv
    bs % upar(irp_abar)                  = state % abar
    bs % upar(irp_zbar)                  = state % zbar
    bs % upar(irp_eta)                   = state % eta
    bs % upar(irp_ye)                    = state % y_e
    bs % upar(irp_dhdY:irp_dhdY+nspec-1) = state % dhdX(:) * aion(:)
    bs % upar(irp_dedY:irp_dedY+nspec-1) = state % dedX(:) * aion(:)

  end subroutine eos_to_bs



  ! Given a burn state, fill the bs_t.

  subroutine burn_to_bs(state, bs)

    !$acc routine seq

    use network, only: nspec, aion
    use eos_type_module, only: eos_t
    use rpar_indices
    use burn_type_module
    use bl_constants_module
    use bs_type_module

    implicit none

    type (burn_t) :: state
    type (bs_t) :: bs

    integer :: i, n

    bs % upar(irp_dens)                           = state % rho / dens_scale
    bs % y(net_itemp)                             = state % T / temp_scale
    bs % y(1:nspec_evolve)                        = state % xn(1:nspec_evolve) / aion(1:nspec_evolve)
    bs % y(net_ienuc)                             = state % e / ener_scale
    bs % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) = state % xn(nspec_evolve+1:nspec) / aion(nspec_evolve+1:nspec)
    bs % upar(irp_cp)                             = state % cp
    bs % upar(irp_cv)                             = state % cv
    bs % upar(irp_abar)                           = state % abar
    bs % upar(irp_zbar)                           = state % zbar
    bs % upar(irp_ye)                             = state % y_e
    bs % upar(irp_eta)                            = state % eta
    bs % upar(irp_dhdY:irp_dhdY+nspec-1)          = state % dhdX(:) * aion(:)
    bs % upar(irp_dedY:irp_dedY+nspec-1)          = state % dedX(:) * aion(:)
    bs % upar(irp_Told)                           = state % T_old
    bs % upar(irp_dcvdt)                          = state % dcvdt
    bs % upar(irp_dcpdt)                          = state % dcpdt

    bs % dydt = state % ydot
    bs % dydt(net_itemp) = bs % dydt(net_itemp) / temp_scale
    bs % dydt(net_ienuc) = bs % dydt(net_ienuc) / ener_scale

    do n = 1, neqs
       bs % jac(n,:) = state % jac(n,:)
    enddo

    bs % jac(net_itemp,:) = bs % jac(net_itemp,:) / temp_scale
    bs % jac(net_ienuc,:) = bs % jac(net_ienuc,:) / ener_scale

    if (state % have_rates) then
       bs % upar(irp_have_rates) = ONE
    else
       bs % upar(irp_have_rates) = -ONE
    endif

    do i = 1, num_rate_groups
       bs % upar(irp_rates+(i-1)*nrates:irp_rates+i*nrates-1) = state % rates(i,:)
    enddo

    if (state % self_heat) then
       bs % upar(irp_self_heat) = ONE
    else
       bs % upar(irp_self_heat) = -ONE
    endif

  end subroutine burn_to_bs



  ! Given a bs_t, set up a burn state.

  subroutine bs_to_burn(bs, state)

    !$acc routine seq

    use network, only: nspec, aion
    use eos_type_module, only: eos_t
    use rpar_indices
    use burn_type_module
    use bl_constants_module
    use bs_type_module

    implicit none

    type (burn_t) :: state
    type (bs_t) :: bs

    integer :: i, n

    state % rho      = bs % upar(irp_dens) * dens_scale
    state % T        = bs % y(net_itemp) * temp_scale
    state % e        = bs % y(net_ienuc) * ener_scale
    state % xn(1:nspec_evolve) = bs % y(1:nspec_evolve) * aion(1:nspec_evolve)
    state % xn(nspec_evolve+1:nspec) = bs % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) * aion(nspec_evolve+1:nspec)
    state % cp       = bs % upar(irp_cp)
    state % cv       = bs % upar(irp_cv)
    state % abar     = bs % upar(irp_abar)
    state % zbar     = bs % upar(irp_zbar)
    state % y_e      = bs % upar(irp_ye)
    state % eta      = bs % upar(irp_eta)
    state % dhdX(:)  = bs % upar(irp_dhdY:irp_dhdY-1+nspec) / aion(:)
    state % dedX(:)  = bs % upar(irp_dedY:irp_dedY-1+nspec) / aion(:)
    state % T_old    = bs % upar(irp_Told)
    state % dcvdt    = bs % upar(irp_dcvdt)
    state % dcpdt    = bs % upar(irp_dcpdt)

    state % ydot = bs % dydt
    state % ydot(net_itemp) = state % ydot(net_itemp) * temp_scale
    state % ydot(net_ienuc) = state % ydot(net_ienuc) * ener_scale

    do n = 1, neqs
       state % jac(n,:) = bs % jac(n,:)
    enddo
    state % jac(net_itemp,:) = state % jac(net_itemp,:) * temp_scale
    state % jac(net_ienuc,:) = state % jac(net_ienuc,:) * ener_scale

    if (bs % upar(irp_have_rates) > ZERO) then
       state % have_rates = .true.
    else
       state % have_rates = .false.
    endif

    do i = 1, num_rate_groups
       state % rates(i,:) = bs % upar(irp_rates+(i-1)*nrates:irp_rates+i*nrates-1)
    enddo

    if (bs % upar(irp_self_heat) > ZERO) then
       state % self_heat = .true.
    else
       state % self_heat = .false.
    endif

  end subroutine bs_to_burn

end module bs_convert_module
