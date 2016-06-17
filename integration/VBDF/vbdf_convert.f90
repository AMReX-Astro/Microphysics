module vbdf_convert_module

  use integration_data, only: temp_scale, dens_scale, ener_scale, inv_temp_scale, inv_dens_scale, inv_ener_scale, aionInv

  implicit none

  public

contains

  ! Given a bdf_ts, set up an EOS state.
  ! We could fill the energy component by storing the initial energy in
  ! ts % upar if we wanted, but we won't do that because (1) if we do call the EOS,
  ! it is always in (rho, T) mode and (2) converting back would imply subtracting
  ! off the nuclear energy from the zone's internal energy, which could lead to
  ! issues from roundoff error if the energy released from burning is small.

  subroutine vbdf_to_eos(state, ts)

    !$acc routine seq

    use network, only: nspec, aion
    use eos_type_module, only: eos_t
    use bdf_type_module
    use burn_type_module
    use rpar_indices

    implicit none

    type (eos_t)  :: state
    type (bdf_ts) :: ts

    state % rho     = ts % upar(irp_dens,1) * dens_scale
    state % T       = ts % y(net_itemp,1) * temp_scale
    state % xn(1:nspec_evolve) = ts % y(1:nspec_evolve,1) * aion(1:nspec_evolve)
    state % xn(nspec_evolve+1:nspec) = ts % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1,1) * aion(nspec_evolve+1:nspec)
    state % cp      = ts % upar(irp_cp,1)
    state % cv      = ts % upar(irp_cv,1)
    state % abar    = ts % upar(irp_abar,1)
    state % zbar    = ts % upar(irp_zbar,1)
    state % eta     = ts % upar(irp_eta,1)
    state % y_e     = ts % upar(irp_ye,1)
    state % dhdX(:) = ts % upar(irp_dhdY:irp_dhdY-1+nspec,1) * aionInv(:)
    state % dedX(:) = ts % upar(irp_dedY:irp_dedY-1+nspec,1) * aionInv(:)

  end subroutine vbdf_to_eos



  ! Given an EOS state, fill the rpar and integration state data.

  subroutine eos_to_vbdf(state, ts)

    !$acc routine seq

    use network, only: nspec, aion
    use eos_type_module, only: eos_t
    use rpar_indices
    use integration_data, only: temp_scale, dens_scale
    use burn_type_module
    use bdf_type_module

    implicit none

    type (eos_t)  :: state
    type (bdf_ts) :: ts

    ts % upar(irp_dens,1)                  = state % rho * inv_dens_scale
    ts % y(net_itemp,1)                    = state % T * inv_temp_scale
    ts % y(1:nspec_evolve,1)               = state % xn(1:nspec_evolve) * aionInv(1:nspec_evolve)
    ts % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1,1) = state % xn(nspec_evolve+1:nspec) * aionInv(nspec_evolve+1:nspec)
    ts % upar(irp_cp,1)                    = state % cp
    ts % upar(irp_cv,1)                    = state % cv
    ts % upar(irp_abar,1)                  = state % abar
    ts % upar(irp_zbar,1)                  = state % zbar
    ts % upar(irp_eta,1)                   = state % eta
    ts % upar(irp_ye,1)                    = state % y_e
    ts % upar(irp_dhdY:irp_dhdY+nspec-1,1) = state % dhdX(:) * aion(:)
    ts % upar(irp_dedY:irp_dedY+nspec-1,1) = state % dedX(:) * aion(:)

  end subroutine eos_to_vbdf



  ! Given a burn state, fill the bdf_ts.

  subroutine burn_to_vbdf(state, ts)

    !$acc routine seq

    use network, only: nspec, aion
    use eos_type_module, only: eos_t
    use rpar_indices
    use burn_type_module
    use bl_constants_module
    use bdf_type_module
    use integration_data

    implicit none

    type (burn_t) :: state
    type (bdf_ts) :: ts

    integer :: n

    ts % upar(irp_dens,1)                           = state % rho * inv_dens_scale
    ts % y(net_itemp,1)                             = state % T * inv_temp_scale
    ts % y(1:nspec_evolve,1)                        = state % xn(1:nspec_evolve) * aionInv(1:nspec_evolve)
    ts % y(net_ienuc,1)                             = state % e * inv_ener_scale
    ts % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1,1) = state % xn(nspec_evolve+1:nspec) * aionInv(nspec_evolve+1:nspec)
    ts % upar(irp_cp,1)                             = state % cp
    ts % upar(irp_cv,1)                             = state % cv
    ts % upar(irp_abar,1)                           = state % abar
    ts % upar(irp_zbar,1)                           = state % zbar
    ts % upar(irp_ye,1)                             = state % y_e
    ts % upar(irp_eta,1)                            = state % eta
    ts % upar(irp_dhdY:irp_dhdY+nspec-1,1)          = state % dhdX(:) * aion(:)
    ts % upar(irp_dedY:irp_dedY+nspec-1,1)          = state % dedX(:) * aion(:)
    ts % upar(irp_Told,1)                           = state % T_old
    ts % upar(irp_dcvdt,1)                          = state % dcvdt
    ts % upar(irp_dcpdt,1)                          = state % dcpdt

    ts % yd(:,1) = state % ydot
    ts % yd(net_itemp,1) = ts % yd(net_itemp,1) * inv_temp_scale
    ts % yd(net_ienuc,1) = ts % yd(net_ienuc,1) * inv_ener_scale

    ts % J(:,:,1) = state % jac
    ts % J(net_itemp,:,1) = ts % J(net_itemp,:,1) * inv_temp_scale
    ts % J(net_ienuc,:,1) = ts % J(net_ienuc,:,1) * inv_ener_scale

    if (state % have_rates) then
       ts % upar(irp_have_rates,1) = ONE
    else
       ts % upar(irp_have_rates,1) = -ONE
    endif

    do n = 1, nrates
       ts % upar(irp_rates+(n-1)*num_rate_groups:irp_rates+n*num_rate_groups-1,1) = state % rates(:,n)
    enddo

    if (state % self_heat) then
       ts % upar(irp_self_heat,1) = ONE
    else
       ts % upar(irp_self_heat,1) = -ONE
    endif

  end subroutine burn_to_vbdf


  ! Given a bdf_ts, set up a burn state.

  subroutine vbdf_to_burn(ts, state)

    !$acc routine seq

    use network, only: nspec, aion
    use eos_type_module, only: eos_t
    use rpar_indices
    use burn_type_module
    use bl_constants_module
    use bdf_type_module
    use integration_data

    implicit none

    type (burn_t) :: state
    type (bdf_ts) :: ts

    integer :: n

    state % rho      = ts % upar(irp_dens,1) * dens_scale
    state % T        = ts % y(net_itemp,1) * temp_scale
    state % e        = ts % y(net_ienuc,1) * ener_scale
    state % xn(1:nspec_evolve) = ts % y(1:nspec_evolve,1) * aion(1:nspec_evolve)
    state % xn(nspec_evolve+1:nspec) = ts % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1,1) * aion(nspec_evolve+1:nspec)
    state % cp       = ts % upar(irp_cp,1)
    state % cv       = ts % upar(irp_cv,1)
    state % abar     = ts % upar(irp_abar,1)
    state % zbar     = ts % upar(irp_zbar,1)
    state % y_e      = ts % upar(irp_ye,1)
    state % eta      = ts % upar(irp_eta,1)
    state % dhdX(:)  = ts % upar(irp_dhdY:irp_dhdY-1+nspec,1) * aionInv(:)
    state % dedX(:)  = ts % upar(irp_dedY:irp_dedY-1+nspec,1) * aionInv(:)
    state % T_old    = ts % upar(irp_Told,1)
    state % dcvdt    = ts % upar(irp_dcvdt,1)
    state % dcpdt    = ts % upar(irp_dcpdt,1)

    state % ydot = ts % yd(:,1)
    state % ydot(net_itemp) = state % ydot(net_itemp) * inv_temp_scale
    state % ydot(net_ienuc) = state % ydot(net_ienuc) * inv_ener_scale

    state % jac = ts % J(:,:,1)
    state % jac(net_itemp,:) = state % jac(net_itemp,:) * inv_temp_scale
    state % jac(net_ienuc,:) = state % jac(net_ienuc,:) * inv_ener_scale

    if (ts % upar(irp_have_rates,1) > ZERO) then
       state % have_rates = .true.
    else
       state % have_rates = .false.
    endif

    do n = 1, nrates
       state % rates(:,n) = ts % upar(irp_rates+(n-1)*num_rate_groups:irp_rates+n*num_rate_groups-1,1)
    enddo

    if (ts % upar(irp_self_heat,1) > ZERO) then
       state % self_heat = .true.
    else
       state % self_heat = .false.
    endif

  end subroutine vbdf_to_burn

end module vbdf_convert_module
