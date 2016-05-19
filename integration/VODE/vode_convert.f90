module vode_convert_module

  use integration_data, only: temp_scale, dens_scale
  use bl_types, only: dp_t

  implicit none

  public

contains

  ! Given an rpar array and the integration state, set up an EOS state.
  ! We could fill the energy component by storing the initial energy in
  ! rpar if we wanted, but we won't do that because (1) if we do call the EOS,
  ! it is always in (rho, T) mode and (2) converting back would imply subtracting
  ! off the nuclear energy from the zone's internal energy, which could lead to
  ! issues from roundoff error if the energy released from burning is small.

  subroutine vode_to_eos(state, y, rpar)

    use network, only: nspec, aion
    use eos_type_module, only: eos_t
    use rpar_indices
    use burn_type_module

    implicit none

    type (eos_t) :: state
    real(dp_t)   :: rpar(n_rpar_comps)
    real(dp_t)   :: y(neqs)

    state % rho     = rpar(irp_dens) * dens_scale
    state % T       = y(net_itemp) * temp_scale
    state % xn(1:nspec_evolve) = y(1:nspec_evolve) * aion(1:nspec_evolve)
    state % xn(nspec_evolve+1:nspec) = rpar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) * aion(nspec_evolve+1:nspec)
    state % cp      = rpar(irp_cp)
    state % cv      = rpar(irp_cv)
    state % abar    = rpar(irp_abar)
    state % zbar    = rpar(irp_zbar)
    state % eta     = rpar(irp_eta)
    state % y_e     = rpar(irp_ye)
    state % dhdX(:) = rpar(irp_dhdY:irp_dhdY-1+nspec) / aion(:)
    state % dedX(:) = rpar(irp_dedY:irp_dedY-1+nspec) / aion(:)

  end subroutine vode_to_eos



  ! Given an EOS state, fill the rpar and integration state data.

  subroutine eos_to_vode(state, y, rpar)

    use network, only: nspec, aion
    use eos_type_module, only: eos_t
    use rpar_indices
    use integration_data, only: temp_scale, dens_scale
    use burn_type_module

    implicit none

    type (eos_t) :: state
    real(dp_t)   :: rpar(n_rpar_comps)
    real(dp_t)   :: y(neqs)

    rpar(irp_dens)                  = state % rho / dens_scale
    y(net_itemp)                    = state % T / temp_scale
    y(1:nspec_evolve)               = state % xn(1:nspec_evolve) / aion(1:nspec_evolve)
    rpar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) = state % xn(nspec_evolve+1:nspec) / aion(nspec_evolve+1:nspec)
    rpar(irp_cp)                    = state % cp
    rpar(irp_cv)                    = state % cv
    rpar(irp_abar)                  = state % abar
    rpar(irp_zbar)                  = state % zbar
    rpar(irp_eta)                   = state % eta
    rpar(irp_ye)                    = state % y_e
    rpar(irp_dhdY:irp_dhdY+nspec-1) = state % dhdX(:) * aion(:)
    rpar(irp_dedY:irp_dedY+nspec-1) = state % dedX(:) * aion(:)

  end subroutine eos_to_vode



  ! Given a burn state, fill the rpar and integration state data.

  subroutine burn_to_vode(state, y, rpar, ydot, jac)

    use network, only: nspec, aion
    use eos_type_module, only: eos_t
    use rpar_indices
    use burn_type_module
    use bl_constants_module
    use integration_data

    implicit none

    type (burn_t) :: state
    real(dp_t)    :: rpar(n_rpar_comps)
    real(dp_t)    :: y(neqs)
    real(dp_t), optional :: ydot(neqs), jac(neqs, neqs)

    integer       :: i

    rpar(irp_dens)                           = state % rho / dens_scale
    y(net_itemp)                             = state % T / temp_scale
    y(1:nspec_evolve)                        = state % xn(1:nspec_evolve) / aion(1:nspec_evolve)
    rpar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) = state % xn(nspec_evolve+1:nspec) / aion(nspec_evolve+1:nspec)
    y(net_ienuc)                             = state % e / ener_scale
    rpar(irp_cp)                             = state % cp
    rpar(irp_cv)                             = state % cv
    rpar(irp_abar)                           = state % abar
    rpar(irp_zbar)                           = state % zbar
    rpar(irp_ye)                             = state % y_e
    rpar(irp_eta)                            = state % eta
    rpar(irp_dhdY:irp_dhdY+nspec-1)          = state % dhdX(:) * aion(:)
    rpar(irp_dedY:irp_dedY+nspec-1)          = state % dedX(:) * aion(:)
    rpar(irp_Told)                           = state % T_old
    rpar(irp_dcvdt)                          = state % dcvdt
    rpar(irp_dcpdt)                          = state % dcpdt

    if (present(ydot)) then
       ydot = state % ydot
       ydot(net_itemp) = ydot(net_itemp) / temp_scale
       ydot(net_ienuc) = ydot(net_ienuc) / ener_scale
    endif

    if (present(jac)) then
       jac = state % jac
       jac(net_itemp,:) = jac(net_itemp,:) / temp_scale
       jac(net_ienuc,:) = jac(net_ienuc,:) / ener_scale
    endif

    if (state % have_rates) then
       rpar(irp_have_rates) = ONE
    else
       rpar(irp_have_rates) = -ONE
    endif

    do i = 1, num_rate_groups
       rpar(irp_rates+(i-1)*nrates:irp_rates+i*nrates-1) = state % rates(i,:)
    enddo

    if (state % self_heat) then
       rpar(irp_self_heat) = ONE
    else
       rpar(irp_self_heat) = -ONE
    endif

  end subroutine burn_to_vode


  ! Given an rpar array and the integration state, set up a burn state.

  subroutine vode_to_burn(y, rpar, state)

    use network, only: nspec, aion
    use eos_type_module, only: eos_t
    use rpar_indices
    use burn_type_module
    use bl_constants_module
    use integration_data

    implicit none

    type (burn_t) :: state
    real(dp_t)    :: rpar(n_rpar_comps)
    real(dp_t)    :: y(neqs)

    integer       :: i

    state % rho      = rpar(irp_dens) * dens_scale
    state % T        = y(net_itemp) * temp_scale
    state % xn(1:nspec_evolve) = y(1:nspec_evolve) * aion(1:nspec_evolve)
    state % xn(nspec_evolve+1:nspec) = rpar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) * aion(nspec_evolve+1:nspec)
    state % e        = y(net_ienuc) * ener_scale
    state % cp       = rpar(irp_cp)
    state % cv       = rpar(irp_cv)
    state % abar     = rpar(irp_abar)
    state % zbar     = rpar(irp_zbar)
    state % y_e      = rpar(irp_ye)
    state % eta      = rpar(irp_eta)
    state % dhdX(:)  = rpar(irp_dhdY:irp_dhdY-1+nspec) / aion(:)
    state % dedX(:)  = rpar(irp_dedY:irp_dedY-1+nspec) / aion(:)
    state % T_old    = rpar(irp_Told)
    state % dcvdt    = rpar(irp_dcvdt)
    state % dcpdt    = rpar(irp_dcpdt)

    if (rpar(irp_have_rates) > ZERO) then
       state % have_rates = .true.
    else
       state % have_rates = .false.
    endif

    do i = 1, num_rate_groups
       state % rates(i,:) = rpar(irp_rates+(i-1)*nrates:irp_rates+i*nrates-1)
    enddo

    if (rpar(irp_self_heat) > ZERO) then
       state % self_heat = .true.
    else
       state % self_heat = .false.
    endif

  end subroutine vode_to_burn

end module vode_convert_module
