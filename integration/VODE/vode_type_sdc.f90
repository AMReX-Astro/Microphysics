module vode_type_module

  use bl_types, only: dp_t
  use rpar_indices, only: n_rpar_comps
  use sdc_type_module, only: SVAR, SVAR_EVOLVE

  implicit none

  private

  integer, parameter :: VODE_NEQS = SVAR_EVOLVE

  public :: VODE_NEQS

contains

  subroutine clean_state(y, rpar)

    real(dp_t) :: y(SVAR_EVOLVE), rpar(n_rpar_comps)

    real (kind=dp_t) :: max_e, ke

    type (eos_t) :: eos_state


    ! Ensure that mass fractions always stay positive.
    y(SFS:SFS+nspec-1) = &
         max(min(y(SFS:SFS+nspec-1), rpar(irp_SRHO)), &
             rpar(irp_SRHO) * 1.d-200)

    ! renormalize abundances as necessary
    if (renormalize_abundances) then
       call renormalize_species(y, rpar)
    endif

    ! Ensure that internal energy never goes above the maximum limit
    ! provided by the EOS. Same for the internal energy implied by the
    ! total energy (which we get by subtracting kinetic energy).
    eos_state % rho = rpar(irp_SRHO) 
    eos_state % T = MAX_TEMP
    eos_state % xn = y(SFS:SFS+nspec-1) / rpar(irp_SRHO)

    call eos(eos_input_rt, eos_state)

    max_e = eos_state % e

    y(SEINT) = min(rpar(irp_SRHO) * max_e, y(SEINT))

    ke = y(SEDEN) - HALF * sum(rpar(irp_SMX:irp_SMZ)**2) / rpar(irp_SRHO)

    y(SEDEN) = min(rpar(irp_SRHO) * max_e + ke, y(SEDEN))

  end subroutine clean_state


  subroutine fill_unevolved_variables(time, y, rpar)

    real(dp_t), intent(in) :: time
    real(dp_t) :: y(SVAR_EVOLVE), rpar(n_rpar_comps)
    
    ! we are always integrating from t = 0, so there is no offset
    ! time needed here
    rpar(irp_SRHO) = rpar(irp_u_init-1+irp_SRHO) + &
         rpar(irp_ydot_a-1+irp_SRHO) * time

    do m = irp_SMX, irp_SMZ
       rpar(m) = rpar(irp_u_init-1+m) + rpar(irp_ydot_a-1+m) * time
    enddo

  end subroutine fill_unevolved_variables

  subroutine renormalize_species(time, y, rpar)

    real(dp_t), intent(in) :: time
    real(dp_t) :: y(SVAR_EVOLVE), rpar(n_rpar_comps)

    real(dp_t) :: nspec_sum

    ! update rho, rho*u, etc.
    call fill_unevolved_variables(time, y, rpar)

    nspec_sum = sum(y(SFS:SFS-1+nspec)) / rpar(irp_SRHO)

    y(SFS:SFS-1+nspec) = y(SFS:SFS-1+nspec) / nspec_sum

  end subroutine renormalize_species


  ! Given a burn state, fill the rpar and integration state data.

  subroutine burn_to_vode(state, y, rpar, ydot, jac)

    use bl_types, only: dp_t
    use bl_constants_module, only: ONE
    use actual_network, only: nspec, nspec_evolve, aion
    use integration_data, only: aionInv
    use rpar_indices, only: irp_dens, irp_nspec, irp_cp, irp_cv, irp_abar, irp_zbar, &
                            irp_ye, irp_eta, irp_cs, irp_dx, &
                            irp_Told, irp_dcvdt, irp_dcpdt, irp_self_heat, &
                            n_rpar_comps, n_not_evolved
    use burn_type_module, only: SVAR_EVOLVE, burn_t, net_itemp, net_ienuc
    use extern_probin_module, only: integrate_molar_fraction

    implicit none

    type (burn_t) :: state
    real(dp_t)    :: rpar(n_rpar_comps)
    real(dp_t)    :: y(SVAR_EVOLVE)
    real(dp_t), optional :: ydot(SVAR_EVOLVE), jac(SVAR_EVOLVE, SVAR_EVOLVE)

    integer :: n

    rpar(irp_dens) = state % rho
    y(net_itemp) = state % T

    if (integrate_molar_fraction) then
       y(1:nspec_evolve) = state % xn(1:nspec_evolve) * aionInv(1:nspec_evolve)
       rpar(irp_nspec:irp_nspec+n_not_evolved-1) = &
            state % xn(nspec_evolve+1:nspec) * aionInv(nspec_evolve+1:nspec)
    else
       y(1:nspec_evolve) = state % xn(1:nspec_evolve)
       rpar(irp_nspec:irp_nspec+n_not_evolved-1) = &
            state % xn(nspec_evolve+1:nspec)
    endif

    y(net_ienuc)                             = state % e
    rpar(irp_cp)                             = state % cp
    rpar(irp_cv)                             = state % cv
    rpar(irp_abar)                           = state % abar
    rpar(irp_zbar)                           = state % zbar
    rpar(irp_ye)                             = state % y_e
    rpar(irp_eta)                            = state % eta
    rpar(irp_cs)                             = state % cs
    rpar(irp_dx)                             = state % dx

    rpar(irp_Told)                           = state % T_old
    rpar(irp_dcvdt)                          = state % dcvdt
    rpar(irp_dcpdt)                          = state % dcpdt

    if (present(ydot)) then
       ydot = state % ydot
    endif

    if (present(jac)) then
       jac = state % jac
    endif

    if (state % self_heat) then
       rpar(irp_self_heat) = ONE
    else
       rpar(irp_self_heat) = -ONE
    endif

  end subroutine burn_to_vode



  ! Given an rpar array and the integration state, set up a burn state.

  subroutine vode_to_burn(y, rpar, state)

    use bl_types, only: dp_t
    use bl_constants_module, only: ZERO
    use actual_network, only: nspec, nspec_evolve, aion
    use integration_data, only: aionInv
    use rpar_indices, only: irp_dens, irp_nspec, irp_cp, irp_cv, irp_abar, irp_zbar, &
                            irp_ye, irp_eta, irp_cs, irp_dx, &
                            irp_Told, irp_dcvdt, irp_dcpdt, irp_self_heat, &
                            n_rpar_comps, n_not_evolved
    use burn_type_module, only: SVAR_EVOLVE, burn_t, net_itemp, net_ienuc
    use extern_probin_module, only: integrate_molar_fraction

    implicit none

    type (burn_t) :: state
    real(dp_t)    :: rpar(n_rpar_comps)
    real(dp_t)    :: y(SVAR_EVOLVE)

    integer :: n

    state % rho      = rpar(irp_dens)
    state % T        = y(net_itemp)
    state % e        = y(net_ienuc)

    if (integrate_molar_fraction) then
       state % xn(1:nspec_evolve) = y(1:nspec_evolve) * aion(1:nspec_evolve)
       state % xn(nspec_evolve+1:nspec) = &
            rpar(irp_nspec:irp_nspec+n_not_evolved-1) * aion(nspec_evolve+1:nspec)
    else
       state % xn(1:nspec_evolve) = y(1:nspec_evolve)
       state % xn(nspec_evolve+1:nspec) = &
            rpar(irp_nspec:irp_nspec+n_not_evolved-1)
    endif

    state % cp       = rpar(irp_cp)
    state % cv       = rpar(irp_cv)
    state % abar     = rpar(irp_abar)
    state % zbar     = rpar(irp_zbar)
    state % y_e      = rpar(irp_ye)
    state % eta      = rpar(irp_eta)
    state % cs       = rpar(irp_cs)
    state % dx       = rpar(irp_dx)

    state % T_old    = rpar(irp_Told)
    state % dcvdt    = rpar(irp_dcvdt)
    state % dcpdt    = rpar(irp_dcpdt)

    if (rpar(irp_self_heat) > ZERO) then
       state % self_heat = .true.
    else
       state % self_heat = .false.
    endif

  end subroutine vode_to_burn

end module vode_type_module
