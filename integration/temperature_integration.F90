module temperature_integration_module

  implicit none

  logical, save, allocatable :: self_heat

#ifdef AMREX_USE_CUDA
  attributes(managed) :: self_heat
#endif

  public

contains

  ! Sets up the temperature equation. This should be called from
  ! within the actual_rhs routine but is provided here as a convenience
  ! since most networks will use the same temperature ODE.

  subroutine temperature_rhs(state)

    !$acc routine seq

    use amrex_constants_module, only: ZERO, ONE
    use amrex_fort_module, only : rt => amrex_real
    use network, only: nspec
    use burn_type_module
    use jacobian_sparsity_module, only: get_jac_entry, set_jac_entry
    use extern_probin_module, only: do_constant_volume_burn, dT_crit, call_eos_in_rhs

    implicit none

    type (burn_t) :: state
    real(rt) :: cv, cp, cvInv, cpInv

    !$gpu

    if (state % self_heat) then

       ! Set up the temperature ODE.  For constant pressure, Dp/Dt = 0, we
       ! evolve :
       !    dT/dt = (1/c_p) [ Hnuc ]
       !
       ! For constant volume, div{U} = 0, and we evolve:
       !    dT/dt = (1/c_v) [ Hnuc ]
       !
       ! See low Mach paper III, including Eq. A3 for details.
       ! Note that we no longer include the chemical potential (dE/dX or dH/dX)
       ! terms because we believe they analytically should vanish.

       if (do_constant_volume_burn) then

          if (.not. call_eos_in_rhs .and. dT_crit < 1.0d19) then

             cv = state % cv + (state % T - state % T_old) * state % dcvdt

          else

             cv = state % cv

          endif

          cvInv = ONE / cv

          state % ydot(net_itemp) = state % ydot(net_ienuc) * cvInv

       else

          if (.not. call_eos_in_rhs .and. dT_crit < 1.0d19) then

             cp = state % cp + (state % T - state % T_old) * state % dcpdt

          else

             cp = state % cp

          endif

          cpInv = ONE / cp

          state % ydot(net_itemp) = state % ydot(net_ienuc) * cpInv

       endif

    endif

  end subroutine temperature_rhs



  ! Sets up the temperature entries in the Jacobian. This should be called from
  ! within the actual_jac routine but is provided here as a convenience
  ! since most networks will use the same temperature ODE.

  subroutine temperature_jac(state)

    !$acc routine seq

    use amrex_constants_module, only: ZERO, ONE
    use amrex_fort_module, only : rt => amrex_real
    use network, only: nspec
    use burn_type_module
    use jacobian_sparsity_module, only: get_jac_entry, set_jac_entry
    use extern_probin_module, only: do_constant_volume_burn, dT_crit, call_eos_in_rhs

    implicit none

    type (burn_t) :: state

    real(rt) :: scratch, cspec, cspecInv

    integer :: k

    !$gpu

    ! Temperature Jacobian elements

    if (state % self_heat) then

       if (do_constant_volume_burn) then

          if (.not. call_eos_in_rhs .and. dT_crit < 1.0d19) then

             cspec = state % cv + (state % T - state % T_old) * state % dcvdt

          else

             cspec = state % cv

          endif

       else

          if (.not. call_eos_in_rhs .and. dT_crit < 1.0d19) then

             cspec = state % cp + (state % T - state % T_old) * state % dcpdt

          else

             cspec = state % cp

          endif

       endif

       cspecInv = ONE / cspec       

       ! d(itemp)/d(yi)
       
       do k = 1, nspec_evolve
          call get_jac_entry(state, net_ienuc, k, scratch)
          scratch = scratch * cspecInv
          call set_jac_entry(state, net_itemp, k, scratch)
       enddo

       ! d(itemp)/d(temp)

       call get_jac_entry(state, net_ienuc, net_itemp, scratch)
       scratch = scratch * cspecInv
       call set_jac_entry(state, net_itemp, net_itemp, scratch)

       ! d(itemp)/d(enuc)

       scratch = ZERO
       call set_jac_entry(state, net_itemp, net_ienuc, scratch)

    endif

  end subroutine temperature_jac



  subroutine temperature_rhs_init()

    use extern_probin_module, only: burning_mode
    use amrex_error_module, only: amrex_error

    implicit none

    ! Provide a default value, then consult the burning_mode.

    allocate(self_heat)
    self_heat = .true.

    if (burning_mode == 0 .or. burning_mode == 2) then
       self_heat = .false.
    else if (burning_mode == 1 .or. burning_mode == 3) then
       self_heat = .true.
    else
       call amrex_error("Error: unknown burning_mode in temperature_rhs_init()")
    end if

  end subroutine temperature_rhs_init

end module temperature_integration_module
