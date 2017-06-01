module temperature_integration_module

  implicit none

  public

contains

  ! Sets up the temperature equation. This should be called from
  ! within the actual_rhs routine but is provided here as a convenience
  ! since most networks will use the same temperature ODE.

#ifdef CUDA  
  attributes(device) &
#endif
  subroutine temperature_rhs(state)

    !$acc routine seq

    use bl_constants_module, only: ZERO, ONE
    use network, only: nspec
    use burn_type_module
    use managed_probin_module, only: cu_do_constant_volume_burn, cu_dT_crit, cu_call_eos_in_rhs    
    use bl_types, only: dp_t

    implicit none

    type (burn_t) :: state
    real(dp_t) :: cv, cp, cvInv, cpInv

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

       if (cu_do_constant_volume_burn) then

          if (.not. cu_call_eos_in_rhs .and. cu_dT_crit < 1.0d19) then

             cv = state % cv + (state % T - state % T_old) * state % dcvdt

          else

             cv = state % cv

          endif

          cvInv = ONE / cv

          state % ydot(net_itemp) = state % ydot(net_ienuc) * cvInv

       else

          if (.not. cu_call_eos_in_rhs .and. cu_dT_crit < 1.0d19) then

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
#ifdef CUDA
  attributes(device) &
#endif
  subroutine temperature_jac(state)

    !$acc routine seq

    use bl_constants_module, only: ZERO, ONE
    use network, only: nspec
    use burn_type_module
    use managed_probin_module, only: cu_do_constant_volume_burn, cu_dT_crit, cu_call_eos_in_rhs
    use bl_types, only: dp_t

    implicit none

    type (burn_t) :: state

    real(dp_t) :: cp, cv, cpInv, cvInv

    ! Temperature Jacobian elements

    if (state % self_heat) then

       if (cu_do_constant_volume_burn) then

          if (.not. cu_call_eos_in_rhs .and. cu_dT_crit < 1.0d19) then

             cv = state % cv + (state % T - state % T_old) * state % dcvdt

          else

             cv = state % cv

          endif

          cvInv = ONE / cv

          ! d(itemp)/d(yi)

          state % jac(net_itemp,1:nspec_evolve) = state % jac(net_ienuc,1:nspec_evolve) * cvInv

          ! d(itemp)/d(temp)

          state % jac(net_itemp, net_itemp) = state % jac(net_ienuc,net_itemp) * cvInv

          ! d(itemp)/d(enuc)

          state % jac(net_itemp, net_ienuc) = ZERO

       else

          if (.not. cu_call_eos_in_rhs .and. cu_dT_crit < 1.0d19) then

             cp = state % cp + (state % T - state % T_old) * state % dcpdt

          else

             cp = state % cp

          endif

          cpInv = ONE / cp

          ! d(itemp)/d(yi)

          state % jac(net_itemp,1:nspec_evolve) = state % jac(net_ienuc,1:nspec_evolve) * cpInv

          ! d(itemp)/d(temp)

          state % jac(net_itemp,net_itemp) = state % jac(net_ienuc,net_itemp) * cpInv

          ! d(itemp)/d(enuc)

          state % jac(net_itemp, net_ienuc) = ZERO

       endif

    endif

  end subroutine temperature_jac

end module temperature_integration_module
