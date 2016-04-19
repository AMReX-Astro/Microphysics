module temperature_integration_module

  implicit none

  public

contains

  ! Sets up the temperature equation. This should be called from
  ! within the actual_rhs routine but is provided here as a convenience
  ! since most networks will use the same temperature ODE.

  subroutine temperature_rhs(state)

    use bl_constants_module, only: ZERO, ONE
    use network, only: nspec, aion
    use burn_type_module
    use extern_probin_module, only: do_constant_volume_burn

    implicit none

    type (burn_t) :: state

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

          state % ydot(net_itemp) = state % ydot(net_ienuc) / state % cv

       else

          state % ydot(net_itemp) = state % ydot(net_ienuc) / state % cp

       endif

    endif

  end subroutine temperature_rhs



  ! Sets up the temperature entries in the Jacobian. This should be called from
  ! within the actual_jac routine but is provided here as a convenience
  ! since most networks will use the same temperature ODE.

  subroutine temperature_jac(state)

    use bl_constants_module, only: ZERO, ONE
    use network, only: nspec, aion
    use burn_type_module
    use extern_probin_module, only: do_constant_volume_burn

    implicit none

    type (burn_t)    :: state

    double precision :: cpInv, cvInv

    ! Temperature Jacobian elements

    if (state % self_heat) then

       if (do_constant_volume_burn) then

          cvInv = ONE / state % cv

          ! d(itemp)/d(yi)

          state % jac(net_itemp,1:nspec_evolve) = state % jac(net_ienuc,1:nspec_evolve) * cvInv

          ! d(itemp)/d(temp)

          state % jac(net_itemp, net_itemp) = state % jac(net_ienuc,net_itemp) * cvInv

       else

          cpInv = ONE / state % cp

          ! d(itemp)/d(yi)

          state % jac(net_itemp,1:nspec_evolve) = state % jac(net_ienuc,1:nspec_evolve) * cpInv

          ! d(itemp)/d(temp)

          state % jac(net_itemp,net_itemp) = state % jac(net_ienuc,net_itemp) * cpInv

       endif

    endif

  end subroutine temperature_jac

end module temperature_integration_module
