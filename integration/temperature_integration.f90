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
    use extern_probin_module, only: do_constant_volume_burn, use_chemical_potential

    implicit none

    type (burn_t)    :: state

    double precision :: dedt, dhdt
    double precision :: cvInv, cpInv

    if (state % self_heat) then

       ! Set up the temperature ODE.  For constant pressure, Dp/Dt = 0, we
       ! evolve :
       !    dT/dt = (1/c_p) [ -sum_i (xi_i omega_i) + Hnuc]
       !
       ! For constant volume, div{U} = 0, and we evolve:
       !    dT/dt = (1/c_v) [ -sum_i ( {e_x}_i omega_i) + Hnuc]
       !
       ! See paper III, including Eq. A3 for details.

       if (do_constant_volume_burn) then

          cvInv = ONE / state % cv

          state % ydot(net_itemp) = state % ydot(net_ienuc) * cvInv

          if (use_chemical_potential .eq. 1) then

             dedt = sum( state % dedX(1:nspec_evolve) * aion(1:nspec_evolve) * state % ydot(1:nspec_evolve) )

             state % ydot(net_itemp) = state % ydot(net_itemp) - dedt * cvInv

          endif

       else

          cpInv = ONE / state % cp

          state % ydot(net_itemp) = state % ydot(net_ienuc) * cpInv

          if (use_chemical_potential .eq. 1) then

             dhdt = sum( state % dhdX(1:nspec_evolve) * aion(1:nspec_evolve) * state % ydot(1:nspec_evolve) )

             state % ydot(net_itemp) = state % ydot(net_itemp) - dhdt * cpInv

          endif

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
    use extern_probin_module, only: do_constant_volume_burn, use_chemical_potential

    implicit none

    type (burn_t)    :: state

    integer          :: j
    double precision :: dedt, dhdt
    double precision :: cpInv, cvInv

    ! Temperature Jacobian elements

    if (state % self_heat) then

       if (do_constant_volume_burn) then

          cvInv = ONE / state % cv

          ! d(itemp)/d(yi)

          state % jac(net_itemp,1:nspec_evolve) = state % jac(net_ienuc,1:nspec_evolve) * cvInv

          ! d(itemp)/d(temp)

          state % jac(net_itemp, net_itemp) = state % jac(net_ienuc,net_itemp) * cvInv

          if (use_chemical_potential .eq. 1) then

             do j = 1, nspec_evolve
                dedt = sum( state % dedX(1:nspec_evolve) * aion(1:nspec_evolve) * state % jac(1:nspec_evolve,j) )
                state % jac(net_itemp,j) = state % jac(net_itemp,j) - dedt * cvInv
             enddo

             dedt = sum( state % dedX(1:nspec_evolve) * aion(1:nspec_evolve) * state % jac(1:nspec_evolve,net_itemp) )
             state % jac(net_itemp,net_itemp) = state % jac(net_itemp,net_itemp) - dedt * cvInv

          endif

       else

          cpInv = ONE / state % cp

          ! d(itemp)/d(yi)

          state % jac(net_itemp,1:nspec_evolve) = state % jac(net_ienuc,1:nspec_evolve) * cpInv

          ! d(itemp)/d(temp)

          state % jac(net_itemp,net_itemp) = state % jac(net_ienuc,net_itemp) * cpInv

          if (use_chemical_potential .eq. 1) then

             do j = 1, nspec_evolve
                dhdt = sum( state % dhdX(1:nspec_evolve) * aion(1:nspec_evolve) * state % jac(1:nspec_evolve,j) )
                state % jac(net_itemp,j) = state % jac(net_itemp,j) - dhdt * cpInv
             enddo

             dhdt = sum ( state % dhdX(1:nspec_evolve) * aion(1:nspec_evolve) * state % jac(1:nspec_evolve,net_itemp) )
             state % jac(net_itemp,net_itemp) = state % jac(net_itemp,net_itemp) - dhdt * cpInv

          endif

       endif

    endif

  end subroutine temperature_jac

end module temperature_integration_module
