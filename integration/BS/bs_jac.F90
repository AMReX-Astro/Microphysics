module jac_module

  implicit none

contains

  subroutine jac(bs)

    !$acc routine seq

    use network, only: aion, aion_inv, nspec_evolve
    use bl_types, only: dp_t
    use bl_constants_module, only: ZERO, ONE
    use actual_rhs_module, only: actual_jac
    use numerical_jac_module, only: numerical_jac
    use extern_probin_module, only: jacobian, burning_mode, burning_mode_factor, &
                                    integrate_temperature, integrate_energy
    use burn_type_module, only: burn_t, net_ienuc, net_itemp
    use bs_type_module, only: bs_t, bs_to_burn, burn_to_bs
    use rpar_indices, only: irp_y_init, irp_t_sound

    implicit none

    type (bs_t) :: bs

    real(dp_t) :: limit_factor, t_sound, t_enuc

    integer :: n

    ! Initialize the Jacobian to zero.
    bs % burn_s % jac(:,:) = ZERO

    ! Call the specific network routine to get the Jacobian.

    call bs_to_burn(bs)

    if (jacobian == 1) then

       call actual_jac(bs % burn_s)

       ! We integrate X, not Y
       do n = 1, nspec_evolve
          bs % burn_s % jac(n,:) = bs % burn_s % jac(n,:) * aion(n)
          bs % burn_s % jac(:,n) = bs % burn_s % jac(:,n) * aionInv(n)
       enddo

       ! Allow temperature and energy integration to be disabled.
       if (.not. integrate_temperature) then
          bs % burn_s % jac(net_itemp,:) = ZERO
       endif

       if (.not. integrate_energy) then
          bs % burn_s % jac(net_ienuc,:) = ZERO
       endif

       ! For burning_mode == 3, limit the rates.
       ! Note that we are limiting with respect to the initial zone energy.

       if (burning_mode == 3) then
          t_enuc = bs % upar(irp_y_init + net_ienuc - 1) / &
               max(abs(bs % burn_s % ydot(net_ienuc)), 1.d-50)
          t_sound = bs % upar(irp_t_sound)

          limit_factor = min(1.0d0, burning_mode_factor * t_enuc / t_sound)

          bs % burn_s % jac = limit_factor * bs % burn_s % jac

       endif

    else

       call numerical_jac(bs % burn_s)

    endif

    call burn_to_bs(bs)

    ! Increment the evaluation counter.

    bs % burn_s % n_jac = bs % burn_s % n_jac + 1

  end subroutine jac

end module jac_module
