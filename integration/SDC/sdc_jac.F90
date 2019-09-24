module jac_module

  implicit none

contains

  subroutine jac(bs)

    !$acc routine seq

    use network, only: aion, aion_inv, nspec_evolve
    use amrex_constants_module, only: ZERO, ONE
    use amrex_fort_module, only : rt => amrex_real
    use network_rhs_module, only: network_jac
    use numerical_jac_module, only: numerical_jac
    use extern_probin_module, only: jacobian, integrate_temperature, integrate_energy, react_boost
    use burn_type_module, only: burn_t, net_ienuc, net_itemp
    use sdc_type_module, only: sdc_t, sdc_to_burn, burn_to_sdc
    use sdc__rpar_indices, only: irp_y_init

    implicit none

    type (bs_t) :: bs

    integer :: n

    ! Initialize the Jacobian to zero.
    bs % burn_s % jac(:,:) = ZERO

    ! Call the specific network routine to get the Jacobian.

    call bs_to_burn(bs)

    if (jacobian == 1) then

       call network_jac(bs % burn_s)

       ! We integrate X, not Y
       do n = 1, nspec_evolve
          bs % burn_s % jac(n,:) = bs % burn_s % jac(n,:) * aion(n)
          bs % burn_s % jac(:,n) = bs % burn_s % jac(:,n) * aion_inv(n)
       enddo

       ! Allow temperature and energy integration to be disabled.
       if (.not. integrate_temperature) then
          bs % burn_s % jac(net_itemp,:) = ZERO
       endif

       if (.not. integrate_energy) then
          bs % burn_s % jac(net_ienuc,:) = ZERO
       endif

    else

       call numerical_jac(bs % burn_s)

    endif

    ! apply fudge factor:
    if (react_boost > ZERO) then
       bs % burn_s % jac(:,:) = react_boost * bs % burn_s % jac(:,:)
    endif

    call burn_to_bs(bs)

    ! Increment the evaluation counter.

    bs % burn_s % n_jac = bs % burn_s % n_jac + 1

  end subroutine jac

end module jac_module
