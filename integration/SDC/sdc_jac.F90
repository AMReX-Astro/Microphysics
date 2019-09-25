module jac_module

  implicit none

contains

  subroutine jac(sdc)

    !$acc routine seq

    use network, only: aion, aion_inv, nspec_evolve
    use amrex_constants_module, only: ZERO, ONE
    use amrex_fort_module, only : rt => amrex_real
    use network_rhs_module, only: network_jac
    use numerical_jac_module, only: numerical_jac
    use extern_probin_module, only: jacobian, integrate_temperature, integrate_energy, react_boost
    use burn_type_module, only: burn_t, net_ienuc, net_itemp
    use sdc_type_module, only: sdc_t, sdc_to_burn, burn_to_sdc
    use sdc_rpar_indices, only: irp_y_init

    implicit none

    type (sdc_t) :: sdc

    integer :: n

    ! Initialize the Jacobian to zero.
    sdc % burn_s % jac(:,:) = ZERO

    ! Call the specific network routine to get the Jacobian.

    call sdc_to_burn(sdc)

    if (jacobian == 1) then

       call network_jac(sdc % burn_s)

       ! We integrate X, not Y
       do n = 1, nspec_evolve
          sdc % burn_s % jac(n,:) = sdc % burn_s % jac(n,:) * aion(n)
          sdc % burn_s % jac(:,n) = sdc % burn_s % jac(:,n) * aion_inv(n)
       enddo

       ! Allow temperature and energy integration to be disabled.
       if (.not. integrate_temperature) then
          sdc % burn_s % jac(net_itemp,:) = ZERO
       endif

       if (.not. integrate_energy) then
          sdc % burn_s % jac(net_ienuc,:) = ZERO
       endif

    else

       call numerical_jac(sdc % burn_s)

    endif

    ! apply fudge factor:
    if (react_boost > ZERO) then
       sdc % burn_s % jac(:,:) = react_boost * sdc % burn_s % jac(:,:)
    endif

    call burn_to_sdc(sdc)

    ! Increment the evaluation counter.

    sdc % burn_s % n_jac = sdc % burn_s % n_jac + 1

  end subroutine jac

end module jac_module
