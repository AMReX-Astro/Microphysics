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
    use burn_type_module, only: burn_t, net_ienuc, net_itemp, neqs
    use bs_type_module, only: bs_t, bs_to_burn, burn_to_bs
    use bs_rpar_indices, only: irp_y_init, irp_t0

    implicit none

    type (bs_t) :: bs
    real(rt) :: J(neqs, neqs)

    integer :: n

    ! Initialize the Jacacobian to zero.
    bs % jac(:,:) = ZERO

    ! Call the specific network routine to get the Jacobian.

    call bs_to_burn(bs)

    if (jacobian == 1) then

       call network_jac(bs % burn_s, bs % jac, bs % upar(irp_t0))

       ! We integrate X, not Y
       do n = 1, nspec_evolve
          bs % jac(n,:) = bs % jac(n,:) * aion(n)
          bs % jac(:,n) = bs % jac(:,n) * aion_inv(n)
       enddo

       ! Allow temperature and energy integration to be disabled.
       if (.not. integrate_temperature) then
          bs % jac(net_itemp,:) = ZERO
       endif

       if (.not. integrate_energy) then
          bs % jac(net_ienuc,:) = ZERO
       endif

    else

       call numerical_jac(bs % burn_s, bs % jac)

    endif

    ! apply fudge factor:
    if (react_boost > ZERO) then
       bs % jac(:,:) = react_boost * bs % jac(:,:)
    endif

    call burn_to_bs(bs)

    ! Increment the evaluation counter.

    bs % burn_s % n_jac = bs % burn_s % n_jac + 1

  end subroutine jac

end module jac_module
