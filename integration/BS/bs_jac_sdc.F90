module bs_jac_module

  implicit none

contains

  subroutine bs_jac(bs)

    !$acc routine seq

    use amrex_constants_module, only: ZERO
    use network_rhs_module, only: network_jac
    use numerical_jac_module, only: numerical_jac
    use extern_probin_module, only: jacobian
    use burn_type_module, only: burn_t
    use bs_type_module, only: bs_t, bs_to_burn, jac_to_bs
    use bs_rpar_indices, only: irp_t0

    implicit none

    type (bs_t) :: bs
    type (burn_t) :: burn

    ! Initialize the Jacobian to zero.

    bs % jac = ZERO

    ! Call the specific network routine to get the Jacobian.

    if (jacobian == 1) then

       burn % jac = ZERO

       call bs_to_burn(bs, burn)

       call network_jac(burn, bs % u(irp_t0))

       call jac_to_bs(bs, burn)

    else

       call numerical_jac(bs)

    endif

    ! Increment the evaluation counter.

    bs % n_jac = bs % n_jac + 1

  end subroutine bs_jac

end module bs_jac_module
