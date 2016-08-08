module jac_module

  implicit none

contains

  subroutine jac(bs)

    !$acc routine seq

    use bl_constants_module, only: ZERO
    use actual_rhs_module, only: actual_jac
    use numerical_jac_module, only: numerical_jac
    use extern_probin_module, only: jacobian
    use burn_type_module, only: burn_t
    use bs_type_module, only: bs_t, bs_to_burn, jac_to_bs

    implicit none

    type (bs_t) :: bs
    type (burn_t) :: burn

    ! Initialize the Jacobian to zero.

    bs % jac = ZERO

    ! Call the specific network routine to get the Jacobian.

    if (jacobian == 1) then

       burn % jac = ZERO

       call bs_to_burn(bs, burn)

       call actual_jac(burn)

       call jac_to_bs(bs, burn)

    else

       call numerical_jac(bs)

    endif

    ! Increment the evaluation counter.

    bs % burn_s % n_jac = bs % burn_s % n_jac + 1

  end subroutine jac

end module jac_module
