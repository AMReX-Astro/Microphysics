module rhs_module

contains

  ! The rhs routine provides the right-hand-side for the BS solver.
  ! This is a generic interface that calls the specific RHS routine in the
  ! network you're actually using.

  subroutine f_rhs(bs)

    !$acc routine seq

    use burn_type_module, only: burn_t
    use bl_constants_module, only: ZERO
    use actual_rhs_module, only: actual_rhs
    use extern_probin_module, only: renormalize_abundances
    use bs_type_module, only: bs_t, clean_state, renormalize_species, rhs_to_bs, bs_to_burn

    implicit none

    type (bs_t) :: bs
    type (burn_t) :: burn

    ! Initialize the RHS to zero.

    burn % ydot = ZERO
    bs % ydot = ZERO

    ! Fix the state as necessary.

    call clean_state(bs)

    ! Renormalize abundances as necessary.

    if (renormalize_abundances) then
       call renormalize_species(bs)
    endif

    ! Tell the burn_t what it needs to know about
    ! the current integration state. This also does
    ! an EOS call to fill the data needed for the RHS.

    call bs_to_burn(bs, burn)

    ! Indicate that we don't yet have valid rates.

    burn % have_rates = .false.

    ! Call the specific network routine to get its RHS.

    call actual_rhs(burn)

    ! Feed the network evaluation into the integration.

    call rhs_to_bs(bs, burn)

    ! Increment the evaluation counter.

    bs % n_rhs = bs % n_rhs + 1

  end subroutine f_rhs



  ! Analytical Jacobian

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

       ! Indicate that we don't yet have valid rates.

       burn % have_rates = .false.

       call bs_to_burn(bs, burn)

       call actual_jac(burn)

       call jac_to_bs(bs, burn)

    else

       call numerical_jac(bs)

    endif

    ! Increment the evaluation counter.

    bs % n_jac = bs % n_jac + 1

  end subroutine jac

end module rhs_module
