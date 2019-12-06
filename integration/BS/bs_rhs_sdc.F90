module bs_rhs_module

contains

  ! The rhs routine provides the right-hand-side for the BS solver.
  ! This is a generic interface that calls the specific RHS routine in the
  ! network you're actually using.

  subroutine f_bs_rhs(bs)

    !$acc routine seq

    use burn_type_module, only: burn_t
    use network_rhs_module, only: network_rhs
    use bs_type_module, only: bs_t, clean_state, rhs_to_bs, bs_to_burn
    use bs_rpar_indices, only: irp_t0

    implicit none

    type (bs_t) :: bs
    type (burn_t) :: burn

    ! Initialize the RHS to zero.

    burn % ydot = ZERO
    bs % ydot = ZERO

    ! Fix the state as necessary.

    call clean_state(bs)

    ! Tell the burn_t what it needs to know about
    ! the current integration state. This also does
    ! an EOS call to fill the data needed for the RHS.

    call bs_to_burn(bs, burn)

    ! Call the specific network routine to get its RHS.

    call network_rhs(burn, bs % u(irp_t0))

    ! Feed the network evaluation into the integration.

    call rhs_to_bs(bs, burn)

    ! Increment the evaluation counter.

    bs % n_rhs = bs % n_rhs + 1

  end subroutine f_bs_rhs

end module bs_rhs_module
