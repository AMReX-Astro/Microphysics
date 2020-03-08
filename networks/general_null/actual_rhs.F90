module actual_rhs_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

  subroutine actual_rhs_init()

    implicit none

    ! Do nothing in this RHS.

  end subroutine actual_rhs_init



  subroutine actual_rhs(state, ydot)

    use burn_type_module, only: burn_t, neqs
    use amrex_constants_module, only: ZERO

    implicit none

    type (burn_t) :: state
    real(rt), intent(inout) :: ydot(neqs)

    !$gpu

    ! Do nothing in this RHS.

    ydot = ZERO

  end subroutine actual_rhs



  subroutine actual_jac(state, jac)

    use burn_type_module, only: burn_t, njrows, njcols
    use amrex_constants_module, only: ZERO

    implicit none

    type (burn_t) :: state
    real(rt), intent(inout) :: jac(njrows, njcols)

    !$gpu

    ! Do nothing in this RHS.

    jac(:,:) = ZERO

  end subroutine actual_jac



  subroutine update_unevolved_species(state)

    use burn_type_module, only: burn_t

    implicit none

    !$gpu

    type (burn_t)    :: state

  end subroutine update_unevolved_species

end module actual_rhs_module
