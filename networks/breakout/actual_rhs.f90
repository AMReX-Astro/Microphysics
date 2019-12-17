module actual_rhs_module

  use burn_type_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

  subroutine actual_rhs_init()

    implicit none

    ! Do nothing in this RHS.

  end subroutine actual_rhs_init



  subroutine actual_rhs(state, ydot)

    implicit none

    type (burn_t), intent(in) :: state
    real(rt)        , intent(inout) :: ydot(neqs)

    ! Do nothing in this RHS.

    ydot = ZERO

  end subroutine actual_rhs



  subroutine actual_jac(state, jac)

    implicit none

    type (burn_t), intent(in) :: state
    real(rt)        , intent(inout) :: jac(njrows, njcols)

    ! Do nothing in this RHS.

    state % jac(:,:) = ZERO

  end subroutine actual_jac

  subroutine update_unevolved_species(state)

    implicit none

    type (burn_t)    :: state

  end subroutine update_unevolved_species

end module actual_rhs_module
