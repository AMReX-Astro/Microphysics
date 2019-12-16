module actual_rhs_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

  subroutine actual_rhs_init()

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    ! Do nothing in this RHS.

  end subroutine actual_rhs_init



  subroutine actual_rhs(state)

    use burn_type_module, only: burn_t
    use amrex_constants_module, only: ZERO

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (burn_t) :: state

    ! Do nothing in this RHS.

    state % ydot = ZERO

  end subroutine actual_rhs



  subroutine actual_jac(state)

    use burn_type_module, only: burn_t
    use amrex_constants_module, only: ZERO

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (burn_t) :: state

    ! Do nothing in this RHS.

    state % jac(:,:) = ZERO

  end subroutine actual_jac



  subroutine update_unevolved_species(state)

    use burn_type_module, only: burn_t

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (burn_t)    :: state

  end subroutine update_unevolved_species

end module actual_rhs_module
