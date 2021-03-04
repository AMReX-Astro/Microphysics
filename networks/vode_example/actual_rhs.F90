module actual_rhs_module

  use amrex_constants_module
  use network
  use burn_type_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

  subroutine actual_rhs_init()

  end subroutine actual_rhs_init


  subroutine actual_rhs(state, ydot)

    implicit none

    type (burn_t), intent(in)    :: state
    real(rt)        , intent(inout) :: ydot(neqs)

  end subroutine actual_rhs


  subroutine actual_jac(state, jac)

    implicit none

    type (burn_t), intent(in)    :: state
    real(rt)        , intent(inout) :: jac(njrows, njcols)

  end subroutine actual_jac

end module actual_rhs_module
