module actual_burner_module

  use amrex_constants_module
  use network
  use eos_type_module
  use burn_type_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

  subroutine actual_burner_init()

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    ! Do nothing in this burner.

  end subroutine actual_burner_init

  subroutine actual_burner(state_in, state_out, dt, time)

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    type (burn_t),    intent(in)    :: state_in
    type (burn_t),    intent(inout) :: state_out
    real(rt)        , intent(in)    :: dt, time

    ! Do nothing in this burner.
    
  end subroutine actual_burner

end module actual_burner_module
