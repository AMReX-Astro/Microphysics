module actual_burner_module

  use eos_type_module

contains

  subroutine actual_burner_init()

    implicit none

    ! Do nothing in this burner.

  end subroutine actual_burner_init



  subroutine actual_burner(state_in, state_out, dt, time)

    implicit none

    type (eos_t), intent(in)    :: state_in
    type (eos_t), intent(inout) :: state_out
    double precision, intent(in)       :: dt, time

    ! Do nothing in this burner.

  end subroutine actual_burner

end module actual_burner_module
