module actual_burner_module

  use bl_types
  use bl_constants_module
  use network
  use eos_type_module
  use actual_burner_data
  
contains

  subroutine actual_burner_init()

    ! Do nothing in this burner.

  end subroutine actual_burner_init



  subroutine actual_burner(state_in, state_out, dt, time)

    implicit none

    type (eos_t),     intent(in)    :: state_in
    type (eos_t),     intent(inout) :: state_out
    double precision, intent(in)    :: dt, time

    ! Do nothing in this burner.
    
  end subroutine actual_burner

end module actual_burner_module
