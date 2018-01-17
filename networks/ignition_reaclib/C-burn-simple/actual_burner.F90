module actual_burner_module

  use network

  implicit none

contains

  subroutine actual_burner_init()
    use reaclib_rates, only: init_reaclib, net_screening_init
    use integrator_module, only: integrator_init

    implicit none

    call integrator_init()

    call init_reaclib()
    call net_screening_init()    
  end subroutine actual_burner_init

  subroutine actual_burner_finalize
    use reaclib_rates, only: term_reaclib, net_screening_finalize

    implicit none
    
    call term_reaclib()
    call net_screening_finalize()    
  end subroutine actual_burner_finalize

#ifdef CUDA
  attributes(device) &
#endif       
  subroutine actual_burner(state_in, state_out, dt, time)

    !$acc routine seq

    use integrator_module, only: integrator
    use burn_type_module, only: burn_t
    use bl_types, only: dp_t

    implicit none

    type (burn_t),    intent(in   ) :: state_in
    type (burn_t),    intent(inout) :: state_out
    real(kind=dp_t),  intent(in   ) :: dt, time

    call integrator(state_in, state_out, dt, time)
  end subroutine actual_burner

end module actual_burner_module
