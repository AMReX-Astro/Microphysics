module integrator_module

  implicit none

  public
  
contains

  subroutine integrator_init()

    use actual_integrator_module, only: actual_integrator_init

    implicit none

    call actual_integrator_init()

  end subroutine integrator_init

#ifdef CUDA
  attributes(device) &
#endif       
  subroutine integrator(state_in, state_out, dt, time)

    !$acc routine seq

    use actual_integrator_module, only: actual_integrator
    use burn_type_module, only: burn_t
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    type (burn_t),  intent(in   ) :: state_in
    type (burn_t),  intent(inout) :: state_out
    real(rt),     intent(in   ) :: dt, time

    call actual_integrator(state_in, state_out, dt, time)

  end subroutine integrator

end module integrator_module
