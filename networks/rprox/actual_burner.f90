! This module contains a version of Wallace and Woosley's (ApJS 45,389
! (1981)) rprox reaction network burner.

module actual_burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module
  use eos_type_module
  use network
  use burn_type_module

contains

  subroutine actual_burner_init()

    use integrator_module, only: integrator_init

    implicit none

    call integrator_init()

  end subroutine actual_burner_init



  subroutine actual_burner(state_in, state_out, dt, time)  

    use integrator_module, only: integrator

    implicit none

    type (burn_t),    intent(in   ) :: state_in
    type (burn_t),    intent(inout) :: state_out
    double precision, intent(in   ) :: dt, time    

    double precision :: T9

    T9 = state_in % T / 10.0**9

    ! Only burn if 0.2 < T9 < 2.5 or X(H1) > 0.05.
    ! The last restriction is a kludge based on the last paragraph of WW81.

    if ((T9 .gt. 0.2d0 .and. T9 .lt. 2.5d0) .or. state_in % xn(ih1) > 0.05d0) then

       ! Call the integration routine.

       call integrator(state_in, state_out, dt, time)

    endif
    
  end subroutine actual_burner

end module actual_burner_module
