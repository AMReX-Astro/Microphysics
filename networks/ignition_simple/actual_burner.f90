module actual_burner_module

  use bl_types
  use bl_constants_module
  use network
  use burn_type_module
  use actual_burner_data

  implicit none

  ! Conversion factor for the nuclear energy generation rate.

  double precision, parameter, private :: avo = 6.0221417930d23
  double precision, parameter, private :: c_light = 2.99792458d10
  double precision, parameter, private :: enuc_conv2 = -avo*c_light*c_light

contains

  subroutine actual_burner_init()

    use integration_module, only: integration_init

    implicit none

    call integration_init()

  end subroutine actual_burner_init



  subroutine actual_burner(state_in, state_out, dt, time)

    use integration_module, only: do_burn

    implicit none

    type (burn_t),    intent(in   ) :: state_in
    type (burn_t),    intent(inout) :: state_out
    double precision, intent(in   ) :: dt, time

    call do_burn(state_in, state_out, dt, time)

  end subroutine actual_burner



  ! Computes the instantaneous energy generation rate

  subroutine ener_gener_rate(dydt, enuc)

    implicit none

    double precision :: dydt(nspec), enuc

    ! This is basically e = m c**2

    enuc = sum(dydt(:) * mion(:)) * enuc_conv2

  end subroutine ener_gener_rate

end module actual_burner_module
