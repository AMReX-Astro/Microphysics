module actual_burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module
  use network
  use burn_type_module

  implicit none

contains

  subroutine actual_burner(state_in, state_out, dt, time)

    use integrator_module, only: integrator

    implicit none

    type (burn_t),       intent(in   ) :: state_in
    type (burn_t),       intent(inout) :: state_out
    double precision,    intent(in   ) :: dt, time

    call integrator(state_in, state_out, dt, time)

  end subroutine actual_burner



  subroutine actual_burner_init()

    use integrator_module, only: integrator_init
    use rates_module, only: rates_init
    use screening_module, only: screening_init
    use integration_data, only: ener_scale

    implicit none

    ener_scale = c_light * c_light

    call integrator_init()

    call rates_init()

    call set_up_screening_factors()

    call screening_init()

  end subroutine actual_burner_init


  subroutine set_up_screening_factors()
    ! Compute and store the more expensive screening factors  

    use screening_module, only: add_screening_factor
    use network, only: aion, zion

    implicit none

    ! note: we need to set these up in the same order that we evaluate the
    ! rates in actual_rhs.f90 (yes, it's ugly)
    call add_screening_factor(zion(ihe4),aion(ihe4),zion(ihe4),aion(ihe4))
    call add_screening_factor(zion(ihe4),aion(ihe4),4.0d0,8.0d0)
    call add_screening_factor(zion(ic12),aion(ic12),zion(ihe4),aion(ihe4))
    call add_screening_factor(zion(ic12),aion(ic12),zion(ic12),aion(ic12))
    call add_screening_factor(zion(ic12),aion(ic12),zion(io16),aion(io16))
    call add_screening_factor(zion(ic12),aion(ic12),zion(io16),aion(io16))
    call add_screening_factor(zion(io16),aion(io16),zion(ihe4),aion(ihe4))
    call add_screening_factor(zion(ine20),aion(ine20),zion(ihe4),aion(ihe4))
    call add_screening_factor(zion(img24),aion(img24),zion(ihe4),aion(ihe4))
    call add_screening_factor(20.0d0,40.0d0,zion(ihe4),aion(ihe4))

  end subroutine set_up_screening_factors

end module actual_burner_module
