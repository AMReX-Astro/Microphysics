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



  ! Compute and store the more expensive screening factors  

  subroutine set_up_screening_factors()

    use screening_module, only: add_screening_factor
    use network, only: aion, zion

    implicit none

    ! note: it is critical that these are called in the exact order
    ! that the screening calls are done in the RHS routine, since we
    ! use that order in the screening

    call add_screening_factor(zion(ihe4),aion(ihe4),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(zion(ihe4),aion(ihe4),4.0d0,8.0d0)
    
    call add_screening_factor(zion(ic12),aion(ic12),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(zion(ic12),aion(ic12),zion(ic12),aion(ic12))
    
    call add_screening_factor(zion(ic12),aion(ic12),zion(io16),aion(io16))
    
    call add_screening_factor(zion(io16),aion(io16),zion(io16),aion(io16))
    
    call add_screening_factor(zion(io16),aion(io16),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(zion(ine20),aion(ine20),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(zion(img24),aion(img24),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(13.0d0,27.0d0,1.0d0,1.0d0)
    
    call add_screening_factor(zion(isi28),aion(isi28),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(15.0d0,31.0d0,1.0d0,1.0d0)
    
    call add_screening_factor(zion(is32),aion(is32),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(17.0d0,35.0d0,1.0d0,1.0d0)
    
    call add_screening_factor(zion(iar36),aion(iar36),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(19.0d0,39.0d0,1.0d0,1.0d0)
    
    call add_screening_factor(zion(ica40),aion(ica40),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(21.0d0,43.0d0,1.0d0,1.0d0)
    
    call add_screening_factor(zion(iti44),aion(iti44),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(23.0d0,47.0d0,1.0d0,1.0d0)
    
    call add_screening_factor(zion(icr48),aion(icr48),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(25.0d0,51.0d0,1.0d0,1.0d0)
    
    call add_screening_factor(zion(ife52),aion(ife52),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(27.0d0,55.0d0,1.0d0,1.0d0)
    
    call add_screening_factor(zion(ife54),aion(ife54),1.0d0,1.0d0)
    
    call add_screening_factor(zion(ife54),aion(ife54),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(zion(ife56),aion(ife56),1.0d0,1.0d0)
    
    call add_screening_factor(1.0d0,2.0d0,zion(ih1),aion(ih1))
    
    call add_screening_factor(zion(ih1),aion(ih1),zion(ih1),aion(ih1))
    
    call add_screening_factor(zion(ihe3),aion(ihe3),zion(ihe3),aion(ihe3))
    
    call add_screening_factor(zion(ihe3),aion(ihe3),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(zion(ic12),aion(ic12),zion(ih1),aion(ih1))
    
    call add_screening_factor(zion(in14),aion(in14),zion(ih1),aion(ih1))
    
    call add_screening_factor(zion(io16),aion(io16),zion(ih1),aion(ih1))
    
    call add_screening_factor(zion(in14),aion(in14),zion(ihe4),aion(ihe4))

  end subroutine set_up_screening_factors

end module actual_burner_module
