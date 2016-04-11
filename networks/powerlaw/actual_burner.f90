module actual_burner_module

  use bl_types
  use bl_constants_module
  use network
  use burn_type_module
  use actual_burner_data
  use extern_probin_module, only: specific_q_burn

  implicit none

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

    double precision :: dydt(nspec_evolve), enuc

    ! This is basically e = m c**2

    enuc = sum(dydt(:) * aion(1:nspec_evolve) * ebin(1:nspec_evolve))

  end subroutine ener_gener_rate

end module actual_burner_module
