module actual_burner_module

  use burn_type_module
  use actual_network
  use actual_network_data

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



  subroutine ener_gener_rate(dydt, ebin, enuc)

    implicit none

    double precision :: dydt(nspec), ebin(nspec), enuc

    enuc = sum(dydt(:) * aion(:) * ebin(:))

  end subroutine ener_gener_rate

end module actual_burner_module
