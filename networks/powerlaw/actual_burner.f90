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

    use rpar_indices

    implicit none

    call init_rpar_indices(nrates, nspec)

  end subroutine actual_burner_init



  subroutine actual_burner(state_in, state_out, dt, time)

    use vode_module, only: vode_burner

    implicit none

    type (burn_t),    intent(in   ) :: state_in
    type (burn_t),    intent(inout) :: state_out
    double precision, intent(in   ) :: dt, time

    call vode_burner(state_in, state_out, dt, time)

  end subroutine actual_burner



  ! Computes the instantaneous energy generation rate

  subroutine ener_gener_rate(dydt, enuc)

    implicit none

    double precision :: dydt(nspec), enuc

    ! This is basically e = m c**2

    enuc = sum(dydt(:) * aion(:) * ebin(:))

  end subroutine ener_gener_rate

end module actual_burner_module
