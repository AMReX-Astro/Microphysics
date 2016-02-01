module actual_burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use rpar_indices
  use eos_module
  use eos_data_module
  use eos_type_module
  use network
  use extern_probin_module, only: specific_q_burn

  implicit none

contains

  subroutine actual_burner_init()

    implicit none

    call init_rpar_indices(nrates, nspec)

  end subroutine actual_burner_init



  subroutine actual_burner(state_in, state_out, dt, time)

    use vode_module, only: vode_burner

    implicit none

    type (eos_t),     intent(in   ) :: state_in
    type (eos_t),     intent(inout) :: state_out
    double precision, intent(in   ) :: dt, time

    call vode_burner(state_in, state_out, dt, time)

  end subroutine actual_burner

end module actual_burner_module
