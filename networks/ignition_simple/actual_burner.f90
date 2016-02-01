module actual_burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module
  use eos_data_module
  use eos_type_module
  use burn_type_module
  use network
  use actual_burner_data

  implicit none

contains

  subroutine actual_burner_init()

    use rpar_indices

    implicit none

    call init_rpar_indices(nrates, nspec)

    irp_rate      = get_next_rpar_index(1)
    irp_dratedt   = get_next_rpar_index(1)
    irp_sc1212    = get_next_rpar_index(1)
    irp_dsc1212dt = get_next_rpar_index(1)
    irp_xc12tmp   = get_next_rpar_index(1)    

  end subroutine actual_burner_init



  subroutine actual_burner(state_in, state_out, dt, time)

    use vode_module, only: vode_burner

    implicit none

    type (burn_t),    intent(in   ) :: state_in
    type (burn_t),    intent(inout) :: state_out
    double precision, intent(in   ) :: dt, time

    call vode_burner(state_in, state_out, dt, time)

  end subroutine actual_burner

end module actual_burner_module
