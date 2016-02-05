module integration_module

  implicit none

  public

contains

  subroutine do_burn(state_in, state_out, dt, time)

    use extern_probin_module, only: integrator
    use vode_module, only: vode_burner
    use bl_error_module, only: bl_error
    use burn_type_module, only: burn_t

    implicit none

    type (burn_t),       intent(in   ) :: state_in
    type (burn_t),       intent(inout) :: state_out
    double precision,    intent(in   ) :: dt, time

    if (integrator == 0) then

       call vode_burner(state_in, state_out, dt, time)

    else

       call bl_error("Error: Unrecognized choice of integrator in do_burn.")

    endif

  end subroutine do_burn

end module integration_module
