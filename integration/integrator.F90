module integrator_module

  implicit none

  public

contains

  subroutine integrator_init()

    use vode_integrator_module, only: vode_integrator_init
    !use vode90_integrator_module, only: vode90_integrator_init
    use bs_integrator_module, only: bs_integrator_init
    use vbdf_integrator_module, only: vbdf_integrator_init

    implicit none

    call vode_integrator_init()
    !call vode90_integrator_init()
    call bs_integrator_init()
    call vbdf_integrator_init()

  end subroutine integrator_init



  subroutine integrator(state_in, state_out, dt, time)

    !$acc routine seq

    use vode_integrator_module, only: vode_integrator
    !use vode90_integrator_module, only: vode90_integrator
    use bs_integrator_module, only: bs_integrator
    use vbdf_integrator_module, only: vbdf_integrator
    use bl_error_module, only: bl_error
    use bl_constants_module, only: ZERO, ONE
    use burn_type_module, only: burn_t
    use bl_types, only: dp_t
    use integration_data, only: integration_status_t
    use extern_probin_module, only: rtol_spec, rtol_temp, rtol_enuc, &
                                    atol_spec, atol_temp, atol_enuc, &
                                    retry_burn, retry_burn_factor, retry_burn_max_change

    implicit none

    type (burn_t),  intent(in   ) :: state_in
    type (burn_t),  intent(inout) :: state_out
    real(dp_t),     intent(in   ) :: dt, time

    type (integration_status_t) :: status
    real(dp_t) :: retry_change_factor

#if (INTEGRATOR == VODE || INTEGRATOR == BS)

    retry_change_factor = ONE

    status % atol_spec = atol_spec
    status % rtol_spec = rtol_spec

    status % atol_temp = atol_temp
    status % rtol_temp = rtol_temp

    status % atol_enuc = atol_enuc
    status % rtol_enuc = rtol_enuc

    do while (retry_change_factor <= retry_burn_max_change)

       status % integration_complete = .true.

#if INTEGRATOR == VODE
       call vode_integrator(state_in, state_out, dt, time, status)
#elif INTEGRATOR == BS
       call bs_integrator(state_in, state_out, dt, time, status)
#else
       call bl_error("Unknown integrator.")
#endif

       if (status % integration_complete) exit

       ! If we got here, the integration failed, try the next burner.

       status % integration_complete = .true.

#if INTEGRATOR == VODE
       print *, "Retrying burn with BS integrator"
       call bs_integrator(state_in, state_out, dt, time, status)
#elif INTEGRATOR == BS
       print *, "Retrying burn with VODE integrator"
       call vode_integrator(state_in, state_out, dt, time, status)
#else
       call bl_error("Unknown integrator.")
#endif

       if (status % integration_complete) exit

       ! If we got here, all integrations failed, loosen the tolerances.

       if (.not. retry_burn) then

          call bl_error("ERROR in burner: integration failed")

       else

          print *, "Retrying burn with looser tolerances"

          retry_change_factor = retry_change_factor * retry_burn_factor

          status % atol_spec = status % atol_spec * retry_burn_factor
          status % rtol_spec = status % rtol_spec * retry_burn_factor

          status % atol_temp = status % atol_temp * retry_burn_factor
          status % rtol_temp = status % rtol_temp * retry_burn_factor

          status % atol_enuc = status % atol_enuc * retry_burn_factor
          status % rtol_enuc = status % rtol_enuc * retry_burn_factor

       end if

    end do

    if (.not. status % integration_complete) then

       call bl_error("ERROR in burner: integration failed")

    endif

#elif (INTEGRATOR == VBDF)

    call vbdf_integrator(state_in, state_out, dt, time)

! #elif (INTEGRATOR == VODE90)

!     call vode90_integrator(state_in, state_out, dt, time)

#else

    call bl_error("Unknown burner.")

#endif

  end subroutine integrator

end module integrator_module
