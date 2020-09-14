module integrator_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  subroutine integrator_init()

    use vode_integrator_module, only: vode_integrator_init

#ifdef NONAKA_PLOT
    use nonaka_plot_module, only: nonaka_init
#endif

    implicit none

    call vode_integrator_init()

#ifdef NONAKA_PLOT
    call nonaka_init()
#endif

  end subroutine integrator_init



  subroutine integrator(state_in, state_out, dt, time)

    !$acc routine seq

    use vode_integrator_module, only: vode_integrator
#ifndef AMREX_USE_CUDA
    use amrex_error_module, only: amrex_error
#endif
    use amrex_constants_module, only: ZERO, ONE
    use integration_data, only: integration_status_t
    use sdc_type_module, only: sdc_t
    use extern_probin_module, only: rtol_spec, rtol_temp, rtol_enuc, &
                                    atol_spec, atol_temp, atol_enuc, &
                                    abort_on_failure, &
                                    retry_burn, retry_burn_factor, retry_burn_max_change

    implicit none

    type (sdc_t),  intent(in   ) :: state_in
    type (sdc_t),  intent(inout) :: state_out
    real(rt),    intent(in   ) :: dt, time

    type (integration_status_t) :: status
    real(rt) :: retry_change_factor
    integer :: current_integrator

    !$gpu

    retry_change_factor = ONE

    status % atol_spec = atol_spec
    status % rtol_spec = rtol_spec

    status % atol_temp = atol_temp
    status % rtol_temp = rtol_temp

    status % atol_enuc = atol_enuc
    status % rtol_enuc = rtol_enuc

    ! Loop until we either succeed at the integration, or run
    ! out of tolerance loosening to try.

    do

       call vode_integrator(state_in, state_out, dt, time, status)

       if (state_out % success) exit

       if (.not. retry_burn) exit

       ! If we got here, the integration failed; loosen the tolerances.

       if (retry_change_factor < retry_burn_max_change) then

          print *, "Retrying burn with looser tolerances in zone ", state_in % i, state_in % j, state_in % k

          retry_change_factor = retry_change_factor * retry_burn_factor

          status % atol_spec = status % atol_spec * retry_burn_factor
          status % rtol_spec = status % rtol_spec * retry_burn_factor

          status % atol_temp = status % atol_temp * retry_burn_factor
          status % rtol_temp = status % rtol_temp * retry_burn_factor

          status % atol_enuc = status % atol_enuc * retry_burn_factor
          status % rtol_enuc = status % rtol_enuc * retry_burn_factor

          print *, "New tolerance loosening factor = ", retry_change_factor

       end if

       ! No need to do the next integrator if we have already succeeded.

       if (state_out % success) exit

       if (.not. retry_burn) exit

    end do

    ! If we get to this point and have not succeded, all available integrators have
    ! failed at all available tolerances; we must either abort, or return to the
    ! driver routine calling the burner, and attempt some other approach such as
    ! subcycling the main advance.

    if (.not. state_out % success) then

       if (abort_on_failure) then
#ifndef CUDA
          call amrex_error("ERROR in burner: integration failed")
#else
          stop
#endif
       end if

    endif

  end subroutine integrator

end module integrator_module
