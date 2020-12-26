module integrator_module

  use amrex_error_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  subroutine integrator_init()

#ifdef INTEGRATOR_HAS_FORTRAN_IMPLEMENTATION
    use integrator_scaling_module, only: integrator_scaling_init
    use vode_integrator_module, only: vode_integrator_init
#else
    use actual_integrator_module, only: actual_integrator_init
#endif
    use temperature_integration_module, only: temperature_rhs_init

#ifdef NONAKA_PLOT
    use nonaka_plot_module, only: nonaka_init
#endif

    implicit none

    call integrator_scaling_init()
#if (INTEGRATOR == 0)
    call vode_integrator_init()
#else
    call actual_integrator_init()
#endif
    call temperature_rhs_init()

#ifdef NONAKA_PLOT
    call nonaka_init()
#endif

  end subroutine integrator_init



  subroutine integrator(state_in, state_out, dt, time)

    !$acc routine seq

#ifdef INTEGRATOR_HAS_FORTRAN_IMPLEMENTATION
#if (INTEGRATOR == 0)
    use vode_integrator_module, only: vode_integrator
#else
    use actual_integrator_module, only: actual_integrator
#endif
    use amrex_constants_module, only: ZERO, ONE
    use integration_data, only: integration_status_t
    use extern_probin_module, only: rtol_spec, rtol_temp, rtol_enuc, &
                                    atol_spec, atol_temp, atol_enuc, &
                                    abort_on_failure, &
                                    retry_burn, retry_burn_factor, retry_burn_max_change
#endif
    use burn_type_module, only: burn_t

    implicit none

    type (burn_t),  intent(in   ) :: state_in
    type (burn_t),  intent(inout) :: state_out
    real(rt),     intent(in   ) :: dt, time

    !$gpu

#ifdef INTEGRATOR_HAS_FORTRAN_IMPLEMENTATION
    type (integration_status_t) :: status
    real(rt) :: retry_change_factor
    integer :: current_integrator

    ! Integrate.  If the burn fails, we loosen the tolerance and try
    ! again, and keep doing so until we succed.  If we keep failing
    ! and hit the threshold for how much tolerance loosening we will
    ! accept, then we abort.

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

       else
          ! we've looosened as much as we are willing to, so bail
          exit
       end if

    end do

    ! If we get to this point and have not succeded, we must either
    ! abort, or return to the driver routine calling the burner, and
    ! attempt some other approach such as subcycling the main advance.

    if (.not. state_out % success) then

       if (abort_on_failure) then
#ifndef CUDA
          call amrex_error("ERROR in burner: integration failed")
#else
          stop
#endif
       end if

    endif

#else

    call actual_integrator(state_in, state_out, dt, time)

#endif

  end subroutine integrator

end module integrator_module
