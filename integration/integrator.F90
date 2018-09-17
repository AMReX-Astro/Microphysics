module integrator_module

  use amrex_error_module

  implicit none

  public

contains

  subroutine integrator_init()

#if ((INTEGRATOR == 0 || INTEGRATOR == 1) && !defined(CUDA))
    use vode_integrator_module, only: vode_integrator_init
    use bs_integrator_module, only: bs_integrator_init
#else
    use actual_integrator_module, only: actual_integrator_init
#endif

    implicit none

#if ((INTEGRATOR == 0 || INTEGRATOR == 1) && !defined(CUDA))
    call vode_integrator_init()
    call bs_integrator_init()
#else
    call actual_integrator_init()
#endif

  end subroutine integrator_init



  subroutine integrator(state_in, state_out, dt, time)

    !$acc routine seq

#if ((INTEGRATOR == 0 || INTEGRATOR == 1) && !defined(CUDA))
    use vode_integrator_module, only: vode_integrator
    use bs_integrator_module, only: bs_integrator
#else
    use actual_integrator_module, only: actual_integrator
#endif
    use amrex_error_module, only: amrex_error
    use amrex_fort_module, only : rt => amrex_real
    use amrex_constants_module, only: ZERO, ONE
    use burn_type_module, only: burn_t
    use integration_data, only: integration_status_t
    use extern_probin_module, only: rtol_spec, rtol_temp, rtol_enuc, &
                                    atol_spec, atol_temp, atol_enuc, &
                                    abort_on_failure, &
                                    retry_burn, retry_burn_factor, retry_burn_max_change

    implicit none

    type (burn_t),  intent(in   ) :: state_in
    type (burn_t),  intent(inout) :: state_out
    real(rt),     intent(in   ) :: dt, time

    !$gpu

#if ((INTEGRATOR == 0 || INTEGRATOR == 1) && !defined(CUDA))
    type (integration_status_t) :: status
    real(rt) :: retry_change_factor
    integer :: current_integrator

    ! Loop through all available integrators. Our strategy will be to
    ! try the default integrator first. If the burn fails, we loosen
    ! the tolerance and try again, and keep doing so until we succed.
    ! If we keep failing and hit the threshold for how much tolerance
    ! loosening we will accept, then we restart with the original tolerance
    ! and switch to another integrator.

    do current_integrator = 0, 1

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

#if (INTEGRATOR == 0)
          if (current_integrator == 0) then
             call vode_integrator(state_in, state_out, dt, time, status)
          else if (current_integrator == 1) then
             call bs_integrator(state_in, state_out, dt, time, status)
          endif
#elif (INTEGRATOR == 1)
          if (current_integrator == 0) then
             call bs_integrator(state_in, state_out, dt, time, status)
          else if (current_integrator == 1) then
             call vode_integrator(state_in, state_out, dt, time, status)
          endif
#endif

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

             if (current_integrator < 1) then

#if (INTEGRATOR == 0)
                print *, "Retrying burn with BS integrator"
#elif (INTEGRATOR == 1)
                print *, "Retrying burn with VODE integrator"
#endif

             end if

             ! Switch to the next integrator (if there is one).

             exit

          end if

       end do

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
          call amrex_error("ERROR in burner: integration failed")
       end if

    endif

#else

    call actual_integrator(state_in, state_out, dt, time)

#endif

  end subroutine integrator

end module integrator_module
