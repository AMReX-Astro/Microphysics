! Common variables and routines for burners
! that use VBDF for their integration.

module actual_integrator_module

  use eos_module
  use network
  use rpar_indices
  use burn_type_module
  use bl_types
  use bdf_type_module
  use bdf

  implicit none

  real(kind=dp_t), parameter, private :: SMALL = 1.d-30


contains

  subroutine actual_integrator_init()

    use bdf, only: init_pascal

    implicit none

    call init_pascal()

  end subroutine actual_integrator_init



  ! Main interface

  subroutine actual_integrator(state_in, state_out, dt, time)

    !$acc routine seq

    use rpar_indices
    use extern_probin_module, only: jacobian, burner_verbose, &
                                    rtol_spec, rtol_temp, rtol_enuc, &
                                    atol_spec, atol_temp, atol_enuc, &
                                    burning_mode, retry_burn, &
                                    retry_burn_factor, retry_burn_max_change, &
                                    call_eos_in_rhs, dT_crit
    use integration_data, only: temp_scale, ener_scale
    use actual_rhs_module, only : update_unevolved_species

    implicit none

    ! Input arguments

    type (burn_t), intent(in   ) :: state_in
    type (burn_t), intent(inout) :: state_out
    real(dp_t),    intent(in   ) :: dt, time
    
    real(dp_t) :: dt_init
    logical, parameter :: RESET = .true.  !.true. means we want to initialize the bdf_ts object
    logical, parameter :: REUSE = .false. !.false. means don't reuse the Jacobian
    real(kind=dp_t), parameter :: DT0 = 1.0d-9 !Initial dt to be used in getting from 
                                               !t to tout.  Also arbitrary,
                                               !multiple values should be
                                               !explored.

    ! Local variables
    integer :: n, i, j, ierr

    real(kind=dp_t) :: atol(neqs), rtol(neqs)   ! input state, abs and rel tolerances
    real(kind=dp_t) :: y0(neqs,bdf_npt), y1(neqs,bdf_npt)
    real(kind=dp_t) :: t0, t1, enuc, dX

    type (eos_t)  :: eos_state_in, eos_state_out, eos_state_temp
    type (bdf_ts) :: ts

    real(dp_t) :: sum
    real(dp_t) :: retry_change_factor

    double precision :: ener_offset

    call bdf_ts_build(ts)

    ! Set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.
    !
    ! **NOTE** if you reduce these tolerances, you probably will need
    ! to (a) decrease dT_crit, (b) increase the maximum number of
    ! steps allowed.

    atol(1:nspec_evolve) = atol_spec ! mass fractions
    atol(net_itemp)      = atol_temp ! temperature
    atol(net_ienuc)      = atol_enuc ! energy generated

    rtol(1:nspec_evolve) = rtol_spec ! mass fractions
    rtol(net_itemp)      = rtol_temp ! temperature
    rtol(net_ienuc)      = rtol_enuc ! energy generated

    ts % atol = atol
    ts % rtol = rtol

    ! Initialize the integration time.

    t0 = ZERO
    t1 = t0 + dt

    ! Convert our input burn state into an EOS type.

    call burn_to_eos(state_in, eos_state_in)

    ! We assume that (rho, T) coming in are valid, do an EOS call
    ! to fill the rest of the thermodynamic variables.

    call eos(eos_input_rt, eos_state_in)

    ! Convert the EOS state data into the form VBDF expects.

    call eos_to_vbdf(eos_state_in, ts)

    ! Get the internal energy e that is consistent with this T.
    ! We will start the zone with this energy and subtract it at
    ! the end. This helps a lot with convergence, rather than
    ! starting e off at zero.

    ener_offset = eos_state_in % e

    ts % y(net_ienuc,1) = ener_offset / ener_scale

    ! Pass through whether we are doing self-heating.

    if (burning_mode == 0 .or. burning_mode == 2) then
       ts % upar(irp_self_heat,1) = -ONE
    else if (burning_mode == 1 .or. burning_mode == 3) then
       ts % upar(irp_self_heat,1) = ONE
    else
#ifndef ACC
       call bl_error("Error: unknown burning_mode in actual_integrator.f90.")
#endif
    endif

    ! Copy in the zone size.

    ts % upar(irp_dx,1) = state_in % dx

    ! Set the sound crossing time.

    ts % upar(irp_t_sound,1) = state_in % dx / eos_state_in % cs

    ! Set the time offset

    ts % upar(irp_t0,1) = time


    ! If we are using the dT_crit functionality and therefore doing a linear
    ! interpolation of the specific heat in between EOS calls, do a second
    ! EOS call here to establish an initial slope.

    ts % upar(irp_Told,1) = eos_state_in % T

    if (dT_crit < 1.0d19) then

       eos_state_temp = eos_state_in
       eos_state_temp % T = eos_state_in % T * (ONE + sqrt(epsilon(ONE)))

       call eos(eos_input_rt, eos_state_temp)

       ts % upar(irp_dcvdt,1) = (eos_state_temp % cv - eos_state_in % cv) / (eos_state_temp % T - eos_state_in % T)
       ts % upar(irp_dcpdt,1) = (eos_state_temp % cp - eos_state_in % cp) / (eos_state_temp % T - eos_state_in % T)

    endif

    do n = 1, neqs
       y0(n,1) = ts % y(n,1)
    end do

    ! Save the initial state.

    ts % upar(irp_y_init:irp_y_init + neqs - 1, 1) = y0(:,1)

    ! get the initial timestep estimate
    call initial_timestep(ts, t0, t1, dt_init)

    ! Call the integration routine.


    call bdf_advance(ts, y0, t0, y1, t1, dt_init, RESET, REUSE, ierr, .true.)

    do n = 1, neqs
       ts % y(n,1) = y1(n,1)
    end do

    ! If we are using hybrid burning and the energy release was negative (or we failed),
    ! re-run this in self-heating mode.

    if ( burning_mode == 2 .and. (ts % y(net_ienuc,1) * ener_scale - ener_offset < ZERO .or. ierr /= BDF_ERR_SUCCESS) ) then

       ts % upar(irp_self_heat,1) = ONE

       call eos_to_vbdf(eos_state_in, ts)

       ts % y(net_ienuc,1) = ener_offset / ener_scale

       ts % upar(irp_Told,1) = eos_state_in % T

       if (dT_crit < 1.0d19) then
          ts % upar(irp_dcvdt,1) = (eos_state_temp % cv - eos_state_in % cv) / (eos_state_temp % T - eos_state_in % T)
          ts % upar(irp_dcpdt,1) = (eos_state_temp % cp - eos_state_in % cp) / (eos_state_temp % T - eos_state_in % T)
       endif

       do n = 1, neqs
          y0(n,1) = ts % y(n,1)
       end do
       call bdf_advance(ts, y0, t0, y1, t1, DT0, RESET, REUSE, ierr, .true.)
       do n = 1, neqs
          ts % y(n,1) = y1(n,1)
       end do

    endif

    ! If we still failed, print out the current state of the integration.

    if (ierr /= BDF_ERR_SUCCESS) then
#ifndef ACC
       print *, 'ERROR: integration failed in net'
       print *, 'ierr = ', ierr
       print *, 'time = ', ts % t
       print *, 'dens = ', state_in % rho
       print *, 'temp start = ', state_in % T
       print *, 'xn start = ', state_in % xn
       print *, 'temp current = ', ts % y(net_itemp,1) * temp_scale
       print *, 'xn current = ', ts % y(1:nspec_evolve,1) * aion(1:nspec_evolve), &
            ts % upar(irp_nspec:irp_nspec+n_not_evolved-1,1) * aion(nspec_evolve+1:)
       print *, 'energy generated = ', ts % y(net_ienuc,1) * ener_scale - ener_offset
#endif
       if (.not. retry_burn) then
#ifndef ACC
          call bl_error("ERROR in burner: integration failed")
#endif
       else

          print *, 'Retrying burn with looser tolerances'

          retry_change_factor = ONE

          do while (ierr /= BDF_ERR_SUCCESS .and. retry_change_factor <= retry_burn_max_change)

             retry_change_factor = retry_change_factor * retry_burn_factor

             ts % atol = atol * retry_burn_factor
             ts % rtol = rtol * retry_burn_factor

             call eos_to_vbdf(eos_state_in, ts)

             ts % y(net_ienuc,1) = ener_offset / ener_scale

             ts % upar(irp_Told,1) = eos_state_in % T

             if (dT_crit < 1.0d19) then
                ts % upar(irp_dcvdt,1) = (eos_state_temp % cv - eos_state_in % cv) / (eos_state_temp % T - eos_state_in % T)
                ts % upar(irp_dcpdt,1) = (eos_state_temp % cp - eos_state_in % cp) / (eos_state_temp % T - eos_state_in % T)
             endif

             do n = 1, neqs
                y0(n,1) = ts % y(n,1)
             end do
             call bdf_advance(ts, y0, t0, y1, t1, DT0, RESET, REUSE, ierr, .true.)
             do n = 1, neqs
                ts % y(n,1) = y1(n,1)
             end do

          enddo

          if (retry_change_factor > retry_burn_max_change .and. ierr /= BDF_ERR_SUCCESS) then
#ifndef ACC
             call bl_error("ERROR in burner: integration failed")
#endif
          endif

       endif

    endif

    ! Subtract the energy offset.

    ts % y(net_ienuc,1) = ts % y(net_ienuc,1) * ener_scale - ener_offset

    ! Store the final data, and then normalize abundances.
    call vbdf_to_burn(ts, state_out)

    if (nspec_evolve < nspec) then
       call update_unevolved_species(state_out)
    endif

    call normalize_abundances_burn(state_out)

    ! store diagnostics
    state_out % n_rhs = ts % nfe
    state_out % n_jac = ts % nje

    if (burner_verbose) then

       ! Print out some integration statistics, if desired.

       print *, 'integration summary: '
       print *, 'dens: ', state_out % rho
       print *, ' temp: ', state_out % T
       print *, ' energy released: ', ts % y(net_ienuc,1) * ener_scale - ener_offset
       print *, 'number of steps taken: ', ts % n
       print *, 'number of f evaluations: ', ts % nfe

    endif

  end subroutine actual_integrator


  subroutine initial_timestep(ts, t0, t1, dt)

    ! this is a version of the timestep estimation algorithm used by
    ! VODE

    ! note, we only do this for the 1st vector in the solution

    !$acc routine seq

    use rhs_module

    type (bdf_ts), intent(inout) :: ts
    real (dp_t), intent(in) :: t0, t1
    real (dp_t), intent(out) :: dt
    type (bdf_ts) :: ts_temp
    real(kind=dp_t) :: h, h_old, hL, hU, ddydtt(neqs), eps, ewt(neqs), yddnorm
    integer :: n

    ts_temp = ts

    eps = maxval(ts % rtol)

    ! Initial lower and upper bounds on the timestep

    hL = 100.0d0 * epsilon(ONE) * max(abs(t0), abs(t1))
    hU = 0.1d0 * abs(t1 - t0)

    ! Initial guess for the iteration

    h = sqrt(hL * hU)
    h_old = 10.0 * h

    ! Iterate on ddydtt = (RHS(t + h, y + h * dydt) - dydt) / h
    call rhs(ts)

    do n = 1, 4

       h_old = h

       ! Get the error weighting -- this is similar to VODE's dewset
       ! routine

       ewt = eps * abs(ts % y(:,1)) + SMALL

       ! Construct the trial point.

       ts_temp % t = ts % t + h
       ts_temp % y(:,1) = ts % y(:,1) + h * ts % yd(:,1)


       ! Call the RHS, then estimate the finite difference.
       call rhs(ts_temp)
       ddydtt = (ts_temp % yd(:,1) - ts % yd(:,1)) / h

       yddnorm = sqrt( sum( (ddydtt*ewt)**2 ) / neqs )

       if (yddnorm*hU*hU > TWO) then
          h = sqrt(TWO / yddnorm)
       else
          h = sqrt(h * hU)
       endif

       if (h_old < TWO * h .and. h_old > HALF * h) exit

    enddo

    ! Save the final timestep, with a bias factor.
    dt = min(max(h/TWO, hL), hU)

  end subroutine initial_timestep

end module actual_integrator_module
