! Common variables and routines for burners
! that use VBDF for their integration.

module actual_integrator_module

  use microphysics_type_module, only: rt, HALF, TWO
  use amrex_error_module
  use eos_type_module
  use eos_module
  use network
  use vbdf_rpar_indices
  use burn_type_module
  use bdf_type_module
  use bdf

  implicit none

  real(rt), parameter, private :: SMALL = 1.e-30_rt


contains

  subroutine actual_integrator_init()

    use bdf, only: init_pascal

    implicit none

    call init_pascal()

  end subroutine actual_integrator_init



  ! Main interface

  subroutine actual_integrator(state_in, state_out, dt, time)

    !$acc routine seq

    use vbdf_rpar_indices
    use extern_probin_module, only: burner_verbose, &
                                    reuse_jac, &
                                    rtol_spec, rtol_temp, rtol_enuc, &
                                    atol_spec, atol_temp, atol_enuc, &
                                    burning_mode, burning_mode_factor, &
                                    retry_burn, retry_burn_factor, retry_burn_max_change, &
                                    call_eos_in_rhs, dT_crit
    use actual_rhs_module, only : update_unevolved_species
    use temperature_integration_module, only: self_heat

    implicit none

    ! Input arguments

    type (burn_t), intent(in   ) :: state_in
    type (burn_t), intent(inout) :: state_out
    real(rt),    intent(in   ) :: dt, time
    
    real(rt) :: dt_init
    logical, parameter :: RESET = .true.  !.true. means we want to initialize the bdf_ts object

    ! Local variables
    integer :: n, ierr

    real(rt) :: atol(neqs), rtol(neqs)   ! input state, abs and rel tolerances
    real(rt) :: y0(neqs,bdf_npt), y1(neqs,bdf_npt)
    real(rt) :: t0, t1

    type (eos_t)  :: eos_state_in, eos_state_temp
    type (bdf_ts) :: ts

    real(rt) :: retry_change_factor

    real(rt) :: ener_offset
    real(rt) :: edot, t_enuc, t_sound, limit_factor

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

    ts % y(net_ienuc,1) = ener_offset

    ! Pass through whether we are doing self-heating.

    if (self_heat) then
       ts % upar(irp_self_heat,1) = ONE
    else
       ts % upar(irp_self_heat,1) = -ONE
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

    if (dT_crit < 1.0e19_rt) then

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


    call bdf_advance(ts, y0, t0, y1, t1, dt_init, &
                     RESET, reuse_jac, ierr, .true.)

    do n = 1, neqs
       ts % y(n,1) = y1(n,1)
    end do

    ! If we are using hybrid burning and the energy release was
    ! negative (or we failed), re-run this in self-heating mode.

    if ( burning_mode == 2 .and. &
         (ts % y(net_ienuc,1) - ener_offset < ZERO .or. &
          ierr /= BDF_ERR_SUCCESS) ) then

       ts % upar(irp_self_heat,1) = ONE

       call eos_to_vbdf(eos_state_in, ts)

       ts % y(net_ienuc,1) = ener_offset

       ts % upar(irp_Told,1) = eos_state_in % T

       if (dT_crit < 1.0e19_rt) then
          ts % upar(irp_dcvdt,1) = (eos_state_temp % cv - eos_state_in % cv) / (eos_state_temp % T - eos_state_in % T)
          ts % upar(irp_dcpdt,1) = (eos_state_temp % cp - eos_state_in % cp) / (eos_state_temp % T - eos_state_in % T)
       endif

       do n = 1, neqs
          y0(n,1) = ts % y(n,1)
       end do
       call bdf_advance(ts, y0, t0, y1, t1, dt_init, RESET, &
                        reuse_jac, ierr, .true.)
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
       print *, 'temp current = ', ts % y(net_itemp,1)
       print *, 'xn current = ', ts % y(1:nspec_evolve,1), &
            ts % upar(irp_nspec:irp_nspec+n_not_evolved-1,1)
       print *, 'energy generated = ', ts % y(net_ienuc,1) - ener_offset
#endif
       if (.not. retry_burn) then
#ifndef ACC
          call amrex_error("ERROR in burner: integration failed")
#endif
       else

          print *, 'Retrying burn with looser tolerances'

          retry_change_factor = ONE

          do while (ierr /= BDF_ERR_SUCCESS .and. retry_change_factor <= retry_burn_max_change)

             retry_change_factor = retry_change_factor * retry_burn_factor

             ts % atol = atol * retry_burn_factor
             ts % rtol = rtol * retry_burn_factor

             call eos_to_vbdf(eos_state_in, ts)

             ts % y(net_ienuc,1) = ener_offset

             ts % upar(irp_Told,1) = eos_state_in % T

             if (dT_crit < 1.0e19_rt) then
                ts % upar(irp_dcvdt,1) = (eos_state_temp % cv - eos_state_in % cv) / (eos_state_temp % T - eos_state_in % T)
                ts % upar(irp_dcpdt,1) = (eos_state_temp % cp - eos_state_in % cp) / (eos_state_temp % T - eos_state_in % T)
             endif

             do n = 1, neqs
                y0(n,1) = ts % y(n,1)
             end do
             call bdf_advance(ts, y0, t0, y1, t1, dt_init, &
                              RESET, reuse_jac, ierr, .true.)
             do n = 1, neqs
                ts % y(n,1) = y1(n,1)
             end do

          enddo

          if (retry_change_factor > retry_burn_max_change .and. ierr /= BDF_ERR_SUCCESS) then
#ifndef ACC
             call amrex_error("ERROR in burner: integration failed")
#endif
          endif

       endif

    endif

    ! Subtract the energy offset.

    ts % y(net_ienuc,1) = ts % y(net_ienuc,1) - ener_offset

    ! Store the final data, and then normalize abundances.
    call vbdf_to_burn(ts, state_out)

    if (nspec_evolve < nspec) then
       call update_unevolved_species(state_out)
    endif

    ! For burning_mode == 3, limit the burning.

    if (burning_mode == 3) then

       t_enuc = eos_state_in % e / max(abs(state_out % e - state_in % e) / max(dt, 1.e-50_rt), 1.e-50_rt)
       t_sound = state_in % dx / eos_state_in % cs

       limit_factor = min(1.0e0_rt, burning_mode_factor * t_enuc / t_sound)

       state_out % e = state_in % e + limit_factor * (state_out % e - state_in % e)
       state_out % xn(:) = state_in % xn(:) + limit_factor * (state_out % xn(:) - state_in % xn(:))

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
       print *, ' energy released: ', ts % y(net_ienuc,1) - ener_offset
       print *, 'number of steps taken: ', ts % n
       print *, 'number of f evaluations: ', ts % nfe

    endif

  end subroutine actual_integrator


  subroutine initial_timestep(ts, t0, t1, dt)

    ! this is a version of the timestep estimation algorithm used by
    ! VODE

    !$acc routine seq

    use rhs_module

    type (bdf_ts), intent(inout) :: ts
    real (rt), intent(in) :: t0, t1
    real (rt), intent(out) :: dt
    type (bdf_ts) :: ts_temp
    real (rt) :: h, h_old, hL, hU, ddydtt(neqs), eps, ewt(neqs), yddnorm
    integer :: n

    ts_temp = ts

    eps = maxval(ts % rtol)

    ! Initial lower and upper bounds on the timestep

    hL = 100.0e0_rt * epsilon(ONE) * max(abs(t0), abs(t1))
    hU = 0.1e0_rt * abs(t1 - t0)

    ! Initial guess for the iteration

    h = sqrt(hL * hU)
    h_old = 10.0_rt * h

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
