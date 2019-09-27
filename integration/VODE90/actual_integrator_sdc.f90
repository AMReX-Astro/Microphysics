! Common variables and routines for burners
! that use VODE for their integration.

module actual_integrator_module

  use amrex_error_module
  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module

  use sdc_type_module
  use vode_type_module
  use vode_rhs_module

  implicit none

  ! Our problem is stiff, so tell ODEPACK that. 21 means stiff, jacobian
  ! function is supplied; 22 means stiff, figure out my jacobian through
  ! differencing.

  integer, parameter :: MF_ANALYTIC_JAC = 21, MF_NUMERICAL_JAC = 22

  ! Tolerance parameters:
  !
  !  itol specifies whether to use an single absolute tolerance for
  !  all variables (1), or to pass an array of absolute tolerances, one
  !  for each variable with a scalar relative tol (2), a scalar absolute
  !  and array of relative tolerances (3), or arrays for both (4).
  !
  !  The error is determined as e(i) = rtol*abs(y(i)) + atol, and must
  !  be > 0.  Since we have some compositions that may be 0 initially,
  !  we will specify both an absolute and a relative tolerance.
  !
  ! We will use arrays for both the absolute and relative tolerances,
  ! since we want to be easier on the temperature than the species.

  integer, parameter :: ITOL = 4

  ! We want to do a normal computation, and get the output values of y(t)
  ! after stepping though dt.

  integer, PARAMETER :: ITASK = 1

  ! We will override the maximum number of steps, so turn on the
  ! optional arguments flag.

  integer, parameter :: IOPT = 1

  ! Declare a real work array of size 22 + 9*NEQ + 2*NEQ**2 and an
  ! integer work array of size 30 + NEQ. These are VODE constants
  ! that depend on the integration mode we're using -- see dvode.f.

  integer, parameter :: LRW = 22 + 9*VODE_NEQS + 2*VODE_NEQS**2
  integer, parameter :: LIW = 30 + VODE_NEQS

contains

  subroutine actual_integrator_init()

  end subroutine actual_integrator_init


  subroutine actual_integrator(state_in, state_out, dt, time)

    use vode_rpar_indices
    use extern_probin_module, only: jacobian, burner_verbose, &
                                    rtol_spec, rtol_temp, rtol_enuc, &
                                    atol_spec, atol_temp, atol_enuc, &
                                    burning_mode, retry_burn, &
                                    retry_burn_factor, retry_burn_max_change, &
                                    call_eos_in_rhs, dT_crit

    ! Input arguments

    type (sdc_t), intent(in   ) :: state_in
    type (sdc_t), intent(inout) :: state_out
    real(rt),    intent(in   ) :: dt, time

    ! Local variables

    real(rt) :: local_time

    ! Work arrays

    real(rt) :: y(VODE_NEQS)
    real(rt) :: atol(VODE_NEQS), rtol(VODE_NEQS)
    real(rt) :: rwork(LRW)
    integer    :: iwork(LIW)
    real(rt) :: rpar(n_rpar_comps)

    integer :: MF_JAC

    ! istate determines the state of the calculation.  A value of 1 meeans
    ! this is the first call to the problem -- this is what we will want.

    integer :: istate

    integer :: ipar

    real(rt) :: sum
    real(rt) :: retry_change_factor

    if (jacobian == 1) then ! Analytical
       MF_JAC = MF_ANALYTIC_JAC
    else if (jacobian == 2) then ! Numerical
       MF_JAC = MF_NUMERICAL_JAC
    else
       call amrex_error("Error: unknown Jacobian mode in actual_integrator.f90.")
    endif

    ! Set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.
    !
    ! **NOTE** if you reduce these tolerances, you probably will need
    ! to (a) decrease dT_crit, (b) increase the maximum number of
    ! steps allowed.

    atol(SFS:SFS-1+nspec) = atol_spec ! mass fractions
    atol(SEDEN)           = atol_enuc ! total energy
    atol(SEINT)           = atol_enuc ! internal energy

    rtol(SFS:SFS-1+nspec) = rtol_spec ! mass fractions
    rtol(SEDEN)           = rtol_enuc ! total energy
    rtol(SEINT)           = rtol_enuc ! internal energy

    ! We want VODE to re-initialize each time we call it.

    istate = 1

    ! Initialize work arrays to zero.

    rwork(:) = ZERO
    iwork(:) = 0

    ! Set the maximum number of steps allowed (the VODE default is 500).

    iwork(6) = 150000

    ! Disable printing of messages about T + H == T unless we are in verbose mode.

    if (burner_verbose) then
       iwork(7) = 1
    else
       iwork(7) = 0
    endif

    ! Initialize the integration time.

    local_time = ZERO

    ! Convert our input sdc state into the form VODE expects

    call sdc_to_vode(state_in, y, rpar)


    ! this is not used but we set it to prevent accessing uninitialzed
    ! data in common routines with the non-SDC integrator
    rpar(irp_self_heat) = -ONE

    ! Set the time offset -- this converts between the local integration
    ! time and the simulation time
    rpar(irp_t0) = time


    ! Call the integration routine.
    call dvode(f_rhs, VODE_NEQS, y, local_time, local_time + dt, &
               ITOL, rtol, atol, ITASK, &
               istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_JAC, rpar, ipar)


    ! Store the final data
    call vode_to_sdc(time, y, rpar, state_out)

    ! get the number of RHS calls and jac evaluations from the VODE
    ! work arrays
    state_out % n_rhs = iwork(12)
    state_out % n_jac = iwork(13)

  end subroutine actual_integrator

end module actual_integrator_module
