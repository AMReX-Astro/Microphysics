! This module contains a version of Wallace and Woosley's (ApJS 45,389
! (1981)) rprox reaction network burner.
!
!

module actual_burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module
  use eos_type_module
  use network

  double precision, parameter :: T2T9 = 1.0d-9

  ! Set the number of independent variables -- this should be
  ! temperature + the number of species which participate
  ! in the evolution equations.
  ! 
  integer, parameter :: NEQ = 1 + nspec

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
  
  integer, parameter :: LRW = 22 + 9*NEQ + 2*NEQ**2
  integer, parameter :: LIW = 30 + NEQ

  ! Should we print out diagnostic output after the solve?
  
  logical, parameter :: verbose = .false.    
  
contains

  subroutine actual_burner(state_in, state_out, dt, time)  

    use rpar_indices

    implicit none

    type (eos_t),     intent(in   ) :: state_in
    type (eos_t),     intent(inout) :: state_out
    double precision, intent(in   ) :: dt, time    

    integer :: j

    double precision :: local_time
    
    ! Work arrays
    
    double precision, pointer :: y(:)
    double precision, pointer :: atol(:), rtol(:)
    double precision, pointer :: rwork(:)
    integer,          pointer :: iwork(:)    

    ! istate determines the state of the calculation.  A value of 1 meeans
    ! this is the first call to the problem -- this is what we will want.
    
    integer :: istate
    
    double precision, pointer :: rpar(:)
    integer :: ipar

    double precision :: T9

    EXTERNAL jac, f_rhs    

    ! Allocate storage for work arrays and for rpar, the
    ! array through which we communicate with VODE from the
    ! RHS routines.
    
    call bl_allocate(y, 1, NEQ)
    call bl_allocate(atol, 1, NEQ)
    call bl_allocate(rtol, 1, NEQ)
    call bl_allocate(rwork, 1, LRW)
    call bl_allocate(iwork, 1, LIW)
    call bl_allocate(rpar, 1, n_rpar_comps)


    ! Set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.  
    !
    ! **NOTE** if you reduce these tolerances, you probably will need
    ! to (a) decrease dT_crit, (b) increase the maximum number of 
    ! steps allowed.

    atol(1:nspec) = 1.d-11  ! mass fractions
    atol(nspec+1) = 1.d-8   ! temperature

    rtol(1:nspec) = 1.d-12  ! mass fractions
    rtol(nspec+1) = 1.d-8   ! temperature

    ! We want VODE to re-initialize each time we call it.
    
    istate = 1

    ! Initialize work arrays to zero.       

    rwork(:) = ZERO
    iwork(:) = 0

    ! Set the maximum number of steps allowed (the VODE default is 500).

    iwork(6) = 15000

    ! Initialize the integration time.

    local_time = ZERO

    ! Temperature (in units of 10**9).

    T9 = state_in % T * T2T9

    ! Abundances are the first nspec values and temperature (T9) is the last.

    y(1:nspec) = state_in % xn(:) / aion(:)
    y(nspec+1) = T9

    ! Set the thermodynamics that are passed via rpar to the RHS routine.

    rpar(irp_dens) = state_in % rho
    rpar(irp_T9_eos) = T9
    rpar(irp_dTcrit) = 1.0d15
    rpar(irp_cp) = state_in % cp
    rpar(irp_dhdX:irp_dhdX+nspec-1) = state_in % dhdX(:)

    ! Only burn if 0.2 < T9 < 2.5 or X(H1) > 0.05.
    ! The last restriction is a kludge based on the last paragraph of WW81.

    if ((T9 .gt. 0.2d0 .and. T9 .lt. 2.5d0) .or. state_in % xn(ih1) > 0.05d0) then

       ! Call the integration routine.

       call dvode(f_rhs, NEQ, y, time, dt, ITOL, rtol, atol, ITASK, &
            istate, IOPT, rwork, LRW, iwork, LIW, jac, &
                                !MF_NUMERICAL_JAC, &
            MF_ANALYTIC_JAC,&
            rpar, ipar)

       if (istate < 0) then
          print *, 'ERROR: integration failed in net'
          print *, 'initial T   = ', state_in % T
          print *, 'initial rho = ', state_in % rho
          print *, 'iniitial X  = ', state_in % xn(:)
          print *, 'istate = ', istate
          print *, 'time = ', local_time
          print *, 'dt = ', dt
          print *, 'output state = ', y
          call bl_error("ERROR in burner: integration failed")
       endif

       ! Store the new mass fractions.
       state_out % xn(:) = max(smallx, min(ONE, y(1:nspec) * aion(:)))

       ! Update the energy.
       state_out % e = state_in % e + sum((state_in % xn(:) - state_out % xn(:)) * ebin(:))

       if (verbose) then
          ! print out some integration statistics, if desired
          print *, 'integration summary: '
          print *, 'dens: ', state_in % rho, ' temp: ', state_in % T
          print *, 'number of steps taken: ', iwork(11)
          print *, 'number of f evaluations: ', iwork(12)
       endif

    endif
    
  end subroutine actual_burner


end module actual_burner_module
