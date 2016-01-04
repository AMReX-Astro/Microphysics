module actual_burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module
  use network

  implicit none
  
  ! Set the number of independent variables -- this should be
  ! temperature, enuc + the number of species which participate
  ! in the evolution equations.
  ! 
  integer, parameter :: NEQ = 2 + nspec

  ! Our problem is stiff, so tell ODEPACK that. 21 means stiff, jacobian 
  ! function is supplied; 22 means stiff, figure out my jacobian through 
  ! differencing.

  integer, parameter :: MF_ANALYTIC_JAC = 21, MF_NUMERICAL_JAC = 22
  integer, parameter :: MF_JAC = MF_NUMERICAL_JAC
  
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
    use network_indices
    use mempool_module, only: bl_allocate, bl_deallocate
    
    implicit none

    type (eos_t),        intent(in   ) :: state_in
    type (eos_t),        intent(inout) :: state_out
    double precision,    intent(in   ) :: dt, time
    
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

    double precision :: sum

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
    
    atol(1:nspec)   = 1.d-12 ! mass fractions
    atol(net_itemp) = 1.d-8  ! temperature
    atol(net_ienuc) = 1.d-8  ! energy generated
    
    rtol(1:nspec)   = 1.d-12 ! mass fractions
    rtol(net_itemp) = 1.d-5  ! temperature
    rtol(net_ienuc) = 1.d-10 ! energy generated

    ! We want VODE to re-initialize each time we call it.
    
    istate = 1

    ! Initialize work arrays to zero.

    rwork(:) = ZERO
    iwork(:) = 0

    ! Set the maximum number of steps allowed (the VODE default is 500).

    iwork(6) = 150000

    ! Initialize the integration time.

    local_time = ZERO

    ! Abundances are the first nspec values and temperature and energy are the last.

    y(1:nspec)   = state_in % xn(:)
    y(net_itemp) = state_in % T
    y(net_ienuc) = ZERO

    ! Density, specific heat at constant pressure (c_p), and dhdX are needed
    ! in the righthand side routine, so we will pass these in to those routines
    ! via the rpar functionality in VODE.

    rpar(irp_dens) = state_in % rho
    rpar(irp_cp)   = state_in % cp
    rpar(irp_cv)   = state_in % cv

    ! The following thermodynamic derivatives are calculated in the EOS.

    ! dhdX = dh/dX |_{p,T}

    rpar(irp_dhdX:irp_dhdX-1+nspec) = state_in % dhdX(:)

    ! dedX = de/dX |_{rho, T}

    rpar(irp_dedX:irp_dedX-1+nspec) = state_in % dEdX(:)

    ! This is just used to make sure everything is happy.

    rpar(irp_smallx) = smallx



    ! Call the integration routine.

    call dvode(f_rhs, NEQ, y, local_time, dt, ITOL, rtol, atol, ITASK, &
         istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_JAC,&
         rpar, ipar)



    if (istate < 0) then
       print *, 'ERROR: integration failed in net'
       print *, 'istate = ', istate
       print *, 'time = ', local_time
       print *, 'dens = ', state_in % rho
       print *, 'temp start = ', state_in % T
       print *, 'xn start = ', state_in % xn
       print *, 'temp current = ', y(net_itemp)
       print *, 'xn current = ', y(1:nspec)
       call bl_error("ERROR in burner: integration failed")
    endif

    ! Store the new mass fractions.
    state_out % xn(:) = max(smallx, min(ONE, y(1:nspec)))    


    ! Energy was integrated in the system -- we use this integrated
    ! energy which contains both the reaction energy release and
    ! neutrino losses. The final energy is the initial energy
    ! plus this energy release. Note that we get a new temperature too,
    ! but we will discard it and the main burner module will do an EOS
    ! call to get a final temperature consistent with this new energy.

    state_out % e = state_in % e + y(net_ienuc)

    if (verbose) then

       ! Print out some integration statistics, if desired.

       print *, 'integration summary: '
       print *, 'dens: ', state_out % rho, ' temp: ', state_out % T
       print *, 'number of steps taken: ', iwork(11)
       print *, 'number of f evaluations: ', iwork(12)

    endif

    call bl_deallocate(y)
    call bl_deallocate(atol)
    call bl_deallocate(rtol)
    call bl_deallocate(rwork)
    call bl_deallocate(iwork)
    call bl_deallocate(rpar)
    
  end subroutine actual_burner

end module actual_burner_module
