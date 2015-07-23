module specific_burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module
  use eos_data_module
  use eos_type_module
  use network

  implicit none

  logical, parameter :: verbose = .false.

  ! set the number of independent variables -- this should be temperature
  ! + the number of species
  integer, parameter :: NEQ = 1 + nspec_advance

  ! our problem is stiff, tell ODEPACK that. 21 means stiff, jacobian 
  ! function is supplied, 22 means stiff, figure out my jacobian through 
  ! differencing
  integer, parameter :: MF_ANALYTIC_JAC = 21, MF_NUMERICAL_JAC = 22

  ! we want to do a normal computation, and get the output values of y(t)
  ! after stepping though dt
  integer, PARAMETER :: ITASK = 1

  ! tolerance parameters:
  !
  !  itol specifies whether to use an single absolute tolerance for
  !  all variables (1), or to pass an array of absolute tolerances, one
  !  for each variable with a scalar relative tol (2), a scalar absolute
  !  and array of relative tolerances (3), or arrays for both (4)
  !  
  !  The error is determined as e(i) = rtol*abs(y(i)) + atol, and must
  !  be > 0.  Since we have some compositions that may be 0 initially,
  !  we will specify both an absolute and a relative tolerance.
  !
  ! We will use arrays for both the absolute and relative tolerances, 
  ! since we want to be easier on the temperature than the species
  integer, parameter :: ITOL = 4

  ! declare a real work array of size 22 + 9*NEQ + 2*NEQ**2 and an
  ! integer work array of since 30 + NEQ
  integer, parameter :: LRW = 22 + 9*NEQ + 2*NEQ**2

  integer, parameter :: LIW = 30 + NEQ

  ! we will override the maximum number of steps, so turn on the 
  ! optional arguments flag
  integer, parameter :: IOPT = 1

contains

  subroutine specific_burner(state_in, state_out, dt, time)

    use burner_aux_module, only : dens_pass, c_p_pass, dhdx_pass, X_O16_pass

    implicit none

    type (eos_t_vector) :: state_in, state_out
    double precision    :: dt, time

    double precision :: enuc

    ! allocate storage for the input state
    double precision, dimension(NEQ) :: y

    double precision, dimension(NEQ) :: atol, rtol

    double precision :: integration_time

    ! istate determines the state of the calculation.  A value of 1 meeans
    ! this is the first call to the problem -- this is what we will want.
    ! Note, istate is changed over the course of the calculation, so it
    ! cannot be a parameter
    integer :: istate
    
    double precision, dimension(LRW) :: rwork
    
    integer, dimension(LIW) :: iwork
  
    double precision :: rpar
    integer :: ipar

    integer :: j

    EXTERNAL jac, f_rhs
    
    ! set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.  
    atol(1:nspec_advance) = 1.d-12    ! mass fractions
    atol(nspec_advance+1) = 1.d-8     ! temperature
       
    rtol(1:nspec_advance) = 1.d-12    ! mass fractions
    rtol(nspec_advance+1) = 1.d-5     ! temperature

    do j = 1, state_in % N

       ! we want VODE to re-initialize each time we call it
       istate = 1

       rwork(:) = ZERO
       iwork(:) = 0

       ! set the maximum number of steps allowed (the VODE default is 500)
       iwork(6) = 15000

       ! initialize the integration time
       integration_time = ZERO

       ! abundances are the first nspec_advance values and temperature is the last
       y(ic12) = state_in % xn(j,ic12)
       y(nspec_advance+1) = state_in % T(j)

       ! Density, specific heat at constant pressure, c_p, and dhdX are needed
       ! in the righthand side routine, so we will pass these in through the
       ! burner_aux module.

       dens_pass    = state_in % rho(j)
       c_p_pass     = state_in % cp(j)
       dhdx_pass(:) = state_in % dhdX(j,:)
       X_O16_pass   = state_in % xn(j,io16)

       ! call the integration routine
       call dvode(f_rhs, NEQ, y, integration_time, dt, ITOL, rtol, atol, ITASK, &
                  istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_ANALYTIC_JAC, &
                  rpar, ipar)

       if (istate < 0) then
          print *, 'ERROR: integration failed in net'
          print *, 'istate = ', istate
          print *, 'time = ', integration_time
          call bl_error("ERROR in burner: integration failed")
       endif


       ! store the new mass fractions -- note, we discard the temperature
       ! here and instead compute the energy release from the binding
       ! energy -- make sure that they are positive
       state_out % xn(j,ic12)  = max(y(ic12), ZERO)
       state_out % xn(j,io16)  = state_in % xn(j,io16)
       state_out % xn(j,img24) = ONE - state_out % xn(j,ic12) - state_out % xn(j,io16)

       ! compute the energy release and update the enthalpy.  Our convention
       ! is that the binding energies are negative, so the energy release is
       ! - sum_k { (Xout(k) - Xin(k)) ebin(k) }
       !
       ! since this version of the network only evolves C12, we can
       ! compute the energy release easily
       enuc = (ebin(img24) - ebin(ic12))*(state_out % xn(j,ic12) - state_in % xn(j,ic12))

       state_out % e(j) = state_in % e(j) + enuc

       if (verbose) then

          ! print out some integration statistics, if desired
          print *, 'integration summary: '
          print *, 'dens: ', state_out % rho(j), ' temp: ', state_out % T(j)
          print *, 'number of steps taken: ', iwork(11)
          print *, 'number of f evaluations: ', iwork(12)

       endif

    enddo

  end subroutine specific_burner

end module specific_burner_module
