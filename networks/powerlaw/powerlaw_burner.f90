module actual_burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use rpar_indices
  use eos_module
  use eos_data_module
  use eos_type_module
  use network
  use extern_probin_module, only: specific_q_burn

  implicit none

  logical, parameter :: verbose = .false.

  ! set the number of independent variables -- this should be temperature
  ! + the number of species
  integer, parameter :: NEQ = nspec

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

  subroutine actual_burner(state_in, state_out, dt, time)    
    
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
  
    double precision :: rpar(n_rpar_comps)
    integer :: ipar

    integer :: j

    EXTERNAL jac, f_rhs

    ! set the tolerances. 
    atol(1:nspec) = 1.d-12    ! mass fractions

    rtol(1:nspec) = 1.d-12    ! mass fractions

    do j = 1, state_in % N

       ! we want VODE to re-initialize each time we call it
       istate = 1

       rwork(:) = ZERO
       iwork(:) = 0


       ! set the maximum number of steps allowed (the VODE default is 500)
       iwork(6) = 15000

       ! initialize the integration time
       integration_time = ZERO

       y(:) = state_in % xn(j,:)

       rpar(irp_dens) = state_in % rho(j)
       rpar(irp_temp) = state_in % T(j)

       ! call the integration routine
       call dvode(f_rhs, NEQ, y, integration_time, dt, ITOL, rtol, atol, ITASK, &
                  istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_NUMERICAL_JAC, &
                  rpar, ipar)

       if (istate < 0) then
          print *, 'ERROR: integration failed in net'
          print *, 'istate = ', istate
          print *, 'time = ', integration_time
          call bl_error("ERROR in burner: integration failed")
       endif

       ! store the new mass fractions -- make sure that they are positive
       state_out % xn(j,ifuel_)= max(y(ifuel_), ZERO)
       state_out % xn(j,iash_) = min(y(iash_), ONE)

       ! compute the energy release from the change in fuel mass fractions.
       enuc = -specific_q_burn*(state_out % xn(j,ifuel_) - state_in % xn(j,ifuel_))

       state_out % e(j) = state_in % e(j) + enuc

       if (verbose) then

          ! print out some integration statistics, if desired
          print *, 'integration summary: '
          print *, 'dens: ', state_in % rho(j), ' temp: ', state_in % T(j)
          print *, 'number of steps taken: ', iwork(11)
          print *, 'number of f evaluations: ', iwork(12)
       endif

    enddo

  end subroutine actual_burner

end module actual_burner_module
