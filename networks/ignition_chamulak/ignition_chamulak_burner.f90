module actual_burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module, only: eos_input_rt, eos
  use eos_type_module
  use network

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

  subroutine actual_burner(state_in, state_out, dt, time)

    use rpar_indices
    
    implicit none

    type (eos_t_vector) :: state_in, state_out
    double precision    :: dt, time

    double precision :: enuc, local_time

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



    ic12 = network_species_index("carbon-12")
    io16 = network_species_index("oxygen-16")
    iash = network_species_index("ash")


    
    ! set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.
    atol(1:nspec_advance) = 1.d-12    ! mass fractions
    atol(nspec_advance+1) = 1.d-8     ! temperature

    rtol(1:nspec_advance) = 1.d-12    ! mass fractions
    rtol(nspec_advance+1) = 1.d-5     ! temperature

    do j = 1, state_in % N

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
       
       y(1:nspec) = state_in % xn(j,:)
       y(nspec+1) = state_in % T(j)

       rpar(irp_dens) = state_in % rho(j)
       rpar(irp_cp)   = state_in % cp(j)
       rpar(irp_dhdX:irp_dhdX-1+nspec) = state_in % dhdX(j,:)
       rpar(irp_o16)  = state_in % xn(j,iO16)



       ! Call the integration routine.
       
       call dvode(f_rhs, NEQ, y, local_time, dt, ITOL, rtol, atol, ITASK, &
            istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_NUMERICAL_JAC,&
            rpar, ipar)



       if (istate < 0) then
          print *, 'ERROR: integration failed in net'
          print *, 'istate = ', istate
          print *, 'time = ', local_time
          call bl_error("ERROR in burner: integration failed")
       endif

       ! Store the new mass fractions.
       state_out % xn(j,iC12) = max(y(iC12), ZERO)
       state_out % xn(j,iO16) = state_in % xn(j,iO16)
       state_out % xn(j,iash) = ONE - state_out % xn(j,iC12) - state_out % xn(j,iO16)

       ! Energy was integrated in the system -- we use this integrated
       ! energy which contains both the reaction energy release and
       ! neutrino losses. The final energy is the initial energy
       ! plus this energy release. Note that we get a new temperature too,
       ! but we will discard it and the main burner module will do an EOS
       ! call to get a final temperature consistent with this new energy.

       enuc = get_ebin_value(dens) * (state_out % xn(j,iC12) - state_in % xn(j,iC12))

       state_out % e(j) = state_in % e(j) + enuc

       if (verbose) then

          ! Print out some integration statistics, if desired.
          
          print *, 'integration summary: '
          print *, 'dens: ', state_out % rho(j), ' temp: ', state_out % T(j)
          print *, 'number of steps taken: ', iwork(11)
          print *, 'number of f evaluations: ', iwork(12)
          
       endif

    enddo
    
    call bl_deallocate(y)
    call bl_deallocate(atol)
    call bl_deallocate(rtol)
    call bl_deallocate(rwork)
    call bl_deallocate(iwork)
    call bl_deallocate(rpar)

  end subroutine actual_burner

end module actual_burner_module
