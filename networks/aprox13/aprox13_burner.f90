! This module contains the triple-alpha reaction network burner.
!
! Given the initial state of the system and the right-hand-side (f_rhs.f90) of
! a linear system of ODEs, this routine calls VODE to get the updated 
! composition, enthalpy and rho_omegadot.
! The temperature is evolved (in an isobaric, self-heating formalism) 
! concurrently with the species to prevent instabilities - see 
! Muller A&A 162, 103-108 (1986) for a disussion on the formation of 
! instabilities.  However, at the end of the call, the changes to temperature 
! are discarded.
!
! This burner provides an explicit Jacobian matrix to the DVODE solver.

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

  ! Tolerance parameters:
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

  ! We want to do a normal computation, and get the output values of y(t)
  ! after stepping though dt.
  integer, PARAMETER :: ITASK = 1

  ! We will override the maximum number of steps, so turn on the 
  ! optional arguments flag
  integer, parameter :: IOPT = 1

  ! Declare a real work array of size 22 + 9*NEQ + 2*NEQ**2 and an
  ! integer work array of since 30 + NEQ
  integer, parameter :: LRW = 22 + 9*NEQ + 2*NEQ**2

  integer, parameter :: LIW = 30 + NEQ

contains

  subroutine actual_burner(state_in, state_out, dt, time)

    use rpar_indices
    use network_indices

    implicit none

    type (eos_t_vector), intent(in   ) :: state_in
    type (eos_t_vector), intent(inout) :: state_out
    double precision,    intent(in   ) :: dt, time
    
    integer :: j

    double precision :: enuc, dX, local_time

    logical, parameter :: verbose = .false.

    ! Allocate storage for the input state

    double precision, dimension(NEQ) :: y
    
    double precision, dimension(NEQ) :: atol, rtol

    double precision, dimension(LRW) :: rwork  

    integer, dimension(LIW) :: iwork

    ! istate determines the state of the calculation.  A value of 1 meeans
    ! this is the first call to the problem -- this is what we will want.
    ! Note, istate is changed over the course of the calculation, so it
    ! cannot be a parameter
    integer :: istate
    
    double precision, allocatable :: rpar(:)
    integer :: ipar

    double precision :: sum

    EXTERNAL jac, f_rhs
    
    ! Allocate storage for rpar -- the scratch array passed into the
    ! rhs and jacobian routines
    allocate(rpar(n_rpar_comps))

    ! Set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.  
    !
    ! **NOTE** if you reduce these tolerances, you probably will need
    ! to (a) decrease dT_crit, (b) increase the maximum number of 
    ! steps allowed.
    atol(1:nspec) = 1.e-12_dp_t    ! mass fractions
    atol(net_itemp) = 1.e-8_dp_t   ! temperature
    atol(net_ienuc) = 1.e-8_dp_t   ! energy generated
    
    rtol(1:nspec) = 1.e-12_dp_t    ! mass fractions
    rtol(net_itemp) = 1.e-5_dp_t   ! temperature
    rtol(net_ienuc) = 1.e-10_dp_t  ! energy generated

    do j = 1, state_in % N

        ! We want VODE to re-initialize each time we call it
       istate = 1

       rwork(:) = ZERO
       iwork(:) = 0

       ! Set the maximum number of steps allowed (the VODE default is 500)
       iwork(6) = 150000

       ! Initialize the integration time
       local_time = ZERO

       ! Abundances are the first nspec values and temperature and energy are the last
       y(1:nspec)   = state_in % xn(j,:)
       y(net_itemp) = state_in % T(j)
       y(net_ienuc) = ZERO

       ! Density, specific heat at constant pressure, c_p, and dhdX are needed
       ! in the righthand side routine, so we will pass these in to those routines
       ! via the rpar functionality in VODE.

       rpar(irp_dens) = state_in % rho(j)
       rpar(irp_cp)   = state_in % cp(j)
       rpar(irp_cv)   = state_in % cv(j)

       ! The following thermodynamic derivatives are calculated in the EOS.

       ! dhdX = dh/dX |_{p,T}
       rpar(irp_dhdX:irp_dhdX-1+nspec) = state_in % dhdX(j,:)

       ! dedX = de/dX |_{rho, T}
       rpar(irp_dedX:irp_dedX-1+nspec) = state_in % dEdX(j,:)

       ! This is just used to make sure everything is happy
       rpar(irp_smallx) = smallx



       ! Call the integration routine
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
       state_out % xn(j,:) = max(smallx, min(ONE, y(1:nspec)))    


       ! Energy was integrated in the system -- we use this integrated
       ! energy which contains both the reaction energy release and
       ! neutrino losses. The final energy is the initial energy
       ! plus this energy release. Note that we get a new temperature too,
       ! but we will discard it and the main burner module will do an EOS
       ! call to get a final temperature consistent with this new energy.

       state_out % e(j) = state_in % e(j) + y(net_ienuc)

       if (verbose) then
          ! Print out some integration statistics, if desired
          print *, 'integration summary: '
          print *, 'dens: ', state_out % rho(j), ' temp: ', state_out % T(j)
          print *, 'number of steps taken: ', iwork(11)
          print *, 'number of f evaluations: ', iwork(12)
       endif

    enddo

  end subroutine actual_burner

end module actual_burner_module
