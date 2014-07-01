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
!

module burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module, only: eos_input_rt, eos
  use eos_type_module
  use network

  private
  public :: burner
  
contains

  subroutine burner(dens, temp, Xin, dt, Xout, rho_omegadot, rho_Hnuc)

    ! outputs:
    !   Xout are the mass fractions after burning through timestep dt
    !   rho_omegadot = rho dX/dt
    !   rho_Hnuc = - sum_k q_k rho_omegadot_k  [erg / cm^3 / s]

    use rpar_indices
    use network_indices

    implicit none

    real(kind=dp_t), intent(in   ) :: dens, temp, Xin(nspec), dt
    real(kind=dp_t), intent(  out) :: Xout(nspec), rho_omegadot(nspec), rho_Hnuc
  
    integer :: n
    real(kind=dp_t) :: enuc, dX

    logical, parameter :: verbose = .false.

    ! set the number of independent variables -- this should be
    ! temperature, enuc + the number of species which participate
    ! in the evolution equations
    ! 
    integer, parameter :: NEQ = 2 + nspec

    ! allocate storage for the input state
    real(kind=dp_t), dimension(NEQ) :: y

    ! our problem is stiff, tell ODEPACK that. 21 means stiff, jacobian 
    ! function is supplied, 22 means stiff, figure out my jacobian through 
    ! differencing
    integer, parameter :: MF_ANALYTIC_JAC = 21, MF_NUMERICAL_JAC = 22


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
    real(kind=dp_t), dimension(NEQ) :: atol, rtol


    real(kind=dp_t) :: time
    

    ! we want to do a normal computation, and get the output values of y(t)
    ! after stepping though dt
    integer, PARAMETER :: ITASK = 1
  

    ! istate determines the state of the calculation.  A value of 1 meeans
    ! this is the first call to the problem -- this is what we will want.
    ! Note, istate is changed over the course of the calculation, so it
    ! cannot be a parameter
    integer :: istate


    ! we will override the maximum number of steps, so turn on the 
    ! optional arguments flag
    integer, parameter :: IOPT = 1
    
    ! declare a real work array of size 22 + 9*NEQ + 2*NEQ**2 and an
    ! integer work array of since 30 + NEQ
    integer, parameter :: LRW = 22 + 9*NEQ + 2*NEQ**2
    real(kind=dp_t), dimension(LRW) :: rwork
    
    integer, parameter :: LIW = 30 + NEQ
    integer, dimension(LIW) :: iwork
    
    real(kind=dp_t), allocatable :: rpar(:)
    integer :: ipar

    type (eos_t) :: eos_state

    real(kind=dp_t) :: sum

    EXTERNAL jac, f_rhs
    
    integer :: i

    real(kind=dp_t) :: smallx=1.e-12

    logical, save :: firstCall = .true.

    if (firstCall) then

       if (.NOT. network_initialized) then
          call bl_error("ERROR in burner: must initialize network first")
       endif

       
       firstCall = .false.
    endif

    ! allocate storage for rpar -- the scratch array passed into the
    ! rhs and jacobian routines
    allocate(rpar(n_rpar_comps))
    

    ! set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.  
    !
    ! Note: VODE computes the tolerance to compare against for
    ! convergence as (quoting vode.f):
    !
    !      EWT(i) = RTOL*abs(Y(i)) + ATOL(i)  if ITOL = 2.
    !      Thus the local error test passes if, in each component,
    !      either the absolute error is less than ATOL (or ATOL(i)),
    !      or the relative error is less than RTOL.
    !      Use RTOL = 0.0 for pure absolute error control, and
    !      use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
    !      control.  Caution: Actual (global) errors may exceed these
    !      local tolerances, so choose them conservatively.
    !
    atol(1:nspec) = 1.e-10_dp_t    ! mass fractions
    atol(itemp) = 1.e-8_dp_t       ! temperature
    atol(ienuc) = 1.e-8_dp_t       ! energy generated

    rtol(1:nspec) = 1.e-10_dp_t    ! mass fractions
    rtol(itemp) = 1.e-5_dp_t       ! temperature
    rtol(ienuc) = 1.e-6_dp_t      ! energy generated

    ! we want VODE to re-initialize each time we call it
    istate = 1

    rwork(:) = ZERO
    iwork(:) = 0


    ! set the maximum number of steps allowed (the VODE default is 500)
    iwork(6) = 15000


    ! initialize the integration time
    time = ZERO

    ! abundances are the first nspec-1 values and temperature is the last
    y(1:nspec) = Xin(:)
    y(itemp)   = temp
    y(ienuc)   = ZERO
    

    ! set the thermodynamics that are passed via rpar to the RHS routine
    ! we need the specific heat at constant pressure and dhdX |_p.  Take
    ! T, rho, Xin as input
    eos_state%rho   = dens
    eos_state%T     = temp
    eos_state%xn(:) = Xin(:)

    call eos(eos_input_rt, eos_state, .false.)


    ! density, specific heat at constant pressure, c_p, and dhdX are needed
    ! in the righthand side routine, so we will pass these in to those routines
    ! via the rpar functionality in VODE.

    rpar(irp_dens) = dens
    rpar(irp_cp)   = eos_state%cp
    rpar(irp_dhdX:irp_dhdX-1+nspec) = eos_state%dhdX(:)

    ! this is just used to make sure everything is happy
    rpar(irp_smallx) = smallx


    ! call the integration routine
    call dvode(f_rhs, NEQ, y, time, dt, ITOL, rtol, atol, ITASK, &
         istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_NUMERICAL_JAC,&
         rpar, ipar)

    if (istate < 0) then
       print *, 'ERROR: integration failed in net'
       print *, 'istate = ', istate
       print *, 'time = ', time
       call bl_error("ERROR in burner: integration failed")
    endif


    ! store the new mass fractions -- note, we discard the temperature
    ! here.  Make sure that they are positive
    do n = 1, nspec
       Xout(n)  = max(y(n), smallx)
    enddo


    ! enforce sum{X_k} = 1
    sum = ZERO
    do n = 1, nspec
       Xout(n) = max(smallx, min(ONE, Xout(n)))
       sum = sum + Xout(n)
    enddo
    Xout(:) = Xout(:)/sum

    ! compute the energy release.  Our convention is that the binding
    ! energies are negative, so the energy release is
    ! - sum_k { (Xout(k) - Xin(k)) ebin(k) }
    !

    rho_omegadot(:) = ZERO
    do n = 1, nspec
       dX = Xout(n) - Xin(n)
       rho_omegadot(n) = dens * dX / dt
    enddo

    ! energy was integrated in the system -- we use this integrated
    ! energy which contains both the reaction energy release and
    ! neutrino losses
    rho_Hnuc = dens*y(ienuc)/dt

    if (verbose) then

       ! print out some integration statistics, if desired
       print *, 'integration summary: '
       print *, 'dens: ', dens, ' temp: ', temp
       print *, 'number of steps taken: ', iwork(11)
       print *, 'number of f evaluations: ', iwork(12)
    endif

  end subroutine burner

end module burner_module
