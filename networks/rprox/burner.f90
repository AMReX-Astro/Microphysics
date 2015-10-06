! This module contains a version of Wallace and Woosley's (ApJS 45,389
! (1981)) rprox reaction network burner.
!
!

module burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module
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

    implicit none

    real(kind=dp_t), intent(in   ) :: dens, temp, Xin(nspec), dt
    real(kind=dp_t), intent(  out) :: Xout(nspec), rho_omegadot(nspec), rho_Hnuc
    
    type(eos_t) :: eos_state

    integer :: n
    real(kind=dp_t) :: enuc, dX, t9
    real(kind=dp_t), parameter :: T2T9 = 1.0e-9_dp_t

    logical, parameter :: verbose = .false.

    ! set the number of independent variables -- this should be temperature(t9)
    ! + the number of species which participate in the evolution equations
    integer, parameter :: NEQ = nspec + 1

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

    EXTERNAL jac, f_rhs
    
    integer :: i

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
    
    ! set the parameters regarding how often to re-evaluate the 
    ! thermodynamics.  T_eos will always store the temperature
    ! that was used for the last EOS call.  dT_crit is the 
    ! relative temperature change beyond which we need to re-evaluate
    ! the thermodynamics
    !
    ! **NOTE** if the burner is not converging (and the temperatures
    ! are shooting up to unrealistically high values), you likely
    ! need to reduce dT_crit to ensure more frequent EOS calls.
!    T_eos = temp
!    dT_crit = 1.e20_dp_t


    ! set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.  
    !
    ! **NOTE** if you reduce these tolerances, you probably will need
    ! to (a) decrease dT_crit, (b) increase the maximum number of 
    ! steps allowed.
    atol(1:nspec) = 1.e-11_dp_t    ! mass fractions
    atol(nspec+1) = 1.e-8_dp_t     ! temperature

    rtol(1:nspec) = 1.e-12_dp_t    ! mass fractions
    rtol(nspec+1) = 1.e-8_dp_t     ! temperature

    ! we want VODE to re-initialize each time we call it
    istate = 1

    rwork(:) = ZERO
    iwork(:) = 0


    ! set the maximum number of steps allowed (the VODE default is 500)
    iwork(6) = 15000


    ! initialize the integration time
    time = ZERO

    ! rates and things are in terms of t9 and this is more on same
    ! scale as species
    t9 = temp*T2T9

    ! abundances are the first nspec values and temperature (t9)is the last
    y(1:nspec) = Xin(:) / aion(:)
    y(nspec+1) = t9

    ! set the thermodynamics that are passed via rpar to the RHS routine
    rpar(irp_dens) = dens
    rpar(irp_T9_eos) = t9
    rpar(irp_dTcrit) = 1e15_dp_t

    ! only burn if 0.2 < t9 < 2.5 or X(h1) > 0.05
    ! the last restriction is a kludge based on the last paragraph of WW81
    if ((t9 .gt. 0.2_dp_t .and. t9 .lt. 2.5_dp_t) .or. Xin(ih1) > 0.05d0) then

       ! call the EOS to get cp and dhdX for temperature evolution
       eos_state%rho = dens
       eos_state%T   = temp
       eos_state%xn  = Xin

       call eos(eos_input_rt, eos_state, .false.)

       rpar(irp_cp) = eos_state%cp
       rpar(irp_dhdX:irp_dhdX+nspec-1) = eos_state%dhdX

       ! call the integration routine
       call dvode(f_rhs, NEQ, y, time, dt, ITOL, rtol, atol, ITASK, &
            istate, IOPT, rwork, LRW, iwork, LIW, jac, &
            !MF_NUMERICAL_JAC, &
            MF_ANALYTIC_JAC,&
            rpar, ipar)

       if (istate < 0) then
          print *, 'ERROR: integration failed in net'
          print *, 'initial T   = ', temp
          print *, 'initial rho = ', dens
          print *, 'iniitial X  = ', Xin
          print *, 'istate = ', istate
          print *, 'time = ', time
          print *, 'dt = ', dt
          print *, 'output state = ', y
          call bl_error("ERROR in burner: integration failed")
       endif

       ! update composition
       Xout(1:nspec) = y(1:nspec) * aion

       ! normalize
       dX = ZERO
       do n = 1, nspec
          Xout(n) = max(ZERO,Xout(n))
          dX = dX + Xout(n)
       end do
       Xout = Xout / dX

       enuc = -sum((Xout - Xin)*ebin)
       rho_Hnuc = dens*enuc/dt
       rho_omegadot = dens*(Xout - Xin)/dt

       if (verbose) then

          ! print out some integration statistics, if desired
          print *, 'integration summary: '
          print *, 'dens: ', dens, ' temp: ', temp
          print *, 'number of steps taken: ', iwork(11)
          print *, 'number of f evaluations: ', iwork(12)
       endif


    else ! just pass through
       Xout = Xin
       rho_Hnuc = ZERO
       rho_omegadot = ZERO
    endif


  end subroutine burner


end module burner_module
