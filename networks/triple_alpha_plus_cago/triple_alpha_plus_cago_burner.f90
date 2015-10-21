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

module actual_burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module, only: eos_input_rt, eos
  use eos_type_module
  use network

  private
  public :: burner, get_enuc_T_sensitivity
  
contains

  subroutine actual_burner(dens, temp, Xin, dt, Xout, rho_omegadot, rho_Hnuc)

    ! outputs:
    !   Xout are the mass fractions after burning through timestep dt
    !   rho_omegadot = rho dX/dt
    !   rho_Hnuc = - sum_k q_k rho_omegadot_k  [erg / cm^3 / s]

    use rpar_indices

    implicit none

    real(kind=dp_t), intent(in   ) :: dens, temp, Xin(nspec), dt
    real(kind=dp_t), intent(  out) :: Xout(nspec), rho_omegadot(nspec), rho_Hnuc
  
    real(kind=dp_t) :: dT_crit, T_eos

    integer :: n
    real(kind=dp_t) :: enuc, dX

    logical, parameter :: verbose = .false.

    ! set the number of independent variables -- this should be temperature
    ! + the number of species which participate in the evolution equations
    ! 
    integer, parameter :: NEQ = 1 + nevolve

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

    integer, save :: ihe4, ic12, io16, ife56
    integer, save :: ir3a, ircago

    
    real(kind=dp_t), allocatable :: rpar(:)
    integer :: ipar

    type (eos_t) :: eos_state

    real(kind=dp_t) :: sum

    EXTERNAL jac, f_rhs
    
    integer :: i

    logical, save :: firstCall = .true.

    if (firstCall) then

       if (.NOT. network_initialized) then
          call bl_error("ERROR in burner: must initialize network first")
       endif

       ihe4  = network_species_index("helium-4")
       ic12  = network_species_index("carbon-12")
       io16  = network_species_index("oxygen-16")
       ife56 = network_species_index("iron-56")

       ir3a   = network_reaction_index("3agc")
       ircago = network_reaction_index("cago")
       

       if (ihe4 < 0 .or. ic12 < 0 .or. io16 < 0 .or. ife56 < 0 .or. &
           ir3a < 0 .or. ircago < 0) then
          call bl_error("ERROR in burner: species undefined")
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
    T_eos = temp
    dT_crit = 1.e20_dp_t


    ! set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.  
    !
    ! **NOTE** if you reduce these tolerances, you probably will need
    ! to (a) decrease dT_crit, (b) increase the maximum number of 
    ! steps allowed.
    atol(1:nevolve) = 1.e-12_dp_t    ! mass fractions
    atol(NEQ) = 1.e-8_dp_t         ! temperature

    rtol(1:nevolve) = 1.e-12_dp_t    ! mass fractions
    rtol(NEQ) = 1.e-5_dp_t         ! temperature

    ! we want VODE to re-initialize each time we call it
    istate = 1

    rwork(:) = ZERO
    iwork(:) = 0


    ! set the maximum number of steps allowed (the VODE default is 500)
    iwork(6) = 15000


    ! initialize the integration time
    time = ZERO

    ! abundances are the first nspec-1 values and temperature is the last
    y(ihe4)  = Xin(ihe4)
    y(ic12)  = Xin(ic12)
    y(io16)  = Xin(io16)
    y(NEQ)   = temp

    ! set the thermodynamics that are passed via rpar to the RHS routine--
    ! these will be updated in f_rhs if the relative temperature change 
    ! exceeds dT_crit

    ! we need the specific heat at constant pressure and dhdX |_p.  Take
    ! T, rho, Xin as input
    eos_state%rho   = dens
    eos_state%T     = temp
    eos_state%xn(:) = Xin(:)

    call eos(eos_input_rt, eos_state)


    ! density, specific heat at constant pressure, c_p, and dhdX are needed
    ! in the righthand side routine, so we will pass these in to those routines
    ! via the rpar functionality in VODE.
    !
    ! Since evaluating the EOS is expensive, we don't call it for every RHS
    ! call -- T_eos and dT_crit control the frequency (see above)
    rpar(irp_dens) = dens
    rpar(irp_cp)   = eos_state%cp
    rpar(irp_dhdX:irp_dhdX-1+nspec) = eos_state%dhdX(:)
    rpar(irp_Teos) = T_eos
    rpar(irp_Tcrit) = dT_crit
    rpar(irp_Y56) = Xin(ife56)/aion(ife56)


    ! call the integration routine
    call dvode(f_rhs, NEQ, y, time, dt, ITOL, rtol, atol, ITASK, &
         istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_ANALYTIC_JAC,&
         rpar, ipar)

    if (istate < 0) then
       print *, 'ERROR: integration failed in net'
       print *, 'istate = ', istate
       print *, 'time = ', time
       call bl_error("ERROR in burner: integration failed")
    endif


    ! store the new mass fractions -- note, we discard the temperature
    ! here and instead compute the energy release from the binding
    ! energy -- make sure that they are positive
    ! 
    Xout(ihe4)  = max(y(ihe4), ZERO)
    Xout(ic12)  = max(y(ic12), ZERO)
    Xout(io16)  = max(y(io16), ZERO)
    Xout(ife56) = Xin(ife56)


    ! enforce sum{X_k} = 1
    sum = ZERO
    do n = 1, nspec
       Xout(n) = max(ZERO, min(ONE, Xout(n)))
       sum = sum + Xout(n)
    enddo
    Xout(:) = Xout(:)/sum

    ! compute the energy release.  Our convention is that the binding
    ! energies are negative, so the energy release is
    ! - sum_k { (Xout(k) - Xin(k)) ebin(k) }
    !
    enuc = ZERO
    rho_omegadot(:) = ZERO
    do n = 1, nspec
       ! Fe56 doesn't participate in reactions - force this
       if (n == ife56) continue

       dX = Xout(n) - Xin(n)

       enuc = enuc - dX * ebin(n)

       rho_omegadot(n) = dens * dX / dt
    enddo

    rho_Hnuc = dens*enuc/dt

    if (verbose) then

       ! print out some integration statistics, if desired
       print *, 'integration summary: '
       print *, 'dens: ', dens, ' temp: ', temp
       print *, 'number of steps taken: ', iwork(11)
       print *, 'number of f evaluations: ', iwork(12)
    endif


  end subroutine actual_burner




  subroutine get_enuc_T_sensitivity(dens, temp, X, denucdT)
    
    ! Calculate the energy generation rate's temperature sensitivity
    ! Used for diagnostic purposes only

    use rates_module
    use screen_module
    use dydt_module

    implicit none

    real(kind=dp_t), intent(IN   ) :: dens, temp, X(nspec)
    real(kind=dp_t), intent(  OUT) :: denucdT

    real(kind=dp_t) :: ymol(nspec)
    real(kind=dp_t) :: rates(nrat), dratesdt(nrat)
    real(kind=dp_t) :: dXdotdT(nspec)
    integer :: k

    ! get the indices
    if (.not. network_initialized) then
       call bl_error("ERROR in get_enuc_T_sensitivity: must initialize network first")
    endif

    ! calculate ymol
    ymol = X / aion

    ! get the d/dT(dX/dt) info, dydt(dratesdT) gives us this
    call make_rates(temp, dens, rates, dratesdt)
    call screen(temp, dens, ymol, rates, dratesdt)
    call dydt(ymol, dratesdt, dXdotdT)

    ! calculate temperature sensitivity
    denucdT = - sum(dXdotdT*ebin)
    

  end subroutine get_enuc_T_sensitivity

end module actual_burner_module
