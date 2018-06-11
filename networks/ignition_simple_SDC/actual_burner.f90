module acutal_burner_module

  use amrex_constants_module
  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real
  use eos_module
  use network

  private
  public :: actual_burner, actual_burner_init

contains

  subroutine actual_burner_init()

  end subroutine actual_burner_init



  subroutine actual_burner(rhoXin, rhohin, dt, rhoout, rhoXout, rhohout, &
                           sdc_rhoX, sdc_rhoh, p0)
    ! outputs:
    !   rhoout is the updated density
    !   rhoXout are the updated density-weighted species (rho X_k)
    !   rhohout is the updated (rho h)

    use actual_rhs_module, only : f_rhs, jac
    use burner_aux_module, only : sdc_rhoX_pass, sdc_rhoh_pass, p0_pass

    implicit none

    real(rt), intent(in   ) :: rhoXin(nspec), rhohin, dt
    real(rt), intent(  out) :: rhoout, rhoXout(nspec), rhohout
    real(rt), intent(in   ) :: sdc_rhoX(nspec), sdc_rhoh
    real(rt), intent(in   ) :: p0

    integer :: n

    logical, parameter :: verbose = .false.

    ! set the number of independent variables -- this should be (rho h)
    ! + the number of species
    integer, parameter :: NEQ = 1 + nspec
  

    ! allocate storage for the input state
    real(rt), dimension(NEQ) :: y


    ! we will always refer to the species by integer indices that come from
    ! the network module -- this makes things robust to a shuffling of the 
    ! species ordering
    integer, save :: ic12, io16, img24

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
    real(rt), dimension(NEQ) :: atol, rtol


    real(rt) :: time
    

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
    real(rt), dimension(LRW) :: rwork
    
    integer, parameter :: LIW = 30 + NEQ
    integer, dimension(LIW) :: iwork
    
    real(rt) :: rpar
    integer :: ipar
    
    logical, save :: firstCall = .true.

    if (firstCall) then

       if (.NOT. network_initialized) then
          call amrex_error("ERROR in burner: must initialize network first")
       endif
     
       ic12 = network_species_index("carbon-12")
       io16 = network_species_index("oxygen-16")
       img24 = network_species_index("magnesium-24")
       
       if (ic12 < 0 .OR. io16 < 0 .OR. img24 < 0) then
          call amrex_error("ERROR in burner: species undefined")
       endif
       
       firstCall = .false.
    endif

    ! set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.  
    atol(1:nspec) = 1.d-12    ! density-weighted mass fractions
    atol(nspec+1) = 1.d-8     ! enthalpy
       
    rtol(1:nspec) = 1.d-12    ! density-weighted mass fractions
    rtol(nspec+1) = 1.d-8     ! enthalpy
    

    ! we want VODE to re-initialize each time we call it
    istate = 1
    
    rwork(:) = ZERO
    iwork(:) = 0
    
    
    ! set the maximum number of steps allowed (the VODE default is 500)
    iwork(6) = 15000
    
    
    ! initialize the integration time
    time = ZERO
    
    
    ! density-weighted abundances are the first nspec values and rhoh is the last
    y(ic12)  = rhoXin(ic12)
    y(io16)  = rhoXin(io16)
    y(img24) = rhoXin(img24)
    y(nspec+1) = rhohin

    ! sdc source terms (sdc_rhoX and sdc_rhoh) are needed
    ! in the righthand side routine, so we will pass these in through the
    ! burner_aux module.
    !
    sdc_rhoX_pass(:) = sdc_rhoX(:)
    sdc_rhoh_pass = sdc_rhoh
    p0_pass = p0

    ! call the integration routine
    call dvode(f_rhs, NEQ, y, time, dt, ITOL, rtol, atol, ITASK, &
               istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_NUMERICAL_JAC, &
               rpar, ipar)

    if (istate < 0) then
       print *, 'ERROR: integration failed in net'
       print *, 'istate = ', istate
       print *, 'time = ', time
       call amrex_error("ERROR in burner: integration failed")
    endif

    ! store the new mass fractions -- note, we discard the temperature
    ! here and instead compute the energy release from the binding
    ! energy -- make sure that they are positive
    rhoXout(ic12)  = max(y(ic12), ZERO)
    rhoXout(io16)  = max(y(io16), ZERO)
    rhoXout(img24) = max(y(img24),ZERO)

    rhoout = rhoXout(ic12) + rhoXout(io16) + rhoXout(img24)
        
    rhohout = y(nspec+1)

    if (verbose) then
       
       ! print out some integration statistics, if desired
       print *, 'integration summary: '
       print *, 'rhoout: ', rhoout
       print *, 'number of steps taken: ', iwork(11)
       print *, 'number of f evaluations: ', iwork(12)
    endif

  end subroutine burner

end module acutal_burner_module
