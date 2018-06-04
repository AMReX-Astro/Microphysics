module actual_rhs_module

  use bl_types
  use bl_constants_module
  use network
  use bl_error_module
  
  implicit none

contains

  subroutine actual_rhs_init()

    implicit none

  end subroutine actual_rhs_init

  subroutine f_rhs(n, t, y, ydot, rpar, ipar)

    use burner_aux_module, only : sdc_rhoX_pass, sdc_rhoh_pass, p0_pass

    implicit none

    ! our convention is that y(1:nspec) are the density-weighted species
    ! (in the same order as defined in network.f90, and y(nspec+1) is
    ! (rho h).
    integer :: n
    real(kind=dp_t) :: y(n), ydot(n)
    real(kind=dp_t) :: rpar
    integer :: ipar

    real(kind=dp_t) :: t

    real(kind=dp_t) :: rho_Hnuc

    integer, save :: ic12, io16, img24

    logical, save :: firstCall = .true.

    if (firstCall) then
       ic12 = network_species_index("carbon-12")
       io16 = network_species_index("oxygen-16")
       img24 = network_species_index("magnesium-24")

       firstCall = .false.
    end if

    ! get the RHS
    call f_rhs_instantaneous_reaction_rates(n, t, y, ydot, rho_Hnuc, rpar, ipar)

    ! now make ydots refer to rhoX and include the sdc sources
    ydot(ic12)  = ydot(ic12)  + sdc_rhoX_pass(ic12)
    ydot(io16)  = ydot(io16)  + sdc_rhoX_pass(io16)
    ydot(img24) = ydot(img24) + sdc_rhoX_pass(img24)
    ydot(nspec+1) = rho_Hnuc + sdc_rhoh_pass

    return

  end subroutine f_rhs

  subroutine jac(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)

    use bl_types
    use bl_constants_module
    use network

    ! we get the thermodynamic state through the burner_aux module -- we freeze
    ! these to the values are the top of the timestep to avoid costly
    ! EOS calls

    implicit none

    integer        , intent(IN   ) :: neq, ml, mu, nrpd, ipar
    real(kind=dp_t), intent(IN   ) :: y(neq), rpar, t
    real(kind=dp_t), intent(  OUT) :: pd(neq,neq)

    ! initialize
    pd(:,:)  = ZERO

    return
  end subroutine jac

  subroutine f_rhs_instantaneous_reaction_rates(n, t, y, ydot, rho_Hnuc, rpar, ipar)

    use eos_module, only: eos_input_rh, eos
    use eos_type_module
    use probin_module, only: use_tfromp
    use screening_module, only: screenz

    implicit none

    ! our convention is that y(1:nspec) are the density-weighted species
    ! (in the same order as defined in network.f90, and y(nspec+1) is
    ! (rho h).
    integer :: n
    real(kind=dp_t) :: y(n), ydot(n)

    integer :: k
    real(kind=dp_t) :: ymass(nspec)

    real(kind=dp_t) :: rpar
    integer :: ipar

    real(kind=dp_t) :: t

    real(kind=dp_t) :: dens, temp, T9, T9a
    real(kind=dp_t) :: rhoh, rho_Hnuc

    real(kind=dp_t) :: rate
    real(kind=dp_t) :: sc1212, dsc1212dt
    real(kind=dp_t) :: xc12tmp

    real(kind=dp_t), PARAMETER :: &
         one_twelvth = 1.0d0/12.0d0, &
         five_sixths = 5.0d0/ 6.0d0, &
         one_third = 1.0d0/ 3.0d0, &
         two_thirds = 2.0d0/ 3.0d0

    real(kind=dp_t) :: scratch

    real(kind=dp_t) :: a, b

    integer, save :: ic12, io16, img24

    real(kind=dp_t) :: X(nspec)

    logical, save :: firstCall = .true.

    type (eos_t) :: eos_state

    if (firstCall) then
       ic12 = network_species_index("carbon-12")
       io16 = network_species_index("oxygen-16")
       img24 = network_species_index("magnesium-24")

       firstCall = .false.
    end if

    ! define density and the mass fractions
    dens = y(1) + y(2) + y(3)

    X(1:nspec) = y(1:nspec)/dens

    ! enthalpy
    rhoh = y(nspec+1)

    ! compute the temperature from the EOS
    if (use_tfromp) then

       call bl_error("f_rhs_instantaneous_reaction_rates needs use_tfromp=F")

    else

       ! get T from (rho, h, X)
       eos_state%rho   = dens
       eos_state%xn(:) = X(1:nspec)
       eos_state%h     = rhoh/dens

       ! need an initial T guess
       eos_state%T = 1.d9

       call eos(eos_input_rh, eos_state)

       temp = eos_state%T

    endif


    ! call the screening routine
    ! compute the molar fractions -- needed for the screening
    ymass(ic12) = X(1) * aion_inv(ic12)
    ymass(io16) = X(2) * aion_inv(io16)
    ymass(img24) = X(3) * aion_inv(img24)

    call screenz(temp,dens,6.0d0,6.0d0,12.0d0,12.0d0,ymass,sc1212,dsc1212dt)


    ! compute some often used temperature constants
    T9     = temp/1.d9
    T9a    = T9/(1.0d0 + 0.0396d0*T9)

    ! compute the CF88 rate
    scratch    = T9a**one_third

    a       = 4.27d26*T9a**five_sixths*T9**(-1.5d0)
    b       = dexp(-84.165d0/scratch - 2.12d-3*T9*T9*T9)
    rate    = a *  b

    ! The change in number density of C12 is
    ! d(n12)/dt = - 2 * 1/2 (n12)**2 <sigma v>
    !
    ! where <sigma v> is the average of the relative velocity times the cross
    ! section for the reaction, and the factor accounting for the total number
    ! of particle pairs has a 1/2 because we are considering a reaction involving 
    ! identical particles (see Clayton p. 293).  Finally, the -2 means that for
    ! each reaction, we lose 2 carbon nuclei.
    !
    ! The corresponding Mg24 change is
    ! d(n24)/dt = + 1/2 (n12)**2 <sigma v>
    !
    ! note that no factor of 2 appears here, because we create only 1 Mg nuclei.
    !
    ! Switching over to mass fractions, using n = rho X N_A/A, where N_A is
    ! Avagadro's number, and A is the mass number of the nucleon, we get
    !
    ! d(X12)/dt = -2 *1/2 (X12)**2 rho N_A <sigma v> / A12
    !
    ! d(X24)/dt = + 1/2 (X12)**2 rho N_A <sigma v> (A24/A12**2)
    !
    ! these are equal and opposite.
    !
    ! The quantity [N_A <sigma v>] is what is tabulated in Caughlin and Fowler.

    ! we will always refer to the species by integer indices that come from
    ! the network module -- this makes things robust to a shuffling of the 
    ! species ordering


    ! changes in X due to reactions only
    xc12tmp = max(X(ic12),0.d0)
    ydot(ic12)  = -one_twelvth*dens*sc1212*rate*xc12tmp**2
    ydot(io16)  = 0.d0
    ydot(img24) =  one_twelvth*dens*sc1212*rate*xc12tmp**2

    ! compute the enthalpy source from reactions (note: eventually, we
    ! we will to add rho_Hext here too)
    rho_Hnuc = 0.0d0
    do k = 1, nspec
       rho_Hnuc = rho_Hnuc - ebin(k)*dens*ydot(k)
    enddo

    ! now make ydots refer to rhoX
    ydot(ic12)  = dens*ydot(ic12)
    ydot(io16)  = dens*ydot(io16)
    ydot(img24) = dens*ydot(img24)

    return

  end subroutine f_rhs_instantaneous_reaction_rates

end module actual_rhs_module
