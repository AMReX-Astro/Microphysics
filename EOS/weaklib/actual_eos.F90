module actual_eos_module

  use amrex_fort_module, only: rt => amrex_real

  use UnitsModule, only: Gram, Centimeter, Kelvin
  use wlEquationOfStateTableModule, only: EquationOfStateTableType
  use wlEOSInversionModule, only: ComputeTemperatureWith_DEY
  use wlInterpolationModule, only: LogInterpolateSingleVariable, &
                                   LogInterpolateDifferentiateSingleVariable, &
                                   ComputeTempFromIntEnergy_Lookup, &
                                   ComputeTempFromIntEnergy_Bisection, &
                                   ComputeTempFromIntEnergy_Secant, &
                                   ComputeTempFromPressure, &
                                   ComputeTempFromPressure_Bisection

  character (len=64), public :: eos_name = "weaklib"
  type (EquationOfStateTableType), target, public :: eos_table
  type (EquationOfStateTableType), pointer, public :: eos_pointer

  public actual_eos, actual_eos_init, actual_eos_finalize, eos_supports_input_type

  private

  integer :: &
       iD_T, iT_T, iY_T, &
       iP_T, iS_T, iE_T, iMe_T, iMp_T, iMn_T, &
       iXp_T, iXn_T, iXa_T, iXh_T, iGm_T
  integer, dimension(3) :: &
       LogInterp
  real(rt) :: &
       OS_P, OS_S, OS_E, OS_Me, OS_Mp, OS_Mn, &
       OS_Xp, OS_Xn, OS_Xa, OS_Xh, OS_Gm
  real(rt), dimension(:), allocatable :: &
       Ds_T, Ts_T, Ys_T
  real(rt), dimension(:,:,:), allocatable :: &
    Ps_T, Ss_T, Es_T, Mes_T, Mps_T, Mns_T, &
    Xps_T, Xns_T, Xas_T, Xhs_T, Gms_T

#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET &
  !$OMP ( Ds_T, Ts_T, Ys_T, &
  !$OMP   OS_P, OS_S, OS_E, OS_Me, OS_Mp, OS_Mn, OS_Xp, OS_Xn, OS_Xa, OS_Xh, OS_Gm, &
  !$OMP   Ps_T, Ss_T, Es_T, Mes_T, Mps_T, Mns_T, Xps_T, Xns_T, Xas_T, Xhs_T, Gms_T )
#elif defined(WEAKLIB_OACC)
  !$ACC DECLARE CREATE &
  !$ACC ( Ds_T, Ts_T, Ys_T, &
  !$ACC   OS_P, OS_S, OS_E, OS_Me, OS_Mp, OS_Mn, OS_Xp, OS_Xn, OS_Xa, OS_Xh, OS_Gm, &
  !$ACC   Ps_T, Ss_T, Es_T, Mes_T, Mps_T, Mns_T, Xps_T, Xns_T, Xas_T, Xhs_T, Gms_T )
#endif

contains


  function eos_supports_input_type(input) result(supported)

    use eos_type_module

    implicit none

    integer, intent(in) :: input
    logical :: supported

    if (input == eos_input_rt .or. &
        input == eos_input_rp .or. &
        input == eos_input_re) then

       supported = .true.

    else

       supported = .false.

    endif

  end function eos_supports_input_type



  subroutine actual_eos(input, state)

    ! Weaklib EOS
    !
    ! The weaklib tables are indexed by log(density),
    ! log(temperature), and electron fraction.  As such, the usual
    ! 'composition' variable passed to the EOS is the electron fraction.
    !
    ! Make sure you use a network that uses ye as a species!

    use amrex_error_module, only: amrex_error
    use eos_type_module
    use weaklib_type_module, only: weaklib_eos_t, eos_to_weaklib, weaklib_to_eos

    implicit none

    ! Input arguments
    integer,      intent(in   ) :: input
    type (eos_t), intent(inout) :: state

    ! Local state
    type (weaklib_eos_t) :: weaklib_state

    call eos_to_weaklib(state, weaklib_state)

    select case (input)

       !---------------------------------------------------------------------------
       ! dens, temp, and ye are inputs;
       ! this is direct table interpolation, so nothing to do here
       !---------------------------------------------------------------------------

    case (eos_input_rt)

       continue


       !---------------------------------------------------------------------------
       ! dens, enthalpy, and xmass are inputs
       !---------------------------------------------------------------------------

    case (eos_input_rh)

       ! NOT CURRENTLY IMPLEMENTED
       call amrex_error("eos_input_th is not supported")


       !---------------------------------------------------------------------------
       ! temp, pres, and ye are inputs; iterate to find density
       !---------------------------------------------------------------------------

    case (eos_input_tp)

       ! NOT CURRENTLY IMPLEMENTED
       call amrex_error("eos_input_th is not supported")


       !---------------------------------------------------------------------------
       ! dens, pres, and ye are inputs; iterate to find the temperature
       !---------------------------------------------------------------------------

    case (eos_input_rp)

       call ComputeTemperatureFromPressure(weaklib_state)


       !---------------------------------------------------------------------------
       ! dens, energy, and ye are inputs; iterate to find temperature
       !---------------------------------------------------------------------------

    case (eos_input_re)

       call ComputeTemperatureFromSpecificInternalEnergy(weaklib_state)


       !---------------------------------------------------------------------------
       ! pres, entropy, and xmass are inputs
       !---------------------------------------------------------------------------

    case (eos_input_ps)
       ! NOT CURRENTLY IMPLEMENTED
       call amrex_error("eos_input_ps is not supported")


       !---------------------------------------------------------------------------
       ! pres, enthalpy, and xmass are inputs
       !---------------------------------------------------------------------------

    case (eos_input_ph)
       ! NOT CURRENTLY IMPLEMENTED
       call amrex_error("eos_input_ph is not supported")


       !---------------------------------------------------------------------------
       ! temp, enthalpy, and xmass are inputs
       !---------------------------------------------------------------------------

    case (eos_input_th)
       ! NOT CURRENTLY IMPLEMENTED
       call amrex_error("eos_input_th is not supported")


       !---------------------------------------------------------------------------
       ! The EOS input doesn't match any of the available options.
       !---------------------------------------------------------------------------

    case default

       call amrex_error("EOS: invalid input")

    end select

    ! Do a final lookup - by now we should have a consistent density and temperature
    call ApplyEquationOfState(weaklib_state)

    call weaklib_to_eos(weaklib_state, state)

  end subroutine actual_eos




  subroutine actual_eos_init

    use extern_probin_module, only: weaklib_eos_table_name
    use amrex_paralleldescriptor_module, only: amrex_pd_ioprocessor

    use wlIOModuleHDF, only: InitializeHDF, FinalizeHDF
    use wlEOSIOModuleHDF, only: ReadEquationOfStateTableHDF

    implicit none

    if (amrex_pd_ioprocessor()) then
       print *, ''
       print *, "Initializing Weaklib EOS on all MPI ranks from table file: ", trim(weaklib_eos_table_name)
       print *, ''
    endif

    eos_pointer => eos_table

    call InitializeHDF()
    call ReadEquationOfStateTableHDF(eos_table, trim(weaklib_eos_table_name))
    call FinalizeHDF()

    ! --- Thermodynamic State Indices ---

    iD_T = eos_table % TS % Indices % iRho
    iT_T = eos_table % TS % Indices % iT
    iY_T = eos_table % TS % Indices % iYe

    ! --- Thermodynamic States ---

    allocate( Ds_T(eos_table % TS % nPoints(iD_T)) )
    Ds_T = eos_table % TS % States(iD_T) % Values

    allocate( Ts_T(eos_table % TS % nPoints(iT_T)) )
    Ts_T = eos_table % TS % States(iT_T) % Values

    allocate( Ys_T(eos_table % TS % nPoints(iY_T)) )
    Ys_T = eos_table % TS % States(iY_T) % Values

    LogInterp = eos_table % TS % LogInterp

    ! --- Dependent Variables Indices ---

    iP_T  = eos_table % DV % Indices % iPressure
    iS_T  = eos_table % DV % Indices % iEntropyPerBaryon
    iE_T  = eos_table % DV % Indices % iInternalEnergyDensity
    iMe_T = eos_table % DV % Indices % iElectronChemicalPotential
    iMp_T = eos_table % DV % Indices % iProtonChemicalPotential
    iMn_T = eos_table % DV % Indices % iNeutronChemicalPotential
    iXp_T = eos_table % DV % Indices % iProtonMassFraction
    iXn_T = eos_table % DV % Indices % iNeutronMassFraction
    iXa_T = eos_table % DV % Indices % iAlphaMassFraction
    iXh_T = eos_table % DV % Indices % iHeavyMassFraction
    iGm_T = eos_table % DV % Indices % iGamma1

    ! --- Dependent Variables Offsets ---

    OS_P  = eos_table % DV % Offsets(iP_T)
    OS_S  = eos_table % DV % Offsets(iS_T)
    OS_E  = eos_table % DV % Offsets(iE_T)
    OS_Me = eos_table % DV % Offsets(iMe_T)
    OS_Mp = eos_table % DV % Offsets(iMp_T)
    OS_Mn = eos_table % DV % Offsets(iMn_T)
    OS_Xp = eos_table % DV % Offsets(iXp_T)
    OS_Xn = eos_table % DV % Offsets(iXn_T)
    OS_Xa = eos_table % DV % Offsets(iXa_T)
    OS_Xh = eos_table % DV % Offsets(iXh_T)
    OS_Gm = eos_table % DV % Offsets(iGm_T)

    ALLOCATE( Ps_T (1:eos_table % DV % nPoints(1), 1:eos_table % DV % nPoints(2), 1:eos_table % DV % nPoints(3)) )
    ALLOCATE( Ss_T (1:eos_table % DV % nPoints(1), 1:eos_table % DV % nPoints(2), 1:eos_table % DV % nPoints(3)) )
    ALLOCATE( Es_T (1:eos_table % DV % nPoints(1), 1:eos_table % DV % nPoints(2), 1:eos_table % DV % nPoints(3)) )
    ALLOCATE( Mes_T(1:eos_table % DV % nPoints(1), 1:eos_table % DV % nPoints(2), 1:eos_table % DV % nPoints(3)) )
    ALLOCATE( Mps_T(1:eos_table % DV % nPoints(1), 1:eos_table % DV % nPoints(2), 1:eos_table % DV % nPoints(3)) )
    ALLOCATE( Mns_T(1:eos_table % DV % nPoints(1), 1:eos_table % DV % nPoints(2), 1:eos_table % DV % nPoints(3)) )
    ALLOCATE( Xps_T(1:eos_table % DV % nPoints(1), 1:eos_table % DV % nPoints(2), 1:eos_table % DV % nPoints(3)) )
    ALLOCATE( Xns_T(1:eos_table % DV % nPoints(1), 1:eos_table % DV % nPoints(2), 1:eos_table % DV % nPoints(3)) )
    ALLOCATE( Xas_T(1:eos_table % DV % nPoints(1), 1:eos_table % DV % nPoints(2), 1:eos_table % DV % nPoints(3)) )
    ALLOCATE( Xhs_T(1:eos_table % DV % nPoints(1), 1:eos_table % DV % nPoints(2), 1:eos_table % DV % nPoints(3)) )
    ALLOCATE( Gms_T(1:eos_table % DV % nPoints(1), 1:eos_table % DV % nPoints(2), 1:eos_table % DV % nPoints(3)) )

    Ps_T (:,:,:) = eos_table % DV % Variables(iP_T ) % Values(:,:,:)
    Ss_T (:,:,:) = eos_table % DV % Variables(iS_T ) % Values(:,:,:)
    Es_T (:,:,:) = eos_table % DV % Variables(iE_T ) % Values(:,:,:)
    Mes_T(:,:,:) = eos_table % DV % Variables(iMe_T) % Values(:,:,:)
    Mps_T(:,:,:) = eos_table % DV % Variables(iMp_T) % Values(:,:,:)
    Mns_T(:,:,:) = eos_table % DV % Variables(iMn_T) % Values(:,:,:)
    Xps_T(:,:,:) = eos_table % DV % Variables(iXp_T) % Values(:,:,:)
    Xns_T(:,:,:) = eos_table % DV % Variables(iXn_T) % Values(:,:,:)
    Xas_T(:,:,:) = eos_table % DV % Variables(iXa_T) % Values(:,:,:)
    Xhs_T(:,:,:) = eos_table % DV % Variables(iXh_T) % Values(:,:,:)
    Gms_T(:,:,:) = eos_table % DV % Variables(iGm_T) % Values(:,:,:)

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET UPDATE TO &
    !$OMP ( OS_P, OS_S, OS_E, OS_Me, OS_Mp, OS_Mn, OS_Xp, OS_Xn, OS_Xa, OS_Xh, OS_Gm )

    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: Ds_T, Ts_T, Ys_T, &
    !$OMP          Ps_T, Ss_T, Es_T, Mes_T, Mps_T, Mns_T, Xps_T, Xns_T, Xas_T, Xhs_T, Gms_T )
#elif defined(WEAKLIB_OACC)
    !$ACC UPDATE DEVICE &
    !$ACC ( Ds_T, Ts_T, Ys_T, &
    !$ACC   OS_P, OS_S, OS_E, OS_Me, OS_Mp, OS_Mn, OS_Xp, OS_Xn, OS_Xa, OS_Xh, OS_Gm, &
    !$ACC   Ps_T, Ss_T, Es_T, Mes_T, Mps_T, Mns_T, Xps_T, Xns_T, Xas_T, Xhs_T, Gms_T )
#endif

  end subroutine actual_eos_init



  subroutine actual_eos_finalize

    implicit none

#if defined(WEAKLIB_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: Ds_T, Ts_T, Ys_T, &
    !$OMP               Ps_T, Ss_T, Es_T, Mes_T, Mps_T, Mns_T, Xps_T, Xns_T, Xas_T, Xhs_T, Gms_T )
#endif

    deallocate(Ds_T, Ts_T, Ys_T)

    deallocate( Ps_T  )
    deallocate( Ss_T  )
    deallocate( Es_T  )
    deallocate( Mes_T )
    deallocate( Mps_T )
    deallocate( Mns_T )
    deallocate( Xps_T )
    deallocate( Xns_T )
    deallocate( Xas_T )
    deallocate( Xhs_T )
    deallocate( Gms_T )

  end subroutine actual_eos_finalize



  subroutine ApplyEquationOfState(state)

    use weaklib_type_module, only: weaklib_eos_t

    implicit none

    type (weaklib_eos_t), intent(inout) :: state

    ! --- Interpolate Pressure ----------------------------------------

    call ComputeDependentVariable(state, state % pressure, Ps_T, OS_P)

    ! --- Interpolate Entropy Per Baryon ------------------------------

    call ComputeDependentVariable(state, state % entropy_per_baryon, Ss_T, OS_S)

    ! --- Interpolate Specific Internal Energy ------------------------

    call ComputeDependentVariable(state, state % specific_internal_energy, Es_T, OS_E)

    ! --- Interpolate Electron Chemical Potential ---------------------

    call ComputeDependentVariable(state, state % electron_chemical_potential, Mes_T, OS_Me)

    ! --- Interpolate Proton Chemical Potential -----------------------

    call ComputeDependentVariable(state, state % proton_chemical_potential, Mps_T, OS_Mp)

    ! --- Interpolate Neutron Chemical Potential ----------------------

    call ComputeDependentVariable(state, state % neutron_chemical_potential, Mns_T, OS_Mn)

    ! --- Interpolate Proton Mass Fraction ----------------------------

    call ComputeDependentVariable(state, state % proton_mass_fraction, Xps_T, OS_Xp)

    ! --- Interpolate Neutron Mass Fraction ---------------------------

    call ComputeDependentVariable(state, state % neutron_mass_fraction, Xns_T, OS_Xn)

    ! --- Interpolate Alpha Mass Fraction -----------------------------

    call ComputeDependentVariable(state, state % alpha_mass_fraction, Xas_T, OS_Xa)

    ! --- Interpolate Heavy Mass Fraction -----------------------------

    call ComputeDependentVariable(state, state % heavy_mass_fraction, Xhs_T, OS_Xh)

    ! --- Gamma1 ------------------------------------------------------

    call ComputeDependentVariable(state, state % gamma_one, Gms_T, OS_Gm)

  end subroutine ApplyEquationOfState


  subroutine ComputeTemperatureFromSpecificInternalEnergy(state)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    use weaklib_type_module, only: weaklib_eos_t

    implicit none

    type (weaklib_eos_t), intent(inout) :: state
    real(rt) :: d_p, e_p, y_p, t_lookup
    integer :: error

    d_p = state % density
    e_p = state % specific_internal_energy
    y_p = state % electron_fraction

    call ComputeTemperatureWith_DEY &
           ( d_p, e_p, y_p, Ds_T, Ts_T, Ys_T, Es_T, OS_E, t_lookup, &
             Error_Option = error )

    state % temperature = t_lookup

  end subroutine ComputeTemperatureFromSpecificInternalEnergy


  subroutine ComputeTemperatureFromPressure(state)

    use weaklib_type_module, only: weaklib_eos_t

    implicit none

    type (weaklib_eos_t), intent(inout) :: state
    real(rt) :: actual_temperature(1)

    call ComputeTempFromPressure_Bisection( &
         [state % density], [state % pressure], [state % electron_fraction], &
         Ds_T, Ts_T, Ys_T, LogInterp, &
         eos_table % DV % Variables(iP_T) % Values, OS_P, &
         actual_temperature)

    state % temperature = actual_temperature(1)

  end subroutine ComputeTemperatureFromPressure


  subroutine ComputeDependentVariable(state, variable, V_T, OS_V)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

    use weaklib_type_module, only: weaklib_eos_t

    implicit none

    type (weaklib_eos_t), intent(inout) :: state
    real(rt), intent(inout) :: variable
    real(rt), dimension(:,:,:), intent(in)  :: V_T
    real(rt),               intent(in)  :: OS_V
    real(rt) :: d_p, t_p, y_p, v_p

    d_p = state % density
    t_p = state % temperature
    y_p = state % electron_fraction

    CALL LogInterpolateSingleVariable &
           ( d_p, t_p, y_p, Ds_T, Ts_T, Ys_T, OS_V, V_T, v_p )

    variable = v_p

  end subroutine ComputeDependentVariable

end module actual_eos_module
