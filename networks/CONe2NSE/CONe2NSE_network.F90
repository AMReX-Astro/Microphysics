module actual_network

  use fundamental_constants_module, only: m_e, m_n, n_A, c_light, m_p      
  
  integer, parameter :: nspec  = 6
  integer, parameter :: naux   = 11
  integer, parameter :: nrates = 0

  integer, parameter :: UFLAM = 1
  integer, parameter :: UFLDT = 2
  integer, parameter :: UFLSP = 3
  integer, parameter :: UCI   = 4
  integer, parameter :: UNEI  = 5
  integer, parameter :: UPHFA = 6
  integer, parameter :: UPHAQ = 7
  integer, parameter :: UPHQN = 8
  integer, parameter :: UYE   = 9
  integer, parameter :: UDYQN = 10
  integer, parameter :: UDQQN = 11

  integer, parameter :: iHe4 = 1
  integer, parameter :: iC12 = 2
  integer, parameter :: iO16 = 3
  integer, parameter :: iNe22 = 4
  integer, parameter :: iMg24 = 5
  integer, parameter :: iSi28 = 6
  
  double precision, save :: aion(nspec), zion(nspec), nion(nspec)
  double precision, save :: bion(nspec), yiion(nspec), yeion(nspec), qion(nspec)

  character (len=16), save :: ratenames(nrates)

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)
  
  character (len=32) :: network_name = "CONe2NSE"

  logical,           save  :: restart  

  integer, parameter :: infile_unit=2
  integer, parameter :: spec_num=238  ! number of lines in SpeciesList.txt file  

contains

  subroutine actual_network_init()

    use bl_constants_module, only: ONE
    use NSE_data, only: NSE_init

    implicit none

    ! Set species data

    spec_names(iHe4 ) = 'helium-4'
    spec_names(iC12 ) = 'carbon-12'
    spec_names(iO16 ) = 'oxygen-16'
    spec_names(iNe22) = 'neon-22'
    spec_names(iMg24) = 'magnesium-24'
    spec_names(iSi28) = 'silicon-28'

    short_spec_names(iHe4 ) = 'he4'
    short_spec_names(iC12 ) = 'c12'
    short_spec_names(iO16 ) = 'o16'
    short_spec_names(iNe22) = 'ne22'
    short_spec_names(iMg24) = 'mg24'
    short_spec_names(iSi28) = 'si28'

    aion(iHe4 ) = 4.0d0
    aion(iC12 ) = 12.0d0
    aion(iO16 ) = 16.0d0
    aion(iNe22) = 22.0d0
    aion(iMg24) = 24.0d0
    aion(iSi28) = 28.0d0
    
    zion(iHe4 ) = 2.0d0
    zion(iC12 ) = 6.0d0
    zion(iO16 ) = 8.0d0
    zion(iNe22) = 10.0d0
    zion(iMg24) = 12.0d0
    zion(iSi28) = 14.0d0

    bion(iHe4 ) = 28.296d0
    bion(iC12 ) = 92.162d0
    bion(iO16 ) = 127.619d0
    bion(iNe22) = 177.77d0
    bion(iMg24) = 198.257d0
    bion(iSi28) = 236.537d0

    yiion(:) = ONE / aion(:)
    yeion(:) = zion(:) / zion(:)
    qion(:)  = bion(:) / aion(:)

    short_aux_names(UFLAM) = 'UFLAM'
    short_aux_names(UFLDT) = 'UFLDT'
    short_aux_names(UFLSP) = 'UFLSP'
    short_aux_names(UCI  ) = 'UCI  '
    short_aux_names(UNEI ) = 'UNEI '
    short_aux_names(UPHFA) = 'UPHFA'
    short_aux_names(UPHAQ) = 'UPHAQ'
    short_aux_names(UPHQN) = 'UPhQN'
    short_aux_names(UYE  ) = 'UYE  '
    short_aux_names(UDYQN) = 'UDYQN'
    short_aux_names(UDQQN) = 'UDQQN'

    ! Initialize the NSE data

    call NSE_init()

  end subroutine actual_network_init

  !   Returns the properties of the unburned fuel and first-stage ashes based on
  !   the initial abundances as passed in
  !
  ! subroutine argements:
  !
  !   xc12initial    :  initial abundance of C12 in this cell
  !   xne22initial   :  initial abundance of Ne22 in this cell
  !   ye_f, ye_a     :  electron fraction for fuel and ash (electrons/Baryon)
  !   yi_f, yi_a     :  Ion "mole fraction" for fuel and ash (ions/Baryon)
  !   qbar_f, qbar_a :  Average binding energy per baryon for fuel and ash
  !
  !  Dean Townsley 2008
  !  
  
  subroutine paraFuelAshProperties(xc12initial, xne22initial, ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a)

    implicit none
    
    double precision, intent(in)  :: xc12initial, xne22initial
    double precision, intent(out) :: ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a

    ! simple implementation for now
    !  ash is just result of C12->Mg24
    !  likely to eventually be based on some manner of table interpolation
    ye_f = xc12initial*yeion(iC12) + xne22initial*yeion(iNe22) + (1.0-xc12initial-xne22initial)*yeion(iO16)
    ye_a = xc12initial*yeion(iMg24) + xne22initial*yeion(iNe22) + (1.0-xc12initial-xne22initial)*yeion(iO16)

    yi_f = xc12initial*yiion(iC12) + xne22initial*yiion(iNe22) + (1.0-xc12initial-xne22initial)*yiion(iO16)
    yi_a = xc12initial*yiion(iMg24) + xne22initial*yiion(iNe22) + (1.0-xc12initial-xne22initial)*yiion(iO16)

    qbar_f = xc12initial*qion(iC12) + xne22initial*qion(iNe22) + (1.0-xc12initial-xne22initial)*qion(iO16)
    qbar_a = xc12initial*qion(iC12) + xne22initial*qion(iNe22) + (1.0-xc12initial-xne22initial)*qion(iO16)

  end subroutine paraFuelAshProperties

end module actual_network
