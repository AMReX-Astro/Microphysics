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

  logical          :: useBurn
  logical          :: useShockBurn
  double precision :: enucDtFactor

  ! Logical file unit for saving detonation ignition points
  integer, save     :: burn_lun = 31
  character(len=80) :: detIgnFileName 
  
  logical,          save :: thermalReact, autoDDT
  double precision, save :: thermalReactInFlameThreshold

  ! table for detonation ignition points
  double precision, save :: IgnRho, IgnRhoFact, IgnPhfa, IgnDist, IgnRad, IgnSep
  integer,          save :: IgnNum, IgnNumMax, IgnRhoCalib
  
  ! variable for saving processor-local and global neutrino loss energy integrals
  double precision, save :: neutLossThisProcStep, neutLoss

  double precision, allocatable, dimension(:), save :: IgnTime, IgnX, IgnY, IgnZ, IgnR
  
contains

  subroutine actual_network_init()

    use extern_probin_module, only: bn_thermalReact, bn_thermalReactInFlameThreshold, bn_neutLoss, &
                                    bn_detIgnFile, bn_detIgnFileName, bn_autoDDT, &
                                    pbIgnRho, pbIgnRhoCalib, pbIgnRhoFact, &
                                    pbIgnPhfa, pbIgnDist, pbIgnRad, &
                                    pbIgnSep, pbIgnNumMax
    use bl_constants_module, only: ONE
    use NSE_data, only: NSE_init

    implicit none

    character(len=4)        :: isotopeName
    real                    :: abar,zbar,bindEnergy,spinfactor
    integer                 :: nfound, nread

    logical                 :: detIgnFile
    integer                 :: i, istat

!    if (.not. restart) then
       neutLoss = 0.0
!    endif
    neutLossThisProcStep = 0.0

    !-------------
    !  determine whether we are detonating manually or automatically
    detIgnFile     = bn_detIgnFile
    detIgnFileName = bn_detIgnFileName
    autoDDT        = bn_autoDDT

    if (detIgnFile .and. autoDDT) then
       call bl_error("Cannot detonate manually and automatically!")
    endif

    !  read detonation ignition points from file
    if (detIgnFile) then

       open(unit=infile_unit,file=detIgnFileName,status='OLD',iostat=istat)
       if (istat /= 0) call bl_error("Unable to open detonation ignition points file")

       ! eat header
       read(infile_unit,*)
       read(infile_unit,*) IgnNum
       allocate(IgnTime(IgnNum))
       allocate(IgnX(IgnNum))
       allocate(IgnY(IgnNum))
       allocate(IgnZ(IgnNum))
       allocate(IgnR(IgnNum))
       do i = 1, IgnNum
          read(infile_unit,*) IgnTime(i), IgnX(i), IgnY(i), IgnZ(i), IgnR(i)
       enddo
!       if (bn_meshMe == MASTER_PE) then
!          write(*,*)  'read detonation ignition points t x y z'
!          do i = 1, IgnNum
!             write(*,*) IgnTime(i), IgnX(i), IgnY(i), IgnZ(i), IgnR(i)
!          enddo
!       endif
       close(unit=infile_unit)

    else if (autoDDT) then

       IgnRho      = pbIgnRho
       IgnRhoCalib = pbIgnRhoCalib
       IgnRhoFact  = pbIgnRhoFact
       IgnPhfa     = pbIgnPhfa
       IgnDist     = pbIgnDist
       IgnRad      = pbIgnRad
       IgnSep      = pbIgnSep
       IgnNumMax   = pbIgnNumMax

       IgnRho = exp( log(IgnRho) - IgnRhoCalib )
       if (IgnRhoFact >= 1.0)  &
            call bl_error("pbIgnRhoFact must be less than 1")
       IgnRhoFact = 10.e0**( IgnRhoFact * log10(IgnRho) ) 

       ! allocate global arrays for detonation points
       IgnNum = 0
       allocate(IgnTime(IgnNumMax),stat=istat)
       if (istat/=0) call bl_error("Cannot allocate IgnTime in paraburn init")
       allocate(IgnX(IgnNumMax),stat=istat)
       if (istat/=0) call bl_error("Cannot allocate IgnX in paraburn init")
       allocate(IgnY(IgnNumMax),stat=istat)
       if (istat/=0) call bl_error("Cannot allocate IgnY in paraburn init")
       allocate(IgnZ(IgnNumMax),stat=istat)
       if (istat/=0) call bl_error("Cannot allocate IgnZ in paraburn init")

       if (restart) then

          open(unit=infile_unit,file=bn_detIgnFileName,status='unknown',IOSTAT=istat)
          if (istat /= 0) call bl_error("Unable to open detonation ignition points file")

          istat = 0
          do i = 1, IgnNumMax

             read(infile_unit,*,IOSTAT=istat) IgnTime(i), IgnX(i),  &
                  IgnY(i), IgnZ(i)
             if (istat/=0) exit

          enddo

          if (istat > 0) then
             call bl_error("Unable to read detonation ignition points file")
          else if (istat < 0) then !EOF reached
             IgnNum = i - 1
          else
             call bl_error("IgnNumMax is too small to read in previous detonation ignition points")
          endif

          close(unit=infile_unit)
       endif

    else

       IgnNum = 0

    endif

    ! use built-in physical constants

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
