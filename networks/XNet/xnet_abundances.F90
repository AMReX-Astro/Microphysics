!***************************************************************************************************
! abundances.f90 10/18/17
! This file contains modules and subroutines associated with matter composition.
!***************************************************************************************************

Module abundances
  !-------------------------------------------------------------------------------------------------
  ! This module contains the abundances of the nuclear species, at the previous time (yo), current
  ! time (y), and trial time (yt), as well as the time derivatives (ydot), at the trial time. The
  ! size of these arrays is allocated in an external routine, where the initial values are set.
  !-------------------------------------------------------------------------------------------------
  Use xnet_types, Only: dp
  Implicit None
  Real(dp), Allocatable :: ystart(:,:) ! Abundances at start time
  Real(dp), Allocatable :: yo(:,:)     ! Abundances at previous time
  Real(dp), Allocatable :: y(:,:)      ! Abundances at current time
  Real(dp), Allocatable :: yt(:,:)     ! Abundances at trial time
  Real(dp), Allocatable :: ydot(:,:)   ! Abundance time derivatives at trial time
  !$omp threadprivate(yo,y,yt,ydot,ystart)

Contains

  Subroutine load_initial_abundances( inab_file, abund_desc, ierr, mask_in )
    !-----------------------------------------------------------------------------------------------
    ! This routine loads the initial abundances at the start time by reading the initial abundance
    ! file or by generating an NSE composition.
    !-----------------------------------------------------------------------------------------------
    Use controls, Only: lun_ab, lun_diag, idiag, t9nse, nzone, szbatch, nzbatchmx, lzactive
!   Use nse, Only: nse_solve, ynse                                                              !NSE
    Use nuc_number, Only: ny
    Use nuclear_data, Only: nname
    Use thermo_data, Only: nstart, tstart, t9start, rhostart, yestart, nh, th, yeh
    Use xnet_eos, Only: y_moment
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: inab_file(nzone)

    ! Output variables
    Character(80), Intent(out) :: abund_desc(nzbatchmx)
    Integer, Intent(out) :: ierr

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(:)

    ! Local variables
    Real(dp) :: yein, yin(ny)
    Real(dp) :: rdt, dt, dye
    Real(dp) :: ytot, abar, zbar, z2bar, zibar
    Integer :: i, izb, izone
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask => mask_in
    Else
      mask => lzactive
    EndIf

    ! Initialize
    yestart = 0.0
    ystart = 0.0
    ierr = 0

    Do izb = 1, nzbatchmx
      If ( mask(izb) ) Then
        izone = izb + szbatch - 1

        ! Interpolate electron fraction from thermo file
        If ( nstart(izb) > 1 .and. nstart(izb) <= nh(izb) ) Then
          rdt = 1.0 / ( th(nstart(izb),izb) - th(nstart(izb)-1,izb) )
          dt = tstart(izb) - th(nstart(izb)-1,izb)
          dye = yeh(nstart(izb),izb) - yeh(nstart(izb)-1,izb)
          yestart(izb) = dt*rdt*dye + yeh(nstart(izb)-1,izb)
        ElseIf ( nstart(izb) == 1 ) Then
          yestart(izb) = yeh(1,izb)
        Else
          yestart(izb) = yeh(nh(izb),izb)
        EndIf

        ! Read the initial abundance file if not in NSE or if invalid Ye from thermo file
!       If ( t9start(izb) <= t9nse .or. yestart(izb) <= 0.0 .or. yestart(izb) >= 1.0 ) Then     !NSE
          Call read_inab_file(inab_file(izone),abund_desc(izb),yein,yin,ierr)

          ! If Ye is not provided in the initial abundance file explicitly, calculate it here
          If ( yein <= 0.0 .or. yein >= 1.0 ) Call y_moment(yin,yein,ytot,abar,zbar,z2bar,zibar)
          yestart(izb) = yein

          ! Log abundance file and description
          If ( idiag >= 0 ) Then
            Write(lun_diag,"(a)") inab_file(izone)
            Write(lun_diag,"(a)") abund_desc(izb)
          EndIf
!       EndIf                                                                                   !NSE

  ! For high temperatures, use NSE to get initial abundance
!       If ( t9start(izb) > t9nse ) Then                                                        !NSE
!         Write(lun_diag,"(a,es10.4,a,es10.4,a,f5.4)") &                                        !NSE
!           & 'NSE abundances for T9=',t9start(izb),', rho=',rhostart(izb),', Ye=',yestart(izb) !NSE
!         Call nse_solve(rhostart(izb),t9start(izb),yestart(izb))                               !NSE
!         ystart(:,izb) = ynse                                                                  !NSE
!       Else                                                                                    !NSE
          ystart(:,izb) = yin
!       EndIf                                                                                   !NSE

        If ( idiag > 0 ) Write(lun_diag,"(a,i6,a,f6.3,a,es10.3,a,f5.4)") &
          & 'Start',nstart(izb),' T9=',t9start(izb),' Rho=',rhostart(izb),' Ye=',yestart(izb)

        ! Log initial abundance
        If ( idiag >= 0 ) Write(lun_diag,"(5(a6,1es10.3))") (nname(i), ystart(i,izb), i=1,ny)
      EndIf
    EndDo

    Return
  End Subroutine load_initial_abundances

  Subroutine read_inab_file( inab_file, abund_desc, yein, yin, ierr )
    !-----------------------------------------------------------------------------------------------
    ! This routine reads initial abundances from a supplied input file.
    !-----------------------------------------------------------------------------------------------
    Use, Intrinsic :: iso_fortran_env, Only: iostat_end
    Use controls, Only: lun_stderr, lun_ab, lun_diag, idiag, getNewUnit
    Use nuc_number, Only: ny
    Use nuclear_data, Only: aa, zz
    Use xnet_interface, Only: index_from_name, string_lc, xnet_terminate
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Character (*), Intent(in) :: inab_file

    ! Output variables
    Character(80), Intent(out) :: abund_desc
    Real(dp), Intent(out) :: yein
    Real(dp), Intent(out) :: yin(ny)
    Integer, Intent(out) :: ierr

    ! Local variables
    Integer, Parameter :: nread_max = 4 ! Maximum number of abundance entries to read at once
    Character(5) :: char_tmp(nread_max), namein
    Real(dp) :: real_tmp(nread_max), xnet, znet, yext, xext, aext, zext
    Integer :: i, inuc

    ! Initialize
    yein = 0.0
    yin = 0.0
    yext = 0.0

    Open(getNewUnit(lun_ab), file=trim(inab_file), action='read', status='old', iostat=ierr)
    If ( ierr /= 0 ) Then
      Write(lun_stderr,"(2a)") 'Failed to open input file: ',trim(inab_file)
      Return
    EndIf

    Read(lun_ab,*) abund_desc
    Do
      char_tmp = '     '
      real_tmp = 0.0

      ! Read nread_max entries at once, for backwards compatability with multiple entries per line
      Read(lun_ab,*,iostat=ierr) (char_tmp(i), real_tmp(i), i=1,nread_max)
      If ( ierr == iostat_end ) Then
        Exit
      ElseIf ( ierr /= 0 ) Then
        Write(lun_stderr,"(3a,i3,a)") 'Failed while trying to read input file: ',trim(inab_file),' (',ierr,')'
        Call xnet_terminate('Failed while trying to read input file: '//trim(inab_file),ierr)
      EndIf

      ! Process each entry, chceking for special cases and converting names to lower-case
      Do i = 1, nread_max
        namein = adjustl(char_tmp(i))
        If ( len_trim(namein) > 0 .and. real_tmp(i) > 0.0 ) Then
          Call string_lc(namein)

          If ( trim(namein) == 'ye' ) Then
            yein = real_tmp(i)
          Else
            Call index_from_name(namein,inuc)
            If ( inuc < 1 .or. inuc > ny ) Then
              Write(lun_diag,"(3a)") 'Input Nuc: ',namein,' not found'
              yext = yext + real_tmp(i)
            Else
              yin(inuc) = real_tmp(i)
            EndIf
          EndIf
        EndIf
      EndDo
    EndDo

    ! If the number of species isn't divisible by 4, we need to parse the last few entries here
    Do i = 1, nread_max
      namein = adjustl(char_tmp(i))
      If ( len_trim(namein) > 0 .and. real_tmp(i) > 0.0 ) Then
        Call string_lc(namein)

        If ( trim(namein) == 'ye' ) Then
          yein = real_tmp(i)
        Else
          Call index_from_name(namein,inuc)
          If ( inuc < 1 .or. inuc > ny ) Then
            Write(lun_diag,"(3a)") 'Input Nuc: ',namein,' not found'
            yext = yext + real_tmp(i)
          Else
            yin(inuc) = real_tmp(i)
          EndIf
        EndIf
      EndIf
    EndDo
    Close(lun_ab)

    ! Total mass fraction inside network
    xnet = sum(yin*aa)
    znet = sum(yin*zz)
    If ( idiag >= 0 ) Write(lun_diag,"(a,4es15.7)") 'ynet, xnet, anet, znet: ',sum(yin),xnet,xnet,znet

    ! Normalize so total mass fraction is one
    yin = yin / xnet

    ! Calculate properties of matter not in network
    If ( yext > 0.0 ) Then
      xext = 1.0 - xnet
      aext = yext / xext
      zext = ( yein - znet ) * aext / xext
      Write(lun_diag,"(a,4es15.7)") 'yext, xext, aext, zext: ',yext,xext,aext,zext
    EndIf

    Return
  End Subroutine read_inab_file

End Module abundances
