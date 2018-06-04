!***************************************************************************************************
! conditions.f90 10/18/17
! This file contains modules and subroutines associated with the thermodynamic conditions in the
! matter undergoing nucleosynthesis.
!***************************************************************************************************

Module conditions
  !-------------------------------------------------------------------------------------------------
  ! This module contains data on the current time and thermodynamic conditions.
  !-------------------------------------------------------------------------------------------------
  Use xnet_types, Only: dp
  Implicit None
  Real(dp), Allocatable :: t(:)                  ! Time at the beginning of the current timestep
  Real(dp), Allocatable :: tt(:)                 ! Trial time for end of current timestep
  Real(dp), Allocatable :: to(:)                 ! Time at the beginning of the previous timestep
  Real(dp), Allocatable :: tdel(:)               ! Trial duration of timestep
  Real(dp), Allocatable :: tdel_next(:)          ! Integrator estimate of duration for next timestep
  Real(dp), Allocatable :: tdel_old(:)           ! Duration of the previous timestep
  Real(dp), Allocatable :: t9t(:),rhot(:),yet(:) ! Temperature, density, and electron fraction at trial time
  Real(dp), Allocatable :: t9(:),rho(:),ye(:)    ! Temperature, density, and electron fraction at current time
  Real(dp), Allocatable :: t9o(:),rhoo(:),yeo(:) ! Temperature, density, and electron fraction at previous time
  Real(dp), Allocatable :: t9dot(:)              ! Time derivative of temperature at trial time
  Real(dp), Allocatable :: cv(:)                 ! Specific heat at constant volume
  !$omp threadprivate(t,tt,to,tdel,tdel_next,tdel_old,t9t,t9,t9o,t9dot,rhot,rho,rhoo,yet,ye,yeo,cv)

  Integer, Allocatable :: nt(:)    ! Point in thermo trajectory for current time
  Integer, Allocatable :: ntt(:)   ! Point in thermo trajectory for trial time
  Integer, Allocatable :: nto(:)   ! Point in thermo trajectory for previous time
  Integer, Allocatable :: ints(:)  ! Nucleus governing timestep
  Integer, Allocatable :: intso(:) ! Nucleus governing last timestep
  !$omp threadprivate(nto,nt,ntt,ints,intso)

End Module conditions

Module thermo_data
  !-------------------------------------------------------------------------------------------------
  ! This module contains the thermodynamic trajectory which the network follows.
  !-------------------------------------------------------------------------------------------------
  Use xnet_types, Only: dp
  Implicit None
  Integer, Parameter    :: nhmx = 50000 ! The max number of thermo points
  Integer, Allocatable  :: nh(:)        ! The actual number of points
  Integer, Allocatable  :: nstart(:)
  Real(dp), Allocatable :: tstart(:),tstop(:),tdelstart(:),t9start(:),rhostart(:),yestart(:)
  Real(dp), Allocatable :: th(:,:),t9h(:,:),rhoh(:,:),yeh(:,:)
  !$omp threadprivate(nh,nstart,tstart,tstop,tdelstart,t9start,rhostart,yestart,th,t9h,rhoh,yeh)

Contains

  Subroutine read_thermo_file( thermo_file, thermo_desc, ierr, mask_in )
    !-----------------------------------------------------------------------------------------------
    ! Read the thermdynamic trajectory
    !-----------------------------------------------------------------------------------------------
    Use, Intrinsic :: iso_fortran_env, Only: iostat_end
    Use controls, Only: idiag, lun_diag, lun_th, nzone, szbatch, nzbatchmx, lzactive, getNewUnit
!   Use neutrino_data, Only: fluxcms, tmevnu                                                    !NNU
    Use xnet_interface, Only: xnet_terminate
    Implicit None

    ! Input variables
    Character(*), Intent(in) :: thermo_file(nzone)

    ! Output variables
    Character(80), Intent(out) :: thermo_desc(nzbatchmx)
    Integer, Intent(out) :: ierr

    ! Optional variables
    Logical, Optional, Target, Intent(in) :: mask_in(:)

    ! Local variables
    Integer :: n, izb, izone
    Logical, Pointer :: mask(:)

    If ( present(mask_in) ) Then
      mask => mask_in
    Else
      mask => lzactive
    EndIf

    ! Initialize
    tstart = 0.0
    tstop = 0.0
    tdelstart = 0.0
    th = 0.0
    t9h = 0.0
    rhoh = 0.0
    yeh = 0.0
!   fluxcms = 0.0                                                                               !NNU
!   tmevnu = 0.0                                                                                !NNU
    ierr = 0

    Do izb = 1, nzbatchmx
      If ( mask(izb) ) Then
        izone = izb + szbatch - 1

        Open(getNewUnit(lun_th), file=trim(thermo_file(izb)), action='read', status='old', iostat=ierr)
        If ( ierr /= 0 ) Then
          Call xnet_terminate('Failed to open input file: '//trim(thermo_file(izb)))
        EndIf

        Read(lun_th,*) thermo_desc(izb)
        Read(lun_th,*) tstart(izb)
        Read(lun_th,*) tstop(izb)
        Read(lun_th,*) tdelstart(izb)
        Do n = 1, nhmx
          Read(lun_th,*,iostat=ierr) th(n,izb),t9h(n,izb),rhoh(n,izb)                           !NOTNSE !NOTNNU
!         Read(lun_th,*,iostat=ierr) th(n,izb),t9h(n,izb),rhoh(n,izb),yeh(n,izb)                !NSE !NOTNNU
!         Read(lun_th,*,iostat=ierr) th(n,izb),t9h(n,izb),rhoh(n,izb),yeh(n,izb),fluxcms(n,:,izb),tmevnu(n,:,izb)       !NNU

          If ( ierr == iostat_end ) Then
            If ( idiag > 2 ) Write(lun_diag,"(a,i6,a)") 'End of Thermo File Reached after ',n,' records'
            Exit
          ElseIf ( ierr /= 0 ) Then
            Call xnet_terminate('Failed while trying to read input file: '//trim(thermo_file(izb)),ierr)
          EndIf
        EndDo
        nh(izb) = n - 1
        Close(lun_th)

        ! Do not use tdelstart from thermo files
        tdelstart(izb) = min(0.0,tdelstart(izb))

        ! Log thermo description
        If ( idiag >= 0 ) Write(lun_diag,"(a)") thermo_desc(izb)
      EndIf
    EndDo

    ! Convert to appropriate units (CGS, except temperature (GK) and neutrino flux)
!   t9h = t9h * 1.0e-9
!   fluxcms = 1.0e-42 * fluxcms                                                                 !NNU

    Return
  End Subroutine read_thermo_file

End Module thermo_data

Subroutine t9rhofind(kstep,tf,nf,t9f,rhof,mask_in)
  !-------------------------------------------------------------------------------------------------
  ! This routine calculates t9 and rho as a function of time, either via interpolation or from an
  ! analytic expression.
  !-------------------------------------------------------------------------------------------------
  Use thermo_data, Only: nh, th, t9h, rhoh
  Use controls, Only: lun_diag, lun_stdout, nzbatchmx, lzactive
  Use xnet_types, Only: dp
  Implicit None

  ! Input variables
  Integer, Intent(in) :: kstep
  Real(dp), Intent(in) :: tf(:)

  ! Input/Output variables
  Integer, Intent(inout) :: nf(size(tf))
  Real(dp), Intent(inout) :: t9f(size(tf)), rhof(size(tf))

  ! Optional variables
  Logical, Optional, Target, Intent(in) :: mask_in(:)

  ! Local variables
  Real(dp) :: dt, rdt, dt9, drho
  Integer :: n, izb
  Logical, Pointer :: mask(:)

  If ( present(mask_in) ) Then
    mask => mask_in
  Else
    mask => lzactive
  End If

  Do izb = 1, nzbatchmx
    If ( mask(izb) ) Then

      ! For constant conditions (nh = 1), set temperature and density
      If ( nh(izb) == 1 ) Then
        t9f(izb) = t9h(1,izb)
        rhof(izb) = rhoh(1,izb)
        nf(izb) = 1

      ! Otherwise, calculate T9 and rho by interpolation
      Else
        Do n = 1, nh(izb)
          If ( tf(izb) <= th(n,izb) ) Exit
        EndDo
        nf(izb) = n
        If ( n > 1 .and. n <= nh(izb) ) Then
          rdt = 1.0 / (th(n,izb)-th(n-1,izb))
          dt = tf(izb) - th(n-1,izb)
          dt9 = t9h(n,izb) - t9h(n-1,izb)
          drho = rhoh(n,izb) - rhoh(n-1,izb)
          t9f(izb) = dt*rdt*dt9 + t9h(n-1,izb)
          rhof(izb) = dt*rdt*drho + rhoh(n-1,izb)
        ElseIf ( n == 1 ) Then
          t9f(izb) = t9h(1,izb)
          rhof(izb) = rhoh(1,izb)
        Else
          t9f(izb) = t9h(nh(izb),izb)
          rhof(izb) = rhoh(nh(izb),izb)
          Write(lun_stdout,*) 'Time beyond thermodynamic range',tf(izb),' >',th(nh(izb),izb)
        EndIf
      EndIf

      ! Output T9 and rho
!     Write(lun_diag,"(a5,i5,3es12.4)") 'T9rho',kstep,tf(izb),t9f(izb),rhof(izb)
    EndIf
  EndDo

  Return
End Subroutine t9rhofind
