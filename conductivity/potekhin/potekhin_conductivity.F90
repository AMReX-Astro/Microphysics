module actual_conductivity_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt), parameter :: Zmin = 1.0_rt
  real(rt), parameter :: Zmax = 60.0_rt
  real(rt), parameter :: TLGmin = 3.0_rt
  real(rt), parameter :: TLGmax = 9.0_rt
  real(rt), parameter :: RLGmin = -6.0_rt
  real(rt), parameter :: RLGmax = 9.75_rt

  ! NB: These parameters must be consistent with the table "condall06.d"
  integer, parameter :: MAXT = 19
  integer, parameter :: MAXR = 64
  integer, parameter :: MAXZ=15

  real(rt) :: AT(MAXT), AR(MAXR), AZ(MAXZ), AKAP(MAXT,MAXR,MAXZ)

contains

  subroutine actual_conductivity_init()

    implicit none

    integer :: iz, ir, it
    integer :: unit

    real(rt) :: z

    ! read in the conductivity table and initialize any data
    open(newunit=unit, file='condall06.d', status='OLD')

    print*,'Reading thermal conductivity data...'

    read(1,'(A)') ! skip the first line

    do iz = 1, MAXZ
       read(unit,*) Z, (AT(it), it = 1,MAXT)
       AZ(iz) = log10(Z)
       do ir = 1, MAXR
          read(unit,*) AR(ir),(AKAP(it,ir,iz), it=1, MAXT)
       enddo
    enddo
    close(unit=unit)

  end subroutine actual_conductivity_init

  subroutine actual_conductivity(state)

    use eos_type_module, only: eos_t
    use network, only : zion, aion, nspec

    implicit none

    type(eos_t), intent(inout) :: state

  end subroutine actual_conductivity

  subroutine CONINTER(Zion, TLG, RLG, CK, DRK, DTK)
    !
    ! This subroutine interpolates the electron thermal conductivity
    ! from the data file "condall06.d"
    ! Version 23.05.99

    ! Input: Zion - ion charge, TLG - lg(T[K]), RLG - lg(rho[g/cc])
    ! Output: CK - Log_{10} thermal conductivity (kappa) [CGS units]
    ! DRK - d log kappa / d log rho
    ! DTK - d log kappa / d log T

    ! it is also possible to obtain all second derivatives

    implicit none

    real(rt), intent(in) :: Zion, TLG, RLG
    real(rt), intent(out) :: ck, drk, dtk

    real(rt) :: ckt0, ckt0z0, ckt0z1
    real(rt) :: ckt1, ckt1z0, ckt1z1
    real(rt) :: cktm, cktmz0, cktmz1
    real(rt) :: cktp, cktpz0, cktpz1
    real(rt) :: dr2kt0, dr2kt0z0, dr2kt0z1
    real(rt) :: dr2kt1, dr2kt1z0, dr2kt1z1
    real(rt) :: dr2ktm, dr2ktmz0, dr2ktmz1
    real(rt) :: dr2ktp, dr2ktpz0, dr2ktpz1
    real(rt) :: drkt0, drkt0z0, drkt0z1
    real(rt) :: drkt1, drkt1z0, drkt1z1
    real(rt) :: drktm, drktmz0, drktmz1
    real(rt) :: drktp, drktpz0, drktpz1
    real(rt) :: drt2k, drtk, dt2k

    integer :: iz, it, ir
    integer :: irm, irp, itm, itp
    real(rt) :: xr, xt, xz0, xz1, zlg

    ZLG = log10(Zion)

    ! initial guess
    iz = MAXZ/2+1
    it = MAXT/2+1
    ir = MAXR/2+1

    call HUNT(AZ, MAXZ, ZLG, iz)
    if (iz == 0 .or. iz == MAXZ) then
       call amrex_error('CONINTER: Z out of range')
    end if

    call HUNT(AT, MAXT, TLG, it)
    if (it == 0 .or. it == MAXT) then
       call amrex_error('CONINTER: T out of range')
    end if

    call HUNT(AR, MAXR, RLG, ir)
    if (ir == 0 .or. ir == MAXR) then
       call amrex_error('CONINTER: rho out of range')
    end if

    ITM = max(1, it-1)
    ITP = min(MAXT, it+2)
    IRM = max(1, ir-1)
    IRP = min(MAXR, ir+2)

    ! Cubic interpolation in RLG:
    ! Z0:
    call CINTERP3(AR(IRM),AR(IR),AR(IR+1),AR(IRP),RLG,IR,MAXR, &
                  AKAP(ITM,IRM,IZ),AKAP(ITM,IR,IZ), &
                  AKAP(ITM,IR+1,IZ),AKAP(ITM,IRP,IZ), &
                  CKTMZ0,DRKTMZ0,DR2KTMZ0,XR)

    call CINTERP3(AR(IRM),AR(IR),AR(IR+1),AR(IRP),RLG,IR,MAXR, &
                  AKAP(IT,IRM,IZ),AKAP(IT,IR,IZ), &
                  AKAP(IT,IR+1,IZ),AKAP(IT,IRP,IZ), &
                  CKT0Z0,DRKT0Z0,DR2KT0Z0,XR)

    call CINTERP3(AR(IRM),AR(IR),AR(IR+1),AR(IRP),RLG,IR,MAXR, &
                  AKAP(IT+1,IRM,IZ),AKAP(IT+1,IR,IZ), &
                  AKAP(IT+1,IR+1,IZ),AKAP(IT+1,IRP,IZ), &
                  CKT1Z0,DRKT1Z0,DR2KT1Z0,XR)

    call CINTERP3(AR(IRM),AR(IR),AR(IR+1),AR(IRP),RLG,IR,MAXR, &
                  AKAP(ITP,IRM,IZ),AKAP(ITP,IR,IZ), &
                  AKAP(ITP,IR+1,IZ),AKAP(ITP,IRP,IZ), &
                  CKTPZ0,DRKTPZ0,DR2KTPZ0,XR)

    ! Z1:
    call CINTERP3(AR(IRM),AR(IR),AR(IR+1),AR(IRP),RLG,IR,MAXR, &
                  AKAP(ITM,IRM,IZ+1),AKAP(ITM,IR,IZ+1), &
                  AKAP(ITM,IR+1,IZ+1),AKAP(ITM,IRP,IZ+1), &
                  CKTMZ1,DRKTMZ1,DR2KTMZ1,XR)

    call CINTERP3(AR(IRM),AR(IR),AR(IR+1),AR(IRP),RLG,IR,MAXR, &
                  AKAP(IT,IRM,IZ+1),AKAP(IT,IR,IZ+1), &
                  AKAP(IT,IR+1,IZ+1),AKAP(IT,IRP,IZ+1), &
                  CKT0Z1,DRKT0Z1,DR2KT0Z1,XR)

    call CINTERP3(AR(IRM),AR(IR),AR(IR+1),AR(IRP),RLG,IR,MAXR, &
                  AKAP(IT+1,IRM,IZ+1),AKAP(IT+1,IR,IZ+1), &
                  AKAP(IT+1,IR+1,IZ+1),AKAP(IT+1,IRP,IZ+1), &
                  CKT1Z1,DRKT1Z1,DR2KT1Z1,XR)

    call CINTERP3(AR(IRM),AR(IR),AR(IR+1),AR(IRP),RLG,IR,MAXR, &
                  AKAP(ITP,IRM,IZ+1),AKAP(ITP,IR,IZ+1), &
                  AKAP(ITP,IR+1,IZ+1),AKAP(ITP,IRP,IZ+1), &
                  CKTPZ1,DRKTPZ1,DR2KTPZ1,XR)

    ! Linear interpolation in ZLG:
    XZ1=(ZLG-AZ(IZ))/(AZ(IZ+1)-AZ(IZ))
    XZ0=1.-XZ1

    CKTM = XZ0*CKTMZ0 + XZ1*CKTMZ1
    DRKTM = XZ0*DRKTMZ0 + XZ1*DRKTMZ1
    DR2KTM = XZ0*DR2KTMZ0 + XZ1*DR2KTMZ1

    CKT0 = XZ0*CKT0Z0 + XZ1*CKT0Z1
    DRKT0 = XZ0*DRKT0Z0 + XZ1*DRKT0Z1
    DR2KT0 = XZ0*DR2KT0Z0 + XZ1*DR2KT0Z1

    CKT1 = XZ0*CKT1Z0 + XZ1*CKT1Z1
    DRKT1 = XZ0*DRKT1Z0 + XZ1*DRKT1Z1
    DR2KT1 = XZ0*DR2KT1Z0 + XZ1*DR2KT1Z1

    CKTP = XZ0*CKTPZ0 + XZ1*CKTPZ1
    DRKTP = XZ0*DRKTPZ0 + XZ1*DRKTPZ1
    DR2KTP = XZ0*DR2KTPZ0 + XZ1*DR2KTPZ1

    ! Cubic interpolation in TLG:
    call CINTERP3(AT(ITM),AT(IT),AT(IT+1),AT(ITP),TLG,IT,MAXT, &
                  CKTM,CKT0,CKT1,CKTP, & ! input: values of lg kappa
                  CK,DTK,DT2K,XT) ! lg kappa, d lg k / d lg T, d2 lg k / d2 lg T
    call CINTERP3(AT(ITM),AT(IT),AT(IT+1),AT(ITP),TLG,IT,MAXT, &
                  DRKTM,DRKT0,DRKT1,DRKTP, &! input: values of d lg k / d lg rho
                  DRK,DRTK,DRT2K,XT) ! d lg k / d lg rho, d2 lgk/(d lgT d lg rho)
    return

  end subroutine CONINTER

  subroutine CINTERP3(ZM, Z0, Z1, ZP, Z, N0, MXNV, VM, V0, V1, VP, VF, DF, D2, XH)
    ! Given 4 values of Z and 4 values of V, find VF corresponding to 5th Z
    !                                                      Version 23.05.99
    ! Output: VF - interpolated value of function
    !         DF - interpolated derivative
    !         D2 - interpolated second derivative
    !         XH - fraction of the path from N0 to N0+1

    real(rt), intent(in) :: zm, z0, z1, zp
    real(rt), intent(in) :: vm, v0, v1, vp
    real(rt), intent(in) :: z
    integer, intent(in) :: n0, mxnv
    real(rt), intent(out) :: vf, df, d2, xh

    real(rt) :: c2, c3, h, hm, hp
    real(rt) :: v01, v11
    real(rt) :: x

    if (N0 <= 0 .or. N0 >= MXNV) then
       call amrex_error('CINTERP: N0 out of range')
    end if

    X = Z - Z0
    H = Z1 - Z0   ! basic interval
    XH = X/H

    if (N0 .gt. 1) then
       HM = Z0 - ZM  ! left adjoint interval
       V01 = ((V1-V0)/H**2 + (V0-VM)/HM**2) / (1./H+1./HM) ! left derivative
    endif

    if (N0 .lt. MXNV-1) then
       HP = ZP - Z1 ! right adjoint interval
       V11 = ((V1-V0)/H**2 + (VP-V1)/HP**2) / (1./H+1./HP) ! right derivative
    endif

    if (N0 .gt. 1 .and. N0 .lt. MXNV-1) then   ! Cubic interpolation
       C2 = 3.*(V1-V0)-H*(V11+2.*V01)
       C3 = H*(V01+V11)-2.*(V1-V0)
       VF = V0+V01*X+C2*XH**2+C3*XH**3
       DF = V01+(2.*C2*XH+3.*C3*XH**2)/H
       D2 = (2.*C2+6.*C3*XH)/H**2
       return
    endif

    if (N0 .eq. 1) then   ! Quadratic interpolation
       C2 = V0-V1+V11*H
       VF = V1-V11*(H-X)+C2*(1.-XH)**2
       DF = V11-2.*C2*(1.-XH)/H
       D2 = 2.*C2/H**2
    else  ! N0=MXNV-1
       C2 = V1-V0-V01*H
       VF = V0+V01*X+C2*XH**2
       DF = V01+2.*C2*XH/H
       D2 = 2.*C2/H**2
    endif

  end subroutine CINTERP3

  subroutine HUNT(XX,N,X,JLO)
    ! W.H.Press, B.P.Flannery, S.A.Teukolsky, W.T.Vetterling
    ! Numerical Receipes(Cambridge Univ., 1986)
    !     Given an array XX of length N, and given a value X,
    !     returns a value JLO such that X is between XX(JLO) and XX(JLO+1).
    !     XX must be monotonic, either increasing or decreasing. 
    !     JLO=0 or JLO=N is returned to indicate that X is out of range.
    !     JLO on input is taken as the initial guess for JLO on output.


    integer, intent(in) :: n
    real(rt), intent(in) :: xx(n)
    real(rt), intent(in) :: x
    integer, intent(inout) :: jlo

    logical :: ASCND
    logical :: bisect

    integer :: inc
    integer :: jm, jhi

    ASCND = XX(N) .gt. XX(1) ! true if ascending order, false otherwise

    bisect = .false.

    if (JLO .le. 0 .or. JLO .gt. N) then ! Input guess not useful.
       JLO=0
       JHI=N+1
       bisect = .true. ! go immediately to bisection
    endif

    if (.not. bisect) then

       INC=1 ! set the hunting increment
       if (X .ge. XX(JLO) .eqv. ASCND) then ! Hunt up:

         JHI=JLO+INC
          do while (X >= XX(jhi) .eqv. ASCND)

             if (JHI .gt. N) then ! Done hunting, since off end of table
                JHI=N+1
                exit
             end if

             JLO=JHI
             INC=INC+INC
             JHI=JLO+INC
          end do

       else ! Hunt down:
          JHI=JLO

          JLO=JHI-INC
          do while (X < XX(jlo) .eqv. ASCND)

             if (JLO .lt. 1) then ! Done hunting, since off end of table
                JLO=0
                exit
             end if

             JHI=JLO
             INC=INC+INC ! so double the increment
             JLO=JHI-INC

          end do

       endif
    end if

    ! Hunt is done, so begin the final bisection phase:
    do while (JHI - JLO /= 1)
       JM = (JHI+JLO)/2.
       if (X .ge. XX(JM) .eqv. ASCND) then
          JLO=JM
       else
          JHI=JM
       endif
    end do

  end subroutine HUNT

end module actual_conductivity_module
