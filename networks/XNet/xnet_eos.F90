!***************************************************************************************************
! eos_starkiller.f90 10/18/17
! Interface to starkiller
! This file contains routines which calculate EoS quantites needed to calculate screening
! corrections for reaction rates.
!***************************************************************************************************

Module xnet_eos
  Use eos_type_module, Only: eos_t
  Implicit None
  Type(eos_t) :: eos_state
  !$omp threadprivate(eos_state)

Contains

  Subroutine eos_initialize
    !-------------------------------------------------------------------------------------------------
    ! This routine initializes starkiller
    !-------------------------------------------------------------------------------------------------
    Use actual_eos_module, Only: actual_eos_init
    Implicit None
    Return
  End Subroutine eos_initialize

  Subroutine eos_cv(rho,t9,y,cv)
    Use xnet_constants, Only: avn, epmev
    Use controls, Only: idiag, iscrn, lun_diag
    Use nuc_number, Only: ny
    Use xnet_types, Only: dp

    Use actual_eos_module, Only: actual_eos
    Use eos_type_module, Only: eos_input_rt
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9, rho, y(ny)

    ! Ouput variables
    Real(dp), Intent(out) :: cv

    ! Local variables
    Real(dp) :: ye, ytot, abar, zbar, z2bar, zibar

    ! Call the eos if it hasn't already been called for screening
    If ( iscrn <= 0 ) Then

      ! Calculate Ye and other needed moments of the abundance distribution
      Call y_moment(y,ye,ytot,abar,zbar,z2bar,zibar)

      ! Load input variables for the eos
      eos_state%rho = rho
      eos_state%T = t9*1e9
      eos_state%abar = abar
      eos_State%zbar = ye*abar

      Call actual_eos(eos_input_rt,eos_state)
    EndIf

    ! Convert units from ergs/g to MeV/nucleon and K to GK
    cv = eos_state%cv * 1.0d9/epmev/avn

    If ( idiag > 0 ) Write(lun_diag,"(a,5es23.15)") 'CV',t9,rho,eos_state%zbar/eos_state%abar,cv

    Return
  End Subroutine eos_cv

  Subroutine xnet_eos_interface(t9,rho,y,ye,ztilde,zinter,lambda0,gammae,dztildedt9)
    !-------------------------------------------------------------------------------------------------
    ! This routine calls the Helmholtz EOS with the input temperature, density and composition. It
    ! returns the factors needed for screening.
    !-------------------------------------------------------------------------------------------------
    Use xnet_constants, Only: avn, bok, clt, e2, ele_en, emass, hbar, pi, pi2, third, two3rd
    Use controls, Only: idiag, iheat, lun_diag
    Use nuc_number, Only: ny
    Use screening_data, Only: thbim2, twm2bi
    Use xnet_types, Only: dp

    Use actual_eos_module, Only: actual_eos
    Use eos_type_module, Only: eos_input_rt
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: t9, rho, y(ny)

    ! Output variables
    Real(dp), Intent(out) :: ztilde, zinter, lambda0, gammae, ye, dztildedt9

    ! Local variables
    Real(dp) :: ytot, bkt, abar, zbar, z2bar, zibar
    Real(dp) :: etae, sratio, efermkt, rel_ef, efc, ae, dsratiodeta

    ! Calculate Ye and other needed moments of the abundance distribution
    Call y_moment(y,ye,ytot,abar,zbar,z2bar,zibar)

    ! Load input variables for the eos
    eos_state%rho = rho
    eos_state%T = t9*1e9
    eos_state%abar = abar
    eos_State%zbar = ye*abar

    ! Call the eos
    Call actual_eos(eos_input_rt,eos_state)
    etae = eos_state%eta

    ! Calculate electon distribution
    bkt = bok*t9
    rel_ef = hbar * (3.0*pi2*rho*avn*ye)**third / (emass*clt)
    efermkt = ele_en * (sqrt(1.0 + rel_ef**2) - 1.0) / bkt
    efc = 0.5*hbar**2 * (3.0*pi2*rho*avn*ye)**two3rd / (emass*bkt)
    ! Write(lun_diag,"(a4,6es12.5)") 'MUh',bkt,etae,efermkt,efc,rel_ef

    ! Calculate ratio f'/f for electrons (Salpeter, Eq. 24)
    Call salpeter_ratio(etae,sratio,dsratiodeta)
    ztilde = sqrt(z2bar + zbar*sratio)
    If ( iheat > 0 ) Then
      dztildedt9 = 0.5*zbar/ztilde * dsratiodeta*eos_state%detadt*1e9
    Else
      dztildedt9 = 0.0
    EndIf

    ! Calculate plasma quantities
    lambda0 = sqrt(4.0*pi*rho*avn*ytot) * (e2/bkt)**1.5 ! DGC, Eq. 3
    ae = (3.0 / (4.0*pi*avn*rho*ye))**third ! electron-sphere radius
    gammae = e2 / (ae*bkt) ! electron Coulomb coupling parameter
    zinter = zibar / (ztilde**thbim2 * zbar**twm2bi)
    If ( idiag > 0 ) Write(lun_diag,"(a14,9es23.15)") 'Helmholtz EOS', &
      & t9,rho,ye,z2bar,zbar,sratio,ztilde,ztilde*lambda0,gammae

    Return
  End Subroutine xnet_eos_interface

  Subroutine salpeter_ratio(eta,ratio,dratiodeta)
    !-------------------------------------------------------------------------------------------------
    ! This routine calculates the Salpeter (1954) ratio f'/f(eta) needed for electron screening.
    ! eta is the ratio of electron chemical potential to kT.
    !
    ! Calculation uses Fermi function relation d/dx f_(k+1) = (k+1) f_k and the rational function
    ! expansions of Fukushima (2015; AMC 259 708) for the F-D integrals of order 1/2, -1/2, and -3/2.
    !-------------------------------------------------------------------------------------------------
    Use controls, Only: iheat, lun_diag
    Use fd, Only: fdm1h, fd1h, fdm3h
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: eta

    ! Output variables
    Real(dp), Intent(out) :: ratio, dratiodeta

    ! Local variables
    Integer :: i
    Real(dp) :: fermip, fermim
    Real(dp) :: dfmdeta, dfpdeta

    ! Calculate f_(-1/2) and f_(1/2)
    fermim = fdm1h(eta)
    fermip = fd1h(eta)

    ! Evalutate the Salpeter ratio
    ratio = 0.5 * fermim/fermip
    If ( iheat > 0 ) Then
      dfmdeta = -0.5 * fdm3h(eta)
      dfpdeta = +0.5 * fermim
      dratiodeta = ratio * (dfmdeta/fermim - dfpdeta/fermip)
    Else
      dratiodeta = 0.0
    EndIf
    ! Write(lun_diag,"(1x,4es12.4)") eta,ratio,fermim,fermip

    Return
  End Subroutine salpeter_ratio

  Subroutine y_moment(y,ye,ytot,abar,zbar,z2bar,zibar)
    !-------------------------------------------------------------------------------------------------
    ! This routine calculates moments of the abundance distribution for the EOS.
    !-------------------------------------------------------------------------------------------------
    Use controls, Only: idiag, lun_diag
    Use nuc_number, Only: ny
    Use nuclear_data, Only: aa, zz, zz2, zzi
    Use xnet_types, Only: dp
    Implicit None

    ! Input variables
    Real(dp), Intent(in) :: y(ny)

    ! Output variables
    Real(dp), Intent(out) :: ye, ytot, abar, zbar, z2bar, zibar

    ! Local variables
    Real(dp) :: atot, ztot

    ! Calculate abundance moments
    ytot  = sum(y)
    atot  = sum(y * aa)
    ztot  = sum(y * zz)
    abar  = atot / ytot
    zbar  = ztot / ytot
    z2bar = sum(y * zz2) / ytot
    zibar = sum(y * zzi) / ytot
    ye = ztot
    If ( idiag > 0 ) Write(lun_diag,"(a4,6es23.15)") 'YMom',ytot,abar,zbar,z2bar,zibar,ye

    Return
  End Subroutine y_moment
End Module xnet_eos
