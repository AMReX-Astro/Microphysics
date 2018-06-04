!***************************************************************************************************
! screening.f90 10/18/17
! This file contains the routines needed to calculate screening corrections for reaction rates.
!***************************************************************************************************

Subroutine screening(mask_in)
  !-------------------------------------------------------------------------------------------------
  ! This routine calculates the screening factors necessary for XNet. An equation of state,
  ! typically Helmholtz (Timmes & Swesty 1999) is used to determine the electron distribution and
  ! chemical potential.
  !
  ! References:
  ! Weak Screening:
  !    Salpeter (1954) Aust J Phys 7 353.
  ! Intermediate Screening:
  !    DeWitt, Graboske & Cooper (1973) ApJ 181 439.
  !    Graboske, DeWitt, Grossman & Cooper (1973) ApJ 181 457.
  ! Strong Screening:
  !    DeWitt & Slattery (2003) Contrib Plasma Phys 43 279.
  !-------------------------------------------------------------------------------------------------
  Use abundances, Only: yt
  Use conditions, Only: rhot, t9t, yet
  Use controls, Only: idiag, iheat, iscrn, lun_diag, nzbatchmx, szbatch, lzactive
  Use cross_sect_data, Only: n1i, n2i, n3i, nreac
  Use nuclear_data, Only: izmax, nname
  Use screening_data
  Use xnet_constants, Only: third
  Use xnet_eos, Only: xnet_eos_interface
  Use xnet_types, Only: dp
  Implicit None

  ! Optional variables
  Logical, Optional, Target, Intent(in) :: mask_in(:)

  ! Local variables
  Integer :: j, mu, izb, izone
  Real(dp), Dimension(nreac(2)) :: h2w, h2i, h2s, lambda12
  Real(dp), Dimension(nreac(2)) :: dh2wdt9, dh2idt9, dh2sdt9
  Real(dp), Dimension(nreac(3)) :: h3w, h3i, h3s, lambda123
  Real(dp), Dimension(nreac(3)) :: dh3wdt9, dh3idt9, dh3sdt9
  Real(dp), Dimension(0:izmax+2) :: gammaz, fhs, fhi, dfhsdt9
  Real(dp) :: ztilde, zinter, lambda0, gammae, dztildedt9
  Logical, Pointer :: mask(:)

  If ( present(mask_in) ) Then
    mask => mask_in
  Else
    mask => lzactive
  EndIf

  Do izb = 1, nzbatchmx
    If ( mask(izb) ) Then

      ! Call EOS to get plasma quantities
      call xnet_eos_interface(t9t(izb),rhot(izb),yt(:,izb),yet(izb), &
        & ztilde,zinter,lambda0,gammae,dztildedt9)

      ! Loop over 1 reactanct reactions to build screening factors.
      h1(:,izb) = 0.0
      If ( iheat > 0 ) dh1dt9(:,izb) = 0.0

      If ( iscrn > 0 ) Then

        ! Calculate screening energies as a function of Z, for prescriptions that follow this approach
        gammaz(0) = 0.0
        gammaz(1:izmax+2) = gammae * zseq53(1:izmax+2)
        fhi(0) = 0.0
        fhi(1:izmax+2) = kbi * zinter * lambda0**bi * zseqi(1:izmax+2)
        fhs(0) = 0.0
        fhs(1:izmax+2) = + cds(1) * gammaz(1:izmax+2) &
          &              + cds(2) * gammaz(1:izmax+2)**cds(5) / cds(5) &
          &              + cds(3) * log(gammaz(1:izmax+2)) &
          &              + cds(4)
        dfhsdt9(0) = 0.0
        dfhsdt9(1:izmax+2) = + cds(1) * gammaz(1:izmax+2) &
          &                  + cds(2) * gammaz(1:izmax+2)**cds(5) &
          &                  + cds(3) * log(gammaz(1:izmax+2))
        dfhsdt9(1:izmax+2) = -dfhsdt9(1:izmax+2) / t9t(izb)

        ! Weak and intermediate screening factors, Table 4 of Graboske et al. (1973)
        lambda12 = zeta2w * ztilde * lambda0
        h2w = lambda12
!       h2i = kbi * zinter * lambda0**bi * zeta2i
        h2i = fhi(iz2c) - fhi(iz21) - fhi(iz22)

        ! Strong screening from Dewitt & Slattery (2003) using linear mixing.
        h2s = fhs(iz21) + fhs(iz22) - fhs(iz2c)

        ! Select Screening factor for 2 reactant reactions
        Where ( iz21 == 0 .or. iz22 == 0 )
          h2(:,izb) = 0.0
        ElseWhere ( lambda12 < 0.1 )
          h2(:,izb) = h2w
        ElseWhere ( lambda12 < 2.0 )
          h2(:,izb) = h2i
        ElseWhere ( lambda12 > 5.0 )
          h2(:,izb) = h2s
        ElseWhere ( h2i < h2s )
          h2(:,izb) = h2i
        ElseWhere
          h2(:,izb) = h2s
        EndWhere

        ! Weak and intermediate screening factors, Table 4 of Graboske+ (1973)
        lambda123 = zeta3w * ztilde * lambda0
        h3w = lambda123
!       h3i = kbi * zinter * lambda0**bi * zeta3i
        h3i = fhi(iz3c) - fhi(iz31) - fhi(iz32) - fhi(iz33)

        ! Strong screening from Dewitt & Slattery (2003) using linear mixing.
        h3s = fhs(iz31) + fhs(iz32) + fhs(iz33) - fhs(iz3c)

        ! Select screening factor for 3 reactant reactions
        Where ( iz31 == 0 .or. iz32 == 0 .or. iz33 == 0 )
          h3(:,izb) = 0.0
        ElseWhere ( lambda123 < 0.1 )
          h3(:,izb) = h3w
        ElseWhere ( lambda123 < 2.0 )
          h3(:,izb) = h3i
        ElseWhere ( lambda123 > 5.0 )
          h3(:,izb) = h3s
        ElseWhere ( h3i < h3s )
          h3(:,izb) = h3i
        ElseWhere
          h3(:,izb) = h3s
        EndWhere

        If ( iheat > 0 ) Then
          dh2wdt9 = +h2w * (dztildedt9/ztilde - 1.5/t9t(izb))
          dh2idt9 = -h2i * (thbim2*dztildedt9/ztilde + bi*1.5/t9t(izb))
          dh2sdt9 = dfhsdt9(iz21) + dfhsdt9(iz22) - dfhsdt9(iz2c)
          Where ( iz21 == 0 .or. iz22 == 0 )
            dh2dt9(:,izb) = 0.0
          ElseWhere ( lambda12 < 0.1 )
            dh2dt9(:,izb) = dh2wdt9
          ElseWhere ( lambda12 < 2.0 )
            dh2dt9(:,izb) = dh2idt9
          ElseWhere ( lambda12 > 5.0 )
            dh2dt9(:,izb) = dh2sdt9
          ElseWhere ( h2i < h2s )
            dh2dt9(:,izb) = dh2idt9
          ElseWhere
            dh2dt9(:,izb) = dh2sdt9
          EndWhere

          dh3wdt9 = +h3w * (dztildedt9/ztilde - 1.5/t9t(izb))
          dh3idt9 = -h3i * (thbim2*dztildedt9/ztilde + bi*1.5/t9t(izb))
          dh3sdt9 = dfhsdt9(iz31) + dfhsdt9(iz32) + dfhsdt9(iz33) - dfhsdt9(iz3c)
          Where ( iz31 == 0 .or. iz32 == 0 .or. iz33 == 0 )
            dh3dt9(:,izb) = 0.0
          ElseWhere ( lambda123 < 0.1 )
            dh3dt9(:,izb) = dh3wdt9
          ElseWhere ( lambda123 < 2.0 )
            dh3dt9(:,izb) = dh3idt9
          ElseWhere ( lambda123 > 5.0 )
            dh3dt9(:,izb) = dh3sdt9
          ElseWhere ( h3i < h3s )
            dh3dt9(:,izb) = dh3idt9
          ElseWhere
            dh3dt9(:,izb) = dh3sdt9
          EndWhere
        EndIf

        If ( idiag > 3 ) Then
          izone = izb + szbatch - 1
          Write(lun_diag,"(a,i5)") 'SCREEN',izone
          Write(lun_diag,"(3a5,i6,5es23.15)") &
            & ('H2',(nname(n2i(j,mu)),j=1,2),mu,lambda12(mu), h2(mu,izb),h2w(mu),h2i(mu),h2s(mu),mu=1,nreac(2))
          Write(lun_diag,"(4a5,i6,5es23.15)") &
            & ('H3',(nname(n3i(j,mu)),j=1,3),mu,lambda123(mu),h3(mu,izb),h3w(mu),h3i(mu),h3s(mu),mu=1,nreac(3))
          If ( iheat > 0 ) Then
            Write(lun_diag,"(a7,2a5,i6,4es23.15)") &
              & ('dH2/dT9',(nname(n2i(j,mu)),j=1,2),mu,dh2dt9(mu,izb),dh2wdt9(mu),dh2idt9(mu),dh2sdt9(mu),mu=1,nreac(2))
            Write(lun_diag,"(a7,3a5,i6,4es23.15)") &
              & ('dH3/dT9',(nname(n3i(j,mu)),j=1,3),mu,dh3dt9(mu,izb),dh3wdt9(mu),dh3idt9(mu),dh3sdt9(mu),mu=1,nreac(3))
          EndIf
        EndIf
      Else
        h2(:,izb) = 0.0
        h3(:,izb) = 0.0
        If ( iheat > 0 ) Then
          dh2dt9(:,izb) = 0.0
          dh3dt9(:,izb) = 0.0
        EndIf
      EndIf
    Else
      h1(:,izb) = 0.0
      h2(:,izb) = 0.0
      h3(:,izb) = 0.0
      If ( iheat > 0 ) Then
        dh1dt9(:,izb) = 0.0
        dh2dt9(:,izb) = 0.0
        dh3dt9(:,izb) = 0.0
      EndIf
    EndIf
  EndDo

  Return
End Subroutine screening
