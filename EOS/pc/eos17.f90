!! Equation of state for fully ionized electron-ion plasmas (EOS EIP)
! A.Y.Potekhin & G.Chabrier, Contrib. Plasma Phys., 50 (2010) 82, 
!       and references therein
! Please communicate comments/suggestions to Alexander Potekhin:
!                                            palex@astro.ioffe.ru
! Previously distributed versions (obsolete):
!   eos2000, eos2002, eos2004, eos2006, eos2007, eos2009, eos10, eos11,
!     eos13, and eos14.
! Last update: 04.03.21. All updates since 2008 are listed below.
!!   L I S T   O F   S U B R O U T I N E S :
!  MAIN (normally commented-out) - example driving routine.
!  MELANGE9 - for arbitrary ionic mixture, renders total (ion+electron)
!          pressure, internal energy, entropy, heat capacity (all
!          normalized to the ionic ideal-gas values), logarithmic
!          derivatives of pressure over temperature and density.
!  EOSFI8 - nonideal (ion-ion + ion-electron + electron-electron)
!          contributions to the free and internal energies, pressure,
!          entropy, heat capacity, derivatives of pressure over
!          logarithm of temperature and over logarithm of density (all
!          normalized to the ionic ideal-gas values) for one ionic
!          component in a mixture.
!  FITION9 - ion-ion interaction contributions to the free and internal
!          energies, pressure, entropy, heat capacity, derivatives of
!          pressure over logarithms of temperature and density.
!  FSCRliq8 - ion-electron (screening) contributions to the free and
!          internal energies, pressure, entropy, heat capacity,
!          derivatives of pressure over logarithms of temperature and
!          density in the liquid phase for one ionic component in a
!          mixture.
!  FSCRsol8 - ion-electron (screening) contributions to the free and
!          internal energies, pressure, entropy, heat capacity,
!          derivatives of pressure over logarithms of temperature and
!          density for monoionic solid.
!  FHARM12 - harmonic (including static-lattice and zero-point)
!          contributions to the free and internal energies, pressure,
!          entropy, heat capacity, derivatives of pressure over
!          logarithms of temperature and density for solid OCP.
!  HLfit12 - the same as FHARM12, but only for thermal contributions
!  ANHARM8 - anharmonic contributions to the free and internal energies,
!          pressure, entropy, heat capacity, derivatives of pressure
!          over logarithms of temperature and density for solid OCP.
!  CORMIX - correction to the linear mixing rule for the Coulomb
!          contributions to the thermodynamic functions in the liquid.
!  ELECT11 - for an ideal electron gas of arbitrary degeneracy and
!          relativity at given temperature and electron chemical
!          potential, renders number density (in atomic units), free
!          energy, pressure, internal energy, entropy, heat capacity 
!          (normalized to the electron ideal-gas values), logarithmic
!          derivatives of pressure over temperature and density.
!  EXCOR7 - electron-electron (exchange-correlation) contributions to
!          the free and internal energies, pressure, entropy, heat
!          capacity, derivatives of pressure over logarithm of
!          temperature and over logarithm of density (all normalized
!          to the classical electron ideal-gas values).
!  FERINV7 - inverse non-relativistic Fermi integrals of orders -1/2,
!          1/2, 3/2, 5/2, and their first and second derivatives.
!  BLIN9 - relativistic Fermi-Dirac integrals of orders 1/2, 3/2, 5/2,
!          and their first, second, and some third derivatives.
!  CHEMFIT7 - electron chemical potential at given density and
!          temperature, and its first derivatives over density and
!          temperature and the second derivative over temperature.
!!   I M P R O V E M E N T S   S I N C E   2 0 0 8 :
!  FHARM8 uses a fit HLfit8 to the thermal free energy of the harmonic
!   Coulomb lattice, which is more accurate than its predecessor FHARM7.
!   Resulting corrections amount up to 20% for the ion heat capacity.
!   Accordingly, S/R D3fit and FthCHA7 deleted (not used anymore).
!  BLIN7 upgraded to BLIN8:
!      - cleaned (a never-reached if-else branch deleted);
!      - Sommerfeld (high-\chi) expansion improved;
!      - some third derivatives added.
!  CORMIX added (and MELANGE7 upgraded to MELANGE8 accordingly).
!  ANHARM7 upgraded to ANHARM8, more consistent with Wigner-Kirkwood.
!  Since the T- and rho-dependences of individual Z values in a mixture
!    are not considered, the corresponding inputs (AYLR, AYLT) are
!    excluded from MELANGE8 (and EOSFI7 changed to EOSFI8 accordingly).
!  ELECT7 upgraded to ELECT9 (high-degeneracy behaviour is improved)
!!   P O S T - P U B L I C A T I O N    (2 0 1 0 +)   IMPROVEMENTS :
!  ELECT9 upgraded (smooth match of two fits at chi >> 1)
!  BLIN8 replaced by BLIN9 - smooth fit interfaces at chi=0.6 and 14.
!  MELANGE8 replaced by MELANGE9 - slightly modified input/output
! 08.08.11 - corrected mistake (unlikely to have an effect) in CHEMFIT7
! 16.11.11 - ELECT9 upgraded to ELECT11 (additional output)
! 20.04.12 - FHARM8 and HLfit8 upgraded to FHARM12 and HLfit12:
!   output of HLfit12 does not include zero-point vibr., but provides U1
! 22.12.12 - MELANGE9 now includes a correction to the linear mixing
!   rule (LMR) for the Madelung energy in the random bcc multi-ion
!   lattice.
! 14.05.13 - an accidental error in programming the newly introduced
!   correction to the LMR is fixed.
! 20.05.13 - calculation of the Wigner-Kirkwood quantum diffraction term
!   for the liquid plasma is moved from EOSFI8 into MELANGE9.
! 10.12.14 - slight cleaning of the text (no effect on the results)
! 28.05.15 - an accidental error in Wigner-Kirkwood entropy correction
!   is fixed (it was in the line "Stot=Stot+FWK*DENSI" since 20.05.13)
! 29.08.15 - eliminated underflow of exp(-THETA) in CHEMFIT7
! 10.08.16 - modified criteria to avoid accuracy loss (round-off errors)
! 07.02.17 - included possibility to switch off the WK (Wigner) terms
! 27.05.17 - safeguard against Zion < 1 is added in FSCRsol8;
!   safeguard against huge (-CHI) values is added in ELECT11.
! 27.01.19 - safeguard against X1=0 in CORMIX.
! 18.04.20 - corrected Wigner-Kirkwood term for heat capacity.
! 04.03.21 - corrected SUBFERMJ: defined parameter EPS (was undefined).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                           MAIN program:               Version 02.06.09
! This driving routine allows one to compile and run this code "as is".
! In practice, however, one usually needs to link subroutines from this
! file to another (external) code, therefore the MAIN program is
! normally commented-out.
      program main
      implicit none
      double precision, parameter :: UN_T6 = .3157746
      integer, parameter :: NMIX = 2
      double precision :: AY(NMIX), AZion(NMIX), ACMI(NMIX)
      double precision :: RHO, RHOlg, T, Tlg, T6, Tnk, TEMP, DENS
      double precision :: Zmean, CMImean, Z2mean, GAMI, P
      double precision :: CHI, TPT, TEGRAD, PRADnkT
      double precision :: PnkT, UNkT, SNk, CV, CHIR, CHIT
      integer :: LIQSOL
      double precision :: x, diff, max_diff, T_arr(3), rho_arr(2)
      integer :: i, j
      AZion(1) = 6.0d0
      AZion(2) = 8.0d0
      ACMI(1) = 12.0d0
      ACMI(2) = 16.0d0
      AY(1) = 0.6d0
      AY(2) = 0.4d0
      T_arr(1) = 1.d9
      T_arr(2) = 5.d9
      T_arr(3) = 1.d6
      rho_arr(1) = 1.d7
      rho_arr(2) = 5.d9

      max_diff = 0.0d0

      do j = 1, 1
         do i = 1, 3
            print *, "iter ", i, j
            T = T_arr(i)
            RHO = RHO_arr(j)
            RHOlg=dlog10(RHO)
            Tlg=dlog10(T)
            T6=10.d0**(Tlg-6.d0)
            RHO=10.d0**RHOlg
            TEMP=T6/UN_T6 ! T [au]
            call MELANGE9(AY,AZion,ACMI,RHO,TEMP, & ! input
                 PRADnkT, & ! additional output - radiative pressure
                 DENS,Zmean,CMImean,Z2mean,GAMI,CHI,TPT,LIQSOL, & ! output param.
                 PnkT,UNkT,SNk,CV,CHIR,CHIT) ! output dimensionless TD functions
            Tnk=8.31447d13/CMImean*RHO*T6 ! n_i kT [erg/cc]
            P=PnkT*Tnk/1.d12 ! P [Mbar]
            TEGRAD=CHIT/(CHIT**2+CHIR*CV/PnkT) ! from Maxwell relat.
            !   --------------------   OUTPUT   --------------------------------   *
            ! Here in the output we have:
            ! RHO - mass density in g/cc
            ! P - total pressure in Mbar (i.e. in 1.e12 dyn/cm^2)
            ! PnkT=P/nkT, where n is the number density of ions, T temperature
            ! CV - heat capacity at constant volume, divided by number of ions, /k
            ! CHIT - logarithmic derivative of pressure \chi_T
            ! CHIR - logarithmic derivative of pressure \chi_\rho
            ! UNkT - internal energy divided by NkT, N being the number of ions
            ! SNk - entropy divided by number of ions, /k
            ! GAMI - ionic Coulomb coupling parameter
            ! TPT=T_p/T, where T_p is the ion plasma temperature
            ! CHI - electron chemical potential, divided by kT
            ! LIQSOL = 0 in the liquid state, = 1 in the solid state

            if (i == 1 .and. j == 1) then
               x = 986087830999.01904d0
            else if (i == 2 .and. j == 1) then
               x = 2495983700684.0181d0
            else if (i == 3 .and. j == 1) then
               x = 826241619577.72607d0
            end if

            diff = abs(x - P) / P
            max_diff = max(diff, max_diff)
            print *, "P DIFF", diff

            if (i == 1 .and. j == 1) then
               x = 16.129464056742833d0
            else if (i == 2 .and. j == 1) then
               x = 8.1653739394820484d0
            else if (i == 3 .and. j == 1) then
               x = 13514.855458323951d0
            end if

            diff = abs(x - PnkT) / PnkT
            max_diff = max(diff, max_diff)
            print *, "PnkT DIFF", diff

            if (i == 1 .and. j == 1) then
               x = 8.5451229292858866d0
            else if (i == 2 .and. j == 1) then
               x = 18.539323243568369d0
            else if (i == 3 .and. j == 1) then
               x = 0.73822827392302692d0
            end if

            diff = abs(x - CV) / CV
            max_diff = max(diff, max_diff)
            print *, "CV DIFF", diff

            if (i == 1 .and. j == 1) then
               x = 0.24165606904443493d0
            else if (i == 2 .and. j == 1) then
               x = 0.88747950206022497d0
            else if (i == 3 .and. j == 1) then
               x = 2.7120648074179433d-5
            end if

            diff = abs(x - CHIT) / CHIT
            max_diff = max(diff, max_diff)
            print *, "CHIT DIFF", diff

            if (i == 1 .and. j == 1) then
               x = 1.3370085960654023d0
            else if (i == 2 .and. j == 1) then
               x = 1.0433031714423413d0
            else if (i == 3 .and. j == 1) then
               x = 1.4524787201645497d0
            end if

            diff = abs(x - CHIR) / CHIR
            max_diff = max(diff, max_diff)
            print *, "CHIR DIFF", diff

            if (i == 1 .and. j == 1) then
               x = 30.712489657322770d0
            else if (i == 2 .and. j == 1) then
               x = 18.110542903803580d0
            else if (i == 3 .and. j == 1) then
               x = 25265.106328521317d0
            end if

            diff = abs(x - UNkT) / UNkT
            max_diff = max(diff, max_diff)
            print *, "UNkT DIFF", diff

            if (i == 1 .and. j == 1) then
               x = 23.797925638433309d0
            else if (i == 2 .and. j == 1) then
               x = 45.817442265862802d0
            else if (i == 3 .and. j == 1) then
               x = 1.0215909624032917d0
            end if

            diff = abs(x - SNk) / SNk
            max_diff = max(diff, max_diff)
            print *, "SNk DIFF", diff

            if (i == 1 .and. j == 1) then
               x = 0.96111630472601972d0
            else if (i == 2 .and. j == 1) then
               x = 0.19172836887561015d0
            else if (i == 3 .and. j == 1) then
               x = 960.24524371490861d0
            end if

            diff = abs(x - GAMI) / GAMI
            max_diff = max(diff, max_diff)
            print *, "GAMI DIFF", diff

            if (i == 1 .and. j == 1) then
               x = 1.2400526419152945d-2
            else if (i == 2 .and. j == 1) then
               x = 2.4705336474828152d-3
            else if (i == 3 .and. j == 1) then
               x = 12.383672318439324d0
            end if

            diff = abs(x - TPT) / TPT
            max_diff = max(diff, max_diff)
            print *, "TPT DIFF", diff

            if (i == 1 .and. j == 1) then
               x = 5.5745494145734744d0
            else if (i == 2 .and. j == 1) then
               x = -0.43436266588208006d0
            else if (i == 3 .and. j == 1) then
               x = 5894.2025691009021d0
            end if

            diff = abs(x - CHI) / CHI
            max_diff = max(diff, max_diff)
            print *, "CHI DIFF", diff

            if (i == 1 .and. j == 1) then
               x = 0
            else if (i == 2 .and. j == 1) then
               x = 0
            else if (i == 3 .and. j == 1) then
               x = 1
            end if

            diff = abs(x - LIQSOL)
            max_diff = max(diff, max_diff)
            print *, "LIQSOL DIFF", diff

         end do
      end do

      print *, "max diff = ", max_diff

      end program main
      
      subroutine MELANGE9(AY,AZion,ACMI,RHO,TEMP,PRADnkT, &
          DENS,Zmean,CMImean,Z2mean,GAMImean,CHI,TPT,LIQSOL, &
          PnkT,UNkT,SNk,CV,CHIR,CHIT)
!                                                       Version 18.04.20
! Difference from v.10.12.14: included switch-off of WK correction
! Stems from MELANGE8 v.26.12.09.
! Difference: output PRADnkT instead of input KRAD
! + EOS of fully ionized electron-ion plasma mixture.     
! Limitations:
! (a) inapplicable in the regimes of
!      (1) bound-state formation,
!      (2) quantum liquid,
!      (3) presence of positrons;
! (b) for the case of a composition gradually depending on RHO or TEMP,
!  second-order functions (CV,CHIR,CHIT in output) should not be trusted
! Choice of the liquid or solid regime - criterion GAMI [because the
!     choice based on comparison of total (non-OCP) free energies can be
!     sometimes dangerous because of the fit uncertainties ("Local field
!     correction" in solid and quantum effects in liquid are unknown)].
! Input: AY - their partial number densities,
!        AZion and ACMI - their charge and mass numbers,
!        RHO - total mass density [g/cc]
!        TEMP - temperature [in a.u.=2Ryd=3.1577e5 K].
! NB: instead of RHO, a true input is CHI, defined below
!     Hence, disagreement between RHO and DENS is the fit error (<0.4%)
! Output:
!         AY - rescaled so that to sum up to 1 and resorted (by AZion)
!         AZion - resorted in ascending order
!         ACMI - resorted in agreement with AZion
!         DENS - electron number density [in a.u.=6.7483346e24 cm^{-3}]
!         Zmean=<Z>, CMImean=<A> - mean ion charge and mass numbers,
!         Z2mean=<Z^2> - mean-square ion charge number
!         GAMImean - effective ion-ion Coulomb coupling constant
!         CHI = mu_e/kT, where mu_e is the electron chem.potential
!         TPT - effective ionic quantum parameter (T_p/T)
!         LIQSOL=0/1 for liquid/solid
!         SNk - dimensionless entropy per 1 ion
!         UNkT - internal energy per kT per ion
!         PnkT - pressure / n_i kT, where n_i is the ion number density
!         PRADnkT - radiative pressure / n_i kT
!         CV - heat capacity per ion, div. by Boltzmann const.
!         CHIR - inverse compressibility -(d ln P / d ln V)_T ("\chi_r")
!         CHIT = (d ln P / d ln T)_V ("\chi_T")
      !implicit double precision (A-H), double precision (O-Z)
      implicit none
      save
      integer, parameter :: NMIX = 2

      double precision, intent(in) :: RHO, TEMP
      double precision, intent(in) :: AY(NMIX), AZion(NMIX), ACMI(NMIX)
      double precision, intent(inout) :: DENS, Zmean, Z2mean, GAMImean
      double precision, intent(inout) :: CHI, TPT
      integer, intent(inout) :: LIQSOL
      double precision, intent(inout) :: SNk, UnkT, PnkT, PRADnkT
      double precision, intent(inout) :: CV, CHIR, CHIT

      double precision, parameter :: CWK = 1.d0 ! Turn on Wigner corrections
      double precision, parameter :: TINY = 1.d-7
      double precision, parameter :: PI = 3.141592653d0
      double precision, parameter :: C53 = 5.d0/3.d0
      double precision, parameter :: C13 = 1.d0/3.d0
      double precision, parameter :: AUM=1822.888d0      ! a.m.u./m_e
      double precision, parameter :: GAMIMELT=175. ! OCP value of Gamma_i for melting
      double precision, parameter :: RSIMELT=140. ! ion density parameter of quantum melting
      double precision, parameter :: RAD=2.554d-7 ! Radiation constant (=4\sigma/c) (in a.u.)
      double precision :: Z52, Z53, Z73, Z321, CMImean, CMI
      double precision :: Zion, Z13, X, X1, X2
      double precision :: UWK, UINTRAD, UMIX, UINTE, UINT, UEid, UC2,UC1
      double precision :: CHIRE, CHITE, CTP, CV1, CV2, CVE, CVMIX, CVtot
      double precision :: DeltaG, DENSI, DNI, DTE, FC1, FC2, FEid, FMIX
      double precision :: DlnDH, DlnDT, DlnDHH, DlnDHT, DlnDTT
      double precision :: FWK, GAME, GAMI
      integer :: i, ix, j
      double precision :: PC1, PC2, PDLR, PDLT, PDR1, PDR2, PDRMIX
      double precision :: PDT1, PDT2, PDTMIX, PEid, PMIX, PRESS, PRESSE
      double precision :: PRESSI, PRESSRAD, PRI, RS, RSI, RZ, SC1, SC2
      double precision :: SEid, Stot, TPT2
      interface
         subroutine chemfit(dens, temp, chi) bind(C, name='chemfit')
           implicit none
           double precision, intent(in), value :: dens, temp
           double precision, intent(inout) :: chi
         end subroutine chemfit
         subroutine elect11(TEMP,CHI, &
              DENS,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE, &
              DlnDH,DlnDT,DlnDHH,DlnDTT,DlnDHT) bind(C, name="elect11")
           implicit none
           double precision, intent(in), value :: TEMP,CHI
           double precision :: DENS,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE, &
                DlnDH,DlnDT,DlnDHH,DlnDTT,DlnDHT
         end subroutine elect11
      end interface
      if (RHO.lt.1.e-19.or.RHO.gt.1.e15) then
         print *, 'MELANGE: RHO out of range'
         stop
      end if
      ! Calculation of average values:
      Zmean=0.
      Z2mean=0.
      Z52=0.
      Z53=0.
      Z73=0.
      Z321=0. ! corr.26.12.09
      CMImean=0.
      do IX=1,NMIX
         Zmean=Zmean+AY(IX)*AZion(IX)
         Z2mean=Z2mean+AY(IX)*AZion(IX)**2
         Z13=AZion(IX)**C13
         Z53=Z53+AY(IX)*Z13**5 
         Z73=Z73+AY(IX)*Z13**7
         Z52=Z52+AY(IX)*dsqrt(AZion(IX))**5
         Z321=Z321+AY(IX)*AZion(IX)*dsqrt(AZion(IX)+1.d0)**3 ! 26.12.09
         CMImean=CMImean+AY(IX)*ACMI(IX)
      enddo
      ! (0) Photons:
      UINTRAD=RAD*TEMP**4
      PRESSRAD=UINTRAD/3.
      ! (1) ideal electron gas (including relativity and degeneracy)
      DENS=RHO/11.20587*Zmean/CMImean ! number density of electrons [au]
      call CHEMFIT(DENS,TEMP,CHI)
      ! NB: CHI can be used as true input instead of RHO or DENS
      call ELECT11(TEMP,CHI, &
       DENS,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE, &
       DlnDH,DlnDT,DlnDHH,DlnDTT,DlnDHT)
      ! NB: at this point DENS is redefined (the difference can be ~0.1%)
      DTE=DENS*TEMP
      PRESSE=PEid*DTE ! P_e [a.u.]
      UINTE=UEid*DTE ! U_e / V [a.u.]
      ! (2) non-ideal Coulomb EIP
      RS=(.75d0/PI/DENS)**C13 ! r_s - electron density parameter
      RSI=RS*CMImean*Z73*AUM ! R_S - ion density parameter
      GAME=1.d0/RS/TEMP ! electron Coulomb parameter Gamma_e
      GAMImean=Z53*GAME   ! effective Gamma_i - ion Coulomb parameter
      if (GAMImean.lt.GAMIMELT.or.RSI.lt.RSIMELT) then
         LIQSOL=0 ! liquid regime
      else
         LIQSOL=1 ! solid regime
      endif
      ! Calculate partial thermodynamic quantities and combine them together:
      UINT=UINTE
      PRESS=PRESSE
      CVtot=CVE*DENS
      Stot=SEid*DENS
      PDLT=PRESSE*CHITE ! d P_e[a.u.] / d ln T
      PDLR=PRESSE*CHIRE ! d P_e[a.u.] / d ln\rho
      DENSI=DENS/Zmean ! number density of all ions
      PRESSI=DENSI*TEMP ! ideal-ions total pressure (normalization)
      TPT2=0.
      CTP=4.d0*PI/AUM/TEMP**2 ! common coefficient for TPT2.10.12.14
      ! Add Coulomb+xc nonideal contributions, and ideal free energy:
      do IX=1,NMIX
        if (AY(IX).ge.TINY) then
         Zion=AZion(IX)
         CMI=ACMI(IX)
         GAMI=Zion**C53*GAME ! Gamma_i for given ion species
         DNI=DENSI*AY(IX) ! number density of ions of given type
         PRI=DNI*TEMP ! = ideal-ions partial pressure (normalization)
         call EOSFI8(LIQSOL,CMI,Zion,RS,GAMI, &
          FC1,UC1,PC1,SC1,CV1,PDT1,PDR1, &
          FC2,UC2,PC2,SC2,CV2,PDT2,PDR2)
         ! First-order TD functions:
         UINT=UINT+UC2*PRI ! internal energy density (e+i+Coul.)
         Stot=Stot+DNI*(SC2-dlog(AY(IX))) !entropy per unit volume[a.u.]
         PRESS=PRESS+PC2*PRI ! pressure (e+i+Coul.) [a.u.]
         ! Second-order functions (they take into account compositional changes):
         CVtot=CVtot+DNI*CV2 ! C_V (e+i+Coul.)/ V (optim.10.12.14)
         PDLT=PDLT+PRI*PDT2 ! d P / d ln T
         PDLR=PDLR+PRI*PDR2 ! d P / d ln\rho
         TPT2=TPT2+CTP*DNI/ACMI(IX)*AZion(IX)**2 ! opt.10.12.14
       end if
      enddo ! next IX
      ! Wigner-Kirkwood perturbative correction for liquid:
      TPT=dsqrt(TPT2) ! effective T_p/T - ion quantum parameter
      ! (in the case of a mixture, this estimate is crude)
      if (LIQSOL.eq.0) then
         FWK=TPT2/24.d0*CWK ! Wigner-Kirkwood (quantum diffr.) term
        if (FWK.gt..7.and.CWK.gt.0.) then
           print*,'MELANGE9: strong quantum effects in liquid!'
           read(*,'(A)')
        endif
         UWK=2.d0*FWK
         UINT=UINT+UWK*PRESSI
         Stot=Stot+FWK*DENSI ! corrected 28.05.15
         PRESS=PRESS+FWK*PRESSI
         CVtot=CVtot-UWK*DENSI ! corrected 18.04.20
         PDLT=PDLT-FWK*PRESSI
         PDLR=PDLR+UWK*PRESSI
      endif
      ! Corrections to the linear mixing rule:
      if (LIQSOL.eq.0) then ! liquid phase
         call CORMIX(RS,GAME,Zmean,Z2mean,Z52,Z53,Z321, &
          FMIX,UMIX,PMIX,CVMIX,PDTMIX,PDRMIX)
      else ! solid phase (only Madelung contribution) [22.12.12]
         FMIX=0.
        do I=1,NMIX
          do J=I+1,NMIX
             RZ=AZion(J)/AZion(I)
             X2=AY(J)/(AY(I)+AY(J))
             X1=dim(1.d0,X2)
             if (X1.lt.TINY) then
                cycle ! 27.01.19
             end if
             if (X2.lt.TINY) then
                cycle
             end if
             X=X2/RZ+(1.d0-1.d0/RZ)*X2**RZ
             GAMI=AZion(I)**C53*GAME ! Gamma_i corrected 14.05.13
             DeltaG=.012*(1.d0-1.d0/RZ**2)*(X1+X2*RZ**C53)
             DeltaG=DeltaG*X/X2*dim(1.d0,X)/X1
             FMIX=FMIX+AY(I)*AY(J)*GAMI*DeltaG
          enddo
        enddo
         UMIX=FMIX
         PMIX=FMIX/3.d0
         CVMIX=0.
         PDTMIX=0.
         PDRMIX=FMIX/2.25d0
      endif
      UINT=UINT+UMIX*PRESSI
      Stot=Stot+DENSI*(UMIX-FMIX)
      PRESS=PRESS+PMIX*PRESSI
      CVtot=CVtot+DENSI*CVMIX
      PDLT=PDLT+PRESSI*PDTMIX
      PDLR=PDLR+PRESSI*PDRMIX
      ! First-order:
      PRADnkT=PRESSRAD/PRESSI ! radiative pressure / n_i k T
      PnkT=PRESS/PRESSI ! P / n_i k T
      UNkT=UINT/PRESSI ! U / N_i k T
      SNk=Stot/DENSI ! S / N_i k
      ! Second-order:
      CV=CVtot/DENSI ! C_V per ion
      CHIR=PDLR/PRESS ! d ln P / d ln\rho
      CHIT=PDLT/PRESS ! d ln P / d ln T
      return
      end

      subroutine EOSFI8(LIQSOL,CMI,Zion,RS,GAMI, &
          FC1,UC1,PC1,SC1,CV1,PDT1,PDR1, &
          FC2,UC2,PC2,SC2,CV2,PDT2,PDR2)
!                                                       Version 16.09.08
!                 call FHARM8 has been replaced by call FHARM12 27.04.12
!                           Wigner-Kirkwood correction excluded 20.05.13
!                                               slight cleaning 10.12.14
! Non-ideal parts of thermodynamic functions in the fully ionized plasma
! Stems from EOSFI5 and EOSFI05 v.04.10.05
! Input: LIQSOL=0/1(liquid/solid), 
!        Zion,CMI - ion charge and mass numbers,
!        RS=r_s (electronic density parameter),
!        GAMI=Gamma_i (ion coupling),
! Output: FC1 and UC1 - non-ideal "ii+ie+ee" contribution to the 
!         free and internal energies (per ion per kT),
!         PC1 - analogous contribution to pressure divided by (n_i kT),
!         CV1 - "ii+ie+ee" heat capacity per ion [units of k]
!         PDT1=(1/n_i kT)*(d P_C/d ln T)_V
!         PDR1=(1/n_i kT)*(d P_C/d ln\rho)_T
! FC2,UC2,PC2,SC2,CV2 - analogous to FC1,UC1,PC1,SC1,CV1, but including
! the part corresponding to the ideal ion gas. This is useful for 
! preventing accuracy loss in some cases (e.g., when SC2 << SC1).
! FC2 does not take into account the entropy of mixing S_{mix}: in a
! mixture, S_{mix}/(N_i k) has to be added externally (see MELANGE9).
! FC2 does not take into account the ion spin degeneracy either.
! When needed, the spin term must be added to the entropy externally.
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(C53=5.d0/3.d0,C76=7.d0/6.d0) ! TINY excl.10.12.14
      parameter (AUM=1822.888d0) ! a.m.u/m_e
      interface
         subroutine excor7(RS,GAME,FXC,UXC,PXC,CVXC,SXC,PDTXC,PDRXC) bind(C, name="excor7")
           implicit none
           double precision, intent(in), value :: RS, GAME
           double precision :: FXC,UXC,PXC,CVXC,SXC,PDTXC,PDRXC
         end subroutine excor7
      end interface

      if (LIQSOL.ne.1.and.LIQSOL.ne.0) then
         print *, 'EOSFI8: invalid LIQSOL'
         stop
      end if
      if (CMI.le..1) then
         print *, 'EOSFI8: too small CMI'
         stop
      end if
      if (Zion.le..1) then
         print *, 'EOSFI8: too small Zion'
         stop
      end if
      if (RS.le..0) then
         print *, 'EOSFI8: invalid RS'
         stop
      end if
      if (GAMI.le..0) then
         print *, 'EOSFI8: invalid GAMI'
         stop
      end if
      GAME=GAMI/Zion**C53
      call EXCOR7(RS,GAME,FXC,UXC,PXC,CVXC,SXC,PDTXC,PDRXC) ! "ee"("xc")
! Calculate "ii" part:
      COTPT=dsqrt(3.d0/AUM/CMI)/Zion**C76 ! auxiliary coefficient
      TPT=GAMI/dsqrt(RS)*COTPT            ! = T_p/T in the OCP
      FidION=1.5*dlog(TPT**2/GAMI)-1.323515
! 1.3235=1+0.5*ln(6/pi); FidION = F_{id.ion gas}/(N_i kT), but without
! the term x_i ln x_i = -S_{mix}/(N_i k).
      if (LIQSOL.eq.0) then               ! liquid
         call FITION9(GAMI, &
          FION,UION,PION,CVii,PDTii,PDRii)
         FItot=FION+FidION
         UItot=UION+1.5
         PItot=PION+1.d0
         CVItot=CVii+1.5d0
         SCItot=UItot-FItot
         PDTi=PDTii+1.d0
         PDRi=PDRii+1.d0
      else                                  ! solid
         call FHARM12(GAMI,TPT, &
          Fharm,Uharm,Pharm,CVharm,Sharm,PDTharm,PDRharm) ! harm."ii"
         call ANHARM8(GAMI,TPT,Fah,Uah,Pah,CVah,PDTah,PDRah) ! anharm.
         FItot=Fharm+Fah
         FION=FItot-FidION
         UItot=Uharm+Uah
         UION=UItot-1.5d0 ! minus 1.5=ideal-gas, in order to get "ii"
         PItot=Pharm+Pah
         PION=PItot-1.d0 ! minus 1=ideal-gas
         PDTi=PDTharm+PDTah
         PDRi=PDRharm+PDRah
         PDTii=PDTi-1.d0 ! minus 1=ideal-gas
         PDRii=PDRi-1.d0 ! minus 1=ideal-gas
         CVItot=CVharm+CVah
         SCItot=Sharm+Uah-Fah
         CVii=CVItot-1.5d0 ! minus 1.5=ideal-gas
      endif
! Calculate "ie" part:
      if (LIQSOL.eq.1) then
         call FSCRsol8(RS,GAMI,Zion,TPT, &
          FSCR,USCR,PSCR,S_SCR,CVSCR,PDTSCR,PDRSCR)
      else
         call FSCRliq8(RS,GAME,Zion, &
          FSCR,USCR,PSCR,CVSCR,PDTSCR,PDRSCR)
         S_SCR=USCR-FSCR
      endif
! Total excess quantities ("ii"+"ie"+"ee", per ion):
      FC0=FSCR+Zion*FXC
      UC0=USCR+Zion*UXC
      PC0=PSCR+Zion*PXC
      SC0=S_SCR+Zion*SXC
      CV0=CVSCR+Zion*CVXC
      PDT0=PDTSCR+Zion*PDTXC
      PDR0=PDRSCR+Zion*PDRXC
      FC1=FION+FC0
      UC1=UION+UC0
      PC1=PION+PC0
      SC1=(UION-FION)+SC0
      CV1=CVii+CV0
      PDT1=PDTii+PDT0
      PDR1=PDRii+PDR0
! Total excess + ideal-ion quantities
      FC2=FItot+FC0
      UC2=UItot+UC0
      PC2=PItot+PC0
      SC2=SCItot+SC0
      CV2=CVItot+CV0
      PDT2=PDTi+PDT0
      PDR2=PDRi+PDR0
      return
      end

! ==================  ELECTRON-ION COULOMB LIQUID  =================== !
      subroutine FITION9(GAMI,FION,UION,PION,CVii,PDTii,PDRii)
!                                                       Version 11.09.08
! Dummy argument Zion is deleted in 2009.
! Non-ideal contributions to thermodynamic functions of classical OCP.
!   Stems from FITION00 v.24.05.00.
! Input: GAMI - ion coupling parameter
! Output: FION - ii free energy / N_i kT
!         UION - ii internal energy / N_i kT
!         PION - ii pressure / n_i kT
!         CVii - ii heat capacity / N_i k
!         PDTii = PION + d(PION)/d ln T = (1/N_i kT)*(d P_{ii}/d ln T)
!         PDRii = PION + d(PION)/d ln\rho
!   Parameters adjusted to Caillol (1999).
      implicit double precision (A-H),double precision (O-Z)
      save
      parameter (A1=-.907347d0,A2=.62849d0,C1=.004500d0,G1=170.0, &
       C2=-8.4d-5,G2=.0037,SQ32=.8660254038d0) ! SQ32=sqrt(3)/2
      A3=-SQ32-A1/dsqrt(A2)
      F0=A1*(dsqrt(GAMI*(A2+GAMI))- &
          A2*dlog(dsqrt(GAMI/A2)+dsqrt(1.+GAMI/A2)))+ &
          2.*A3*(dsqrt(GAMI)-datan(dsqrt(GAMI)))
      U0=dsqrt(GAMI)**3*(A1/dsqrt(A2+GAMI)+A3/(1.d0+GAMI))
!   This is the zeroth approximation. Correction:
      UION=U0+C1*GAMI**2/(G1+GAMI)+C2*GAMI**2/(G2+GAMI**2)
      FION=F0+C1*(GAMI-G1*dlog(1.d0+GAMI/G1))+ &
        C2/2.*dlog(1.d0+GAMI**2/G2)
      CVii=-0.5*dsqrt(GAMI)**3*(A1*A2/dsqrt(A2+GAMI)**3+ &
       A3*(1.d0-GAMI)/(1.d0+GAMI)**2) - &
       GAMI**2*(C1*G1/(G1+GAMI)**2+C2*(G2-GAMI**2)/(G2+GAMI**2)**2)
      PION=UION/3.
      PDRii=(4.*UION-CVii)/9. ! p_{ii} + d p_{ii} / d ln\rho
      PDTii=CVii/3. ! p_{ii} + d p_{ii} / d ln T
      return
      end

      subroutine FSCRliq8(RS,GAME,Zion, &
          FSCR,USCR,PSCR,CVSCR,PDTSCR,PDRSCR) ! fit to the el.-ion scr.
!                                                       Version 11.09.08
!                                                       cleaned 16.06.09
! Stems from FSCRliq7 v. 09.06.07. Included a check for RS=0.
!   INPUT: RS - density parameter, GAME - electron Coulomb parameter,
!          Zion - ion charge number,
!   OUTPUT: FSCR - screening (e-i) free energy per kT per 1 ion,
!           USCR - internal energy per kT per 1 ion (screen.contrib.)
!           PSCR - pressure divided by (n_i kT) (screen.contrib.)
!           CVSCR - heat capacity per 1 ion (screen.contrib.)
!           PDTSCR,PDRSCR = PSCR + d PSCR / d ln(T,\rho)
      implicit double precision(A-H),double precision(O-Z)
      save
      parameter(XRS=.0140047,TINY=1.d-19)
      if (RS.lt.0.) then
         print *, 'FSCRliq8: RS < 0'
         stop
      end if
      if (RS.lt.TINY) then
         FSCR=0.
         USCR=0.
         PSCR=0.
         CVSCR=0.
         PDTSCR=0.
         PDRSCR=0.
         return
      endif
      SQG=sqrt(GAME)
      SQR=sqrt(RS)
      SQZ1=dsqrt(1.+Zion)
      SQZ=dsqrt(Zion)
      CDH0=Zion/1.73205 ! 1.73205=sqrt(3.)
      CDH=CDH0*(SQZ1**3-SQZ**3-1.)
      SQG=sqrt(GAME)
      ZLN=dlog(Zion)
      Z13=exp(ZLN/3.) ! Zion**(1./3.)
      X=XRS/RS ! relativity parameter
      CTF=Zion**2*.2513*(Z13-1.+.2/sqrt(Z13))
! Thomas-Fermi constant; .2513=(18/175)(12/\pi)^{2/3}
      P01=1.11*exp(.475*ZLN)
      P03=0.2+0.078*ZLN**2
      PTX=1.16+.08*ZLN
      TX=GAME**PTX
      TXDG=PTX*TX/GAME
      TXDGG=(PTX-1.)*TXDG/GAME
      TY1=1./(1.d-3*Zion**2+2.*GAME)
      TY1DG=-2.*TY1**2
      TY1DGG=-4.*TY1*TY1DG
      TY2=1.+6.*RS**2
      TY2DX=-12.*RS**2/X
      TY2DXX=-3.*TY2DX/X
      TY=RS**3/TY2*(1.+TY1)
      TYX=3./X+TY2DX/TY2
      TYDX=-TY*TYX
      TYDG=RS**3*TY1DG/TY2
      P1=(Zion-1.)/9.
      COR1=1.+P1*TY
      COR1DX=P1*TYDX
      COR1DG=P1*TYDG
      COR1DXX=P1*(TY*(3./X**2+(TY2DX/TY2)**2-TY2DXX/TY2)-TYDX*TYX)
      COR1DGG=P1*RS**3*TY1DGG/TY2
      COR1DXG=-P1*TYDG*TYX
      U0=.78*sqrt(GAME/Zion)*RS**3
      U0DX=-3.*U0/X
      U0DG=.5*U0/GAME
      U0DXX=-4.*U0DX/X
      U0DGG=-.5*U0DG/GAME
      U0DXG=-3.*U0DG/X
      D0DG=Zion**3
      D0=GAME*D0DG+21.*RS**3
      D0DX=-63.*RS**3/X
      D0DXX=252.*RS**3/X**2
      COR0=1.+U0/D0
      COR0DX=(U0DX-U0*D0DX/D0)/D0
      COR0DG=(U0DG-U0*D0DG/D0)/D0
      COR0DXX=(U0DXX-(2.*U0DX*D0DX+U0*D0DXX)/D0+2.*(D0DX/D0)**2)/D0
      COR0DGG=(U0DGG-2.*U0DG*D0DG/D0+2.*U0*(D0DG/D0)**2)/D0
      COR0DXG=(U0DXG-(U0DX*D0DG+U0DG*D0DX)/D0+2.*U0*D0DX*D0DG/D0**2)/D0
! Relativism:
      RELE=dsqrt(1.d0+X**2)
      Q1=.18/dsqrt(dsqrt(Zion))
      Q2=.2+.37/dsqrt(Zion)
      H1U=1.+X**2/5.
      H1D=1.+Q1*X+Q2*X**2
      H1=H1U/H1D
      H1X=.4*X/H1U-(Q1+2.*Q2*X)/H1D
      H1DX=H1*H1X
      H1DXX=H1DX*H1X+ &
       H1*(.4/H1U-(.4*X/H1U)**2-2.*Q2/H1D+((Q1+2.*Q2*X)/H1D)**2)
      UP=CDH*SQG+P01*CTF*TX*COR0*H1
      UPDX=P01*CTF*TX*(COR0DX*H1+COR0*H1DX)
      UPDG=.5*CDH/SQG+P01*CTF*(TXDG*COR0+TX*COR0DG)*H1
      UPDXX=P01*CTF*TX*(COR0DXX*H1+2.*COR0DX*H1DX+COR0*H1DXX)
      UPDGG=-.25*CDH/(SQG*GAME)+ &
       P01*CTF*(TXDGG*COR0+2.*TXDG*COR0DG+TX*COR0DGG)*H1
      UPDXG=P01*CTF*(TXDG*(COR0DX*H1+COR0*H1DX)+ &
       TX*(COR0DXG*H1+COR0DG*H1DX))
      DN1=P03*SQG+P01/RS*TX*COR1
      DN1DX=P01*TX*(COR1/XRS+COR1DX/RS)
      DN1DG=.5*P03/SQG+P01/RS*(TXDG*COR1+TX*COR1DG)
      DN1DXX=P01*TX/XRS*(2.*COR1DX+X*COR1DXX)
      DN1DGG=-.25*P03/(GAME*SQG)+ &
       P01/RS*(TXDGG*COR1+2.*TXDG*COR1DG+TX*COR1DGG)
      DN1DXG=P01*(TXDG*(COR1/XRS+COR1DX/RS)+TX*(COR1DG/XRS+COR1DXG/RS))
      DN=1.+DN1/RELE
      DNDX=DN1DX/RELE-X*DN1/RELE**3
      DNDXX=(DN1DXX-((2.*X*DN1DX+DN1)-3.*X**2*DN1/RELE**2)/RELE**2)/RELE
      DNDG=DN1DG/RELE
      DNDGG=DN1DGG/RELE
      DNDXG=DN1DXG/RELE-X*DN1DG/RELE**3
      FSCR=-UP/DN*GAME
      FX=(UP*DNDX/DN-UPDX)/DN
      FXDG=((UPDG*DNDX+UPDX*DNDG+UP*DNDXG-2.*UP*DNDX*DNDG/DN)/DN- &
       UPDXG)/DN
      FDX=FX*GAME
      FG=(UP*DNDG/DN-UPDG)/DN
      FDG=FG*GAME-UP/DN
      FDGDH=SQG*DNDG/DN**2 ! d FDG / d CDH
      FDXX=((UP*DNDXX+2.*(UPDX*DNDX-UP*DNDX**2/DN))/DN-UPDXX)/DN*GAME
      FDGG=2.*FG+GAME*((2.*DNDG*(UPDG-UP*DNDG/DN)+UP*DNDGG)/DN-UPDGG)/DN
      FDXG=FX+GAME*FXDG
      USCR=GAME*FDG
      CVSCR=-GAME**2*FDGG
      PSCR=(X*FDX+GAME*FDG)/3.
      PDTSCR=-GAME**2*(X*FXDG+FDGG)/3.
      PDRSCR=(12.*PSCR+X**2*FDXX+2.*X*GAME*FDXG+GAME**2*FDGG)/9.
      return
      end

! ==============   SUBROUTINES FOR THE SOLID STATE   ================= !
      subroutine FSCRsol8(RS,GAMI,ZNUCL,TPT, &
          FSCR,USCR,PSCR,S_SCR,CVSCR,PDTSCR,PDRSCR)
!                                                       Version 28.05.08
!                    undefined zero variable Q1DXG is wiped out 21.06.10
!                                 accuracy-loss safeguard added 10.08.16
!                              safequard against Zion < 1 added 27.05.17
! Fit to the el.-ion screening in bcc or fcc Coulomb solid
! Stems from FSCRsol8 v.09.06.07. Included a check for RS=0.
!   INPUT: RS - el. density parameter, GAMI - ion coupling parameter,
!          ZNUCL - ion charge, TPT=T_p/T - ion quantum parameter
!   OUTPUT: FSCR - screening (e-i) free energy per kT per 1 ion,
!           USCR - internal energy per kT per 1 ion (screen.contrib.)
!           PSCR - pressure divided by (n_i kT) (screen.contrib.)
!           S_SCR - screening entropy contribution / (N_i k)
!           CVSCR - heat capacity per 1 ion (screen.contrib.)
!           PDTSCR,PDRSCR = PSCR + d PSCR / d ln(T,\rho)
      implicit double precision(A-H),double precision(O-Z)
      save
      dimension AP(4) ! parameters of the fit
      parameter (C13=1.d0/3.d0,ENAT=2.7182818285d0,TINY=1.d-19)
      data AP/1.1866,.684,17.9,41.5/,PX/.205/ ! for bcc lattice
      if (RS.lt.0.) then
         print *, 'FSCRliq8: RS < 0'
         stop
      end if
      if (RS.lt.TINY) then
         FSCR=0.
         USCR=0.
         PSCR=0.
         S_SCR=0.
         CVSCR=0.
         PDTSCR=0.
         PDRSCR=0.
         return
      endif
      Zion=ZNUCL
      if (Zion.lt.1.d0) then ! 27.05.17
         print*,'FSCRsol8 WARNING: Z =',Zion,' < 1.'
         Zion=1.d0
      endif
      XSR=.0140047/RS ! relativity parameter
      Z13=Zion**C13
      P1=.00352*(1.-AP(1)/Zion**.267+.27/Zion)
      P2=1.d0+2.25/Z13* &
      (1.+AP(2)*Zion**5+.222*Zion**6)/(1.+.222*Zion**6)
      ZLN=dlog(Zion)
      Finf=sqrt(P2/XSR**2+1.)*Z13**2*P1 ! The TF limit
      FinfX=-P2/((P2+XSR**2)*XSR)
      FinfDX=Finf*FinfX
      FinfDXX=FinfDX*FinfX-FinfDX*(P2+3.*XSR**2)/((P2+XSR**2)*XSR)
      R1=AP(4)/(1.+ZLN)
      R2=.395*ZLN+.347/Zion/sqrt(Zion)
      R3=1.d0/(1.d0+ZLN*sqrt(ZLN)*.01+.097/Zion**2)
      Q1U=R1+AP(3)*XSR**2
      Q1D=1.d0+R2*XSR**2
      Q1=Q1U/Q1D
      Q1X=2.*XSR*(AP(3)/Q1U-R2/Q1D)
      Q1XDX=Q1X/XSR+4.*XSR**2*((R2/Q1D)**2-(AP(3)/Q1U)**2)
      Q1DX=Q1*Q1X
      Q1DXX=Q1DX*Q1X+Q1*Q1XDX
! New quantum factor, in order to suppress CVSCR at TPT >> 1
      if (TPT.lt.6./PX) then
         Y0=(PX*TPT)**2
         Y0DX=Y0/XSR
         Y0DG=2.*Y0/GAMI
         Y0DXX=0.
         Y0DGG=Y0DG/GAMI
         Y0DXG=Y0DG/XSR
         Y1=dexp(Y0)
         Y1DX=Y1*Y0DX
         Y1DG=Y1*Y0DG
         Y1DXX=Y1*(Y0DX**2+Y0DXX)
         Y1DGG=Y1*(Y0DG**2+Y0DGG)
         Y1DXG=Y1*(Y0DX*Y0DG+Y0DXG)
         SA=1.d0+Y1
         SUPA=dlog(SA)
         SUPADX=Y1DX/SA
         SUPADG=Y1DG/SA
         SUPADXX=(Y1DXX-Y1DX**2/SA)/SA
         SUPADGG=(Y1DGG-Y1DG**2/SA)/SA
         SUPADXG=(Y1DXG-Y1DX*Y1DG/SA)/SA
         EM2=ENAT-2.d0
         SB=ENAT-EM2/Y1
         SUPB=dlog(SB)
         EM2Y1=EM2/(Y1**2*SB)
         SUPBDX=EM2Y1*Y1DX
         SUPBDG=EM2Y1*Y1DG
         SUPBDXX=EM2Y1*(Y1DXX-2.d0*Y1DX**2/Y1-Y1DX*SUPBDX)
         SUPBDGG=EM2Y1*(Y1DGG-2.d0*Y1DG**2/Y1-Y1DG*SUPBDG)
         SUPBDXG=EM2Y1*(Y1DXG-2.d0*Y1DX*Y1DG/Y1-Y1DG*SUPBDX)
         SUP=dsqrt(SUPA/SUPB)
         SUPX=.5d0*(SUPADX/SUPA-SUPBDX/SUPB)
         SUPDX=SUP*SUPX
         SUPG=.5d0*(SUPADG/SUPA-SUPBDG/SUPB)
         SUPDG=SUP*SUPG
         SUPDXX=SUPDX*SUPX+ &
          SUP*.5d0*(SUPADXX/SUPA-(SUPADX/SUPA)**2- &
                  SUPBDXX/SUPB+(SUPBDX/SUPB)**2)
         SUPDGG=SUPDG*SUPG+ &
          SUP*.5d0*(SUPADGG/SUPA-(SUPADG/SUPA)**2- &
                  SUPBDGG/SUPB+(SUPBDG/SUPB)**2)
         SUPDXG=SUPDX*SUPG+ &
          SUP*.5d0*((SUPADXG-SUPADX*SUPADG/SUPA)/SUPA- &
                  (SUPBDXG-SUPBDX*SUPBDG/SUPB)/SUPB)
      else
         SUP=PX*TPT
         SUPDX=.5d0*PX*TPT/XSR
         SUPDG=PX*TPT/GAMI
         SUPDXX=-.5d0*SUPDX/XSR
         SUPDGG=0.
         SUPDXG=SUPDX/GAMI
      endif
      GR3=(GAMI/SUP)**R3
      GR3X=-R3*SUPDX/SUP
      GR3DX=GR3*GR3X
      GR3DXX=GR3DX*GR3X-R3*GR3*(SUPDXX/SUP-(SUPDX/SUP)**2)
      GR3G=R3*(1.d0/GAMI-SUPDG/SUP)
      GR3DG=GR3*GR3G
      GR3DGG=GR3DG*GR3G+GR3*R3*((SUPDG/SUP)**2-SUPDGG/SUP-1.d0/GAMI**2)
      GR3DXG=GR3DG*GR3X+GR3*R3*(SUPDX*SUPDG/SUP**2-SUPDXG/SUP)
      W=1.d0+Q1/GR3
      WDX=Q1DX/GR3-Q1*GR3DX/GR3**2
      WDG=-Q1*GR3DG/GR3**2
      WDXX=Q1DXX/GR3- &
        (2.d0*Q1DX*GR3DX+Q1*(GR3DXX-2.d0*GR3DX**2/GR3))/GR3**2
      WDGG=Q1*(2.d0*GR3DG**2/GR3-GR3DGG)/GR3**2
      WDXG=-(Q1DX*GR3DG+Q1*(GR3DXG-2.d0*GR3DX*GR3DG/GR3))/GR3**2
      FSCR=-GAMI*Finf*W
      FDX=-GAMI*(FinfDX*W+Finf*WDX)
      FDXX=-GAMI*(FinfDXX*W+2.d0*FinfDX*WDX+Finf*WDXX)
      FDG=-Finf*W-GAMI*Finf*WDG
      FDGG=-2.d0*Finf*WDG-GAMI*Finf*WDGG
      if (dabs(FDGG).lt.TINY) FDGG=0. ! 10.08.16: roundoff err.safeguard
      FDXG=-FinfDX*W-Finf*WDX-GAMI*(FinfDX*WDG+Finf*WDXG)
      S_SCR=-GAMI**2*Finf*WDG
      USCR=S_SCR+FSCR
      CVSCR=-GAMI**2*FDGG
      PSCR=(XSR*FDX+GAMI*FDG)/3.d0
      PDTSCR=GAMI**2*(XSR*Finf*(FinfX*WDG+WDXG)-FDGG)/3.d0
      PDRSCR=(12.d0*PSCR+XSR**2*FDXX+2.d0*XSR*GAMI*FDXG+ &
       GAMI**2*FDGG)/9.d0
      return
      end

      subroutine ANHARM8(GAMI,TPT,Fah,Uah,Pah,CVah,PDTah,PDRah)
! ANHARMONIC free energy                                Version 27.07.07
!                                                       cleaned 16.06.09
! Stems from ANHARM8b. Difference: AC=0., B1=.12 (.1217 - over accuracy)
! Input: GAMI - ionic Gamma, TPT=Tp/T - ionic quantum parameter
! Output: anharm.free en. Fah=F_{AH}/(N_i kT), internal energy Uah,
!   pressure Pah=P_{AH}/(n_i kT), specific heat CVah = C_{V,AH}/(N_i k),
!   PDTah = Pah + d Pah / d ln T, PDRah = Pah + d Pah / d ln\rho
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(NM=3)
      dimension AA(NM)
      data AA/10.9,247.,1.765d5/ ! Farouki & Hamaguchi'93
      data B1/.12/ ! coeff.at \eta^2/\Gamma at T=0
      CK=B1/AA(1) ! fit coefficient
      TPT2=TPT**2
      TPT4=TPT2**2
      TQ=B1*TPT2/GAMI ! quantum dependence
      TK2=CK*TPT2
      SUP=dexp(-TK2) ! suppress.factor of class.anharmonicity
      Fah=0.
      Uah=0.
      Pah=0.
      CVah=0.
      PDTah=0.
      PDRah=0.
      SUPGN=SUP
      do N=1,NM
         CN=N
         SUPGN=SUPGN/GAMI ! SUP/Gamma^n
         ACN=AA(N)
         Fah=Fah-ACN/CN*SUPGN
         Uah=Uah+(ACN*(1.+2.*TK2/CN))*SUPGN
         PN=AA(N)/3.+TK2*AA(N)/CN
         Pah=Pah+PN*SUPGN
         CVah=CVah+((CN+1.)*AA(N)+(4.-2./CN)*AA(N)*TK2+ &
           4.*AA(N)*CK**2/CN*TPT4)*SUPGN
         PDTah=PDTah+(PN*(1.+CN+2.*TK2)-2./CN*AA(N)*TK2)*SUPGN
         PDRah=PDRah+(PN*(1.-CN/3.-TK2)+AA(N)/CN*TK2)*SUPGN
      enddo
      Fah=Fah-TQ
      Uah=Uah-TQ
      Pah=Pah-TQ/1.5
      PDRah=PDRah-TQ/4.5
      return
      end

      subroutine FHARM12(GAMI,TPT, &
        Fharm,Uharm,Pharm,CVth,Sth,PDTharm,PDRharm)
! Thermodynamic functions of a harmonic crystal, incl.stat.Coul.lattice
! 
!                                                       Version 27.04.12
! Stems from FHARM8 v.15.02.08
! Replaced HLfit8 with HLfit12: rearranged output.
! Input: GAMI - ionic Gamma, TPT=T_{p,i}/T
! Output: Fharm=F/(N_i T), Uharm=U/(N_i T), Pharm=P/(n_i T),
! CVth=C_V/N_i, Sharm=S/N_i
! PDTharm = Pharm + d Pharm / d ln T, PDRharm = Pharm + d Pharm/d ln\rho
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(CM=.895929256d0) ! Madelung
      call HLfit12(TPT,F,U,CVth,Sth,U1,CW,1)
      U0=-CM*GAMI ! perfect lattice
      E0=1.5d0*U1*TPT ! zero-point energy
      Uth=U+E0
      Fth=F+E0
      Uharm=U0+Uth
      Fharm=U0+Fth
      Pharm=U0/3.d0+Uth/2.d0
      PDTharm=.5d0*CVth
      PDRharm=U0/2.25d0+.75d0*Uth-.25d0*CVth
      return
      end

      subroutine HLfit12(eta,F,U,CV,S,U1,CW,LATTICE)
!                                                       Version 24.04.12
! Stems from HLfit8 v.03.12.08;
!   differences: E0 excluded from  U and F;
!   U1 and d(CV)/d\ln(T) are added on the output.
! Fit to thermal part of the thermodynamic functions.
! Baiko, Potekhin, & Yakovlev (2001).
! Zero-point lattice quantum energy 1.5u_1\eta EXCLUDED (unlike HLfit8).
! Input: eta=Tp/T, LATTICE=1 for bcc, 2 for fcc
! Output: F and U (normalized to NkT) - due to phonon excitations,
!   CV and S (normalized to Nk) in the HL model,
!   U1 - the 1st phonon moment,
!   CW=d(CV)/d\ln(T)
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(EPS=1.d-5,TINY=1.d-99)
      if (LATTICE.eq.1) then ! bcc lattice
         CLM=-2.49389d0 ! 3*ln<\omega/\omega_p>
         U1=.5113875d0
         ALPHA=.265764d0
         BETA=.334547d0
         GAMMA=.932446d0
         A1=.1839d0
         A2=.593586d0
         A3=.0054814d0
         A4=5.01813d-4
         A6=3.9247d-7
         A8=5.8356d-11
         B0=261.66d0
         B2=7.07997d0
         B4=.0409484d0
         B5=.000397355d0
         B6=5.11148d-5
         B7=2.19749d-6
         C9=.004757014d0
         C11=.0047770935d0
      elseif (LATTICE.eq.2) then ! fcc lattice
         CLM=-2.45373d0
         U1=.513194d0
         ALPHA=.257591d0
         BETA=.365284d0
         GAMMA=.9167070d0
         A1=.0
         A2=.532535d0
         A3=.0
         A4=3.76545d-4
         A6=2.63013d-7
         A8=6.6318d-11
         B0=303.20d0
         B2=7.7255d0
         B4=.0439597d0
         B5=.000114295d0
         B6=5.63434d-5
         B7=1.36488d-6
         C9=.00492387d0
         C11=.00437506d0
      else
         print *, 'HLfit: unknown lattice type'
         stop
      endif
      if (eta.gt.1./EPS) then ! asymptote of Eq.(13) of BPY'01
         U=3./(C11*eta**3)
         F=-U/3.
         CV=4.*U
         S=U-F
         return
      elseif (eta.lt.EPS) then ! Eq.(17) of BPY'01
        if (eta.lt.TINY) then
           print *, 'HLfit: eta is too small'
           stop
        end if
         F=3.*dlog(eta)+CLM-1.5*U1*eta+eta**2/24. 
         U=3.-1.5*U1*eta+eta**2/12.
         CV=3.-eta**2/12.
         S=U-F
         return
      endif
      eta2=eta**2
      eta3=eta2*eta
      eta4=eta3*eta
      eta5=eta4*eta
      eta6=eta5*eta
      eta7=eta6*eta
      eta8=eta7*eta
      B9=A6*C9
      B11=A8*C11
      UP=1.+A1*eta+A2*eta2+A3*eta3+A4*eta4+A6*eta6+A8*eta8
      DN=B0+B2*eta2+B4*eta4+B5*eta5+B6*eta6+ &
       B7*eta7+eta8*(B9*eta+B11*eta3)
      EA=dexp(-ALPHA*eta)
      EB=dexp(-BETA*eta)
      EG=dexp(-GAMMA*eta)
      F=dlog(1.d0-EA)+dlog(1.d0-EB)+dlog(1.-EG)-UP/DN ! F_{thermal}/NT
      UP1=A1+ &
      2.*A2*eta+3.*A3*eta2+4.*A4*eta3+6.*A6*eta5+8.*A8*eta7
      UP2=2.*A2+6.*A3*eta+12.*A4*eta2+30.*A6*eta4+56.*A8*eta6
      UP3=6.*A3+24.*A4*eta+120.*A6*eta3+336*A8*eta5
      DN1=2.*B2*eta+4.*B4*eta3+5.*B5*eta4+6.*B6*eta5+ &
       7.*B7*eta6+eta8*(9.*B9+11.*B11*eta2)
      DN2=2.*B2+12.*B4*eta2+20.*B5*eta3+30.*B6*eta4+ &
       42.*B7*eta5+72.*B9*eta7+110.*B11*eta8*eta
      DN3=24.*B4*eta+60.*B5*eta2+120.*B6*eta3+ &
       210.*B7*eta4+504.*B9*eta6+990.*B11*eta8
      DF1=ALPHA*EA/(1.d0-EA)+BETA*EB/(1.d0-EB)+GAMMA*EG/(1.d0-EG)- &
       (UP1*DN-DN1*UP)/DN**2 ! int.en./NT/eta = df/d\eta
      DF2=ALPHA**2*EA/(1.d0-EA)**2+BETA**2*EB/(1.d0-EB)**2+ &
       GAMMA**2*EG/(1.d0-EG)**2+ &
       ((UP2*DN-DN2*UP)*DN-2.*(UP1*DN-DN1*UP)*DN1)/DN**3 ! -d2f/d\eta^2
      U=DF1*eta
      CV=DF2*eta2
      DF3=-ALPHA**3*EA/(1.d0-EA)**3*(1.+EA)- &
       BETA**3*EB/(1.d0-EB)**3*(1.+EB)- &
       GAMMA**3*EG/(1.d0-EG)**3*(1.+EG)+ &
       UP3/DN-(3.*UP2*DN1+3.*UP1*DN2+UP*DN3)/DN**2+ &
       6.*DN1*(UP1*DN1+UP*DN2)/DN**3-6.*UP*DN1**3/DN**4 ! -d3f/d\eta^3
      CW=-2.*CV-eta3*DF3
      S=U-F
      return
      end

      subroutine CORMIX(RS,GAME,Zmean,Z2mean,Z52,Z53,Z321, &
       FMIX,UMIX,PMIX,CVMIX,PDTMIX,PDRMIX)
!                                                       Version 02.07.09
! Correction to the linear mixing rule for moderate to small Gamma
! Input: RS=r_s (if RS=0, then OCP, otherwise EIP)
!        GAME=\Gamma_e
!        Zmean=<Z> (average Z of all ions, without electrons)
!        Z2mean=<Z^2>, Z52=<Z^2.5>, Z53=<Z^{5/3}>, Z321=<Z(Z+1)^1.5>
! Output: FMIX=\Delta f - corr.to the reduced free energy f=F/N_{ion}kT
!         UMIX=\Delta u - corr.to the reduced internal energy u
!         PMIX=\Delta u - corr.to the reduced pressure P=P/n_{ion}kT
!         CVMIX=\Delta c - corr.to the reduced heat capacity c_V
!         PDTMIX=(1/n_{ion}kT)d\Delta P / d ln T
!               = \Delta p +  d \Delta p / d ln T
!         PDRMIX=(1/n_{ion}kT)d\Delta P / d ln n_e
! (composition is assumed fixed: Zmean,Z2mean,Z52,Z53=constant)
      implicit double precision (A-H), double precision (O-Z)
      parameter (TINY=1.d-9)
      GAMImean=GAME*Z53
      if (RS.lt.TINY) then ! OCP
         Dif0=Z52-dsqrt(Z2mean**3/Zmean)
      else
         Dif0=Z321-dsqrt((Z2mean+Zmean)**3/Zmean)
      endif
      DifR=Dif0/Z52
      DifFDH=Dif0*GAME*sqrt(GAME/3.) ! F_DH - F_LM(DH)
      D=Z2mean/Zmean**2
      if (dabs(D-1.d0).lt.TINY) then ! no correction
         FMIX=0.
         UMIX=0.
         PMIX=0.
         CVMIX=0.
         PDTMIX=0.
         PDRMIX=0.
         return
      endif
      P3=D**(-0.2)
      D0=(2.6*DifR+14.*DifR**3)/(1.d0-P3)
      GP=D0*GAMImean**P3
      FMIX0=DifFDH/(1.+GP)
      Q=D**2*.0117
      R=1.5/P3-1.
      GQ=Q*GP
      FMIX=FMIX0/(1.+GQ)**R
      G=1.5-P3*GP/(1.+GP)-R*P3*GQ/(1.+GQ)
      UMIX=FMIX*G
      PMIX=UMIX/3.d0
      GDG=-P3**2*(GP/(1.d0+GP)**2+R*GQ/(1.d0+GQ)**2) ! d G /d ln Gamma
      UDG=UMIX*G+FMIX*GDG ! d u_mix /d ln Gamma
      CVMIX=UMIX-UDG
      PDTMIX=PMIX-UDG/3.
      PDRMIX=PMIX+UDG/9.
      return
      end
