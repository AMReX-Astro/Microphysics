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
         subroutine CORMIX(RS,GAME,Zmean,Z2mean,Z52,Z53,Z321, &
              FMIX,UMIX,PMIX,CVMIX,PDTMIX,PDRMIX) bind(C, name="cormix")
           implicit none
           double precision, value :: RS,GAME,Zmean,Z2mean,Z52,Z53,Z321
           double precision :: FMIX,UMIX,PMIX,CVMIX,PDTMIX,PDRMIX
         end subroutine CORMIX
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
         subroutine fscrsol8(RS,GAMI,Zion,TPT, &
              FSCR,USCR,PSCR,S_SCR,CVSCR,PDTSCR,PDRSCR) bind(C, name="fscrsol8")
           implicit none
           double precision, intent(in), value :: RS, GAMI, Zion, TPT
           double precision :: FSCR,USCR,PSCR,S_SCR,CVSCR,PDTSCR,PDRSCR
         end subroutine fscrsol8
         subroutine anharm8(GAMI,TPT,Fah,Uah,Pah,CVah,PDTah,PDRah) bind(C, name="anharm8")
           implicit none
           double precision, intent(in), value :: GAMI,TPT
           double precision :: Fah,Uah,Pah,CVah,PDTah,PDRah
         end subroutine anharm8
         subroutine fharm12(GAMI,TPT, &
              Fharm,Uharm,Pharm,CVharm,Sharm,PDTharm,PDRharm) bind(C, name="fharm12")
           implicit none
           double precision, intent(in), value :: GAMI,TPT
           double precision :: Fharm,Uharm,Pharm,CVharm,Sharm,PDTharm,PDRharm
         end subroutine fharm12
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
