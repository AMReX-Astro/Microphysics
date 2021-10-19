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
         subroutine EOSFI8(LIQSOL,CMI,Zion,RS,GAMI, &
              FC1,UC1,PC1,SC1,CV1,PDT1,PDR1, &
              FC2,UC2,PC2,SC2,CV2,PDT2,PDR2) bind(C, name="eosfi8")
           implicit none
           integer, value :: LIQSOL
           double precision, value :: CMI,Zion,RS,GAMI
           double precision :: FC1,UC1,PC1,SC1,CV1,PDT1,PDR1, &
                FC2,UC2,PC2,SC2,CV2,PDT2,PDR2
         end subroutine EOSFI8
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
