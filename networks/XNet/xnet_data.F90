!***************************************************************************************************
! data.f90 10/18/17
! This file contains the nuclear and reaction data structures and the subroutines which read in the
! data and allocate the arrays.
!***************************************************************************************************

Module nuc_number
  !-------------------------------------------------------------------------------------------------
  ! This module contains ny, the number of nuclear species whose abundances are evolved by the
  ! network. The value of ny is read in by read_reaction_data or read_nuclear_data.
  !-------------------------------------------------------------------------------------------------
  Implicit None
  Integer :: ny
End Module nuc_number

Module nuclear_data
  !-------------------------------------------------------------------------------------------------
  ! This module contains the essential data for each included species. Their array sizes are set in
  ! the routine read_nuclear_data.
  !-------------------------------------------------------------------------------------------------
  Use xnet_types, Only: dp
  Implicit None
  Character(5), Allocatable :: nname(:)            ! Nuclei names (e.g. he4) (nname(0) is placeholder for non-nuclei)
  Real(dp), Allocatable     :: aa(:), zz(:), nn(:) ! Mass (aa), proton (zz), and neutron (nn) numbers
  Real(dp), Allocatable     :: be(:), mex(:)       ! Binding energy and mass excess [MeV c^{-2}]
  Real(dp), Allocatable     :: mm(:)               ! Mass of nuclei [g]
  Integer, Allocatable      :: ia(:), iz(:), in(:) ! Integer copies of mass, proton, and neutron numbers
  Integer                   :: inmin, inmax        ! Min and max neutron numbers
  Integer                   :: izmin, izmax        ! Min and max proton numbers

  ! Commonly used powers of Z
  Real(dp), Allocatable     :: zz2(:), zz53(:), zzi(:) ! zz^2, zz^{5/3}, and zz^{3b-1}

  ! Get the neutron and proton mass excess (and masses) from netwinv for consistency
  Integer  :: ineut, iprot        ! Indices for neutron and proton in network
  Real(dp) :: mex_n ! = 8.0713179 ! Neutron mass excess [MeV c^{-2}]
  Real(dp) :: mex_p ! = 7.2889848 ! Proton mass excess [MeV c^{-2}]
End Module nuclear_data

Module part_funct_data
  !-------------------------------------------------------------------------------------------------
  ! This module contains the nuclear partition function data. g is the interpolation data (ng,ny),
  ! gg is the current parition function, and angm is the J. gg(0) and angm(0) are placeholders for
  ! non-nuclei. The array sizes are set in read_nuclear_data.
  !-------------------------------------------------------------------------------------------------
  Use xnet_types, Only: dp
  Implicit None
  Integer, Parameter    :: ng = 24      ! Number of grid points for partition function data
  Real(dp)              :: t9i(ng)      ! Temperature grid for partition function data
  Real(dp), Allocatable :: g(:,:)       ! Partition function data
  Real(dp), Allocatable :: gg(:,:)      ! Interpolated partition function
  Real(dp), Allocatable :: angm(:)      ! Angular momentum
  Real(dp), Allocatable :: dlngdt9(:,:) ! d(ln(partition functions))/dT9
  !$omp threadprivate(gg,dlngdt9)

End Module part_funct_data

Module cross_sect_data
  !-------------------------------------------------------------------------------------------------
  ! This module contains the data needed to calculate the cross sections.
  ! The csect variables are the results, velocity integrated cross section * density dependence.
  !-------------------------------------------------------------------------------------------------
  Use xnet_types, Only: dp
  Implicit None
  Integer               :: nreac(3)                        ! # of reactions with i reactants
  Integer, Allocatable  :: n1i(:,:), n2i(:,:), n3i(:,:)    ! List of nuclei affected by each reaction
  Real(dp), Allocatable :: rc1(:,:), rc2(:,:), rc3(:,:)    ! REACLIB parameters
  Real(dp), Allocatable :: csect1(:,:), csect2(:,:), csect3(:,:)
  Real(dp), Allocatable :: dcsect1dt9(:,:), dcsect2dt9(:,:), dcsect3dt9(:,:)
  !$omp threadprivate(csect1,csect2,csect3,dcsect1dt9,dcsect2dt9,dcsect3dt9)

  ! Reaction flags to indicate variations in how the rate is calculated
  Integer, Allocatable :: iwk1(:), iwk2(:), iwk3(:)    ! Weak Reaction
  Integer, Allocatable :: ires1(:), ires2(:), ires3(:) ! Resonant Reaction
  Integer, Allocatable :: irev1(:), irev2(:), irev3(:) ! Reverse Reaction

  ! Additional rate or energy factors
  Real(dp), Allocatable :: q1(:), q2(:), q3(:)          ! Reaction Q values
  Real(dp), Allocatable :: rpf1(:), rpf2(:)             ! Partition function ratios for detailed balance
  Real(dp), Allocatable :: dlnrpf1dt9(:), dlnrpf2dt9(:) ! d(ln(partition function ratios))/dT9
  !$omp threadprivate(rpf1,rpf2,dlnrpf1dt9,dlnrpf2dt9)

  !-------------------------------------------------------------------------------------------------
  ! Aside from the REACLIB formated data, this dataset includes pointers to sets of external
  ! reaction rates, in the form of indices encoded in rc{1,2,3}(1).
  !-------------------------------------------------------------------------------------------------

  ! FFN are electron and positron capture rates encoded in the data format of Fuller, Fowler & Neuman (1982,1985).
  Integer              :: nffn    ! The number of FFN reactions
  Integer, Allocatable :: iffn(:) ! Pointers to FFN list

  ! NNU are neutrino and antineutrino capture rates from Zinner & Langanke, implemented by Carla Froehlich.
  Integer              :: nnnu    ! The number of NNU reactions
  Integer, Allocatable :: innu(:) ! Pointers to NNU list
End Module cross_sect_data

Module screening_data
  !-------------------------------------------------------------------------------------------------
  ! This module contains data used to calculate the screening corrections that appear in the
  ! reaction rate exponents.
  !-------------------------------------------------------------------------------------------------
  Use xnet_types, Only: dp
  Implicit None
  Real(dp), Allocatable :: h1(:,:), h2(:,:), h3(:,:)             ! Screening factors
  Real(dp), Allocatable :: dh1dt9(:,:), dh2dt9(:,:), dh3dt9(:,:) ! d(screening factors)/dT9
  !$omp threadprivate(h1,h2,h3,dh1dt9,dh2dt9,dh3dt9)

  ! Proton (charge) numbers for individual reactants in 2- and 3-reactant reactions
  Real(dp), Allocatable :: z21(:), z22(:), z31(:), z32(:), z33(:)
  Integer, Allocatable  :: iz21(:), iz22(:), iz31(:), iz32(:), iz33(:)

  ! Composite proton (charge) numbers for reactants in 2- and 3-reactant reactions
  Real(dp), Allocatable :: z2c(:), z3c(:)
  Integer, Allocatable  :: iz2c(:), iz3c(:)

  ! Screening factors from Table 4 of Graboske+ (1973)
  Real(dp), Parameter   :: bw = 1.0,  kbw = 0.5            ! Weak screening parameters
  Real(dp), Parameter   :: bi = 0.86, kbi = 0.38           ! Intermediate screening parmaeters
  Real(dp), Parameter   :: bip1 = 1.86                     ! bi + 1
  Real(dp), Parameter   :: thbim1 = 1.58                   ! 3*bi - 1
  Real(dp), Parameter   :: thbim2 = 0.58                   ! 3*bi - 2
  Real(dp), Parameter   :: twm2bi = 0.28                   ! 2 - 2*bi
  Real(dp), Allocatable :: zeta2w(:), zeta2i(:), zeta3w(:), zeta3i(:) ! Reaction charge parameter

  ! Strong screening fitting coefficients from DeWitt & Slattery (2003), Eq. 4:
  !   f(gamma) = a*gamma + (1/s)*b*gamma^s + c*ln(gamma) + d
  Real(dp), Parameter   :: cds(5) = (/ -0.899172, 0.602249, -0.274823, -1.401915, 0.3230064 /)

  ! Other data useful for calculating screening factors
  Real(dp), Allocatable :: zseq(:), zseq53(:), zseqi(:)    ! Sequence of numbers spanning the range of charge numbers in screening
End Module screening_data

Module reac_rate_data
  !-------------------------------------------------------------------------------------------------
  ! This module contains the data necessary to calculate the reaction rates and to map to each
  ! species those reactions which affect it.
  !-------------------------------------------------------------------------------------------------
  Use xnet_types, Only: dp
  Implicit None
  Integer               :: nan(3)                          ! Size of extended reaction->nuclei arrays
  Integer, Allocatable  :: la(:,:), le(:,:)                ! Extents in extended reaction arrays for each reactant
  Integer, Allocatable  :: mu1(:), mu2(:), mu3(:)          ! Index mapping rates to extended arrays

  ! nij(k) is the jth reactant in i-reactant reactions for reaction k in the extended arrays
  Integer, Allocatable  :: n11(:), n21(:), n22(:), n31(:), n32(:), n33(:)

  ! Factors to indicate creation/destruction of species and avoid double-counting for identical reactants
  Real(dp), Allocatable :: a1(:), a2(:), a3(:)

  ! Reaction rates after folding cross-sections in with counting factors
  Real(dp), Allocatable :: b1(:,:), b2(:,:), b3(:,:)       ! Coefficiencts of the Y terms in Eq. 10 of Hix & Meyer (2006)
  !$omp threadprivate(b1,b2,b3)

End Module reac_rate_data

Subroutine read_nuclear_data(data_dir,data_desc)
  !-------------------------------------------------------------------------------------------------
  ! This routine reads, from the file netwinv, the nuclei included along with the nuclear data which
  ! will be needed for later calculations. This data includes the atomic number, the number of
  ! protons and neutrons, and the binding energy (calculated from the tabulated mass excess). Also
  ! the tabulations of the partition functions, g, are read in for later interpolation. Once the set
  ! of nuclear data is read in, it is assigned to the proper nuclei.
  !-------------------------------------------------------------------------------------------------
  Use controls, Only: iheat, nzbatchmx, getNewUnit
  Use nuc_number, Only: ny
  Use nuclear_data
  Use part_funct_data, Only: ng, t9i, g, gg, angm, dlngdt9
  Use screening_data, Only: thbim1
  Use xnet_constants, Only: avn, m_e, m_n, m_p, m_u, five3rd
  Use xnet_types, Only: dp
  Implicit None

  ! Input variables
  Character(*), Intent(in) :: data_dir ! Network data directory

  ! Output variables
  Character(80), Intent(out) :: data_desc ! Brief network description

  ! Local variables
  Real(dp) :: spin ! Ground state spin
  Integer :: lun_desc, lun_winv
  Integer :: i, j, ierr
  Integer :: it9i(ng) ! Raw partition function temperature grid from file
  Character(5) :: nam

  ! Read in the data description
  Open(getNewUnit(lun_desc), file=trim(data_dir)//"/net_desc", status='old', iostat=ierr)
  If ( ierr == 0 ) Then
    Read(lun_desc,"(a80)") data_desc
    Close(lun_desc)
  Else
    data_desc = ""
  EndIf

  ! Read the size of the network and partition function data
  Open(getNewUnit(lun_winv), file=trim(data_dir)//"/netwinv", status='old')
  Read(lun_winv,"(i5)") ny

  ! Read in the partition function temperature grid, and fix endpoints
  Read(lun_winv,"(24i3)") (it9i(i),i=1,ng)
  t9i = real(it9i,dp)
  t9i(1:ng-1) = 0.01*t9i(1:ng-1) ; t9i(ng) = 0.1*t9i(ng)

  ! Skip the list of nuclei, use names from each data entry instead
  Do i = 1, ny
    Read(lun_winv,"(a5)") nam
  EndDo

  ! Set size of nuclear data arrays and read in nuclear data and partition function interpolation table.
  ! nname(0), gg(0) and angm(0) are placeholders for non-nuclei.
  Allocate (nname(0:ny))
  Allocate (aa(ny),zz(ny),nn(ny),be(ny),mex(ny),mm(ny),ia(ny),iz(ny),in(ny))
  Allocate (zz2(ny),zz53(ny),zzi(ny))
  Allocate (g(ng,ny),angm(0:ny))
  nname(0) = ' === '
  angm(0) = 0.0
  ineut = 0 ; iprot = 0
  Do i = 1, ny
    Read(lun_winv,*) nname(i), aa(i), iz(i), in(i), spin, mex(i)
    Read(lun_winv,*) (g(j,i), j=1,ng)
    If ( iz(i) == 0 .and. in(i) == 1 ) ineut = i
    If ( iz(i) == 1 .and. in(i) == 0 ) iprot = i
    ia(i) = nint(aa(i))
    zz(i) = real(iz(i),dp)
    nn(i) = real(in(i),dp)
    angm(i) = 2.0*spin + 1.0
  EndDo
  Close(lun_winv)
  inmin = minval(in) ; inmax = maxval(in)
  izmin = minval(iz) ; izmax = maxval(iz)
  zz2 = zz*zz        ; zz53 = zz**five3rd ; zzi = zz**thbim1

  ! For consistency, use neutron and proton mass excess from netwinv to
  ! calculate binding energies. Otherwise, use CODATA recommended 2014 values.
  If ( ineut > 0 ) Then
    mex_n = mex(ineut)
  Else
    mex_n = m_n - m_u
  EndIf
  If ( iprot > 0 ) Then
    mex_p = mex(iprot)
  Else
    mex_p = m_p + m_e - m_u
  EndIf
  be(:) = mex_n*nn(:) + mex_p*zz(:) - mex(:)

  ! Uncomment the commented end of the line below to use the actual mass instead of A*m_u
  mm(:) = aa(:) / avn! + mex(:)*epmev/(clt*clt)
! mm(:) = zz(:)*(m_p+m_e) + nn(:)*m_n - be(:)*epmev/(clt*clt)

  ! Allocate threadprivate arrays
  !$omp parallel
  Allocate (gg(0:ny,nzbatchmx))
  If(iheat>0) Allocate (dlngdt9(0:ny,nzbatchmx))
  !$omp end parallel

  Return
End Subroutine read_nuclear_data

Subroutine read_reaction_data(data_dir)
  !-------------------------------------------------------------------------------------------------
  ! This routine reads in the necessary reaction data.
  !-------------------------------------------------------------------------------------------------
  Use controls, Only: iheat, iscrn, lun_stderr, nzbatchmx, getNewUnit
  Use cross_sect_data
  Use nuc_number, Only: ny
  Use nuclear_data, Only: izmax, nname, zz
  Use reac_rate_data
  Use screening_data
  Use xnet_constants, Only: five3rd
  Use xnet_interface, Only: read_ffn_data, xnet_terminate
  Use xnet_types, Only: dp
  Implicit None

  ! Input variables
  Character(*), Intent(in) :: data_dir

  ! Local variables
  Integer :: i, j, n, l
  Integer :: nr1, nr2, nr3
  Integer :: lun_s3, lun_s4

  ! Read in nuclear set and numbers of reactions
  Open(getNewUnit(lun_s4), file=trim(data_dir)//"/nets4", form='unformatted', status='old')
  Read(lun_s4) ny
  Read(lun_s4) (nname(i), i=1,ny)
  Read(lun_s4) nffn, nnnu
  Read(lun_s4) (nreac(i), i=1,3)

  ! Allocate reaction arrays
  nr1 = nreac(1) ; nr2 = nreac(2) ; nr3 = nreac(3)
  Allocate (iffn(nr1),innu(nr1))
  Allocate (n1i(4,nr1),iwk1(nr1),ires1(nr1),irev1(nr1),rc1(7,nr1),q1(nr1))
  Allocate (n2i(5,nr2),iwk2(nr2),ires2(nr2),irev2(nr2),rc2(7,nr2),q2(nr2))
  Allocate (n3i(6,nr3),iwk3(nr3),ires3(nr3),irev3(nr3),rc3(7,nr3),q3(nr3))

  ! If there are FFN rates, read in the FFN data and set FFN array sizes
  If ( nffn > 0 ) Call read_ffn_data(nffn,data_dir)                                             !FFN

  ! If there are NNU rates, read in the NNU data and set NNU array sizes
! If ( nnnu > 0 ) Call read_nnu_data(nnnu,data_dir)                                             !NNU

  ! Read in the reaction cross section data
  Open(getNewUnit(lun_s3), file=trim(data_dir)//"/nets3", form='unformatted', status='old')

  ! Read in reaction arrays for 1 reactant reactions
  Do j = 1, nr1
    Read(lun_s3) n, (n1i(l,j), l=1,4), iwk1(j), ires1(j), irev1(j), (rc1(l,j), l=1,7), q1(j)
    If ( n /= j ) Then
      Write(lun_stderr,*) 'Error in nets3, 1',j,n
      Call xnet_terminate('Error in nets3, 1')
    EndIf
  EndDo

  ! Calculate pointers to non-REACLIB data
  Where ( iwk1 == 2 .or. iwk1 == 3 ) ! FFN reaction
    iffn = nint(rc1(1,:))
    innu = 0
  ElseWhere ( iwk1 == 7 .or. iwk1 == 8 ) ! NNU reaction
    iffn = 0
    innu = nint(rc1(1,:))
  ElseWhere
    iffn = 0 ; innu = 0
  EndWhere

  ! Read in reaction arrays for 2 reactant reactions
  Do j = 1, nr2
    Read(lun_s3) n, (n2i(l,j), l=1,5), iwk2(j), ires2(j), irev2(j), (rc2(l,j), l=1,7), q2(j)
    If ( n /= j ) Then
      Write(lun_stderr,*) 'Error in nets3, 2',j,n
      Call xnet_terminate('Error in nets3, 2')
    EndIf
  EndDo

  ! Read in reaction arrays for 3 reactant reactions
  Do j = 1, nr3
    Read(lun_s3) n, (n3i(l,j), l=1,6), iwk3(j), ires3(j), irev3(j), (rc3(l,j), l=1,7), q3(j)
    If ( n /= j ) Then
      Write(lun_stderr,*) 'Error in nets3, 3',j,n
      Call xnet_terminate('Error in nets3, 3')
    EndIf
  EndDo

  ! Allocate and read extents in extended reaction arrays for each reactant
  Allocate (la(3,ny),le(3,ny))
  Do i = 1, ny
    Read(lun_s4) n, (la(j,i), le(j,i), j=1,3)
    If ( n /= i ) Then
      Write(lun_stderr,*) 'Error in nets4',i,n
      Call xnet_terminate('Error in nets4')
    EndIf
  EndDo
  Close(lun_s4)

  ! Allocate and read extended reaction->nuclei arrays linking nuclei to the reactions which affect them
  nan(:) = le(:,ny)
  Allocate (mu1(nan(1)),a1(nan(1)),n11(nan(1)))
  Allocate (mu2(nan(2)),a2(nan(2)),n21(nan(2)),n22(nan(2)))
  Allocate (mu3(nan(3)),a3(nan(3)),n31(nan(3)),n32(nan(3)),n33(nan(3)))
  Do j = 1, nan(1)
    Read(lun_s3) a1(j), mu1(j)
  EndDo
  n11 = n1i(1,mu1)
  Do j = 1, nan(2)
    Read(lun_s3) a2(j), mu2(j)
  EndDo
  n21 = n2i(1,mu2) ; n22 = n2i(2,mu2)
  Do j = 1, nan(3)
    Read(lun_s3) a3(j), mu3(j)
  EndDo
  n31 = n3i(1,mu3) ; n32 = n3i(2,mu3) ; n33 = n3i(3,mu3)
  Close(lun_s3)

  ! Allocate and initialize screening arrays
  If ( iscrn > 0 ) Then
    Allocate (zseq(0:izmax+2),zseq53(0:izmax+2),zseqi(0:izmax+2))
    zseq = (/ (real(i,dp), i=0,izmax+2) /)
    zseq53 = zseq**five3rd
    zseqi = zseq**bip1

    ! 2-reactant screening terms
    Allocate (z21(nr2),z22(nr2),iz21(nr2),iz22(nr2),iz2c(nr2),z2c(nr2),zeta2w(nr2),zeta2i(nr2))
    z21 = zz(n2i(1,:)) ; z22 = zz(n2i(2,:))
    iz21 = nint(z21)   ; iz22 = nint(z22)
    iz2c = iz21 + iz22
    z2c = real(iz2c,dp)
    zeta2w = z21*z22
    zeta2i = z2c**bip1 - z21**bip1 - z22**bip1

    ! 3-reactant screening terms
    Allocate (z31(nr3),z32(nr3),z33(nr3),iz31(nr3),iz32(nr3),iz33(nr3),iz3c(nr3),z3c(nr3),zeta3w(nr3),zeta3i(nr3))
    z31 = zz(n3i(1,:)) ; z32 = zz(n3i(2,:)) ; z33 = zz(n3i(3,:))
    iz31 = nint(z31)   ; iz32 = nint(z32)   ; iz33 = nint(z33)
    iz3c = iz31 + iz32 + iz33
    z3c = real(iz3c,dp)
    zeta3w = z31*z32 + z31*z33 + z32*z33
    zeta3i = z3c**bip1 - z31**bip1 - z32**bip1 - z33**bip1
  EndIf

  ! Allocate threadprivate arrays
  !$omp parallel default(shared)
  Allocate (csect1(nr1,nzbatchmx),h1(nr1,nzbatchmx),rpf1(nr1))
  Allocate (csect2(nr2,nzbatchmx),h2(nr2,nzbatchmx),rpf2(nr2))
  Allocate (csect3(nr3,nzbatchmx),h3(nr3,nzbatchmx))
  Allocate (b1(nan(1),nzbatchmx))
  Allocate (b2(nan(2),nzbatchmx))
  Allocate (b3(nan(3),nzbatchmx))
  If ( iheat > 0 ) Then
    Allocate (dcsect1dt9(nr1,nzbatchmx),dh1dt9(nr1,nzbatchmx),dlnrpf1dt9(nr1))
    Allocate (dcsect2dt9(nr2,nzbatchmx),dh2dt9(nr2,nzbatchmx),dlnrpf2dt9(nr2))
    Allocate (dcsect3dt9(nr3,nzbatchmx),dh3dt9(nr3,nzbatchmx))
  EndIf
  !$omp end parallel

  Return
End Subroutine read_reaction_data
