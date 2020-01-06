module actual_network

  use physical_constants, only: ERG_PER_MeV
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  public

  character (len=32), parameter :: network_name = "pynucastro"

  real(rt), parameter :: avo = 6.0221417930d23
  real(rt), parameter :: c_light = 2.99792458d10
  real(rt), parameter :: enuc_conv2 = -avo*c_light*c_light

  real(rt), parameter :: ev2erg  = 1.60217648740d-12
  real(rt), parameter :: mev2erg = ev2erg*1.0d6
  real(rt), parameter :: mev2gr  = mev2erg/c_light**2

  real(rt), parameter :: mass_neutron  = 1.67492721184d-24
  real(rt), parameter :: mass_proton   = 1.67262163783d-24
  real(rt), parameter :: mass_electron = 9.10938215450d-28

  integer, parameter :: nrates = 19
  integer, parameter :: num_rate_groups = 4

  ! Evolution and auxiliary
  integer, parameter :: nspec_evolve = 13
  integer, parameter :: naux  = 0

  ! Number of nuclear species in the network
  integer, parameter :: nspec = 13

  ! Number of reaclib rates
  integer, parameter :: nrat_reaclib = 19
  integer, parameter :: number_reaclib_sets = 48

  ! Number of tabular rates
  integer, parameter :: nrat_tabular = 0

  ! Binding Energies Per Nucleon (MeV)
  real(rt) :: ebind_per_nucleon(nspec)

  ! aion: Nucleon mass number A
  ! zion: Nucleon atomic number Z
  ! nion: Nucleon neutron number N
  ! bion: Binding Energies (ergs)

  ! Nuclides
  integer, parameter :: jp   = 1
  integer, parameter :: jhe4   = 2
  integer, parameter :: jc12   = 3
  integer, parameter :: jc13   = 4
  integer, parameter :: jn13   = 5
  integer, parameter :: jn14   = 6
  integer, parameter :: jn15   = 7
  integer, parameter :: jo14   = 8
  integer, parameter :: jo15   = 9
  integer, parameter :: jo16   = 10
  integer, parameter :: jo17   = 11
  integer, parameter :: jf17   = 12
  integer, parameter :: jf18   = 13

  ! Reactions
  integer, parameter :: k_n13__c13__weak__wc12   = 1
  integer, parameter :: k_o14__n14__weak__wc12   = 2
  integer, parameter :: k_o15__n15__weak__wc12   = 3
  integer, parameter :: k_f17__o17__weak__wc12   = 4
  integer, parameter :: k_p_c12__n13   = 5
  integer, parameter :: k_he4_c12__o16   = 6
  integer, parameter :: k_p_c13__n14   = 7
  integer, parameter :: k_p_n13__o14   = 8
  integer, parameter :: k_p_n14__o15   = 9
  integer, parameter :: k_he4_n14__f18   = 10
  integer, parameter :: k_p_n15__o16   = 11
  integer, parameter :: k_p_o16__f17   = 12
  integer, parameter :: k_p_o17__f18   = 13
  integer, parameter :: k_he4_n13__p_o16   = 14
  integer, parameter :: k_p_n15__he4_c12   = 15
  integer, parameter :: k_he4_o14__p_f17   = 16
  integer, parameter :: k_p_o17__he4_n14   = 17
  integer, parameter :: k_p_f18__he4_o15   = 18
  integer, parameter :: k_he4_he4_he4__c12   = 19

  ! reactvec indices
  integer, parameter :: i_rate        = 1
  integer, parameter :: i_drate_dt    = 2
  integer, parameter :: i_scor        = 3
  integer, parameter :: i_dscor_dt    = 4
  integer, parameter :: i_dqweak      = 5
  integer, parameter :: i_epart       = 6

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  real(rt), allocatable, save :: aion(:), zion(:), bion(:)
  real(rt), allocatable, save :: nion(:), mion(:), wion(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: aion, zion, bion, nion, mion, wion
#endif

  !$acc declare create(aion, zion, bion, nion, mion, wion)

#ifdef REACT_SPARSE_JACOBIAN
  ! Shape of Jacobian in Compressed Sparse Row format
  integer, parameter   :: NETWORK_SPARSE_JAC_NNZ = 109
  integer, allocatable :: csr_jac_col_index(:), csr_jac_row_count(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: csr_jac_col_index, csr_jac_row_count
#endif
#endif

contains

  subroutine actual_network_init()

    implicit none

    integer :: i

    ! Allocate ion info arrays
    allocate(aion(nspec))
    allocate(zion(nspec))
    allocate(bion(nspec))
    allocate(nion(nspec))
    allocate(mion(nspec))
    allocate(wion(nspec))

    spec_names(jp)   = "hydrogen-1"
    spec_names(jhe4)   = "helium-4"
    spec_names(jc12)   = "carbon-12"
    spec_names(jc13)   = "carbon-13"
    spec_names(jn13)   = "nitrogen-13"
    spec_names(jn14)   = "nitrogen-14"
    spec_names(jn15)   = "nitrogen-15"
    spec_names(jo14)   = "oxygen-14"
    spec_names(jo15)   = "oxygen-15"
    spec_names(jo16)   = "oxygen-16"
    spec_names(jo17)   = "oxygen-17"
    spec_names(jf17)   = "fluorine-17"
    spec_names(jf18)   = "fluorine-18"

    short_spec_names(jp)   = "h1"
    short_spec_names(jhe4)   = "he4"
    short_spec_names(jc12)   = "c12"
    short_spec_names(jc13)   = "c13"
    short_spec_names(jn13)   = "n13"
    short_spec_names(jn14)   = "n14"
    short_spec_names(jn15)   = "n15"
    short_spec_names(jo14)   = "o14"
    short_spec_names(jo15)   = "o15"
    short_spec_names(jo16)   = "o16"
    short_spec_names(jo17)   = "o17"
    short_spec_names(jf17)   = "f17"
    short_spec_names(jf18)   = "f18"

    ebind_per_nucleon(jp)   = 0.00000000000000d+00
    ebind_per_nucleon(jhe4)   = 7.07391500000000d+00
    ebind_per_nucleon(jc12)   = 7.68014400000000d+00
    ebind_per_nucleon(jc13)   = 7.46984900000000d+00
    ebind_per_nucleon(jn13)   = 7.23886300000000d+00
    ebind_per_nucleon(jn14)   = 7.47561400000000d+00
    ebind_per_nucleon(jn15)   = 7.69946000000000d+00
    ebind_per_nucleon(jo14)   = 7.05227800000000d+00
    ebind_per_nucleon(jo15)   = 7.46369200000000d+00
    ebind_per_nucleon(jo16)   = 7.97620600000000d+00
    ebind_per_nucleon(jo17)   = 7.75072800000000d+00
    ebind_per_nucleon(jf17)   = 7.54232800000000d+00
    ebind_per_nucleon(jf18)   = 7.63163800000000d+00

    aion(jp)   = 1.00000000000000d+00
    aion(jhe4)   = 4.00000000000000d+00
    aion(jc12)   = 1.20000000000000d+01
    aion(jc13)   = 1.30000000000000d+01
    aion(jn13)   = 1.30000000000000d+01
    aion(jn14)   = 1.40000000000000d+01
    aion(jn15)   = 1.50000000000000d+01
    aion(jo14)   = 1.40000000000000d+01
    aion(jo15)   = 1.50000000000000d+01
    aion(jo16)   = 1.60000000000000d+01
    aion(jo17)   = 1.70000000000000d+01
    aion(jf17)   = 1.70000000000000d+01
    aion(jf18)   = 1.80000000000000d+01

    zion(jp)   = 1.00000000000000d+00
    zion(jhe4)   = 2.00000000000000d+00
    zion(jc12)   = 6.00000000000000d+00
    zion(jc13)   = 6.00000000000000d+00
    zion(jn13)   = 7.00000000000000d+00
    zion(jn14)   = 7.00000000000000d+00
    zion(jn15)   = 7.00000000000000d+00
    zion(jo14)   = 8.00000000000000d+00
    zion(jo15)   = 8.00000000000000d+00
    zion(jo16)   = 8.00000000000000d+00
    zion(jo17)   = 8.00000000000000d+00
    zion(jf17)   = 9.00000000000000d+00
    zion(jf18)   = 9.00000000000000d+00

    nion(jp)   = 0.00000000000000d+00
    nion(jhe4)   = 2.00000000000000d+00
    nion(jc12)   = 6.00000000000000d+00
    nion(jc13)   = 7.00000000000000d+00
    nion(jn13)   = 6.00000000000000d+00
    nion(jn14)   = 7.00000000000000d+00
    nion(jn15)   = 8.00000000000000d+00
    nion(jo14)   = 6.00000000000000d+00
    nion(jo15)   = 7.00000000000000d+00
    nion(jo16)   = 8.00000000000000d+00
    nion(jo17)   = 9.00000000000000d+00
    nion(jf17)   = 8.00000000000000d+00
    nion(jf18)   = 9.00000000000000d+00

    do i = 1, nspec
       bion(i) = ebind_per_nucleon(i) * aion(i) * ERG_PER_MeV
    end do

    ! Set the mass
    mion(:) = nion(:) * mass_neutron + zion(:) * (mass_proton + mass_electron) &
         - bion(:)/(c_light**2)

    ! Molar mass
    wion(:) = avo * mion(:)

    ! Common approximation
    !wion(:) = aion(:)

    !$acc update device(aion, zion, bion, nion, mion, wion)

#ifdef REACT_SPARSE_JACOBIAN
    ! Set CSR format metadata for Jacobian
    allocate(csr_jac_col_index(NETWORK_SPARSE_JAC_NNZ))
    allocate(csr_jac_row_count(nspec_evolve + 3)) ! neq + 1

    csr_jac_col_index = [ &
      1, &
      2, &
      3, &
      4, &
      5, &
      6, &
      7, &
      8, &
      10, &
      11, &
      13, &
      14, &
      1, &
      2, &
      3, &
      5, &
      6, &
      7, &
      8, &
      11, &
      13, &
      14, &
      1, &
      2, &
      3, &
      7, &
      14, &
      1, &
      4, &
      5, &
      14, &
      1, &
      2, &
      3, &
      5, &
      14, &
      1, &
      2, &
      4, &
      6, &
      8, &
      11, &
      14, &
      1, &
      7, &
      9, &
      14, &
      1, &
      2, &
      5, &
      8, &
      14, &
      1, &
      6, &
      9, &
      13, &
      14, &
      1, &
      2, &
      3, &
      5, &
      7, &
      10, &
      14, &
      1, &
      11, &
      12, &
      14, &
      1, &
      2, &
      8, &
      10, &
      12, &
      14, &
      1, &
      2, &
      6, &
      11, &
      13, &
      14, &
      1, &
      2, &
      3, &
      4, &
      5, &
      6, &
      7, &
      8, &
      9, &
      10, &
      11, &
      12, &
      13, &
      14, &
      1, &
      2, &
      3, &
      4, &
      5, &
      6, &
      7, &
      8, &
      9, &
      10, &
      11, &
      12, &
      13, &
      14, &
      15  ]

    csr_jac_row_count = [ &
      1, &
      13, &
      23, &
      28, &
      32, &
      37, &
      44, &
      48, &
      53, &
      58, &
      65, &
      69, &
      75, &
      81, &
      95, &
      110  ]
#endif

  end subroutine actual_network_init


  subroutine actual_network_finalize()
    ! Deallocate storage arrays

    if (allocated(aion)) then
       deallocate(aion)
    endif

    if (allocated(zion)) then
       deallocate(zion)
    endif

    if (allocated(bion)) then
       deallocate(bion)
    endif

    if (allocated(nion)) then
       deallocate(nion)
    endif

    if (allocated(mion)) then
       deallocate(mion)
    endif

    if (allocated(wion)) then
       deallocate(wion)
    endif

#ifdef REACT_SPARSE_JACOBIAN
    if (allocated(csr_jac_col_index)) then
       deallocate(csr_jac_col_index)
    endif

    if (allocated(csr_jac_row_count)) then
       deallocate(csr_jac_row_count)
    endif
#endif

  end subroutine actual_network_finalize


  subroutine ener_gener_rate(dydt, enuc)
    ! Computes the instantaneous energy generation rate

    !$acc routine seq

    implicit none

    real(rt) :: dydt(nspec), enuc

    !$gpu

    ! This is basically e = m c**2

    enuc = sum(dydt(:) * mion(:)) * enuc_conv2

  end subroutine ener_gener_rate

end module actual_network
