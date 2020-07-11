module actual_network

  use network_properties
  use physical_constants, only: ERG_PER_MeV
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  public

  character (len=32), parameter :: network_name = "pynucastro"

  real(rt), parameter :: avo = 6.0221417930e23_rt
  real(rt), parameter :: c_light = 2.99792458e10_rt
  real(rt), parameter :: enuc_conv2 = -avo*c_light*c_light

  real(rt), parameter :: ev2erg  = 1.60217648740e-12_rt
  real(rt), parameter :: mev2erg = ev2erg * 1.0e6_rt
  real(rt), parameter :: mev2gr  = mev2erg / c_light**2

  real(rt), parameter :: mass_neutron  = 1.67492721184e-24_rt
  real(rt), parameter :: mass_proton   = 1.67262163783e-24_rt
  real(rt), parameter :: mass_electron = 9.10938215450e-28_rt

  integer, parameter :: nrates = 18


  ! For each rate, we need: rate, drate/dT, screening, dscreening/dT
  integer, parameter :: num_rate_groups = 4

  ! Number of reaclib rates
  integer, parameter :: nrat_reaclib = 14
  integer, parameter :: number_reaclib_sets = 34

  ! Number of tabular rates
  integer, parameter :: nrat_tabular = 4

  ! Binding Energies Per Nucleon (MeV)
  real(rt) :: ebind_per_nucleon(nspec)

  ! bion: Binding Energies (ergs)

  ! Nuclides
  integer, parameter :: jp   = 1
  integer, parameter :: jhe4   = 2
  integer, parameter :: jo16   = 3
  integer, parameter :: jo20   = 4
  integer, parameter :: jf20   = 5
  integer, parameter :: jne20   = 6
  integer, parameter :: jmg24   = 7
  integer, parameter :: jal27   = 8
  integer, parameter :: jsi28   = 9
  integer, parameter :: jp31   = 10
  integer, parameter :: js32   = 11

  ! Reactions
  integer, parameter :: k_ne20__he4_o16   = 1
  integer, parameter :: k_he4_o16__ne20   = 2
  integer, parameter :: k_he4_ne20__mg24   = 3
  integer, parameter :: k_he4_mg24__si28   = 4
  integer, parameter :: k_p_al27__si28   = 5
  integer, parameter :: k_he4_al27__p31   = 6
  integer, parameter :: k_he4_si28__s32   = 7
  integer, parameter :: k_p_p31__s32   = 8
  integer, parameter :: k_o16_o16__p_p31   = 9
  integer, parameter :: k_o16_o16__he4_si28   = 10
  integer, parameter :: k_he4_mg24__p_al27   = 11
  integer, parameter :: k_p_al27__he4_mg24   = 12
  integer, parameter :: k_he4_si28__p_p31   = 13
  integer, parameter :: k_p_p31__he4_si28   = 14
  integer, parameter :: k_f20__o20   = 15
  integer, parameter :: k_ne20__f20   = 16
  integer, parameter :: k_o20__f20   = 17
  integer, parameter :: k_f20__ne20   = 18

  real(rt), allocatable, save :: bion(:), mion(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: bion, mion
#endif

  !$acc declare create(bion, mion)

#ifdef REACT_SPARSE_JACOBIAN
  ! Shape of Jacobian in Compressed Sparse Row format
  integer, parameter   :: NETWORK_SPARSE_JAC_NNZ = 90
  integer, allocatable :: csr_jac_col_index(:), csr_jac_row_count(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: csr_jac_col_index, csr_jac_row_count
#endif
#endif

contains

  subroutine actual_network_init()

    implicit none

    integer :: i

    call network_properties_init()

    ! Allocate ion info arrays
    allocate(bion(nspec))
    allocate(mion(nspec))

    ebind_per_nucleon(jp)   = 0.00000000000000e+00_rt
    ebind_per_nucleon(jhe4)   = 7.07391500000000e+00_rt
    ebind_per_nucleon(jo16)   = 7.97620600000000e+00_rt
    ebind_per_nucleon(jo20)   = 7.56857000000000e+00_rt
    ebind_per_nucleon(jf20)   = 7.72013400000000e+00_rt
    ebind_per_nucleon(jne20)   = 8.03224000000000e+00_rt
    ebind_per_nucleon(jmg24)   = 8.26070900000000e+00_rt
    ebind_per_nucleon(jal27)   = 8.33155300000000e+00_rt
    ebind_per_nucleon(jsi28)   = 8.44774400000000e+00_rt
    ebind_per_nucleon(jp31)   = 8.48116700000000e+00_rt
    ebind_per_nucleon(js32)   = 8.49312900000000e+00_rt

    do i = 1, nspec
       bion(i) = ebind_per_nucleon(i) * aion(i) * ERG_PER_MeV
    end do

    ! Set the mass
    mion(:) = nion(:) * mass_neutron + zion(:) * (mass_proton + mass_electron) &
         - bion(:)/(c_light**2)


    !$acc update device(bion, mion)

#ifdef REACT_SPARSE_JACOBIAN
    ! Set CSR format metadata for Jacobian
    allocate(csr_jac_col_index(NETWORK_SPARSE_JAC_NNZ))
    allocate(csr_jac_row_count(nspec_evolve + 3)) ! neq + 1

    csr_jac_col_index = [ &
      1, &
      2, &
      3, &
      7, &
      8, &
      9, &
      10, &
      12, &
      1, &
      2, &
      3, &
      6, &
      7, &
      8, &
      9, &
      10, &
      12, &
      2, &
      3, &
      6, &
      12, &
      4, &
      5, &
      12, &
      4, &
      5, &
      6, &
      12, &
      2, &
      3, &
      5, &
      6, &
      12, &
      1, &
      2, &
      6, &
      7, &
      8, &
      12, &
      1, &
      2, &
      7, &
      8, &
      12, &
      1, &
      2, &
      3, &
      7, &
      8, &
      9, &
      10, &
      12, &
      1, &
      2, &
      3, &
      8, &
      9, &
      10, &
      12, &
      1, &
      2, &
      9, &
      10, &
      11, &
      12, &
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
      13  ]

    csr_jac_row_count = [ &
      1, &
      9, &
      18, &
      22, &
      25, &
      29, &
      34, &
      40, &
      45, &
      53, &
      60, &
      66, &
      78, &
      91  ]
#endif

  end subroutine actual_network_init


  subroutine actual_network_finalize()
    ! Deallocate storage arrays

    if (allocated(bion)) then
       deallocate(bion)
    endif

    if (allocated(mion)) then
       deallocate(mion)
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
