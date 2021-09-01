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

  integer, parameter :: nrates = 31


  ! For each rate, we need: rate, drate/dT, screening, dscreening/dT
  integer, parameter :: num_rate_groups = 4

  ! Number of reaclib rates
  integer, parameter :: nrat_reaclib = 31
  integer, parameter :: number_reaclib_sets = 64

  ! Number of tabular rates
  integer, parameter :: nrat_tabular = 0

  ! Binding Energies Per Nucleon (MeV)
  real(rt) :: ebind_per_nucleon(nspec)

  ! bion: Binding Energies (ergs)

  ! Nuclides
  integer, parameter :: jp   = 1
  integer, parameter :: jd   = 2
  integer, parameter :: jhe3   = 3
  integer, parameter :: jhe4   = 4
  integer, parameter :: jbe7   = 5
  integer, parameter :: jb8   = 6
  integer, parameter :: jc12   = 7
  integer, parameter :: jc13   = 8
  integer, parameter :: jn13   = 9
  integer, parameter :: jn14   = 10
  integer, parameter :: jn15   = 11
  integer, parameter :: jo14   = 12
  integer, parameter :: jo15   = 13
  integer, parameter :: jo16   = 14
  integer, parameter :: jo17   = 15
  integer, parameter :: jf17   = 16
  integer, parameter :: jf18   = 17

  ! Reactions
  integer, parameter :: k_n13__c13__weak__wc12   = 1
  integer, parameter :: k_o14__n14__weak__wc12   = 2
  integer, parameter :: k_o15__n15__weak__wc12   = 3
  integer, parameter :: k_f17__o17__weak__wc12   = 4
  integer, parameter :: k_b8__he4_he4__weak__wc12   = 5
  integer, parameter :: k_p_p__d__weak__bet_pos_   = 6
  integer, parameter :: k_p_p__d__weak__electron_capture   = 7
  integer, parameter :: k_p_d__he3   = 8
  integer, parameter :: k_d_d__he4   = 9
  integer, parameter :: k_p_he3__he4__weak__bet_pos_   = 10
  integer, parameter :: k_he4_he3__be7   = 11
  integer, parameter :: k_p_be7__b8   = 12
  integer, parameter :: k_p_c12__n13   = 13
  integer, parameter :: k_he4_c12__o16   = 14
  integer, parameter :: k_p_c13__n14   = 15
  integer, parameter :: k_p_n13__o14   = 16
  integer, parameter :: k_p_n14__o15   = 17
  integer, parameter :: k_he4_n14__f18   = 18
  integer, parameter :: k_p_n15__o16   = 19
  integer, parameter :: k_p_o16__f17   = 20
  integer, parameter :: k_p_o17__f18   = 21
  integer, parameter :: k_d_he3__p_he4   = 22
  integer, parameter :: k_he4_n13__p_o16   = 23
  integer, parameter :: k_p_n15__he4_c12   = 24
  integer, parameter :: k_he4_o14__p_f17   = 25
  integer, parameter :: k_p_o17__he4_n14   = 26
  integer, parameter :: k_p_f18__he4_o15   = 27
  integer, parameter :: k_he3_he3__p_p_he4   = 28
  integer, parameter :: k_d_be7__p_he4_he4   = 29
  integer, parameter :: k_he3_be7__p_p_he4_he4   = 30
  integer, parameter :: k_he4_he4_he4__c12   = 31

  real(rt), allocatable, save :: bion(:), mion(:)

#ifdef REACT_SPARSE_JACOBIAN
  ! Shape of Jacobian in Compressed Sparse Row format
  integer, parameter   :: NETWORK_SPARSE_JAC_NNZ = 145
  integer, allocatable :: csr_jac_col_index(:), csr_jac_row_count(:)
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
    ebind_per_nucleon(jd)   = 1.11228300000000e+00_rt
    ebind_per_nucleon(jhe3)   = 2.57268000000000e+00_rt
    ebind_per_nucleon(jhe4)   = 7.07391500000000e+00_rt
    ebind_per_nucleon(jbe7)   = 5.37154800000000e+00_rt
    ebind_per_nucleon(jb8)   = 4.71715500000000e+00_rt
    ebind_per_nucleon(jc12)   = 7.68014400000000e+00_rt
    ebind_per_nucleon(jc13)   = 7.46984900000000e+00_rt
    ebind_per_nucleon(jn13)   = 7.23886300000000e+00_rt
    ebind_per_nucleon(jn14)   = 7.47561400000000e+00_rt
    ebind_per_nucleon(jn15)   = 7.69946000000000e+00_rt
    ebind_per_nucleon(jo14)   = 7.05227800000000e+00_rt
    ebind_per_nucleon(jo15)   = 7.46369200000000e+00_rt
    ebind_per_nucleon(jo16)   = 7.97620600000000e+00_rt
    ebind_per_nucleon(jo17)   = 7.75072800000000e+00_rt
    ebind_per_nucleon(jf17)   = 7.54232800000000e+00_rt
    ebind_per_nucleon(jf18)   = 7.63163800000000e+00_rt

    do i = 1, nspec
       bion(i) = ebind_per_nucleon(i) * aion(i) * ERG_PER_MeV
    end do

    ! Set the mass
    mion(:) = nion(:) * mass_neutron + zion(:) * (mass_proton + mass_electron) &
         - bion(:)/(c_light**2)


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
      7, &
      8, &
      9, &
      10, &
      11, &
      12, &
      14, &
      15, &
      17, &
      18, &
      1, &
      2, &
      3, &
      5, &
      18, &
      1, &
      2, &
      3, &
      4, &
      5, &
      18, &
      1, &
      2, &
      3, &
      4, &
      5, &
      6, &
      7, &
      9, &
      10, &
      11, &
      12, &
      15, &
      17, &
      18, &
      1, &
      2, &
      3, &
      4, &
      5, &
      18, &
      1, &
      5, &
      6, &
      18, &
      1, &
      4, &
      7, &
      11, &
      18, &
      1, &
      8, &
      9, &
      18, &
      1, &
      4, &
      7, &
      9, &
      18, &
      1, &
      4, &
      8, &
      10, &
      12, &
      15, &
      18, &
      1, &
      11, &
      13, &
      18, &
      1, &
      4, &
      9, &
      12, &
      18, &
      1, &
      10, &
      13, &
      17, &
      18, &
      1, &
      4, &
      7, &
      9, &
      11, &
      14, &
      18, &
      1, &
      15, &
      16, &
      18, &
      1, &
      4, &
      12, &
      14, &
      16, &
      18, &
      1, &
      4, &
      10, &
      15, &
      17, &
      18, &
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
      15, &
      16, &
      17, &
      18, &
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
      15, &
      16, &
      17, &
      18, &
      19  ]

    csr_jac_row_count = [ &
      1, &
      16, &
      21, &
      27, &
      41, &
      47, &
      51, &
      56, &
      60, &
      65, &
      72, &
      76, &
      81, &
      86, &
      93, &
      97, &
      103, &
      109, &
      127, &
      146  ]
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

    implicit none

    real(rt) :: dydt(nspec), enuc

    ! This is basically e = m c**2

    enuc = sum(dydt(:) * mion(:)) * enuc_conv2

  end subroutine ener_gener_rate

end module actual_network
