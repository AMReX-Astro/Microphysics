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

  integer, parameter :: nrates = 96


  ! For each rate, we need: rate, drate/dT, screening, dscreening/dT
  integer, parameter :: num_rate_groups = 4

  ! Number of reaclib rates
  integer, parameter :: nrat_reaclib = 96
  integer, parameter :: number_reaclib_sets = 158

  ! Number of tabular rates
  integer, parameter :: nrat_tabular = 0

  ! Binding Energies Per Nucleon (MeV)
  real(rt) :: ebind_per_nucleon(nspec)

  ! bion: Binding Energies (ergs)

  ! Nuclides
  integer, parameter :: jp   = 1
  integer, parameter :: jhe4   = 2
  integer, parameter :: jc12   = 3
  integer, parameter :: jc14   = 4
  integer, parameter :: jn13   = 5
  integer, parameter :: jn14   = 6
  integer, parameter :: jo16   = 7
  integer, parameter :: jo18   = 8
  integer, parameter :: jf18   = 9
  integer, parameter :: jne20   = 10
  integer, parameter :: jne21   = 11
  integer, parameter :: jmg24   = 12
  integer, parameter :: jal27   = 13
  integer, parameter :: jsi28   = 14
  integer, parameter :: jp31   = 15
  integer, parameter :: js32   = 16
  integer, parameter :: jcl35   = 17
  integer, parameter :: jar36   = 18
  integer, parameter :: jk39   = 19
  integer, parameter :: jca40   = 20
  integer, parameter :: jsc43   = 21
  integer, parameter :: jti44   = 22
  integer, parameter :: jv47   = 23
  integer, parameter :: jcr48   = 24
  integer, parameter :: jmn51   = 25
  integer, parameter :: jfe52   = 26
  integer, parameter :: jco55   = 27
  integer, parameter :: jni56   = 28

  ! Reactions
  integer, parameter :: k_c14__n14__weak__wc12   = 1
  integer, parameter :: k_f18__o18__weak__wc12   = 2
  integer, parameter :: k_n13__p_c12   = 3
  integer, parameter :: k_o16__he4_c12   = 4
  integer, parameter :: k_o18__he4_c14   = 5
  integer, parameter :: k_f18__he4_n14   = 6
  integer, parameter :: k_ne20__he4_o16   = 7
  integer, parameter :: k_mg24__he4_ne20   = 8
  integer, parameter :: k_si28__p_al27   = 9
  integer, parameter :: k_si28__he4_mg24   = 10
  integer, parameter :: k_p31__he4_al27   = 11
  integer, parameter :: k_s32__p_p31   = 12
  integer, parameter :: k_s32__he4_si28   = 13
  integer, parameter :: k_cl35__he4_p31   = 14
  integer, parameter :: k_ar36__p_cl35   = 15
  integer, parameter :: k_ar36__he4_s32   = 16
  integer, parameter :: k_k39__he4_cl35   = 17
  integer, parameter :: k_ca40__p_k39   = 18
  integer, parameter :: k_ca40__he4_ar36   = 19
  integer, parameter :: k_sc43__he4_k39   = 20
  integer, parameter :: k_ti44__p_sc43   = 21
  integer, parameter :: k_ti44__he4_ca40   = 22
  integer, parameter :: k_v47__he4_sc43   = 23
  integer, parameter :: k_cr48__p_v47   = 24
  integer, parameter :: k_cr48__he4_ti44   = 25
  integer, parameter :: k_mn51__he4_v47   = 26
  integer, parameter :: k_fe52__p_mn51   = 27
  integer, parameter :: k_fe52__he4_cr48   = 28
  integer, parameter :: k_co55__he4_mn51   = 29
  integer, parameter :: k_ni56__p_co55   = 30
  integer, parameter :: k_ni56__he4_fe52   = 31
  integer, parameter :: k_c12__he4_he4_he4   = 32
  integer, parameter :: k_p_c12__n13   = 33
  integer, parameter :: k_he4_c12__o16   = 34
  integer, parameter :: k_he4_c14__o18   = 35
  integer, parameter :: k_he4_n14__f18   = 36
  integer, parameter :: k_he4_o16__ne20   = 37
  integer, parameter :: k_he4_ne20__mg24   = 38
  integer, parameter :: k_he4_mg24__si28   = 39
  integer, parameter :: k_p_al27__si28   = 40
  integer, parameter :: k_he4_al27__p31   = 41
  integer, parameter :: k_he4_si28__s32   = 42
  integer, parameter :: k_p_p31__s32   = 43
  integer, parameter :: k_he4_p31__cl35   = 44
  integer, parameter :: k_he4_s32__ar36   = 45
  integer, parameter :: k_p_cl35__ar36   = 46
  integer, parameter :: k_he4_cl35__k39   = 47
  integer, parameter :: k_he4_ar36__ca40   = 48
  integer, parameter :: k_p_k39__ca40   = 49
  integer, parameter :: k_he4_k39__sc43   = 50
  integer, parameter :: k_he4_ca40__ti44   = 51
  integer, parameter :: k_p_sc43__ti44   = 52
  integer, parameter :: k_he4_sc43__v47   = 53
  integer, parameter :: k_he4_ti44__cr48   = 54
  integer, parameter :: k_p_v47__cr48   = 55
  integer, parameter :: k_he4_v47__mn51   = 56
  integer, parameter :: k_he4_cr48__fe52   = 57
  integer, parameter :: k_p_mn51__fe52   = 58
  integer, parameter :: k_he4_mn51__co55   = 59
  integer, parameter :: k_he4_fe52__ni56   = 60
  integer, parameter :: k_p_co55__ni56   = 61
  integer, parameter :: k_c12_c12__he4_ne20   = 62
  integer, parameter :: k_he4_n13__p_o16   = 63
  integer, parameter :: k_p_o16__he4_n13   = 64
  integer, parameter :: k_c12_o16__p_al27   = 65
  integer, parameter :: k_c12_o16__he4_mg24   = 66
  integer, parameter :: k_o16_o16__p_p31   = 67
  integer, parameter :: k_o16_o16__he4_si28   = 68
  integer, parameter :: k_he4_f18__p_ne21   = 69
  integer, parameter :: k_he4_ne20__c12_c12   = 70
  integer, parameter :: k_c12_ne20__p_p31   = 71
  integer, parameter :: k_c12_ne20__he4_si28   = 72
  integer, parameter :: k_p_ne21__he4_f18   = 73
  integer, parameter :: k_he4_mg24__p_al27   = 74
  integer, parameter :: k_he4_mg24__c12_o16   = 75
  integer, parameter :: k_p_al27__he4_mg24   = 76
  integer, parameter :: k_p_al27__c12_o16   = 77
  integer, parameter :: k_he4_si28__p_p31   = 78
  integer, parameter :: k_he4_si28__c12_ne20   = 79
  integer, parameter :: k_he4_si28__o16_o16   = 80
  integer, parameter :: k_p_p31__he4_si28   = 81
  integer, parameter :: k_p_p31__c12_ne20   = 82
  integer, parameter :: k_p_p31__o16_o16   = 83
  integer, parameter :: k_he4_s32__p_cl35   = 84
  integer, parameter :: k_p_cl35__he4_s32   = 85
  integer, parameter :: k_he4_ar36__p_k39   = 86
  integer, parameter :: k_p_k39__he4_ar36   = 87
  integer, parameter :: k_he4_ca40__p_sc43   = 88
  integer, parameter :: k_p_sc43__he4_ca40   = 89
  integer, parameter :: k_he4_ti44__p_v47   = 90
  integer, parameter :: k_p_v47__he4_ti44   = 91
  integer, parameter :: k_he4_cr48__p_mn51   = 92
  integer, parameter :: k_p_mn51__he4_cr48   = 93
  integer, parameter :: k_he4_fe52__p_co55   = 94
  integer, parameter :: k_p_co55__he4_fe52   = 95
  integer, parameter :: k_he4_he4_he4__c12   = 96

  real(rt), allocatable, save :: bion(:), mion(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: bion, mion
#endif

  !$acc declare create(bion, mion)

#ifdef REACT_SPARSE_JACOBIAN
  ! Shape of Jacobian in Compressed Sparse Row format
  integer, parameter   :: NETWORK_SPARSE_JAC_NNZ = 317
  integer, allocatable :: csr_jac_col_index(:), csr_jac_row_count(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: csr_jac_col_index, csr_jac_row_count
#endif
#endif

contains

  subroutine actual_network_init()

    implicit none

    integer :: i

    ! initialize the autogenerated network properties
    call network_properties_init()

    ! Allocate ion info arrays
    allocate(bion(nspec))
    allocate(mion(nspec))

    ebind_per_nucleon(jp)   = 0.00000000000000e+00_rt
    ebind_per_nucleon(jhe4)   = 7.07391500000000e+00_rt
    ebind_per_nucleon(jc12)   = 7.68014400000000e+00_rt
    ebind_per_nucleon(jc14)   = 7.52031900000000e+00_rt
    ebind_per_nucleon(jn13)   = 7.23886300000000e+00_rt
    ebind_per_nucleon(jn14)   = 7.47561400000000e+00_rt
    ebind_per_nucleon(jo16)   = 7.97620600000000e+00_rt
    ebind_per_nucleon(jo18)   = 7.76709700000000e+00_rt
    ebind_per_nucleon(jf18)   = 7.63163800000000e+00_rt
    ebind_per_nucleon(jne20)   = 8.03224000000000e+00_rt
    ebind_per_nucleon(jne21)   = 7.97171300000000e+00_rt
    ebind_per_nucleon(jmg24)   = 8.26070900000000e+00_rt
    ebind_per_nucleon(jal27)   = 8.33155300000000e+00_rt
    ebind_per_nucleon(jsi28)   = 8.44774400000000e+00_rt
    ebind_per_nucleon(jp31)   = 8.48116700000000e+00_rt
    ebind_per_nucleon(js32)   = 8.49312900000000e+00_rt
    ebind_per_nucleon(jcl35)   = 8.52027800000000e+00_rt
    ebind_per_nucleon(jar36)   = 8.51990900000000e+00_rt
    ebind_per_nucleon(jk39)   = 8.55702500000000e+00_rt
    ebind_per_nucleon(jca40)   = 8.55130300000000e+00_rt
    ebind_per_nucleon(jsc43)   = 8.53082500000000e+00_rt
    ebind_per_nucleon(jti44)   = 8.53352000000000e+00_rt
    ebind_per_nucleon(jv47)   = 8.58222500000000e+00_rt
    ebind_per_nucleon(jcr48)   = 8.57226900000000e+00_rt
    ebind_per_nucleon(jmn51)   = 8.63377200000000e+00_rt
    ebind_per_nucleon(jfe52)   = 8.60957400000000e+00_rt
    ebind_per_nucleon(jco55)   = 8.66961800000000e+00_rt
    ebind_per_nucleon(jni56)   = 8.64277900000000e+00_rt

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
      5, &
      7, &
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
      19, &
      20, &
      21, &
      22, &
      23, &
      24, &
      25, &
      26, &
      27, &
      28, &
      29, &
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
      19, &
      20, &
      21, &
      22, &
      23, &
      24, &
      25, &
      26, &
      27, &
      28, &
      29, &
      1, &
      2, &
      3, &
      5, &
      7, &
      10, &
      12, &
      13, &
      14, &
      15, &
      29, &
      2, &
      4, &
      8, &
      29, &
      1, &
      2, &
      3, &
      5, &
      7, &
      29, &
      2, &
      4, &
      6, &
      9, &
      29, &
      1, &
      2, &
      3, &
      5, &
      7, &
      10, &
      12, &
      13, &
      14, &
      15, &
      29, &
      2, &
      4, &
      8, &
      9, &
      29, &
      1, &
      2, &
      6, &
      9, &
      11, &
      29, &
      1, &
      2, &
      3, &
      7, &
      10, &
      12, &
      14, &
      15, &
      29, &
      1, &
      2, &
      9, &
      11, &
      29, &
      1, &
      2, &
      3, &
      7, &
      10, &
      12, &
      13, &
      14, &
      29, &
      1, &
      2, &
      3, &
      7, &
      12, &
      13, &
      14, &
      15, &
      29, &
      1, &
      2, &
      3, &
      7, &
      10, &
      12, &
      13, &
      14, &
      15, &
      16, &
      29, &
      1, &
      2, &
      3, &
      7, &
      10, &
      13, &
      14, &
      15, &
      16, &
      17, &
      29, &
      1, &
      2, &
      14, &
      15, &
      16, &
      17, &
      18, &
      29, &
      1, &
      2, &
      15, &
      16, &
      17, &
      18, &
      19, &
      29, &
      1, &
      2, &
      16, &
      17, &
      18, &
      19, &
      20, &
      29, &
      1, &
      2, &
      17, &
      18, &
      19, &
      20, &
      21, &
      29, &
      1, &
      2, &
      18, &
      19, &
      20, &
      21, &
      22, &
      29, &
      1, &
      2, &
      19, &
      20, &
      21, &
      22, &
      23, &
      29, &
      1, &
      2, &
      20, &
      21, &
      22, &
      23, &
      24, &
      29, &
      1, &
      2, &
      21, &
      22, &
      23, &
      24, &
      25, &
      29, &
      1, &
      2, &
      22, &
      23, &
      24, &
      25, &
      26, &
      29, &
      1, &
      2, &
      23, &
      24, &
      25, &
      26, &
      27, &
      29, &
      1, &
      2, &
      24, &
      25, &
      26, &
      27, &
      28, &
      29, &
      1, &
      2, &
      25, &
      26, &
      27, &
      28, &
      29, &
      1, &
      2, &
      26, &
      27, &
      28, &
      29, &
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
      19, &
      20, &
      21, &
      22, &
      23, &
      24, &
      25, &
      26, &
      27, &
      28, &
      29, &
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
      19, &
      20, &
      21, &
      22, &
      23, &
      24, &
      25, &
      26, &
      27, &
      28, &
      29, &
      30  ]

    csr_jac_row_count = [ &
      1, &
      27, &
      56, &
      67, &
      71, &
      77, &
      82, &
      93, &
      98, &
      104, &
      113, &
      118, &
      127, &
      136, &
      147, &
      158, &
      166, &
      174, &
      182, &
      190, &
      198, &
      206, &
      214, &
      222, &
      230, &
      238, &
      246, &
      253, &
      259, &
      288, &
      318  ]
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
