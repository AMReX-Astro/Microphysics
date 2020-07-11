module actual_network

  use network_properties
  use amrex_fort_module, only : rt => amrex_real

  implicit none


  real(rt)        , parameter, private :: clight = 2.99792458e10_rt
  real(rt)        , parameter, private :: ev2erg  = 1.60217648740e-12_rt
  real(rt)        , parameter, private :: mev2erg = ev2erg*1.0e6_rt
  real(rt)        , parameter, private :: mev2gr  = mev2erg/clight**2

  character (len=32), parameter :: network_name = "ignition_simple"

  real(rt)        , parameter, private :: mn = 1.67492721184e-24_rt
  real(rt)        , parameter, private :: mp = 1.67262163783e-24_rt
  real(rt)        , parameter, private :: me = 9.1093821545e-28_rt

  integer, parameter :: ic12  = 1
  integer, parameter :: io16  = 2
  integer, parameter :: img24 = 3

  real(rt)        , allocatable :: bion(:), mion(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: bion, mion
#endif

  !$acc declare create(bion, mion)

  integer, parameter :: nrates = 1
  integer, parameter :: num_rate_groups = 4

  ! Conversion factor for the nuclear energy generation rate.

  real(rt)        , parameter :: avo = 6.0221417930e23_rt
  real(rt)        , parameter :: enuc_conv2 = -avo*clight*clight

#ifdef REACT_SPARSE_JACOBIAN
  ! Shape of Jacobian in Compressed Sparse Row format
  integer, parameter :: NETWORK_SPARSE_JAC_NNZ = 7
  integer, allocatable :: csr_jac_col_index(:), csr_jac_row_count(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: csr_jac_col_index, csr_jac_row_count
#endif
#endif

contains

  subroutine actual_network_init

    implicit none

    call network_properties_init()

    allocate(bion(nspec))
    allocate(mion(nspec))

    ! Binding energies per nucleus in MeV
    bion(ic12)  = 92.16294e0_rt
    bion(io16)  = 127.62093e0_rt
    bion(img24) = 198.2579e0_rt

    ! Set the mass
    mion(:) = nion(:) * mn + zion(:) * (mp + me) - bion(:) * mev2gr

#ifdef REACT_SPARSE_JACOBIAN
    ! Set CSR format metadata for Jacobian
    allocate(csr_jac_col_index(NETWORK_SPARSE_JAC_NNZ))
    allocate(csr_jac_row_count(nspec + 3)) ! neq + 1

    csr_jac_col_index = [1, 2, 1, 2, 1, 2, 3]
    csr_jac_row_count = [1, 3, 5, 8]
#endif

    !$acc update device(nion, mion)

  end subroutine actual_network_init


  subroutine actual_network_finalize

    implicit none

    call network_properties_finalize()

    if (allocated(bion)) then
       deallocate(bion)
    endif
    if (allocated(mion)) then
       deallocate(mion)
    endif

  end subroutine actual_network_finalize

end module actual_network
