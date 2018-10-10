module actual_network

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  double precision, parameter, private :: clight = 2.99792458d10
  double precision, parameter, private :: ev2erg  = 1.60217648740d-12
  double precision, parameter, private :: mev2erg = ev2erg*1.0d6
  double precision, parameter, private :: mev2gr  = mev2erg/clight**2

  double precision, parameter, private :: mn = 1.67492721184d-24
  double precision, parameter, private :: mp = 1.67262163783d-24
  double precision, parameter, private :: me = 9.1093821545d-28

  integer, parameter :: nspec = 3
  integer, parameter :: nspec_evolve = 1
  integer, parameter :: naux  = 0

  integer, parameter :: ic12  = 1
  integer, parameter :: io16  = 2
  integer, parameter :: img24 = 3

  character (len=16), save :: spec_names(nspec) 
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  double precision, allocatable :: aion(:), zion(:), nion(:)
  double precision, allocatable :: bion(:), mion(:), wion(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: aion, zion, nion, bion, mion, wion
#endif

  !$acc declare create(aion, zion, bion, nion, mion, wion)

  integer, parameter :: nrates = 1
  integer, parameter :: num_rate_groups = 4

  ! Conversion factor for the nuclear energy generation rate.

  double precision, parameter :: avo = 6.0221417930d23
  double precision, parameter :: enuc_conv2 = -avo*clight*clight

  ! Shape of Jacobian in Compressed Sparse Row format
  integer, parameter :: NETWORK_CSR_JAC_NNZ = 7
  integer, allocatable :: csr_jac_col_index(:), csr_jac_row_count(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: csr_jac_col_index, csr_jac_row_count
#endif

contains

  subroutine actual_network_init

    implicit none

    spec_names(ic12)  = "carbon-12"
    spec_names(io16)  = "oxygen-16"
    spec_names(img24) = "magnesium-24"

    short_spec_names(ic12)  = "C12"
    short_spec_names(io16)  = "O16"
    short_spec_names(img24) = "Mg24"

    allocate(aion(nspec))
    allocate(zion(nspec))
    allocate(nion(nspec))
    allocate(bion(nspec))
    allocate(mion(nspec))
    allocate(wion(nspec))
    
    aion(ic12)  = 12.0d0
    aion(io16)  = 16.0d0
    aion(img24) = 24.0d0

    zion(ic12)  = 6.0d0
    zion(io16)  = 8.0d0
    zion(img24) = 12.0d0

    ! Binding energies per nucleus in MeV
    bion(ic12)  = 92.16294d0
    bion(io16)  = 127.62093d0
    bion(img24) = 198.2579d0

    ! Set the number of neutrons
    nion(:) = aion(:) - zion(:)

    ! Set the mass
    mion(:) = nion(:) * mn + zion(:) * (mp + me) - bion(:) * mev2gr

    ! Molar mass
    wion(:) = avo * mion(:)

    ! Common approximation
    wion(:) = aion(:)

    ! Set CSR format metadata for Jacobian
    allocate(csr_jac_col_index(NETWORK_CSR_JAC_NNZ))
    allocate(csr_jac_row_count(nspec_evolve + 3)) ! neq + 1

    csr_jac_col_index = [1, 2, 1, 2, 1, 2, 3]
    csr_jac_row_count = [1, 3, 5, 8]

    !$acc update device(aion, zion, bion, nion, mion, wion)

  end subroutine actual_network_init


  subroutine actual_network_finalize

    implicit none

    if (allocated(aion)) then
       deallocate(aion)
    endif
    if (allocated(zion)) then
       deallocate(zion)
    endif
    if (allocated(nion)) then
       deallocate(nion)
    endif
    if (allocated(bion)) then
       deallocate(bion)
    endif
    if (allocated(mion)) then
       deallocate(mion)
    endif
    if (allocated(wion)) then
       deallocate(wion)
    endif

  end subroutine actual_network_finalize

end module actual_network
