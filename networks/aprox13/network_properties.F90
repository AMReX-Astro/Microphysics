! An automatically file of network properties.  This provides the properties
! of a set of non-reacting species.
!
! nspec            -- the number of species
! naux             -- the number of auxiliary variables
!
! aion             -- atomic number
! zion             -- proton number
!
! spec_names       -- the name of the isotope
! short_spec_names -- an abbreviated form of the species name
!
! aux_names        -- the name of the auxiliary variable
! short_aux_names  -- an abbreviated form of the auxiliary variable


module network_properties

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, parameter :: nspec = 13
  integer, parameter :: naux =  0

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len=16), save :: aux_names(naux)
  character (len= 5), save :: short_aux_names(naux)

  double precision, allocatable, save :: aion(:), zion(:), nion(:)

  !$acc declare create(aion, zion, nion)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: aion, zion, nion
#endif

contains

  subroutine network_properties_init

    spec_names(1) = "helium-4"
    spec_names(2) = "carbon-12"
    spec_names(3) = "oxygen-16"
    spec_names(4) = "neon-20"
    spec_names(5) = "magnesium-24"
    spec_names(6) = "silicon-28"
    spec_names(7) = "sulfur-32"
    spec_names(8) = "argon-36"
    spec_names(9) = "calcium-40"
    spec_names(10) = "titanium-44"
    spec_names(11) = "chromium-48"
    spec_names(12) = "iron-52"
    spec_names(13) = "nickel-56"

    short_spec_names(1) = "He4"
    short_spec_names(2) = "C12"
    short_spec_names(3) = "O16"
    short_spec_names(4) = "Ne20"
    short_spec_names(5) = "Mg24"
    short_spec_names(6) = "Si28"
    short_spec_names(7) = "S32"
    short_spec_names(8) = "Ar36"
    short_spec_names(9) = "Ca40"
    short_spec_names(10) = "Ti44"
    short_spec_names(11) = "Cr48"
    short_spec_names(12) = "Fe52"
    short_spec_names(13) = "Ni56"

    allocate(aion(nspec))
    allocate(zion(nspec))
    allocate(nion(nspec))

    aion(1) = 4.0_rt
    aion(2) = 12.0_rt
    aion(3) = 16.0_rt
    aion(4) = 20.0_rt
    aion(5) = 24.0_rt
    aion(6) = 28.0_rt
    aion(7) = 32.0_rt
    aion(8) = 36.0_rt
    aion(9) = 40.0_rt
    aion(10) = 44.0_rt
    aion(11) = 48.0_rt
    aion(12) = 52.0_rt
    aion(13) = 56.0_rt

    zion(1) = 2.0_rt
    zion(2) = 6.0_rt
    zion(3) = 8.0_rt
    zion(4) = 10.0_rt
    zion(5) = 12.0_rt
    zion(6) = 14.0_rt
    zion(7) = 16.0_rt
    zion(8) = 18.0_rt
    zion(9) = 20.0_rt
    zion(10) = 22.0_rt
    zion(11) = 24.0_rt
    zion(12) = 26.0_rt
    zion(13) = 28.0_rt

    ! Set the number of neutrons
    nion(:) = aion(:) - zion(:)



    !$acc update device(aion, zion)

  end subroutine network_properties_init



  subroutine network_properties_finalize

    implicit none

    if (allocated(aion)) then
       deallocate(aion)
    end if

    if (allocated(zion)) then
       deallocate(zion)
    end if

    if (allocated(nion)) then
       deallocate(nion)
    end if

  end subroutine network_properties_finalize

end module network_properties
