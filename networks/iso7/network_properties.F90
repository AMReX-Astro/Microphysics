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

  integer, parameter :: nspec = 7
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
    spec_names(7) = "nickel-56"

    short_spec_names(1) = "He4"
    short_spec_names(2) = "C12"
    short_spec_names(3) = "O16"
    short_spec_names(4) = "Ne20"
    short_spec_names(5) = "Mg24"
    short_spec_names(6) = "Si28"
    short_spec_names(7) = "Ni56"

    allocate(aion(nspec))
    allocate(zion(nspec))
    allocate(nion(nspec))

    aion(1) = 4.0
    aion(2) = 12.0
    aion(3) = 16.0
    aion(4) = 20.0
    aion(5) = 24.0
    aion(6) = 28.0
    aion(7) = 56.0

    zion(1) = 2.0
    zion(2) = 6.0
    zion(3) = 8.0
    zion(4) = 10.0
    zion(5) = 12.0
    zion(6) = 14.0
    zion(7) = 28.0

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
