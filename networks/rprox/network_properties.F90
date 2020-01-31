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

  integer, parameter :: nspec = 10
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

    spec_names(1) = "carbon-12"
    spec_names(2) = "oxygen-14"
    spec_names(3) = "oxygen-15"
    spec_names(4) = "oxygen-16"
    spec_names(5) = "flourine-17"
    spec_names(6) = "magnesium-22"
    spec_names(7) = "sulfur-30"
    spec_names(8) = "nickel-56"
    spec_names(9) = "hydrogen-1"
    spec_names(10) = "helium-4"

    short_spec_names(1) = "C12"
    short_spec_names(2) = "O14"
    short_spec_names(3) = "O15"
    short_spec_names(4) = "O16"
    short_spec_names(5) = "F17"
    short_spec_names(6) = "Mg22"
    short_spec_names(7) = "S30"
    short_spec_names(8) = "Ni56"
    short_spec_names(9) = "H1"
    short_spec_names(10) = "He4"

    allocate(aion(nspec))
    allocate(zion(nspec))
    allocate(nion(nspec))

    aion(1) = 12.0
    aion(2) = 14.0
    aion(3) = 15.0
    aion(4) = 16.0
    aion(5) = 17.0
    aion(6) = 22.0
    aion(7) = 30.0
    aion(8) = 56.0
    aion(9) = 1.0
    aion(10) = 4.0

    zion(1) = 6.0
    zion(2) = 8.0
    zion(3) = 8.0
    zion(4) = 8.0
    zion(5) = 9.0
    zion(6) = 12.0
    zion(7) = 16.0
    zion(8) = 28.0
    zion(9) = 1.0
    zion(10) = 2.0

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
