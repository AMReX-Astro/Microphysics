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

  integer, parameter :: nspec = 3
  integer, parameter :: nspec_evolve = 3
  integer, parameter :: naux =  0

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len=16), save :: aux_names(naux)
  character (len= 5), save :: short_aux_names(naux)

  double precision, allocatable, save :: aion(:), zion(:)

  !$acc declare create(aion, zion)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: aion, zion
#endif

contains

  subroutine network_properties_init

    spec_names(1) = "fuel"
    spec_names(2) = "ash"
    spec_names(3) = "inert"

    short_spec_names(1) = "fuel"
    short_spec_names(2) = "ash"
    short_spec_names(3) = "inert"

    allocate(aion(nspec))
    allocate(zion(nspec))

    aion(1) = 2.0
    aion(2) = 4.0
    aion(3) = 8.0

    zion(1) = 1.0
    zion(2) = 2.0
    zion(3) = 4.0


    !$acc update device(aion, zion)

  end subroutine network_properties_init



  subroutine network_properties_finalize

    implicit none

    deallocate(aion)
    deallocate(zion)

  end subroutine network_properties_finalize

end module network_properties
