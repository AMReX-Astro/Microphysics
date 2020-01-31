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

    spec_names(1) = "hydrogen-1"
    spec_names(2) = "helium-4"
    spec_names(3) = "oxygen-14"
    spec_names(4) = "oxygen-15"
    spec_names(5) = "neon-18"
    spec_names(6) = "silicon-25"
    spec_names(7) = "iron-56"

    short_spec_names(1) = "H1"
    short_spec_names(2) = "He4"
    short_spec_names(3) = "O14"
    short_spec_names(4) = "O15"
    short_spec_names(5) = "Ne18"
    short_spec_names(6) = "Si25"
    short_spec_names(7) = "Fi56"

    allocate(aion(nspec))
    allocate(zion(nspec))
    allocate(nion(nspec))

    aion(1) = 1.0
    aion(2) = 4.0
    aion(3) = 14.0
    aion(4) = 15.0
    aion(5) = 18.0
    aion(6) = 25.0
    aion(7) = 56.0

    zion(1) = 1.0
    zion(2) = 2.0
    zion(3) = 8.0
    zion(4) = 8.0
    zion(5) = 10.0
    zion(6) = 14.0
    zion(7) = 26.0

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
