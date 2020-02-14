! An automatically generated file of network properties.  This provides the properties
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

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, parameter :: nspec = 19
  integer, parameter :: naux =  0

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len=16), save :: aux_names(naux)
  character (len= 5), save :: short_aux_names(naux)

  real(rt), allocatable, save :: aion(:), zion(:), nion(:)

  !$acc declare create(aion, zion, nion)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: aion, zion, nion
#endif

contains

  subroutine network_properties_init

    spec_names(1) = "hydrogen-1"
    spec_names(2) = "helium-3"
    spec_names(3) = "helium-4"
    spec_names(4) = "carbon-12"
    spec_names(5) = "nitrogen-14"
    spec_names(6) = "oxygen-16"
    spec_names(7) = "neon-20"
    spec_names(8) = "magnesium-24"
    spec_names(9) = "silicon-28"
    spec_names(10) = "sulfur-32"
    spec_names(11) = "argon-36"
    spec_names(12) = "calcium-40"
    spec_names(13) = "titanium-44"
    spec_names(14) = "chromium-48"
    spec_names(15) = "iron-52"
    spec_names(16) = "iron-54"
    spec_names(17) = "nickel-56"
    spec_names(18) = "neutron"
    spec_names(19) = "proton"

    short_spec_names(1) = "H1"
    short_spec_names(2) = "He3"
    short_spec_names(3) = "He4"
    short_spec_names(4) = "C12"
    short_spec_names(5) = "N14"
    short_spec_names(6) = "O16"
    short_spec_names(7) = "Ne20"
    short_spec_names(8) = "Mg24"
    short_spec_names(9) = "Si28"
    short_spec_names(10) = "S32"
    short_spec_names(11) = "Ar36"
    short_spec_names(12) = "Ca40"
    short_spec_names(13) = "Ti44"
    short_spec_names(14) = "Cr48"
    short_spec_names(15) = "Fe52"
    short_spec_names(16) = "Fe54"
    short_spec_names(17) = "Ni56"
    short_spec_names(18) = "n"
    short_spec_names(19) = "p"

    allocate(aion(nspec))
    allocate(zion(nspec))
    allocate(nion(nspec))

    aion(1) = 1.0_rt
    aion(2) = 3.0_rt
    aion(3) = 4.0_rt
    aion(4) = 12.0_rt
    aion(5) = 14.0_rt
    aion(6) = 16.0_rt
    aion(7) = 20.0_rt
    aion(8) = 24.0_rt
    aion(9) = 28.0_rt
    aion(10) = 32.0_rt
    aion(11) = 36.0_rt
    aion(12) = 40.0_rt
    aion(13) = 44.0_rt
    aion(14) = 48.0_rt
    aion(15) = 52.0_rt
    aion(16) = 54.0_rt
    aion(17) = 56.0_rt
    aion(18) = 1.0_rt
    aion(19) = 1.0_rt

    zion(1) = 1.0_rt
    zion(2) = 2.0_rt
    zion(3) = 2.0_rt
    zion(4) = 6.0_rt
    zion(5) = 7.0_rt
    zion(6) = 8.0_rt
    zion(7) = 10.0_rt
    zion(8) = 12.0_rt
    zion(9) = 14.0_rt
    zion(10) = 16.0_rt
    zion(11) = 18.0_rt
    zion(12) = 20.0_rt
    zion(13) = 22.0_rt
    zion(14) = 24.0_rt
    zion(15) = 26.0_rt
    zion(16) = 26.0_rt
    zion(17) = 28.0_rt
    zion(18) = 0.0_rt
    zion(19) = 1.0_rt

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
