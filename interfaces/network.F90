! the network module provides the information about the species we are
! advecting:
!
! This module contains the following routines:
!
!  network_init          -- initialize the isotope properties
!
!  network_species_index -- return the index of the species given its name
!
!  network_finalize      -- do any network cleanup

module network

  use amrex_fort_module, only : rt => amrex_real
  use network_properties

  implicit none

  logical :: network_initialized = .false.

  ! temporary hack until we eliminate nspec_evolve from application codes
  integer, parameter :: nspec_evolve = nspec

contains

  subroutine network_init()

    use amrex_error_module, only : amrex_error
    use amrex_constants_module, only : ONE

    implicit none

    call network_properties_init()

    if ( nspec .le. 0 ) then
       call amrex_error("Network cannot have a negative number of species.")
    endif

    if ( naux .lt. 0 ) then
       call amrex_error("Network cannot have a negative number of auxiliary variables.")
    endif

    network_initialized = .true.

  end subroutine network_init


  function network_species_index(name) result(r)

    character(len=*) :: name
    integer :: r, n

    r = -1

    do n = 1, nspec
       if (name == spec_names(n) .or. name == short_spec_names(n)) then
          r = n
          return
       endif
    enddo

  end function network_species_index


  function network_aux_index(name) result(r)

    character(len=*) :: name
    integer :: r, n

    r = -1

    do n = 1, naux
       if (name == aux_names(n) .or. name == short_aux_names(n)) then
          r = n
          return
       endif
    enddo

  end function network_aux_index


  function get_network_species_name(index) result(name)

    character(len=128) :: name
    integer :: index

    if (index < 1 .or. index > nspec) then
       name = ""
    else
       name = spec_names(index)
    endif

  end function get_network_species_name


  function get_network_short_species_name(index) result(name)

    character(len=128) :: name
    integer :: index

    if (index < 1 .or. index > nspec) then
       name = ""
    else
       name = short_spec_names(index)
    endif

  end function get_network_short_species_name


  subroutine network_finalize()
    implicit none

  end subroutine network_finalize

end module network
