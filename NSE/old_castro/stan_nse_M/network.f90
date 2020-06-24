! the network module provides the information about the species we are
! advecting: 
!
! nspec      -- the number of species
!
! aion       -- atomic number
! zion       -- proton number
! eion       -- nuclear binding energy (in erg/g)
!
! spec_names -- the name of the isotope
! short_spec_names -- the abbreviated name of the isotope
!
! This module contains two routines:
!
!  network_init()        -- initialize the isotope properties
!
!  network_species_index -- return the index of the species given its name
!

module network

  use bl_types

  implicit none

  integer, parameter :: nspec = 7
  integer, parameter :: naux  = 3

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)

  character (len=16), save ::  aux_names(naux)
  character (len= 5), save :: short_aux_names(naux)

  real(kind=dp_t), save :: aion(nspec+naux), zion(nspec+naux), ebin(nspec+naux)

  logical, save :: network_initialized = .false.

contains
  
  subroutine network_init()

    integer :: ic12, io16, ihe4, ine20, img24, isi28, ife56
    integer :: iye, ibea, iabar

    ! integer keys -- for convinence.  In all other places, we will find
    ! these by querying based on species name using network_species_index
    ic12  = 1
    io16  = 2
    ihe4  = 3
    ine20 = 4
    img24 = 5
    isi28 = 6
    ife56 = 7

    iye   = 1
    ibea  = 2
    iabar = 3

    spec_names(ic12)  = "carbon-12"
    short_spec_names(ic12)  = "C12"

    spec_names(io16)  = "oxygen-16"
    short_spec_names(io16)  = "O16"

    spec_names(ihe4)  = "helium-4"
    short_spec_names(ihe4)  = "He4"

    spec_names(ine20)  = "neon-20"
    short_spec_names(ine20)  = "Ne20"

    spec_names(img24)  = "magnesium-24"
    short_spec_names(img24)  = "Mg24"

    spec_names(isi28)  = "silicon-28"
    short_spec_names(isi28)  = "Si28"

    spec_names(ife56)  = "iron-56"
    short_spec_names(ife56)  = "Fe56"


    aux_names(iye)  = "Ye"
    short_aux_names(iye)  = "Ye"

    aux_names(ibea)  = "BeA"
    short_aux_names(ibea)  = "BeA"

    aux_names(iabar) = "Abar"   
    short_aux_names(iabar) = "Abar" 

    
    network_initialized = .true.

  end subroutine network_init

  function network_species_index(name) result(r)

    character(len=*) :: name
    integer :: r, n

    r = -1

    do n = 1, nspec
       if (trim(name) == trim(spec_names(n))) then
          r = n
          exit
       endif
    enddo

    do n = 1, naux
       if (trim(name) == trim(aux_names(n))) then
          r = n+nspec
          exit
       endif
    enddo

    return
  end function network_species_index

end module network
