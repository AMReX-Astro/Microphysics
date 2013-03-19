! the network module provides the information about the species we are
! advecting: 
!
! nspec      -- the number of species
!
! aion       -- atomic number
! zion       -- proton number
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
  use bl_constants_module

  implicit none

  ! conversion factors for the nuclear energy generation rate detlap
  ! is the mass excess of the proton in amu detlan is the mass excess
  ! of the neutron in amu -- these come from the original const.dek
  ! from the public network
  real (kind=dp_t), parameter :: avo     = 6.0221417930d23
  real (kind=dp_t), parameter :: clight  = 2.99792458d10
  real (kind=dp_t), parameter :: kerg    = 1.380650424d-16
  real (kind=dp_t), parameter :: ev2erg  = 1.60217648740d-12
  real (kind=dp_t), parameter :: kev     = kerg/ev2erg
  real (kind=dp_t), parameter :: amu     = 1.66053878283d-24

  real (kind=dp_t), parameter :: mn      = 1.67492721184d-24
  real (kind=dp_t), parameter :: mp      = 1.67262163783d-24
  real (kind=dp_t), parameter :: me      = 9.1093821545d-28

  real (kind=dp_t), parameter :: deltap     = 7.288969d0
  real (kind=dp_t), parameter :: deltan     = 8.071323d0

  real (kind=dp_t), parameter :: enuc_conv  = ev2erg*1.0d6*avo
  real (kind=dp_t), parameter :: enuc_conv2 = -avo*clight*clight

  real(kind=dp_t), parameter :: mev2erg = ev2erg*1.0d6
  real(kind=dp_t), parameter :: mev2gr  = mev2erg/clight**2

  real(kind=dp_t), parameter :: pi = M_PI

  character (len=*), parameter :: network_name = "aprox13"

  integer, parameter :: nspec = 13
  integer, parameter :: naux  = 0
  integer, parameter :: nrat = 67

  character (len=16), save :: spec_names(nspec) 
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  real(kind=dp_t), save :: aion(nspec), zion(nspec), nion(nspec)
  real(kind=dp_t), save :: bion(nspec), mion(nspec), wion(nspec)

  character (len=16), save :: ratnam(nrat)


  logical, save :: network_initialized = .false.

contains
  
  subroutine network_init()

    use network_indices

    integer :: i

    ! the following comes directly from init_aprox13

    short_spec_names(ihe4)  = 'he4'
    short_spec_names(ic12)  = 'c12'
    short_spec_names(io16)  = 'o16'
    short_spec_names(ine20) = 'ne20'
    short_spec_names(img24) = 'mg24'
    short_spec_names(isi28) = 'si28'
    short_spec_names(is32)  = 's32'
    short_spec_names(iar36) = 'ar36'
    short_spec_names(ica40) = 'ca40'
    short_spec_names(iti44) = 'ti44'
    short_spec_names(icr48) = 'cr48'
    short_spec_names(ife52) = 'fe52'
    short_spec_names(ini56) = 'ni56'

    spec_names(:) = short_spec_names(:)

    ! set the number of nucleons in the element
    aion(ihe4)  = 4.0d0
    aion(ic12)  = 12.0d0
    aion(io16)  = 16.0d0
    aion(ine20) = 20.0d0
    aion(img24) = 24.0d0
    aion(isi28) = 28.0d0
    aion(is32)  = 32.0d0
    aion(iar36) = 36.0d0
    aion(ica40) = 40.0d0
    aion(iti44) = 44.0d0
    aion(icr48) = 48.0d0
    aion(ife52) = 52.0d0
    aion(ini56) = 56.0d0

    ! set the number of protons in the element
    zion(ihe4)  = 2.0d0
    zion(ic12)  = 6.0d0
    zion(io16)  = 8.0d0
    zion(ine20) = 10.0d0
    zion(img24) = 12.0d0
    zion(isi28) = 14.0d0
    zion(is32)  = 16.0d0
    zion(iar36) = 18.0d0
    zion(ica40) = 20.0d0
    zion(iti44) = 22.0d0
    zion(icr48) = 24.0d0
    zion(ife52) = 26.0d0
    zion(ini56) = 28.0d0

    ! set the binding energy of the element (MeV)
    bion(ihe4)  =  28.29603d0
    bion(ic12)  =  92.16294d0
    bion(io16)  = 127.62093d0
    bion(ine20) = 160.64788d0
    bion(img24) = 198.25790d0
    bion(isi28) = 236.53790d0
    bion(is32)  = 271.78250d0
    bion(iar36) = 306.72020d0
    bion(ica40) = 342.05680d0
    bion(iti44) = 375.47720d0
    bion(icr48) = 411.46900d0
    bion(ife52) = 447.70800d0
    bion(ini56) = 484.00300d0

    do i=1, nspec
       ! set the number of neutrons and mass
       nion(i) = aion(i) - zion(i)

       mion(i) = nion(i)*mn + zion(i)*(mp+me) - bion(i)*mev2gr

       ! molar mass
       wion(i) = avo * mion(i)

       ! a common approximation
       wion(i) = aion(i)
    enddo
    
    ! set the names of the reaction rates 
    ratnam(ir3a)   = 'r3a  '
    ratnam(irg3a)  = 'rg3a '
    ratnam(ircag)  = 'rcag '
    ratnam(ir1212) = 'r1212'
    ratnam(ir1216) = 'r1216'
    ratnam(ir1616) = 'r1616'
    ratnam(iroga)  = 'roga '
    ratnam(iroag)  = 'roag '
    ratnam(irnega) = 'rnega'
    ratnam(irneag) = 'rneag'
    ratnam(irmgga) = 'rmgga'
    ratnam(irmgag) = 'rmgag'
    ratnam(irsiga) = 'rsiga'
    ratnam(irmgap) = 'rmgap'
    ratnam(iralpa) = 'ralpa'
    ratnam(iralpg) = 'ralpg'
    ratnam(irsigp) = 'rsigp'
    ratnam(irsiag) = 'rsiag'
    ratnam(irsga)  = 'rsga '
    ratnam(irsiap) = 'rsiap'
    ratnam(irppa)  = 'rppa '
    ratnam(irppg)  = 'rppg '
    ratnam(irsgp)  = 'rsgp '
    ratnam(irsag)  = 'rsag '
    ratnam(irarga) = 'rarga'
    ratnam(irsap)  = 'rsap '
    ratnam(irclpa) = 'rclpa'
    ratnam(irclpg) = 'rclpg'
    ratnam(irargp) = 'rargp'
    ratnam(irarag) = 'rarag'
    ratnam(ircaga) = 'rcaga'
    ratnam(irarap) = 'rarap'
    ratnam(irkpa)  = 'rkpa '
    ratnam(irkpg)  = 'rkpg '
    ratnam(ircagp) = 'rcagp'
    ratnam(ircaag) = 'rcaag'
    ratnam(irtiga) = 'rtiga'
    ratnam(ircaap) = 'rcaap'
    ratnam(irscpa) = 'rscpa'
    ratnam(irscpg) = 'rscpg'
    ratnam(irtigp) = 'rtigp'
    ratnam(irtiag) = 'rtiag'
    ratnam(ircrga) = 'rcrga'
    ratnam(irtiap) = 'rtiap'
    ratnam(irvpa)  = 'rvpa '
    ratnam(irvpg)  = 'rvpg '
    ratnam(ircrgp) = 'rcrgp'
    ratnam(ircrag) = 'rcrag'
    ratnam(irfega) = 'rfega'
    ratnam(ircrap) = 'rcrap'
    ratnam(irmnpa) = 'rmnpa'
    ratnam(irmnpg) = 'rmnpg'
    ratnam(irfegp) = 'rfegp'
    ratnam(irfeag) = 'rfeag'
    ratnam(irniga) = 'rniga'
    ratnam(irfeap) = 'rfeap'
    ratnam(ircopa) = 'rcopa'
    ratnam(ircopg) = 'rcopg'
    ratnam(irnigp) = 'rnigp'

    ratnam(irr1)   = 'r1   '
    ratnam(irs1)   = 's1   '
    ratnam(irt1)   = 't1   '
    ratnam(iru1)   = 'u1   '
    ratnam(irv1)   = 'v1   '
    ratnam(irw1)   = 'w1   '
    ratnam(irx1)   = 'x1   '
    ratnam(iry1)   = 'y1   '

    network_initialized = .true.

  end subroutine network_init


  function network_species_index(name) result(r)

    character(len=*) :: name
    integer :: r, n

    r = -1

    do n = 1, nspec
       if (trim(name) == trim(spec_names(n)) .or. &
           trim(name) == trim(short_spec_names(n))) then
          r = n
          exit
       endif
    enddo
    return
  end function network_species_index


  subroutine ener_gener_rate(dydt,enuc)

    ! computes the instantaneous energy generation rate

    ! declare the pass
    double precision dydt(nspec), enuc
    
    ! local variables
    integer          i
    
    ! instantaneous energy generation rate

    ! this form misses n <-> p differences 
    
    ! enuc = 0.0d0 
    ! do i=1,ionmax 
    !   enuc = enuc + dydt(i) * bion(i) 
    ! enddo 
    ! enuc = enuc * enuc_conv


    ! this form gets the n <-> p differences 
    
    ! enuc = 0.0d0 
    ! do i=1,ionmax
    !      enuc = enuc + dydt(i) * (bion(i) - zion(i)*deltap - nion(i)*deltan)
    ! enddo 
    ! enuc = enuc * enuc_conv
    
    ! this form is closest to e = m c**2 and gives the same results as
    ! the form above
    
    enuc = 0.0d0
    do i=1,nspec
       enuc = enuc + dydt(i) * mion(i)
    enddo
    enuc = enuc * enuc_conv2
    
    return
  end subroutine ener_gener_rate


  subroutine network_finalize()

  end subroutine network_finalize

end module network
