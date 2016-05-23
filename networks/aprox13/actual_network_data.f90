module actual_network_data

  implicit none

  integer, parameter :: nspec = 13
  integer, parameter :: nspec_evolve = 13
  integer, parameter :: naux  = 0

  integer, parameter :: ihe4  = 1
  integer, parameter :: ic12  = 2
  integer, parameter :: io16  = 3
  integer, parameter :: ine20 = 4
  integer, parameter :: img24 = 5
  integer, parameter :: isi28 = 6
  integer, parameter :: is32  = 7
  integer, parameter :: iar36 = 8
  integer, parameter :: ica40 = 9
  integer, parameter :: iti44 = 10
  integer, parameter :: icr48 = 11
  integer, parameter :: ife52 = 12
  integer, parameter :: ini56 = 13

  double precision :: aion(nspec), zion(nspec), nion(nspec)
  double precision :: bion(nspec), mion(nspec), wion(nspec)

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  character (len=32), save :: network_name = "aprox13"

  ! Some fundamental physical constants

  double precision, parameter :: avo = 6.0221417930d23
  double precision, parameter :: c_light = 2.99792458d10

  double precision, parameter :: ev2erg  = 1.60217648740d-12
  double precision, parameter :: mev2erg = ev2erg*1.0d6
  double precision, parameter :: mev2gr  = mev2erg/c_light**2

  double precision, parameter :: mn = 1.67492721184d-24
  double precision, parameter :: mp = 1.67262163783d-24
  double precision, parameter :: me = 9.1093821545d-28

  ! Conversion factor for the nuclear energy generation rate.

  double precision, parameter :: enuc_conv2 = -avo*c_light*c_light

  !$acc declare create(aion, zion, nion, bion, mion, wion)

end module actual_network_data
