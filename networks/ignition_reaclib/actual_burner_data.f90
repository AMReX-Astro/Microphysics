module actual_burner_data

  implicit none

  integer, parameter :: nrates = 4

  ! Conversion factor for the nuclear energy generation rate.

  double precision, parameter :: avo = 6.0221417930d23
  double precision, parameter :: c_light = 2.99792458d10
  double precision, parameter :: enuc_conv2 = -avo*c_light*c_light

  !$acc declare create(avo, c_light, enuc_conv2)

end module actual_burner_data
