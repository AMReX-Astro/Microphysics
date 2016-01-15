module network_indices

  implicit none

  ! ids of the species
  
  integer, parameter :: ihe4  = 1
  integer, parameter :: ic12  = 2
  integer, parameter :: io16  = 3
  integer, parameter :: ine20 = 4
  integer, parameter :: img24 = 5
  integer, parameter :: isi28 = 6
  integer, parameter :: ini56 = 7

  ! set the id numbers of the reaction rates

  integer, parameter :: ircag   = 1
  integer, parameter :: iroga   = 2
  integer, parameter :: ir3a    = 3
  integer, parameter :: irg3a   = 4
  integer, parameter :: ir1212  = 5
  integer, parameter :: ir1216  = 6
  integer, parameter :: ir1616  = 7
  integer, parameter :: iroag   = 8
  integer, parameter :: irnega  = 9
  integer, parameter :: irneag  = 10
  integer, parameter :: irmgga  = 11
  integer, parameter :: irmgag  = 12
  integer, parameter :: irsiga  = 13
  integer, parameter :: ircaag  = 14
  integer, parameter :: irtiga  = 15

  integer, parameter :: irsi2ni = 16
  integer, parameter :: irni2si = 17

  ! rpar data

  integer :: irp_dydt
  integer :: irp_rates
  integer :: irp_drdy1
  integer :: irp_drdy2

end module network_indices
