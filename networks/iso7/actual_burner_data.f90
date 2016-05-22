module actual_burner_data

  integer, parameter :: nrates  = 17

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

  character (len=20), save :: ratenames(nrates)

  integer, save :: irp_dydt
  integer, save :: irp_rates
  integer, save :: irp_drdy1
  integer, save :: irp_drdy2

contains

  ! Computes the instantaneous energy generation rate

  subroutine ener_gener_rate(dydt, enuc)

    use actual_network_data, only: nspec, mion, enuc_conv2

    implicit none

    double precision :: dydt(nspec), enuc

    ! This is basically e = m c**2

    enuc = sum(dydt(:) * mion(:)) * enuc_conv2

  end subroutine ener_gener_rate

end module actual_burner_data
