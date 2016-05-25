module actual_burner_data

  implicit none

  integer, parameter :: nrates = 1

  ! Conversion factor for the nuclear energy generation rate.

  double precision, parameter :: avo = 6.0221417930d23
  double precision, parameter :: c_light = 2.99792458d10
  double precision, parameter :: enuc_conv2 = -avo*c_light*c_light

contains

  ! Computes the instantaneous energy generation rate

  subroutine ener_gener_rate(dydt, enuc)

    use network

    !$acc routine seq

    implicit none

    double precision :: dydt(nspec_evolve), enuc

    ! This is basically e = m c**2

    ! Note that since we don't explicitly evolve Mg24
    ! in this network, we need to explicitly add its
    ! contribution in this routine.

    enuc = dydt(ic12) * (mion(img24) - mion(ic12)) * enuc_conv2

  end subroutine ener_gener_rate

end module actual_burner_data
