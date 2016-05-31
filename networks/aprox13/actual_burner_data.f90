module actual_burner_data

  use actual_network_data, only: nrates

  implicit none

  ! Table interpolation data

  double precision, parameter :: tab_tlo = 6.0d0, tab_thi = 10.0d0
  integer, parameter :: tab_per_decade = 500
  integer, parameter :: nrattab = int(tab_thi - tab_tlo) * tab_per_decade + 1
  integer, parameter :: tab_imax = int(tab_thi - tab_tlo) * tab_per_decade + 1
  double precision, parameter :: tab_tstp = (tab_thi - tab_tlo) / dble(tab_imax - 1)

  double precision :: rattab(nrates, nrattab)
  double precision :: drattabdt(nrates, nrattab)
  double precision :: drattabdd(nrates, nrattab)
  double precision :: ttab(nrattab)

  !$acc declare create(rattab, drattabdt, drattabdd, ttab)

contains

  ! Computes the instantaneous energy generation rate

  subroutine ener_gener_rate(dydt, enuc)

    !$acc routine seq

    use actual_network_data, only: nspec, mion, enuc_conv2

    implicit none

    double precision :: dydt(nspec), enuc

    ! This is basically e = m c**2

    enuc = sum(dydt(:) * mion(:)) * enuc_conv2

  end subroutine ener_gener_rate

end module actual_burner_data
