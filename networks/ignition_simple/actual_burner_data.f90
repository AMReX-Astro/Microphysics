module actual_burner_data

  implicit none

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
    ! contribution in this routine. We can factor out
    ! the common factor of dydt(ic12), we just need to
    ! account for a factor of aion(ic12) / aion(img24)
    ! for the second term to make the expression work.

    enuc = dydt(ic12) * (mion(ic12) - mion(img24) / 2) * enuc_conv2

  end subroutine ener_gener_rate

end module actual_burner_data
