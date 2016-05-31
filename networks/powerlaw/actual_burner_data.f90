module actual_burner_data

  implicit none

contains

  ! Computes the instantaneous energy generation rate

  subroutine ener_gener_rate(dydt, enuc)

    use network

    implicit none

    double precision :: dydt(nspec_evolve), enuc

    ! This is basically e = m c**2

    enuc = sum(dydt(:) * aion(1:nspec_evolve) * ebin(1:nspec_evolve))

  end subroutine ener_gener_rate

end module actual_burner_data
