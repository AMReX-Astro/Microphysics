module actual_burner_data

  implicit none

contains

  subroutine ener_gener_rate(dydt, enuc)

    use network

    implicit none

    double precision :: dydt(nspec), enuc

    enuc = sum(dydt(:) * ebin(:))

  end subroutine ener_gener_rate

end module actual_burner_data
