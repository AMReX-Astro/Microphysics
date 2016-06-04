! This module contains the dydt routine which change in species with respect 
! to time based on the reaction rates.
!
! Note that the rates are determined in terms of molar fractions and as such
! must be converted back to mass fractions before exiting.
!
module dydt_module

  use bl_types
  use bl_constants_module
  use network

  implicit none

contains 

  subroutine dydt(ymol, rates, ydot)

    double precision, intent(IN   ) :: ymol(nspec), rates(nrates)
    double precision, intent(  OUT) :: ydot(nspec)

    ydot(:) = ZERO

    ! He4 reactions
    ydot(ihe4_) = - THREE * ymol(ihe4_)**THREE * rates(ir3a_)              &
                 - ONE * ymol(ic12_) * ymol(ihe4_) * rates(ircago_)

    ! C12 reactions
    ydot(ic12_) =  ONE * ymol(ihe4_)**THREE * rates(ir3a_)                &
                 - ONE * ymol(ic12_) * ymol(ihe4_) * rates(ircago_)

    ! O16 reactions
    ydot(io16_) =  ONE * ymol(ic12_) * ymol(ihe4_) * rates(ircago_)

  end subroutine dydt


end module dydt_module
