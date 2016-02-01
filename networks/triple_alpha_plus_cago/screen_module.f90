! This module contains the screen routine which applies screening corrections
! to the reaction rates.  The input composition must be in terms of molar 
! fractions.  The triple-alpha reaction goes through a compound nucleus 
! channel and therefore screening must be applied to both reactions.
!
! A call is made to screenz which lives in screen.f to apply the screening
! to a single reaction.
module screen_module

  use bl_types
  use bl_constants_module
  use network

  implicit none

contains

  subroutine screen(temp, dens, ymol, rates, dratesdt)

    use screening_module, only: screenz

    real(kind=dp_t), intent(IN   ) :: temp, dens, ymol(nspec)
    real(kind=dp_t), intent(INOUT) :: rates(nrat)
    real(kind=dp_t), intent(INOUT) :: dratesdt(nrat)

    real(kind=dp_t) :: scorr1, dscorr1dt
    real(kind=dp_t) :: scorr2, dscorr2dt
    real(kind=dp_t) ::  scorr,  dscorrdt

    real(kind=dp_t) :: rates_in(nrat), dratesdt_in(nrat)

    
    rates_in    = rates
    dratesdt_in = dratesdt
    
    ! triple alpha going through compound nucleus channel
    call screenz(temp, dens,              &
         zion(ihe4_), zion(ihe4_),          &
         aion(ihe4_), aion(ihe4_),          &
         ymol, aion, zion, nspec, scorr1, dscorr1dt)

    call screenz(temp, dens,              &
         zion(ihe4_), FOUR,                &
         aion(ihe4_), EIGHT,               &
         ymol, aion, zion, nspec, scorr2, dscorr2dt)

    scorr    = scorr1 * scorr2
    dscorrdt = dscorr1dt * scorr2 + scorr1 * dscorr2dt

    rates(ir3a_)    = rates_in(ir3a_) * scorr
    dratesdt(ir3a_) = dratesdt_in(ir3a_) * scorr +                        &
                     rates_in(ir3a_) * dscorrdt

    ! C12 + alpha --> O16
    call screenz(temp, dens,     &
         zion(ic12_), zion(ihe4_), &
         aion(ic12_), aion(ihe4_), &
         ymol, aion, zion, nspec, scorr, dscorrdt)

    rates(ircago_)    = rates_in(ircago_) * scorr
    dratesdt(ircago_) = dratesdt_in(ircago_) * scorr + &
                        rates_in(ircago_) * dscorrdt

    return

  end subroutine screen

end module screen_module
