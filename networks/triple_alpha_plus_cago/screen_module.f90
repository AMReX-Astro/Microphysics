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

    !$acc routine seq
    
    use screening_module, only: screenz

    real(kind=dp_t), intent(IN   ) :: temp, dens, ymol(nspec)
    real(kind=dp_t), intent(INOUT) :: rates(nrates)
    real(kind=dp_t), intent(INOUT) :: dratesdt(nrates)

    real(kind=dp_t) :: scorr1, dscorr1dt
    real(kind=dp_t) :: scorr2, dscorr2dt
    real(kind=dp_t) ::  scorr,  dscorrdt

    real(kind=dp_t) :: rates_in(nrates), dratesdt_in(nrates)

    
    rates_in    = rates
    dratesdt_in = dratesdt
    
    ! triple alpha going through compound nucleus channel
    call screenz(temp, dens,              &
         zion(ihe4), zion(ihe4),          &
         aion(ihe4), aion(ihe4),          &
         ymol, scorr1, dscorr1dt)

    call screenz(temp, dens,              &
         zion(ihe4), FOUR,                &
         aion(ihe4), EIGHT,               &
         ymol, scorr2, dscorr2dt)

    scorr    = scorr1 * scorr2
    dscorrdt = dscorr1dt * scorr2 + scorr1 * dscorr2dt

    rates(ir3a)    = rates_in(ir3a) * scorr
    dratesdt(ir3a) = dratesdt_in(ir3a) * scorr +                        &
                     rates_in(ir3a) * dscorrdt

    ! C12 + alpha --> O16
    call screenz(temp, dens,     &
         zion(ic12), zion(ihe4), &
         aion(ic12), aion(ihe4), &
         ymol, scorr, dscorrdt)

    rates(ircago)    = rates_in(ircago) * scorr
    dratesdt(ircago) = dratesdt_in(ircago) * scorr + &
                        rates_in(ircago) * dscorrdt

    return

  end subroutine screen

end module screen_module
