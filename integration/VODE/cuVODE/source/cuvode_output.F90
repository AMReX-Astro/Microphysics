module cuvode_output_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, parameter :: LUNIT = 6

contains

  subroutine xerrwd(MSG, NMES, NERR, LEVEL, NI, I1, I2, NR, R1, R2)

    implicit none

    character (len=80), intent(in) :: MSG
    real(rt), intent(in) :: R1, R2
    integer,  intent(in) :: NMES, NERR, LEVEL, NI, I1, I2, NR

    character (len=50), parameter :: fmt10 = "(1X,A)"
    character (len=50), parameter :: fmt20 = "(6X,'In above message,  I1 =',I10)"
    character (len=50), parameter :: fmt30 = "(6X,'In above message,  I1 =',I10,3X,'I2 =',I10)"
    character (len=50), parameter :: fmt40 = "(6X,'In above message,  R1 =',D21.13)"
    character (len=50), parameter :: fmt50 = "(6X,'In above,  R1 =',D21.13,3X,'R2 =',D21.13)"

    WRITE (LUNIT, fmt10) MSG

    select case(NI)
    case(1)
       WRITE (LUNIT, fmt20) I1
    case(2)
       WRITE (LUNIT, fmt30) I1, I2
    end select

    select case(NR)
    case(1)
       WRITE (LUNIT, fmt40) R1
    case(2)
       WRITE (LUNIT, fmt50) R1, R2
    end select

    IF (LEVEL .NE. 2) RETURN

    STOP

  end subroutine xerrwd
  
end module cuvode_output_module
