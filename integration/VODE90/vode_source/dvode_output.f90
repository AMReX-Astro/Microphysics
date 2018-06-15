module dvode_output_module

  use amrex_fort_module, only : rt => amrex_real

  use dvode_constants_module
  
  implicit none

contains

  function iumach() result(ium_return)
    ! ***BEGIN PROLOGUE  IUMACH
    ! ***PURPOSE  Provide standard output unit number.
    ! ***CATEGORY  R1
    ! ***TYPE      INTEGER (IUMACH-I)
    ! ***KEYWORDS  MACHINE CONSTANTS
    ! ***AUTHOR  Hindmarsh, Alan C., (LLNL)
    ! ***DESCRIPTION
    !  *Usage:
    !         INTEGER  LOUT, IUMACH
    !         LOUT = IUMACH()
    ! 
    !  *Function Return Values:
    !      LOUT : the standard logical unit for Fortran output.
    ! 
    ! ***REFERENCES  (NONE)
    ! ***ROUTINES CALLED  (NONE)
    ! ***REVISION HISTORY  (YYMMDD)
    !    930915  DATE WRITTEN
    !    930922  Made user-callable, and other cosmetic changes. (FNF)
    ! ***END PROLOGUE  IUMACH
    ! 
    ! *Internal Notes:
    !   The built-in value of 6 is standard on a wide range of Fortran
    !   systems.  This may be machine-dependent.
    ! **End
    ! ***FIRST EXECUTABLE STATEMENT  IUMACH
    integer, parameter :: ium = 6
    integer :: ium_return
    ium_return = ium
    return
  end function iumach
  
  function ixsav(IPAR, IVALUE, ISET) result(ixs)
    ! ***BEGIN PROLOGUE  IXSAV
    ! ***SUBSIDIARY
    ! ***PURPOSE  Save and recall error message control parameters.
    ! ***CATEGORY  R3C
    ! ***TYPE      ALL (IXSAV-A)
    ! ***AUTHOR  Hindmarsh, Alan C., (LLNL)
    ! ***DESCRIPTION
    ! 
    !   IXSAV saves and recalls one of two error message parameters:
    !     LUNIT, the logical unit number to which messages are printed, and
    !     MESFLG, the message print flag.
    !   This is a modification of the SLATEC library routine J4SAVE.
    ! 
    !   Saved local variables..
    !    LUNIT  = Logical unit number for messages.  The default is obtained
    !             by a call to IUMACH (may be machine-dependent).
    !    MESFLG = Print control flag..
    !             1 means print all messages (the default).
    !             0 means no printing.
    ! 
    !   On input..
    !     IPAR   = Parameter indicator (1 for LUNIT, 2 for MESFLG).
    !     IVALUE = The value to be set for the parameter, if ISET = .TRUE.
    !     ISET   = Logical flag to indicate whether to read or write.
    !              If ISET = .TRUE., the parameter will be given
    !              the value IVALUE.  If ISET = .FALSE., the parameter
    !              will be unchanged, and IVALUE is a dummy argument.
    ! 
    !   On return..
    !     IXSAV = The (old) value of the parameter.
    ! 
    ! ***SEE ALSO  XERRWD, XERRWV
    ! ***ROUTINES CALLED  IUMACH
    ! ***REVISION HISTORY  (YYMMDD)
    !    921118  DATE WRITTEN
    !    930329  Modified prologue to SLATEC format. (FNF)
    !    930915  Added IUMACH call to get default output unit.  (ACH)
    !    930922  Minor cosmetic changes. (FNF)
    !    010425  Type declaration for IUMACH added. (ACH)
    ! ***END PROLOGUE  IXSAV
    ! 
    !  Subroutines called by IXSAV.. None
    !  Function routine called by IXSAV.. IUMACH
    ! -----------------------------------------------------------------------
    ! **End
    
    logical :: ISET
    integer :: IPAR, IVALUE
    integer :: ixs
    integer, save :: LUNIT = -1
    integer, save :: MESFLG = 1

    select case(IPAR)
    case(1)
       IF (LUNIT .EQ. -1) LUNIT = IUMACH()
       ixs = LUNIT
       IF (ISET) LUNIT = IVALUE
    case(2)
       ixs = MESFLG
       IF (ISET) MESFLG = IVALUE
    case default
       ixs = 0
    end select
  end function ixsav

  subroutine xsetf(MFLAG)
    ! ***BEGIN PROLOGUE  XSETF
    ! ***PURPOSE  Reset the error print control flag.
    ! ***CATEGORY  R3A
    ! ***TYPE      ALL (XSETF-A)
    ! ***KEYWORDS  ERROR CONTROL
    ! ***AUTHOR  Hindmarsh, Alan C., (LLNL)
    ! ***DESCRIPTION
    ! 
    !    XSETF sets the error print control flag to MFLAG:
    !       MFLAG=1 means print all messages (the default).
    !       MFLAG=0 means no printing.
    ! 
    ! ***SEE ALSO  XERRWD, XERRWV
    ! ***REFERENCES  (NONE)
    ! ***ROUTINES CALLED  IXSAV
    ! ***REVISION HISTORY  (YYMMDD)
    !    921118  DATE WRITTEN
    !    930329  Added SLATEC format prologue. (FNF)
    !    930407  Corrected SEE ALSO section. (FNF)
    !    930922  Made user-callable, and other cosmetic changes. (FNF)
    ! ***END PROLOGUE  XSETF
    ! 
    !  Subroutines called by XSETF.. None
    !  Function routine called by XSETF.. IXSAV
    ! -----------------------------------------------------------------------
    ! **End
    integer :: MFLAG, JUNK

    IF (MFLAG .EQ. 0 .OR. MFLAG .EQ. 1) JUNK = IXSAV(2,MFLAG,.TRUE.)
    return
  end subroutine xsetf

  subroutine xsetun(LUN)
    ! ***BEGIN PROLOGUE  XSETUN
    ! ***PURPOSE  Reset the logical unit number for error messages.
    ! ***CATEGORY  R3B
    ! ***TYPE      ALL (XSETUN-A)
    ! ***KEYWORDS  ERROR CONTROL
    ! ***DESCRIPTION
    ! 
    !    XSETUN sets the logical unit number for error messages to LUN.
    ! 
    ! ***AUTHOR  Hindmarsh, Alan C., (LLNL)
    ! ***SEE ALSO  XERRWD, XERRWV
    ! ***REFERENCES  (NONE)
    ! ***ROUTINES CALLED  IXSAV
    ! ***REVISION HISTORY  (YYMMDD)
    !    921118  DATE WRITTEN
    !    930329  Added SLATEC format prologue. (FNF)
    !    930407  Corrected SEE ALSO section. (FNF)
    !    930922  Made user-callable, and other cosmetic changes. (FNF)
    ! ***END PROLOGUE  XSETUN
    ! 
    !  Subroutines called by XSETUN.. None
    !  Function routine called by XSETUN.. IXSAV
    ! -----------------------------------------------------------------------
    ! **End
    integer :: LUN, JUNK

    IF (LUN .GT. 0) JUNK = IXSAV (1,LUN,.TRUE.)
    return
  end subroutine xsetun

  subroutine xerrwd(MSG, NMES, NERR, LEVEL, NI, I1, I2, NR, R1, R2)
    ! ***BEGIN PROLOGUE  XERRWD
    ! ***SUBSIDIARY
    ! ***PURPOSE  Write error message with values.
    ! ***CATEGORY  R3C
    ! ***TYPE      DOUBLE PRECISION (XERRWV-S, XERRWD-D)
    ! ***AUTHOR  Hindmarsh, Alan C., (LLNL)
    ! ***DESCRIPTION
    ! 
    !   Subroutines XERRWD, XSETF, XSETUN, and the function routine IXSAV,
    !   as given here, constitute a simplified version of the SLATEC error
    !   handling package.
    ! 
    !   All arguments are input arguments.
    ! 
    !   MSG    = The message (character array).
    !   NMES   = The length of MSG (number of characters).
    !   NERR   = The error number (not used).
    !   LEVEL  = The error level..
    !            0 or 1 means recoverable (control returns to caller).
    !            2 means fatal (run is aborted--see note below).
    !   NI     = Number of integers (0, 1, or 2) to be printed with message.
    !   I1,I2  = Integers to be printed, depending on NI.
    !   NR     = Number of reals (0, 1, or 2) to be printed with message.
    !   R1,R2  = Reals to be printed, depending on NR.
    ! 
    !   Note..  this routine is machine-dependent and specialized for use
    !   in limited context, in the following ways..
    !   1. The argument MSG is assumed to be of type CHARACTER, and
    !      the message is printed with a format of (1X,A).
    !   2. The message is assumed to take only one line.
    !      Multi-line messages are generated by repeated calls.
    !   3. If LEVEL = 2, control passes to the statement   STOP
    !      to abort the run.  This statement may be machine-dependent.
    !   4. R1 and R2 are assumed to be in double precision and are printed
    !      in D21.13 format.
    ! 
    ! ***ROUTINES CALLED  IXSAV
    ! ***REVISION HISTORY  (YYMMDD)
    !    920831  DATE WRITTEN
    !    921118  Replaced MFLGSV/LUNSAV by IXSAV. (ACH)
    !    930329  Modified prologue to SLATEC format. (FNF)
    !    930407  Changed MSG from CHARACTER*1 array to variable. (FNF)
    !    930922  Minor cosmetic change. (FNF)
    ! ***END PROLOGUE  XERRWD
    ! 
    ! *Internal Notes:
    ! 
    !  For a different default logical unit number, IXSAV (or a subsidiary
    !  routine that it calls) will need to be modified.
    !  For a different run-abort command, change the statement following
    !  statement 100 at the end.
    ! -----------------------------------------------------------------------
    !  Subroutines called by XERRWD.. None
    !  Function routine called by XERRWD.. IXSAV
    ! -----------------------------------------------------------------------
    ! **End
    real(rt) :: R1, R2
    integer    :: NMES, NERR, LEVEL, NI, I1, I2, NR
    integer    :: LUNIT, MESFLG
    character (len=80) :: MSG
    character (len=50), parameter :: fmt10 = "(1X,A)"
    character (len=50), parameter :: fmt20 = "(6X,'In above message,  I1 =',I10)"
    character (len=50), parameter :: fmt30 = "(6X,'In above message,  I1 =',I10,3X,'I2 =',I10)"
    character (len=50), parameter :: fmt40 = "(6X,'In above message,  R1 =',D21.13)"
    character (len=50), parameter :: fmt50 = "(6X,'In above,  R1 =',D21.13,3X,'R2 =',D21.13)"

    LUNIT = IXSAV(1, 0, .FALSE.)
    MESFLG = IXSAV(2, 0, .FALSE.)
    if (MESFLG .NE. 0) then
       WRITE (LUNIT,fmt10)  MSG
       select case(NI)
       case(1)
          WRITE (LUNIT, fmt20) I1
       case(2)
          WRITE (LUNIT, fmt30) I1,I2
       end select
       select case(NR)
       case(1)
          WRITE (LUNIT, fmt40) R1
       case(2)
          WRITE (LUNIT, fmt50) R1,R2
       end select
    end if
    IF (LEVEL .NE. 2) RETURN
    STOP
  end subroutine xerrwd
  
end module dvode_output_module
