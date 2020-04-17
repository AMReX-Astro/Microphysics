!! Uniform method of printing error/warning messages.  Error messages
!! print a message to the default output unit and then call the PARALLEL_ABORT
!! to the process (or set of processes if parallel)

module bl_error_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

  subroutine bl_error(str)
    character(len=*), intent(in), optional :: str
    if ( present(str) ) then
       write(*,*) "BOXLIB ERROR: ", str
    else
       write(*,*) "BOXLIB ERROR"
    end if
  end subroutine bl_error
end module bl_error_module
