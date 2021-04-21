
module actual_eos_module

  implicit none

  character (len=64), parameter :: eos_name = "tillotson" 
 
contains

  subroutine actual_eos_init

    implicit none

  end subroutine actual_eos_init


  subroutine is_input_valid(input, valid)

    implicit none

    integer, intent(in) :: input
    logical, intent(out) :: valid

    valid = .false.

  end subroutine is_input_valid


  subroutine actual_eos(input, state)

    use amrex_error_module, only: amrex_error
    use eos_type_module, only: eos_t

    implicit none

    integer,      intent(in   ) :: input
    type (eos_t), intent(inout) :: state

    ! Not implemented in Fortran

    call amrex_error("No Fortran implementation for this EOS")

  end subroutine actual_eos

  subroutine actual_eos_finalize
    
    implicit none
  
  end subroutine actual_eos_finalize

end module actual_eos_module
