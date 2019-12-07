!
! A module to provide integer indices into the various storage arrays
! for accessing the different variables by name.
!
module variables

  use amrex_fort_module, only : rt => amrex_real

  use network, only: nspec, spec_names

  use actual_conductivity_module, only : cond_name

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, parameter :: MAX_NAME_LEN=20

  type plot_t
     integer :: irho = -1
     integer :: itemp = -1
     integer :: ih = -1
     integer :: ie = -1
     integer :: ip = -1
     integer :: is = -1
     integer :: iconductivity = -1
     integer :: ispec = -1

     integer :: n_plot_comps = 0

     character(len=MAX_NAME_LEN), allocatable :: names(:)

   contains
     procedure :: next_index => get_next_plot_index

  end type plot_t

  type(plot_t), allocatable :: p

#ifdef AMREX_USE_CUDA
  attributes(managed) :: p
#endif

contains

  function get_next_plot_index(this, num) result (next)

    ! return the next starting index for a plotfile quantity, and
    ! increment the counter of plotfile quantities, n_plot_comps, by
    ! num

    use amrex_fort_module, only : rt => amrex_real
    class(plot_t), intent(inout) :: this
    integer, intent(in) :: num
    integer :: next

    next = this%n_plot_comps + 1
    this%n_plot_comps = this%n_plot_comps + num

    return
  end function get_next_plot_index

  subroutine init_variables() bind(C, name="init_variables")

    use amrex_fort_module, only : rt => amrex_real
    integer :: n

    allocate(p)

    ! variable information
    p % irho      = p % next_index(1)
    p % itemp     = p % next_index(1)
    p % ih        = p % next_index(1)
    p % ie        = p % next_index(1)
    p % ip        = p % next_index(1)
    p % is        = p % next_index(1)
    p % iconductivity   = p % next_index(1)
    p % ispec     = p % next_index(nspec)

    allocate(p%names(p%n_plot_comps))

    p % names(p % irho) = "density"
    p % names(p % itemp) = "temperature"
    p % names(p % ih) = "specific_enthalpy"
    p % names(p % ie) = "specific_energy"
    p % names(p % ip) = "pressure"
    p % names(p % is) = "specific_entropy"
    p % names(p % iconductivity)   = "conductivity"
    do n = 0, nspec-1
       p % names(p % ispec + n) = "X_" // adjustl(trim(spec_names(n+1)))
    enddo

  end subroutine init_variables

  subroutine get_ncomp(ncomp_in) bind(C, name="get_ncomp")

    use amrex_fort_module, only : rt => amrex_real
    integer, intent(inout) :: ncomp_in

    ncomp_in = p % n_plot_comps

  end subroutine get_ncomp

  subroutine get_name_len(nlen_in) bind(C, name="get_name_len")

    use amrex_fort_module, only : rt => amrex_real
    integer, intent(inout) :: nlen_in

    nlen_in = MAX_NAME_LEN

  end subroutine get_name_len

  subroutine get_var_name(cstring, idx) bind(C, name="get_var_name")

    use iso_c_binding

    use amrex_fort_module, only : rt => amrex_real
    implicit none
    type(c_ptr), intent(inout) :: cstring
    integer, intent(in) :: idx

    ! include space for the NULL termination
    character(MAX_NAME_LEN+1), pointer :: fstring
    integer :: slen

    allocate(fstring)

    ! C++ is 0-based, so add 1 to the idx
    fstring = p % names(idx+1)
    slen = len_trim(fstring)
    fstring(slen+1:slen+1) = c_null_char

    cstring = c_loc(fstring)

  end subroutine get_var_name

  subroutine get_cond_len(nlen_in) bind(C, name="get_cond_len")

    use amrex_fort_module, only : rt => amrex_real
    integer, intent(inout) :: nlen_in

    nlen_in = len(cond_name)

  end subroutine get_cond_len

  subroutine get_cond_name(cond_string) bind(C, name="get_cond_name")

    use iso_c_binding

    use amrex_fort_module, only : rt => amrex_real
    implicit none
    type(c_ptr), intent(inout) :: cond_string

    ! include space for the NULL termination
    character(len(cond_name)+1), pointer :: fstring
    integer :: slen

    allocate(fstring)

    fstring = cond_name
    slen = len_trim(fstring)
    fstring(slen+1:slen+1) = c_null_char

    cond_string = c_loc(fstring)

  end subroutine get_cond_name

end module variables
