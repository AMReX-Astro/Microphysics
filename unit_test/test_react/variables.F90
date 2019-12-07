!
! A module to provide integer indices into the various storage arrays
! for accessing the different variables by name.
!
module variables

  use amrex_fort_module, only : rt => amrex_real

  use network, only: nspec, spec_names

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, parameter :: MAX_NAME_LEN=20

  type plot_t
     integer :: irho = -1
     integer :: itemp = -1
     integer :: ispec = -1
     integer :: ispec_old = -1
     integer :: irodot = -1
     integer :: irho_Hnuc = -1

     integer :: n_plot_comps = 0

     character(len=MAX_NAME_LEN), allocatable :: names(:)

   contains
     procedure :: next_index => get_next_plot_index

  end type plot_t

  type(plot_t), allocatable :: p

#if defined(AMREX_USE_CUDA)
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

    ! variable information
    allocate(p)

    p % irho      = p % next_index(1)
    p % itemp     = p % next_index(1)
    p % ispec     = p % next_index(nspec)
    p % ispec_old = p % next_index(nspec)
    p % irodot    = p % next_index(nspec)
    p % irho_Hnuc = p % next_index(1)

    allocate(p%names(p%n_plot_comps))

    p % names(p % irho) = "density"
    p % names(p % itemp) = "temperature"
    do n = 0, nspec-1
       p % names(p % ispec + n) = "Xnew_" // adjustl(trim(spec_names(n+1)))
       p % names(p % ispec_old + n) = "Xold_" // adjustl(trim(spec_names(n+1)))
       p % names(p % irodot + n) = "wdot_" // adjustl(trim(spec_names(n+1)))
    enddo
    p % names(p % irho_Hnuc) = "rho_Hnuc"

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
    integer :: len

    allocate(fstring)

    ! C++ is 0-based, so add 1 to the idx
    fstring = p % names(idx+1)
    len = len_trim(fstring)
    fstring(len+1:len+1) = c_null_char

    cstring = c_loc(fstring)

  end subroutine get_var_name

  subroutine finalize_variables()

    deallocate(p%names)

  end subroutine finalize_variables

end module variables

