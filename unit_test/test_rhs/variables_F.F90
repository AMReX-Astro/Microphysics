!
! A module to provide integer indices into the various storage arrays
! for accessing the different variables by name.
!
module variables

  use amrex_fort_module, only : rt => amrex_real

  use network, only: nspec, spec_names
  use burn_type_module, only: neqs

  implicit none

  integer, parameter :: MAX_NAME_LEN=50

  type plot_t
     integer :: irho = -1
     integer :: itemp = -1
     integer :: ispec = -1
     integer :: ispec_old = -1
     integer :: itemp_dot = -1
     integer :: ienuc_dot = -1
     integer :: ijac = -1

     integer :: n_plot_comps = 0

   contains
     procedure :: next_index => get_next_plot_index

  end type plot_t

  type(plot_t), allocatable :: p

#if defined(AMREX_USE_CUDA) && defined(AMREX_USE_GPU_PRAGMA)
  attributes(managed) :: p
#endif

contains

  function get_next_plot_index(this, num) result (next)

    ! return the next starting index for a plotfile quantity, and
    ! increment the counter of plotfile quantities, n_plot_comps, by
    ! num

    class(plot_t), intent(inout) :: this
    integer, intent(in) :: num
    integer :: next

    next = this%n_plot_comps + 1
    this%n_plot_comps = this%n_plot_comps + num

    return
  end function get_next_plot_index

  subroutine init_variables_F() bind(C, name="init_variables_F")

    integer :: n, i, j

    ! variable information
    allocate(p)

    p % irho      = p % next_index(1)
    p % itemp     = p % next_index(1)
    p % ispec     = p % next_index(nspec)
    p % ispec_old = p % next_index(nspec)
    p % itemp_dot = p % next_index(1)
    p % ienuc_dot = p % next_index(1)
    p % ijac      = p % next_index(neqs * neqs)

  end subroutine init_variables_F

  subroutine finalize_variables()

  end subroutine finalize_variables

end module variables

