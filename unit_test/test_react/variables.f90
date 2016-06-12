!
! A module to provide integer indices into the various storage arrays
! for accessing the different variables by name.
!
module variables

  use bl_types
  use network, only: nspec

  implicit none

  type plot_t
     integer :: irho = -1
     integer :: itemp = -1
     integer :: ispec = -1
     integer :: irodot = -1
     integer :: irho_Hnuc = -1

     integer :: n_plot_comps = 0

     character(len=20), allocatable :: names(:)

   contains
     procedure :: next_index => get_next_plot_index

  end type plot_t

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

  subroutine init_variables(p)

    type(plot_t), intent(inout) :: p

    ! variable information
    p % irho      = p % next_index(1)
    p % itemp     = p % next_index(1)
    p % ispec     = p % next_index(nspec)
    p % irodot    = p % next_index(nspec)
    p % irho_Hnuc = p % next_index(1)

    allocate(p % names(p % n_plot_comps))

    p % names(p % irho) = "density"
    p % names(p % itemp) = "temperature"
    do n = 0, nspec-1
       p % names(p % ispec + n) = "X_" // adjustl(trim(spec_names(n+1)))
       p % names(p % irodot + n) = "wdot_" // adjustl(trim(spec_names(n+1)))
    enddo
    p % names(p % irho_Hnuc) = "rho_Hnuc"

  end subroutine init_plot_variables


end module variables
