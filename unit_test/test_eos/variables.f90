!
! A module to provide integer indices into the various storage arrays
! for accessing the different variables by name.
!
module variables

  use bl_types
  use network, only: nspec, spec_names

  implicit none

  type plot_t
     integer :: irho = -1
     integer :: itemp = -1
     integer :: ih = -1
     integer :: ie = -1
     integer :: ip = -1
     integer :: is = -1
     integer :: ispec = -1

     integer :: ierr_T_eos_rh = -1
     integer :: ierr_rho_eos_tp = -1
     integer :: ierr_T_eos_rp = -1
     integer :: ierr_T_eos_re = -1
     integer :: ierr_rho_eos_ps = -1
     integer :: ierr_T_eos_ps = -1
     integer :: ierr_rho_eos_ph = -1
     integer :: ierr_T_eos_ph = -1
     integer :: ierr_rho_eos_th = -1

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

    integer :: n

    ! variable information
    p % irho      = p % next_index(1)
    p % itemp     = p % next_index(1)
    p % ih        = p % next_index(1)
    p % ie        = p % next_index(1)
    p % ip        = p % next_index(1)
    p % is        = p % next_index(1)
    p % ispec     = p % next_index(nspec)

    p % ierr_T_eos_rh   = p % next_index(1)
    p % ierr_rho_eos_tp = p % next_index(1)
    p % ierr_T_eos_rp = p % next_index(1)
    p % ierr_T_eos_re = p % next_index(1)
    p % ierr_rho_eos_ps = p % next_index(1)
    p % ierr_T_eos_ps = p % next_index(1)
    p % ierr_rho_eos_ph = p % next_index(1)
    p % ierr_T_eos_ph = p % next_index(1)
    p % ierr_rho_eos_th = p % next_index(1)

    allocate(p%names(p%n_plot_comps))

    p % names(p % irho) = "density"
    p % names(p % itemp) = "temperature"
    p % names(p % ih) = "specific_enthalpy"
    p % names(p % ie) = "specific_energy"
    p % names(p % ip) = "pressure"
    p % names(p % ip) = "specific_entropy"
    do n = 0, nspec-1
       p % names(p % ispec + n) = "X_" // adjustl(trim(spec_names(n+1)))
    enddo

    p % names(p % ierr_T_eos_rh)   = "err_T_eos_rh"
    p % names(p % ierr_rho_eos_tp) = "err_rho_eos_tp"
    p % names(p % ierr_T_eos_rp)   = "err_T_eos_rp"
    p % names(p % ierr_T_eos_re)   = "err_T_eos_re"
    p % names(p % ierr_rho_eos_ps) = "err_rho_eos_ps"
    p % names(p % ierr_T_eos_ps)   = "err_T_eos_ps"
    p % names(p % ierr_rho_eos_ph) = "err_rho_eos_ph"
    p % names(p % ierr_T_eos_ph)   = "err_T_eos_ph"
    p % names(p % ierr_rho_eos_th) = "err_rho_eos_th"

  end subroutine init_variables

end module variables
