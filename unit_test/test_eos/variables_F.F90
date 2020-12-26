!
! A module to provide integer indices into the various storage arrays
! for accessing the different variables by name.
!
module variables

  use amrex_fort_module, only : rt => amrex_real

  use network, only: nspec, spec_names

  use actual_eos_module, only : eos_name

  implicit none

  integer, parameter :: MAX_NAME_LEN=20

  type plot_t_f
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

     integer :: icv = -1
     integer :: icp = -1
     integer :: ine = -1
     integer :: inp = -1
     integer :: ieta = -1
     integer :: ipele = -1
     integer :: ippos = -1
     integer :: imu = -1
     integer :: imue = -1
     integer :: iye = -1
     integer :: idpdt = -1
     integer :: idpdr = -1
     integer :: idedt = -1
     integer :: idedr = -1
     integer :: idsdt = -1
     integer :: idsdr = -1
     integer :: idhdt = -1
     integer :: idhdr = -1
     integer :: idpdx = -1
     integer :: idedx = -1
     integer :: idhdx = -1
     integer :: igam1 = -1
     integer :: ics = -1
     integer :: iabar = -1
     integer :: izbar = -1
     integer :: idpda = -1
     integer :: idpdz = -1
     integer :: ideda = -1
     integer :: idedz = -1
     integer :: idpde = -1
     integer :: idpdre = -1

     integer :: n_plot_comps = 0

   contains
     procedure :: next_index => get_next_plot_index

  end type plot_t_f

  type(plot_t_f), allocatable :: p

#if defined(AMREX_USE_CUDA) && defined(AMREX_USE_GPU_PRAGMA)
  attributes(managed) :: p
#endif

contains

  function get_next_plot_index(this, num) result (next)

    ! return the next starting index for a plotfile quantity, and
    ! increment the counter of plotfile quantities, n_plot_comps, by
    ! num

    class(plot_t_f), intent(inout) :: this
    integer, intent(in) :: num
    integer :: next

    next = this%n_plot_comps + 1
    this%n_plot_comps = this%n_plot_comps + num

    return
  end function get_next_plot_index

  subroutine init_variables_F() bind(C, name="init_variables_F")

    integer :: n

    allocate(p)

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

    p % icv = p % next_index(1)
    p % icp = p % next_index(1)
    p % ine = p % next_index(1)
    p % inp = p % next_index(1)
    p % ieta = p % next_index(1)
    p % ipele = p % next_index(1)
    p % ippos = p % next_index(1)
    p % imu = p % next_index(1)
    p % imue = p % next_index(1)
    p % iye = p % next_index(1)
    p % idpdt = p % next_index(1)
    p % idpdr = p % next_index(1)
    p % idedt = p % next_index(1)
    p % idedr = p % next_index(1)
    p % idsdt = p % next_index(1)
    p % idsdr = p % next_index(1)
    p % idhdt = p % next_index(1)
    p % idhdr = p % next_index(1)
    p % idpdx = p % next_index(nspec)
    p % idedx = p % next_index(nspec)
    p % idhdx = p % next_index(nspec)
    p % igam1 = p % next_index(1)
    p % ics = p % next_index(1)
    p % iabar = p % next_index(1)
    p % izbar = p % next_index(1)
    p % idpda = p % next_index(1)
    p % idpdz = p % next_index(1)
    p % ideda = p % next_index(1)
    p % idedz = p % next_index(1)
    p % idpde = p % next_index(1)
    p % idpdre = p % next_index(1)

  end subroutine init_variables_F

  subroutine get_ncomp(ncomp_in) bind(C, name="get_ncomp")

    integer, intent(inout) :: ncomp_in

    ncomp_in = p % n_plot_comps

  end subroutine get_ncomp

end module variables
