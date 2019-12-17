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

    class(plot_t), intent(inout) :: this
    integer, intent(in) :: num
    integer :: next

    next = this%n_plot_comps + 1
    this%n_plot_comps = this%n_plot_comps + num

    return
  end function get_next_plot_index

  subroutine init_variables() bind(C, name="init_variables")

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

    allocate(p%names(p%n_plot_comps))

    p % names(p % irho) = "density"
    p % names(p % itemp) = "temperature"
    p % names(p % ih) = "specific_enthalpy"
    p % names(p % ie) = "specific_energy"
    p % names(p % ip) = "pressure"
    p % names(p % is) = "specific_entropy"
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

    p % names(p % icv) = "c_v"
    p % names(p % icp) = "c_p"
    p % names(p % ine) = "n_e"
    p % names(p % inp) = "n_p"
    p % names(p % ieta) = "eta"
    p % names(p % ipele) = "p_ele"
    p % names(p % ippos) = "p_pos"
    p % names(p % imu) = "mu"
    p % names(p % imue) = "mu_e"
    p % names(p % iye) = "Y_e"
    p % names(p % idpdt) = "dp_dT"
    p % names(p % idpdr) = "dp_drho"
    p % names(p % idedt) = "de_dT"
    p % names(p % idedr) = "de_drho"
    p % names(p % idsdt) = "ds_dT"
    p % names(p % idsdr) = "ds_drho"
    p % names(p % idhdt) = "dh_dT"
    p % names(p % idhdr) = "dh_drho"
    do n = 0, nspec-1
       p % names(p % idpdx + n) = "dp_dX_" // adjustl(trim(spec_names(n+1)))
       p % names(p % idedx + n) = "de_dX_" // adjustl(trim(spec_names(n+1)))
       p % names(p % idhdx + n) = "dh_dX_" // adjustl(trim(spec_names(n+1)))
    enddo
    p % names(p % igam1) = "Gamma_1"
    p % names(p % ics) = "soundspeed"
    p % names(p % iabar) = "Abar"
    p % names(p % izbar) = "Zbar"
    p % names(p % idpda) = "dp_dA"
    p % names(p % idpdz) = "dp_dZ"
    p % names(p % ideda) = "dp_dA"
    p % names(p % idedz) = "de_dZ"
    p % names(p % idpde) = "dp_de_rho"
    p % names(p % idpdre) = "dp_drho_e"

  end subroutine init_variables

  subroutine get_ncomp(ncomp_in) bind(C, name="get_ncomp")

    integer, intent(inout) :: ncomp_in

    ncomp_in = p % n_plot_comps

  end subroutine get_ncomp

  subroutine get_name_len(nlen_in) bind(C, name="get_name_len")

    integer, intent(inout) :: nlen_in

    nlen_in = MAX_NAME_LEN

  end subroutine get_name_len

  subroutine get_var_name(cstring, idx) bind(C, name="get_var_name")

    use iso_c_binding

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

  subroutine get_eos_len(nlen_in) bind(C, name="get_eos_len")

    integer, intent(inout) :: nlen_in

    nlen_in = len(eos_name)

  end subroutine get_eos_len

  subroutine get_eos_name(eos_string) bind(C, name="get_eos_name")

    use iso_c_binding

    implicit none
    type(c_ptr), intent(inout) :: eos_string

    ! include space for the NULL termination
    character(len(eos_name)+1), pointer :: fstring
    integer :: slen

    allocate(fstring)

    fstring = eos_name
    slen = len_trim(fstring)
    fstring(slen+1:slen+1) = c_null_char

    eos_string = c_loc(fstring)

  end subroutine get_eos_name

end module variables
