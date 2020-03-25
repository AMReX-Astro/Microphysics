!
! A module to provide integer indices into the various storage arrays
! for accessing the different variables by name.
!
module variables

  use amrex_fort_module, only : rt => amrex_real

  use network, only: nspec, spec_names

  use actual_conductivity_module, only : cond_name

  implicit none

  integer, parameter :: MAX_NAME_LEN=20

  type plot_t
     integer :: irho = -1
     integer :: itemp = -1
     integer :: ispec = -1
     integer :: iscn_he4_he4 = -1
     integer :: iscn_he4_be8 = -1
     integer :: iscn_c12_he4 = -1
     integer :: iscn_c12_c12 = -1
     integer :: iscn_c12_o16 = -1
     integer :: iscn_o16_o16 = -1
     integer :: iscn_o16_he4 = -1
     integer :: iscn_ne20_he4 = -1
     integer :: iscn_mg24_he4 = -1
     integer :: iscn_al27_p = -1
     integer :: iscn_si28_he4 = -1
     integer :: iscn_p31_p = -1
     integer :: iscn_s32_he4 = -1
     integer :: iscn_cl35_p = -1
     integer :: iscn_ar36_he4 = -1
     integer :: iscn_k39_p = -1
     integer :: iscn_ca40_he4 = -1
     integer :: iscn_sc43_p = -1
     integer :: iscn_ti44_he4 = -1
     integer :: iscn_v47_p = -1
     integer :: iscn_cr48_he4 = -1
     integer :: iscn_mn51_p = -1
     integer :: iscn_fe52_he4 = -1
     integer :: iscn_co55_p = -1
     integer :: iscn_fe54_p = -1
     integer :: iscn_fe54_he4 = -1
     integer :: iscn_fe56_p = -1
     integer :: iscn_d_p = -1
     integer :: iscn_p_p = -1
     integer :: iscn_he3_he3 = -1
     integer :: iscn_he3_he4 = -1
     integer :: iscn_c12_p = -1
     integer :: iscn_n14_p = -1
     integer :: iscn_o16_p = -1
     integer :: iscn_n14_he4 = -1

     integer :: iscn_he4_he4_dt = -1
     integer :: iscn_he4_be8_dt = -1
     integer :: iscn_c12_he4_dt = -1
     integer :: iscn_c12_c12_dt = -1
     integer :: iscn_c12_o16_dt = -1
     integer :: iscn_o16_o16_dt = -1
     integer :: iscn_o16_he4_dt = -1
     integer :: iscn_ne20_he4_dt = -1
     integer :: iscn_mg24_he4_dt = -1
     integer :: iscn_al27_p_dt = -1
     integer :: iscn_si28_he4_dt = -1
     integer :: iscn_p31_p_dt = -1
     integer :: iscn_s32_he4_dt = -1
     integer :: iscn_cl35_p_dt = -1
     integer :: iscn_ar36_he4_dt = -1
     integer :: iscn_k39_p_dt = -1
     integer :: iscn_ca40_he4_dt = -1
     integer :: iscn_sc43_p_dt = -1
     integer :: iscn_ti44_he4_dt = -1
     integer :: iscn_v47_p_dt = -1
     integer :: iscn_cr48_he4_dt = -1
     integer :: iscn_mn51_p_dt = -1
     integer :: iscn_fe52_he4_dt = -1
     integer :: iscn_co55_p_dt = -1
     integer :: iscn_fe54_p_dt = -1
     integer :: iscn_fe54_he4_dt = -1
     integer :: iscn_fe56_p_dt = -1
     integer :: iscn_d_p_dt = -1
     integer :: iscn_p_p_dt = -1
     integer :: iscn_he3_he3_dt = -1
     integer :: iscn_he3_he4_dt = -1
     integer :: iscn_c12_p_dt = -1
     integer :: iscn_n14_p_dt = -1
     integer :: iscn_o16_p_dt = -1
     integer :: iscn_n14_he4_dt = -1

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
    p % ispec     = p % next_index(nspec)

    p % iscn_he4_he4 = p % next_index(1)
    p % iscn_he4_be8 = p % next_index(1)
    p % iscn_c12_he4 = p % next_index(1)
    p % iscn_c12_c12 = p % next_index(1)
    p % iscn_c12_o16 = p % next_index(1)
    p % iscn_o16_o16 = p % next_index(1)
    p % iscn_o16_he4 = p % next_index(1)
    p % iscn_ne20_he4 = p % next_index(1)
    p % iscn_mg24_he4 = p % next_index(1)
    p % iscn_al27_p = p % next_index(1)
    p % iscn_si28_he4 = p % next_index(1)
    p % iscn_p31_p = p % next_index(1)
    p % iscn_s32_he4 = p % next_index(1)
    p % iscn_cl35_p = p % next_index(1)
    p % iscn_ar36_he4 = p % next_index(1)
    p % iscn_k39_p = p % next_index(1)
    p % iscn_ca40_he4 = p % next_index(1)
    p % iscn_sc43_p = p % next_index(1)
    p % iscn_ti44_he4 = p % next_index(1)
    p % iscn_v47_p = p % next_index(1)
    p % iscn_cr48_he4 = p % next_index(1)
    p % iscn_mn51_p = p % next_index(1)
    p % iscn_fe52_he4 = p % next_index(1)
    p % iscn_co55_p = p % next_index(1)
    p % iscn_fe54_p = p % next_index(1)
    p % iscn_fe54_he4 = p % next_index(1)
    p % iscn_fe56_p = p % next_index(1)
    p % iscn_d_p = p % next_index(1)
    p % iscn_p_p = p % next_index(1)
    p % iscn_he3_he3 = p % next_index(1)
    p % iscn_he3_he4 = p % next_index(1)
    p % iscn_c12_p = p % next_index(1)
    p % iscn_n14_p = p % next_index(1)
    p % iscn_o16_p = p % next_index(1)
    p % iscn_n14_he4 = p % next_index(1)

    p % iscn_he4_he4_dt = p % next_index(1)
    p % iscn_he4_be8_dt = p % next_index(1)
    p % iscn_c12_he4_dt = p % next_index(1)
    p % iscn_c12_c12_dt = p % next_index(1)
    p % iscn_c12_o16_dt = p % next_index(1)
    p % iscn_o16_o16_dt = p % next_index(1)
    p % iscn_o16_he4_dt = p % next_index(1)
    p % iscn_ne20_he4_dt = p % next_index(1)
    p % iscn_mg24_he4_dt = p % next_index(1)
    p % iscn_al27_p_dt = p % next_index(1)
    p % iscn_si28_he4_dt = p % next_index(1)
    p % iscn_p31_p_dt = p % next_index(1)
    p % iscn_s32_he4_dt = p % next_index(1)
    p % iscn_cl35_p_dt = p % next_index(1)
    p % iscn_ar36_he4_dt = p % next_index(1)
    p % iscn_k39_p_dt = p % next_index(1)
    p % iscn_ca40_he4_dt = p % next_index(1)
    p % iscn_sc43_p_dt = p % next_index(1)
    p % iscn_ti44_he4_dt = p % next_index(1)
    p % iscn_v47_p_dt = p % next_index(1)
    p % iscn_cr48_he4_dt = p % next_index(1)
    p % iscn_mn51_p_dt = p % next_index(1)
    p % iscn_fe52_he4_dt = p % next_index(1)
    p % iscn_co55_p_dt = p % next_index(1)
    p % iscn_fe54_p_dt = p % next_index(1)
    p % iscn_fe54_he4_dt = p % next_index(1)
    p % iscn_fe56_p_dt = p % next_index(1)
    p % iscn_d_p_dt = p % next_index(1)
    p % iscn_p_p_dt = p % next_index(1)
    p % iscn_he3_he3_dt = p % next_index(1)
    p % iscn_he3_he4_dt = p % next_index(1)
    p % iscn_c12_p_dt = p % next_index(1)
    p % iscn_n14_p_dt = p % next_index(1)
    p % iscn_o16_p_dt = p % next_index(1)
    p % iscn_n14_he4_dt = p % next_index(1)


    allocate(p%names(p%n_plot_comps))

    p % names(p % irho) = "density"
    p % names(p % itemp) = "temperature"
    do n = 0, nspec-1
       p % names(p % ispec + n) = "X_" // adjustl(trim(spec_names(n+1)))
    enddo
    p % names(p % iscn_he4_he4) = "scn_he4_he4"
    p % names(p % iscn_he4_be8) = "scn_he4_be8"
    p % names(p % iscn_c12_he4) = "scn_c12_he4"
    p % names(p % iscn_c12_c12) = "scn_c12_c12"
    p % names(p % iscn_c12_o16) = "scn_c12_o16"
    p % names(p % iscn_o16_o16) = "scn_o16_o16"
    p % names(p % iscn_o16_he4) = "scn_o16_he4"
    p % names(p % iscn_ne20_he4) = "scn_ne20_he4"
    p % names(p % iscn_mg24_he4) = "scn_mg24_he4"
    p % names(p % iscn_al27_p) = "scn_al27_p"
    p % names(p % iscn_si28_he4) = "scn_si28_he4"
    p % names(p % iscn_p31_p) = "scn_p31_p"
    p % names(p % iscn_s32_he4) = "scn_s32_he4"
    p % names(p % iscn_cl35_p) = "scn_cl35_p"
    p % names(p % iscn_ar36_he4) = "scn_ar36_he4"
    p % names(p % iscn_k39_p) = "scn_k39_p"
    p % names(p % iscn_ca40_he4) = "scn_ca40_he4"
    p % names(p % iscn_sc43_p) = "scn_sc43_p"
    p % names(p % iscn_ti44_he4) = "scn_ti44_he4"
    p % names(p % iscn_v47_p) = "scn_v47_p"
    p % names(p % iscn_cr48_he4) = "scn_cr48_he4"
    p % names(p % iscn_mn51_p) = "scn_mn51_p"
    p % names(p % iscn_fe52_he4) = "scn_fe52_he4"
    p % names(p % iscn_co55_p) = "scn_co55_p"
    p % names(p % iscn_fe54_p) = "scn_fe54_p"
    p % names(p % iscn_fe54_he4) = "scn_fe54_he4"
    p % names(p % iscn_fe56_p) = "scn_fe56_p"
    p % names(p % iscn_d_p) = "scn_d_p"
    p % names(p % iscn_p_p) = "scn_p_p"
    p % names(p % iscn_he3_he3) = "scn_he3_he3"
    p % names(p % iscn_he3_he4) = "scn_he3_he4"
    p % names(p % iscn_c12_p) = "scn_c12_p"
    p % names(p % iscn_n14_p) = "scn_n14_p"
    p % names(p % iscn_o16_p) = "scn_o16_p"
    p % names(p % iscn_n14_he4) = "scn_n14_he4"

    p % names(p % iscn_he4_he4_dt) = "scn_he4_he4_dt"
    p % names(p % iscn_he4_be8_dt) = "scn_he4_be8_dt"
    p % names(p % iscn_c12_he4_dt) = "scn_c12_he4_dt"
    p % names(p % iscn_c12_c12_dt) = "scn_c12_c12_dt"
    p % names(p % iscn_c12_o16_dt) = "scn_c12_o16_dt"
    p % names(p % iscn_o16_o16_dt) = "scn_o16_o16_dt"
    p % names(p % iscn_o16_he4_dt) = "scn_o16_he4_dt"
    p % names(p % iscn_ne20_he4_dt) = "scn_ne20_he4_dt"
    p % names(p % iscn_mg24_he4_dt) = "scn_mg24_he4_dt"
    p % names(p % iscn_al27_p_dt) = "scn_al27_p_dt"
    p % names(p % iscn_si28_he4_dt) = "scn_si28_he4_dt"
    p % names(p % iscn_p31_p_dt) = "scn_p31_p_dt"
    p % names(p % iscn_s32_he4_dt) = "scn_s32_he4_dt"
    p % names(p % iscn_cl35_p_dt) = "scn_cl35_p_dt"
    p % names(p % iscn_ar36_he4_dt) = "scn_ar36_he4_dt"
    p % names(p % iscn_k39_p_dt) = "scn_k39_p_dt"
    p % names(p % iscn_ca40_he4_dt) = "scn_ca40_he4_dt"
    p % names(p % iscn_sc43_p_dt) = "scn_sc43_p_dt"
    p % names(p % iscn_ti44_he4_dt) = "scn_ti44_he4_dt"
    p % names(p % iscn_v47_p_dt) = "scn_v47_p_dt"
    p % names(p % iscn_cr48_he4_dt) = "scn_cr48_he4_dt"
    p % names(p % iscn_mn51_p_dt) = "scn_mn51_p_dt"
    p % names(p % iscn_fe52_he4_dt) = "scn_fe52_he4_dt"
    p % names(p % iscn_co55_p_dt) = "scn_co55_p_dt"
    p % names(p % iscn_fe54_p_dt) = "scn_fe54_p_dt"
    p % names(p % iscn_fe54_he4_dt) = "scn_fe54_he4_dt"
    p % names(p % iscn_fe56_p_dt) = "scn_fe56_p_dt"
    p % names(p % iscn_d_p_dt) = "scn_d_p_dt"
    p % names(p % iscn_p_p_dt) = "scn_p_p_dt"
    p % names(p % iscn_he3_he3_dt) = "scn_he3_he3_dt"
    p % names(p % iscn_he3_he4_dt) = "scn_he3_he4_dt"
    p % names(p % iscn_c12_p_dt) = "scn_c12_p_dt"
    p % names(p % iscn_n14_p_dt) = "scn_n14_p_dt"
    p % names(p % iscn_o16_p_dt) = "scn_o16_p_dt"
    p % names(p % iscn_n14_he4_dt) = "scn_n14_he4_dt"

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

  subroutine get_cond_len(nlen_in) bind(C, name="get_cond_len")

    integer, intent(inout) :: nlen_in

    nlen_in = len(cond_name)

  end subroutine get_cond_len

  subroutine get_cond_name(cond_string) bind(C, name="get_cond_name")

    use iso_c_binding

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
