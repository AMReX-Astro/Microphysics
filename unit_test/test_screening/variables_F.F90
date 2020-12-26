!
! A module to provide integer indices into the various storage arrays
! for accessing the different variables by name.
!
module variables

  use amrex_fort_module, only : rt => amrex_real

  use network, only: nspec, spec_names

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

  end subroutine init_variables_F

end module variables
