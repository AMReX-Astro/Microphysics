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
        integer :: ini56 = -1
        integer :: ic12ag = -1
        integer :: ic12ag_deboer17 = -1
        integer :: itriplealf = -1
        integer :: ic12c12 = -1
        integer :: ic12o16 = -1
        integer :: io16o16 = -1
        integer :: io16ag = -1
        integer :: ine20ag = -1
        integer :: img24ag = -1
        integer :: img24ap = -1
        integer :: ial27pg = -1
        integer :: ial27pg_old = -1
        integer :: isi28ag = -1
        integer :: isi28ap = -1
        integer :: ip31pg = -1
        integer :: is32ag = -1
        integer :: is32ap = -1
        integer :: icl35pg = -1
        integer :: iar36ag = -1
        integer :: iar36ap = -1
        integer :: ik39pg = -1
        integer :: ica40ag = -1
        integer :: ica40ap = -1
        integer :: isc43pg = -1
        integer :: iti44ag = -1
        integer :: iti44ap = -1
        integer :: iv47pg = -1
        integer :: icr48ag = -1
        integer :: icr48ap = -1
        integer :: imn51pg = -1
        integer :: ife52ag = -1
        integer :: ife52ap = -1
        integer :: ico55pg = -1
        integer :: ipp = -1
        integer :: ipng = -1
        integer :: idpg = -1
        integer :: ihe3ng = -1
        integer :: ihe3he3 = -1
        integer :: ihe3he4 = -1
        integer :: ic12pg = -1
        integer :: in14pg = -1
        integer :: in15pg = -1
        integer :: in15pa = -1
        integer :: io16pg = -1
        integer :: in14ag = -1
        integer :: ife52ng = -1
        integer :: ife53ng = -1
        integer :: ife54ng = -1
        integer :: ife54pg = -1
        integer :: ife54ap = -1
        integer :: ife55ng = -1
        integer :: ife56pg = -1
        integer :: ilanganke = -1
        integer :: iecapnuc = -1

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

        integer :: i, n

        ! number of rates, not include langanke and ecapnuc
        integer, parameter :: n_rates = 52
        integer :: first_rate_index

        allocate(p)

        ! variable information
        p % irho = p % next_index(1)
        p % itemp = p % next_index(1)
        p % ini56 = p % next_index(1)

        p % ic12ag = p % next_index(4)
        p % ic12ag_deboer17 = p % next_index(4)
        p % itriplealf = p % next_index(4)
        p % ic12c12 = p % next_index(4)
        p % ic12o16 = p % next_index(4)
        p % io16o16 = p % next_index(4)
        p % io16ag = p % next_index(4)
        p % ine20ag = p % next_index(4)
        p % img24ag = p % next_index(4)
        p % img24ap = p % next_index(4)
        p % ial27pg = p % next_index(4)
        p % ial27pg_old = p % next_index(4)
        p % isi28ag = p % next_index(4)
        p % isi28ap = p % next_index(4)
        p % ip31pg = p % next_index(4)
        p % is32ag = p % next_index(4)
        p % is32ap = p % next_index(4)
        p % icl35pg = p % next_index(4)
        p % iar36ag = p % next_index(4)
        p % iar36ap = p % next_index(4)
        p % ik39pg = p % next_index(4)
        p % ica40ag = p % next_index(4)
        p % ica40ap = p % next_index(4)
        p % isc43pg = p % next_index(4)
        p % iti44ag = p % next_index(4)
        p % iti44ap = p % next_index(4)
        p % iv47pg = p % next_index(4)
        p % icr48ag = p % next_index(4)
        p % icr48ap = p % next_index(4)
        p % imn51pg = p % next_index(4)
        p % ife52ag = p % next_index(4)
        p % ife52ap = p % next_index(4)
        p % ico55pg = p % next_index(4)
        p % ipp = p % next_index(4)
        p % ipng = p % next_index(4)
        p % idpg = p % next_index(4)
        p % ihe3ng = p % next_index(4)
        p % ihe3he3= p % next_index(4)
        p % ihe3he4 = p % next_index(4)
        p % ic12pg = p % next_index(4)
        p % in14pg = p % next_index(4)
        p % in15pg = p % next_index(4)
        p % in15pa = p % next_index(4)
        p % io16pg = p % next_index(4)
        p % in14ag = p % next_index(4)
        p % ife52ng = p % next_index(4)
        p % ife53ng = p % next_index(4)
        p % ife54ng = p % next_index(4)
        p % ife54pg = p % next_index(4)
        p % ife54ap = p % next_index(4)
        p % ife55ng = p % next_index(4)
        p % ife56pg = p % next_index(4)

        ! langanke and ecapnuc are different so not included in n_tests 
        p % ilanganke = p % next_index(2)
        p % iecapnuc = p % next_index(4)

    end subroutine init_variables_F

end module variables
