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

        allocate(p%names(p%n_plot_comps))

        p % names(p % irho) = "density"
        p % names(p % itemp) = "temperature"
        p % names(p % ini56) = "X(ni56)"

        p % names(p % ic12ag) = "c12ag"
        p % names(p % ic12ag_deboer17) = "c12ag_deboer17"
        p % names(p % itriplealf) = "triplealf"
        p % names(p % ic12c12) = "c12c12"
        p % names(p % ic12o16) = "c12o16"
        p % names(p % io16o16) = "o16o16"
        p % names(p % io16ag) = "o16ag"
        p % names(p % ine20ag) = "ne20ag"
        p % names(p % img24ag) = "mg24ag"
        p % names(p % img24ap) = "mg24ap"
        p % names(p % ial27pg) = "al27pg"
        p % names(p % ial27pg_old) = "al27pg_old"
        p % names(p % isi28ag) = "si28ag"
        p % names(p % isi28ap) = "si28ap"
        p % names(p % ip31pg) = "p31pg"
        p % names(p % is32ag) = "s32ag"
        p % names(p % is32ap) = "s32ap"
        p % names(p % icl35pg) = "cl35pg"
        p % names(p % iar36ag) = "ar36ag"
        p % names(p % iar36ap) = "ar36ap"
        p % names(p % ik39pg) = "k39pg"
        p % names(p % ica40ag) = "ca40ag"
        p % names(p % ica40ap) = "ca40ap"
        p % names(p % isc43pg) = "sc43pg"
        p % names(p % iti44ag) = "ti44ag"
        p % names(p % iti44ap) = "ti44ap"
        p % names(p % iv47pg) = "v47pg"
        p % names(p % icr48ag) = "cr48ag"
        p % names(p % icr48ap) = "cr48ap"
        p % names(p % imn51pg) = "mn51pg"
        p % names(p % ife52ag) = "fe52ag"
        p % names(p % ife52ap) = "fe52ap"
        p % names(p % ico55pg) = "co55pg"
        p % names(p % ipp) = "pp"
        p % names(p % ipng) = "png"
        p % names(p % idpg) = "dpg"
        p % names(p % ihe3ng) = "he3ng"
        p % names(p % ihe3he3) = "he3he3"
        p % names(p % ihe3he4) = "he3he4"
        p % names(p % ic12pg) = "c12pg"
        p % names(p % in14pg) = "n14pg"
        p % names(p % in15pg) = "n15pg"
        p % names(p % in15pa) = "n15pa"
        p % names(p % io16pg) = "o16pg"
        p % names(p % in14ag) = "n14ag"
        p % names(p % ife52ng) = "fe52ng"
        p % names(p % ife53ng) = "fe53ng"
        p % names(p % ife54ng) = "fe54ng"
        p % names(p % ife54pg) = "fe54pg"
        p % names(p % ife54ap) = "fe54ap"
        p % names(p % ife55ng) = "fe55ng"
        p % names(p % ife56pg) = "fe56pg"

        first_rate_index = p % ic12ag

        do n = 0, n_rates-1
            do i = 1, 3
                p % names(first_rate_index + 4*n + i) = trim(p % names(first_rate_index + 4*n))
            enddo

            p % names(first_rate_index + 4*n) = trim(p % names(first_rate_index + 4*n)) // trim("_fr")
            p % names(first_rate_index + 4*n + 1) = trim(p % names(first_rate_index + 4*n + 1)) // trim("_dfrdt")
            p % names(first_rate_index + 4*n + 2) = trim(p % names(first_rate_index + 4*n + 2)) // trim("_rr")
            p % names(first_rate_index + 4*n + 3) = trim(p % names(first_rate_index + 4*n + 3)) // trim("_drrdt")
        enddo

        p % names(p % ilanganke) = "langanke_rn56ec"
        p % names(p % ilanganke+1) = "langanke_sn56ec"

        p % names(p % iecapnuc) = "ecapnuc_rpen"
        p % names(p % iecapnuc+1) = "ecapnuc_rnep"
        p % names(p % iecapnuc+2) = "ecapnuc_spenc"
        p % names(p % iecapnuc+3) = "ecapnuc_snepc"

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
