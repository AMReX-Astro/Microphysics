module eos_type_module

  use bl_types, only: dp_t
  use network, only: nspec, naux

  implicit none

  integer, parameter :: eos_input_rt = 1  ! rho, T are inputs
  integer, parameter :: eos_input_rh = 2  ! rho, h are inputs
  integer, parameter :: eos_input_tp = 3  ! T, p are inputs
  integer, parameter :: eos_input_rp = 4  ! rho, p are inputs
  integer, parameter :: eos_input_re = 5  ! rho, e are inputs
  integer, parameter :: eos_input_ps = 6  ! p, s are inputs
  integer, parameter :: eos_input_ph = 7  ! p, h are inputs
  integer, parameter :: eos_input_th = 8  ! T, h are inputs

  ! these are used to allow for a generic interface to the 
  ! root finding
  integer, parameter :: itemp = 1
  integer, parameter :: idens = 2
  integer, parameter :: iener = 3
  integer, parameter :: ienth = 4
  integer, parameter :: ientr = 5
  integer, parameter :: ipres = 6

  ! error codes
  integer, parameter :: ierr_general         = 1
  integer, parameter :: ierr_input           = 2
  integer, parameter :: ierr_iter_conv       = 3
  integer, parameter :: ierr_neg_e           = 4
  integer, parameter :: ierr_neg_p           = 5
  integer, parameter :: ierr_neg_h           = 6
  integer, parameter :: ierr_neg_s           = 7
  integer, parameter :: ierr_iter_var        = 8
  integer, parameter :: ierr_init            = 9
  integer, parameter :: ierr_init_xn         = 10
  integer, parameter :: ierr_out_of_bounds   = 11
  integer, parameter :: ierr_not_implemented = 12

  ! Smallest possible temperature and density permitted by the user.

  real(kind=dp_t), save :: smallt = 1.d-200
  real(kind=dp_t), save :: smalld = 1.d-200

  real(kind=dp_t), save :: smallx = 1.d-200

  ! Minimum and maximum temperature, density, and ye permitted by the EOS.

  real(kind=dp_t), save :: mintemp = 1.d-200
  real(kind=dp_t), save :: maxtemp = 1.d200
  real(kind=dp_t), save :: mindens = 1.d-200
  real(kind=dp_t), save :: maxdens = 1.d200
  real(kind=dp_t), save :: minye   = 1.d-200
  real(kind=dp_t), save :: maxye   = 1.d0 + 1.d-12

  ! A generic structure holding thermodynamic quantities and their derivatives,
  ! plus some other quantities of interest.

  ! rho      -- mass density (g/cm**3)
  ! T        -- temperature (K)
  ! xn       -- the mass fractions of the individual isotopes
  ! p        -- the pressure (dyn/cm**2)
  ! h        -- the enthalpy (erg/g)
  ! e        -- the internal energy (erg/g)
  ! s        -- the entropy (erg/g/K)
  ! c_v      -- specific heat at constant volume
  ! c_p      -- specific heat at constant pressure
  ! ne       -- number density of electrons + positrons
  ! np       -- number density of positrons only
  ! eta      -- degeneracy parameter
  ! pele     -- electron pressure + positron pressure
  ! ppos     -- position pressure only
  ! mu       -- mean molecular weight
  ! mu_e     -- mean number of nucleons per electron
  ! y_e      -- electron fraction == 1 / mu_e
  ! dPdT     -- d pressure/ d temperature
  ! dPdr     -- d pressure/ d density
  ! dedT     -- d energy/ d temperature
  ! dedr     -- d energy/ d density
  ! dsdT     -- d entropy/ d temperature
  ! dsdr     -- d entropy/ d density
  ! dhdT     -- d enthalpy/ d temperature
  ! dhdr     -- d enthalpy/ d density
  ! dPdX     -- d pressure / d xmass
  ! dhdX     -- d enthalpy / d xmass at constant pressure
  ! gam1     -- first adiabatic index (d log P/ d log rho) |_s
  ! cs       -- sound speed
  ! abar     -- average atomic number ( sum_k {X_k} ) / ( sum_k {X_k/A_k} )
  ! zbar     -- average proton number ( sum_k {Z_k X_k/ A_k} ) / ( sum_k {X_k/A_k} )
  ! dpdA     -- d pressure/ d abar
  ! dpdZ     -- d pressure/ d zbar
  ! dedA     -- d energy/ d abar
  ! dedZ     -- d energy/ d zbar

  ! Initialize the main quantities to an unphysical number
  ! so that we know if the user forgot to initialize them
  ! when calling the EOS in a particular mode.

  real(kind=dp_t), parameter :: init_num  = -1.0d200
  real(kind=dp_t), parameter :: init_test = -1.0d199

  type :: eos_t

    real(kind=dp_t) :: rho         = init_num
    real(kind=dp_t) :: T           = init_num
    real(kind=dp_t) :: p           = init_num
    real(kind=dp_t) :: e           = init_num
    real(kind=dp_t) :: h           = init_num
    real(kind=dp_t) :: s           = init_num
    real(kind=dp_t) :: dpdT        != init_num
    real(kind=dp_t) :: dpdr        != init_num
    real(kind=dp_t) :: dedT        != init_num
    real(kind=dp_t) :: dedr        != init_num
    real(kind=dp_t) :: dhdT        != init_num
    real(kind=dp_t) :: dhdr        != init_num
    real(kind=dp_t) :: dsdT        != init_num
    real(kind=dp_t) :: dsdr        != init_num
    real(kind=dp_t) :: dpde        != init_num
    real(kind=dp_t) :: dpdr_e      != init_num

    real(kind=dp_t) :: xn(nspec)   = init_num
    real(kind=dp_t) :: aux(naux)   = init_num
    real(kind=dp_t) :: cv          != init_num
    real(kind=dp_t) :: cp          != init_num
    real(kind=dp_t) :: xne         != init_num
    real(kind=dp_t) :: xnp         != init_num
    real(kind=dp_t) :: eta         != init_num
    real(kind=dp_t) :: pele        != init_num
    real(kind=dp_t) :: ppos        != init_num
    real(kind=dp_t) :: mu          != init_num
    real(kind=dp_t) :: mu_e        != init_num
    real(kind=dp_t) :: y_e         != init_num
    real(kind=dp_t) :: dedX(nspec) != init_num
    real(kind=dp_t) :: dpdX(nspec) != init_num
    real(kind=dp_t) :: dhdX(nspec) != init_num
    real(kind=dp_t) :: gam1        != init_num
    real(kind=dp_t) :: cs          != init_num

    real(kind=dp_t) :: abar        != init_num
    real(kind=dp_t) :: zbar        != init_num
    real(kind=dp_t) :: dpdA        != init_num

    real(kind=dp_t) :: dpdZ        != init_num
    real(kind=dp_t) :: dedA        != init_num
    real(kind=dp_t) :: dedZ        != init_num

    real(kind=dp_t) :: smallt
    real(kind=dp_t) :: smalld

    logical :: reset                != .false.
    logical :: check_small          != .true.
    logical :: check_inputs         != .true.

  end type eos_t

  !$acc declare create(smallt, smalld, smallx)
  !$acc declare create(mintemp, maxtemp, mindens, maxdens, minye, maxye)

contains

  ! Given a set of mass fractions, calculate quantities that depend
  ! on the composition like abar and zbar.

  subroutine composition(state)

    !$acc routine seq

    use bl_constants_module, only: ONE
    use network, only: aion, zion

    implicit none

    type (eos_t), intent(inout) :: state

    ! Calculate abar, the mean nucleon number,
    ! zbar, the mean proton number,
    ! mu, the mean molecular weight,
    ! mu_e, the mean number of nucleons per electron, and
    ! y_e, the electron fraction.

    state % mu_e = ONE / (sum(state % xn(:) * zion(:) / aion(:)))
    state % y_e = ONE / state % mu_e

    state % abar = ONE / (sum(state % xn(:) / aion(:)))
    state % zbar = state % abar / state % mu_e

  end subroutine composition



  ! Compute thermodynamic derivatives with respect to xn(:)

  subroutine composition_derivatives(state)

    !$acc routine seq

    use bl_constants_module, only: ZERO
    use network, only: aion, zion

    implicit none

    type (eos_t), intent(inout) :: state

    state % dpdX(:) = state % dpdA * (state % abar/aion(:)) * (aion(:) - state % abar) &
                    + state % dpdZ * (state % abar/aion(:))   &
                                   * (zion(:) - state % zbar)

    state % dEdX(:) = state % dedA * (state % abar/aion(:))   &
                                   * (aion(:) - state % abar) &
                    + state % dedZ * (state % abar/aion(:))   &
                                   * (zion(:) - state % zbar)

    if (state % dPdr .ne. ZERO) then

       state % dhdX(:) = state % dedX(:) &
                       + (state % p / state % rho**2 - state % dedr) &
                       *  state % dPdX(:) / state % dPdr

    endif

  end subroutine composition_derivatives



  ! Normalize the mass fractions: they must be individually positive
  ! and less than one, and they must all sum to unity.

  subroutine normalize_abundances(state)

    !$acc routine seq

    use bl_constants_module, only: ONE
    use extern_probin_module, only: small_x

    implicit none

    type (eos_t), intent(inout) :: state

    state % xn(:) = max(small_x, min(ONE, state % xn(:)))

    state % xn(:) = state % xn(:) / sum(state % xn(:))

  end subroutine normalize_abundances

end module eos_type_module
