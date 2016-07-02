module bs_type_module

  use bl_types, only: dp_t
  use burn_type_module, only: neqs
  use rpar_indices, only: n_rpar_comps

  implicit none

  ! BS parameters -- see the discussion in 16.4
  integer, parameter :: KMAXX = 7
  integer :: nseq(KMAXX+1)

  ! error codes
  integer, parameter :: IERR_NONE = 0
  integer, parameter :: IERR_DT_TOO_SMALL = -100
  integer, parameter :: IERR_TOO_MANY_STEPS = -101
  integer, parameter :: IERR_DT_UNDERFLOW = -102
  integer, parameter :: IERR_NO_CONVERGENCE = -103

  integer, parameter :: IERR_LU_DECOMPOSITION_ERROR = -200

  real(kind=dp_t), parameter :: S1 = 0.25_dp_t
  real(kind=dp_t), parameter :: S2 = 0.7_dp_t

  real(kind=dp_t), parameter :: RED_BIG_FACTOR = 0.7_dp_t
  real(kind=dp_t), parameter :: RED_SMALL_FACTOR = 1.e-5_dp_t
  real(kind=dp_t), parameter :: SCALMX = 0.1_dp_t

  type bs_t
     logical :: first
     real(kind=dp_t) :: eps_old
     real(kind=dp_t) :: dt_did
     real(kind=dp_t) :: dt_next
     real(kind=dp_t) :: a(KMAXX+1)
     real(kind=dp_t) :: alpha(KMAXX, KMAXX)
     real(kind=dp_t) :: t_new
     integer :: kmax
     integer :: kopt

     real(kind=dp_t) :: y(neqs), dydt(neqs), jac(neqs, neqs)
     real(kind=dp_t) :: atol(neqs), rtol(neqs)
     real(kind=dp_t) :: upar(n_rpar_comps)
     real(kind=dp_t) :: t, dt, tmax
     integer         :: n
     integer         :: n_rhs, n_jac

     integer :: i, j, k

  end type bs_t

  !$acc declare create(nseq)

contains

  subroutine clean_state(state)

    use bl_constants_module, only: ONE
    use extern_probin_module, only: renormalize_abundances
    use actual_network, only: aion, nspec, nspec_evolve
    use integration_data, only: aionInv, temp_scale
    use burn_type_module, only: net_itemp

    implicit none

    type (bs_t) :: state

    ! Ensure that mass fractions always stay positive.

    state % y(1:nspec_evolve) = max(state % y(1:nspec_evolve) * aion(1:nspec_evolve), 1.d-30) * aionInv(1:nspec_evolve)
    state % y(1:nspec_evolve) = min(state % y(1:nspec_evolve) * aion(1:nspec_evolve), ONE) * aionInv(1:nspec_evolve)

    ! Ensure that the temperature always stays within reasonable limits.

    state % y(net_itemp) = min(1.0d11 / temp_scale, max(state % y(net_itemp), 1.0d4 / temp_scale))

  end subroutine clean_state



  subroutine renormalize_species(state)

    use actual_network, only: aion, nspec, nspec_evolve
    use rpar_indices, only: irp_nspec

    implicit none

    type (bs_t) :: state

    real(dp_t) :: nspec_sum

    nspec_sum = sum(state % y(1:nspec_evolve) * aion(1:nspec_evolve)) + &
                sum(state % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) * aion(nspec_evolve+1:nspec))

    state % y(1:nspec_evolve) = state % y(1:nspec_evolve) / nspec_sum
    state % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) = state % upar(irp_nspec:irp_nspec+nspec-nspec_evolve-1) / nspec_sum

  end subroutine renormalize_species

end module bs_type_module
