module bs_type_module

  use bl_types
  use burn_type_module
  use rpar_indices

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
  end type bs_t

  !$acc declare create(nseq)

end module bs_type_module
