module bdf_type_module

  use bl_types
  use burn_type_module, only: neqs
  use rpar_indices, only: n_rpar_comps

  implicit none

  integer, parameter :: bdf_npt = 1         ! number of points
  integer, parameter :: bdf_max_order = 5   ! maximum order (1 to 6)
  integer :: A(0:bdf_max_order, 0:bdf_max_order) ! pascal matrix, shared by all ts
  !$acc declare create(A)

  !
  ! bdf time-stepper
  !
  type :: bdf_ts

     integer  :: npt
     integer  :: neq
     integer  :: max_order
     integer  :: max_steps                  ! maximum allowable number of steps
     integer  :: max_iters                  ! maximum allowable number of newton iterations
     integer  :: verbose                    ! verbosity level
     real(dp_t) :: dt_min                   ! minimum allowable step-size
     real(dp_t) :: eta_min                  ! minimum allowable step-size shrink factor
     real(dp_t) :: eta_max                  ! maximum allowable step-size growth factor
     real(dp_t) :: eta_thresh               ! step-size growth threshold
     integer  :: max_j_age                  ! maximum age of Jacobian
     integer  :: max_p_age                  ! maximum age of newton iteration matrix

     logical  :: debug
     integer  :: dump_unit

     real(dp_t) :: rtol(neqs)               ! relative tolerances
     real(dp_t) :: atol(neqs)               ! absolute tolerances

     ! state
     real(dp_t) :: t                        ! current time
     real(dp_t) :: t1                       ! final time
     real(dp_t) :: dt                       ! current time step
     real(dp_t) :: dt_nwt                   ! dt used when building newton iteration matrix
     integer  :: k                          ! current order
     integer  :: n                          ! current step
     integer  :: j_age                      ! age of Jacobian
     integer  :: p_age                      ! age of newton iteration matrix
     integer  :: k_age                      ! number of steps taken at current order
     real(dp_t) :: tq(-1:2)                 ! error coefficients (test quality)
     real(dp_t) :: tq2save
     logical  :: refactor

     real(dp_t) :: J(neqs,neqs,bdf_npt)               ! Jacobian matrix
     real(dp_t) :: P(neqs,neqs,bdf_npt)               ! Newton iteration matrix
     real(dp_t) :: z(neqs,bdf_npt,0:bdf_max_order)    ! Nordsieck histroy array, indexed as (dof, p, n)
     real(dp_t) :: z0(neqs,bdf_npt,0:bdf_max_order)   ! Nordsieck predictor array
     real(dp_t) :: h(0:bdf_max_order)                 ! time steps, h = [ h_n, h_{n-1}, ..., h_{n-k} ]
     real(dp_t) :: l(0:bdf_max_order)                 ! predictor/corrector update coefficients
     real(dp_t) :: shift(0:bdf_max_order)             ! scratch array to hold shifted arrays
     real(dp_t) :: upar(n_rpar_comps,bdf_npt)         ! array of user parameters (passed to
                                                      !    user's Jacobian and f)
     real(dp_t) :: y(neqs,bdf_npt)                    ! current y
     real(dp_t) :: yd(neqs,bdf_npt)                   ! current \dot{y}
     real(dp_t) :: rhs(neqs,bdf_npt)                  ! solver rhs
     real(dp_t) :: e(neqs,bdf_npt)                    ! accumulated correction
     real(dp_t) :: e1(neqs,bdf_npt)                   ! accumulated correction, previous step
     real(dp_t) :: ewt(neqs,bdf_npt)                  ! cached error weights
     real(dp_t) :: b(neqs,bdf_npt)                    ! solver work space
     integer    :: ipvt(neqs,bdf_npt)                 ! pivots (neq,npts)

     ! counters
     integer :: nfe                         ! number of function evaluations
     integer :: nje                         ! number of Jacobian evaluations
     integer :: nlu                         ! number of factorizations
     integer :: nit                         ! number of non-linear solver iterations
     integer :: nse                         ! number of non-linear solver errors
     integer :: ncse                        ! number of consecutive non-linear solver errors
     integer :: ncit                        ! number of current non-linear solver iterations
     integer :: ncdtmin                     ! number of consecutive times we tried to shrink beyond the minimum time step

  end type bdf_ts

end module bdf_type_module
