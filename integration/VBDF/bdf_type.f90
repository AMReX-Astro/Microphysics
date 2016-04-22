module bdf_type_module

  use bl_types

  implicit none

  !
  ! bdf time-stepper
  !
  type :: bdf_ts

     integer  :: neq                        ! number of equations (degrees of freedom) per point
     integer  :: npt                        ! number of points
     integer  :: max_order                  ! maximum order (1 to 6)
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

     real(dp_t), allocatable :: rtol(:)         ! relative tolerances
     real(dp_t), allocatable :: atol(:)         ! absolute tolerances

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
     !real(dp_t) :: temp_data
     !real(dp_t) :: temp_data(3,1,0:3) !z-like, if max_order=3
     !real(dp_t) :: temp_data(0:3)     !l-like, if max_order=3
     real(dp_t) :: temp_data(2,1)
     logical  :: refactor

     real(dp_t), allocatable :: J(:,:,:)        ! Jacobian matrix
     real(dp_t), allocatable :: P(:,:,:)        ! Newton iteration matrix
     real(dp_t), allocatable :: z(:,:,:)        ! Nordsieck histroy array, indexed as (dof, p, n)
     real(dp_t), allocatable :: z0(:,:,:)       ! Nordsieck predictor array
     real(dp_t), allocatable :: h(:)            ! time steps, h = [ h_n, h_{n-1}, ..., h_{n-k} ]
     real(dp_t), allocatable :: l(:)            ! predictor/corrector update coefficients
     real(dp_t), allocatable :: shift(:)        ! scratch array to hold shifted arrays
     real(dp_t), allocatable :: upar(:,:)       ! array of user parameters (passed to
     !    user's Jacobian and f)
     real(dp_t), allocatable :: y(:,:)          ! current y
     real(dp_t), allocatable :: yd(:,:)         ! current \dot{y}
     real(dp_t), allocatable :: rhs(:,:)        ! solver rhs
     real(dp_t), allocatable :: e(:,:)          ! accumulated correction
     real(dp_t), allocatable :: e1(:,:)         ! accumulated correction, previous step
     real(dp_t), allocatable :: ewt(:,:)        ! cached error weights
     real(dp_t), allocatable :: b(:,:)          ! solver work space
     integer,    allocatable :: ipvt(:,:)         ! pivots (neq,npts)
     integer,    allocatable :: A(:,:)            ! pascal matrix

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
