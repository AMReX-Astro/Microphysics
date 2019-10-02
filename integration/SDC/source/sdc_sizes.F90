module sdc_sizes_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  ! these are parameters for the SDC 4th order Radau method
  integer, parameter :: SDC_NODES = 4
  integer, parameter :: SDC_MAX_ITERATIONS = 4

  real(rt), parameter :: dt_sdc(0:SDC_NODES-1) = [0.0d0, (4.0d0 - sqrt(6.0d0))/10.0d0, &
                                   (4.0d0 + sqrt(6.0d0))/10.0d0, 1.0d0]
  real(rt), parameter :: node_weights(0:SDC_NODES-1) = [0.0d0, (16.0d0 - sqrt(6.0d0))/36.0d0, &
                                        (16.0d0 + sqrt(6.0d0))/36.0d0, 1.0d0/9.0d0]

  integer, parameter :: NEWTON_SUCCESS = 0
  integer, parameter :: SINGULAR_MATRIX = -1
  integer, parameter :: CONVERGENCE_FAILURE = -2

  type ::sdc_diag_t
     integer :: count
     integer :: retries
     integer :: newton_solver_calls
     integer :: newton_iterations

  end type sdc_diag_t

end module sdc_sizes_module
