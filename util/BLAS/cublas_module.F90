module cublas_module

  ! The cublas_module provides Fortran interfaces for the C cuBLAS functions
  ! Above each interface subroutine is listed the C-style argument list
  ! for the corresponding cuBLAS function.
  implicit none

#ifdef CUDA
  interface cuda_blas

     ! cublasStatus_t cublasDaxpy(cublasHandle_t handle, int n,
     ! const double          *alpha,
     ! const double          *x, int incx,
     ! double                *y, int incy)
     attributes(device) &
     subroutine cublasDaxpy(n, alpha, x, incx, y, incy) bind(C,name='cublasDaxpy')
       use iso_c_binding
       implicit none
       integer(c_int), value  :: n, incx, incy
       real(c_double), value  :: alpha
       real(c_double), pointer :: x(:), y(:)
     end subroutine cublasDaxpy

     ! cublasStatus_t cublasDcopy(cublasHandle_t handle, int n,
     ! const double          *x, int incx,
     ! double                *y, int incy)
     attributes(device) &
     subroutine cublasDcopy(n, x, incx, y, incy) bind(C,name='cublasDcopy')
       use iso_c_binding
       implicit none
       integer(c_int), value  :: n, incx, incy
       real(c_double), pointer :: x(:), y(:)
     end subroutine cublasDcopy

     ! cublasStatus_t cublasDdot (cublasHandle_t handle, int n,
     ! const double          *x, int incx,
     ! const double          *y, int incy,
     ! double          *result)
     attributes(device) &
     subroutine cublasDdot(n, x, incx, y, incy, res) bind(C,name='cublasDdot')
       use iso_c_binding
       implicit none
       integer(c_int), value  :: n, incx, incy
       real(c_double), pointer :: x(:), y(:), res(:)
     end subroutine cublasDdot

     ! cublasStatus_t cublasDgemm(cublasHandle_t handle,
     ! cublasOperation_t transa, cublasOperation_t transb,
     ! int m, int n, int k,
     ! const double          *alpha,
     ! const double          *A, int lda,
     ! const double          *B, int ldb,
     ! const double          *beta,
     ! double          *C, int ldc)
     attributes(device) &
     subroutine cublasDgemm(cta, ctb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) bind(C,name='cublasDgemm')
       use iso_c_binding
       implicit none
       character(1,c_char), value :: cta, ctb
       integer(c_int), value      :: m, n, k, lda, ldb, ldc
       real(c_double), value      :: alpha, beta
       real(c_double), pointer :: A(:,:), B(:,:), C(:,:)
     end subroutine cublasDgemm

     ! cublasStatus_t  cublasDscal(cublasHandle_t handle, int n,
     ! const double          *alpha,
     ! double          *x, int incx)
     attributes(device) &
     subroutine cublasDscal(n, alpha, x, incx) bind(C,name='cublasDscal')
       use iso_c_binding, only: c_int, c_double
       use cudafor, only: c_devptr
       implicit none
       integer(c_int), value  :: n, incx
       real(c_double), value  :: alpha
       type(c_devptr), device :: x
     end subroutine cublasDscal

     ! cublasStatus_t cublasIdamax(cublasHandle_t handle, int n,
     ! const double *x, int incx, int *result)
     attributes(device) &
     subroutine cublasIdamax(n, x, incx, res) bind(C,name='cublasIdamax')
       use iso_c_binding
       implicit none
       integer(c_int), value  :: n, incx
       integer(c_int)         :: res
       real(c_double), pointer :: x(:)
     end subroutine cublasIdamax
     
  end interface cuda_blas
#endif
end module cublas_module
