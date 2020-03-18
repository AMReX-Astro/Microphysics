module linpack_module

  use cuvode_parameters_module, only: VODE_NEQS

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif  
  subroutine dgesl(a, ipvt, b)

    implicit none

    integer, parameter :: lda = VODE_NEQS
    integer, parameter :: n = VODE_NEQS

    integer, intent(in) :: ipvt(n)
    real(rt)        , intent(in) :: a(lda,n)
    real(rt)        , intent(inout) :: b(n)

    real(rt) :: t
    integer  :: k, kb, l, nm1

    !$gpu

    nm1 = n - 1

    ! solve a * x = b
    ! first solve l * y = b
    if (nm1 >= 1) then
       do k = 1, nm1
          l = ipvt(k)
          t = b(l)
          if (l /= k) then
             b(l) = b(k)
             b(k) = t
          end if
          b(k+1:n) = b(k+1:n) + t * a(k+1:n,k)
       end do
    end if

    ! now solve u * x = y
    do kb = 1, n
       k = n + 1 - kb
       b(k) = b(k)/a(k,k)
       t = -b(k)
       b(1:k-1) = b(1:k-1) + t * a(1:k-1,k)
    end do

  end subroutine dgesl

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  subroutine dgefa (a, ipvt, info)

    integer, parameter :: lda = VODE_NEQS
    integer, parameter :: n = VODE_NEQS
    integer ipvt(n), info
    real(rt)         a(lda, n)
    ! 
    !      dgefa factors a real(rt)         matrix by gaussian elimination.
    ! 
    !      dgefa is usually called by dgeco, but it can be called
    !      directly with a saving in time if  rcond  is not needed.
    !      (time for dgeco) = (1 + 9/n)*(time for dgefa) .
    ! 
    !      on entry
    ! 
    !         a       real(rt)        (lda, n)
    !                 the matrix to be factored.
    ! 
    !         lda     integer
    !                 the leading dimension of the array  a .
    ! 
    !         n       integer
    !                 the order of the matrix  a .
    ! 
    !      on return
    ! 
    !         a       an upper triangular matrix and the multipliers
    !                 which were used to obtain it.
    !                 the factorization can be written  a = l*u  where
    !                 l  is a product of permutation and unit lower
    !                 triangular matrices and  u  is upper triangular.
    ! 
    !         ipvt    integer(n)
    !                 an integer array of pivot indices.
    ! 
    !         info    integer
    !                 = 0  normal value.
    !                 = k  if  u(k,k) .eq. 0.0 .  this is not an error
    !                      condition for this subroutine, but it does
    !                      indicate that dgesl or dgedi will divide by zero
    !                      if called.  use  rcond  in dgeco for a reliable
    !                      indication of singularity.
    ! 
    !      linpack. this version dated 08/14/78 .
    !      cleve moler, university of new mexico, argonne national lab.
    ! 
    !      subroutines and functions
    ! 
    !      blas vdaxpy
    ! 
    !      internal variables
    ! 
    real(rt)         t
    integer j,k,kp1,l,nm1

    !$gpu

    ! gaussian elimination with partial pivoting

    info = 0
    nm1 = n - 1

    if (nm1 >= 1) then

       do k = 1, nm1

          kp1 = k + 1

          ! find l = pivot index
          l = idamax(n-k+1,a(k:n,k)) + k - 1
          ipvt(k) = l

          ! zero pivot implies this column already triangularized
          if (a(l,k) /= 0.0e0_rt) then

             ! interchange if necessary
             if (l /= k) then
                t = a(l,k)
                a(l,k) = a(k,k)
                a(k,k) = t
             end if

             ! compute multipliers
             t = -1.0e0_rt / a(k,k)
             a(k+1:n,k) = a(k+1:n,k) * t

             ! row elimination with column indexing
             do j = kp1, n
                t = a(l,j)
                if (l /= k) then
                   a(l,j) = a(k,j)
                   a(k,j) = t
                end if
                a(k+1:n,j) = a(k+1:n,j) + t * a(k+1:n,k)
             end do

          else

             info = k

          end if

       end do

    end if

    ipvt(n) = n
    if (a(n,n) .eq. 0.0e0_rt) info = n

  end subroutine dgefa

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  function idamax(N, x) result(index)

    implicit none

    integer, intent(in) :: N
    real(rt)        , intent(in) :: x(N)

    real(rt)         :: dmax
    integer :: i, index

    !$gpu

    index = 1

    dmax = abs(x(1))
    DO i = 2, N
       IF (abs(X(i)) .le. dmax) cycle
       index = i
       dmax = abs(x(i))
    end do

  end function idamax

end module linpack_module
