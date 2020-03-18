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
    integer, parameter :: job = 0

    integer, intent(in) :: ipvt(n)
    real(rt)        , intent(in) :: a(lda,n)
    real(rt)        , intent(inout) :: b(n)

    ! 
    !      dgesl solves the real(rt)         system
    !      a * x = b  or  trans(a) * x = b
    !      using the factors computed by dgeco or dgefa.
    ! 
    !      on entry
    ! 
    !         a       real(rt)        (lda, n)
    !                 the output from dgeco or dgefa.
    ! 
    !         lda     integer
    !                 the leading dimension of the array  a .
    ! 
    !         n       integer
    !                 the order of the matrix  a .
    ! 
    !         ipvt    integer(n)
    !                 the pivot array from dgeco or dgefa.
    ! 
    !         b       real(rt)        (n)
    !                 the right hand side array.
    ! 
    !         job     integer
    !                 = 0         to solve  a*x = b ,
    !                 = nonzero   to solve  trans(a)*x = b  where
    !                             trans(a)  is the transpose.
    ! 
    !      on return
    ! 
    !         b       the solution array  x .
    ! 
    !      error condition
    ! 
    !         a division by zero will occur if the input factor contains a
    !         zero on the diagonal.  technically this indicates singularity
    !         but it is often caused by improper arguments or improper
    !         setting of lda .  it will not occur if the subroutines are
    !         called correctly and if dgeco has set rcond .gt. 0.0
    !         or dgefa has set info .eq. 0 .
    ! 
    !      to compute  inverse(a) * c  where  c  is a matrix
    !      with  p  columns
    !            call dgeco(a,lda,n,ipvt,rcond,z)
    !            if (rcond is too small) go to ...
    !            do 10 j = 1, p
    !               call dgesl(a,lda,n,ipvt,c(1,j),0)
    !         10 continue
    ! 
    !      linpack. this version dated 08/14/78 .
    !      cleve moler, university of new mexico, argonne national lab.
    ! 
    !      internal variables

    real(rt)         t
    integer k,kb,l,nm1

    !$gpu

    nm1 = n - 1
    if (job == 0) then
       ! job = 0, solve a * x = b
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
       enddo

    else

       ! job = nonzero, solve  trans(a) * x = b
       ! first solve trans(u) * y = b
       do k = 1, n
          t = sum(a(1:k-1,k) * b(1:k-1))
          b(k) = (b(k) - t)/a(k,k)
       end do

       ! now solve trans(l) * x = y
       if (nm1 >= 1) then
          do kb = 1, nm1
             k = n - kb
             b(k) = b(k) + sum(a(k+1:n,k) * b(k+1:n))
             l = ipvt(k)
             if (l /= k) then
                t = b(l)
                b(l) = b(k)
                b(k) = t
             end if
          end do
       end if

    end if

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
