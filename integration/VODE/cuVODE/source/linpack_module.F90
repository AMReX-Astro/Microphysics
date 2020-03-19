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

    integer,  intent(in   ) :: ipvt(VODE_NEQS)
    real(rt), intent(in   ) :: a(VODE_NEQS,VODE_NEQS)
    real(rt), intent(inout) :: b(VODE_NEQS)

    real(rt) :: t
    integer  :: k, kb, l, nm1

    !$gpu

    nm1 = VODE_NEQS - 1

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
          b(k+1:VODE_NEQS) = b(k+1:VODE_NEQS) + t * a(k+1:VODE_NEQS,k)
       end do
    end if

    ! now solve u * x = y
    do kb = 1, VODE_NEQS
       k = VODE_NEQS + 1 - kb
       b(k) = b(k) / a(k,k)
       t = -b(k)
       b(1:k-1) = b(1:k-1) + t * a(1:k-1,k)
    end do

  end subroutine dgesl

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  subroutine dgefa (a, ipvt, info)

    real(rt), intent(inout) :: a(VODE_NEQS, VODE_NEQS)
    integer,  intent(inout) :: ipvt(VODE_NEQS), info

    ! dgefa factors a matrix by gaussian elimination.
    ! a is returned in the form a = l * u where
    ! l is a product of permutation and unit lower
    ! triangular matrices and u is upper triangular.

    real(rt) :: t
    integer  :: j, k, kp1, l, nm1

    !$gpu

    ! gaussian elimination with partial pivoting

    info = 0
    nm1 = VODE_NEQS - 1

    if (nm1 >= 1) then

       do k = 1, nm1

          kp1 = k + 1

          ! find l = pivot index
          l = idamax(VODE_NEQS-k+1, a(k:VODE_NEQS,k)) + k - 1
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
             a(k+1:VODE_NEQS,k) = a(k+1:VODE_NEQS,k) * t

             ! row elimination with column indexing
             do j = kp1, VODE_NEQS
                t = a(l,j)
                if (l /= k) then
                   a(l,j) = a(k,j)
                   a(k,j) = t
                end if
                a(k+1:VODE_NEQS,j) = a(k+1:VODE_NEQS,j) + t * a(k+1:VODE_NEQS,k)
             end do

          else

             info = k

          end if

       end do

    end if

    ipvt(VODE_NEQS) = VODE_NEQS
    if (a(VODE_NEQS,VODE_NEQS) .eq. 0.0e0_rt) info = VODE_NEQS

  end subroutine dgefa

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  function idamax(N, x) result(index)

    implicit none

    integer,  intent(in) :: N
    real(rt), intent(in) :: x(N)

    real(rt) :: dmax
    integer  :: i, index

    !$gpu

    index = 1

    dmax = abs(x(1))
    do i = 2, N
       if (abs(x(i)) .le. dmax) cycle
       index = i
       dmax = abs(x(i))
    end do

  end function idamax

end module linpack_module
