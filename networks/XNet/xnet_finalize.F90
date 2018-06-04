subroutine bn_xnetFinalize()
  !use Driver_interface, ONLY : Driver_getComm, Driver_getMype, Driver_getNumProcs
  use bn_xnetData, ONLY : xnet_writeTimers
  use timers, ONLY : timer_burner, timer_xnet, timer_tstep, timer_nraph, timer_deriv, &
     timer_jacob, timer_solve, timer_csect, timer_scrn
  use xnet_interface, ONLY : jacobian_finalize

  use mpi
  !$ use omp_lib

  implicit none

#include "constants.h"

  integer, parameter :: ntimers = 9

  real(8), allocatable, dimension(:,:,:) :: t_ALL
  real(8), allocatable, dimension(:,:) :: t_ALL_omp_min, t_ALL_omp_max, t_ALL_omp_avg
  real(8), allocatable, dimension(:) :: t_ALL_min, t_ALL_max, t_ALL_avg
  real(8), allocatable, dimension(:,:) :: t_MPI

  integer :: ierr, i, j, k
  integer :: bn_globalMe, bn_globalNumProcs, bn_globalComm, mythread, nthreads

  character(2) :: cstrLen, cntimers, cstrLenM9
  character(12) :: headerFormat
  character(16) :: allFormat
  character(19) :: ompFormat
  character(19) :: mpiFormat

  integer, parameter :: strLen = 14
  character(strLen), parameter :: tHeader(ntimers+2) = [ character(strLen) :: &
    "           MPI", &
    "           OMP", &
    "      t_burner", &
    "        t_xnet", &
    "       t_tstep", &
    "       t_nraph", &
    "       t_deriv", &
    "       t_jacob", &
    "       t_solve", &
    "       t_csect", &
    "        t_scrn" ]

  if (xnet_writeTimers) then

    !TODO: replace these with Microphysics routines
    !call Driver_getMype(GLOBAL_COMM,bn_globalMe)
    !call Driver_getNumProcs(GLOBAL_COMM,bn_globalNumProcs)
    !call Driver_getComm(GLOBAL_COMM,bn_globalComm)

     write(cstrLen,'(i2)') strLen
     write(cntimers,'(i2)') ntimers
     write(cstrLenM9,'(i2)') max(strLen-9,1)
     headerFormat = '(2A'//trim(adjustl(cstrLen))//','//trim(adjustl(cntimers))//'A'//trim(adjustl(cstrLen))//')'
     allFormat = '(2I'//trim(adjustl(cstrLen))//','//trim(adjustl(cntimers))//'ES'//trim(adjustl(cstrLen))//'.'// &
                 trim(adjustl(cstrLenM9))//')'
     ompFormat = '(I'//trim(adjustl(cstrLen))//',A'//trim(adjustl(cstrLen))//','//trim(adjustl(cntimers))//'ES'// &
                 trim(adjustl(cstrLen))//'.'//trim(adjustl(cstrLenM9))//')'
     mpiFormat = '('//trim(adjustl(cstrLen))//'x,A'//trim(adjustl(cstrLen))//','//trim(adjustl(cntimers))//'ES'// &
                 trim(adjustl(cstrLen))//'.'//trim(adjustl(cstrLenM9))//')'

     nthreads = 1
     !$omp parallel default(shared)
     !$omp master
     !$ nthreads = omp_get_num_threads()
     !$omp end master
     !$omp end parallel

     allocate(t_MPI(ntimers,nthreads))
     allocate(t_ALL(ntimers,nthreads,bn_globalNumProcs))
     allocate(t_ALL_omp_min(ntimers,bn_globalNumProcs))
     allocate(t_ALL_omp_max(ntimers,bn_globalNumProcs))
     allocate(t_ALL_omp_avg(ntimers,bn_globalNumProcs))
     allocate(t_ALL_min(ntimers))
     allocate(t_ALL_max(ntimers))
     allocate(t_ALL_avg(ntimers))

     mythread = 0
     !$omp parallel default(shared) private(mythread)
     !$ mythread = omp_get_thread_num()
     t_MPI(1,mythread+1) = timer_burner
     t_MPI(2,mythread+1) = timer_xnet
     t_MPI(3,mythread+1) = timer_tstep
     t_MPI(4,mythread+1) = timer_nraph
     t_MPI(5,mythread+1) = timer_deriv
     t_MPI(6,mythread+1) = timer_jacob
     t_MPI(7,mythread+1) = timer_solve
     t_MPI(8,mythread+1) = timer_csect
     t_MPI(9,mythread+1) = timer_scrn
     !$omp end parallel

     call MPI_GATHER(t_MPI, ntimers*nthreads, MPI_DOUBLE_PRECISION, t_ALL, ntimers*nthreads, MPI_DOUBLE_PRECISION, 0, bn_GlobalComm, ierr)

     if (bn_globalMe == 0) then
        t_ALL_omp_min(:,:) = minval(t_ALL,dim=2)
        t_ALL_omp_max(:,:) = maxval(t_ALL,dim=2)
        t_ALL_omp_avg(:,:) = sum(t_ALL,dim=2) / nthreads
        t_ALL_min(:) = minval(t_ALL_omp_min,dim=2)
        t_ALL_max(:) = maxval(t_ALL_omp_max,dim=2)
        t_ALL_avg(:) = sum(t_ALL_omp_max,dim=2) / bn_globalNumProcs
        write(*,*)
        write(*,headerFormat) (tHeader(k),k=1,ntimers+2)
        do i = 1, bn_globalNumProcs
           write(*,'(A)') repeat("-",strLen*(ntimers+2))
           do j = 1, nthreads
              write(*,allFormat) i, j, (t_ALL(k,j,i),k=1,ntimers)
           end do
           write(*,'(2('//trim(cstrLen)//'x),A)') repeat("-",strLen*ntimers)
           write(*,ompFormat) i, 'MIN:', (t_ALL_omp_min(k,i),k=1,ntimers)
           write(*,ompFormat) i, 'MAX:', (t_ALL_omp_max(k,i),k=1,ntimers)
           write(*,ompFormat) i, 'AVG:', (t_ALL_omp_avg(k,i),k=1,ntimers)
        end do
        write(*,'('//trim(cstrLen)//'x,A)') repeat("-",strLen*(ntimers+1))
        write(*,headerFormat) (tHeader(k),k=1,ntimers+2)
        write(*,'('//trim(cstrLen)//'x,A)') repeat("-",strLen*(ntimers+1))
        write(*,mpiFormat) 'MIN:', (t_ALL_min(k),k=1,ntimers)
        write(*,mpiFormat) 'MAX:', (t_ALL_max(k),k=1,ntimers)
        write(*,mpiFormat) 'AVG:', (t_ALL_avg(k),k=1,ntimers)
     end if

     call MPI_BARRIER(bn_GlobalComm,ierr)

     deallocate(t_ALL)
     deallocate(t_ALL_omp_min)
     deallocate(t_ALL_omp_max)
     deallocate(t_ALL_omp_avg)
     deallocate(t_ALL_min)
     deallocate(t_ALL_max)
     deallocate(t_ALL_avg)
     deallocate(t_MPI)

  end if

  call jacobian_finalize

  return

end subroutine bn_xnetFinalize
