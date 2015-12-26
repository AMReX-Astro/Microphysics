! Aaron Jackson 2009
!
! This subroutine communicates whether any process has found ignition
! conditions.
! If processes have, the master process determines the order of communication.
! Then, sequentially, each process with ignition conditions broadcasts the
! corresponding detonation points to everyone.
! After this communication, each process throws out detonation points that
! are too close to one another, preferentially keeping detonations further
! away from the center of the star

subroutine bn_paraAllSpark( ignition_coords, det_num, &
                            det_xCoord, det_yCoord, det_zCoord )

  use Burn_data, ONLY : pbIgnRad, pbIgnDist, pbIgnSep, pbIgnNum, pbIgnNumMax, &
                        pbIgnX, pbIgnY, pbIgnZ, &
                        bn_meshMe, bn_meshComm, bn_meshNumProcs

  implicit none

 include 'mpif.h'
#include "constants.h"  
#include "Flash.h"

  real, allocatable, dimension(:), intent(inout) :: ignition_coords
  integer, intent(inout) :: det_num
  real, allocatable, dimension(:), intent(out) :: det_xCoord, det_yCoord, det_zCoord

  logical, allocatable, dimension(:) :: detMask
  integer, allocatable, dimension(:) :: bcast_order, displ
  integer :: bcastNum, bcastNumMax, ignNum
  integer :: i, j, k, l, ii, jj, kk, ll, istat
  real :: dist

  real, allocatable, dimension(:) :: recv_buf
  real, allocatable, dimension(:) :: radius

  allocate(bcast_order(0:bn_meshNumProcs-1),stat=istat)
  if (istat/=0) call Driver_abortFlash("Unable to allocate bcast_order")

  ! gather all ignition_conditions first
  call MPI_AllGather( det_num, 1, MPI_INTEGER,     &
                      bcast_order, 1, MPI_INTEGER, &
                      bn_meshComm, istat )
  if (istat/=0) &
     call Driver_abortFlash("Error occured while calling MPI_AllGather")

  det_num = SUM(bcast_order)

  ! if there are no ignition conditions, then we are done
  if (det_num .eq. 0) return

  allocate(recv_buf(det_num*NDIM),stat=istat)
  if (istat/=0) call Driver_abortFlash("Unable to allocate recv_buf")
  allocate(displ(0:bn_meshNumProcs-1),stat=istat)
  if (istat/=0) call Driver_abortFlash("Unable to allocate displ")


  displ(0) = 0
  bcast_order(0) = NDIM * bcast_order(0)
  do i = 1, bn_meshNumProcs-1
     displ(i) = displ(i-1) + bcast_order(i-1)
     bcast_order(i) = NDIM * bcast_order(i)
  enddo

  call MPI_AllGatherV( ignition_coords,      &
                       bcast_order(bn_meshMe),&
                       MPI_DOUBLE_PRECISION, &
                       recv_buf,             &
                       bcast_order,          &
                       displ,                &
                       MPI_DOUBLE_PRECISION, &
                       bn_meshComm,       &
                       istat )
  if (istat/=0)  &
     call Driver_abortFlash("Error occured while calling MPI_AllGatherV")

  if ( allocated(ignition_coords) ) deallocate(ignition_coords)
  deallocate(bcast_order)
  deallocate(displ)

  allocate(radius(det_num),stat=istat)
  if (istat/=0) call Driver_abortFlash("Unable to allocate radius")

  ! otherwise we need to allocate space for the detonation coordinates
  allocate(det_xCoord(det_num),stat=istat)
  if (istat/=0) call Driver_abortFlash("Unable to allocate det_xCoord")

#if NDIM >= 2
  allocate(det_yCoord(det_num),stat=istat)
  if (istat/=0) call Driver_abortFlash("Unable to allocate det_yCoord")
#endif

#if NDIM > 2
  allocate(det_zCoord(det_num),stat=istat)
  if (istat/=0) call Driver_abortFlash("Unable to allocate det_zCoord")
#endif

  do i = 1, det_num 

     det_xCoord(i) = recv_buf( (i-1)*NDIM + 1)
     radius(i) = det_xCoord(i)**2
#if NDIM >= 2
     det_yCoord(i) = recv_buf( (i-1)*NDIM + 2)
     radius(i) = radius(i) + det_yCoord(i)**2
#endif
#if NDIM > 2
     det_zCoord(i) = recv_buf( (i-1)*NDIM + 3)
     radius(i) = radius(i) + det_zCoord(i)**2
#endif
     radius(i) = sqrt( radius(i) )

  enddo

  deallocate(recv_buf)

  allocate(detMask(det_num),stat=istat)
  if (istat/=0) call Driver_abortFlash("Unable to allocate detMask")

  detMask(:) = .false.
  ! now determine if any detonation points overlap
  do l = 1, det_num

     ! check with other detonation points found this time step
     do ll = l+1, det_num

        ! make sure we are not comparing the same point and that it has
        ! not been previously removed
        if ( .not. (detMask(l) .or. detMask(ll)) ) then

           dist = ( det_xCoord(l) - det_xCoord(ll) )**2
#if NDIM >= 2
           dist = dist + ( det_yCoord(l) - det_yCoord(ll) )**2
#endif
#if NDIM > 2
           dist = dist + ( det_zCoord(l) - det_zCoord(ll) )**2
#endif
           dist = sqrt( dist )

           if ( dist <= pbIgnSep ) then

              if (radius(l) >= radius(ll)) then
                 ! get rid of 2nd radius point
                 detMask(ll) = .true.
              else
                 ! get rid of 1st radius point
                 detMask(l) = .true.
              endif

           endif

        endif

     enddo

  enddo

  deallocate(radius)

  ! now adjust and move arrays
  if (any(detMask)) then

     l = 1
     do while (l <= det_num)

        if (detMask(l)) then
           ! remove this point
           det_num = det_num - 1

           do ll = l, det_num

              det_xCoord(ll) = det_xCoord(ll+1)
#if NDIM >= 2
              det_yCoord(ll) = det_yCoord(ll+1)
#endif
#if NDIM > 2
              det_zCoord(ll) = det_zCoord(ll+1)
#endif
              detMask(ll) = detMask(ll+1)

           enddo

           ! move l back 1
           l = l - 1

        endif

        l = l + 1

     enddo

  endif

  deallocate(detMask)

  return
  
end subroutine bn_paraAllSpark
