! Aaron Jackson 2009

subroutine bn_paraAllIgnite( ignition_conditions, det_num )

  use Burn_data, ONLY: bn_meshComm

  implicit none

 include 'mpif.h'
#include "constants.h"  
#include "Flash.h"

  integer, intent(in) :: det_num
  logical, allocatable, dimension(:), intent(inout) :: ignition_conditions

  logical, dimension(det_num) :: send_buf
  integer :: ierr


  send_buf = ignition_conditions

  call MPI_AllReduce( send_buf,            &
                      ignition_conditions, &
                      det_num,             &
                      MPI_LOGICAL,         &
                      MPI_LAND,            &
                      bn_meshComm,         &
                      ierr )

  return
  
end subroutine bn_paraAllIgnite
