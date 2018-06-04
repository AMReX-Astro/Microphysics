!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/XNet/bn_initNetwork
!!
!! NAME
!!
!!  bn_initNetwork
!!
!!
!! SYNOPSIS
!!
!!  bn_initNetwork()
!!
!!
!! DESCRIPTION
!!
!!  Initializes variables used in XNet.
!!  Called from Burn_init.
!!
!!***


subroutine bn_initNetwork
  use Burn_data, ONLY: aion, bion, zion, aioninv, zionsq, ionam
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use bn_xnetData, ONLY : xnet_myid, xnet_nproc, mythread, nthread, &
                          xnet_data_dir, xnet_nzbatchmx, xnet_iweak0, &
                          xnet_iscrn, xnet_iprocess, xnet_isolv, xnet_kstmx, &
                          xnet_kitmx, xnet_ijac, xnet_iconvc, xnet_changemx, &
                          xnet_yacc, xnet_tolm, xnet_tolc, xnet_ymin, &
                          xnet_tdel_maxmult, xnet_iheat, xnet_changemxt, &
                          xnet_tolt9, xnet_idiag, xnet_itsout, &
                          xnet_writeTimers, xnet_data_desc, xnet_aa, xnet_zz, &
                          xnet_be, xnet_nname
  use bn_interface, ONLY : bn_xnetInit

  !$ use omp_lib

  implicit none

#include "Flash.h"
#include "constants.h"

  character(len=5) :: tmp_name
  integer :: i

  ! Initialize MPI/OpenMP identifiers
  call Driver_getMype(GLOBAL_COMM,xnet_myid)
  call Driver_getNumProcs(GLOBAL_COMM,xnet_nproc)
  mythread = 1
  nthread = 1
  !$omp parallel default(shared)
  !$ mythread = omp_get_thread_num()
  !$omp master
  !$ nthread = omp_get_num_threads()
  !$omp end master
  !$omp end parallel

  ! Read and distribute control file data (get from FLASH runtime parameters instead)
  !call control_bcast(xnet_data_dir)

  ! XNet/REACLIB Data directory
  call RuntimeParameters_get('xnet_data_dir',xnet_data_dir)

  ! XNet controls
  call RuntimeParameters_get('xnet_nzbatchmx',xnet_nzbatchmx)
  call RuntimeParameters_get('xnet_iweak',xnet_iweak0)
  call RuntimeParameters_get('xnet_iscrn',xnet_iscrn)
  call RuntimeParameters_get('xnet_iprocess',xnet_iprocess)
  call RuntimeParameters_get('xnet_isolv',xnet_isolv)
  call RuntimeParameters_get('xnet_kstmx',xnet_kstmx)
  call RuntimeParameters_get('xnet_kitmx',xnet_kitmx)
  call RuntimeParameters_get('xnet_ijac',xnet_ijac)
  call RuntimeParameters_get('xnet_iconvc',xnet_iconvc)
  call RuntimeParameters_get('xnet_changemx',xnet_changemx)
  call RuntimeParameters_get('xnet_yacc',xnet_yacc)
  call RuntimeParameters_get('xnet_tolm',xnet_tolm)
  call RuntimeParameters_get('xnet_tolc',xnet_tolc)
  call RuntimeParameters_get('xnet_ymin',xnet_ymin)
  call RuntimeParameters_get('xnet_tdel_maxmult',xnet_tdel_maxmult)
  call RuntimeParameters_get('xnet_iheat',xnet_iheat)
  call RuntimeParameters_get('xnet_changemxt',xnet_changemxt)
  call RuntimeParameters_get('xnet_tolt9',xnet_tolt9)
  call RuntimeParameters_get('xnet_idiag',xnet_idiag)
  !call RuntimeParameters_get('xnet_itsout',xnet_itsout)
  xnet_itsout = 0

  call RuntimeParameters_get('xnet_writeTimers',xnet_writeTimers)

  ! Read, broadcast, and allocate XNet data
  call bn_xnetInit(xnet_data_dir,xnet_data_desc)

  ! Some bn_* routines get nuclear data from Burn_data, so copy from XNet data
  do i = 1, NSPECIES
    tmp_name = adjustl(xnet_nname(i))
    ionam(i) = tmp_name(1:4)
    aion(i) = xnet_aa(i)
    zion(i) = xnet_zz(i)
    bion(i) = xnet_be(i)
    aioninv(i) = 1.0e0 / aion(i)
    zionsq(i) = zion(i) * zion(i)
  end do

  return
end subroutine bn_initNetwork
