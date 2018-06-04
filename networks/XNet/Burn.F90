!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/XNet/Burn
!!
!! NAME
!!
!!  Burn
!!
!!
!! SYNOPSIS
!!
!!   call Burn ( integer, intent(IN)    :: blockCount, 
!!               integer(:), intent(IN) :: blockList, 
!!               real, intent(IN)       ::  dt  )    
!!
!! DESCRIPTION
!!
!!  Apply burner to all blocks in specified list.
!!
!! ARGUMENTS
!!
!!   blockCount -- dimension of blockList
!!   blockList -- array of blocks which should receive burning
!!   dt  --       passed to the internal bn_burner module  
!!
!! PARAMETERS
!!
!!  useBurn -- Boolean, True.  Turns on burning module
!!  useBurnTable -- Boolean, False.  Controls the generation of reaction rates.
!!                TRUE interpolates from a stored table; FALSE generates them
!!                analytically.
!!  useShockBurn -- Boolean, FALSE.  Controls whether burning is allowed inside
!!                a regime experiencing shocks
!!  algebra -- Integer, 1, [1,2].  Controls choice of linear algebra package used
!!                for matrix solution.  1=Ma28 sparse package, 2=Gift hardwired package.
!!  odeStepper -- Integer, 1, [1,2].  Controls time integration routines.
!!                1=Bader-Deuflhard variable order, 2=Rosenbrock 4th order
!!  nuclearTempMin/Max -- Real, 1.1E+8/1.0E+12.  Minimum and maximum temperature
!!                ranges where burning can occur
!!  nuclearDensMin/Max -- Real, 1.0E-10/1.0E+14.  Minimum and maximum density range
!!                where burning can occur.
!!  nuclearNI56Max -- Real, 1.0.  Maximum mass fraction of nickel where burning
!!                can occur.
!!  enucDtFactor -- Real, 1.0E+30.  Timestep limiter.See Burn_computeDt for details.              
!!
!! NOTES
!!
!!  The burning unit adds a new mesh variable ENUC_VAR which is the nuclear energy 
!!             generation rate
!!
!!***

!!REORDER(4): solnData

#include "Flash.h"

subroutine Burn (  blockCount, blockList, dt  )    

  use bn_interface, ONLY : bn_burner   
  use bn_xnetData, ONLY : xnet_myid, xnet_nzbatchmx, xnet_inuc2unk
  use Burn_data, ONLY : bn_nuclearTempMin, bn_nuclearTempMax, bn_nuclearDensMin, &
       &   bn_nuclearDensMax, bn_nuclearNI56Max, bn_useShockBurn, &
       &   bn_smallx, bn_useBurn, ionam
  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_interface, ONLY : Eos_wrapped
  use Grid_interface, ONLY : Grid_fillGuardCells, &
       Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_getBlkPtr, &
       Grid_releaseBlkPtr
  use Hydro_interface, ONLY : Hydro_detectShock
  use Simulation_interface, ONLY : Simulation_mapStrToInt
  use Timers_interface, ONLY : Timers_start, Timers_stop
#ifdef FLASH_GRID_PARAMESH
  use tree, ONLY : bflags
#endif

  !$ use omp_lib

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Eos.h"

  !args
  integer, intent(in) :: blockCount
  integer, intent(in) :: blockList(blockCount)
  real,    intent(in) :: dt

  ! locals
  integer :: blockID, thisBlock
  real, pointer, dimension(:,:,:,:) :: solnData
  real, allocatable, dimension(:)   :: xCoord, yCoord, zCoord
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer :: iSize, jSize, kSize, iSizeGC, jSizeGC, kSizeGC

  logical :: okBurnTemp, okBurnDens, okBurnShock, okBurnNickel
  logical :: getGuardCells = .true.

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC, &
                  GRID_JLO_GC:GRID_JHI_GC, &
                  GRID_KLO_GC:GRID_KHI_GC) :: shock
  real,    dimension(NSPECIES,NXB,NYB,NZB,blockCount), target :: xIn, xOut
  real,    dimension(NXB,NYB,NZB,blockCount),          target :: sdot, tmp, rho
  logical, dimension(NXB,NYB,NZB,blockCount),          target :: burnedZone
#else
  real,    allocatable :: shock(:,:,:)
  real,    allocatable, target :: xIn(:,:,:,:,:), xOut(:,:,:,:,:)
  real,    allocatable, target :: sdot(:,:,:,:), tmp(:,:,:,:), rho(:,:,:,:)
  logical, allocatable, target :: burnedZone(:,:,:,:)
#endif
  real,    pointer :: xIn_batch(:,:,:), xOut_batch(:,:,:)
  real,    pointer :: sdot_batch(:,:), tmp_batch(:,:), rho_batch(:,:)
  logical, pointer :: burnedZone_batch(:,:)

  integer :: nzones, batchCount
  integer, dimension(:), allocatable :: sumBurn_TS_batch
  integer, dimension(blockCount) :: batch_lo, batch_hi
  integer, dimension(blockCount) :: sumBurn_TS

  real :: ei, ek, enuc
  integer :: i, j, k, m, n, ii, jj, kk, mm, nn

  ! ----------------------- check if burning is requested in runtime parameters -------
  if (.not. bn_useBurn) return

  !---------------------------------------------------------------------------------

  ! start the timer ticking
  call Timers_start("burn")

  call Timers_start("burn_top")

  if (.NOT. bn_useShockBurn) then
     call Grid_fillGuardCells(CENTER, ALLDIR)
  endif

#ifndef FIXEDBLOCKSIZE
  call Grid_getBlkIndexLimits(1,blkLimits,blkLimitsGC)
  allocate(xIn(NSPECIES,iSize,jSize,kSize,blockCount))
  allocate(xOut(NSPECIES,iSize,jSize,kSize,blockCount))
  allocate(sdot(iSize,jSize,kSize,blockCount))
  allocate(tmp(iSize,jSize,kSize,blockCount))
  allocate(rho(iSize,jSize,kSize,blockCount))
  allocate(burnedZone(iSize,jSize,kSize,blockCount))
#endif

  burnedZone = .FALSE.
  nzones = 0

  ! loop over list of blocks passed in
  do thisBlock = 1, blockCount

     blockID = blockList(thisBlock)

     ! get dimensions/limits and coordinates
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     iSizeGC = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     jSizeGC = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     kSizeGC = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
     iSize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
     jSize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
     kSize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1

#ifndef FIXEDBLOCKSIZE
     allocate(shock(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
                    blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
                    blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
#endif

     ! identify the range of batches in each block (use floor/ceil in case of overlap)
     batch_lo(thisBlock) = nzones / xnet_nzbatchmx + 1
     nzones = nzones + iSize * jSize * kSize
     batch_hi(thisBlock) = (nzones + xnet_nzbatchmx - 1) / xnet_nzbatchmx

     ! allocate space for dimensions
     allocate(xCoord(iSizeGC))
     allocate(yCoord(jSizeGC))
     allocate(zCoord(kSizeGC))

     call Grid_getCellCoords(IAXIS,blockID,CENTER,getGuardCells,xCoord,iSizeGC)
     call Grid_getCellCoords(JAXIS,blockID,CENTER,getGuardCells,yCoord,jSizeGC)
     call Grid_getCellCoords(KAXIS,blockID,CENTER,getGuardCells,zCoord,kSizeGC)

     ! Get a pointer to solution data 
     call Grid_getBlkPtr(blockID,solnData)

     if (.NOT. bn_useShockBurn) then
        call Hydro_detectShock(solnData, shock, blkLimits, blkLimitsGC, (/0,0,0/), &
             xCoord,yCoord,zCoord)
     else
        shock(:,:,:) = 0
     endif

     solnData(NMPI_VAR,:,:,:) = xnet_myid

     !$omp parallel do &
     !$omp collapse(3) &
     !$omp default(shared) &
     !$omp private(k,kk,j,jj,i,ii,okBurnTemp,okBurnDens,okBurnShock,okBurnNickel)
#ifndef FIXEDBLOCKSIZE
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do k = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do k = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              kk = k - blkLimits(LOW,KAXIS) + 1
              jj = j - blkLimits(LOW,JAXIS) + 1
              ii = i - blkLimits(LOW,IAXIS) + 1
#else
     do k = GRID_KLO, GRID_KHI
        do j = GRID_JLO, GRID_JHI
           do i = GRID_ILO, GRID_IHI
              kk = k - GRID_KLO + 1
              jj = j - GRID_JLO + 1
              ii = i - GRID_ILO + 1
#endif

              tmp(ii,jj,kk,thisBlock)  = solnData(TEMP_VAR,i,j,k)
              rho(ii,jj,kk,thisBlock)  = solnData(DENS_VAR,i,j,k)
              sdot(ii,jj,kk,thisBlock) = 0.0e0

              ! Map the solution data into the order required by bn_burner
              xIn(1:NSPECIES,ii,jj,kk,thisBlock) = solnData(xnet_inuc2unk,i,j,k)

              okBurnTemp = .FALSE.
              okBurnDens = .FALSE.
              okBurnShock = .FALSE.
              okBurnNickel = .FALSE.

              okBurnTemp = (tmp(ii,jj,kk,thisBlock) >= bn_nuclearTempMin .AND. tmp(ii,jj,kk,thisBlock) <= bn_nuclearTempMax)
              okBurnDens = (rho(ii,jj,kk,thisBlock) >= bn_nuclearDensMin .AND. rho(ii,jj,kk,thisBlock) <= bn_nuclearDensMax)
              okBurnShock = (shock(i,j,k) == 0.0 .OR. (shock(i,j,k) == 1.0 .AND. bn_useShockBurn))

              if (okBurnTemp .AND. okBurnDens .AND. okBurnShock) then

                 if (NI56_SPEC /= NONEXISTENT) then
                    okBurnNickel = (solnData(NI56_SPEC,i,j,k) <  bn_nuclearNI56Max)
                 else    ! nickel is not even a species in this simulation, so we'll always burn
                    okBurnNickel = .TRUE.
                 endif

                 if (okBurnNickel) then
                    burnedZone(ii,jj,kk,thisBlock) = .TRUE.
                 endif

              endif

           end do
        end do
     end do
     !$omp end parallel do

     call Grid_releaseBlkPtr(blockID,solnData)
     deallocate(xCoord)
     deallocate(yCoord)
     deallocate(zCoord)

#ifndef FIXEDBLOCKSIZE
     deallocate(shock)
#endif

  end do

  call Timers_stop("burn_top")

  call Timers_start("burn_middle")

  ! get number of batches needed for all local zones (round up)
  batchCount = (nzones + xnet_nzbatchmx - 1) / xnet_nzbatchmx

  ! reshape all local zone data arrays into batches
  tmp_batch (1:xnet_nzbatchmx,1:batchCount) => tmp (:,:,:,:)
  rho_batch (1:xnet_nzbatchmx,1:batchCount) => rho (:,:,:,:)
  sdot_batch(1:xnet_nzbatchmx,1:batchCount) => sdot(:,:,:,:)
  xIn_batch (1:NSPECIES,1:xnet_nzbatchmx,1:batchCount) => xIn (:,:,:,:,:)
  xOut_batch(1:NSPECIES,1:xnet_nzbatchmx,1:batchCount) => xOut(:,:,:,:,:)
  burnedZone_batch(1:xnet_nzbatchmx,1:batchCount) => burnedZone(:,:,:,:)

  allocate(sumBurn_TS_batch(batchCount))

  !$omp parallel do &
  !$omp schedule(runtime) &
  !$omp default(shared)
  do m = 1, batchCount
     ! Do the actual burn
     call bn_burner(dt, tmp_batch(:,m), rho_batch(:,m), xIn_batch(:,:,m), &
          xOut_batch(:,:,m), sdot_batch(:,m), burnedZone_batch(:,m), sumBurn_TS_batch(m))
  end do
  !$omp end parallel do

  ! get maximum sumBurn_TS_batch per block
  sumBurn_TS = 0
  do thisBlock = 1, blockCount
     sumBurn_TS(thisBlock) = maxval( sumBurn_TS_batch(batch_lo(thisBlock):batch_hi(thisBlock)) )
  end do

  deallocate(sumBurn_TS_batch)

  call Timers_stop("burn_middle")

  call Timers_start("burn_bottom")

  ! loop over list of blocks passed in
  do thisBlock = 1, blockCount

     blockID = blockList(thisBlock)

     ! get dimensions/limits and coordinates
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ! Get a pointer to solution data 
     call Grid_getBlkPtr(blockID,solnData)

     ! Now put updated local data arrays back into unk through solnData pointer
     !$omp parallel do &
     !$omp collapse(3) &
     !$omp default(shared) &
     !$omp private(k,kk,j,jj,i,ii,ei,ek,enuc)
#ifndef FIXEDBLOCKSIZE
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do k = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do k = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              kk = k - blkLimits(LOW,KAXIS) + 1
              jj = j - blkLimits(LOW,JAXIS) + 1
              ii = i - blkLimits(LOW,IAXIS) + 1
#else
     do k = GRID_KLO, GRID_KHI
        do j = GRID_JLO, GRID_JHI
           do i = GRID_ILO, GRID_IHI
              kk = k - GRID_KLO + 1
              jj = j - GRID_JLO + 1
              ii = i - GRID_ILO + 1
#endif

              ! Map the solution data into the order required by bn_burner
              solnData(xnet_inuc2unk,i,j,k) = xOut(1:NSPECIES,ii,jj,kk,thisBlock)

              !  NOTE should probably do something here with eintSwitch for consistency
              !  LBR will settle for simply using internal energy!
              ! kinetic energy
              ek = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 +  &
                 solnData(VELY_VAR,i,j,k)**2 +  &
                 solnData(VELZ_VAR,i,j,k)**2)

              ! internal energy, add on nuclear rate*timestep
              enuc = dt*sdot(ii,jj,kk,thisBlock)
              ei = solnData(ENER_VAR,i,j,k) + enuc - ek

#ifdef EINT_VAR
              solnData(EINT_VAR,i,j,k) = ei
#endif
              solnData(ENER_VAR,i,j,k) = ei + ek
#ifdef EELE_VAR
              solnData(EELE_VAR,i,j,k) = solnData(EELE_VAR,i,j,k) + enuc
#endif
              solnData(ENUC_VAR,i,j,k) = sdot(ii,jj,kk,thisBlock)

           end do
        end do
     end do
     !$omp end parallel do
   
#ifdef FLASH_GRID_PARAMESH
     bflags(1,blockID) = sumBurn_TS(thisBlock)
#endif
     solnData(MTSB_VAR,:,:,:) = sumBurn_TS(thisBlock)

     ! we've altered the EI, let's equilabrate
     call Timers_start("eos")
     if (any(burnedZone(:,:,:,thisBlock))) then

#ifdef FLASH_UHD_3T
        call Eos_wrapped(MODE_DENS_EI_GATHER,blkLimits,blockID) ! modified for 3T
#else 
        call Eos_wrapped(MODE_DENS_EI,blkLimits,blockID)
#endif

     end if
     call Timers_stop("eos")

     call Grid_releaseBlkPtr(blockID,solnData)

  end do

  call Timers_stop("burn_bottom")

#ifndef FIXEDBLOCKSIZE
  deallocate(xIn)
  deallocate(xOut)
  deallocate(sdot)
  deallocate(tmp)
  deallocate(rho)
  deallocate(burnedZone)
#endif

  call Timers_stop("burn")

  return
end subroutine Burn
