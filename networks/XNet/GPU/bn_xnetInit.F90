!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/XNet/GPU/bn_xnetInit
!!
!! NAME
!!
!!  bn_xnetInit
!!
!!
!! SYNOPSIS
!!
!!  bn_xnetInit ( character(*), intent(IN) :: data_dir,
!!                character(80), intent(OUT) :: data_desc )
!!
!! DESCRIPTION
!!
!!  Calls XNet routines to read and broadcast XNet data 
!!  Called from bn_initNetwork.
!!
!! ARGUMENTS
!!
!!   data_dir -- nuclear data directory containing REACLIB-formatted data
!!   data_desc -- brief description of  network
!!
!!***

subroutine bn_xnetInit(data_dir,data_desc)
  use abundances, ONLY : ystart, yo, y, yt, ydot
  Use conditions, ONLY : t, tt, to, tdel, tdel_next, tdel_old, t9t, rhot, yet, &
    t9, rho, ye, t9o, rhoo, yeo, t9dot, cv, nt, ntt, nto, ints, intso
  use controls, ONLY : idiag, iheat, iscrn, iprocess, iweak, &
    kmon, ktot, nzbatchmx, lzactive, myid, nproc, mythread, nthread, &
    lun_diag, lun_ev, lun_ts, lun_stdout, getNewUnit
  use gpu_controls, ONLY : gpu_init
  use nuc_number, ONLY : ny
  Use thermo_data, ONLY : nstart, tstart, tstop, tdelstart, t9start, rhostart, &
    yestart, nh, th, t9h, rhoh, yeh, nhmx
  use xnet_eos, ONLY : eos_initialize
  use xnet_interface, ONLY : flux_init, name_ordered, net_preprocess
  use xnet_interface_mpi, ONLY : control_bcast, jacobian_bcast, match_bcast, netdata_bcast

  implicit none

#include "Flash.h"
#include "constants.h"

  ! Input variables
  character(*), intent(in) :: data_dir

  ! Output variables
  character(80), intent(out) :: data_desc

  ! Local variables
  character(80) :: diag_file
  character(80), parameter :: diag_file_base = 'xnet_diag'

  call gpu_init

  !$omp parallel default(shared) private(diag_file)

  ! Open diagnositic output file
  if ( idiag > 0 ) then
    diag_file = trim(diag_file_base)
    call name_ordered(diag_file,myid,nproc)
    if ( nthread > 1 ) call name_ordered(diag_file,mythread,nthread)
    open(getNewUnit(lun_diag), file=diag_file)
  else
    lun_diag = lun_stdout
  endif

  ! Allocate control arrays
  allocate (lzactive(nzbatchmx))
  allocate (iweak(nzbatchmx),lun_ev(nzbatchmx),lun_ts(nzbatchmx))
  allocate (kmon(2,nzbatchmx),ktot(2,nzbatchmx))

  !$omp end parallel

  if ( iprocess > 0 .and. myid == MASTER_PE ) call net_preprocess( lun_stdout, data_dir, data_dir )

  ! Read and distribute nuclear and reaction data
  call netdata_bcast(data_dir,data_desc)

  ! Read in and distribute data for reaction matching and fluxes
  call match_bcast(data_dir)

  ! Initialize flux tracking
  call flux_init

  ! BCast jacobian data
  call jacobian_bcast(data_dir)

  ! Initialize EoS for screening or self-heating
  If ( iscrn > 0 .or. iheat > 0 ) Call eos_initialize

  !$omp parallel default(shared)

  ! Allocate abundance arrays
  allocate (y(ny,nzbatchmx),yo(ny,nzbatchmx),yt(ny,nzbatchmx),ydot(ny,nzbatchmx),ystart(ny,nzbatchmx))

  ! Allocate conditions arrays
  allocate (t(nzbatchmx),tt(nzbatchmx),to(nzbatchmx), &
    &       tdel(nzbatchmx),tdel_next(nzbatchmx),tdel_old(nzbatchmx), &
    &       t9(nzbatchmx),t9t(nzbatchmx),t9o(nzbatchmx), &
    &       rho(nzbatchmx),rhot(nzbatchmx),rhoo(nzbatchmx), &
    &       ye(nzbatchmx),yet(nzbatchmx),yeo(nzbatchmx), &
    &       nt(nzbatchmx),ntt(nzbatchmx),nto(nzbatchmx), &
    &       ints(nzbatchmx),intso(nzbatchmx), &
    &       t9dot(nzbatchmx),cv(nzbatchmx))

  ! Allocate thermo history arrays
  allocate (nh(nzbatchmx),nstart(nzbatchmx), &
    &       tstart(nzbatchmx),tstop(nzbatchmx),tdelstart(nzbatchmx), &
    &       t9start(nzbatchmx),rhostart(nzbatchmx),yestart(nzbatchmx), &
    &       th(nhmx,nzbatchmx),t9h(nhmx,nzbatchmx),rhoh(nhmx,nzbatchmx),yeh(nhmx,nzbatchmx))

  !$omp end parallel

  return
end subroutine bn_xnetInit