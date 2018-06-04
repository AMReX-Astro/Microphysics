module actual_network

   use bl_types, only: dp_t
   use xnet_constants, only: avn, clt, epmev, m_n, m_p, m_e

   implicit none

   real(dp_t), parameter, private :: clight = clt
   real(dp_t), parameter, private :: mev2erg = epmev
   real(dp_t), parameter, private :: ev2erg  = mev2erg*1.0d-6
   real(dp_t), parameter, private :: mev2gr  = mev2erg/clight**2

   real(dp_t), parameter, private :: mn = m_n*mev2gr
   real(dp_t), parameter, private :: mp = m_p*mev2gr
   real(dp_t), parameter, private :: me = m_e*mev2gr

   integer :: nspec
   integer :: nspec_evolve
   integer :: naux

   character (len=16), allocatable :: spec_names(:) 
   character (len= 5), allocatable :: short_spec_names(:)
   character (len= 5), allocatable :: short_aux_names(:)

   character (len=32), save :: network_name = "xnet"

   real(dp_t), allocatable :: aion(:), zion(:), bion(:)
   real(dp_t), allocatable :: nion(:), mion(:), wion(:)

   !$acc declare create(aion, zion, bion, nion, mion, wion)

   integer :: nrates
   integer :: num_rate_groups

   ! Conversion factor for the nuclear energy generation rate.

   real(dp_t), parameter :: avo = avn
   real(dp_t), parameter :: enuc_conv2 = -avo*clight*clight

contains

   subroutine actual_network_init

      use controls, only: mythread, nthread, myid, nproc, nzbatchmx, iweak0, &
         iscrn, iprocess, isolv, kstmx, kitmx, ijac, iconvc, changemx, yacc, &
         tolm, tolc, ymin, tdel_maxmult, iheat, changemxt, tolt9, idiag, itsout
      use extern_probin_module, only: xnet_data_dir, xnet_nzbatchmx, xnet_iweak0, &
         xnet_iscrn, xnet_iprocess, xnet_isolv, xnet_kstmx, &
         xnet_kitmx, xnet_ijac, xnet_iconvc, xnet_changemx, &
         xnet_yacc, xnet_tolm, xnet_tolc, xnet_ymin, &
         xnet_tdel_maxmult, xnet_iheat, xnet_changemxt, &
         xnet_tolt9, xnet_idiag, xnet_itsout
      use nuclear_data, only: aa, zz, nn, be, mm, nname
      use nuc_number, only: ny
      !use parallel
      use xnet_interface, only: xnet_init

      use mpi
      !$use omp_lib

      implicit none

      integer :: i, ierr
      character (len=80) :: data_desc

      ! Initialize MPI/OpenMP identifiers
      ! TODO: get these from parallel module?
      call mpi_comm_rank( MPI_COMM_WORLD, myid, ierr)
      call mpi_comm_size( MPI_COMM_WORLD, nproc, ierr)
      mythread = 1
      nthread = 1
      !$omp parallel default(shared)
      !$ mythread = omp_get_thread_num()
      !$omp master
      !$ nthread = omp_get_num_threads()
      !$omp end master
      !$omp end parallel

      nzbatchmx = xnet_nzbatchmx
      iweak0 = xnet_iweak0
      iscrn = xnet_iscrn
      iprocess = xnet_iprocess
      isolv = xnet_isolv
      kstmx = xnet_kstmx
      kitmx = xnet_kitmx
      ijac = xnet_ijac
      iconvc = xnet_iconvc
      changemx = xnet_changemx
      yacc = xnet_yacc
      tolm = xnet_tolm
      tolc = xnet_tolc
      ymin = xnet_ymin
      tdel_maxmult = xnet_tdel_maxmult
      iheat = xnet_iheat
      changemxt = xnet_changemxt
      tolt9 = xnet_tolt9
      idiag = xnet_idiag
      itsout = xnet_itsout

      call xnet_init(xnet_data_dir,data_desc)

      nspec = ny
      nspec_evolve = ny
      naux = 0

      allocate(spec_names(nspec))
      allocate(short_spec_names(nspec))
      allocate(short_aux_names(naux))

      allocate(aion(nspec))
      allocate(zion(nspec))
      allocate(nion(nspec))
      allocate(bion(nspec))
      allocate(mion(nspec))
      allocate(wion(nspec))

      do i = 1, nspec
         spec_names(i) = nname(i)
         short_spec_names(i) = nname(i)
      end do
      aion(:) = aa(:)
      zion(:) = zz(:)
      nion(:) = nn(:)
      bion(:) = be(:)
      mion(:) = mm(:)
      wion(:) = aion(:)

      !$acc update device(aion, zion, bion, nion, mion, wion)

      nrates = sum(nreac(:))
      num_rate_groups = 1

   end subroutine actual_network_init

   subroutine actual_network_finalize

      use xnet_interface, only: xnet_finalize

      implicit none

      call xnet_finalize

   end subroutine actual_network_finalize

end module actual_network
