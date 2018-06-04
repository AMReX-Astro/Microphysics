!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/XNet/bn_xnetData
!!
!! NAME
!!  
!!  bn_xnetData
!!
!!
!! SYNOPSIS
!! 
!!  use bn_xnetData
!!
!! DESCRIPTION
!!
!!  Contains variables required by the XNet nuclear network
!!
!!***

module bn_xnetData

  ! Rename XNet runtime controls to match FLASH naming conventions
  use controls, ONLY : &
    xnet_iweak0 => iweak0, &             ! Saves input iweak flag
    xnet_iscrn => iscrn, &               ! If =0, screening is ignored
    xnet_iprocess => iprocess, &         ! If >0, then run network pre-processing (slows calculation)
    xnet_nzbatchmx => nzbatchmx, &       ! Maximum number of zones in a batch
    xnet_isolv => isolv, &               ! Sets the integration method (1=BE, 2=BD)
    xnet_kstmx => kstmx, &               ! Max # of timesteps before exit
    xnet_kitmx => kitmx, &               ! Max # of iterations or substeps within a timestep
    xnet_ijac => ijac, &                 ! Rebuild the Jacobian every ijac iterations after the first
    xnet_iconvc => iconvc, &             ! Controls type of convergence condition (0=mass)
    xnet_changemx => changemx, &         ! Relative abundance change used to guess timestep
    xnet_yacc => yacc, &                 ! Min abundance required to be accuracte, used in timestep determination.
    xnet_tolc => tolc, &                 ! The iterative convergence test limit
    xnet_tolm => tolm, &                 ! Max network mass error
    xnet_ymin => ymin, &                 ! Minimum abundance, y < ymin = 0
    xnet_tdel_maxmult => tdel_maxmult, & ! new timestep <= tdel_maxmult * previous timestep
    xnet_iheat => iheat, &               ! If >0, then implicitly couple network to temperature (self-heating)
    xnet_changemxt => changemxt, &       ! Relative temperature change used to guess timestep
    xnet_tolt9 => tolt9, &               ! The iterative convergence test limit
    xnet_idiag => idiag, &               ! Sets level of diagnostic output
    xnet_itsout => itsout, &             ! Sets level of time series output
    xnet_lun_diag => lun_diag, &         ! Logical units for per-thread diagnostic output file
    xnet_myid => myid, &                 ! Global MPI rank ID
    xnet_nproc => nproc, &               ! Number of MPI ranks
    mythread, &                          ! OpenMP thread ID
    nthread                              ! Number of OpenMP threads
  
  use nuclear_data, ONLY : &
    xnet_aa => aa, &
    xnet_zz => zz, &
    xnet_nn => nn, &
    xnet_be => be, &
    xnet_nname => nname

  use xnet_eos, ONLY : &
    xnet_inuc2unk => inuc2unk
  
  ! Other controls needed for XNet to operate in FLASH
  character(80) :: xnet_data_dir         ! Nuclear data directory
  character(80) :: xnet_data_desc        ! Brief network description

  logical :: xnet_writeTimers

end module bn_xnetData
