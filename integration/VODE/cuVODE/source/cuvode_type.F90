module cuvode_types_module

  use amrex_fort_module, only: rt => amrex_real
  use cuvode_parameters_module, only: VODE_NEQS
  use vode_rpar_indices, only: n_rpar_comps

  implicit none

  real(rt), parameter :: UROUND = epsilon(1.0_rt)
  real(rt), parameter :: SRUR = sqrt(UROUND)
  real(rt), parameter :: CCMXJ = 0.2e0_rt
  real(rt), parameter :: HMIN = 0.0_rt
  real(rt), parameter :: HMXI = 0.0_rt

  ! For the backward differentiation formula (BDF) integration
  ! the maximum order should be no greater than 5.
  integer, parameter :: VODE_MAXORD = 5
  integer, parameter :: VODE_LMAX = VODE_MAXORD + 1

  ! How many timesteps should pass before refreshing the Jacobian
  integer, parameter :: max_steps_between_jacobian_evals = 50

  ! Type dvode_t contains the integration solution and control variables
  type :: dvode_t
     real(rt) :: HU
     real(rt) :: ACNRM, CONP, CRATE, DRC, EL(13)
     real(rt) :: ETA, ETAMAX, H, HNEW, HSCAL, PRL1
     real(rt) :: RC, RL1, TAU(13), TQ(5), TN
     integer  :: NFE, NJE, NST
     integer  :: ICF, IPUP, JCUR, JSTART, JSV, KFLAG
     integer  :: L
     integer  :: MXSTEP
     integer  :: NEWH, NEWQ, NQ, NQNYH, NQWAIT, NSLJ
     integer  :: NSLP

     ! Tolerances
     real(rt) :: RTOL(VODE_NEQS), ATOL(VODE_NEQS)

     ! Real parameters
     real(rt) :: RPAR(n_rpar_comps)

     ! State flag
     integer :: ISTATE

     ! Local time and integration end time
     real(rt) :: T, TOUT

     ! Integration array
     real(rt) :: Y(VODE_NEQS)

     ! Jacobian
     real(rt) :: jac(VODE_NEQS*VODE_NEQS)

     ! Saved Jacobian
     real(rt) :: jac_save(VODE_NEQS*VODE_NEQS)

     real(rt) :: yh(VODE_NEQS, VODE_LMAX)
     real(rt) :: ewt(VODE_NEQS)
     real(rt) :: savf(VODE_NEQS)
     real(rt) :: acor(VODE_NEQS)

     ! Jacobian method
     integer  :: jacobian

  end type dvode_t

contains

#ifndef AMREX_USE_CUDA
  subroutine print_state(dvode_state)
    use amrex_fort_module, only : rt => amrex_real
    type(dvode_t) :: dvode_state

    write(*,*) 'HU = ', dvode_state % HU
    write(*,*) 'ACNRM = ', dvode_state % ACNRM
    write(*,*) 'CONP = ', dvode_state % CONP
    write(*,*) 'CRATE = ', dvode_state % CRATE
    write(*,*) 'DRC = ', dvode_state % DRC
    write(*,*) 'EL(1) = ', dvode_state % EL(1)
    write(*,*) 'EL(2) = ', dvode_state % EL(2)
    write(*,*) 'EL(3) = ', dvode_state % EL(3)
    write(*,*) 'EL(4) = ', dvode_state % EL(4)
    write(*,*) 'EL(5) = ', dvode_state % EL(5)
    write(*,*) 'EL(6) = ', dvode_state % EL(6)
    write(*,*) 'EL(7) = ', dvode_state % EL(7)
    write(*,*) 'EL(8) = ', dvode_state % EL(8)
    write(*,*) 'EL(9) = ', dvode_state % EL(9)
    write(*,*) 'EL(10) = ', dvode_state % EL(10)
    write(*,*) 'EL(11) = ', dvode_state % EL(11)
    write(*,*) 'EL(12) = ', dvode_state % EL(12)
    write(*,*) 'EL(13) = ', dvode_state % EL(13)
    write(*,*) 'ETA = ', dvode_state % ETA
    write(*,*) 'ETAMAX = ', dvode_state % ETAMAX
    write(*,*) 'H = ', dvode_state % H
    write(*,*) 'HNEW = ', dvode_state % HNEW
    write(*,*) 'HSCAL = ', dvode_state % HSCAL
    write(*,*) 'PRL1 = ', dvode_state % PRL1
    write(*,*) 'RC = ', dvode_state % RC
    write(*,*) 'RL1 = ', dvode_state % RL1
    write(*,*) 'TAU(1) = ', dvode_state % TAU(1)
    write(*,*) 'TAU(2) = ', dvode_state % TAU(2)
    write(*,*) 'TAU(3) = ', dvode_state % TAU(3)
    write(*,*) 'TAU(4) = ', dvode_state % TAU(4)
    write(*,*) 'TAU(5) = ', dvode_state % TAU(5)
    write(*,*) 'TAU(6) = ', dvode_state % TAU(6)
    write(*,*) 'TAU(7) = ', dvode_state % TAU(7)
    write(*,*) 'TAU(8) = ', dvode_state % TAU(8)
    write(*,*) 'TAU(9) = ', dvode_state % TAU(9)
    write(*,*) 'TAU(10) = ', dvode_state % TAU(10)
    write(*,*) 'TAU(11) = ', dvode_state % TAU(11)
    write(*,*) 'TAU(12) = ', dvode_state % TAU(12)
    write(*,*) 'TAU(13) = ', dvode_state % TAU(13)
    write(*,*) 'TQ(1) = ', dvode_state % TQ(1)
    write(*,*) 'TQ(2) = ', dvode_state % TQ(2)
    write(*,*) 'TQ(3) = ', dvode_state % TQ(3)
    write(*,*) 'TQ(4) = ', dvode_state % TQ(4)
    write(*,*) 'TQ(5) = ', dvode_state % TQ(5)
    write(*,*) 'TN = ', dvode_state % TN
    write(*,*) 'NFE = ', dvode_state % NFE
    write(*,*) 'NJE = ', dvode_state % NJE
    write(*,*) 'NST = ', dvode_state % NST
    write(*,*) 'ICF = ', dvode_state % ICF
    write(*,*) 'IPUP = ', dvode_state % IPUP
    write(*,*) 'JCUR = ', dvode_state % JCUR
    write(*,*) 'JSTART = ', dvode_state % JSTART
    write(*,*) 'JSV = ', dvode_state % JSV
    write(*,*) 'KFLAG = ', dvode_state % KFLAG
    write(*,*) 'L = ', dvode_state % L
    write(*,*) 'MXSTEP = ', dvode_state % MXSTEP
    write(*,*) 'NEWH = ', dvode_state % NEWH
    write(*,*) 'NEWQ = ', dvode_state % NEWQ
    write(*,*) 'NQ = ', dvode_state % NQ
    write(*,*) 'NQNYH = ', dvode_state % NQNYH
    write(*,*) 'NQWAIT = ', dvode_state % NQWAIT
    write(*,*) 'NSLJ = ', dvode_state % NSLJ
    write(*,*) 'NSLP = ', dvode_state % NSLP
  end subroutine print_state
#endif

end module cuvode_types_module
