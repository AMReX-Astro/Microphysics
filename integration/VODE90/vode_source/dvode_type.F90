module dvode_type_module

  use bl_types, only: dp_t
  use vode_parameters_module, only: VODE_NEQS
  use rpar_indices, only: n_rpar_comps

  use dvode_constants_module
  
  implicit none

  type :: dvode_t
     ! Variables previously in common blocks
     real(dp_t) :: HU
     real(dp_t) :: ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13)
     real(dp_t) :: ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1
     real(dp_t) :: RC, RL1, TAU(13), TQ(5), TN, UROUND
     integer    :: NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST     
     integer    :: ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH
     integer    :: L, LENWM
     integer    :: LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP
     integer    :: NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ
     integer    :: NSLP

     ! Tolerances
     real(dp_t) :: RTOL(VODE_NEQS), ATOL(VODE_NEQS)

     ! Real parameters
     real(dp_t) :: RPAR(n_rpar_comps)

     ! State flag
     integer    :: ISTATE

     ! Local time and integration end time
     real(dp_t) :: T, TOUT

     ! Integration vector
     real(dp_t) :: Y(VODE_NEQS)
  end type dvode_t

contains

#ifndef CUDA  
  subroutine print_state(dvode_state)
    type(dvode_t) :: dvode_state

    write(*,*) 'HU = ', dvode_state % HU
    write(*,*) 'ACNRM = ', dvode_state % ACNRM
    write(*,*) 'CCMXJ = ', dvode_state % CCMXJ
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
    write(*,*) 'HMIN = ', dvode_state % HMIN
    write(*,*) 'HMXI = ', dvode_state % HMXI
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
    write(*,*) 'UROUND = ', dvode_state % UROUND
    write(*,*) 'NCFN = ', dvode_state % NCFN
    write(*,*) 'NETF = ', dvode_state % NETF
    write(*,*) 'NFE = ', dvode_state % NFE
    write(*,*) 'NJE = ', dvode_state % NJE
    write(*,*) 'NLU = ', dvode_state % NLU
    write(*,*) 'NNI = ', dvode_state % NNI
    write(*,*) 'NQU = ', dvode_state % NQU
    write(*,*) 'NST = ', dvode_state % NST
    write(*,*) 'ICF = ', dvode_state % ICF
    write(*,*) 'INIT = ', dvode_state % INIT
    write(*,*) 'IPUP = ', dvode_state % IPUP
    write(*,*) 'JCUR = ', dvode_state % JCUR
    write(*,*) 'JSTART = ', dvode_state % JSTART
    write(*,*) 'JSV = ', dvode_state % JSV
    write(*,*) 'KFLAG = ', dvode_state % KFLAG
    write(*,*) 'KUTH = ', dvode_state % KUTH
    write(*,*) 'L = ', dvode_state % L
    write(*,*) 'LENWM = ', dvode_state % LENWM
    write(*,*) 'LOCJS = ', dvode_state % LOCJS
    write(*,*) 'METH = ', dvode_state % METH
    write(*,*) 'MITER = ', dvode_state % MITER
    write(*,*) 'MSBJ = ', dvode_state % MSBJ
    write(*,*) 'MXHNIL = ', dvode_state % MXHNIL
    write(*,*) 'MXSTEP = ', dvode_state % MXSTEP
    write(*,*) 'NEWH = ', dvode_state % NEWH
    write(*,*) 'NEWQ = ', dvode_state % NEWQ
    write(*,*) 'NHNIL = ', dvode_state % NHNIL
    write(*,*) 'NQ = ', dvode_state % NQ
    write(*,*) 'NQNYH = ', dvode_state % NQNYH
    write(*,*) 'NQWAIT = ', dvode_state % NQWAIT
    write(*,*) 'NSLJ = ', dvode_state % NSLJ
    write(*,*) 'NSLP = ', dvode_state % NSLP
  end subroutine print_state
#endif
  
end module dvode_type_module
