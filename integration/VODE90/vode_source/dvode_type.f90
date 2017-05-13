module dvode_type_module

  implicit none

  type :: dvode_t
     real(dp_t) :: HU
     real(dp_t) :: ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13)
     real(dp_t) :: ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1
     real(dp_t) :: RC, RL1, TAU(13), TQ(5), TN, UROUND
     integer    :: NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST     
     integer    :: ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH
     integer    :: L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM
     integer    :: LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP
     integer    :: N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ
     integer    :: NSLP, NYH
  end type dvode_t

end module dvode_type_module
