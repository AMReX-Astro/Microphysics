module dvode_type_module

  use rpar_indices, only: n_rpar_comps, n_ipar_comps

  implicit none

  type :: dvode_t
     ! Variables previously in common blocks
     real(dp_t) :: HU
     real(dp_t) :: ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13)
     real(dp_t) :: ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1
     real(dp_t) :: RC, RL1, TAU(13), TQ(5), TN, UROUND
     integer    :: NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST     
     integer    :: ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH
     integer    :: L, LMAX, LIWM
     integer    :: LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP
     integer    :: N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ
     integer    :: NSLP, NYH
     integer    :: LYH, LEWT, LACOR, LSAVF, LWM
     
     ! ! Variables passed around the integration work routines
     ! integer    :: NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, LIW
     ! integer    :: IPAR(n_ipar_comps)
     ! integer, allocatable :: IWORK(:)
     ! real(dp_t) :: T, TOUT, RPAR(n_rpar_comps)
     ! real(dp_t), allocatable :: Y(:), RTOL(:), ATOL(:), RWORK(:)
  end type dvode_t

! contains

!   subroutine fill_state(dvode_state, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, &
!        ISTATE, IOPT, RWORK, LRW, IWORK, LIW, RPAR, IPAR)
!     type(dvode_t), intent(inout) :: dvode_state
!     integer, intent(in) :: NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, LIW
!     integer, intent(in) :: IWORK(LIW), IPAR(n_ipar_comps)
!     real(dp_t), intent(in) :: T, TOUT
!     real(dp_t), intent(in) :: Y(NEQ), RTOL(NEQ), ATOL(NEQ), RWORK(LRW), RPAR(n_rpar_comps)

!     dvode_state % NEQ = NEQ
!     dvode_state % ITOL = ITOL
!     dvode_state % ITASK = ITASK
!     dvode_state % ISTATE = ISTATE
!     dvode_state % IOPT = IOPT
!     dvode_state % LRW = LRW
!     dvode_state % LIW = LIW
!     dvode_state % T = T
!     dvode_state % TOUT = TOUT
!     dvode_state % IPAR(:) = IPAR(:)
!     dvode_state % RPAR(:) = RPAR(:)
    
!     allocate(dvode_state % IWORK(LIW))
!     dvode_state % IWORK(:) = IWORK(:)
!     allocate(dvode_state % Y(NEQ))
!     dvode_state % Y(:) = Y(:)
!     allocate(dvode_state % RTOL(NEQ))
!     dvode_state % RTOL(:) = RTOL(:)
!     allocate(dvode_state % ATOL(NEQ))
!     dvode_state % ATOL(:) = ATOL(:)
!     allocate(dvode_state % RWORK(LRW))
!     dvode_state % RWORK(:) = RWORK(:)
!   end subroutine fill_state

!   subroutine destroy_state(dvode_state)
!     type(dvode_t), intent(inout) :: dvode_state

!     deallocate(dvode_state % IWORK(LIW))
!     deallocate(dvode_state % Y(NEQ))
!     deallocate(dvode_state % RTOL(NEQ))
!     deallocate(dvode_state % ATOL(NEQ))
!     deallocate(dvode_state % RWORK(LRW))
!   end subroutine destroy_state
  
end module dvode_type_module
