module cuvode_dvjust_module

  use cuvode_parameters_module, only: VODE_LMAX, VODE_NEQS, VODE_MAXORD
  use cuvode_types_module, only: dvode_t
  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  subroutine dvjust(IORD, vstate)

    !$acc routine seq
    
    ! -----------------------------------------------------------------------
    !  This subroutine adjusts the YH array on reduction of order,
    !  and also when the order is increased.
    !  Communication with DVJUST uses the following:
    !  IORD  = An integer flag used to indicate an order
    !          increase (IORD = +1) or an order decrease (IORD = -1).
    !  HSCAL = Step size H used in scaling of Nordsieck array YH.
    !          (If IORD = +1, DVJUST assumes that HSCAL = TAU(1).)
    !  See References 1 and 2 for details.
    ! -----------------------------------------------------------------------

    implicit none

    ! Declare arguments
    type(dvode_t), intent(inout) :: vstate
    integer,       intent(in   ) :: IORD

    ! Declare local variables
    real(rt) :: ALPH0, ALPH1, HSUM, PROD, T1, XI,XIOLD
    integer    :: I, IBACK, J, JP1, LP1, NQM1, NQM2, NQP1

    !$gpu

    IF ((vstate % NQ .EQ. 2) .AND. (IORD .NE. 1)) RETURN

    NQM1 = vstate % NQ - 1
    NQM2 = vstate % NQ - 2

    ! -----------------------------------------------------------------------
    !  Stiff option...
    !  Check to see if the order is being increased or decreased.
    ! -----------------------------------------------------------------------

    IF (IORD /= 1) then

       ! Order decrease. ------------------------------------------------------
       do J = 1, VODE_LMAX
          vstate % EL(J) = 0.0_rt
       end do
       vstate % EL(3) = 1.0_rt
       HSUM = 0.0_rt
       do J = 1,NQM2
          ! Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
          HSUM = HSUM + vstate % TAU(J)
          XI = HSUM/vstate % HSCAL
          JP1 = J + 1
          do IBACK = 1, JP1
             I = (J + 4) - IBACK
             vstate % EL(I) = vstate % EL(I) * XI + vstate % EL(I-1)
          end do
       end do
       ! Subtract correction terms from YH array. -----------------------------
       do J = 3,vstate % NQ
          do I = 1, VODE_NEQS
             vstate % YH(I,J) = vstate % YH(I,J) - &
                  vstate % YH(I,vstate % L) * vstate % EL(J)
          end do
       end do

    else

       ! Order increase. ------------------------------------------------------
       do J = 1, VODE_LMAX
          vstate % EL(J) = 0.0_rt
       end do
       vstate % EL(3) = 1.0_rt
       ALPH0 = -1.0_rt
       ALPH1 = 1.0_rt
       PROD = 1.0_rt
       XIOLD = 1.0_rt
       HSUM = vstate % HSCAL
       IF (vstate % NQ /= 1) then
          do J = 1, NQM1
             ! Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
             JP1 = J + 1
             HSUM = HSUM + vstate % TAU(JP1)
             XI = HSUM/vstate % HSCAL
             PROD = PROD*XI
             ALPH0 = ALPH0 - 1.0_rt/REAL(JP1)
             ALPH1 = ALPH1 + 1.0_rt/XI
             do IBACK = 1, JP1
                I = (J + 4) - IBACK
                vstate % EL(I) = vstate % EL(I) * XIOLD + vstate % EL(I-1)
             end do
             XIOLD = XI
          end do
       end IF

       T1 = (-ALPH0 - ALPH1)/PROD
       ! Load column L + 1 in YH array. ---------------------------------------
       LP1 = vstate % L + 1
       do I = 1, VODE_NEQS
          vstate % YH(I,LP1) = T1 * vstate % YH(I,VODE_LMAX)
       end do
       ! Add correction terms to YH array. ------------------------------------
       NQP1 = vstate % NQ + 1
       do J = 3, NQP1
          vstate % YH(:,J) = vstate % YH(:,J) + vstate % EL(J) * vstate % YH(:,LP1)
       end do

    end if
    
    return

  end subroutine dvjust

end module cuvode_dvjust_module
