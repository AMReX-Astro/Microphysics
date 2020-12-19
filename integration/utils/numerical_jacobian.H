#ifndef NUMERICAL_JACOBIAN_H
#define NUMERICAL_JACOBIAN_H

#include <network.H>
#include <burn_type.H>
#include <extern_parameters.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void numerical_jac(burn_t& state, JacNetArray2D& jac)
{


    // the choice of eps should be ~ sqrt(eps), where eps is machine epsilon.
    // this balances truncation vs. roundoff error in the differencing

    constexpr Real eps = 1.e-8_rt;

    for (int i = 1; i <= neqs; i++) {
        for (int j = 1; j <= neqs; j++) {
            jac(i,j) = 0.0_rt;
        }
    }

    if (centered_diff_jac) {
        burn_t state_delp = state;
        burn_t state_delm = state;

        YdotNetArray1D ydotp;
        YdotNetArray1D ydotm;

        // species derivatives

        for (int n = 1; n <= NumSpec; n++) {
            // perturb species

            for (int q = 0; q < NumSpec; q++) {
                state_delp.xn[q] = state.xn[q];
            }

            // xn is 0-based
            state_delp.xn[n-1] *= (1.0_rt + eps);

            actual_rhs(state_delp, ydotp);

            // We integrate X, so convert from the Y we got back from the RHS

            for (int q = 1; q <= NumSpec; q++) {
                // note ydotp is 1-based, but aion is 0-based
                ydotp(q) *= aion[q-1];
            }

            for (int q = 0; q < NumSpec; q++) {
                state_delm.xn[q] = state.xn[q];
            }

            // xn is 0-based
            state_delm.xn[n-1] *= (1.0_rt - eps);

            actual_rhs(state_delm, ydotm);

            for (int q = 1; q <= NumSpec; q++) {
                // note ydotm is 1-based, but aion is 0-based
                ydotm(q) *= aion[q-1];
            }

            for (int m = 1; m <= neqs; m++) {
                jac(m,n) = 0.5_rt * (ydotp(m) - ydotm(m)) / (eps * state.xn[n-1]);
            }
        }


        // temperature derivative

        for (int q = 0; q < NumSpec; q++) {
            state_delp.xn[q] = state.xn[q];
        }
        state_delp.T *= (1.0_rt + eps);

        actual_rhs(state_delp, ydotp);

        for (int q = 1; q <= NumSpec; q++) {
            ydotp(q) *= aion[q-1];
        }
%%%%%%%

       state_delm % xn = state % xn
       state_delm % T  = state % T * (ONE - eps)

       call actual_rhs(state_delm, ydotm)

       ydotm(1:nspec) = ydotm(1:nspec) * aion(1:nspec)

       do m = 1, neqs
          scratch = HALF*(ydotp(m) - ydotm(m)) / (eps * state % T)
          call set_jac_entry(jac, m, net_itemp, scratch)
       enddo

       ! energy derivatives -- these are all 0! (yay!)
       scratch = ZERO
       do m = 1, neqs
          call set_jac_entry(jac, m, net_ienuc, scratch)
       enddo

    else
       call copy_burn_t(state_delp, state)
       call copy_burn_t(state_delm, state)   ! note: delm here is actually just the input

       ! default
       call actual_rhs(state_delm, ydotm)

       ydotm(1:nspec) = ydotm(1:nspec) * aion(1:nspec)

       ! species derivatives
       do n = 1, nspec
          ! perturb species -- we send in X, but ydot is in terms of dY/dt, not dX/dt
          state_delp % xn = state % xn

          h = eps * abs(state % xn(n))
          if (h == 0) then
             h = eps
          endif

          state_delp % xn(n) = state % xn(n) + h

          call actual_rhs(state_delp, ydotp)

          ! We integrate X, so convert from the Y we got back from the RHS

          ydotp(1:nspec) = ydotp(1:nspec) * aion(1:nspec)

          do m = 1, neqs
             scratch = (ydotp(m) - ydotm(m)) / h
             call set_jac_entry(jac, m, n, scratch)
          enddo
       enddo

       ! temperature derivative
       state_delp % xn = state % xn

       h = eps * abs(state % T)
       if (h == 0) then
          h = eps
       endif

       state_delp % T = state % T + h

       call actual_rhs(state_delp, ydotp)

       ydotp(1:nspec) = ydotp(1:nspec) * aion(1:nspec)

       do m = 1, neqs
          scratch = (ydotp(m) - ydotm(m)) / h
          call set_jac_entry(jac, m, net_itemp, scratch)
       enddo

       ! energy derivatives -- these are all 0! (yay!)
       scratch = ZERO
       do m = 1, neqs
          call set_jac_entry(jac, m, net_ienuc, scratch)
       enddo

    endif

  end subroutine numerical_jac

#endif
