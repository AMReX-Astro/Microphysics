#ifndef actual_rhs_H
#define actual_rhs_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <extern_parameters.H>
#include <actual_network.H>
#include <jacobian_utilities.H>
#include <screen.H>

using namespace amrex;
using namespace ArrayUtil;


struct rate_eval_t {
    Real unscreened_rates(NumGroups, NumRates);
    Real screened_rates(NumRates);
};


// Indices into rate groups in the rate_eval_t type
const int i_rate        = 1;
const int i_drate_dt    = 2;
const int i_scor        = 3;
const int i_dscor_dt    = 4;



AMREX_GPU_HOST_DEVICE AMREX_INLINE
void evaluate_rates(const burn_t& state, rate_eval_t& rate_eval) {


    // create molar fractions

    Array1D<Real, 1, NumSpec> Y;
    for (int n = 1; n <= NumSpec; ++n) {
        Y(n) = state.xn[n-1] * aion_inv[n-1];
    }


    // Calculate Reaclib rates

    plasma_state_t pstate;
    fill_plasma_state(pstate, state.T, state.rho, Y);

//////

    do i = 1, nrat_reaclib
       call reaclib_evaluate(pstate, state % T, i, rate, drate_dt)
       rate_eval % unscreened_rates(i_rate, i) = rate
       rate_eval % unscreened_rates(i_drate_dt, i) = drate_dt
    end do

    ! Evaluate screening factors
    if (screen_reaclib) then

      call screen5(pstate, 1, scor, dscor_dt, dscor_dd)

      call screen5(pstate, 2, scor2, dscor2_dt, dscor2_dd)
      rate_eval % unscreened_rates(i_scor,1) = scor * scor2
      rate_eval % unscreened_rates(i_dscor_dt,1) = scor * dscor2_dt + dscor_dt * scor2


      call screen5(pstate, 3, scor, dscor_dt, dscor_dd)
      rate_eval % unscreened_rates(i_scor,2) = scor
      rate_eval % unscreened_rates(i_dscor_dt,2) = dscor_dt

    end if


    ! Compute screened rates
    rate_eval % screened_rates = rate_eval % unscreened_rates(i_rate, :) * &
                                 rate_eval % unscreened_rates(i_scor, :)

  end subroutine evaluate_rates


  subroutine actual_rhs(state, ydot)
    
    use extern_probin_module, only: do_constant_volume_burn, disable_thermal_neutrinos
    use burn_type_module, only: net_itemp, net_ienuc, neqs
    use sneut_module, only: sneut5
    use temperature_integration_module, only: temperature_rhs

    implicit none

    type(burn_t), intent(in) :: state
    real(rt), intent(inout) :: ydot(neqs)

    type(rate_eval_t) :: rate_eval
    real(rt) :: Y(nspec), ydot_nuc(nspec)
    integer :: i, j
    real(rt) :: rhoy, ye, enuc
    real(rt) :: sneut, dsneutdt, dsneutdd, snuda, snudz

    ! Set molar abundances
    Y(:) = state % xn(:) * aion_inv(:)

    call evaluate_rates(state, rate_eval)

    call rhs_nuc(state, ydot_nuc, Y, rate_eval % screened_rates)
    ydot(1:nspec) = ydot_nuc

    ! ion binding energy contributions
    call ener_gener_rate(ydot_nuc, enuc)

    ! include reaction neutrino losses (non-thermal)

    ! Get the thermal neutrino losses
    if (.not. disable_thermal_neutrinos) then
       call sneut5(state % T, state % rho, state % abar, state % zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)
    else
       sneut = ZERO
    end if

    ! Append the energy equation (this is erg/g/s)
    ydot(net_ienuc) = enuc - sneut

    ! Append the temperature equation
    call temperature_rhs(state, ydot)

  end subroutine actual_rhs


  subroutine rhs_nuc(state, ydot_nuc, Y, screened_rates)

    implicit none

    type (burn_t), intent(in) :: state
    real(rt), intent(out) :: ydot_nuc(nspec)
    real(rt), intent(in)  :: Y(nspec)
    real(rt), intent(in)  :: screened_rates(nrates)



    ydot_nuc(jhe4) = ( &
      -screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4)*state % rho - 0.5e0_rt* &
      screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**3*state % rho**2 &
       )

    ydot_nuc(jc12) = ( &
      -screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4)*state % rho + &
      0.16666666666666667e0_rt*screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**3* &
      state % rho**2 &
       )

    ydot_nuc(jo16) = ( &
      screened_rates(k_he4_c12__o16)*Y(jc12)*Y(jhe4)*state % rho &
       )


  end subroutine rhs_nuc


  subroutine actual_jac(state, jac)

    use burn_type_module, only: net_itemp, net_ienuc, neqs, njrows, njcols
    use extern_probin_module, only: disable_thermal_neutrinos
    use sneut_module, only: sneut5
    use temperature_integration_module, only: temperature_jac
    use jacobian_sparsity_module, only: get_jac_entry, set_jac_entry, set_jac_zero

    implicit none
    
    type(burn_t), intent(in) :: state
    real(rt), intent(inout) :: jac(njrows, njcols)

    type(rate_eval_t) :: rate_eval
    real(rt) :: screened_rates_dt(nrates)
    real(rt) :: Y(nspec), yderivs(nspec)
    real(rt) :: ye, rhoy, b1, scratch
    real(rt) :: sneut, dsneutdt, dsneutdd, snuda, snudz
    integer  :: j, k

    ! Set molar abundances
    Y(:) = state % xn(:) * aion_inv(:)
    
    call evaluate_rates(state, rate_eval)

    ! Zero out the Jacobian
    call set_jac_zero(jac)

    ! Species Jacobian elements with respect to other species
    call jac_nuc(state, jac, Y, rate_eval % screened_rates)

    ! Evaluate the species Jacobian elements with respect to temperature by
    ! calling the RHS using the temperature derivative of the screened rate
    screened_rates_dt = rate_eval % unscreened_rates(i_rate, :) * &
                        rate_eval % unscreened_rates(i_dscor_dt, :) + &
                        rate_eval % unscreened_rates(i_drate_dt, :) * &
                        rate_eval % unscreened_rates(i_scor, :)

    call rhs_nuc(state, yderivs, Y, screened_rates_dt)

    do k = 1, nspec
       call set_jac_entry(jac, k, net_itemp, yderivs(k))
    enddo

    ! Energy generation rate Jacobian elements with respect to species
    do j = 1, nspec
       do k = 1, nspec
          call get_jac_entry(jac, k, j, yderivs(k))
       enddo
       call ener_gener_rate(yderivs, scratch)
       call set_jac_entry(jac, net_ienuc, j, scratch)
    enddo

    ! Account for the thermal neutrino losses
    if (.not. disable_thermal_neutrinos) then
       call sneut5(state % T, state % rho, state % abar, state % zbar, sneut, dsneutdt, dsneutdd, snuda, snudz)

       do j = 1, nspec
          b1 = (-state % abar * state % abar * snuda + (zion(j) - state % zbar) * state % abar * snudz)
          call get_jac_entry(jac, net_ienuc, j, scratch)
          scratch = scratch - b1
          call set_jac_entry(jac, net_ienuc, j, scratch)
       enddo
    endif

    ! Energy generation rate Jacobian element with respect to temperature
    do k = 1, nspec
       call get_jac_entry(jac, k, net_itemp, yderivs(k))
    enddo
    call ener_gener_rate(yderivs, scratch)
    if (.not. disable_thermal_neutrinos) then
       scratch = scratch - dsneutdt
    endif
    call set_jac_entry(jac, net_ienuc, net_itemp, scratch)

    ! Temperature Jacobian elements
    call temperature_jac(state, jac)

  end subroutine actual_jac


  subroutine jac_nuc(state, jac, Y, screened_rates)

    use jacobian_sparsity_module, only: set_jac_entry

    implicit none

    type(burn_t), intent(in) :: state
    real(rt), intent(inout) :: jac(njrows, njcols)

    real(rt), intent(in)  :: Y(nspec)
    real(rt), intent(in)  :: screened_rates(nrates)
    real(rt) :: scratch



    scratch = (&
      -screened_rates(k_he4_c12__o16)*Y(jc12)*state % rho - 1.5e0_rt* &
      screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**2*state % rho**2 &
       )
    call set_jac_entry(jac, jhe4, jhe4, scratch)

    scratch = (&
      -screened_rates(k_he4_c12__o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jhe4, jc12, scratch)

    scratch = (&
      -screened_rates(k_he4_c12__o16)*Y(jc12)*state % rho + 0.5e0_rt* &
      screened_rates(k_he4_he4_he4__c12)*Y(jhe4)**2*state % rho**2 &
       )
    call set_jac_entry(jac, jc12, jhe4, scratch)

    scratch = (&
      -screened_rates(k_he4_c12__o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jc12, jc12, scratch)

    scratch = (&
      screened_rates(k_he4_c12__o16)*Y(jc12)*state % rho &
       )
    call set_jac_entry(jac, jo16, jhe4, scratch)

    scratch = (&
      screened_rates(k_he4_c12__o16)*Y(jhe4)*state % rho &
       )
    call set_jac_entry(jac, jo16, jc12, scratch)


  end subroutine jac_nuc

AMREX_INLINE
void actual_rhs_init () {

    init_reaclib();

    screening_init();

    set_up_screening_factors();

}


#endif
