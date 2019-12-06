module net_utils

   use microphysics_type_module

   implicit none

   integer, save :: handle_net

   integer, save :: solver_choice, decsol_choice, num_reactions, species

   integer, pointer, save :: which_rates(:)

   real(rt), dimension(:), pointer, save :: x_previous, peak_abundance, &
          peak_time, d_dxdt_dRho, d_dxdt_dT, dxdt, xin, x_initial

   real(rt), dimension(:,:), pointer, save :: d_dxdt_dx

   real(rt), save :: burn_ergs, burn_ergs_total, burn_neu_total, &
          burn_rho, burn_temp, data_output_min_t, eps_neu_total, &
          eps_nuc_total, eta, theta_e_for_graboske_et_al

   integer, save :: screening_mode

   real(rt), dimension(:), pointer, save :: category_factors, rate_factors

end module net_utils

