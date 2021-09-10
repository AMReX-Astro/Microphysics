n_cell = 16

do_cxx = 1

prefix = react_aprox19_sdc_

unit_test.dens_min   = 1.e5
unit_test.dens_max   = 5.e8
unit_test.temp_min   = 5.e7
unit_test.temp_max   = 5.e9

unit_test.C_nse = 0.2
unit_test.T_nse = 5.e8
unit_test.rho_nse = 3.e7

# NSE doesn't really care about the timestep
unit_test.tmax = 1.e-9

unit_test.primary_species_1 = helium-4
unit_test.primary_species_2 = iron-54

integrator.jacobian = 1


