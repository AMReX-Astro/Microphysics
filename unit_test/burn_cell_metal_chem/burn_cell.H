#ifndef BURN_CELL_H
#define BURN_CELL_H

#include <fstream>
#include <iostream>

#include <burner.H>
#include <eos.H>
#include <extern_parameters.H>
#include <network.H>

amrex::Real grav_constant = C::Gconst;

AMREX_INLINE
auto burn_cell_c() -> int {

  using namespace unit_test_rp;

  burn_t state;

  amrex::Real numdens[NumSpec] = {-1.0};
  int n;
  int x;
  int nn;

  amrex::Real metal = network_rp::metallicity;

  for (n = 1; n <= NumSpec; ++n) {
    switch (n) {

    case 1:
      numdens[n - 1] = primary_species_1;
      break;
    case 2:
      numdens[n - 1] = primary_species_2;
      break;
    case 3:
      numdens[n - 1] = primary_species_3;
      break;
    case 4:
      numdens[n - 1] = primary_species_4;
      break;
    case 5:
      numdens[n - 1] = primary_species_5;
      break;
    case 6:
      numdens[n - 1] = primary_species_6;
      break;
    case 7:
      numdens[n - 1] = primary_species_7;
      break;
    case 8:
      numdens[n - 1] = primary_species_8;
      break;
    case 9:
      numdens[n - 1] = primary_species_9;
      break;
    case 10:
      numdens[n - 1] = primary_species_10;
      break;
    case 11:
      numdens[n - 1] = primary_species_11;
      break;
    case 12:
      numdens[n - 1] = primary_species_12;
      break;
    case 13:
      numdens[n - 1] = primary_species_13;
      break;
    case 14:
      numdens[n - 1] = primary_species_14;
      break;
    case 15:
      numdens[n - 1] = primary_species_15;
      break;
    case 16:
      numdens[n - 1] = primary_species_16;
      break;
    case 17:
      numdens[n - 1] = primary_species_17*metal;
      break;
    case 18:
      numdens[n - 1] = primary_species_18;
      break;
    case 19:
      numdens[n - 1] = primary_species_19;
      break;
    case 20:
      numdens[n - 1] = primary_species_20;
      break;
    case 21:
      numdens[n - 1] = primary_species_21;
      break;
    case 22:
      numdens[n - 1] = primary_species_22;
      break;
    case 23:
      numdens[n - 1] = primary_species_23*metal;
      break;
    case 24:
      numdens[n - 1] = primary_species_24;
      break;
    case 25:
      numdens[n - 1] = primary_species_25;
      break;
    case 26:
      numdens[n - 1] = primary_species_26;
      break;
    case 27:
      numdens[n - 1] = primary_species_27;
      break;
    case 28:
      numdens[n - 1] = primary_species_28;
      break;
    case 29:
      numdens[n - 1] = primary_species_29;
      break;
    case 30:
      numdens[n - 1] = primary_species_30;
      break;
    case 31:
      numdens[n - 1] = primary_species_31;
      break;
    case 32:
      numdens[n - 1] = primary_species_32;
      break;
    case 33:
      numdens[n - 1] = primary_species_33;
      break;
    case 34:
      numdens[n - 1] = primary_species_34;
      break;
    }
  }

  //scale number densities by initial ninit
  for (n = 0; n < NumSpec; ++n) {
    numdens[n] *= ninit;
  }

  //if metallicity is 0, reset metal number densities to 0
  if (metal == 0) {
    for (n = 0; n < NumSpec; ++n) {
      if ((n < 2) || (n > 15)) {
        state.xn[n] = 0.0;
      }
    }
  }


  // Echo initial conditions at burn and fill burn state input

  std::cout << "Redshift: " << network_rp::redshift << std::endl;
  std::cout << "Metallicity: " << metal << std::endl;
  std::cout << "Dust2gas Ratio: " << network_rp::dust2gas_ratio << std::endl;
  std::cout << " " << std::endl;

  std::cout << "Maximum Time (s): " << tmax << std::endl;
  std::cout << "State Temperature (K): " << temperature << std::endl;
  for (n = 0; n < NumSpec; ++n) {
    std::cout << "Number Density input (" << short_spec_names_cxx[n]
              << "): " << numdens[n] << std::endl;
  }

  state.T = temperature;
  // set initial Tdust to CMB
  state.aux[0] = 2.73*(1.0 + network_rp::redshift);

  // find the density in g/cm^3
  amrex::Real rhotot = 0.0_rt;
  amrex::Real sum_numdens = 0.0_rt;
  for (n = 0; n < NumSpec; ++n) {
    state.xn[n] = numdens[n];
    rhotot += state.xn[n] * spmasses[n]; // spmasses contains the masses of all
                                         // species, defined in EOS
    sum_numdens += state.xn[n];
  }

  state.rho = rhotot;
  std::cout << "rho: " << rhotot << std::endl;

  // call the EOS to set initial internal energy e
  eos(eos_input_rt, state);

  std::cout << "initial eint: " << state.e << std::endl;

  // name of output file
  std::ofstream state_over_time("state_over_time.txt");

  // save the initial state -- we'll use this to determine
  // how much things changed over the entire burn
  burn_t state_in = state;

  // output the data in columns, one line per timestep
  state_over_time << std::setw(10) << "# Time";
  state_over_time << std::setw(15) << "NumberDensity";
  state_over_time << std::setw(15) << "Density";
  state_over_time << std::setw(15) << "Tgas";
  state_over_time << std::setw(15) << "Tdust";
  for (x = 0; x < NumSpec; ++x) {
    const std::string &element = short_spec_names_cxx[x];
    state_over_time << std::setw(15) << element;
  }
  state_over_time << std::endl;

  amrex::Real t = 0.0;

  state_over_time << std::setw(10) << t;
  state_over_time << std::setw(15) << sum_numdens;
  state_over_time << std::setw(15) << state.rho;
  state_over_time << std::setw(15) << state.T;
  state_over_time << std::setw(15) << state.aux[0];
  for (x = 0; x < NumSpec; ++x) {
    state_over_time << std::setw(15) << state.xn[x];
  }
  state_over_time << std::endl;

  // loop over steps, burn, and output the current state
  // the loop below is similar to that used in krome and GPUAstroChem
  amrex::Real dd = rhotot;
  amrex::Real dd1 = 0.0_rt;

  for (n = 0; n < nsteps; n++) {

    dd1 = dd;

    amrex::Real rhotmp = 0.0_rt;

    for (nn = 0; nn < NumSpec; ++nn) {
      rhotmp += state.xn[nn] * spmasses[nn];
    }

    // find the freefall time
    amrex::Real tff = std::sqrt(M_PI * 3.0 / (32.0 * rhotmp * grav_constant));
    amrex::Real dt = tff_reduc * tff;
    // scale the density
    dd += dt * (dd / tff);

    //std::cout << "step params " << dd << ", " << tff << ", " << tff_reduc << ", " << rhotmp << ", " << state.xn << std::endl;

    // stop the test if dt is very small
    if (dt < 10) {
      break;
    }

    // stop the test if we have reached very high densities
    if (dd > 2e-6) {
      break;
    }

    // scale the number densities
    for (nn = 0; nn < NumSpec; ++nn) {
      state.xn[nn] *= dd / dd1;
    }

    // input the scaled density in burn state
    state.rho *= dd / dd1;

    //std::cout << "before burn: " << state.rho << ", " << state.T << ", " << state.xn << ", " << state.e << std::endl;

    integrator_rp::ode_max_dt = dt*1e-1;
    // do the actual integration
    burner(state, dt);

    //std::cout << "after burn: " << state << std::endl;

    // ensure positivity and normalize
    //amrex::Real inmfracs[NumSpec] = {-1.0};
    //amrex::Real insum = 0.0_rt;
    //for (int nn = 0; nn < NumSpec; ++nn) {
    //  state.xn[nn] = amrex::max(state.xn[nn], small_x);
    //  inmfracs[nn] = spmasses[nn] * state.xn[nn] / state.rho;
    //  insum += inmfracs[nn];
    //}

    //for (int nn = 0; nn < NumSpec; ++nn) {
    //  inmfracs[nn] /= insum;
    //  // update the number densities with conserved mass fractions
    //  state.xn[nn] = inmfracs[nn] * state.rho / spmasses[nn];
    //}

    // update the number density of electrons due to charge conservation
    state.xn[2] = -state.xn[5] - state.xn[9] + state.xn[3] + state.xn[6] + state.xn[8] +
                  2.0 * state.xn[13] + state.xn[11] + state.xn[14] + state.xn[16] + state.xn[21] +
                  state.xn[24] + state.xn[26] + state.xn[28] + state.xn[29] + state.xn[31];

    //std::cout << "conserved elec: " << state.xn[2] << std::endl;

    // reconserve mass fractions post charge conservation
    //insum = 0;
    //for (int nn = 0; nn < NumSpec; ++nn) {
    //  state.xn[nn] = amrex::max(state.xn[nn], small_x);
    //  inmfracs[nn] = spmasses[nn] * state.xn[nn] / state.rho;
    //  insum += inmfracs[nn];
    //}

    //for (int nn = 0; nn < NumSpec; ++nn) {
    //  inmfracs[nn] /= insum;
    //  // update the number densities with conserved mass fractions
    //  state.xn[nn] = inmfracs[nn] * state.rho / spmasses[nn];
    //}

    // get the updated T
    eos(eos_input_re, state);

    t += dt;

    // get number density
    // make abundance 0 for all metals if metallicity is 0
    for (nn = 0; nn < NumSpec; ++nn) {
      if (metal == 0 && ((nn < 2) || (nn > 15))) {
        state.xn[nn] = 0.0;
      }
      sum_numdens += state.xn[nn];
    }

    state_over_time << std::setw(10) << t;
    state_over_time << std::setw(15) << sum_numdens;
    state_over_time << std::setw(15) << state.rho;
    state_over_time << std::setw(12) << state.T;
    state_over_time << std::setw(12) << state.aux[0];
    for (x = 0; x < NumSpec; ++x) {
      state_over_time << std::setw(15) << state.xn[x];
    }
    state_over_time << std::endl;
    std::cout << "step " << n << " done with dt = " << dt << std::endl;

  }
  state_over_time.close();

  // output diagnostics to the terminal
  std::cout << "------------------------------------" << std::endl;
  std::cout << "successful? " << state.success << std::endl;

  std::cout << "------------------------------------" << std::endl;
  std::cout << "T initial = " << state_in.T << std::endl;
  std::cout << "T final =   " << state.T << std::endl;
  std::cout << "Tdust initial = " << state_in.aux[0] << std::endl;
  std::cout << "Tdust final =   " << state.aux[0] << std::endl;
  std::cout << "Eint initial = " << state_in.e << std::endl;
  std::cout << "Eint final =   " << state.e << std::endl;
  std::cout << "rho initial = " << state_in.rho << std::endl;
  std::cout << "rho final =   " << state.rho << std::endl;

  std::cout << "------------------------------------" << std::endl;
  std::cout << "New number densities: " << std::endl;
  for (n = 0; n < NumSpec; ++n) {
    std::cout << state.xn[n] << std::endl;
  }

  return state.success;
}
#endif
