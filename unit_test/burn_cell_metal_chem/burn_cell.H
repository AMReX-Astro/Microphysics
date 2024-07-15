#ifndef BURN_CELL_H
#define BURN_CELL_H

#include <fstream>
#include <iostream>

#include <burner.H>
#include <eos.H>
#include <extern_parameters.H>
#include <network.H>

amrex::Real grav_constant = 6.674e-8;

AMREX_INLINE
auto burn_cell_c() -> int {

  using namespace unit_test_rp;

  burn_t state;

  amrex::Real numdens[NumSpec] = {-1.0};

  for (int n = 1; n <= NumSpec; ++n) {
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
      numdens[n - 1] = primary_species_17;
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
      numdens[n - 1] = primary_species_23;
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
    case 35:
      numdens[n - 1] = primary_species_35;
      break;
    case 36:
      numdens[n - 1] = primary_species_36;
      break;
    case 37:
      numdens[n - 1] = primary_species_37;
      break;
    case 38:
      numdens[n - 1] = primary_species_38;
      break;
    case 39:
      numdens[n - 1] = primary_species_39;
      break;
    }
  }

  // Echo initial conditions at burn and fill burn state input

  std::cout << "Maximum Time (s): " << tmax << std::endl;
  std::cout << "State Temperature (K): " << temperature << std::endl;
  for (int n = 0; n < NumSpec; ++n) {
    std::cout << "Number Density input (" << short_spec_names_cxx[n]
              << "): " << numdens[n] << std::endl;
  }

  state.T = temperature;

  // find the density in g/cm^3
  amrex::Real rhotot = 0.0_rt;
  for (int n = 0; n < NumSpec; ++n) {
    state.xn[n] = numdens[n];
    rhotot += state.xn[n] * spmasses[n]; // spmasses contains the masses of all
                                         // species, defined in EOS
  }

  state.rho = rhotot;

  // normalize -- just in case

  amrex::Real mfracs[NumSpec] = {-1.0};
  amrex::Real msum = 0.0_rt;
  for (int n = 0; n < NumSpec; ++n) {
    mfracs[n] = state.xn[n] * spmasses[n] / rhotot;
    msum += mfracs[n];
  }

  for (int n = 0; n < NumSpec; ++n) {
    mfracs[n] /= msum;
    // use the normalized mass fractions to obtain the corresponding number
    // densities
    state.xn[n] = mfracs[n] * rhotot / spmasses[n];
  }

#ifdef DEBUG
  for (int n = 0; n < NumSpec; ++n) {
    std::cout << "Mass fractions (" << short_spec_names_cxx[n]
              << "): " << mfracs[n] << std::endl;
    std::cout << "Number Density conserved (" << short_spec_names_cxx[n]
              << "): " << state.xn[n] << std::endl;
  }
#endif

  // call the EOS to set initial internal energy e
  eos(eos_input_rt, state);

  // name of output file
  std::ofstream state_over_time("state_over_time.txt");

  // save the initial state -- we'll use this to determine
  // how much things changed over the entire burn
  burn_t state_in = state;

  // output the data in columns, one line per timestep
  state_over_time << std::setw(10) << "# Time";
  state_over_time << std::setw(15) << "Density";
  state_over_time << std::setw(15) << "Temperature";
  for (int x = 0; x < NumSpec; ++x) {
    const std::string &element = short_spec_names_cxx[x];
    state_over_time << std::setw(15) << element;
  }
  state_over_time << std::endl;

  amrex::Real t = 0.0;

  state_over_time << std::setw(10) << t;
  state_over_time << std::setw(15) << state.rho;
  state_over_time << std::setw(15) << state.T;
  for (int x = 0; x < NumSpec; ++x) {
    state_over_time << std::setw(15) << state.xn[x];
  }
  state_over_time << std::endl;

  // loop over steps, burn, and output the current state
  // the loop below is similar to that used in krome and pynucastro
  amrex::Real dd = rhotot;
  amrex::Real dd1 = 0.0_rt;

  for (int n = 0; n < nsteps; n++) {

    dd1 = dd;

    amrex::Real rhotmp = 0.0_rt;

    for (int nn = 0; nn < NumSpec; ++nn) {
      rhotmp += state.xn[nn] * spmasses[nn];
    }

    // find the freefall time
    amrex::Real tff = std::sqrt(M_PI * 3.0 / (32.0 * rhotmp * grav_constant));
    amrex::Real dt = tff_reduc * tff;
    // scale the density
    dd += dt * (dd / tff);

#ifdef DEBUG
    std::cout << "step params " << dd << ", " << tff << ", " << rhotmp
              << std::endl;
#endif

    // stop the test if dt is very small
    if (dt < 10) {
      break;
    }

    // stop the test if we have reached very high densities
    if (dd > 2e-6) {
      break;
    }

    std::cout << "step " << n << " done" << std::endl;

    // scale the number densities
    for (int nn = 0; nn < NumSpec; ++nn) {
      state.xn[nn] *= dd / dd1;
    }

    // input the scaled density in burn state
    state.rho *= dd / dd1;

    // do the actual integration
    burner(state, dt);

    // ensure positivity and normalize
    amrex::Real inmfracs[NumSpec] = {-1.0};
    amrex::Real insum = 0.0_rt;
    for (int nn = 0; nn < NumSpec; ++nn) {
      state.xn[nn] = amrex::max(state.xn[nn], small_x);
      inmfracs[nn] = spmasses[nn] * state.xn[nn] / state.rho;
      insum += inmfracs[nn];
    }

    for (int nn = 0; nn < NumSpec; ++nn) {
      inmfracs[nn] /= insum;
      // update the number densities with conserved mass fractions
      state.xn[nn] = inmfracs[nn] * state.rho / spmasses[nn];
    }

    // update the number density of electrons due to charge conservation
    state.xn[0] = -state.xn[3] - state.xn[7] - state.xn[17] - state.xn[25] +
                  state.xn[1] + state.xn[6] + state.xn[13] + state.xn[4] + state.xn[9] + 
                  2.0 * state.xn[12] + state.xn[15] + state.xn[23] + state.xn[34] + state.xn[35] +
                  state.xn[11] + state.xn[18] + state.xn[20] + state.xn[22] + state.xn[32] +
                  state.xn[26] + state.xn[28] + state.xn[30] + state.xn[37];

    // reconserve mass fractions post charge conservation
    insum = 0;
    for (int nn = 0; nn < NumSpec; ++nn) {
      state.xn[nn] = amrex::max(state.xn[nn], small_x);
      inmfracs[nn] = spmasses[nn] * state.xn[nn] / state.rho;
      insum += inmfracs[nn];
    }

    for (int nn = 0; nn < NumSpec; ++nn) {
      inmfracs[nn] /= insum;
      // update the number densities with conserved mass fractions
      state.xn[nn] = inmfracs[nn] * state.rho / spmasses[nn];
    }

    // get the updated T
    eos(eos_input_re, state);

    t += dt;

    state_over_time << std::setw(10) << t;
    state_over_time << std::setw(15) << state.rho;
    state_over_time << std::setw(12) << state.T;
    for (int x = 0; x < NumSpec; ++x) {
      state_over_time << std::setw(15) << state.xn[x];
    }
    state_over_time << std::endl;
  }
  state_over_time.close();

  // output diagnostics to the terminal
  std::cout << "------------------------------------" << std::endl;
  std::cout << "successful? " << state.success << std::endl;

  std::cout << "------------------------------------" << std::endl;
  std::cout << "T initial = " << state_in.T << std::endl;
  std::cout << "T final =   " << state.T << std::endl;
  std::cout << "Eint initial = " << state_in.e << std::endl;
  std::cout << "Eint final =   " << state.e << std::endl;
  std::cout << "rho initial = " << state_in.rho << std::endl;
  std::cout << "rho final =   " << state.rho << std::endl;

  std::cout << "------------------------------------" << std::endl;
  std::cout << "New number densities: " << std::endl;
  for (int n = 0; n < NumSpec; ++n) {
    std::cout << state.xn[n] << std::endl;
  }

  return state.success;
}
#endif
