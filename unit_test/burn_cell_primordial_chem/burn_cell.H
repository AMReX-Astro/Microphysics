#ifndef BURN_CELL_H
#define BURN_CELL_H

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include <AMReX_GpuAsyncArray.H>
#include <AMReX_Gpu.H>
#include <actual_network.H>
#include <burner.H>
#include <eos.H>
#include <extern_parameters.H>
#include <network.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
auto balance_charge_local(burn_t& state) -> void {
  state.xn[0] = -state.xn[3] - state.xn[7] + state.xn[1] + state.xn[12] +
                state.xn[6] + state.xn[4] + state.xn[9] +
                2.0_rt * state.xn[11];
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
auto collapse_perturbation_factor(const amrex::Long cell_idx,
                                  const int step) -> amrex::Real {
  std::uint64_t x =
      static_cast<std::uint64_t>(cell_idx) + 0x9e3779b97f4a7c15ULL;
  x ^= static_cast<std::uint64_t>(step) + 0xbf58476d1ce4e5b9ULL +
       (x << 6U) + (x >> 2U);
  x ^= x >> 30U;
  x *= 0xbf58476d1ce4e5b9ULL;
  x ^= x >> 27U;
  x *= 0x94d049bb133111ebULL;
  x ^= x >> 31U;

  constexpr amrex::Real inv_two_to_53 =
      1.0_rt / 9007199254740992.0_rt;
  return 0.5_rt + static_cast<amrex::Real>(x >> 11U) * inv_two_to_53;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
auto initialize_burn_state() -> burn_t {

  burn_t state;

  amrex::Real numdens[NumSpec] = {-1.0};

  for (int n = 1; n <= NumSpec; ++n) {
    switch (n) {

    case 1:
      numdens[n - 1] = unit_test_rp::primary_species_1;
      break;
    case 2:
      numdens[n - 1] = unit_test_rp::primary_species_2;
      break;
    case 3:
      numdens[n - 1] = unit_test_rp::primary_species_3;
      break;
    case 4:
      numdens[n - 1] = unit_test_rp::primary_species_4;
      break;
    case 5:
      numdens[n - 1] = unit_test_rp::primary_species_5;
      break;
    case 6:
      numdens[n - 1] = unit_test_rp::primary_species_6;
      break;
    case 7:
      numdens[n - 1] = unit_test_rp::primary_species_7;
      break;
    case 8:
      numdens[n - 1] = unit_test_rp::primary_species_8;
      break;
    case 9:
      numdens[n - 1] = unit_test_rp::primary_species_9;
      break;
    case 10:
      numdens[n - 1] = unit_test_rp::primary_species_10;
      break;
    case 11:
      numdens[n - 1] = unit_test_rp::primary_species_11;
      break;
    case 12:
      numdens[n - 1] = unit_test_rp::primary_species_12;
      break;
    case 13:
      numdens[n - 1] = unit_test_rp::primary_species_13;
      break;
    case 14:
      numdens[n - 1] = unit_test_rp::primary_species_14;
      break;
    }
  }

  state.T = unit_test_rp::temperature;

  // find the density in g/cm^3
  amrex::Real rhotot = 0.0_rt;
  for (int n = 0; n < NumSpec; ++n) {
    state.xn[n] = numdens[n];
    rhotot += state.xn[n] * spmasses[n];
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

#if defined(DEBUG) && !defined(AMREX_USE_GPU)
  for (int n = 0; n < NumSpec; ++n) {
    std::cout << "Mass fractions (" << short_spec_names_cxx[n]
              << "): " << mfracs[n] << std::endl;
    std::cout << "Number Density conserved (" << short_spec_names_cxx[n]
              << "): " << state.xn[n] << std::endl;
  }
#endif

  // call the EOS to set initial internal energy e
  eos(eos_input_rt, state);

  return state;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
auto evolve_burn_cell(burn_t& state, const amrex::Long cell_idx) -> int {

  // loop over steps, burn, and output the current state
  // the loop below is similar to that used in krome and pynucastro
  amrex::Real dd = state.rho;
  amrex::Real dd1 = 0.0_rt;

  for (int n = 0; n < unit_test_rp::nsteps; n++) {

    dd1 = dd;

    amrex::Real rhotmp = 0.0_rt;

    for (int nn = 0; nn < NumSpec; ++nn) {
      rhotmp += state.xn[nn] * spmasses[nn];
    }

    // find the freefall time
    constexpr amrex::Real grav_constant = 6.674e-8_rt;
    amrex::Real tff = std::sqrt(M_PI * 3.0 / (32.0 * rhotmp * grav_constant));
    amrex::Real dt = unit_test_rp::tff_reduc * tff;
    // scale the density
    dd += dt * (dd / tff);

#if defined(DEBUG) && !defined(AMREX_USE_GPU)
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

    // scale the number densities
    for (int nn = 0; nn < NumSpec; ++nn) {
      state.xn[nn] *= dd / dd1;
    }

    // input the scaled density in burn state
    state.rho *= dd / dd1;

    if ((n + 1) % 20 == 0) {
      const amrex::Real perturbation =
          collapse_perturbation_factor(cell_idx, n + 1);
      dd *= perturbation;
      state.rho *= perturbation;
      for (int nn = 0; nn < NumSpec; ++nn) {
        state.xn[nn] *= perturbation;
      }
    }

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
    balance_charge_local(state);

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
  }

  return state.success;
}

AMREX_INLINE
auto evolve_burn_cell_with_history(burn_t state) -> int {

  std::ofstream state_over_time;
  state_over_time.open("state_over_time.txt");

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
  amrex::Real dd = state.rho;
  amrex::Real dd1 = 0.0_rt;

  for (int n = 0; n < unit_test_rp::nsteps; n++) {

    dd1 = dd;

    amrex::Real rhotmp = 0.0_rt;

    for (int nn = 0; nn < NumSpec; ++nn) {
      rhotmp += state.xn[nn] * spmasses[nn];
    }

    // find the freefall time
    constexpr amrex::Real grav_constant = 6.674e-8_rt;
    amrex::Real tff = std::sqrt(M_PI * 3.0 / (32.0 * rhotmp * grav_constant));
    amrex::Real dt = unit_test_rp::tff_reduc * tff;
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

    if ((n + 1) % 20 == 0) {
      const amrex::Real perturbation =
          collapse_perturbation_factor(0, n + 1);
      dd *= perturbation;
      state.rho *= perturbation;
      for (int nn = 0; nn < NumSpec; ++nn) {
        state.xn[nn] *= perturbation;
      }
    }

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
    balance_charge_local(state);

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

AMREX_INLINE
auto burn_cell_c(const int ncell) -> int {

  std::int64_t total_cells = 1;
  total_cells *= ncell;
  total_cells *= ncell;
  total_cells *= ncell;

  std::cout << "integrating cubic grid: " << ncell << " x " << ncell << " x "
            << ncell << " cells (" << total_cells << " total)" << std::endl;

  std::cout << "All cells use the same initial thermodynamic state and "
               "composition"
            << std::endl;

  int num_failed = 0;
  amrex::AsyncArray<int> aa_num_failed(&num_failed, 1);
  int* num_failed_d = aa_num_failed.data();

  amrex::ParallelFor(total_cells, [=] AMREX_GPU_DEVICE (amrex::Long idx) noexcept {
    const int i = static_cast<int>(idx % ncell);
    const int j = static_cast<int>((idx / ncell) % ncell);
    const int k = static_cast<int>(idx / (ncell * ncell));

    burn_t state = initialize_burn_state();
    state.i = i;
    state.j = j;
    state.k = k;

    if (!evolve_burn_cell(state, idx)) {
      amrex::Gpu::Atomic::Add(num_failed_d, 1);
    }
  });

  aa_num_failed.copyToHost(&num_failed, 1);
  amrex::Gpu::synchronize();

  if (num_failed > 0) {
    std::cout << num_failed << " cell integrations failed" << std::endl;
    return 0;
  }

  std::cout << "completed " << total_cells << " cell integrations"
            << std::endl;

  std::cout << "integrating representative cell (0, 0, 0) for diagnostics..."
            << std::endl;
  std::cout << "Maximum Time (s): " << unit_test_rp::tmax << std::endl;
  std::cout << "State Temperature (K): " << unit_test_rp::temperature
            << std::endl;

  burn_t state = initialize_burn_state();
  for (int n = 0; n < NumSpec; ++n) {
    std::cout << "Number Density input (" << short_spec_names_cxx[n]
              << "): " << state.xn[n] << std::endl;
  }

  return evolve_burn_cell_with_history(state);
}

#endif
