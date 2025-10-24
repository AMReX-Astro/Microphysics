#ifndef BURN_CELL_H
#define BURN_CELL_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include <AMReX_GpuContainers.H>
#include <AMReX_GpuDevice.H>
#include <burner.H>
#include <eos.H>
#include <extern_parameters.H>
#include <network.H>
#include <actual_network.H>

amrex::Real grav_constant = 6.674e-8;

namespace {

AMREX_INLINE
void init_burn_state(burn_t &state, amrex::Real density_scale,
                     amrex::Real temp_offset, bool print_initial) {

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

  for (int n = 0; n < NumSpec; ++n) {
    numdens[n] *= density_scale;
  }

  if (print_initial) {
    std::cout << "Maximum Time (s): " << unit_test_rp::tmax << std::endl;
    std::cout << "State Temperature (K): "
              << unit_test_rp::temperature + temp_offset << std::endl;
    for (int n = 0; n < NumSpec; ++n) {
      std::cout << "Number Density input (" << short_spec_names_cxx[n]
                << "): " << numdens[n] << std::endl;
    }
  }

  state = burn_t{};
  state.T = unit_test_rp::temperature + temp_offset;

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

#ifdef DEBUG
  if (print_initial) {
    for (int n = 0; n < NumSpec; ++n) {
      std::cout << "Mass fractions (" << short_spec_names_cxx[n]
                << "): " << mfracs[n] << std::endl;
      std::cout << "Number Density conserved (" << short_spec_names_cxx[n]
                << "): " << state.xn[n] << std::endl;
    }
  }
#endif

  eos(eos_input_rt, state);
}

AMREX_INLINE
void enforce_post_burn_constraints(burn_t &state) {

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
  balance_charge(state);

  // reconserve mass fractions post charge conservation
  insum = 0.0_rt;
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

} // namespace

AMREX_INLINE
auto burn_cell_c() -> int {

  burn_t state{};
  init_burn_state(state, 1.0_rt, 0.0_rt, true);

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
  amrex::Real dd = state.rho;
  amrex::Real dd1 = 0.0_rt;

  for (int n = 0; n < unit_test_rp::nsteps; n++) {

    dd1 = dd;

    amrex::Real rhotmp = 0.0_rt;

    for (int nn = 0; nn < NumSpec; ++nn) {
      rhotmp += state.xn[nn] * spmasses[nn];
    }

    // find the freefall time
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

    // do the actual integration
    {
      state.suppress_failure_output = true;
      amrex::AsyncArray<burn_t> state_device(&state, 1);
      burn_t* state_ptr = state_device.data();

      amrex::ParallelFor(1,
                         [state_ptr, dt] AMREX_GPU_DEVICE(int idx) noexcept {
                           burner(state_ptr[idx], dt);
                         });

      state_device.copyToHost(&state, 1);
      amrex::Gpu::streamSynchronize();
      state.suppress_failure_output = false;
    }

    enforce_post_burn_constraints(state);

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

AMREX_INLINE
auto burn_cell_multi_c() -> int {

  constexpr int n_zones = 128;

  amrex::Gpu::HostVector<burn_t> states_h(n_zones);
  amrex::Gpu::HostVector<amrex::Real> dd_h(n_zones);
  amrex::Gpu::HostVector<amrex::Real> dt_h(n_zones);
  amrex::Gpu::HostVector<amrex::Real> time_h(n_zones);
  std::vector<burn_t> states_before(n_zones);
  std::vector<int> active_h(n_zones, 1);
  int retry_total = 0;

  for (int iz = 0; iz < n_zones; ++iz) {
    amrex::Real temp_offset =
        0.0001_rt * unit_test_rp::temperature *
        (static_cast<amrex::Real>(iz) / static_cast<amrex::Real>(n_zones - 1) -
         0.5_rt);
    init_burn_state(states_h[iz], 1.0_rt, temp_offset, false);
    dd_h[iz] = states_h[iz].rho;
    dt_h[iz] = 0.0_rt;
    time_h[iz] = 0.0_rt;
  }

  amrex::Gpu::DeviceVector<burn_t> states_d(n_zones);
  amrex::Gpu::DeviceVector<amrex::Real> dt_d(n_zones);

  for (int step = 0; step < unit_test_rp::nsteps; ++step) {
    bool any_active = false;

    for (int iz = 0; iz < n_zones; ++iz) {
      if (active_h[iz] == 0) {
        dt_h[iz] = 0.0_rt;
        continue;
      }

      amrex::Real dd1 = dd_h[iz];
      amrex::Real rhotmp = 0.0_rt;
      for (int nn = 0; nn < NumSpec; ++nn) {
        rhotmp += states_h[iz].xn[nn] * spmasses[nn];
      }
      amrex::Real tff =
          std::sqrt(M_PI * 3.0_rt / (32.0_rt * rhotmp * grav_constant));
      amrex::Real dt = unit_test_rp::tff_reduc * tff;
      amrex::Real dd_new = dd1 + dt * (dd1 / tff);

      if (dt < 10.0_rt || dd_new > 2.0e-6_rt) {
        active_h[iz] = 0;
        dt_h[iz] = 0.0_rt;
        continue;
      }

      amrex::Real ratio = dd_new / dd1;
      for (int nn = 0; nn < NumSpec; ++nn) {
        states_h[iz].xn[nn] *= ratio;
      }
      states_h[iz].rho *= ratio;

      states_h[iz].suppress_failure_output = true;

      states_before[iz] = states_h[iz];

      dd_h[iz] = dd_new;
      dt_h[iz] = dt;
      any_active = true;
    }

    if (!any_active) {
      break;
    }

    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, states_h.begin(),
                          states_h.end(), states_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, dt_h.begin(),
                          dt_h.end(), dt_d.begin());
    amrex::Gpu::streamSynchronize();

    burn_t *states_ptr = states_d.data();
    const amrex::Real *dt_ptr = dt_d.data();

    amrex::ParallelFor(
        n_zones, [states_ptr, dt_ptr] AMREX_GPU_DEVICE(int idx) noexcept {
          amrex::Real dt = dt_ptr[idx];
          if (dt <= 0.0_rt) {
            return;
          }
          burner(states_ptr[idx], dt);
        });
    amrex::Gpu::streamSynchronize();

    amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, states_d.begin(),
                          states_d.end(), states_h.begin());
    amrex::Gpu::streamSynchronize();

    const int max_split_depth = 6;
    constexpr amrex::Real min_split_dt = 1.0_rt;
    std::function<bool(burn_t&, amrex::Real, int)> integrate_with_retries;
    integrate_with_retries = [&](burn_t& state, amrex::Real dt,
                                 int depth) -> bool {
      burn_t backup = state;
      state.time = 0.0_rt;
      state.suppress_failure_output = true;
      burner(state, dt);
      if (state.success && state.e > 0.0_rt) {
        state.suppress_failure_output = false;
        return true;
      }
      if (depth >= max_split_depth || 0.5_rt * dt <= min_split_dt) {
        state = backup;
        return false;
      }
      amrex::Real half_dt = 0.5_rt * dt;
      state = backup;
      if (!integrate_with_retries(state, half_dt, depth + 1)) {
        return false;
      }
      return integrate_with_retries(state, half_dt, depth + 1);
    };

    for (int iz = 0; iz < n_zones; ++iz) {
      if (dt_h[iz] <= 0.0_rt) {
        continue;
      }

      auto& state = states_h[iz];

      if (!state.success || state.e <= 0.0_rt) {
        burn_t retry_state = states_before[iz];
        retry_state.suppress_failure_output = false;
        ++retry_total;
        if (integrate_with_retries(retry_state, dt_h[iz], 0)) {
          state = retry_state;
        } else {
          state = states_before[iz];
          state.success = false;
          active_h[iz] = 0;
          dt_h[iz] = 0.0_rt;
          continue;
        }
      }

      state.suppress_failure_output = false;
      enforce_post_burn_constraints(state);
      time_h[iz] += dt_h[iz];
    }
  }

  int failures = 0;
  amrex::Real t_min = std::numeric_limits<amrex::Real>::max();
  amrex::Real t_max = -std::numeric_limits<amrex::Real>::max();
  bool any_success = false;

  for (int iz = 0; iz < n_zones; ++iz) {
    if (!states_h[iz].success || states_h[iz].e <= 0.0_rt) {
      ++failures;
      continue;
    }
    any_success = true;
    t_min = std::min(t_min, states_h[iz].T);
    t_max = std::max(t_max, states_h[iz].T);
  }

  std::cout << "multi-zone burn complete; zones failed = " << failures
            << std::endl;
  if (any_success) {
    std::cout << "multi-zone temperature range: [" << t_min << ", " << t_max
              << "]" << std::endl;
    std::cout << "multi-zone retries applied: " << retry_total << std::endl;
  } else {
    std::cout << "multi-zone burn had no successful zones" << std::endl;
  }

  return (failures == 0);
}
#endif
