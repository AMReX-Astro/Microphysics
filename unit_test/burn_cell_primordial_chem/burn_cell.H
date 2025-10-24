#ifndef BURN_CELL_H
#define BURN_CELL_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
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
  init_burn_state(state, 1.0_rt, 0.0_rt, false);

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

    //    std::cout << "step " << n << " done" << std::endl;

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

  constexpr int column_width = 44;
  auto print_state_row = [&](const std::string &label, amrex::Real initial_value,
                             amrex::Real final_value) {
    std::ostringstream left_stream;
    left_stream << label << " " << std::scientific << std::setprecision(6)
                << initial_value;
    std::ostringstream right_stream;
    right_stream << label << " " << std::scientific << std::setprecision(6)
                 << final_value;
    std::cout << std::left << std::setw(column_width) << left_stream.str()
              << std::left << std::setw(column_width) << right_stream.str()
              << std::endl;
  };

  std::cout << "------------------------------------" << std::endl;
  std::cout << std::left << std::setw(column_width) << "Initial State"
            << std::left << std::setw(column_width) << "Final State"
            << std::endl;
  std::cout << std::string(2 * column_width, '-') << std::endl;
  print_state_row("Temperature (K):", state_in.T, state.T);
  print_state_row("Density (g/cm^3):", state_in.rho, state.rho);
  print_state_row("Eint (erg/g):", state_in.e, state.e);
  for (int n = 0; n < NumSpec; ++n) {
    const std::string label = "n(" + short_spec_names_cxx[n] + "):";
    print_state_row(label, state_in.xn[n], state.xn[n]);
  }
  std::cout << std::string(2 * column_width, '-') << std::endl;
  std::cout << std::right;

  return state.success;
}

AMREX_INLINE
auto burn_cell_multi_c() -> int {

  constexpr int n_zones = 128;
  constexpr amrex::Real rel_tol = 1.0e-10_rt;
  constexpr amrex::Real abs_tol = 1.0e-12_rt;

  struct MultiZoneSummary {
    int failures = 0;
    int retries = 0;
    amrex::Real t_min = std::numeric_limits<amrex::Real>::max();
    amrex::Real t_max = -std::numeric_limits<amrex::Real>::max();
    bool any_success = false;
  };

  auto run_multi_zone = [&](bool use_gpu, std::vector<burn_t> &final_states,
                            std::vector<amrex::Real> &final_times,
                            MultiZoneSummary &summary) -> bool {
    amrex::Gpu::HostVector<burn_t> states_h(n_zones);
    amrex::Gpu::HostVector<amrex::Real> dd_h(n_zones);
    amrex::Gpu::HostVector<amrex::Real> dt_h(n_zones);
    amrex::Gpu::HostVector<amrex::Real> time_h(n_zones);
    std::vector<burn_t> states_before(n_zones);
    std::vector<int> active_h(n_zones, 1);

    summary = MultiZoneSummary{};
    int retry_total = 0;

    for (int iz = 0; iz < n_zones; ++iz) {
      amrex::Real temp_offset =
          0.5_rt * unit_test_rp::temperature *
          (static_cast<amrex::Real>(iz) /
               static_cast<amrex::Real>(n_zones - 1) -
           0.5_rt);
      init_burn_state(states_h[iz], 1.0_rt, temp_offset, false);
      dd_h[iz] = states_h[iz].rho;
      dt_h[iz] = 0.0_rt;
      time_h[iz] = 0.0_rt;
      active_h[iz] = 1;
    }

    amrex::Gpu::DeviceVector<burn_t> states_d;
    amrex::Gpu::DeviceVector<amrex::Real> dt_d;
    if (use_gpu) {
      states_d.resize(n_zones);
      dt_d.resize(n_zones);
    }

    const int max_split_depth = 6;
    constexpr amrex::Real min_split_dt = 1.0_rt;

    auto burn_backend = [&](burn_t &state, amrex::Real dt) {
      state.time = 0.0_rt;
      state.suppress_failure_output = true;

      if (use_gpu) {
        amrex::AsyncArray<burn_t> state_device(&state, 1);
        burn_t *state_ptr = state_device.data();

        amrex::ParallelFor(
            1, [state_ptr, dt] AMREX_GPU_DEVICE(int idx) noexcept {
              burner(state_ptr[idx], dt);
            });

        state_device.copyToHost(&state, 1);
        amrex::Gpu::streamSynchronize();
      } else {
        burner(state, dt);
      }
    };

    std::function<bool(burn_t &, amrex::Real, int)> integrate_with_retries;
    integrate_with_retries = [&](burn_t &state, amrex::Real dt,
                                 int depth) -> bool {
      burn_t backup = state;
      burn_backend(state, dt);
      if (state.success && state.e > 0.0_rt) {
        state.suppress_failure_output = false;
        return true;
      }
      if (depth >= max_split_depth || 0.5_rt * dt <= min_split_dt) {
        state = backup;
        state.suppress_failure_output = false;
        return false;
      }
      amrex::Real half_dt = 0.5_rt * dt;
      state = backup;
      if (!integrate_with_retries(state, half_dt, depth + 1)) {
        state.suppress_failure_output = false;
        return false;
      }
      bool status = integrate_with_retries(state, half_dt, depth + 1);
      state.suppress_failure_output = false;
      return status;
    };

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

      if (use_gpu) {
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
      } else {
        for (int iz = 0; iz < n_zones; ++iz) {
          amrex::Real dt = dt_h[iz];
          if (dt <= 0.0_rt) {
            continue;
          }
          burn_backend(states_h[iz], dt);
        }
      }

      for (int iz = 0; iz < n_zones; ++iz) {
        if (dt_h[iz] <= 0.0_rt) {
          continue;
        }

        auto &state = states_h[iz];

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

    summary.failures = 0;
    summary.t_min = std::numeric_limits<amrex::Real>::max();
    summary.t_max = -std::numeric_limits<amrex::Real>::max();
    summary.any_success = false;
    summary.retries = retry_total;

    for (int iz = 0; iz < n_zones; ++iz) {
      if (!states_h[iz].success || states_h[iz].e <= 0.0_rt) {
        ++summary.failures;
        continue;
      }
      summary.any_success = true;
      summary.t_min = std::min(summary.t_min, states_h[iz].T);
      summary.t_max = std::max(summary.t_max, states_h[iz].T);
    }

    final_states.assign(states_h.begin(), states_h.end());
    final_times.assign(time_h.begin(), time_h.end());

    return (summary.failures == 0);
  };

  std::vector<burn_t> cpu_states;
  std::vector<amrex::Real> cpu_times;
  MultiZoneSummary cpu_summary;
  bool cpu_success = run_multi_zone(false, cpu_states, cpu_times, cpu_summary);

  std::vector<burn_t> gpu_states;
  std::vector<amrex::Real> gpu_times;
  MultiZoneSummary gpu_summary;
  bool gpu_success = run_multi_zone(true, gpu_states, gpu_times, gpu_summary);

  auto print_summary = [&](const std::string &label,
                           const MultiZoneSummary &summary) {
    std::cout << label << " multi-zone burn complete; zones failed = "
              << summary.failures << std::endl;
    if (summary.any_success) {
      std::cout << label << " multi-zone temperature range: ["
                << summary.t_min << ", " << summary.t_max << "]"
                << std::endl;
      std::cout << label << " multi-zone retries applied: "
                << summary.retries << std::endl;
    } else {
      std::cout << label << " multi-zone burn had no successful zones"
                << std::endl;
    }
  };

  print_summary("CPU", cpu_summary);
  print_summary("GPU", gpu_summary);

  auto within_tol = [&](amrex::Real a, amrex::Real b) -> bool {
    amrex::Real diff = std::abs(a - b);
    amrex::Real scale = std::max(std::abs(a), std::abs(b));
    return diff <= (abs_tol + rel_tol * scale);
  };

  struct MismatchStats {
    int count = 0;
    amrex::Real max_abs_diff = 0.0_rt;
    int max_zone = -1;
    amrex::Real cpu_value_at_max = 0.0_rt;
    amrex::Real gpu_value_at_max = 0.0_rt;
  };

  std::map<std::string, MismatchStats> mismatch_summary;

  auto record_mismatch = [&](int zone, const std::string &field,
                             amrex::Real cpu_value, amrex::Real gpu_value) {
    auto &stats = mismatch_summary[field];
    stats.count++;
    amrex::Real diff = std::abs(cpu_value - gpu_value);
    if (stats.max_zone < 0 || diff > stats.max_abs_diff) {
      stats.max_abs_diff = diff;
      stats.max_zone = zone;
      stats.cpu_value_at_max = cpu_value;
      stats.gpu_value_at_max = gpu_value;
    }
  };

  for (int iz = 0; iz < n_zones; ++iz) {
    bool cpu_zone_success =
        cpu_states[iz].success && cpu_states[iz].e > 0.0_rt;
    bool gpu_zone_success =
        gpu_states[iz].success && gpu_states[iz].e > 0.0_rt;

    if (cpu_zone_success != gpu_zone_success) {
      record_mismatch(iz, "success_flag",
                      cpu_zone_success ? 1.0_rt : 0.0_rt,
                      gpu_zone_success ? 1.0_rt : 0.0_rt);
      continue;
    }
    if (!cpu_zone_success) {
      continue;
    }

    if (!within_tol(cpu_times[iz], gpu_times[iz])) {
      record_mismatch(iz, "time", cpu_times[iz], gpu_times[iz]);
    }

    if (!within_tol(cpu_states[iz].rho, gpu_states[iz].rho)) {
      record_mismatch(iz, "rho", cpu_states[iz].rho, gpu_states[iz].rho);
    }
    if (!within_tol(cpu_states[iz].T, gpu_states[iz].T)) {
      record_mismatch(iz, "T", cpu_states[iz].T, gpu_states[iz].T);
    }
    if (!within_tol(cpu_states[iz].e, gpu_states[iz].e)) {
      record_mismatch(iz, "e", cpu_states[iz].e, gpu_states[iz].e);
    }
    for (int n = 0; n < NumSpec; ++n) {
      if (!within_tol(cpu_states[iz].xn[n], gpu_states[iz].xn[n])) {
        std::string label =
            "xn(" + std::string(short_spec_names_cxx[n]) + ")";
        record_mismatch(iz, label, cpu_states[iz].xn[n],
                        gpu_states[iz].xn[n]);
      }
    }
  }

  if (mismatch_summary.empty()) {
    std::cout << "CPU/GPU multi-zone states agree within tolerances"
              << std::endl;
  } else {
    std::cout << "CPU/GPU mismatches detected:" << std::endl;
    for (const auto &entry : mismatch_summary) {
      const auto &field = entry.first;
      const auto &stats = entry.second;
      std::cout << "  field " << field << " mismatched in " << stats.count
                << " zone(s)";
      if (stats.max_zone >= 0) {
        std::cout << "; largest difference at zone " << stats.max_zone
                  << " (CPU = " << stats.cpu_value_at_max
                  << ", GPU = " << stats.gpu_value_at_max << ")";
      }
      std::cout << std::endl;
    }
  }

  return (cpu_success && gpu_success && mismatch_summary.empty());
}
#endif
