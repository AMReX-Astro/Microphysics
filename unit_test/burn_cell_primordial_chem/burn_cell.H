#ifndef BURN_CELL_H
#define BURN_CELL_H

#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <complex>
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
#include <AMReX_GpuLaunch.H>
#include <burner.H>
#include <eos.H>
#include <extern_parameters.H>
#include <actual_rhs.H>
#include <integrator_data.H>
#include <network.H>
#include <actual_network.H>

amrex::Real grav_constant = 6.674e-8;

namespace {

using Complex = std::complex<amrex::Real>;
using ComplexMatrix = std::vector<std::vector<Complex>>;

struct BurnZoneKernel {
  burn_t *states_ptr;
  const amrex::Real *dt_ptr;

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void operator()(int iz) const noexcept {
    amrex::Real dt = dt_ptr[iz];
    if (dt <= 0.0_rt) {
      return;
    }
    burn_t &state = states_ptr[iz];
    state.time = 0.0_rt;
    state.suppress_failure_output = true;
    burner(state, dt);
  }
};

AMREX_INLINE
auto identity_complex_matrix(int n) -> ComplexMatrix {
  ComplexMatrix I(n, std::vector<Complex>(n, Complex(0.0_rt, 0.0_rt)));
  for (int i = 0; i < n; ++i) {
    I[i][i] = Complex(1.0_rt, 0.0_rt);
  }
  return I;
}

AMREX_INLINE
auto matmul(const ComplexMatrix &A, const ComplexMatrix &B) -> ComplexMatrix {
  const int n = static_cast<int>(A.size());
  ComplexMatrix C(n, std::vector<Complex>(n, Complex(0.0_rt, 0.0_rt)));
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < n; ++k) {
      const Complex aik = A[i][k];
      if (std::abs(aik) == 0.0_rt) {
        continue;
      }
      for (int j = 0; j < n; ++j) {
        C[i][j] += aik * B[k][j];
      }
    }
  }
  return C;
}

AMREX_INLINE
void qr_decompose_mgs(const ComplexMatrix &A, ComplexMatrix &Q,
                      ComplexMatrix &R) {
  const int n = static_cast<int>(A.size());
  Q.assign(n, std::vector<Complex>(n, Complex(0.0_rt, 0.0_rt)));
  R.assign(n, std::vector<Complex>(n, Complex(0.0_rt, 0.0_rt)));

  std::vector<Complex> v(n);

  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      v[i] = A[i][j];
    }

    for (int i = 0; i < j; ++i) {
      Complex rij(0.0_rt, 0.0_rt);
      for (int k = 0; k < n; ++k) {
        rij += std::conj(Q[k][i]) * v[k];
      }
      R[i][j] = rij;
      for (int k = 0; k < n; ++k) {
        v[k] -= rij * Q[k][i];
      }
    }

    amrex::Real vnorm2 = 0.0_rt;
    for (int k = 0; k < n; ++k) {
      vnorm2 += std::norm(v[k]);
    }
    const amrex::Real vnorm = std::sqrt(vnorm2);
    R[j][j] = Complex(vnorm, 0.0_rt);

    if (vnorm > 0.0_rt) {
      for (int k = 0; k < n; ++k) {
        Q[k][j] = v[k] / vnorm;
      }
    } else {
      Q[j][j] = Complex(1.0_rt, 0.0_rt);
    }
  }
}

AMREX_INLINE
auto jacobian_eigenvalues_qr(const RArray2D &jac) -> std::vector<Complex> {
  constexpr int n = INT_NEQS;
  ComplexMatrix A(n, std::vector<Complex>(n, Complex(0.0_rt, 0.0_rt)));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      A[i][j] = Complex(static_cast<amrex::Real>(jac(i + 1, j + 1)), 0.0_rt);
    }
  }

  ComplexMatrix Q;
  ComplexMatrix R;
  const auto I = identity_complex_matrix(n);

  constexpr int max_iter = 256;
  constexpr amrex::Real offdiag_tol = 1.0e-12_rt;

  for (int iter = 0; iter < max_iter; ++iter) {
    amrex::Real offdiag_norm = 0.0_rt;
    for (int i = 1; i < n; ++i) {
      for (int j = 0; j < i; ++j) {
        offdiag_norm += std::norm(A[i][j]);
      }
    }
    offdiag_norm = std::sqrt(offdiag_norm);
    if (offdiag_norm < offdiag_tol) {
      break;
    }

    const Complex mu = A[n - 1][n - 1];
    ComplexMatrix shifted = A;
    for (int i = 0; i < n; ++i) {
      shifted[i][i] -= mu * I[i][i];
    }

    qr_decompose_mgs(shifted, Q, R);
    A = matmul(R, Q);
    for (int i = 0; i < n; ++i) {
      A[i][i] += mu * I[i][i];
    }
  }

  std::vector<Complex> eigvals(n, Complex(0.0_rt, 0.0_rt));
  for (int i = 0; i < n; ++i) {
    eigvals[i] = A[i][i];
  }
  return eigvals;
}

struct JacobianStiffness {
  amrex::Real max_abs_real = 0.0_rt;
  amrex::Real min_abs_real = std::numeric_limits<amrex::Real>::infinity();
  amrex::Real ratio = std::numeric_limits<amrex::Real>::infinity();
  bool valid = false;
};

AMREX_INLINE
auto compute_jacobian_stiffness(const burn_t &state_in) -> JacobianStiffness {
  burn_t state = state_in;
  RArray2D jac;
  actual_jac(state, jac);

  const auto eigvals = jacobian_eigenvalues_qr(jac);
  constexpr amrex::Real tiny = 1.0e-300_rt;

  JacobianStiffness result;
  for (const auto &lam : eigvals) {
    const amrex::Real ar = std::abs(lam.real());
    if (!std::isfinite(ar)) {
      continue;
    }
    result.max_abs_real = std::max(result.max_abs_real, ar);
    if (ar > tiny) {
      result.min_abs_real = std::min(result.min_abs_real, ar);
    }
  }

  if (result.max_abs_real > tiny &&
      std::isfinite(result.min_abs_real) &&
      result.min_abs_real > tiny) {
    result.ratio = result.max_abs_real / result.min_abs_real;
    result.valid = true;
  }

  return result;
}

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
auto burn_cell_multi_c(int n_zones, bool enable_gpu,
                       const std::string &compare_final_state_file = "")
    -> int {

  constexpr amrex::Real rel_tol = 1.0e-10_rt;
  constexpr amrex::Real abs_tol = 1.0e-12_rt;

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_zones > 0,
                                   "n_zones must be positive");

  struct MultiZoneSummary {
    int failures = 0;
    int retries = 0;
    bool stiffness_enabled = false;
    amrex::Real t_min = std::numeric_limits<amrex::Real>::max();
    amrex::Real t_max = -std::numeric_limits<amrex::Real>::max();
    amrex::Real stiffness_ratio_min = std::numeric_limits<amrex::Real>::max();
    amrex::Real stiffness_ratio_max = 0.0_rt;
    int stiffness_samples = 0;
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
    summary.stiffness_enabled = unit_test_rp::compute_stiffness_diagnostics;
    int retry_total = 0;

    for (int iz = 0; iz < n_zones; ++iz) {
      amrex::Real temp_offset = 0.0_rt;
      if (n_zones > 1) {
        temp_offset = 0.5_rt * unit_test_rp::temperature *
                      (static_cast<amrex::Real>(iz) /
                           static_cast<amrex::Real>(n_zones - 1) -
                       0.5_rt);
      }
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

    auto precompute_total_steps = [&]() -> int {
      auto dd_pc = dd_h;
      auto active_pc = active_h;
      int steps = 0;

      for (int step = 0; step < unit_test_rp::nsteps; ++step) {
        int any_active_int = 0;
#ifdef AMREX_USE_OMP
#pragma omp parallel for schedule(dynamic, 64) reduction(max : any_active_int)
#endif
        for (int iz = 0; iz < n_zones; ++iz) {
          if (active_pc[iz] == 0) {
            continue;
          }

          amrex::Real dd1 = dd_pc[iz];
          amrex::Real tff =
              std::sqrt(M_PI * 3.0_rt / (32.0_rt * dd1 * grav_constant));
          amrex::Real dt = unit_test_rp::tff_reduc * tff;
          amrex::Real dd_new = dd1 + dt * (dd1 / tff);

          if (dt < 10.0_rt || dd_new > 2.0e-6_rt) {
            active_pc[iz] = 0;
            continue;
          }

          dd_pc[iz] = dd_new;
          any_active_int = 1;
        }

        const bool any_active = (any_active_int != 0);
        if (!any_active) {
          break;
        }
        ++steps;
      }

      return steps;
    };

    const int total_steps = precompute_total_steps();

    for (int step = 0; step < unit_test_rp::nsteps; ++step) {
      int any_active_int = 0;
#ifdef AMREX_USE_OMP
#pragma omp parallel for schedule(dynamic, 64) reduction(max : any_active_int)
#endif
      for (int iz = 0; iz < n_zones; ++iz) {
        if (active_h[iz] == 0) {
          dt_h[iz] = 0.0_rt;
          continue;
        }

        amrex::Real dd1 = dd_h[iz];
        amrex::Real tff =
            std::sqrt(M_PI * 3.0_rt / (32.0_rt * dd1 * grav_constant));
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
        any_active_int = 1;
      }

      const bool any_active = (any_active_int != 0);
      if (!any_active) {
        break;
      }

      if (unit_test_rp::burn_progress_interval > 0 &&
          (step + 1) % unit_test_rp::burn_progress_interval == 0) {
        std::cout << "Burn progress: timestep " << (step + 1) << " / "
                  << total_steps << std::endl;
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
        burn_t *states_ptr = states_h.data();
        const amrex::Real *dt_ptr = dt_h.data();
#ifdef AMREX_USE_OMP
#pragma omp parallel for schedule(dynamic, 32)
#endif
        for (int iz = 0; iz < n_zones; ++iz) {
          amrex::Real dt = dt_ptr[iz];
          if (dt <= 0.0_rt) {
            continue;
          }
          burn_t &state = states_ptr[iz];
          state.time = 0.0_rt;
          state.suppress_failure_output = true;
          burner(state, dt);
        }
      }

      int retry_total_step = 0;
      amrex::Real stiffness_ratio_min_step =
          std::numeric_limits<amrex::Real>::max();
      amrex::Real stiffness_ratio_max_step = 0.0_rt;
      int stiffness_samples_step = 0;

#ifdef AMREX_USE_OMP
#pragma omp parallel for schedule(dynamic, 32) \
    reduction(+ : retry_total_step, stiffness_samples_step) \
    reduction(min : stiffness_ratio_min_step) reduction(max : stiffness_ratio_max_step)
#endif
      for (int iz = 0; iz < n_zones; ++iz) {
        if (dt_h[iz] <= 0.0_rt) {
          continue;
        }

        auto &state = states_h[iz];

        if (!state.success || state.e <= 0.0_rt) {
          burn_t retry_state = states_before[iz];
          retry_state.suppress_failure_output = false;
          ++retry_total_step;
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
        if (summary.stiffness_enabled) {
          const auto stiff = compute_jacobian_stiffness(state);
          if (stiff.valid) {
            stiffness_ratio_min_step =
                std::min(stiffness_ratio_min_step, stiff.ratio);
            stiffness_ratio_max_step =
                std::max(stiffness_ratio_max_step, stiff.ratio);
            ++stiffness_samples_step;
          }
        }
        time_h[iz] += dt_h[iz];
      }

      retry_total += retry_total_step;
      if (summary.stiffness_enabled && stiffness_samples_step > 0) {
        summary.stiffness_ratio_min =
            std::min(summary.stiffness_ratio_min, stiffness_ratio_min_step);
        summary.stiffness_ratio_max =
            std::max(summary.stiffness_ratio_max, stiffness_ratio_max_step);
        summary.stiffness_samples += stiffness_samples_step;
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

  auto print_summary = [&](const std::string &label,
                           const MultiZoneSummary &summary) {
    const std::string prefix = label + " multi-zone";
    std::cout << prefix << " summary:" << std::endl;
    std::cout << "  burn complete             : zones failed = "
              << summary.failures << std::endl;
    if (summary.any_success) {
      std::cout << "  temperature range         : ["
                << summary.t_min << ", " << summary.t_max << "]" << std::endl;
      std::cout << "  retries applied           : "
                << summary.retries << std::endl;
      if (summary.stiffness_enabled) {
        if (summary.stiffness_samples > 0) {
          std::cout << "  stiffness ratio range     : ["
                    << std::scientific << summary.stiffness_ratio_min << ", "
                    << summary.stiffness_ratio_max << "]"
                    << std::defaultfloat << " (samples="
                    << summary.stiffness_samples << ")" << std::endl;
        } else {
          std::cout << "  stiffness ratio range     : unavailable" << std::endl;
        }
      } else {
        std::cout << "  stiffness diagnostics     : disabled" << std::endl;
      }
    } else {
      std::cout << "  status                    : no successful zones"
                << std::endl;
    }
  };

  struct BackendRunResult {
    std::vector<burn_t> states;
    std::vector<amrex::Real> times;
    MultiZoneSummary summary;
    bool success = false;
    std::chrono::duration<double> wall{0.0};
  };

  struct TolLevelResult {
    amrex::Real scale = 1.0_rt;
    BackendRunResult cpu;
    BackendRunResult gpu;
  };

  const amrex::Real base_rtol_spec = integrator_rp::rtol_spec;
  const amrex::Real base_atol_spec = integrator_rp::atol_spec;
  const amrex::Real base_rtol_enuc = integrator_rp::rtol_enuc;
  const amrex::Real base_atol_enuc = integrator_rp::atol_enuc;
  const amrex::Real base_retry_rtol_spec = integrator_rp::retry_rtol_spec;
  const amrex::Real base_retry_atol_spec = integrator_rp::retry_atol_spec;
  const amrex::Real base_retry_rtol_enuc = integrator_rp::retry_rtol_enuc;
  const amrex::Real base_retry_atol_enuc = integrator_rp::retry_atol_enuc;

  const bool compare_against_saved = !compare_final_state_file.empty();
  const bool do_richardson =
      !compare_against_saved && unit_test_rp::richardson_enable;
  const int n_tol_levels = do_richardson ? 3 : 1;
  const std::array<amrex::Real, 3> richardson_scales = {
      0.5_rt, 1.0_rt, 2.0_rt};
  constexpr int richardson_loose_level = 0;
  constexpr int richardson_base_level = 1;
  constexpr int richardson_tight_level = 2;
  const amrex::Real richardson_ratio =
      richardson_scales[richardson_tight_level] /
      richardson_scales[richardson_base_level];
  const bool compare_backend_is_gpu = compare_against_saved && enable_gpu;
  const bool run_gpu_backend =
      compare_against_saved ? compare_backend_is_gpu : enable_gpu;

  auto apply_tolerance_scale = [&](amrex::Real scale) {
    integrator_rp::rtol_spec = base_rtol_spec / scale;
    integrator_rp::atol_spec = base_atol_spec / scale;
    integrator_rp::rtol_enuc = base_rtol_enuc / scale;
    integrator_rp::atol_enuc = base_atol_enuc / scale;

    if (base_retry_rtol_spec > 0.0_rt) {
      integrator_rp::retry_rtol_spec = base_retry_rtol_spec / scale;
    }
    if (base_retry_atol_spec > 0.0_rt) {
      integrator_rp::retry_atol_spec = base_retry_atol_spec / scale;
    }
    if (base_retry_rtol_enuc > 0.0_rt) {
      integrator_rp::retry_rtol_enuc = base_retry_rtol_enuc / scale;
    }
    if (base_retry_atol_enuc > 0.0_rt) {
      integrator_rp::retry_atol_enuc = base_retry_atol_enuc / scale;
    }
  };

  auto restore_base_tolerances = [&]() {
    integrator_rp::rtol_spec = base_rtol_spec;
    integrator_rp::atol_spec = base_atol_spec;
    integrator_rp::rtol_enuc = base_rtol_enuc;
    integrator_rp::atol_enuc = base_atol_enuc;
    integrator_rp::retry_rtol_spec = base_retry_rtol_spec;
    integrator_rp::retry_atol_spec = base_retry_atol_spec;
    integrator_rp::retry_rtol_enuc = base_retry_rtol_enuc;
    integrator_rp::retry_atol_enuc = base_retry_atol_enuc;
  };

  auto format_tol = [](amrex::Real tol) {
    std::ostringstream out;
    out << std::scientific << std::setprecision(10) << tol;
    return out.str();
  };

  std::vector<TolLevelResult> level_results(n_tol_levels);
  bool cpu_success = true;
  bool gpu_success = true;

  for (int ilev = 0; ilev < n_tol_levels; ++ilev) {
    amrex::Real scale = do_richardson ? richardson_scales[ilev] : 1.0_rt;
    level_results[ilev].scale = scale;
    apply_tolerance_scale(scale);

    if (do_richardson) {
      std::cout << "Richardson tolerances (level " << (ilev + 1)
                << ", scale=" << scale << "):" << std::endl;
      std::cout << "  rtol_spec / atol_spec      : "
                << format_tol(integrator_rp::rtol_spec) << " / "
                << format_tol(integrator_rp::atol_spec) << std::endl;
      std::cout << "  rtol_enuc / atol_enuc      : "
                << format_tol(integrator_rp::rtol_enuc) << " / "
                << format_tol(integrator_rp::atol_enuc) << std::endl;
      if (integrator_rp::retry_rtol_spec > 0.0_rt &&
          integrator_rp::retry_atol_spec > 0.0_rt) {
        std::cout << "  retry_rtol_spec / retry_atol_spec : "
                  << format_tol(integrator_rp::retry_rtol_spec) << " / "
                  << format_tol(integrator_rp::retry_atol_spec) << std::endl;
      }
      if (integrator_rp::retry_rtol_enuc > 0.0_rt &&
          integrator_rp::retry_atol_enuc > 0.0_rt) {
        std::cout << "  retry_rtol_enuc / retry_atol_enuc : "
                  << format_tol(integrator_rp::retry_rtol_enuc) << " / "
                  << format_tol(integrator_rp::retry_atol_enuc) << std::endl;
      }
    }

    std::ostringstream level_label;
    if (do_richardson) {
      level_label << " (tol level " << (ilev + 1) << ", scale=" << scale
                  << ")";
    }
    if (!compare_against_saved || !compare_backend_is_gpu) {
      auto cpu_start = std::chrono::steady_clock::now();
      level_results[ilev].cpu.success =
          run_multi_zone(false, level_results[ilev].cpu.states,
                         level_results[ilev].cpu.times,
                         level_results[ilev].cpu.summary);
      auto cpu_end = std::chrono::steady_clock::now();
      level_results[ilev].cpu.wall = cpu_end - cpu_start;
      cpu_success = cpu_success && level_results[ilev].cpu.success;

      const std::string cpu_label = "CPU" + level_label.str();
      print_summary(cpu_label, level_results[ilev].cpu.summary);
      std::cout << "  burn walltime (s)         : "
                << level_results[ilev].cpu.wall.count() << std::endl;
    } else if (ilev == 0) {
      std::cout << "CPU multi-zone burn skipped (saved-state compare mode)"
                << std::endl;
    }

    if (run_gpu_backend) {
      auto gpu_start = std::chrono::steady_clock::now();
      level_results[ilev].gpu.success =
          run_multi_zone(true, level_results[ilev].gpu.states,
                         level_results[ilev].gpu.times,
                         level_results[ilev].gpu.summary);
      auto gpu_end = std::chrono::steady_clock::now();
      level_results[ilev].gpu.wall = gpu_end - gpu_start;
      gpu_success = gpu_success && level_results[ilev].gpu.success;

      const std::string gpu_label = "GPU" + level_label.str();
      print_summary(gpu_label, level_results[ilev].gpu.summary);
      std::cout << "  burn walltime (s)         : "
                << level_results[ilev].gpu.wall.count() << std::endl;
    } else if (ilev == 0) {
      std::cout << "GPU multi-zone burn skipped (no GPU backend)"
                << std::endl;
    }
  }

  restore_base_tolerances();

  const auto &cpu_states = level_results.front().cpu.states;
  const auto &cpu_times = level_results.front().cpu.times;
  const auto &gpu_states = level_results.front().gpu.states;
  const auto &gpu_times = level_results.front().gpu.times;
  const auto &compare_states = compare_backend_is_gpu ? gpu_states : cpu_states;
  const auto &compare_times = compare_backend_is_gpu ? gpu_times : cpu_times;
  const bool compare_backend_success =
      compare_backend_is_gpu ? level_results.front().gpu.success
                             : level_results.front().cpu.success;

  struct RichardsonEstimate {
    bool valid = false;
    amrex::Real p = 0.0_rt;
    amrex::Real err = 0.0_rt;
    amrex::Real delta01 = 0.0_rt;
  };

  auto estimate_richardson = [&](amrex::Real y_loose, amrex::Real y_base,
                                 amrex::Real y_tight) -> RichardsonEstimate {
    RichardsonEstimate est;
    const amrex::Real d01 = y_loose - y_base;
    const amrex::Real d12 = y_base - y_tight;
    est.delta01 = d01;

    constexpr amrex::Real tiny = 1.0e-300_rt;
    const amrex::Real abs_d01 = std::abs(d01);
    const amrex::Real abs_d12 = std::abs(d12);

    if (abs_d01 <= tiny && abs_d12 <= tiny) {
      est.valid = true;
      return est;
    }

    if (abs_d01 <= tiny || abs_d12 <= tiny) {
      est.err = d01;
      return est;
    }

    amrex::Real ratio = abs_d01 / abs_d12;
    if (!(ratio > 0.0_rt) || !std::isfinite(ratio)) {
      est.err = d01;
      return est;
    }

    est.p = std::log(ratio) / std::log(richardson_ratio);
    if (!std::isfinite(est.p)) {
      est.err = d01;
      return est;
    }

    amrex::Real denom = std::pow(richardson_ratio, est.p) - 1.0_rt;
    if (std::abs(denom) <= tiny || !std::isfinite(denom)) {
      est.err = d01;
      return est;
    }

    est.err = d01 / denom;
    est.valid = true;
    return est;
  };

  auto within_tol = [&](amrex::Real a, amrex::Real b) -> bool {
    amrex::Real diff = std::abs(a - b);
    amrex::Real scale = std::max(std::abs(a), std::abs(b));
    return diff <= (abs_tol + rel_tol * scale);
  };

  if (compare_against_saved) {
    struct SavedZoneData {
      bool present = false;
      int zone = -1;
      amrex::Real time = 0.0_rt;
      amrex::Real rho = 0.0_rt;
      amrex::Real T = 0.0_rt;
      amrex::Real e = 0.0_rt;
      int success = 0;
      std::array<amrex::Real, NumSpec> xn{};
      int trunc_valid = 0;
      amrex::Real trunc_err_T = 0.0_rt;
      amrex::Real trunc_err_e = 0.0_rt;
      std::array<amrex::Real, NumSpec> trunc_err_xn{};
    };

    std::ifstream ifs(compare_final_state_file);
    if (!ifs) {
      std::cout << "ERROR: unable to open saved final-state file "
                << compare_final_state_file << std::endl;
      return false;
    }

    std::vector<SavedZoneData> saved(n_zones);
    std::string line;
    int parsed_rows = 0;
    while (std::getline(ifs, line)) {
      if (line.empty() || line[0] == '#') {
        continue;
      }

      std::istringstream iss(line);
      SavedZoneData row;
      if (!(iss >> row.zone >> row.time >> row.rho >> row.T >> row.e >>
            row.success)) {
        std::cout << "ERROR: malformed row in saved-state file: " << line
                  << std::endl;
        return false;
      }
      for (int n = 0; n < NumSpec; ++n) {
        if (!(iss >> row.xn[n])) {
          std::cout << "ERROR: missing species values in saved-state file row: "
                    << line << std::endl;
          return false;
        }
      }
      if (!(iss >> row.trunc_valid >> row.trunc_err_T >> row.trunc_err_e)) {
        std::cout << "ERROR: missing truncation fields in saved-state file row: "
                  << line << std::endl;
        return false;
      }
      for (int n = 0; n < NumSpec; ++n) {
        if (!(iss >> row.trunc_err_xn[n])) {
          std::cout << "ERROR: missing species truncation fields in saved-state "
                       "file row: "
                    << line << std::endl;
          return false;
        }
      }

      if (row.zone < 0 || row.zone >= n_zones) {
        std::cout << "ERROR: zone index " << row.zone
                  << " is out of bounds for n_zones = " << n_zones
                  << std::endl;
        return false;
      }
      row.present = true;
      saved[row.zone] = row;
      ++parsed_rows;
    }

    std::cout << "Loaded " << parsed_rows << " zone rows from "
              << compare_final_state_file << std::endl;

    int missing_rows = 0;
    int zones_failed = 0;
    int zones_passed = 0;
    amrex::Real max_saved_err_units = 0.0_rt;
    int max_saved_err_zone = -1;
    std::string max_saved_err_field = "none";

    for (int iz = 0; iz < n_zones; ++iz) {
      if (!saved[iz].present) {
        ++missing_rows;
        ++zones_failed;
        std::cout << "zone " << iz << " : FAIL (missing from saved file)"
                  << std::endl;
        continue;
      }

      const auto &ref = saved[iz];
      std::vector<std::string> mismatches;

      const bool cur_success =
          compare_states[iz].success && compare_states[iz].e > 0.0_rt;
      const bool ref_success = (ref.success != 0);
      if (cur_success != ref_success) {
        mismatches.push_back("success");
      }

      if (!within_tol(compare_times[iz], ref.time)) {
        mismatches.push_back("time");
      }
      if (!within_tol(compare_states[iz].rho, ref.rho)) {
        mismatches.push_back("rho");
      }

      auto within_saved_err = [&](amrex::Real cur, amrex::Real reference,
                                  amrex::Real saved_err) {
        const amrex::Real threshold = 2.0_rt * std::abs(saved_err);
        return std::abs(cur - reference) <= threshold;
      };
      auto update_saved_err_units = [&](const std::string &field,
                                        amrex::Real cur,
                                        amrex::Real reference,
                                        amrex::Real saved_err) {
        const amrex::Real diff = std::abs(cur - reference);
        amrex::Real units = 0.0_rt;
        if (std::abs(saved_err) > 0.0_rt) {
          units = diff / std::abs(saved_err);
        } else if (diff > 0.0_rt) {
          units = std::numeric_limits<amrex::Real>::infinity();
        }
        if (units > max_saved_err_units) {
          max_saved_err_units = units;
          max_saved_err_zone = iz;
          max_saved_err_field = field;
        }
      };

      if (ref.trunc_valid != 0) {
        update_saved_err_units("T", compare_states[iz].T, ref.T,
                               ref.trunc_err_T);
        if (!within_saved_err(compare_states[iz].T, ref.T, ref.trunc_err_T)) {
          mismatches.push_back("T");
        }
        update_saved_err_units("e", compare_states[iz].e, ref.e,
                               ref.trunc_err_e);
        if (!within_saved_err(compare_states[iz].e, ref.e, ref.trunc_err_e)) {
          mismatches.push_back("e");
        }
        for (int n = 0; n < NumSpec; ++n) {
          const std::string species_field =
              "xn(" + std::string(short_spec_names_cxx[n]) + ")";
          update_saved_err_units(species_field, compare_states[iz].xn[n],
                                 ref.xn[n], ref.trunc_err_xn[n]);
          if (!within_saved_err(compare_states[iz].xn[n], ref.xn[n],
                                ref.trunc_err_xn[n])) {
            mismatches.push_back(species_field);
          }
        }
      } else {
        if (!within_tol(compare_states[iz].T, ref.T)) {
          mismatches.push_back("T");
        }
        if (!within_tol(compare_states[iz].e, ref.e)) {
          mismatches.push_back("e");
        }
        for (int n = 0; n < NumSpec; ++n) {
          if (!within_tol(compare_states[iz].xn[n], ref.xn[n])) {
            mismatches.push_back("xn(" + std::string(short_spec_names_cxx[n]) +
                                 ")");
          }
        }
      }

      if (mismatches.empty()) {
        ++zones_passed;
        std::cout << "zone " << iz << " : PASS" << std::endl;
      } else {
        ++zones_failed;
        std::cout << "zone " << iz << " : FAIL [";
        for (int i = 0; i < static_cast<int>(mismatches.size()); ++i) {
          if (i > 0) {
            std::cout << ", ";
          }
          std::cout << mismatches[i];
        }
        std::cout << "]" << std::endl;
      }
    }

    std::cout << "Saved-state comparison summary: passed=" << zones_passed
              << ", failed=" << zones_failed
              << ", missing_rows=" << missing_rows << std::endl;
    if (max_saved_err_zone >= 0) {
      std::cout << "Largest saved-state difference in truncation-error units: "
                << std::fixed << std::setprecision(3) << max_saved_err_units
                << std::defaultfloat << " (zone " << max_saved_err_zone
                << ", field " << max_saved_err_field << ")" << std::endl;
    }
    return (compare_backend_success && zones_failed == 0);
  }

  auto report_richardson = [&](const std::string &label,
                               const std::vector<burn_t> &l0,
                               const std::vector<burn_t> &l1,
                               const std::vector<burn_t> &l2) {
    if (!do_richardson) {
      return;
    }
    if (l0.size() != l1.size() || l0.size() != l2.size()) {
      std::cout << label
                << " Richardson estimate skipped (inconsistent level sizes)"
                << std::endl;
      return;
    }

    int valid_zones = 0;
    int valid_temp_estimates = 0;
    amrex::Real max_abs_temp_err = 0.0_rt;
    amrex::Real max_abs_energy_err = 0.0_rt;
    amrex::Real max_abs_species_err = 0.0_rt;
    amrex::Real max_rel_temp_err = 0.0_rt;
    amrex::Real max_rel_energy_err = 0.0_rt;
    amrex::Real max_rel_species_err = 0.0_rt;
    int species_rel_count = 0;
    std::vector<int> species_rel_hit_count(NumSpec, 0);

    const amrex::Real trace_rel_cutoff =
        amrex::max(unit_test_rp::richardson_trace_species_cutoff, 0.0_rt);

    auto rel_error = [](amrex::Real err, amrex::Real ref) {
      constexpr amrex::Real tiny = 1.0e-300_rt;
      return std::abs(err) / amrex::max(std::abs(ref), tiny);
    };

    std::cout << label << " Richardson truncation estimate:" << std::endl;
    std::cout << "  levels with scales         : ["
              << richardson_scales[richardson_loose_level] << ", "
              << richardson_scales[richardson_base_level] << ", "
              << richardson_scales[richardson_tight_level] << "]"
              << std::endl;

    for (int iz = 0; iz < n_zones; ++iz) {
      const bool ok0 = l0[iz].success && l0[iz].e > 0.0_rt;
      const bool ok1 = l1[iz].success && l1[iz].e > 0.0_rt;
      const bool ok2 = l2[iz].success && l2[iz].e > 0.0_rt;
      if (!(ok0 && ok1 && ok2)) {
        continue;
      }

      ++valid_zones;

      const RichardsonEstimate t_est =
          estimate_richardson(l0[iz].T, l1[iz].T, l2[iz].T);
      const RichardsonEstimate e_est =
          estimate_richardson(l0[iz].e, l1[iz].e, l2[iz].e);

      if (t_est.valid) {
        ++valid_temp_estimates;
      }
      max_abs_temp_err = std::max(max_abs_temp_err, std::abs(t_est.err));
      max_abs_energy_err = std::max(max_abs_energy_err, std::abs(e_est.err));
      max_rel_temp_err =
          std::max(max_rel_temp_err, rel_error(t_est.err, l1[iz].T));
      max_rel_energy_err =
          std::max(max_rel_energy_err, rel_error(e_est.err, l1[iz].e));

      amrex::Real zone_species_linf = 0.0_rt;
      amrex::Real zone_rel_species_linf = 0.0_rt;
      std::ostringstream zone_species_list;
      bool zone_first_species = true;
      amrex::Real zone_max_xn = 0.0_rt;
      for (int n = 0; n < NumSpec; ++n) {
        zone_max_xn = std::max(zone_max_xn, std::abs(l1[iz].xn[n]));
      }
      const amrex::Real zone_trace_cutoff = trace_rel_cutoff * zone_max_xn;
      for (int n = 0; n < NumSpec; ++n) {
        const RichardsonEstimate xn_est =
            estimate_richardson(l0[iz].xn[n], l1[iz].xn[n], l2[iz].xn[n]);
        zone_species_linf = std::max(zone_species_linf, std::abs(xn_est.err));
        if (std::abs(l1[iz].xn[n]) >= zone_trace_cutoff) {
          zone_rel_species_linf = std::max(
              zone_rel_species_linf, rel_error(xn_est.err, l1[iz].xn[n]));
          ++species_rel_count;
          ++species_rel_hit_count[n];
          if (!zone_first_species) {
            zone_species_list << ", ";
          }
          zone_species_list << short_spec_names_cxx[n];
          zone_first_species = false;
        }
      }
      max_abs_species_err = std::max(max_abs_species_err, zone_species_linf);
      max_rel_species_err =
          std::max(max_rel_species_err, zone_rel_species_linf);

      if (unit_test_rp::richardson_print_per_zone) {
        std::cout << label << " zone " << iz << ": |E_T|=" << std::abs(t_est.err)
                  << ", rel(E_T)=" << rel_error(t_est.err, l2[iz].T)
                  << ", p_T=" << t_est.p << ", |E_e|=" << std::abs(e_est.err)
                  << ", rel(E_e)=" << rel_error(e_est.err, l2[iz].e)
                  << ", p_e=" << e_est.p << ", |E_xn|_inf=" << zone_species_linf
                  << ", rel(E_xn)_inf(non-trace)=" << zone_rel_species_linf
                  << ", non-trace species: ["
                  << (zone_first_species ? "none" : zone_species_list.str())
                  << "]"
                  << std::endl;
      }
    }

    std::cout << label << " Richardson summary:" << std::endl;
    std::cout << "  valid zones                : " << valid_zones << "/"
              << n_zones << std::endl;
    std::cout << "  valid temperature estimates: " << valid_temp_estimates
              << std::endl;
    std::cout << "  max |E_T|                  : " << max_abs_temp_err
              << std::endl;
    std::cout << "  max |E_e|                  : " << max_abs_energy_err
              << std::endl;
    std::cout << "  max |E_xn|_inf             : " << max_abs_species_err
              << std::endl;
    std::cout << "  max rel(E_T)               : " << max_rel_temp_err
              << std::endl;
    std::cout << "  max rel(E_e)               : " << max_rel_energy_err
              << std::endl;
    std::cout << "  max rel(E_xn)_inf(non-trace): " << max_rel_species_err
              << std::endl;
    std::cout << "  trace rel cutoff           : " << trace_rel_cutoff
              << std::endl;
    std::cout << "  species count (relative)   : " << species_rel_count
              << std::endl;

    std::ostringstream considered_species;
    bool first_species = true;
    for (int n = 0; n < NumSpec; ++n) {
      if (species_rel_hit_count[n] > 0) {
        if (!first_species) {
          considered_species << ", ";
        }
        considered_species << short_spec_names_cxx[n];
        first_species = false;
      }
    }
    std::cout << label << " non-trace species used in relative metric:"
              << std::endl;
    std::cout << "  species                    : ["
              << (first_species ? "none" : considered_species.str()) << "]"
              << std::endl;
  };

  if (do_richardson && n_tol_levels == 3) {
    report_richardson("CPU", level_results[richardson_loose_level].cpu.states,
                      level_results[richardson_base_level].cpu.states,
                      level_results[richardson_tight_level].cpu.states);
  }

  struct ZoneTruncationEstimate {
    bool valid = false;
    amrex::Real temp_err = 0.0_rt;
    amrex::Real energy_err = 0.0_rt;
    std::array<amrex::Real, NumSpec> species_err{};
  };

  std::vector<ZoneTruncationEstimate> cpu_zone_truncation;
  if (do_richardson && n_tol_levels == 3) {
    cpu_zone_truncation.resize(n_zones);
    const auto &cpu_l0 = level_results[richardson_loose_level].cpu.states;
    const auto &cpu_l1 = level_results[richardson_base_level].cpu.states;
    const auto &cpu_l2 = level_results[richardson_tight_level].cpu.states;

    for (int iz = 0; iz < n_zones; ++iz) {
      auto &est = cpu_zone_truncation[iz];
      const bool ok0 = cpu_l0[iz].success && cpu_l0[iz].e > 0.0_rt;
      const bool ok1 = cpu_l1[iz].success && cpu_l1[iz].e > 0.0_rt;
      const bool ok2 = cpu_l2[iz].success && cpu_l2[iz].e > 0.0_rt;
      if (!(ok0 && ok1 && ok2)) {
        continue;
      }

      est.valid = true;
      est.temp_err =
          estimate_richardson(cpu_l0[iz].T, cpu_l1[iz].T, cpu_l2[iz].T).err;
      est.energy_err =
          estimate_richardson(cpu_l0[iz].e, cpu_l1[iz].e, cpu_l2[iz].e).err;
      for (int n = 0; n < NumSpec; ++n) {
        est.species_err[n] =
            estimate_richardson(cpu_l0[iz].xn[n], cpu_l1[iz].xn[n],
                                cpu_l2[iz].xn[n])
                .err;
      }
    }
  }

  {
    const std::string output_file = "burn_cell_final_state.txt";
    std::ofstream ofs(output_file);
    if (!ofs) {
      std::cout << "WARNING: unable to open " << output_file
                << " for writing final state output" << std::endl;
    } else {
      ofs << std::scientific << std::setprecision(16);
      ofs << "# burn_cell_primordial_chem final state (CPU, tolerance level 1)\n";
      ofs << "# columns: zone time rho T e success";
      for (int n = 0; n < NumSpec; ++n) {
        ofs << " xn_" << short_spec_names_cxx[n];
      }
      ofs << " truncation_valid truncation_err_T truncation_err_e";
      for (int n = 0; n < NumSpec; ++n) {
        ofs << " truncation_err_xn_" << short_spec_names_cxx[n];
      }
      ofs << "\n";

      for (int iz = 0; iz < n_zones; ++iz) {
        ofs << iz << " " << cpu_times[iz] << " " << cpu_states[iz].rho << " "
            << cpu_states[iz].T << " " << cpu_states[iz].e << " "
            << (cpu_states[iz].success ? 1 : 0);
        for (int n = 0; n < NumSpec; ++n) {
          ofs << " " << cpu_states[iz].xn[n];
        }

        if (!cpu_zone_truncation.empty() && cpu_zone_truncation[iz].valid) {
          ofs << " 1 " << cpu_zone_truncation[iz].temp_err << " "
              << cpu_zone_truncation[iz].energy_err;
          for (int n = 0; n < NumSpec; ++n) {
            ofs << " " << cpu_zone_truncation[iz].species_err[n];
          }
        } else {
          ofs << " 0 0.0 0.0";
          for (int n = 0; n < NumSpec; ++n) {
            ofs << " 0.0";
          }
        }
        ofs << "\n";
      }
      std::cout << "Wrote final state and truncation estimates to "
                << output_file << std::endl;
    }
  }

  const bool use_richardson_gpu_failure =
      enable_gpu && do_richardson && (n_tol_levels == 3);

  struct ZoneRichardsonThresholds {
    bool valid = false;
    amrex::Real temp_err = 0.0_rt;
    amrex::Real energy_err = 0.0_rt;
    std::array<amrex::Real, NumSpec> species_err{};
    std::array<int, NumSpec> species_non_trace{};
  };

  std::vector<ZoneRichardsonThresholds> cpu_zone_thresholds;
  if (use_richardson_gpu_failure) {
    cpu_zone_thresholds.resize(n_zones);
    const auto &cpu_l0 = level_results[richardson_loose_level].cpu.states;
    const auto &cpu_l1 = level_results[richardson_base_level].cpu.states;
    const auto &cpu_l2 = level_results[richardson_tight_level].cpu.states;
    const amrex::Real trace_rel_cutoff =
        amrex::max(unit_test_rp::richardson_trace_species_cutoff, 0.0_rt);

    for (int iz = 0; iz < n_zones; ++iz) {
      auto &thr = cpu_zone_thresholds[iz];
      const bool ok0 = cpu_l0[iz].success && cpu_l0[iz].e > 0.0_rt;
      const bool ok1 = cpu_l1[iz].success && cpu_l1[iz].e > 0.0_rt;
      const bool ok2 = cpu_l2[iz].success && cpu_l2[iz].e > 0.0_rt;
      if (!(ok0 && ok1 && ok2)) {
        continue;
      }

      thr.valid = true;
      thr.temp_err =
          std::abs(estimate_richardson(cpu_l0[iz].T, cpu_l1[iz].T, cpu_l2[iz].T)
                       .err);
      thr.energy_err =
          std::abs(estimate_richardson(cpu_l0[iz].e, cpu_l1[iz].e, cpu_l2[iz].e)
                       .err);

      amrex::Real zone_max_xn = 0.0_rt;
      for (int n = 0; n < NumSpec; ++n) {
        zone_max_xn = std::max(zone_max_xn, std::abs(cpu_l1[iz].xn[n]));
      }
      const amrex::Real zone_trace_cutoff = trace_rel_cutoff * zone_max_xn;

      for (int n = 0; n < NumSpec; ++n) {
        thr.species_err[n] = std::abs(estimate_richardson(
                                          cpu_l0[iz].xn[n], cpu_l1[iz].xn[n],
                                          cpu_l2[iz].xn[n])
                                          .err);
        thr.species_non_trace[n] =
            (std::abs(cpu_l1[iz].xn[n]) >= zone_trace_cutoff) ? 1 : 0;
      }
    }
  }

  struct MismatchStats {
    int count = 0;
    amrex::Real max_abs_diff = 0.0_rt;
    int max_zone = -1;
    amrex::Real cpu_value_at_max = 0.0_rt;
    amrex::Real gpu_value_at_max = 0.0_rt;
  };

  struct GpuFailureSummary {
    int zones_failed = 0;
    int temp_failures = 0;
    int energy_failures = 0;
    int species_failures = 0;
    int missing_richardson = 0;
  };

  std::map<std::string, MismatchStats> mismatch_summary;
  bool gpu_failure_detected = false;
  GpuFailureSummary gpu_failure_summary;

  struct NormAccumulator {
    amrex::Real diff_l1 = 0.0_rt;
    amrex::Real diff_l2_sq = 0.0_rt;
    amrex::Real diff_linf = 0.0_rt;
    amrex::Real ref_l1 = 0.0_rt;
    amrex::Real ref_l2_sq = 0.0_rt;
    amrex::Real ref_linf = 0.0_rt;
    void record(amrex::Real diff, amrex::Real ref) {
      amrex::Real abs_diff = std::abs(diff);
      amrex::Real abs_ref = std::abs(ref);
      diff_l1 += abs_diff;
      diff_l2_sq += diff * diff;
      diff_linf = std::max(diff_linf, abs_diff);
      ref_l1 += abs_ref;
      ref_l2_sq += ref * ref;
      ref_linf = std::max(ref_linf, abs_ref);
    }
  };

  std::map<std::string, NormAccumulator> field_norms;

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

  auto record_difference = [&](const std::string &field, amrex::Real cpu_value,
                               amrex::Real gpu_value, bool track_norm) {
    if (!track_norm) {
      return;
    }
    field_norms[field].record(cpu_value - gpu_value, cpu_value);
  };

  if (enable_gpu) {
    for (int iz = 0; iz < n_zones; ++iz) {
      bool cpu_zone_success =
          cpu_states[iz].success && cpu_states[iz].e > 0.0_rt;
      bool gpu_zone_success =
          gpu_states[iz].success && gpu_states[iz].e > 0.0_rt;

      if (cpu_zone_success != gpu_zone_success) {
        record_mismatch(iz, "success_flag",
                        cpu_zone_success ? 1.0_rt : 0.0_rt,
                        gpu_zone_success ? 1.0_rt : 0.0_rt);
        record_difference("success_flag",
                          cpu_zone_success ? 1.0_rt : 0.0_rt,
                          gpu_zone_success ? 1.0_rt : 0.0_rt, false);
        if (use_richardson_gpu_failure) {
          gpu_failure_detected = true;
          ++gpu_failure_summary.zones_failed;
        }
        continue;
      }
      if (!cpu_zone_success) {
        continue;
      }

      if (!within_tol(cpu_times[iz], gpu_times[iz])) {
        record_mismatch(iz, "time", cpu_times[iz], gpu_times[iz]);
      }
      record_difference("time", cpu_times[iz], gpu_times[iz], false);

      if (!within_tol(cpu_states[iz].rho, gpu_states[iz].rho)) {
        record_mismatch(iz, "rho", cpu_states[iz].rho, gpu_states[iz].rho);
      }
      record_difference("rho", cpu_states[iz].rho, gpu_states[iz].rho, false);

      bool zone_failed_by_richardson = false;
      if (use_richardson_gpu_failure) {
        const auto &thr = cpu_zone_thresholds[iz];
        if (!thr.valid) {
          record_mismatch(iz, "richardson_unavailable", 1.0_rt, 0.0_rt);
          zone_failed_by_richardson = true;
          ++gpu_failure_summary.missing_richardson;
        } else {
          if (std::abs(cpu_states[iz].T - gpu_states[iz].T) >
              2.0_rt * thr.temp_err) {
            record_mismatch(iz, "T_richardson_fail", cpu_states[iz].T,
                            gpu_states[iz].T);
            zone_failed_by_richardson = true;
            ++gpu_failure_summary.temp_failures;
          }
          if (std::abs(cpu_states[iz].e - gpu_states[iz].e) >
              2.0_rt * thr.energy_err) {
            record_mismatch(iz, "e_richardson_fail", cpu_states[iz].e,
                            gpu_states[iz].e);
            zone_failed_by_richardson = true;
            ++gpu_failure_summary.energy_failures;
          }
          for (int n = 0; n < NumSpec; ++n) {
            if (thr.species_non_trace[n] == 0) {
              continue;
            }
            if (std::abs(cpu_states[iz].xn[n] - gpu_states[iz].xn[n]) >
                2.0_rt * thr.species_err[n]) {
              std::string label = "xn_richardson_fail(" +
                                  std::string(short_spec_names_cxx[n]) + ")";
              record_mismatch(iz, label, cpu_states[iz].xn[n],
                              gpu_states[iz].xn[n]);
              zone_failed_by_richardson = true;
              ++gpu_failure_summary.species_failures;
            }
          }
        }
      } else {
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

      record_difference("temperature", cpu_states[iz].T, gpu_states[iz].T,
                        true);
      record_difference("e", cpu_states[iz].e, gpu_states[iz].e, false);
      for (int n = 0; n < NumSpec; ++n) {
        std::string norm_label =
            "number_density(" + std::string(short_spec_names_cxx[n]) + ")";
        record_difference(norm_label, cpu_states[iz].xn[n],
                          gpu_states[iz].xn[n], true);
      }

      if (use_richardson_gpu_failure && zone_failed_by_richardson) {
        gpu_failure_detected = true;
        ++gpu_failure_summary.zones_failed;
      }
    }
  }

  if (!field_norms.empty()) {
    std::cout << "CPU/GPU multi-zone difference norms:" << std::endl;
    auto safe_norm_ratio = [](amrex::Real num, amrex::Real denom) {
      if (denom > 0.0_rt) {
        return num / denom;
      }
      return (num == 0.0_rt) ? 0.0_rt
                             : std::numeric_limits<amrex::Real>::infinity();
    };
    auto print_norm = [&](const std::string &label,
                          const NormAccumulator &acc) {
      amrex::Real diff_l2 = std::sqrt(acc.diff_l2_sq);
      amrex::Real ref_l2 = std::sqrt(acc.ref_l2_sq);
      std::cout << "  " << label
                << ": L1 = " << safe_norm_ratio(acc.diff_l1, acc.ref_l1)
                << ", L2 = " << safe_norm_ratio(diff_l2, ref_l2)
                << ", L_inf = "
                << safe_norm_ratio(acc.diff_linf, acc.ref_linf) << std::endl;
    };
    auto temp_it = field_norms.find("temperature");
    if (temp_it != field_norms.end()) {
      print_norm("temperature", temp_it->second);
    }
    for (int n = 0; n < NumSpec; ++n) {
      std::string label =
          "number_density(" + std::string(short_spec_names_cxx[n]) + ")";
      auto it = field_norms.find(label);
      if (it != field_norms.end()) {
        print_norm(label, it->second);
      }
    }
  } else if (enable_gpu) {
    std::cout << "No comparable CPU/GPU multi-zone states to compute norms"
              << std::endl;
  }

  if (mismatch_summary.empty() && enable_gpu) {
    if (use_richardson_gpu_failure) {
      std::cout << "CPU/GPU multi-zone states pass Richardson-based GPU "
                   "failure criterion"
                << std::endl;
    } else {
      std::cout << "CPU/GPU multi-zone states agree within tolerances"
                << std::endl;
    }
  } else if (enable_gpu) {
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
    if (use_richardson_gpu_failure) {
      std::cout << "GPU failure summary (|CPU-GPU| > 2*|Richardson error|): "
                << "zones_failed=" << gpu_failure_summary.zones_failed
                << ", temp_failures=" << gpu_failure_summary.temp_failures
                << ", energy_failures=" << gpu_failure_summary.energy_failures
                << ", species_failures=" << gpu_failure_summary.species_failures
                << ", missing_richardson="
                << gpu_failure_summary.missing_richardson << std::endl;
    }
  }

  if (!enable_gpu) {
    return cpu_success;
  }
  if (use_richardson_gpu_failure) {
    return (cpu_success && gpu_success && !gpu_failure_detected);
  }
  return (cpu_success && gpu_success && mismatch_summary.empty());
}
#endif
