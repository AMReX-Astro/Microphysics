#include <cmath>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include <AMReX.H>
#include <AMReX_ParmParse.H>

#ifdef AMREX_USE_CUDA
#include <cuda_runtime.h>
#endif

#include <actual_network.H>
#include <burner.H>
#include <eos.H>
#include <extern_parameters.H>
#include <network.H>
#include <primordial_rodas5p_hopper.H>
#include <unit_test.H>

using namespace amrex;

namespace {

void check_cuda (const cudaError_t err, const char* what)
{
#ifdef AMREX_USE_CUDA
    if (err != cudaSuccess) {
        amrex::Abort(std::string(what) + ": " + cudaGetErrorString(err));
    }
#else
    amrex::ignore_unused(err, what);
#endif
}

burn_t make_state ()
{
    burn_t state{};

    const amrex::Real numdens[NumSpec] = {
        unit_test_rp::primary_species_1,
        unit_test_rp::primary_species_2,
        unit_test_rp::primary_species_3,
        unit_test_rp::primary_species_4,
        unit_test_rp::primary_species_5,
        unit_test_rp::primary_species_6,
        unit_test_rp::primary_species_7,
        unit_test_rp::primary_species_8,
        unit_test_rp::primary_species_9,
        unit_test_rp::primary_species_10,
        unit_test_rp::primary_species_11,
        unit_test_rp::primary_species_12,
        unit_test_rp::primary_species_13,
        unit_test_rp::primary_species_14
    };

    state.T = unit_test_rp::temperature;

    amrex::Real rhotot = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        state.xn[n] = numdens[n];
        rhotot += state.xn[n] * spmasses[n];
    }

    state.rho = rhotot;

    amrex::Real msum = 0.0_rt;
    amrex::Real mfracs[NumSpec]{};
    for (int n = 0; n < NumSpec; ++n) {
        mfracs[n] = state.xn[n] * spmasses[n] / rhotot;
        msum += mfracs[n];
    }

    for (int n = 0; n < NumSpec; ++n) {
        mfracs[n] /= msum;
        state.xn[n] = mfracs[n] * rhotot / spmasses[n];
    }

    eos(eos_input_rt, state);

    return state;
}

bool close_enough (const amrex::Real a, const amrex::Real b)
{
    constexpr amrex::Real rel_tol = 1.0e-10_rt;
    constexpr amrex::Real abs_tol = 1.0e-30_rt;
    return std::abs(a - b) <= abs_tol + rel_tol * amrex::max(std::abs(a), std::abs(b));
}

bool compare_state (const burn_t& ref, const burn_t& got)
{
    bool ok = true;

    ok = ok && ref.success == got.success;
    ok = ok && ref.error_code == got.error_code;
    ok = ok && ref.n_rhs == got.n_rhs;
    ok = ok && ref.n_jac == got.n_jac;
    ok = ok && ref.n_step == got.n_step;
    ok = ok && close_enough(ref.T, got.T);
    ok = ok && close_enough(ref.e, got.e);
    ok = ok && close_enough(ref.time, got.time);

    for (int n = 0; n < NumSpec; ++n) {
        ok = ok && close_enough(ref.xn[n], got.xn[n]);
    }

    if (! ok) {
        std::cout << "reference: success=" << ref.success
                  << " error_code=" << ref.error_code
                  << " n_rhs=" << ref.n_rhs
                  << " n_jac=" << ref.n_jac
                  << " n_step=" << ref.n_step
                  << " T=" << ref.T
                  << " e=" << ref.e
                  << " time=" << ref.time << std::endl;
        std::cout << "hopper:    success=" << got.success
                  << " error_code=" << got.error_code
                  << " n_rhs=" << got.n_rhs
                  << " n_jac=" << got.n_jac
                  << " n_step=" << got.n_step
                  << " T=" << got.T
                  << " e=" << got.e
                  << " time=" << got.time << std::endl;
        for (int n = 0; n < NumSpec; ++n) {
            if (! close_enough(ref.xn[n], got.xn[n])) {
                std::cout << "species " << n << " ref=" << ref.xn[n]
                          << " hopper=" << got.xn[n] << std::endl;
            }
        }
    }

    return ok;
}

} // namespace

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    int success = 0;

    {
        ParmParse const pp("unit_test");
        std::string const run_prefix = "burn_cell_primordial_chem_";
        std::string input_run_prefix;
        pp.query("run_prefix", input_run_prefix);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(run_prefix == input_run_prefix,
                                         "input file is missing or incorrect!");

        init_unit_test();
        eos_init(unit_test_rp::small_temp, unit_test_rp::small_dens);
        network_init();

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(integrator_rp::rosenbrock_tableau == 0,
                                         "Hopper prototype requires Rodas5P (integrator.rosenbrock_tableau = 0)");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(integrator_rp::jacobian == 1,
                                         "Hopper prototype requires analytic Jacobian (integrator.jacobian = 1)");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(integrator_rp::use_number_densities,
                                         "Hopper prototype requires number-density mode");

        constexpr int nstates = primordial_rodas5p_hopper::zones_per_block;
        const amrex::Real dt = unit_test_rp::tmax / 1000.0_rt;

        std::vector<burn_t> reference(nstates);
        std::vector<burn_t> hopper(nstates);

        for (int i = 0; i < nstates; ++i) {
            reference[i] = make_state();
            hopper[i] = reference[i];
            reference[i].i = i;
            hopper[i].i = i;
            burner(reference[i], dt);
        }

#ifdef AMREX_USE_CUDA
        burn_t* d_states = nullptr;
        check_cuda(cudaMalloc(&d_states, nstates * sizeof(burn_t)), "cudaMalloc(d_states)");
        check_cuda(cudaMemcpy(d_states, hopper.data(), nstates * sizeof(burn_t), cudaMemcpyHostToDevice),
                   "cudaMemcpy H2D states");
        primordial_rodas5p_hopper::burn_device(d_states, dt, nstates);
        check_cuda(cudaDeviceSynchronize(), "cudaDeviceSynchronize");
        check_cuda(cudaMemcpy(hopper.data(), d_states, nstates * sizeof(burn_t), cudaMemcpyDeviceToHost),
                   "cudaMemcpy D2H states");
        check_cuda(cudaFree(d_states), "cudaFree(d_states)");
#else
        amrex::Abort("Hopper prototype test requires AMREX_USE_CUDA");
#endif

        success = 1;
        for (int i = 0; i < nstates; ++i) {
            if (! compare_state(reference[i], hopper[i])) {
                success = 0;
            }
        }

        std::cout << "Hopper prototype shared bytes = "
                  << primordial_rodas5p_hopper::shared_bytes() << std::endl;
        std::cout << "Hopper prototype compare success = " << success << std::endl;
    }

    amrex::Finalize();

    return success ? 0 : 1;
}
