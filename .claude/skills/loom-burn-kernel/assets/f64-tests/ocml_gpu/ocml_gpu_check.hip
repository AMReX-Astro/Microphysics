// GPU bit-for-bit check of the transcribed OCML f64 recipes
// (../../kernel/f64_math.H: exp, log, cbrt; ../../kernel/f64_math_extra.H: pow, log10,
// tanh) against hipcc's own OCML on AMDGPU (gfx90a). See README.md.
//
// For every input the kernels compute three results:
//   recipe  chemistry::f64_math::<fn>, compiled as a __device__ function with
//           device_shim/loomcxx/scalar.h standing in for the Loom ops;
//   stock   std::<fn> in device code, exactly what the stock Microphysics HIP
//           build calls. The clang HIP headers map pow/log10/tanh/log/cbrt
//           to __ocml_<fn>_f64; upstream clang 23 maps exp(double) to
//           __builtin_exp (llvm.exp.f64), which the AMDGPU backend may
//           expand itself instead of calling OCML;
//   direct  __ocml_<fn>_f64 called by name (informational; for pow the
//           compiler may merge it with the stock call, both being the same
//           pure function).
// The host compares recipe vs stock bit patterns per sweep.
//
// Build: see build_and_run.sh (hipcc -std=c++20 -O3 -ffp-contract=off
// --offload-arch=gfx90a -I device_shim -I ../../kernel).
#include <hip/hip_runtime.h>

#define DEVICE __device__
#include "f64_math_extra.H" // also includes f64_math.H

#include <chrono>
#include <cinttypes>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <random>
#include <string>
#include <vector>

namespace fm = chemistry::f64_math;

#define HIP_CHECK(call)                                                              \
    do {                                                                             \
        hipError_t err_ = (call);                                                    \
        if (err_ != hipSuccess) {                                                    \
            std::fprintf(stderr, "HIP error %s at %s:%d: %s\n", #call, __FILE__,     \
                         __LINE__, hipGetErrorString(err_));                         \
            std::exit(2);                                                            \
        }                                                                            \
    } while (0)

// ---------------------------------------------------------------- device side

enum Fn : int { F_EXP, F_LOG, F_CBRT, F_LOG10, F_TANH };

template <int F>
__global__ void k_unary(const double* __restrict__ x, double* __restrict__ r,
                        double* __restrict__ s, double* __restrict__ d, int n) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    const double a = x[i];
    if constexpr (F == F_EXP) {
        r[i] = fm::exp(a); s[i] = std::exp(a); d[i] = __ocml_exp_f64(a);
    } else if constexpr (F == F_LOG) {
        r[i] = fm::log(a); s[i] = std::log(a); d[i] = __ocml_log_f64(a);
    } else if constexpr (F == F_CBRT) {
        r[i] = fm::cbrt(a); s[i] = std::cbrt(a); d[i] = __ocml_cbrt_f64(a);
    } else if constexpr (F == F_LOG10) {
        r[i] = fm::log10(a); s[i] = std::log10(a); d[i] = __ocml_log10_f64(a);
    } else {
        r[i] = fm::tanh(a); s[i] = std::tanh(a); d[i] = __ocml_tanh_f64(a);
    }
}

// pow with both arguments loaded at run time.
__global__ void k_pow(const double* __restrict__ x, const double* __restrict__ y,
                      double* __restrict__ r, double* __restrict__ s, double* __restrict__ d,
                      int n) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    const double a = x[i], b = y[i];
    r[i] = fm::pow(a, b);
    s[i] = std::pow(a, b);
    d[i] = __ocml_pow_f64(a, b);
}

// pow with a compile-time constant exponent, as in the rate and EOS code
// (e.g. std::pow(plasg, 0.25e0_rt), std::pow(y_12, 0.14_rt)). After OCML is
// inlined the compiler may fold y-dependent parts; this checks that the
// folding changes no bits on either side.
constexpr int kNumConstExps = 16;
__host__ __device__ constexpr double const_exp(int i) {
    switch (i) {
    case 0: return 2.0 / 3.0;
    case 1: return 1.0 / 3.0;
    case 2: return -3.0 / 2.0;
    case 3: return -2.0 / 3.0;
    case 4: return 5.0 / 6.0;
    case 5: return 0.53;
    case 6: return 3.263;
    case 7: return 9.227;
    case 8: return 0.25;
    case 9: return 1.9885;
    case 10: return 0.10;
    case 11: return 5.0 / 3.0;
    case 12: return 5.0 / 12.0;
    case 13: return 0.14;
    case 14: return 0.19;
    default: return 1.0;
    }
}
constexpr int kConstOneIndex = 15; // const_exp(15) == 1.0

template <int I>
__global__ void k_pow_const_y(const double* __restrict__ x, double* __restrict__ r,
                              double* __restrict__ s, double* __restrict__ d, int n) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    constexpr double y = const_exp(I);
    const double a = x[i];
    r[i] = fm::pow(a, y);
    s[i] = std::pow(a, y);
    d[i] = __ocml_pow_f64(a, y);
}

// std::pow(10.0e0_rt, t), as in the rates/ code and the helmholtz table setup.
__global__ void k_pow_base10(const double* __restrict__ t, double* __restrict__ r,
                             double* __restrict__ s, double* __restrict__ d, int n) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    const double b = t[i];
    r[i] = fm::pow(10.0, b);
    s[i] = std::pow(10.0, b);
    d[i] = __ocml_pow_f64(10.0, b);
}

// VODE (vode_dvstep.H): std::pow(x, 1.0_rt / FLOTL) with FLOTL = (Real)L, and
// std::pow(x, vstate.L) with L a short (HIP's __HIP_OVERLOAD2 template turns
// it into pow((double)x, (double)L)). The recipe is called as the Loom
// translation calls it: pow(x, 1.0 / (Real)L) and pow(x, (Real)L).
template <bool Root>
__global__ void k_pow_vode(const double* __restrict__ x, const short* __restrict__ L,
                           double* __restrict__ r, double* __restrict__ s,
                           double* __restrict__ d, int n) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    const double a = x[i];
    const short l = L[i];
    if constexpr (Root) {
        const double flotl = (double)l;
        r[i] = fm::pow(a, 1.0 / flotl);
        s[i] = std::pow(a, 1.0 / flotl);
        d[i] = __ocml_pow_f64(a, 1.0 / flotl);
    } else {
        r[i] = fm::pow(a, (double)l);
        s[i] = std::pow(a, l);
        d[i] = __ocml_pow_f64(a, (double)l);
    }
}

// Diagnostic for question (a): IEEE 1.0 / x (what the recipes use) against
// an emulation of OCML's MATH_FAST_RCP, assumed to be v_rcp_f64 followed by
// two fma Newton steps. The exact OCML macro body is not reproduced here; if
// this reports differences, compare with ocml/src/opts.h / ep.h of the
// installed device-libs. The recipe-vs-stock sweeps are the real test.
__global__ void k_rcp(const double* __restrict__ x, double* __restrict__ ieee,
                      double* __restrict__ newton, int n) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    const double a = x[i];
    ieee[i] = 1.0 / a;
    double rr = __builtin_amdgcn_rcp(a);
    double e = __builtin_fma(-a, rr, 1.0);
    rr = __builtin_fma(e, rr, rr);
    e = __builtin_fma(-a, rr, 1.0);
    rr = __builtin_fma(e, rr, rr);
    newton[i] = rr;
}

// ------------------------------------------------------------------ host side

static uint64_t bits_of(double x) {
    uint64_t u;
    std::memcpy(&u, &x, 8);
    return u;
}
static double from_bits(uint64_t u) {
    double x;
    std::memcpy(&x, &u, 8);
    return x;
}
static int64_t ordered(double x) {
    int64_t i = (int64_t)bits_of(x);
    return i < 0 ? (int64_t)0x8000000000000000LL - i : i;
}
static uint64_t ulp_distance(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) {
        return (std::isnan(a) && std::isnan(b)) ? 0 : UINT64_MAX;
    }
    if (a == b) return 0; // +0 vs -0 is reported as a bit difference with 0 ulp
    int64_t d = ordered(a) - ordered(b);
    return d < 0 ? (uint64_t)(-d) : (uint64_t)d;
}

// Device buffers, sized once.
static size_t g_cap = 0;
static double *d_x = nullptr, *d_y = nullptr, *d_r = nullptr, *d_s = nullptr, *d_d = nullptr;
static short* d_L = nullptr;
static double g_gpu_ms = 0.0;

static void ensure_capacity(size_t n) {
    if (n <= g_cap) return;
    if (g_cap) {
        HIP_CHECK(hipFree(d_x)); HIP_CHECK(hipFree(d_y)); HIP_CHECK(hipFree(d_r));
        HIP_CHECK(hipFree(d_s)); HIP_CHECK(hipFree(d_d)); HIP_CHECK(hipFree(d_L));
    }
    HIP_CHECK(hipMalloc((void**)&d_x, n * sizeof(double)));
    HIP_CHECK(hipMalloc((void**)&d_y, n * sizeof(double)));
    HIP_CHECK(hipMalloc((void**)&d_r, n * sizeof(double)));
    HIP_CHECK(hipMalloc((void**)&d_s, n * sizeof(double)));
    HIP_CHECK(hipMalloc((void**)&d_d, n * sizeof(double)));
    HIP_CHECK(hipMalloc((void**)&d_L, n * sizeof(short)));
    g_cap = n;
}

static dim3 grid_for(int n) { return dim3((unsigned)((n + 255) / 256)); }

struct Outputs {
    std::vector<double> r, s, d;
};

// Runs `launch` (which launches one or more kernels on d_x/d_y/d_L -> d_r,
// d_s, d_d) and copies the three result arrays back.
static Outputs run_gpu(const std::vector<double>& x, const std::vector<double>* y,
                       const std::vector<short>* L, const std::function<void(int)>& launch) {
    const int n = (int)x.size();
    ensure_capacity(x.size());
    HIP_CHECK(hipMemcpy(d_x, x.data(), n * sizeof(double), hipMemcpyHostToDevice));
    if (y) HIP_CHECK(hipMemcpy(d_y, y->data(), n * sizeof(double), hipMemcpyHostToDevice));
    if (L) HIP_CHECK(hipMemcpy(d_L, L->data(), n * sizeof(short), hipMemcpyHostToDevice));
    hipEvent_t e0, e1;
    HIP_CHECK(hipEventCreate(&e0));
    HIP_CHECK(hipEventCreate(&e1));
    HIP_CHECK(hipEventRecord(e0, 0));
    launch(n);
    HIP_CHECK(hipGetLastError());
    HIP_CHECK(hipEventRecord(e1, 0));
    HIP_CHECK(hipEventSynchronize(e1));
    float ms = 0.0f;
    HIP_CHECK(hipEventElapsedTime(&ms, e0, e1));
    g_gpu_ms += ms;
    HIP_CHECK(hipEventDestroy(e0));
    HIP_CHECK(hipEventDestroy(e1));
    Outputs o;
    o.r.resize(n); o.s.resize(n); o.d.resize(n);
    HIP_CHECK(hipMemcpy(o.r.data(), d_r, n * sizeof(double), hipMemcpyDeviceToHost));
    HIP_CHECK(hipMemcpy(o.s.data(), d_s, n * sizeof(double), hipMemcpyDeviceToHost));
    HIP_CHECK(hipMemcpy(o.d.data(), d_d, n * sizeof(double), hipMemcpyDeviceToHost));
    return o;
}

static int g_fail = 0;            // some sweep has a non-NaN bit difference
static uint64_t g_nan_payload = 0; // both NaN but different bits
static uint64_t g_stock_direct = 0;

struct Example {
    uint64_t ulp;
    double x, y, r, s;
};

static void print_val(const char* tag, double v) {
    std::printf(" %s=%a(0x%016" PRIx64 ")", tag, v, bits_of(v));
}

// Compares recipe vs stock for one sweep and prints the OCMLGPU line.
static void compare(const std::string& name, const std::vector<double>& x,
                    const std::vector<double>* y, const Outputs& o) {
    const size_t n = x.size();
    uint64_t identical = 0, differ = 0, nan_payload = 0, stock_direct = 0, max_ulp = 0;
    uint64_t recipe_direct = 0;
    uint64_t ulp_hist[3] = {0, 0, 0}; // 1 ulp (or 0 ulp: signed zero), 2 ulp, >2
    const int kKeep = 4;
    std::vector<Example> worst;
    for (size_t i = 0; i < n; ++i) {
        const double r = o.r[i], s = o.s[i];
        if (bits_of(o.s[i]) != bits_of(o.d[i])) ++stock_direct;
        if (bits_of(o.r[i]) != bits_of(o.d[i])) ++recipe_direct;
        if (bits_of(r) == bits_of(s)) {
            ++identical;
            continue;
        }
        if (std::isnan(r) && std::isnan(s)) {
            ++nan_payload;
            if (worst.size() < (size_t)kKeep && max_ulp == 0) {
                worst.push_back({0, x[i], y ? (*y)[i] : 0.0, r, s});
            }
            continue;
        }
        ++differ;
        const uint64_t u = ulp_distance(r, s);
        ulp_hist[u <= 1 ? 0 : (u == 2 ? 1 : 2)]++;
        if (u > max_ulp) max_ulp = u;
        // Keep the kKeep largest-ulp examples (first seen wins ties).
        Example ex{u, x[i], y ? (*y)[i] : 0.0, r, s};
        if (worst.size() < (size_t)kKeep) {
            worst.push_back(ex);
        } else {
            size_t m = 0;
            for (size_t k = 1; k < worst.size(); ++k)
                if (worst[k].ulp < worst[m].ulp) m = k;
            if (ex.ulp > worst[m].ulp) worst[m] = ex;
        }
    }
    char ulp_str[32];
    if (max_ulp == UINT64_MAX) std::snprintf(ulp_str, sizeof ulp_str, "nan_vs_number");
    else std::snprintf(ulp_str, sizeof ulp_str, "%" PRIu64, max_ulp);
    std::printf("OCMLGPU %-32s n=%-8zu bit_identical=%-8" PRIu64 " differ=%-6" PRIu64
                " max_ulp=%s [<=1:%" PRIu64 " 2:%" PRIu64 " >2:%" PRIu64 "]"
                " nan_payload_diff=%" PRIu64 " stock_vs_direct_differ=%" PRIu64
                " recipe_vs_direct_differ=%" PRIu64 "\n",
                name.c_str(), n, identical, differ, ulp_str, ulp_hist[0], ulp_hist[1],
                ulp_hist[2], nan_payload, stock_direct, recipe_direct);
    for (const Example& ex : worst) {
        std::printf("    example ulp=%" PRIu64 ":", ex.ulp == UINT64_MAX ? 0 : ex.ulp);
        print_val("x", ex.x);
        if (y) print_val("y", ex.y);
        print_val("recipe", ex.r);
        print_val("stock", ex.s);
        std::printf("\n");
    }
    if (differ) g_fail = 1;
    g_nan_payload += nan_payload;
    g_stock_direct += stock_direct;
}

// ---------------------------------------------------------------- the sweeps

static std::mt19937_64 rng(20261018); // same seed as ../test_ocml.cpp
static double U(double a, double b) { return std::uniform_real_distribution<double>(a, b)(rng); }
static double log_uniform(double lo, double hi) {
    return std::exp2(U(std::log2(lo), std::log2(hi)));
}
static double sgn() { return U(0, 1) < 0.5 ? -1.0 : 1.0; }
static uint64_t no_nan_bits(uint64_t b) {
    if (((b >> 52) & 0x7ff) == 0x7ff) b &= ~(1ULL << 62);
    return b;
}

template <int F>
static void unary_sweep(const char* fname, const char* sname, long n,
                        const std::function<double()>& gen) {
    std::vector<double> x(n);
    for (long i = 0; i < n; ++i) x[i] = gen();
    Outputs o = run_gpu(x, nullptr, nullptr, [](int m) {
        hipLaunchKernelGGL(k_unary<F>, grid_for(m), dim3(256), 0, 0, d_x, d_r, d_s, d_d, m);
    });
    compare(std::string(fname) + "/" + sname, x, nullptr, o);
}

template <int F>
static void unary_values(const char* fname, const char* sname, const std::vector<double>& x) {
    Outputs o = run_gpu(x, nullptr, nullptr, [](int m) {
        hipLaunchKernelGGL(k_unary<F>, grid_for(m), dim3(256), 0, 0, d_x, d_r, d_s, d_d, m);
    });
    compare(std::string(fname) + "/" + sname, x, nullptr, o);
}

static void launch_pow(int m) {
    hipLaunchKernelGGL(k_pow, grid_for(m), dim3(256), 0, 0, d_x, d_y, d_r, d_s, d_d, m);
}

static void pow_values(const char* sname, const std::vector<double>& x,
                       const std::vector<double>& y) {
    Outputs o = run_gpu(x, &y, nullptr, launch_pow);
    compare(std::string("pow/") + sname, x, &y, o);
}

static void pow_sweep(const char* sname, long n,
                      const std::function<void(double&, double&)>& gen) {
    std::vector<double> x(n), y(n);
    for (long i = 0; i < n; ++i) gen(x[i], y[i]);
    pow_values(sname, x, y);
}

template <int I>
static void launch_const_y(int which, double* xp, double* rp, double* sp, double* dp, int m) {
    if (which == I) {
        hipLaunchKernelGGL(k_pow_const_y<I>, grid_for(m), dim3(256), 0, 0, xp, rp, sp, dp, m);
    } else if constexpr (I + 1 < kNumConstExps) {
        launch_const_y<I + 1>(which, xp, rp, sp, dp, m);
    }
}

static const std::vector<double> kSpecialUnary = {
    0.0, -0.0, 0x1p-1074, -0x1p-1074, 0x1p-1022, -0x1p-1022, 0x1.ffffffffffffep-1023, 1.0,
    -1.0, 10.0, 0.1, 1e-300, 1e300, 0x1.fffffffffffffp1023, -0x1.fffffffffffffp1023, INFINITY,
    -INFINITY, NAN, -NAN, from_bits(0x7ff0000000000001ULL) /* sNaN */,
    from_bits(0x7ff8000000000123ULL), from_bits(0xfff4000000000000ULL), -2.0, 0x1p-27,
    0x1.fffffffffffffp-28, 19.0625, 0x1.3000000000001p4, -19.0625, 30.0, -30.0, 32.0, 33.0,
    710.0, -710.0, 709.782712893384, 709.7827128933841, 709.79, -708.3964185322641,
    -745.1332191019411, -745.1332191019412, -745.14, 1024.0, -1075.0, 1025.0, -1076.0, 2.0, 0.5,
    8.0, 27.0, -27.0, 1e-310, 3.0, 0x1.0000000000001p0, 0x1.fffffffffffffp-1,
    0x1.5555555555555p-1, 0x1.5555555555556p-1, 1000.0, 1e22, 1e-22, 1e15, 1e16, 2.718281828459045};

int main(int argc, char** argv) {
    const long N = argc > 1 ? std::atol(argv[1]) : (1L << 20);
    const auto t0 = std::chrono::steady_clock::now();

    hipDeviceProp_t prop;
    HIP_CHECK(hipGetDeviceProperties(&prop, 0));
    std::printf("device: %s (%s), samples per random sweep: %ld\n", prop.name,
                prop.gcnArchName, N);
    ensure_capacity((size_t)N > 8192 ? (size_t)N : 8192);

    // ------------------------------------------------------------------ pow
    std::printf("== pow ==\n");
    {
        const double inf = INFINITY, nan = NAN;
        const double sub = 0x1p-1074, big = 0x1.fffffffffffffp1023;
        const std::vector<double> xs = {0.0, -0.0, sub, -sub, 0x1p-1022, 0.5, -0.5, 1.0, -1.0,
                                        2.0, -2.0, 3.0, -3.0, 2.5, -2.5, 10.0,
                                        0x1.0000000000001p0, 0x1.fffffffffffffp-1, big, -big,
                                        inf, -inf, nan, -nan};
        const std::vector<double> ys = {0.0, -0.0, 0.5, -0.5, 1.0, -1.0, 2.0, -2.0, 3.0, -3.0,
                                        2.5, -2.5, 0x1p53, -0x1p53, 0x1p53 + 2,
                                        0x1.fffffffffffffp52, -0x1.fffffffffffffp52, 1023.5,
                                        1e300, -1e300, 0x1p-1074, -0x1p-1074, 1075.0, -1075.0,
                                        1024.0, inf, -inf, nan, -nan};
        std::vector<double> x, y;
        for (double a : xs)
            for (double b : ys) { x.push_back(a); y.push_back(b); }
        pow_values("special_grid", x, y);
    }
    pow_sweep("wide", N, [](double& x, double& y) {
        x = log_uniform(0x1p-1074, 0x1p1023) * (U(0, 1) < 0.5 ? 1.0 : 2.0);
        y = U(-1, 1) * log_uniform(1e-3, 1e3);
    });
    pow_sweep("bits_random", N, [](double& x, double& y) {
        x = from_bits(rng() & 0x7fefffffffffffffULL);
        y = from_bits(no_nan_bits(rng()));
    });
    pow_sweep("bits_random_signed_x", N, [](double& x, double& y) {
        x = from_bits(no_nan_bits(rng()));
        y = from_bits(no_nan_bits(rng()));
    });
    pow_sweep("near_one_base", N, [](double& x, double& y) {
        x = 1.0 + sgn() * std::exp2(-U(1, 53));
        y = sgn() * log_uniform(1.0, 1e16);
    });
    pow_sweep("integer_y_signed_x", N, [](double& x, double& y) {
        x = sgn() * log_uniform(1e-6, 1e6);
        y = (double)std::uniform_int_distribution<int>(-60, 60)(rng);
    });
    pow_sweep("half_integer_y", N, [](double& x, double& y) {
        x = log_uniform(1e-30, 1e30);
        y = (double)std::uniform_int_distribution<int>(-20, 20)(rng) + 0.5;
    });
    pow_sweep("rates_T9_1e-3..1e6", N, [](double& x, double& y) {
        x = log_uniform(1e-3, 1e6);
        y = const_exp(std::uniform_int_distribution<int>(0, kNumConstExps - 2)(rng));
    });
    {
        // The same exponents as compile-time constants (one kernel per value).
        const long chunk = N / kNumConstExps > 0 ? N / kNumConstExps : 1;
        std::vector<double> x(chunk * kNumConstExps), y(x.size());
        for (int c = 0; c < kNumConstExps; ++c)
            for (long i = 0; i < chunk; ++i) {
                x[c * chunk + i] = log_uniform(1e-3, 1e6);
                y[c * chunk + i] = const_exp(c);
            }
        Outputs o = run_gpu(x, nullptr, nullptr, [chunk](int) {
            for (int c = 0; c < kNumConstExps; ++c) {
                const long off = c * chunk;
                launch_const_y<0>(c, d_x + off, d_r + off, d_s + off, d_d + off, (int)chunk);
            }
        });
        compare("pow/rates_T9_const_y", x, &y, o);
    }
    pow_sweep("rates_triplealf_args", N, [](double& x, double& y) {
        const double t9 = log_uniform(1e-3, 1e6);
        const bool first = U(0, 1) < 0.5;
        x = first ? 0.025 * (1.0 / t9) : t9 / 0.025;
        y = first ? 3.263 : 9.227;
    });
    {
        // VODE: x^(1/L) and x^L, L = 1..6, with L a short as in vode_dvstep.H.
        std::vector<double> x(N), y(N);
        std::vector<short> L(N);
        for (long i = 0; i < N; ++i) {
            x[i] = log_uniform(1e-14, 1e8);
            L[i] = (short)std::uniform_int_distribution<int>(1, 6)(rng);
            y[i] = 1.0 / (double)L[i];
        }
        Outputs o = run_gpu(x, nullptr, &L, [](int m) {
            hipLaunchKernelGGL(k_pow_vode<true>, grid_for(m), dim3(256), 0, 0, d_x, d_L, d_r,
                               d_s, d_d, m);
        });
        compare("pow/vode_root_x^(1/L)", x, &y, o);
        for (long i = 0; i < N; ++i) {
            x[i] = log_uniform(1e-3, 1e3);
            y[i] = (double)L[i];
        }
        o = run_gpu(x, nullptr, &L, [](int m) {
            hipLaunchKernelGGL(k_pow_vode<false>, grid_for(m), dim3(256), 0, 0, d_x, d_L, d_r,
                               d_s, d_d, m);
        });
        compare("pow/vode_int_x^L_short", x, &y, o);
    }
    {
        std::vector<double> t(N), ten(N, 10.0);
        for (long i = 0; i < N; ++i) t[i] = U(-12, 15);
        pow_values("pow10_runtime_base", ten, t);
        Outputs o = run_gpu(t, nullptr, nullptr, [](int m) {
            hipLaunchKernelGGL(k_pow_base10, grid_for(m), dim3(256), 0, 0, d_x, d_r, d_s, d_d,
                               m);
        });
        compare("pow/pow10_const_base", t, nullptr, o);
    }
    pow_sweep("overflow_underflow_edges", N, [](double& x, double& y) {
        x = log_uniform(0.5, 4.0);
        const double target = (U(0, 1) < 0.5) ? U(1015, 1030) : -U(1015, 1080);
        y = target / std::log2(x);
    });
    pow_sweep("subnormal_base", N, [](double& x, double& y) {
        x = from_bits(rng() & 0x000fffffffffffffULL);
        y = U(-1, 1) * log_uniform(1e-3, 2.0);
    });

    // pow(x, 1) == x? (OCML has no y == 1 special case.)
    {
        const long m = N < 200000 ? N : 200000;
        std::vector<double> x(m), one(m, 1.0);
        for (long i = 0; i < m; ++i) x[i] = log_uniform(1e-300, 1e300);
        Outputs o = run_gpu(x, &one, nullptr, launch_pow);
        Outputs oc = run_gpu(x, nullptr, nullptr, [](int k) {
            hipLaunchKernelGGL(k_pow_const_y<kConstOneIndex>, grid_for(k), dim3(256), 0, 0, d_x,
                               d_r, d_s, d_d, k);
        });
        long rr = 0, ss = 0, dd = 0, rc = 0, sc = 0;
        for (long i = 0; i < m; ++i) {
            rr += bits_of(o.r[i]) == bits_of(x[i]);
            ss += bits_of(o.s[i]) == bits_of(x[i]);
            dd += bits_of(o.d[i]) == bits_of(x[i]);
            rc += bits_of(oc.r[i]) == bits_of(x[i]);
            sc += bits_of(oc.s[i]) == bits_of(x[i]);
        }
        std::printf("OCMLGPU pow(x,1)==x n=%ld  runtime y: recipe=%ld stock=%ld direct=%ld"
                    "  const y: recipe=%ld stock=%ld  (not equal: recipe %.4f%%, stock %.4f%%)\n",
                    m, rr, ss, dd, rc, sc, 100.0 * (m - rr) / m, 100.0 * (m - ss) / m);
        compare("pow/x^1_runtime_y", x, &one, o);
        compare("pow/x^1_const_y", x, &one, oc);
        std::vector<double> three = {3.0};
        std::vector<double> y1 = {1.0};
        Outputs o3 = run_gpu(three, &y1, nullptr, launch_pow);
        std::printf("    pow(3,1): recipe=%a stock=%a\n", o3.r[0], o3.s[0]);
    }

    // ---------------------------------------------------------- log10, tanh
    std::printf("== log10 ==\n");
    unary_values<F_LOG10>("log10", "special_values", kSpecialUnary);
    unary_sweep<F_LOG10>("log10", "bits_random_positive", N,
                         [] { return from_bits(rng() & 0x7fefffffffffffffULL); });
    unary_sweep<F_LOG10>("log10", "near_one", N,
                         [] { return 1.0 + sgn() * std::exp2(-U(1, 53)); });
    unary_sweep<F_LOG10>("log10", "T_and_rho_1e-5..1e16", N,
                         [] { return log_uniform(1e-5, 1e16); });
    unary_sweep<F_LOG10>("log10", "subnormal", N / 8 + 1,
                         [] { return from_bits(rng() & 0x000fffffffffffffULL); });
    {
        std::vector<double> x;
        for (int k = -22; k <= 22; ++k) x.push_back(std::pow(10.0, k)); // exact on host
        unary_values<F_LOG10>("log10", "exact_powers_of_ten", x);
    }

    std::printf("== tanh ==\n");
    unary_values<F_TANH>("tanh", "special_values", kSpecialUnary);
    unary_sweep<F_TANH>("tanh", "-25..25", N, [] { return U(-25, 25); });
    unary_sweep<F_TANH>("tanh", "small_2^-40..1", N, [] { return sgn() * std::exp2(-U(0, 40)); });
    unary_sweep<F_TANH>("tanh", "bits_random", N, [] { return from_bits(no_nan_bits(rng())); });

    // -------------------------------------------------------- exp, log, cbrt
    std::printf("== exp ==\n");
    unary_values<F_EXP>("exp", "special_values", kSpecialUnary);
    unary_sweep<F_EXP>("exp", "-745..710", N, [] { return U(-745.2, 709.8); });
    unary_sweep<F_EXP>("exp", "-50..50", N, [] { return U(-50, 50); });
    unary_sweep<F_EXP>("exp", "small_2^-60..1", N, [] { return sgn() * std::exp2(-U(0, 60)); });
    unary_sweep<F_EXP>("exp", "overflow_subnormal_edges", N,
                       [] { return U(0, 1) < 0.5 ? U(700, 712) : U(-760, -700); });
    unary_sweep<F_EXP>("exp", "bits_random", N, [] { return from_bits(no_nan_bits(rng())); });

    std::printf("== log ==\n");
    unary_values<F_LOG>("log", "special_values", kSpecialUnary);
    unary_sweep<F_LOG>("log", "bits_random_positive", N,
                       [] { return from_bits(rng() & 0x7fefffffffffffffULL); });
    unary_sweep<F_LOG>("log", "near_one", N, [] { return 1.0 + sgn() * std::exp2(-U(1, 53)); });
    unary_sweep<F_LOG>("log", "1e-5..1e16", N, [] { return log_uniform(1e-5, 1e16); });
    unary_sweep<F_LOG>("log", "subnormal", N / 8 + 1,
                       [] { return from_bits(rng() & 0x000fffffffffffffULL); });

    std::printf("== cbrt ==\n");
    unary_values<F_CBRT>("cbrt", "special_values", kSpecialUnary);
    unary_sweep<F_CBRT>("cbrt", "bits_random", N, [] { return from_bits(no_nan_bits(rng())); });
    unary_sweep<F_CBRT>("cbrt", "1e-30..1e30_signed", N,
                        [] { return sgn() * log_uniform(1e-30, 1e30); });
    unary_sweep<F_CBRT>("cbrt", "0.5..2", N, [] { return U(0.5, 2.0); });
    unary_sweep<F_CBRT>("cbrt", "subnormal", N / 8 + 1,
                        [] { return from_bits(rng() & 0x800fffffffffffffULL); });
    {
        std::vector<double> x;
        for (long k = 1; k <= 100000; ++k) x.push_back((double)(k * k * k));
        unary_values<F_CBRT>("cbrt", "exact_cubes_1..1e5", x);
    }

    // ------------------------------------- diagnostic: MATH_FAST_RCP vs IEEE
    {
        std::printf("== diagnostic: IEEE 1/x vs v_rcp_f64 + 2 Newton steps (assumed "
                    "MATH_FAST_RCP form; not part of PASS) ==\n");
        const char* names[] = {"1+m_m_in_[2/3,4/3)", "tanh_e.hi_[1,e^32]", "random_normal"};
        for (int k = 0; k < 3; ++k) {
            std::vector<double> x(N);
            for (long i = 0; i < N; ++i) {
                if (k == 0) x[i] = 1.0 + U(2.0 / 3.0, 4.0 / 3.0);
                else if (k == 1) x[i] = std::exp(U(0, 32));
                else x[i] = from_bits((rng() & 0x7fdfffffffffffffULL) | 0x0010000000000000ULL);
            }
            Outputs o = run_gpu(x, nullptr, nullptr, [](int m) {
                hipLaunchKernelGGL(k_rcp, grid_for(m), dim3(256), 0, 0, d_x, d_r, d_s, m);
            });
            long diff = 0, shown = 0;
            for (long i = 0; i < N; ++i) {
                if (bits_of(o.r[i]) != bits_of(o.s[i])) {
                    ++diff;
                    if (shown++ < 3) {
                        std::printf("    rcp example:");
                        print_val("x", x[i]);
                        print_val("ieee", o.r[i]);
                        print_val("newton", o.s[i]);
                        std::printf("\n");
                    }
                }
            }
            std::printf("OCMLGPU-DIAG rcp/%s n=%ld differ=%ld\n", names[k], N, diff);
        }
    }

    const double wall =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("GPU kernel time %.1f ms, wall %.1f s\n", g_gpu_ms, wall);
    std::printf("total nan_payload_diff=%" PRIu64 " stock_vs_direct_differ=%" PRIu64 "\n",
                g_nan_payload, g_stock_direct);
    if (!g_fail && g_nan_payload == 0) {
        std::printf("PASS: every sweep is bit-identical (recipe vs stock)\n");
        return 0;
    }
    if (!g_fail) {
        std::printf("FAIL: only NaN payloads differ (%" PRIu64 " results); all non-NaN results "
                    "bit-identical\n", g_nan_payload);
        return 1;
    }
    std::printf("FAIL: some sweep has non-NaN bit differences (see differ= above)\n");
    return 1;
}
