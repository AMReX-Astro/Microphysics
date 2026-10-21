// Host validation of the OCML f64 pow/log10/tanh transcriptions in
// f64_math_extra.H. See README.md in this directory for how to build and run.
//
// For every sweep it reports the ulp distance of the recipe from the correctly
// rounded result (MPFR, 53 bits, IEEE subnormals), and of the macOS libm from
// the same reference, plus the distance between the recipe and libm. Special
// cases (zeros, infinities, NaNs, +-1, domain errors, overflow, underflow) must
// agree with libm exactly (any NaN matches any NaN).
#define DEVICE inline
#include "f64_math_extra.H"

#include <mpfr.h>

#include <cinttypes>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>
#include <random>
#include <string>
#include <vector>

namespace fm = chemistry::f64_math;

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

// Monotone map of doubles onto integers; ulp distance is the difference.
static int64_t ordered(double x) {
    int64_t i = (int64_t)bits_of(x);
    return i < 0 ? (int64_t)0x8000000000000000LL - i : i;
}
static uint64_t ulp_distance(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) {
        return (std::isnan(a) && std::isnan(b)) ? 0 : UINT64_MAX;
    }
    if (a == b) return 0; // +0 == -0
    int64_t d = ordered(a) - ordered(b);
    return d < 0 ? (uint64_t)(-d) : (uint64_t)d;
}

// Correctly rounded references.
static double cr_pow(double x, double y) {
    mpfr_t a, b, r;
    mpfr_inits2(53, a, b, r, (mpfr_ptr)0);
    mpfr_set_d(a, x, MPFR_RNDN);
    mpfr_set_d(b, y, MPFR_RNDN);
    int t = mpfr_pow(r, a, b, MPFR_RNDN);
    t = mpfr_subnormalize(r, t, MPFR_RNDN);
    double z = mpfr_get_d(r, MPFR_RNDN);
    mpfr_clears(a, b, r, (mpfr_ptr)0);
    return z;
}
static double cr_unary(double x, int (*f)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t)) {
    mpfr_t a, r;
    mpfr_inits2(53, a, r, (mpfr_ptr)0);
    mpfr_set_d(a, x, MPFR_RNDN);
    int t = f(r, a, MPFR_RNDN);
    t = mpfr_subnormalize(r, t, MPFR_RNDN);
    double z = mpfr_get_d(r, MPFR_RNDN);
    mpfr_clears(a, r, (mpfr_ptr)0);
    return z;
}

struct Stats {
    std::string name;
    uint64_t n = 0;
    uint64_t hist[4] = {0, 0, 0, 0}; // recipe vs CR: 0, 1, 2, >=3 ulp
    uint64_t max_ulp = 0;
    double worst_x = 0, worst_y = 0, worst_got = 0, worst_ref = 0;
    uint64_t libm_max_ulp = 0;   // libm vs CR
    uint64_t libm_exact = 0;     // libm == CR
    uint64_t vs_libm_max = 0;    // recipe vs libm
    uint64_t vs_libm_equal = 0;  // recipe == libm (bits)
    uint64_t special_mismatch = 0;
};

static uint64_t g_hash = 1469598103934665603ULL; // FNV-1a of all outputs
static void hash_add(double v) {
    uint64_t u = bits_of(v);
    for (int i = 0; i < 8; ++i) {
        g_hash ^= (u >> (8 * i)) & 0xff;
        g_hash *= 1099511628211ULL;
    }
}

static bool is_special_result(double v) {
    return std::isnan(v) || std::isinf(v) || v == 0.0;
}

static void record(Stats& s, double x, double y, double got, double libm, double ref) {
    hash_add(got);
    ++s.n;
    if (is_special_result(ref) || is_special_result(libm) || is_special_result(got)) {
        // Overflow, underflow to zero, NaN: must agree with libm exactly.
        bool same = (std::isnan(got) && std::isnan(libm)) || bits_of(got) == bits_of(libm);
        if (!same) {
            if (s.special_mismatch < 10) {
                std::printf("  SPECIAL MISMATCH %s x=%a y=%a got=%a libm=%a cr=%a\n",
                            s.name.c_str(), x, y, got, libm, ref);
            }
            ++s.special_mismatch;
        }
    }
    uint64_t u = ulp_distance(got, ref);
    s.hist[u >= 3 ? 3 : u]++;
    if (u > s.max_ulp || (s.n == 1)) {
        if (u > s.max_ulp) s.max_ulp = u;
        s.worst_x = x; s.worst_y = y; s.worst_got = got; s.worst_ref = ref;
    }
    uint64_t ul = ulp_distance(libm, ref);
    if (ul > s.libm_max_ulp) s.libm_max_ulp = ul;
    if (ul == 0) ++s.libm_exact;
    uint64_t ug = ulp_distance(got, libm);
    if (ug > s.vs_libm_max) s.vs_libm_max = ug;
    if (bits_of(got) == bits_of(libm) || (std::isnan(got) && std::isnan(libm))) ++s.vs_libm_equal;
}

static void print(const Stats& s) {
    std::printf("%-28s n=%-8" PRIu64 " max_ulp=%-3" PRIu64
                " [0:%.5f 1:%.5f 2:%.6f >=3:%" PRIu64 "]  libm_max=%" PRIu64
                " libm_cr=%.5f  vs_libm_max=%" PRIu64 " vs_libm_eq=%.5f  special_mismatch=%" PRIu64 "\n",
                s.name.c_str(), s.n, s.max_ulp, double(s.hist[0]) / s.n,
                double(s.hist[1]) / s.n, double(s.hist[2]) / s.n, s.hist[3],
                s.libm_max_ulp, double(s.libm_exact) / s.n, s.vs_libm_max,
                double(s.vs_libm_equal) / s.n, s.special_mismatch);
    if (s.max_ulp > 0) {
        std::printf("    worst: x=%a y=%a got=%a cr=%a\n", s.worst_x, s.worst_y,
                    s.worst_got, s.worst_ref);
    }
}

static std::mt19937_64 rng(20261018);
static double U(double a, double b) {
    return std::uniform_real_distribution<double>(a, b)(rng);
}
static double log_uniform(double lo, double hi) {
    return std::exp2(U(std::log2(lo), std::log2(hi)));
}

static int g_fail = 0;

static void pow_sweep(const char* name, long n, const std::function<void(double&, double&)>& gen) {
    Stats s;
    s.name = std::string("pow/") + name;
    for (long i = 0; i < n; ++i) {
        double x, y;
        gen(x, y);
        record(s, x, y, fm::pow(x, y), std::pow(x, y), cr_pow(x, y));
    }
    print(s);
    if (s.max_ulp > 2 || s.special_mismatch) g_fail = 1;
}

static void unary_sweep(const char* fname, const char* name, long n, double (*recipe)(double),
                        double (*libm)(double), int (*mp)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t),
                        const std::function<double()>& gen) {
    Stats s;
    s.name = std::string(fname) + "/" + name;
    for (long i = 0; i < n; ++i) {
        double x = gen();
        record(s, x, 0.0, recipe(x), libm(x), cr_unary(x, mp));
    }
    print(s);
    if (s.max_ulp > 2 || s.special_mismatch) g_fail = 1;
}

static double r_exp(double x) { return fm::exp(x); }
static double r_log(double x) { return fm::log(x); }
static double r_cbrt(double x) { return fm::cbrt(x); }
static double l_exp(double x) { return std::exp(x); }
static double l_log(double x) { return std::log(x); }
static double l_cbrt(double x) { return std::cbrt(x); }
static double r_log10(double x) { return fm::log10(x); }
static double r_tanh(double x) { return fm::tanh(x); }
static double l_log10(double x) { return std::log10(x); }
static double l_tanh(double x) { return std::tanh(x); }

// Special-value grid: every result must match libm exactly where C/IEEE
// define it (NaN matches NaN); other grid points go into the ulp stats.
static void pow_special_grid() {
    const double inf = INFINITY, nan = NAN;
    const double sub = 0x1p-1074, big = 0x1.fffffffffffffp1023;
    std::vector<double> xs = {0.0, -0.0, sub, -sub, 0x1p-1022, 0.5, -0.5, 1.0, -1.0,
                              2.0, -2.0, 3.0, -3.0, 2.5, -2.5, 10.0, 0x1.0000000000001p0,
                              0x1.fffffffffffffp-1, big, -big, inf, -inf, nan, -nan};
    std::vector<double> ys = {0.0, -0.0, 0.5, -0.5, 1.0, -1.0, 2.0, -2.0, 3.0, -3.0,
                              2.5, -2.5, 0x1p53, -0x1p53, 0x1p53 + 2, 0x1.fffffffffffffp52,
                              -0x1.fffffffffffffp52, 1023.5, 1e300, -1e300, 0x1p-1074, -0x1p-1074,
                              1075.0, -1075.0, 1024.0, inf, -inf, nan, -nan};
    Stats s;
    s.name = "pow/special_grid";
    uint64_t exact_checked = 0, exact_mismatch = 0;
    for (double x : xs) {
        for (double y : ys) {
            double got = fm::pow(x, y), lib = std::pow(x, y), ref = cr_pow(x, y);
            record(s, x, y, got, lib, ref);
            bool defined = std::isnan(x) || std::isnan(y) || std::isinf(x) || std::isinf(y) ||
                           x == 0.0 || y == 0.0 || x == 1.0 || std::fabs(x) == 1.0 ||
                           (x < 0.0 && std::trunc(y) != y) || is_special_result(lib);
            if (defined) {
                ++exact_checked;
                bool same = (std::isnan(got) && std::isnan(lib)) || bits_of(got) == bits_of(lib);
                if (!same) {
                    ++exact_mismatch;
                    std::printf("  GRID MISMATCH x=%a y=%a got=%a libm=%a\n", x, y, got, lib);
                }
            }
        }
    }
    print(s);
    std::printf("    special grid: %" PRIu64 " C-defined cases checked bitwise, %" PRIu64
                " mismatches\n", exact_checked, exact_mismatch);
    if (exact_mismatch || s.max_ulp > 2 || s.special_mismatch) g_fail = 1;
}

// C/Annex F special cases (NaN, inf, zero results or arguments, +-1) must
// match libm bit for bit. Other listed points (e.g. tiny |x| where tanh(x)
// rounds to x) are judged against MPFR, not libm: glibc's tanh is 1 ulp off
// at 0x1.fffffffffffffp-28 while the recipe is correctly rounded there
// (Codex E2E run, 2026-10-19). A libm-only difference is reported, not failed.
static void unary_special(const char* fname, double (*recipe)(double), double (*libm)(double),
                          int (*mp)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t)) {
    const double inf = INFINITY, nan = NAN;
    std::vector<double> xs = {0.0, -0.0, 0x1p-1074, -0x1p-1074, 0x1p-1022, 1.0, -1.0, 10.0,
                              0.1, 1e-300, 1e300, 0x1.fffffffffffffp1023, inf, -inf, nan, -nan,
                              -2.0, 0x1p-27, 0x1.fffffffffffffp-28, 19.0625, 0x1.3000000000001p4,
                              -19.0625, 30.0, -30.0, 710.0, -710.0};
    uint64_t mismatches = 0;
    for (double x : xs) {
        double got = recipe(x), lib = libm(x);
        hash_add(got);
        bool annex_f = std::isnan(lib) || std::isinf(lib) || lib == 0.0 || x == 0.0 ||
                       std::isinf(x) || std::isnan(x) || std::fabs(lib) == 1.0;
        bool tiny = std::fabs(x) < 0x1p-27;
        bool same = (std::isnan(got) && std::isnan(lib)) || bits_of(got) == bits_of(lib);
        bool nan_sign_ok = !std::isnan(got) || std::signbit(got) == std::signbit(x) ||
                           std::strcmp(fname, "log10") == 0;
        bool ok = nan_sign_ok;
        if (annex_f) {
            ok = ok && same;
        } else if (tiny) {
            const double cr = cr_unary(x, mp);
            ok = ok && bits_of(got) == bits_of(cr);
            if (bits_of(lib) != bits_of(cr)) {
                std::printf("  %s note: libm is not correctly rounded at x=%a (libm=%a cr=%a got=%a)\n",
                            fname, x, lib, cr, got);
            }
        }
        if (!ok) {
            ++mismatches;
            std::printf("  %s SPECIAL MISMATCH x=%a got=%a libm=%a\n", fname, x, got, lib);
        }
    }
    std::printf("%s special values: %zu checked, %" PRIu64 " mismatches\n", fname, xs.size(),
                mismatches);
    if (mismatches) g_fail = 1;
}

int main(int argc, char** argv) {
    long N = argc > 1 ? std::atol(argv[1]) : 1000000;
    mpfr_set_emin(-1073);
    mpfr_set_emax(1024);

    std::printf("== pow ==\n");
    pow_special_grid();
    pow_sweep("wide", N, [](double& x, double& y) {
        x = log_uniform(0x1p-1074, 0x1p1023) * (U(0, 1) < 0.5 ? 1.0 : 2.0);
        y = U(-1, 1) * log_uniform(1e-3, 1e3);
    });
    pow_sweep("bits_random", N, [](double& x, double& y) {
        uint64_t a = rng() & 0x7fefffffffffffffULL, b = rng();
        if (((b >> 52) & 0x7ff) == 0x7ff) b &= ~(1ULL << 62);
        x = from_bits(a);
        y = from_bits(b);
    });
    pow_sweep("near_one_base", N, [](double& x, double& y) {
        x = 1.0 + (U(0, 1) < 0.5 ? -1.0 : 1.0) * std::exp2(-U(1, 53));
        y = (U(0, 1) < 0.5 ? -1.0 : 1.0) * log_uniform(1.0, 1e16);
    });
    pow_sweep("integer_y_signed_x", N, [](double& x, double& y) {
        x = (U(0, 1) < 0.5 ? -1.0 : 1.0) * log_uniform(1e-6, 1e6);
        y = (double)std::uniform_int_distribution<int>(-60, 60)(rng);
    });
    pow_sweep("half_integer_y", N, [](double& x, double& y) {
        x = log_uniform(1e-30, 1e30);
        y = (double)std::uniform_int_distribution<int>(-20, 20)(rng) + 0.5;
    });
    {
        const double exps[] = {2.0 / 3.0, 1.0 / 3.0, -3.0 / 2.0, -2.0 / 3.0, 5.0 / 6.0, 0.53,
                               3.263, 9.227, 0.25, 1.9885, 0.10, 5.0 / 3.0, 5.0 / 12.0};
        pow_sweep("rates_T9_1e-3..1e6", N, [&](double& x, double& y) {
            x = log_uniform(1e-3, 1e6);
            y = exps[std::uniform_int_distribution<int>(0, 12)(rng)];
        });
        pow_sweep("rates_triplealf_args", N, [&](double& x, double& y) {
            double t9 = log_uniform(1e-3, 1e6);
            bool first = U(0, 1) < 0.5;
            x = first ? 0.025 * (1.0 / t9) : t9 / 0.025;
            y = first ? 3.263 : 9.227;
        });
    }
    pow_sweep("vode_eta_root", N, [](double& x, double& y) {
        x = log_uniform(1e-14, 1e8);
        y = 1.0 / (double)std::uniform_int_distribution<int>(1, 13)(rng);
    });
    pow_sweep("vode_h_ratio_int", N, [](double& x, double& y) {
        x = log_uniform(1e-3, 1e3);
        y = (double)std::uniform_int_distribution<int>(1, 13)(rng);
    });
    pow_sweep("pow10_table", N, [](double& x, double& y) {
        x = 10.0;
        y = U(-12, 15);
    });
    pow_sweep("overflow_underflow_edges", N, [](double& x, double& y) {
        x = log_uniform(0.5, 4.0);
        double target = (U(0, 1) < 0.5) ? U(1015, 1030) : -U(1015, 1080);
        y = target / std::log2(x);
    });

    {
        // OCML does not special-case y == 1: report how often pow(x, 1) != x
        // (informational; VODE calls pow(x, 1.0 / L) and pow(x, L) with L = 1).
        long n1 = 0, m = 200000;
        for (long i = 0; i < m; ++i) {
            double x = log_uniform(1e-300, 1e300);
            if (fm::pow(x, 1.0) != x) ++n1;
        }
        std::printf("pow(x, 1.0) != x for %ld of %ld random x (e.g. pow(3, 1) = %a)\n", n1, m,
                    fm::pow(3.0, 1.0));
    }

    std::printf("== log10 ==\n");
    unary_special("log10", r_log10, l_log10, mpfr_log10);
    unary_sweep("log10", "bits_random_positive", N, r_log10, l_log10, mpfr_log10, [] {
        return from_bits(rng() & 0x7fefffffffffffffULL);
    });
    unary_sweep("log10", "near_one", N, r_log10, l_log10, mpfr_log10, [] {
        return 1.0 + (U(0, 1) < 0.5 ? -1.0 : 1.0) * std::exp2(-U(1, 53));
    });
    unary_sweep("log10", "T_and_rho_1e-5..1e16", N, r_log10, l_log10, mpfr_log10, [] {
        return log_uniform(1e-5, 1e16);
    });
    {
        // Exact powers of ten: does log10(10^k) return exactly k?
        int not_exact = 0, n = 0;
        for (int k = -22; k <= 22; ++k) {
            double x = std::pow(10.0, k); // exact for |k| <= 22
            ++n;
            if (fm::log10(x) != (double)k) {
                ++not_exact;
                std::printf("  log10(1e%d) = %a (libm %a)\n", k, fm::log10(x), std::log10(x));
            }
        }
        std::printf("log10 exact powers of ten: %d/%d not exactly k\n", not_exact, n);
    }

    std::printf("== tanh ==\n");
    unary_special("tanh", r_tanh, l_tanh, mpfr_tanh);
    unary_sweep("tanh", "-25..25", N, r_tanh, l_tanh, mpfr_tanh, [] { return U(-25, 25); });
    unary_sweep("tanh", "small_2^-40..1", N, r_tanh, l_tanh, mpfr_tanh, [] {
        return (U(0, 1) < 0.5 ? -1.0 : 1.0) * std::exp2(-U(0, 40));
    });
    unary_sweep("tanh", "bits_random", N, r_tanh, l_tanh, mpfr_tanh, [] {
        uint64_t b = rng();
        if (((b >> 52) & 0x7ff) == 0x7ff) b &= ~(1ULL << 62);
        return from_bits(b);
    });

    // The base recipes from f64_math.H (exp, log, cbrt), swept directly; above
    // they are only exercised inside pow/log10/tanh.
    std::printf("== exp / log / cbrt (f64_math.H) ==\n");
    unary_sweep("exp", "-745..710", N, r_exp, l_exp, mpfr_exp, [] { return U(-745, 710); });
    unary_sweep("exp", "-1..1", N, r_exp, l_exp, mpfr_exp, [] { return U(-1, 1); });
    unary_sweep("log", "bits_random_positive", N, r_log, l_log, mpfr_log, [] {
        return from_bits(rng() & 0x7fefffffffffffffULL);
    });
    unary_sweep("log", "near_one", N, r_log, l_log, mpfr_log, [] {
        return 1.0 + (U(0, 1) < 0.5 ? -1.0 : 1.0) * std::exp2(-U(1, 53));
    });
    unary_sweep("cbrt", "bits_random", N, r_cbrt, l_cbrt, mpfr_cbrt, [] {
        uint64_t b = rng();
        if (((b >> 52) & 0x7ff) == 0x7ff) b &= ~(1ULL << 62);
        return from_bits(b);
    });
    unary_sweep("cbrt", "1e-30..1e30", N, r_cbrt, l_cbrt, mpfr_cbrt, [] {
        return log_uniform(1e-30, 1e30);
    });

    std::printf("output hash: %016" PRIx64 "\n", g_hash);
    std::printf("%s\n", g_fail ? "FAIL" : "PASS");
    return g_fail;
}
