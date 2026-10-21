// Standalone burn reproducer (see README.md).
//
// One HIP program, no Microphysics or AMReX: it compiles the port's kernel
// tree (burn_kernel.cpp and the headers next to it) with hipcc as an
// ordinary __global__ kernel (-DMPL_HIP_STANDALONE), reads a case file
// written by test_react (loom_case_format.H: BurnParams, helmholtz tables,
// initial state, stock reference result), burns every zone on the GPU and
// compares the result with the reference bit for bit. Optionally it also
// launches a Loom code object of the same kernel root on the same inputs
// (--hsaco), so hipcc, Loom and the stock Microphysics build can be compared
// from one process on one GPU.
//
// Comparison rule: every double is compared by its bit pattern, except that
// any NaN matches any NaN; n_rhs, n_step and the per-zone success flag must
// match exactly.
//
// Exit status: 0 if every comparison that ran is bit-exact and every
// repeated run is identical, 1 on a mismatch, 2 on a usage or runtime error.
#include <hip/hip_runtime.h>

#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include "loom_case_format.H"
#include "burn_kernel.cpp"

#if !defined(MPL_HIP_STANDALONE)
#error "build with -DMPL_HIP_STANDALONE (see Makefile)"
#endif

#if defined(MPL_HIP_MATH_RECIPES)
#define STANDALONE_MATH "recipes (the Loom path's f64 OCML transcriptions)"
#else
#define STANDALONE_MATH "ocml (hipcc's device libm, as in the stock HIP build)"
#endif

namespace {

[[noreturn]] void die(const char* fmt, ...)
{
    std::va_list ap;
    va_start(ap, fmt);
    std::fprintf(stderr, "burn_standalone: ");
    std::vfprintf(stderr, fmt, ap);
    std::fprintf(stderr, "\n");
    va_end(ap);
    std::exit(2);
}

#define CHECK(call)                                                                   \
    do {                                                                              \
        hipError_t check_err_ = (call);                                               \
        if (check_err_ != hipSuccess) {                                               \
            die("%s:%d: %s failed: %s", __FILE__, __LINE__, #call,                    \
                hipGetErrorString(check_err_));                                       \
        }                                                                             \
    } while (0)

// ---------------------------------------------------------------------------
// case file

struct Case {
    LoomCaseHeader h{};
    std::vector<unsigned char> params;
    std::vector<double> tables[LOOM_CASE_NTABLES];
    std::vector<double> state_in;
    std::vector<double> ref_state;      // empty if no reference
    std::vector<int> ref_n_rhs;
    std::vector<int> ref_ok;

    long ncells() const { return static_cast<long>(h.nx) * h.ny * h.nz; }
    bool has_ref() const { return h.ref_kind != LOOM_CASE_REF_NONE; }
};

const char* ref_kind_name(std::int32_t k)
{
    switch (k) {
    case LOOM_CASE_REF_NONE: return "none";
    case LOOM_CASE_REF_STOCK_CPU: return "stock Microphysics (CPU)";
    case LOOM_CASE_REF_STOCK_HIP: return "stock Microphysics (HIP)";
    case LOOM_CASE_REF_STANDALONE_HIP: return "standalone hipcc kernel";
    case LOOM_CASE_REF_STANDALONE_LOOM: return "Loom code object";
    default: return "unknown";
    }
}

template <typename T>
void read_exact(std::FILE* fp, std::vector<T>& v, std::size_t n, const char* what, const char* path)
{
    v.resize(n);
    if (n != 0 && std::fread(v.data(), sizeof(T), n, fp) != n) {
        die("%s: truncated while reading %s", path, what);
    }
}

// table lengths the kernel was compiled for (eos.H: f[jmax][imax][9], ...)
std::uint64_t expected_table_len(int t)
{
    const std::uint64_t imax = mpl::helmholtz::imax;
    const std::uint64_t jmax = mpl::helmholtz::jmax;
    switch (t) {
    case 0: return jmax * imax * 9;          // f
    case 1: case 2: case 3: return jmax * imax * 4;  // dpdf, ef, xf
    case 4: return imax;                     // d
    case 5: case 6: case 7: case 8: case 9: return jmax;  // t, dt_sav, dt2_sav, dti_sav, dt2i_sav
    default: return imax;                    // dd_sav, dd2_sav, ddi_sav, dd2i_sav
    }
}

Case read_case(const char* path)
{
    std::FILE* fp = std::fopen(path, "rb");
    if (fp == nullptr) { die("cannot open %s", path); }
    Case c;
    if (std::fread(&c.h, sizeof c.h, 1, fp) != 1) { die("%s: too short for a case header", path); }
    if (std::memcmp(c.h.magic, LOOM_CASE_MAGIC, sizeof c.h.magic) != 0) { die("%s: not a case file (bad magic)", path); }
    if (c.h.version != LOOM_CASE_VERSION || c.h.header_bytes != sizeof(LoomCaseHeader)) {
        die("%s: case format version %u / header %u bytes, this program reads version %u / %zu bytes",
            path, c.h.version, c.h.header_bytes, LOOM_CASE_VERSION, sizeof(LoomCaseHeader));
    }
    c.h.description[sizeof c.h.description - 1] = '\0';
    if (c.h.burn_params_bytes != sizeof(mpl::BurnParams)) {
        die("%s: BurnParams is %" PRIu64 " bytes in the case, %zu in this kernel: "
            "the case was written for a different port or kernel version",
            path, c.h.burn_params_bytes, sizeof(mpl::BurnParams));
    }
    if (c.h.num_spec != mpl::NumSpec) {
        die("%s: NumSpec %d in the case, %d in this kernel", path, c.h.num_spec, mpl::NumSpec);
    }
    if (c.h.nx <= 0 || c.h.ny <= 0 || c.h.nz <= 0 || c.h.ncomp <= 0) { die("%s: bad box or component count", path); }
    const int maxidx = std::max({c.h.irho, c.h.itemp, c.h.ispec_old + c.h.num_spec - 1, c.h.ispec + c.h.num_spec - 1,
                                 c.h.irodot + c.h.num_spec - 1, c.h.irho_hnuc});
    if (maxidx >= c.h.ncomp) { die("%s: component index %d out of range (ncomp %d)", path, maxidx, c.h.ncomp); }

    read_exact(fp, c.params, c.h.burn_params_bytes, "BurnParams", path);
    for (int t = 0; t < LOOM_CASE_NTABLES; ++t) {
        if (c.h.table_len[t] != expected_table_len(t)) {
            die("%s: table %d has %" PRIu64 " entries, the kernel expects %" PRIu64,
                path, t, c.h.table_len[t], expected_table_len(t));
        }
        read_exact(fp, c.tables[t], c.h.table_len[t], "a helmholtz table", path);
    }
    const std::size_t n = static_cast<std::size_t>(c.ncells());
    read_exact(fp, c.state_in, n * c.h.ncomp, "state_in", path);
    if (c.has_ref()) {
        read_exact(fp, c.ref_state, n * c.h.ncomp, "ref_state", path);
        read_exact(fp, c.ref_n_rhs, 2 * n, "ref_n_rhs", path);
        read_exact(fp, c.ref_ok, n, "ref_success", path);
    }
    unsigned char extra = 0;
    if (std::fread(&extra, 1, 1, fp) != 0) { die("%s: trailing bytes after the last section", path); }
    std::fclose(fp);
    return c;
}

template <typename T>
void write_exact(std::FILE* fp, const T* p, std::size_t n, const char* path)
{
    if (n != 0 && std::fwrite(p, sizeof(T), n, fp) != n) { die("write to %s failed", path); }
}

void write_case(const char* path, const Case& c)
{
    std::FILE* fp = std::fopen(path, "wb");
    if (fp == nullptr) { die("cannot open %s for writing", path); }
    write_exact(fp, &c.h, 1, path);
    write_exact(fp, c.params.data(), c.params.size(), path);
    for (int t = 0; t < LOOM_CASE_NTABLES; ++t) { write_exact(fp, c.tables[t].data(), c.tables[t].size(), path); }
    write_exact(fp, c.state_in.data(), c.state_in.size(), path);
    if (c.has_ref()) {
        write_exact(fp, c.ref_state.data(), c.ref_state.size(), path);
        write_exact(fp, c.ref_n_rhs.data(), c.ref_n_rhs.size(), path);
        write_exact(fp, c.ref_ok.data(), c.ref_ok.size(), path);
    }
    if (std::fclose(fp) != 0) { die("closing %s failed", path); }
}

// The one-zone case of global zone (i, j, k): same params and tables, the
// zone's inputs and (if present) its reference result.
Case zone_case(const Case& c, int i, int j, int k)
{
    const int ii = i - c.h.lo_i;
    const int jj = j - c.h.lo_j;
    const int kk = k - c.h.lo_k;
    if (ii < 0 || ii >= c.h.nx || jj < 0 || jj >= c.h.ny || kk < 0 || kk >= c.h.nz) {
        die("zone (%d,%d,%d) is outside the box [%d..%d]x[%d..%d]x[%d..%d]", i, j, k,
            c.h.lo_i, c.h.lo_i + c.h.nx - 1, c.h.lo_j, c.h.lo_j + c.h.ny - 1, c.h.lo_k, c.h.lo_k + c.h.nz - 1);
    }
    const long n = c.ncells();
    const long z = ii + static_cast<long>(jj) * c.h.nx + static_cast<long>(kk) * c.h.nx * c.h.ny;
    Case s;
    s.h = c.h;
    s.h.nx = s.h.ny = s.h.nz = 1;
    s.h.lo_i = i;
    s.h.lo_j = j;
    s.h.lo_k = k;
    s.params = c.params;
    for (int t = 0; t < LOOM_CASE_NTABLES; ++t) { s.tables[t] = c.tables[t]; }
    for (int m = 0; m < c.h.ncomp; ++m) { s.state_in.push_back(c.state_in[z + m * n]); }
    if (c.has_ref()) {
        for (int m = 0; m < c.h.ncomp; ++m) { s.ref_state.push_back(c.ref_state[z + m * n]); }
        s.ref_n_rhs = {c.ref_n_rhs[z], c.ref_n_rhs[z + n]};
        s.ref_ok = {c.ref_ok[z]};
        s.h.ref_failures = (c.ref_ok[z] == 1) ? 0 : 1;
    }
    return s;
}

// ---------------------------------------------------------------------------
// burns

struct Result {
    std::vector<double> state;
    std::vector<int> n_rhs;
    std::vector<int> ok;
    float ms{0.0f};
    int failures() const
    {
        return static_cast<int>(std::count_if(ok.begin(), ok.end(), [](int v) { return v != 1; }));
    }
};

// Device copies of one case's inputs, plus output buffers.
struct DeviceCase {
    mpl::BurnParams* params{nullptr};
    double* tables[LOOM_CASE_NTABLES]{};
    double* state{nullptr};
    int* n_rhs{nullptr};
    int* ok{nullptr};

    explicit DeviceCase(const Case& c)
    {
        CHECK(hipMalloc(reinterpret_cast<void**>(&params), sizeof(mpl::BurnParams)));
        CHECK(hipMemcpy(params, c.params.data(), sizeof(mpl::BurnParams), hipMemcpyHostToDevice));
        for (int t = 0; t < LOOM_CASE_NTABLES; ++t) {
            CHECK(hipMalloc(reinterpret_cast<void**>(&tables[t]), c.tables[t].size() * sizeof(double)));
            CHECK(hipMemcpy(tables[t], c.tables[t].data(), c.tables[t].size() * sizeof(double), hipMemcpyHostToDevice));
        }
        CHECK(hipMalloc(reinterpret_cast<void**>(&state), c.state_in.size() * sizeof(double)));
        CHECK(hipMalloc(reinterpret_cast<void**>(&n_rhs), 2 * c.ncells() * sizeof(int)));
        CHECK(hipMalloc(reinterpret_cast<void**>(&ok), c.ncells() * sizeof(int)));
    }
    ~DeviceCase()
    {
        (void) hipFree(ok);
        (void) hipFree(n_rhs);
        (void) hipFree(state);
        for (double* t : tables) { (void) hipFree(t); }
        (void) hipFree(params);
    }
    DeviceCase(const DeviceCase&) = delete;
    DeviceCase& operator=(const DeviceCase&) = delete;

    // fresh inputs for one burn: state_in, n_rhs = -1, ok = -1
    void reset(const Case& c)
    {
        CHECK(hipMemcpy(state, c.state_in.data(), c.state_in.size() * sizeof(double), hipMemcpyHostToDevice));
        CHECK(hipMemset(n_rhs, 0xff, 2 * c.ncells() * sizeof(int)));
        CHECK(hipMemset(ok, 0xff, c.ncells() * sizeof(int)));
    }

    Result fetch(const Case& c, float ms) const
    {
        Result r;
        r.ms = ms;
        r.state.resize(c.state_in.size());
        r.n_rhs.resize(2 * static_cast<std::size_t>(c.ncells()));
        r.ok.resize(static_cast<std::size_t>(c.ncells()));
        CHECK(hipMemcpy(r.state.data(), state, r.state.size() * sizeof(double), hipMemcpyDeviceToHost));
        CHECK(hipMemcpy(r.n_rhs.data(), n_rhs, r.n_rhs.size() * sizeof(int), hipMemcpyDeviceToHost));
        CHECK(hipMemcpy(r.ok.data(), ok, r.ok.size() * sizeof(int), hipMemcpyDeviceToHost));
        return r;
    }
};

// The kernel arguments, each in a variable of the exact parameter type of
// mpl::burn_zones (for hipModuleLaunchKernel, which takes their addresses).
struct KernelArgs {
    const mpl::BurnParams* params;
    const mpl::Real* tab[LOOM_CASE_NTABLES];
    mpl::Real* state;
    long jstride, kstride, nstride;
    int lo_i, lo_j, lo_k, nx, ny, nz;
    int irho, itemp, ispec_old, ispec, irodot, irho_hnuc;
    int* n_rhs;
    long n_rhs_jstride, n_rhs_kstride, n_rhs_nstride;
    mpl::Real dt;
    int* zone_success;

    KernelArgs(const Case& c, const DeviceCase& d)
    {
        params = d.params;
        for (int t = 0; t < LOOM_CASE_NTABLES; ++t) { tab[t] = d.tables[t]; }
        state = d.state;
        jstride = c.h.nx;
        kstride = static_cast<long>(c.h.nx) * c.h.ny;
        nstride = c.ncells();
        lo_i = c.h.lo_i; lo_j = c.h.lo_j; lo_k = c.h.lo_k;
        nx = c.h.nx; ny = c.h.ny; nz = c.h.nz;
        irho = c.h.irho; itemp = c.h.itemp; ispec_old = c.h.ispec_old;
        ispec = c.h.ispec; irodot = c.h.irodot; irho_hnuc = c.h.irho_hnuc;
        n_rhs = d.n_rhs;
        n_rhs_jstride = jstride; n_rhs_kstride = kstride; n_rhs_nstride = nstride;
        dt = c.h.dt;
        zone_success = d.ok;
    }
};

using KernelSig = void (*)(const mpl::BurnParams*,
    const mpl::Real*, const mpl::Real*, const mpl::Real*, const mpl::Real*, const mpl::Real*, const mpl::Real*,
    const mpl::Real*, const mpl::Real*, const mpl::Real*, const mpl::Real*, const mpl::Real*, const mpl::Real*,
    const mpl::Real*, const mpl::Real*,
    mpl::Real*, long, long, long, int, int, int, int, int, int,
    int, int, int, int, int, int,
    int*, long, long, long, mpl::Real, int*);
static_assert(std::is_same_v<decltype(&mpl::burn_zones), KernelSig>,
              "mpl::burn_zones does not have the standard kernel-root signature (assets/root/burn_kernel.cpp)");

unsigned int blocks_for(long ncells) { return static_cast<unsigned int>((ncells + 127) / 128); }

Result run_hip(const Case& c, DeviceCase& d)
{
    d.reset(c);
    KernelArgs a(c, d);
    hipEvent_t start{}, stop{};
    CHECK(hipEventCreate(&start));
    CHECK(hipEventCreate(&stop));
    CHECK(hipEventRecord(start, nullptr));
    hipLaunchKernelGGL(mpl::burn_zones, dim3(blocks_for(c.ncells())), dim3(128), 0, nullptr,
                       a.params, a.tab[0], a.tab[1], a.tab[2], a.tab[3], a.tab[4], a.tab[5], a.tab[6],
                       a.tab[7], a.tab[8], a.tab[9], a.tab[10], a.tab[11], a.tab[12], a.tab[13],
                       a.state, a.jstride, a.kstride, a.nstride, a.lo_i, a.lo_j, a.lo_k, a.nx, a.ny, a.nz,
                       a.irho, a.itemp, a.ispec_old, a.ispec, a.irodot, a.irho_hnuc,
                       a.n_rhs, a.n_rhs_jstride, a.n_rhs_kstride, a.n_rhs_nstride, a.dt, a.zone_success);
    CHECK(hipGetLastError());
    CHECK(hipEventRecord(stop, nullptr));
    CHECK(hipEventSynchronize(stop));
    float ms = 0.0f;
    CHECK(hipEventElapsedTime(&ms, start, stop));
    (void) hipEventDestroy(start);
    (void) hipEventDestroy(stop);
    return d.fetch(c, ms);
}

// A Loom code object of the same root: its workgroup count is fixed at
// loom-compile time, so it is launched with that many blocks (--loom-blocks,
// default ceil(zones / 128) of the case it was compiled for); lanes past the
// case's zone count do nothing.
Result run_loom(const Case& c, DeviceCase& d, hipFunction_t fn, unsigned int blocks)
{
    d.reset(c);
    KernelArgs a(c, d);
    void* args[] = {
        &a.params,
        &a.tab[0], &a.tab[1], &a.tab[2], &a.tab[3], &a.tab[4], &a.tab[5], &a.tab[6],
        &a.tab[7], &a.tab[8], &a.tab[9], &a.tab[10], &a.tab[11], &a.tab[12], &a.tab[13],
        &a.state, &a.jstride, &a.kstride, &a.nstride,
        &a.lo_i, &a.lo_j, &a.lo_k, &a.nx, &a.ny, &a.nz,
        &a.irho, &a.itemp, &a.ispec_old, &a.ispec, &a.irodot, &a.irho_hnuc,
        &a.n_rhs, &a.n_rhs_jstride, &a.n_rhs_kstride, &a.n_rhs_nstride,
        &a.dt, &a.zone_success};
    hipEvent_t start{}, stop{};
    CHECK(hipEventCreate(&start));
    CHECK(hipEventCreate(&stop));
    CHECK(hipEventRecord(start, nullptr));
    CHECK(hipModuleLaunchKernel(fn, blocks, 1, 1, 128, 1, 1, 0, nullptr, args, nullptr));
    CHECK(hipGetLastError());
    CHECK(hipEventRecord(stop, nullptr));
    CHECK(hipEventSynchronize(stop));
    float ms = 0.0f;
    CHECK(hipEventElapsedTime(&ms, start, stop));
    (void) hipEventDestroy(start);
    (void) hipEventDestroy(stop);
    return d.fetch(c, ms);
}

// ---------------------------------------------------------------------------
// comparison

bool same_bits(double a, double b)
{
    if (a != a && b != b) { return true; }     // any NaN matches any NaN
    std::uint64_t x = 0;
    std::uint64_t y = 0;
    std::memcpy(&x, &a, sizeof x);
    std::memcpy(&y, &b, sizeof y);
    return x == y;
}

std::string comp_name(const LoomCaseHeader& h, int m)
{
    if (m == h.irho) { return "density"; }
    if (m == h.itemp) { return "temperature"; }
    if (m == h.irho_hnuc) { return "rho_hnuc"; }
    if (m >= h.ispec_old && m < h.ispec_old + h.num_spec) { return "X_old[" + std::to_string(m - h.ispec_old) + "]"; }
    if (m >= h.ispec && m < h.ispec + h.num_spec) { return "X_new[" + std::to_string(m - h.ispec) + "]"; }
    if (m >= h.irodot && m < h.irodot + h.num_spec) { return "rodot[" + std::to_string(m - h.irodot) + "]"; }
    return "comp" + std::to_string(m);
}

std::uint64_t bits(double x)
{
    std::uint64_t u = 0;
    std::memcpy(&u, &x, sizeof u);
    return u;
}

struct CompWorst {
    long differing{0};
    double max_rel{0.0};
    long zone{-1};
    double a{0.0}, b{0.0};
};

// Compare result B against result A zone by zone; print a summary line and
// the worst zone per differing component. Returns true if bit-exact.
bool compare(const char* label, const char* name_a, const Case& c, const std::vector<double>& sa,
             const std::vector<int>& na, const std::vector<int>& oa, const char* name_b, const Result& b)
{
    const long n = c.ncells();
    std::vector<CompWorst> w(static_cast<std::size_t>(c.h.ncomp));
    long identical = 0, nrhs_diff = 0, nstep_diff = 0, ok_diff = 0, first_bad = -1;
    for (long z = 0; z < n; ++z) {
        bool same = true;
        for (int m = 0; m < c.h.ncomp; ++m) {
            const double x = sa[z + m * n];
            const double y = b.state[z + m * n];
            if (same_bits(x, y)) { continue; }
            same = false;
            CompWorst& cw = w[m];
            ++cw.differing;
            double rel = (x != 0.0) ? std::abs(x - y) / std::abs(x) : std::numeric_limits<double>::infinity();
            if (!std::isfinite(rel)) { rel = std::numeric_limits<double>::infinity(); }
            if (cw.zone < 0 || rel > cw.max_rel) { cw.max_rel = rel; cw.zone = z; cw.a = x; cw.b = y; }
        }
        if (na[z] != b.n_rhs[z]) { ++nrhs_diff; same = false; }
        if (na[z + n] != b.n_rhs[z + n]) { ++nstep_diff; same = false; }
        if (oa[z] != b.ok[z]) { ++ok_diff; same = false; }
        if (same) { ++identical; } else if (first_bad < 0) { first_bad = z; }
    }
    const int fa = static_cast<int>(std::count_if(oa.begin(), oa.end(), [](int v) { return v != 1; }));
    const bool exact = identical == n;
    std::printf("COMPARE %s a=%s b=%s zones=%ld identical_zones=%ld nrhs_diff_zones=%ld nstep_diff_zones=%ld "
                "success_diff_zones=%ld failures_a=%d failures_b=%d bit_exact=%d\n",
                label, name_a, name_b, n, identical, nrhs_diff, nstep_diff, ok_diff, fa, b.failures(), exact ? 1 : 0);
    auto ijk = [&](long z, int& i, int& j, int& k) {
        i = c.h.lo_i + static_cast<int>(z % c.h.nx);
        j = c.h.lo_j + static_cast<int>((z / c.h.nx) % c.h.ny);
        k = c.h.lo_k + static_cast<int>(z / (static_cast<long>(c.h.nx) * c.h.ny));
    };
    for (int m = 0; m < c.h.ncomp; ++m) {
        const CompWorst& cw = w[m];
        if (cw.differing == 0) { continue; }
        int i = 0, j = 0, k = 0;
        ijk(cw.zone, i, j, k);
        std::printf("  %-14s differing_zones=%ld max_rel=%.3e worst zone (%d,%d,%d) rho=%.17g T=%.17g\n"
                    "                 %s=%.17g [0x%016" PRIx64 "]  %s=%.17g [0x%016" PRIx64 "]\n",
                    comp_name(c.h, m).c_str(), cw.differing, cw.max_rel, i, j, k,
                    c.state_in[cw.zone + c.h.irho * n], c.state_in[cw.zone + c.h.itemp * n],
                    name_a, cw.a, bits(cw.a), name_b, cw.b, bits(cw.b));
    }
    if (first_bad >= 0) {
        int i = 0, j = 0, k = 0;
        ijk(first_bad, i, j, k);
        std::printf("  first differing zone (%d,%d,%d): n_rhs %d/%d n_step %d/%d success %d/%d "
                    "(rerun it alone with --zone %d,%d,%d)\n",
                    i, j, k, na[first_bad], b.n_rhs[first_bad], na[first_bad + n], b.n_rhs[first_bad + n],
                    oa[first_bad], b.ok[first_bad], i, j, k);
    }
    return exact;
}

bool same_result(const Result& a, const Result& b)
{
    if (a.n_rhs != b.n_rhs || a.ok != b.ok || a.state.size() != b.state.size()) { return false; }
    for (std::size_t q = 0; q < a.state.size(); ++q) {
        if (!same_bits(a.state[q], b.state[q])) { return false; }
    }
    return true;
}

void usage()
{
    std::fprintf(stderr,
        "usage: burn_standalone CASE [options]\n"
        "  --repeats N         burns per kernel (default 2; repeats must be bit-identical)\n"
        "  --no-hip            skip the hipcc-compiled kernel\n"
        "  --hsaco FILE        also run this Loom code object of the same root\n"
        "  --kernel NAME       its kernel symbol (default mpl.burn_zones)\n"
        "  --loom-blocks N     workgroups it was compiled for (default ceil(zones/128) of CASE)\n"
        "  --zone I,J,K        burn only this zone (global AMReX index) of CASE\n"
        "  --out FILE          write CASE (or its --zone subset) with the hipcc result as reference\n"
        "  --out-loom FILE     same, with the Loom result as reference\n"
        "  --write-subset FILE write the --zone subset of CASE with its own reference and exit\n"
        "  --print-inputs      print the inputs (rho, T, X) of every zone burned, as hex floats\n");
    std::exit(2);
}

Case with_result(const Case& c, const Result& r, std::int32_t kind, const char* what)
{
    Case o = c;
    o.h.ref_kind = kind;
    o.h.ref_failures = r.failures();
    o.ref_state = r.state;
    o.ref_n_rhs = r.n_rhs;
    o.ref_ok = r.ok;
    std::string desc = std::string(what) + "; inputs from: " + c.h.description;
    std::memset(o.h.description, 0, sizeof o.h.description);
    std::strncpy(o.h.description, desc.c_str(), sizeof o.h.description - 1);
    return o;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2 || argv[1][0] == '-') { usage(); }
    const char* case_path = argv[1];
    int repeats = 2;
    bool run_hipcc = true;
    std::string hsaco, kernel = "mpl.burn_zones", out, out_loom, write_subset;
    long loom_blocks = -1;
    bool have_zone = false, print_inputs = false;
    int zi = 0, zj = 0, zk = 0;
    for (int a = 2; a < argc; ++a) {
        const std::string opt = argv[a];
        auto value = [&]() -> const char* {
            if (a + 1 >= argc) { usage(); }
            return argv[++a];
        };
        if (opt == "--repeats") { repeats = std::max(1, std::atoi(value())); }
        else if (opt == "--no-hip") { run_hipcc = false; }
        else if (opt == "--hsaco") { hsaco = value(); }
        else if (opt == "--kernel") { kernel = value(); }
        else if (opt == "--loom-blocks") { loom_blocks = std::atol(value()); }
        else if (opt == "--zone") {
            if (std::sscanf(value(), "%d,%d,%d", &zi, &zj, &zk) != 3) { usage(); }
            have_zone = true;
        }
        else if (opt == "--out") { out = value(); }
        else if (opt == "--out-loom") { out_loom = value(); }
        else if (opt == "--write-subset") { write_subset = value(); }
        else if (opt == "--print-inputs") { print_inputs = true; }
        else { usage(); }
    }

    const Case full = read_case(case_path);
    const Case c = have_zone ? zone_case(full, zi, zj, zk) : full;
    if (loom_blocks < 0) { loom_blocks = blocks_for(full.ncells()); }

    if (!write_subset.empty()) {
        if (!have_zone) { die("--write-subset needs --zone"); }
        write_case(write_subset.c_str(), c);
        std::printf("wrote %s: zone (%d,%d,%d), reference %s\n", write_subset.c_str(), zi, zj, zk,
                    ref_kind_name(c.h.ref_kind));
        return 0;
    }

    int dev = 0;
    hipDeviceProp_t prop{};
    CHECK(hipGetDevice(&dev));
    CHECK(hipGetDeviceProperties(&prop, dev));
    int runtime = 0;
    CHECK(hipRuntimeGetVersion(&runtime));
    std::printf("burn_standalone: device %s (%s), HIP runtime %d\n", prop.name, prop.gcnArchName, runtime);
    std::printf("  compiled by: %s, math: %s\n", __VERSION__, STANDALONE_MATH);
    std::printf("  case: %s, %d x %d x %d zones at (%d,%d,%d), %d components, NumSpec %d, dt %.17g\n",
                case_path, c.h.nx, c.h.ny, c.h.nz, c.h.lo_i, c.h.lo_j, c.h.lo_k, c.h.ncomp, c.h.num_spec, c.h.dt);
    std::printf("  reference: %s, %d failed zones\n  written by: %s\n", ref_kind_name(c.h.ref_kind),
                c.h.ref_failures, c.h.description);

    if (print_inputs) {
        const long n = c.ncells();
        for (long z = 0; z < n; ++z) {
            std::printf("  zone %ld: rho=%a T=%a X=", z, c.state_in[z + c.h.irho * n], c.state_in[z + c.h.itemp * n]);
            for (int s = 0; s < c.h.num_spec; ++s) {
                std::printf("%s%a", s ? "," : "", c.state_in[z + (c.h.ispec_old + s) * n]);
            }
            std::printf("\n");
        }
    }

    DeviceCase d(c);
    bool all_ok = true;

    std::vector<Result> hip_runs;
    if (run_hipcc) {
        hipFuncAttributes attr{};
        CHECK(hipFuncGetAttributes(&attr, reinterpret_cast<const void*>(&mpl::burn_zones)));
        std::printf("  hipcc kernel: %d VGPRs, %zu B private (scratch) per lane, max %d threads per block\n",
                    attr.numRegs, attr.localSizeBytes, attr.maxThreadsPerBlock);
        for (int r = 0; r < repeats; ++r) {
            hip_runs.push_back(run_hip(c, d));
            std::printf("RUN kind=hipcc run=%d ms=%.3f failures=%d\n", r, hip_runs.back().ms, hip_runs.back().failures());
        }
        bool rep = true;
        for (int r = 1; r < repeats; ++r) { rep = rep && same_result(hip_runs[0], hip_runs[r]); }
        std::printf("REPEAT kind=hipcc runs=%d identical=%d\n", repeats, rep ? 1 : 0);
        all_ok = all_ok && rep;
    }

    std::vector<Result> loom_runs;
    if (!hsaco.empty()) {
        hipModule_t mod{};
        hipFunction_t fn{};
        CHECK(hipModuleLoad(&mod, hsaco.c_str()));
        CHECK(hipModuleGetFunction(&fn, mod, kernel.c_str()));
        std::printf("  Loom code object: %s, kernel %s, %ld workgroups of 128\n", hsaco.c_str(), kernel.c_str(), loom_blocks);
        if (blocks_for(c.ncells()) > static_cast<unsigned int>(loom_blocks)) {
            die("the case has %ld zones, more than %ld workgroups of 128 cover", c.ncells(), loom_blocks);
        }
        for (int r = 0; r < repeats; ++r) {
            loom_runs.push_back(run_loom(c, d, fn, static_cast<unsigned int>(loom_blocks)));
            std::printf("RUN kind=loom run=%d ms=%.3f failures=%d\n", r, loom_runs.back().ms, loom_runs.back().failures());
        }
        bool rep = true;
        for (int r = 1; r < repeats; ++r) { rep = rep && same_result(loom_runs[0], loom_runs[r]); }
        std::printf("REPEAT kind=loom runs=%d identical=%d\n", repeats, rep ? 1 : 0);
        all_ok = all_ok && rep;
        (void) hipModuleUnload(mod);
    }

    if (c.has_ref() && !hip_runs.empty()) {
        all_ok = compare("hipcc_vs_reference", "reference", c, c.ref_state, c.ref_n_rhs, c.ref_ok, "hipcc", hip_runs[0]) && all_ok;
    }
    if (c.has_ref() && !loom_runs.empty()) {
        all_ok = compare("loom_vs_reference", "reference", c, c.ref_state, c.ref_n_rhs, c.ref_ok, "loom", loom_runs[0]) && all_ok;
    }
    if (!hip_runs.empty() && !loom_runs.empty()) {
        all_ok = compare("loom_vs_hipcc", "hipcc", c, hip_runs[0].state, hip_runs[0].n_rhs, hip_runs[0].ok, "loom", loom_runs[0]) && all_ok;
    }
    if (!c.has_ref() && (hip_runs.empty() || loom_runs.empty())) {
        std::printf("(no reference in the case and only one kernel: nothing to compare)\n");
    }

    if (!out.empty()) {
        if (hip_runs.empty()) { die("--out needs the hipcc kernel (drop --no-hip)"); }
        const std::string what = std::string("standalone hipcc kernel (") + STANDALONE_MATH + ") on " + prop.gcnArchName;
        write_case(out.c_str(), with_result(c, hip_runs[0], LOOM_CASE_REF_STANDALONE_HIP, what.c_str()));
        std::printf("wrote %s (reference: hipcc result)\n", out.c_str());
    }
    if (!out_loom.empty()) {
        if (loom_runs.empty()) { die("--out-loom needs --hsaco"); }
        const std::string what = std::string("Loom code object ") + hsaco + " on " + prop.gcnArchName;
        write_case(out_loom.c_str(), with_result(c, loom_runs[0], LOOM_CASE_REF_STANDALONE_LOOM, what.c_str()));
        std::printf("wrote %s (reference: Loom result)\n", out_loom.c_str());
    }

    std::printf("STANDALONE_RESULT %s\n", all_ok ? "PASS" : "FAIL");
    return all_ok ? 0 : 1;
}
