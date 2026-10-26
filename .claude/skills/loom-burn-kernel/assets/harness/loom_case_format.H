// Case file format for the standalone burn reproducer (see README.md).
//
// Written by test_react (loom_case.H, unit_test.loom_dump_case=<file>) and
// by the standalone driver (main.hip --out); read by the standalone driver.
// Self-contained (no AMReX, no Microphysics): this one file is copied as is
// into unit_test/test_react/ and into every standalone bundle.
//
// Layout (native byte order, x86-64 / little-endian; no padding between
// sections):
//
//   LoomCaseHeader
//   BurnParams                    burn_params_bytes bytes, the kernel's
//                                 mpl::BurnParams as the host filled it
//   14 helmholtz tables           table_len[t] doubles each, in the kernel
//                                 argument order (f, dpdf, ef, xf, d, t,
//                                 dt_sav, dt2_sav, dti_sav, dt2i_sav,
//                                 dd_sav, dd2_sav, ddi_sav, dd2i_sav)
//   state_in                      ncomp * nx * ny * nz doubles: the state box
//                                 before the burn, AMReX Array4 order (i
//                                 fastest, then j, k, component)
//   and, when ref_kind != 0, the result of one reference burn of state_in:
//   ref_state                     ncomp * nx * ny * nz doubles (same order)
//   ref_n_rhs                     2 * nx * ny * nz int32 (component 0 =
//                                 n_rhs, 1 = n_step; same order)
//   ref_success                   nx * ny * nz int32, 1 = burn succeeded
//
// A reader must check magic, version and header_bytes, then
// burn_params_bytes == sizeof(mpl::BurnParams), num_spec == mpl::NumSpec and
// the table lengths against the kernel's table dimensions.
#ifndef LOOM_CASE_FORMAT_H
#define LOOM_CASE_FORMAT_H

#include <cstdint>

#define LOOM_CASE_MAGIC "MPLCASE1"
#define LOOM_CASE_VERSION 1u
#define LOOM_CASE_NTABLES 14

// What produced the reference result stored in the file.
enum LoomCaseRefKind : std::int32_t {
    LOOM_CASE_REF_NONE = 0,
    LOOM_CASE_REF_STOCK_CPU = 1,        // stock Microphysics do_react, CPU build
    LOOM_CASE_REF_STOCK_HIP = 2,        // stock Microphysics do_react, HIP build
    LOOM_CASE_REF_STANDALONE_HIP = 3,   // standalone driver, hipcc-compiled kernel
    LOOM_CASE_REF_STANDALONE_LOOM = 4   // standalone driver, Loom code object
};

struct LoomCaseHeader {
    char magic[8];                      // LOOM_CASE_MAGIC, not NUL-terminated
    std::uint32_t version;              // LOOM_CASE_VERSION
    std::uint32_t header_bytes;         // sizeof(LoomCaseHeader)
    std::int32_t num_spec;              // NumSpec of the network
    std::int32_t ncomp;                 // state components (plot_t n_plot_comps)
    std::int32_t nx, ny, nz;            // box extents
    std::int32_t lo_i, lo_j, lo_k;      // box lower corner
    std::int32_t irho, itemp, ispec_old, ispec, irodot, irho_hnuc;  // plot_t indices
    std::int32_t ref_kind;              // LoomCaseRefKind
    std::int32_t ref_failures;          // failed zones in the reference burn
    double dt;                          // burn time (unit_test.tmax)
    std::uint64_t burn_params_bytes;    // sizeof(mpl::BurnParams) of the writer
    std::uint64_t table_len[LOOM_CASE_NTABLES];
    char description[1024];             // free text, NUL-terminated: what wrote
                                        // it (modules, compiler, commits, device)
};

static_assert(sizeof(LoomCaseHeader) == 1232, "LoomCaseHeader layout changed");

#endif
