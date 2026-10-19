# Bringing up a new network

For a network in the `networks/rhs.H` framework (it defines `rhs_data` in
`actual_network.H`; list them with
`grep -l rhs_data networks/*/actual_network.H`). Same EOS (helmholtz), same
integrator (VODE), Strang, neutrinos off as milestone 1 (SKILL.md Phase 0).

For the first port, Phase 3 translation and host validation took about 75
minutes (parallel agents), followed by several rounds of Loom compile fixes. Most of that work was the network-independent modules, now in
`modules/`; what remains per network is mainly the rate functions.

## What is reused and what is per network

| Piece | Reuse | Source |
|---|---|---|
| `support.H`, `f64_math*.H`, `dual.H` | as is | `kernel/` |
| `eos.H` (helmholtz) | as is (network data from `network_data.H`) | `modules/` |
| `screen.H` (screen5, plasma state) | as is | `modules/` |
| `network_rhs.H` (rhs.H machinery: rhs, jac, construct_rate, LU) | as is | `modules/` |
| `vode.H`, `burner.H` | as is | `modules/` |
| `network_data.H`, `burn.H` | **generated** | `scripts/gen_network_data.py` |
| `rates.H` (tfactors + this network's rate functions, `_r` and `_d`) | **hand-translated** | `prompts/rates.md` |
| `network.H` (NetParams, rate dispatch, postprocess_rate, ener_gener_rate) | **hand-translated** | `prompts/network.md` |
| `burn_kernel.cpp` (kernel root `mpl::burn_zones`) | as is | `root/` |
| `loom_eos_params.H` (host `fill_eos_params`, `fill_eos_tables`) | as is | `harness/` |
| `loom_params.H` (the other host fill_* functions) | template; fill the `PORT:` lines | `harness/` |
| validators | network-neutral templates for dual (`harness/validate_dual.H`), EOS, rates, network, integrator and kernel; the only per-network part is `validators/validate_port.H` (parameter variants, rate branch temperatures, direct rate-function checks: fill the `PORT:` lines). `validate_network.H` also holds the `burn_t` layout `static_assert`s that `validate_integrator.H` relies on | `validators/`, `harness/` |
| test_react wiring | `harness/test_react.patch` (GNUmakefile, Make.package, _parameters, main.cpp), plus `loom_launch.H`, `loom_validate.H/.cpp`, `loom_eos_params.H`, `loom_case.H` and `standalone/loom_case_format.H` as is | `harness/`, `standalone/` |
| standalone HIP reproducer | `main.hip`, `Makefile`, `README.md` as is; bundled from the port by `scripts/make_standalone.sh` (no per-network part) | `standalone/` |
| smoke kernels | templates | `smoke/` |

The generic modules include the per-network files by name (`network_data.H`,
`burn.H`, `rates.H`, `network.H`). The hooks `network.H` must define before
including `network_rhs.H` are listed in `network_rhs.H`'s header comment.
`rates.H` (namespace `mpl`, including `support.H`, `dual.H` and
`network_data.H`, which supplies `Rates::` and `NumRates`) must declare
exactly these; `screen.H`, `network.H`, the smoke kernels and the validators
call them:

```cpp
struct RatesParams {                         // scalars only (kernel ABI)
    int  screening_enable_debye_huckel_skip;     // screening_rp::... (bool -> int); screen.H reads it
    Real screening_debye_huckel_skip_threshold;  // screening_rp::...;                screen.H reads it
    // + one field per runtime parameter the rates read, named <namespace>_<name>
};
template <typename T> struct rate_t { T fr; T rr; };          // interfaces/rhs_type.H
template <typename T> struct tf_t { T temp, t9, ..., lnt9; };  // interfaces/tfactors.H, same fields and order
MPL_DEVICE tf_t<Real> get_tfactors_r(Real temp);
MPL_DEVICE tf_t<Dual> get_tfactors_d(Dual temp);
// one pair per stock rate function reached, e.g. rate_c12ag(tf, den, fr, rr):
MPL_DEVICE rate_t<Real> rate_<name>_r(const tf_t<Real>* tf, const Real den[, const RatesParams* params]);
MPL_DEVICE rate_t<Dual> rate_<name>_d(const tf_t<Dual>* tf, const Real den[, const RatesParams* params]);
// the network's RHS::evaluate_analytical_rate<rate>, den = 1.0 inside as in the stock call:
template <int rate> MPL_DEVICE rate_t<Real> evaluate_analytical_rate_r(const tf_t<Real>* tf, const RatesParams* params);
template <int rate> MPL_DEVICE rate_t<Dual> evaluate_analytical_rate_d(const tf_t<Dual>* tf, const RatesParams* params);
```

`params` is passed to a rate function only if it reads one;
`evaluate_analytical_rate_{r,d}` always take it.

## Steps

1. **Configure and inventory** (SKILL.md Phases 0 and 2):
   `uv run --no-project <skill>/scripts/inventory.py --out <notes>/inventory -- NETWORK_DIR=<net> EOS_DIR=helmholtz INTEGRATOR_DIR=VODE USE_NEUTRINOS=FALSE XTRA_CXXFLAGS=-ffp-contract=off`
   then `scripts/specialize.py --inv <notes>/inventory`. Check the "CPU vs HIP
   defines" section and the per-function table for constructs the modules
   don't handle yet (other screen methods, NSE, aux species).
2. **Kernel directory**: `unit_test/test_react/loom/<net>_vode/`. Copy
   `kernel/*`, `modules/*`, `root/*`, `validators/*` and
   `harness/validate_dual.H` (the validators include each other and the
   translated headers from this directory), and `smoke/*` into its
   `tests/import/`, and write the port's `README.md` (configuration,
   translation rules from `prompts/_common.md`, file ownership, build and
   validation commands). Wire test_react if it isn't wired yet: from the repo
   root `git apply <skill>/assets/harness/test_react.patch`,
   and copy `harness/loom_launch.H`, `loom_validate.H`, `loom_validate.cpp`,
   `loom_eos_params.H`, `loom_case.H` and `standalone/loom_case_format.H`
   into `unit_test/test_react/`. Every validate or
   Loom build then passes `LOOM_KERNEL_DIR=loom/<net>_vode` to make; the
   harness headers reach the kernel's `eos.H`, `burner.H` and
   `burn_kernel.cpp` through it.
3. **Generate the network data**:
   `uv run --no-project <skill>/scripts/gen_network_data.py --out <work> --emit-dir <kernel dir> --support-dir <kernel dir> -- COMP=llvm CXX=clang++ CC=clang NETWORK_DIR=<net> EOS_DIR=helmholtz INTEGRATOR_DIR=VODE USE_NEUTRINOS=FALSE`
   It static_asserts every value against the stock headers before writing.
   **Validate the EOS and duals now**, before touching rates or network:
   `USE_LOOM_VALIDATE=TRUE LOOM_KERNEL_DIR=loom/<net>_vode
   XTRA_CPPFLAGS="-DLOOM_VALIDATE_DUAL -DLOOM_VALIDATE_EOS"` (command in
   `references/test-react-integration.md`, "Native validation"). They need
   only the copied modules, `network_data.H` and `loom_eos_params.H`; no
   `rates.H`, `network.H` or `loom_params.H` (measured on aprox13: 600.6M EOS
   checks, 0 mismatches, in that state). Zero mismatches before step 4.
4. **Translate the rates** (the declarations above) with `prompts/_common.md` + `prompts/rates.md`
   (one sub-agent, or yourself). If the inventory shows autodiff shapes `dual.H` doesn't
   document, run `prompts/dual.md` first.
5. **Translate `network.H`** with `prompts/_common.md` + `prompts/network.md`.
   Then write `unit_test/test_react/loom_params.H` from
   `harness/loom_params.H`: fill `fill_rates_params` / `fill_net_params` with
   one line per `RatesParams` / `NetParams` field, and add a
   `loom_require_supported_settings()` abort for every runtime setting the
   translation omits (`network.use_tables` is already rejected). The
   validators and the launcher all include it; never redefine a `fill_*`
   function in a validator. Lint every kernel file except `validate_*`
   (SKILL.md Phase 3): `rates.H` and `network.H` must be clean; the 18
   `warning[no-reference]` hits from the shipped `eos.H` (13) and `vode.H`
   (5) are expected (scalar references the probes showed import).
6. **Validate everything**. Fill the `PORT:` lines of `validate_port.H`
   first (the network's runtime switches as parameter variants, every rate
   branch temperature, direct calls of rate functions that take a density).
   The rates, network, integrator and kernel validators include
   `loom_params.H` (and through it `burner.H`), so they build only once
   `rates.H`, `network.H` and `loom_params.H` exist. AMReX's make does not
   track flag changes, and deleting `loom_validate.o` was not enough in the
   E2E run: give each validator set its own build tree and executable,
   `TMP_BUILD_DIR=tmp_val_<set> EBASE=loomval_<set>` (e.g. `<set>` = `eos`,
   `net`, `burn`). One combined build with
   `USE_LOOM_VALIDATE=TRUE LOOM_KERNEL_DIR=loom/<net>_vode` and
   `XTRA_CPPFLAGS="-DLOOM_VALIDATE_DUAL -DLOOM_VALIDATE_EOS -DLOOM_VALIDATE_RATES -DLOOM_VALIDATE_NETWORK -DLOOM_VALIDATE_INTEGRATOR -DLOOM_VALIDATE_KERNEL"`;
   adapt `validators/validate_integrator.H`'s zone setup only if the network
   needs it. Zero mismatches.
7. **Smoke-compile on the GPU host after each module validates**, not at the
   end: `scripts/loom_smoke.sh launch <name> <kernel dir> tests/import/<smoke>.cpp <roots>`
   with templates from `smoke/` (rates, network, full burn), one at a time:
   `loom_smoke.sh wait <name>` blocks until a run finishes. A fused full
   burn compiles in 9 min (triple_alpha, Rosetta) to ~30 min (aprox13),
   nearly all of it in register repair (see
   `references/footguns.md`); use Loom tools built with the allocator fix
   (`scripts/build_loom_tools.sh` with `patches/alloc-relocation-overlap.patch`).
   SUBRANGE/026 errors: see `MPL_ASSUME` in `references/loom-frontend-rules.md`.
8. **Kernel root, launcher, GPU comparison** (SKILL.md Phases 5-6).

## Not covered yet

- Neutrinos (sneut5: sin/cos OCML recipes and `Dual3` paths are not wired).
- NSE (table or net), SDC burners, screening methods other than screen5.
- pynucastro-generated networks (most of `networks/`): they don't use the
  hand-written rate library in `rates/`; their generated rate code needs its own approach (ideally a
  Loom-compatible pynucastro output mode).
- Integrators other than VODE (`prompts/integrator.md`), EOSes other than
  helmholtz (`prompts/eos.md`).
