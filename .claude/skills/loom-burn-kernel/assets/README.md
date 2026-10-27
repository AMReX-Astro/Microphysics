# Reusable Loom kernel library

Validated pieces every Microphysics → Loom port needs. Copy `kernel/` into
the port's kernel directory (`unit_test/test_react/loom/<network>_<integrator>/`)
and `harness/loom_validate.*` into `unit_test/test_react/`, then adapt.

Snapshot taken 2026-10-18 from the aprox13 + helmholtz + VODE port (branch
`loom-aprox13-vode`, Microphysics `60562512`, vendored autodiff last changed in
`f3ca3777`). When a port improves one of these files, copy it back here and
update this README.

**Start here for a new network: `NEW_NETWORK.md`.**

| Directory | Contents |
|---|---|
| `kernel/` | platform layer, OCML f64 math, POD duals (below) |
| `modules/` | network-independent translated modules (below) |
| `root/` | `burn_kernel.cpp`: the kernel root `mpl::burn_zones` (one lane per zone, reproduces test_react's `do_react`) |
| `validators/` | network-independent validator templates (EOS, rates + screening, network RHS/Jacobian/LU, whole-burn integrator, kernel root vs `do_react`) plus `validate_port.H`, the per-port hooks |
| `harness/` | test_react validation hook, dual validator, `loom_launch.H` (HIP launcher for test_react), `loom_params.H` template |
| `smoke/` | Loom smoke-compile kernels for `scripts/loom_smoke.sh` |
| `prompts/` | agent prompt templates: `_common.md` (prepend to every prompt) + per-task prompts |

## modules/

Network-independent for any `rhs.H`-framework network with helmholtz + VODE,
Strang, screen5, neutrinos off, HIP code path. They include the per-network
files by name (`network_data.H`, `burn.H`, `rates.H`, `network.H`); see
`NEW_NETWORK.md` for the contract.

| File | What | Host (bit-exact vs stock) |
|---|---|---|
| `eos.H` | helmholtz EOS: composition, eos(), reset, override, prad limiter | 600.6M checks, 0 mismatches |
| `screen.H` | screen5, plasma state, generated screen factors | in rates validator: 11.7M checks, 0 mismatches |
| `network_rhs.H` | rhs.H: rhs, jac, construct_rate, intermediate rates, LU (strict `<` bounds, no `MPL_ASSUME`) | 31.1M checks, 0 mismatches |
| `vode.H` | VODE + Strang integrator glue (`aion`/`aion_inv` from `network_data.H`; stock-shaped loops) | whole burns: 11.6M checks, 0 mismatches |
| `root/burn_kernel.cpp` | kernel root | vs stock do_react: 1.86M checks, 0 mismatches (16^3, 32^3, sub-box, prad limiter on) |
| `burner.H` | burner, retry, `BurnParams` | (whole-burn validator) |

At the pin `5e6af0ec1a` (gfx90a, 2026-10-27), aprox13 built from these
modules is bit-exact against the stock HIP reference on the test_react
4096-zone case (`n_rhs`, `n_step` and success included; standalone driver,
case file from the stock build), with 256 VGPR and 70,360 B scratch
(compile 761 s); triple_alpha_plus_cago was not rerun there. The modules
were unchanged by the pin move. VODE's species-loop `MPL_ASSUME`s are still
required there: without them the burn fails SUBRANGE/026.

Status at the previous pin (`fa8eba4554`, gfx90a), for aprox13 and triple_alpha_plus_cago:
- host, all six validators: 0 mismatches (aprox13 948M checks,
  triple_alpha_plus_cago 411M); rerun them after any module change;
- Loom production root (`burn_zones`, 32 workgroups, 4096 zones): bit-exact
  against the stock HIP test_react reference, `n_rhs`, `n_step` and success
  included; 256 VGPR, scratch aprox13 72,232 B, triple_alpha_plus_cago
  42,520 B. A hipcc build of the same kernel tree (OCML or recipe math) is
  bit-exact too. Every zone succeeds in these runs, so the failure exits are
  checked by translation only;
- Loom is about 7.4x slower than hipcc on the same source (256 VGPR and
  42-72 KB of spills vs hipcc's 408-492 registers including AGPRs, 2.4-6 KB);
- two rewrites remained: `MPL_ASSUME` on VODE's species loops with `break`
  (SUBRANGE/026 in the full burn, although probe f17 passes; still required
  at `5e6af0ec1a`) and `actual_log_screen`'s `done` flag (probe f24; no
  longer required at `5e6af0ec1a`).

Update this section when a port moves a module forward.

## Launcher (harness/loom_launch.H) and the per-network loom_params.H

`loom_launch.H` loads the Loom code object with `hipModuleLoad`, fills
`BurnParams` on the host, passes the helmholtz tables by their
`AMREX_GPU_MANAGED` addresses, launches `burn_zones` (one box, zone count
checked against `unit_test.loom_ncells`), and in compare mode runs stock and
Loom burns alternately and aborts unless every zone is bit-identical. It is
compiled with `USE_HIP=TRUE USE_LOOM_BURN=TRUE LOOM_KERNEL_DIR=loom/<port>`;
the patched GNUmakefile (`harness/test_react.patch`) turns `LOOM_KERNEL_DIR`
into `-DLOOM_KERNEL_DIR_PATH` (unquoted; the headers stringize it). It aborts if the device
has no managed memory (the tables are passed by address). Its comparison treats any NaN as equal to
any NaN (`loom_same_bits`; `references/test-react-integration.md`,
"Comparing outputs").

It includes the port's own `unit_test/test_react/loom_params.H`, made from
the `harness/loom_params.H` template (the port fills in the per-network
parts): host functions
`fill_rates_params`, `fill_net_params`,
`fill_integrator_params`, `fill_burn_params` reading the stock runtime
parameters into the kernel's params structs (plus `fill_eos_params` and
`fill_eos_tables` from the network-independent `loom_eos_params.H`, which it
includes), and `loom_require_supported_settings()` aborting on unported
settings (`network.use_tables`). The burn validators include it too; the EOS
validator includes only `loom_eos_params.H`.

`scripts/loom_gpu_compare.sh` runs Phase 6 on the GPU host: production
compile of the root, argument-metadata dump, optional managed-memory probe,
and test_react compare runs: one run with the inputs file by default; add
the prad-limiter run with
`RUNS=";eos.prad_limiter_rho_c=1.0e7 eos.prad_limiter_delta_rho=3.0e6"`.

## validators/, smoke/, prompts/

- `validators/validate_{eos,integrator}.H`: network-independent templates:
  stock-state construction via the stock EOS, the test_react zone grid as
  `do_react` builds it, whole-burn `burn_t` comparison, parameter toggle and
  restore, sensitivity checks, retry/failure variants. Add the network's own
  hard zones (by its species indices) to the integrator template.
- `validators/validate_{rates,network}.H`: network-neutral templates
  (snapshot of the triple_alpha_plus_cago port's validators, 2026-10-18,
  made generic): rates loop over `Rates::NumRates` and the screened pairs,
  network over compositions built from `NumSpec`; `validate_network.H` also
  holds the `burn_t` layout `static_assert`s. Their per-network inputs are
  the `PORT:` hooks in `validators/validate_port.H` (runtime switches,
  branch temperatures, direct rate calls). A port may add checks of its own.
- `validators/validate_kernel.H`: its `dt` + 1 ulp self-check runs on the
  full 16^3 grid (a cold corner doesn't burn for gentle networks).
- `smoke/smoke_{wave1,wave2,wave3,burn}.cpp`: EOS + math, rates + screening,
  network RHS/Jacobian/LU, full burn. Paths assume they live in
  `<kernel dir>/tests/import/`.
- `prompts/`: `_common.md` holds every rule and mid-flight correction from the
  aprox13 port; `rates.md` and `network.md` are the per-network prompts;
  `eos.md`, `integrator.md`, `dual.md` are for a new EOS, integrator, or
  autodiff shapes.

## kernel/

| File | What | Validation status |
|---|---|---|
| `support.H` | Platform layer: `MPL_DEVICE`/`MPL_KERNEL`, `Real`, `mpl::math` bindings (host = the std call stock code makes, Loom = Loom op or OCML recipe), `powi` with AMReX's grouping, `min/max/clamp` with std semantics, `Array1D`/`Array2D` (1-based accessors), `MPL_ASSUME` | Bindings checked through the module validators |
| `f64_math.H` | OCML f64 `exp`, `log`, `cbrt` (from hrx-system `experimental/loom-codegen-comparison`, `7cf21edd85`) | Used by the chemistry campaigns on gfx90a/gfx942; bit-identical to hipcc's OCML on gfx90a |
| `f64_math_extra.H` | OCML f64 `pow`, `log10`, `tanh` (llvm-project `amd/device-libs` `d366fa84f3fd`, same revision) | Host: max 1 ulp vs MPFR in every sweep, all C Annex F special cases match libm. GPU (gfx90a): pow/log10/tanh and `f64_math.H`'s exp/log/cbrt bit-identical to hipcc's OCML over ~1M samples per sweep; only 3 cbrt NaN-payload differences (irrelevant under the NaN rule) |
| `dual.H` | POD `Dual` (= `autodiff::dual`) and `Dual3` (= `dual_array<1,3>`) with free functions that reproduce autodiff's expression-template evaluation bit for bit under `AUTODIFF_STRICT_ASSOCIATIVITY`; its header comment is the translation guide (expression shape → call sequence, ~40 worked examples) | 24.9M checks, 0 mismatches against autodiff (`harness/validate_dual.H`) |

Notes:
- `MPL_ASSUME(cond)` is `__builtin_assume(cond)` under Loom and a no-op on
  the host, so it never changes validated host results. Put
  `MPL_ASSUME(lo <= i && i <= hi);` at the top of each loop body whose index
  into a private array has runtime bounds (VODE order loops); without it,
  compile fails with SUBRANGE/026. Counted loops with compile-time limits
  (the network LU) need none when written with strict `<` bounds. Rules for the
  condition: `references/loom-frontend-rules.md`, SUBRANGE row.
- The Loom side has no f64 exp/log/pow/log10/tanh/sin/cos/cbrt (probes); these
  recipes are the only path. sin/cos are not transcribed yet (only needed
  with neutrinos on).
- OCML `pow` has no `y == 1` special case: `pow(3, 1)` is 1 ulp below 3. The
  GPU check confirmed hipcc's OCML does the same (`pow(x, 1) != x` for the
  same 6.21% of samples on both sides), so it is correct for bit comparison
  against the HIP reference, but it differs from libm.
- `f64_math_extra.H` replaces OCML's `MATH_FAST_RCP` (hardware reciprocal +
  Newton) with IEEE division, as `f64_math.H` does; the GPU check found no
  difference on gfx90a.
- `dual.H` was derived for the reachable aprox13 path. Expression shapes it
  doesn't document (e.g. number^dual in use, `NumberDualDualMulOp`) must be
  added and validated before a translation uses them. If the vendored autodiff
  changes, rerun `validate_dual.H`.
- Never write derivative arithmetic by hand: compound assignment and mixed
  number/dual forms have non-obvious autodiff evaluation orders (see the
  "Warnings for translators" part of the guide).

## harness/

- `loom_validate.H`, `loom_validate.cpp`: the native validation hook for
  test_react (bit comparison helpers; dispatcher for `validate_<module>.H`
  selected by `-DLOOM_VALIDATE_<MODULE>`). Copy `loom_validate.H`,
  `loom_validate.cpp`, `loom_launch.H`, `loom_eos_params.H` and
  `loom_params.H` into `unit_test/test_react/`.
- `loom_eos_params.H` (copy as is): `fill_eos_params` / `fill_eos_tables`,
  network-independent; the only host header the EOS validator needs, so the
  EOS validates before any rates/network translation. `loom_params.H`
  includes it, and it defines `LOOM_KERNEL_PATH(file)` from the unquoted
  `-DLOOM_KERNEL_DIR_PATH` the GNUmakefile sets.
- `test_react.patch`: the whole test_react wiring against stock
  `development` (`git apply <skill>/assets/harness/test_react.patch`
  from the repo root): GNUmakefile `USE_LOOM_VALIDATE` / `USE_LOOM_BURN`
  blocks (`LOOM_KERNEL_DIR` required; the kernel directory goes **last** on the
  include path, after `Make.unit_test`, so stock `eos.H`/`burner.H` keep
  winning), `Make.package` entries, the `unit_test.loom_*` runtime parameters,
  and the two `main.cpp` hooks (validation after `network_init()`, the
  `loom_burn_main` branch around the
  ParReduce). Network-neutral.
- `validate_dual.H` (copy into the kernel directory): the dual validator, primitive tests only (every `dual_*`
  / `dual3_*` function vs autodiff on edge and random values). A port adds
  its own section comparing whole stock rate, screening and rhs.H expressions
  against the guide's call sequences.

The recipe accuracy tests (host MPFR sweep, GPU bit check against hipcc's
OCML) and the standalone HIP reproducer with its case-file writer are in
the `loom-burn-kernel-tools` follow-up.
