# Translating the Microphysics burn path

Microphysics device code is template-heavy C++ written for AMReX. The Loom
importer accepts a C-like subset (see `loom-frontend-rules.md`). The
translation is therefore a *specialization*: one network, one EOS, one
integrator, one set of `-D` flags, written out as plain functions over plain
arrays. Keep that in mind when choosing what to generate and what to write
by hand.

## Inventory limits (Phase 2)

`scripts/inventory.py` generates the configuration's headers into its own
`TMP_BUILD_DIR` (the test's `tmp_build_dir` may hold another network's
headers, so it is neither touched nor trusted). Cover its limits by reading
code:
- It analyzes the CPU build's `-D` set, but also takes the `USE_HIP` build's
  compile line (`--hip-arch`, default gfx90a; a stub `hipconfig`/`hipcc`
  stands in when ROCm isn't installed, since `make -n` never runs them) and
  lists the differing defines in "CPU vs HIP defines". The kernel follows the
  HIP side (Phase 4 makes the CPU reference match it). Also note
  `AMREX_USE_GPU` branches and macros like `HIP_CONSTEXPR` that change
  meaning between the two.
- Reachability is by unqualified name, so it is an upper bound; host-only
  code (printing, `#ifndef AMREX_USE_GPU`) shows up as reachable on CPU.
- Code reachable only through a runtime flag (`use_tables`, `jacobian`) is
  reachable. Decide with the user whether such branches are ported or
  rejected at launch ("Runtime-gated paths" below).
- It does not measure ternary nesting depth or base classes.
- Keep the matcher list in step with Phase 1 results.

## General method

- **Specialize at the source level, not the semantics level.** Resolve
  templates and `if constexpr` branches for the pinned configuration
  (`neqs`, `NumSpec`, `NumAux`, `integrator_rp::jacobian`, `SCREEN_METHOD`,
  `NSE`, `STRANG`). Keep compile-time constants compile-time (`constexpr int`
  trip counts); keep runtime parameters runtime.
- **Statement for statement, same order.** Keep each expression's
  evaluation order and the `+`/`*` grouping the source already has (Loom
  imports with `--fp-model=precise`; contraction policy in SKILL.md). Changing
  it is a non-bit-exact change (user approval, documented).
- **Keep one-to-one traceability.** Each translated function carries a
  comment with its source path and function name. If you outline a helper
  (nested ternaries, early-return loops), name it after its origin.
- **`amrex::Real` → `Real`.** The kernel cannot include AMReX, so define
  `using Real = double;` in the kernel namespace and `static_assert`
  (host side) that `amrex::Real` is `double` (`PRECISION=DOUBLE`).
- **1-based arrays.** `Array1D<T,1,N>` and `Array2D<T,1,N,1,N>` become
  0-based flat arrays; translate index expressions mechanically
  (`a(i)` → `a[i-1]`, `m(i,j)` → `m[(i-1)*N + (j-1)]`). Do this in one place
  (a generator or a macro) because hand-converted indices are the easiest
  bug to introduce and the hardest to see. Note Array2D's storage order and
  keep row-major *consistently* in the translation.

## burn_t and the integrator state

- Keep `burn_t` as a struct with the same fields (drop only fields that are
  `#ifdef`'d out for the configuration). Pass `burn_t*` everywhere a
  `burn_t&` was.
- Integrator types (`dvode_t` in `integration/VODE/vode_type.H`,
  the Rosenbrock type in `integration/Rosenbrock/rosenbrock_type.H`) hold
  `Array1D/Array2D` members and are used through member functions. Turn them
  into structs of flat arrays plus free functions taking a pointer. The
  chemistry kernel's `ScratchView` (a struct of pointers into one scratch
  record) is the pattern that imported cleanly.
- Control flow in the integrators (VODE `dvstep`, `dvnlsd`, `dvjac`; the
  Rosenbrock step/reject loop) uses `return` inside loops and multi-exit
  loops. Rewrite with explicit `done`/`retry` flags and a result code,
  preserving every exit point and its counter updates. Add `MPL_ASSUME`
  bounds on stage/retry counters and order-dependent indices
  (`loom-frontend-rules.md`, SUBRANGE row).
- Keep the counters (`n_rhs`, `n_jac`, `n_step`) exactly as the source
  increments them; test_react writes `n_rhs` and `n_step` out and you will
  compare them.

## Runtime parameters

Every parameter the burn reads (`integrator_rp::rtol_spec`, `atol_*`,
`ode_max_steps`, `jacobian`, `use_jacobian_caching`, `screening_rp::*`,
`network_rp::*`, `eos_rp::*`, `small_temp`, ...) is an
`AMREX_GPU_MANAGED` global in Microphysics. Loom cannot read globals, so:
- Generate a POD `BurnParams` with one field per parameter actually read on
  the burn path (Phase 2 inventory), with the same types.
- The host fills it from the `*_rp::` namespaces after `init_unit_test()`,
  copies it to device memory, and passes a `const BurnParams*`.
- In the kernel, read each field once into a local `const` near the top of the
  root and pass values or the pointer down. Parameters that are constant for
  the run but read in hot loops are worth hoisting (campaign lesson:
  tolerances and controller settings as constants were measurable wins);
  never replace them with literal values, which would specialize on the
  inputs file.

### Runtime-gated paths that are not ported

Some parameters switch on whole alternative code paths that the configured
inputs never use. Porting them adds kernel inputs, tables and untested code.
For these, the first milestone rejects the non-default value at launch: the
host launcher checks the parameter and calls `amrex::Abort` with a message
naming it, the parameter is left out of `BurnParams`, and the kernel keeps
only the default branch. Record each one in the agent note.

| Parameter | Default | Not ported | Decision |
|---|---|---|---|
| `network.use_tables` | 0 | tabulated rates: `ttab`, `rattab` (~9 MB of duals for aprox13), `rate_tab_t::initialize`, `evaluate_tabulated_rate` | user, 2026-10-18 |

Ported although the inputs leave them off: `eos.prad_limiter_rho_c` /
`prad_limiter_delta_rho` (default -1; brings in `tanh`; user, 2026-10-18).
Phase 6 includes a run with the limiter on.

Candidates still to raise with the user when reached: `integrator.jacobian`
= 2 (numerical Jacobian; default 1), `screening.enable_debye_huckel_skip`.

Ask the user before adding a row; a parameter whose non-default value the
inputs file sets, or that a production code relies on, has to be ported.

## Tables and constant data

- **Device tables** (`AMREX_GPU_MANAGED` arrays: helmholtz `f`, `dpdf`, `ef`,
  `xf`, `d`, `t`, `*_sav` in `EOS/helmholtz/actual_eos_data.H`; tabulated
  rates; NSE tables; `rhs.H` rate data) become pointer arguments, one per
  table: a pointer to a struct of pointers does not import (probes c03/c06).
  The kernel root rebuilds the tables struct by value from them (c04), as
  `assets/root/burn_kernel.cpp` does for the 14 helmholtz tables
  (`EosTables`). On moth (MI210, XNACK off) the host passes the
  `AMREX_GPU_MANAGED` table addresses directly; they are device-readable.
  Index arithmetic for `f[jmax][imax][9]`
  must match the C array layout exactly; validate one lookup natively.
- **`constexpr` arrays** (`aion`, `zion`, `bion`, rate coefficient tables,
  partition-function grids) fail at namespace or `static constexpr` member
  scope. Use a `constexpr` function with an if-chain (folds when the index
  is constant), or a local `const Real t[N] = {...}` inside the function that
  uses it. Scalar `constexpr` constants are fine.

## Generated network data (`rhs.H` networks)

For an `rhs.H`-framework network (aprox13/19/21, iso7, ignition_chamulak,
ignition_simple, powerlaw, rprox, triple_alpha_plus_cago) don't hand-write the
network data: `scripts/gen_network_data.py` generates `network_data.H`
(Species/Rates indices, every `rhs_data` field, intermediate-rate indices,
zion/aion/bion/mion, screen factors) and `burn.H` from the stock headers;
`assets/modules/network_rhs.H` is the generic rhs.H machinery. What stays per
network (`rates.H`, `network.H`) is in `assets/NEW_NETWORK.md`.

```bash
uv run --no-project <skill>/scripts/gen_network_data.py \
  --out <work dir> --emit-dir <kernel dir> --support-dir <kernel dir> \
  [--flags <inventory>/flags.txt | -- NETWORK_DIR=<net> USE_NEUTRINOS=FALSE ...]
```

It static_asserts every generated value against the stock constexpr source
before writing anything, and fails (or warns) on configurations the kernel
does not port (NEUTRINOS, NSE, NONAKA_PLOT, non-screen5 screening).

## Autodiff (screening, rates, he-burn networks)

`autodiff::dual` and `admath::` rely on operator overloading, expression
templates and base classes, none of which import. Replace them with POD dual
structs and free functions, passed and returned by value (probes o04/o05).
The skill ships this as `assets/kernel/dual.H` (validated bit for bit against
autodiff; its header is the translation guide): use it rather than writing a
new one, and extend it, with validation, for shapes it doesn't cover:
`struct Dual { double v; double d; };` for `autodiff::dual`, and a fixed-size
gradient array for `dual_array<1, N>` / `GradArray` (screening's composition
derivatives).

The free functions must reproduce the exact operation sequence of
`util/autodiff/autodiff/forward/dual/dual.hpp` as Microphysics builds it,
which is with `AUTODIFF_STRICT_ASSOCIATIVITY` defined
(`util/microphysics_autodiff.H:18`). For example `a * b` is evaluated as
`self = a; self *= b`, i.e. `grad = a.grad * b.val; grad += a.val * b.grad;
val = a.val * b.val`; nested expressions go through a temporary
(`assign(tmp, expr); self *= tmp`), and number * dual scales the gradient.
Read the `assign*` functions for every operation you use (add, sub, mul, div,
inv, pow, exp, log, sqrt, ...) and copy their formulas and order; don't derive
the calculus yourself, or the result won't be bit-exact. Write one host test
per function comparing the POD version against `autodiff::dual` bit for bit.

This is the riskiest part of the translation:
- validate each expanded function natively against the dual version at many
  points (including both screening regimes and the blend region);
- write the derivative accumulators as SSA `const` chains, not `+=` on
  mutable scalars (both a correctness-review aid and, in the campaigns, a
  12% performance difference);
- if a whole function's derivative is forward-differentiated in one big
  sweep, watch the import time (8+ min once); split it if it grows.

## Math

- Bind `sqrt`, `fabs`, `fma`, `floor`, `trunc`, `rint` to `loom::scalar`
  operations (they lower for f64).
- `exp`, `log`, `cbrt`, `pow`, `log10`, `tanh`: use the skill's
  `assets/kernel/f64_math.H` and `f64_math_extra.H` (OCML transcriptions,
  bound in `assets/kernel/support.H`).
- Anything else the burn uses (`sin`/`cos` for sneut5, `expm1`, `log1p`,
  `atan`, `erf`, ...): transcribe the matching OCML
  routine (ROCm-Device-Libs `ocml/src/*D*.cl`) statement for statement in the
  same style, with explicit `fma`, one operation per statement, and the same
  special-case handling. Do **not** substitute identities like
  `pow(x,y) = exp(y*log(x))` or `log10(x) = log(x)/ln10`: they are not the
  functions the HIP reference calls, and the error grows with `|y log x|`.
- `amrex::Math::powi<N>` / integer powers: keep the template recursion and its
  multiplication grouping (`support.H` `powi`).
- `std::min/max`: by-value helpers returning the *first* operand on ties and
  unordered comparisons (not minnum/maxnum).
- Validate every recipe against MPFR/host libm over a sweep that covers the
  argument ranges the burn actually sees (record max ulp error), and bit for
  bit against hipcc's OCML on the GPU (`assets/f64-tests/`, `ocml_gpu/`). The
  shipped exp/log/cbrt/pow/log10/tanh recipes are bit-identical to OCML on
  gfx90a (2026-10-18).
- **Calls with a constant argument** (e.g. `rate_c12ag`'s deboer17 branch:
  `exp(24.1)`, `exp(7.4)`). hipcc may fold them at compile time, correctly
  rounded, while the Loom kernel evaluates the OCML recipe at run time; the
  two can differ in the last bit, which breaks bit-exactness only on that
  path. In Phase 6, list every recipe call whose argument is a compile-time
  constant, check the stock HIP code object's disassembly for each (a call
  or inlined OCML sequence vs a literal), and include a GPU compare run with
  the switch that reaches it turned on. If hipcc folds, the bit-exact fix is
  to use the folded constant on the Loom side too: a non-statement-for-statement
  rewrite that needs the user's approval and a note in the hand-off.

## NSE

If the configuration enables NSE (`NSE_TABLE` or `NSE_NET`), `burner()`
branches on `in_nse()`. Port the check and the NSE update too, or ask the
user to pick a configuration without NSE for the first port. Don't silently
drop the branch: test_react's grid reaches NSE temperatures.

## Things that must not leak into the kernel

- `std::cout`/printing, `amrex::Abort`, `AMREX_ASSERT` (host-only; keep their
  conditions as error codes if they change control flow),
- `#ifndef AMREX_USE_GPU` host-only branches (the Loom path follows the GPU
  branch, since the HIP reference does),
- NONAKA_PLOT or other diagnostics.
