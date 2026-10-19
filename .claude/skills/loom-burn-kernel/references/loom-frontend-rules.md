# Loom C++ frontend: construct → rewrite

Status as of hrx-system `7cf21edd85` (react-kernels, 2026-10-18; compiler
sources identical to public main `4031286f0b`), measured
with this skill's `probes/` kit (results and probe names in
`probes/README.md`). **Re-probe at your commit** (SKILL.md Phase 1):
rows marked "changed recently" or anything older than the commit you use
may be stale.

## Import invocation

Command lines for import and compile are in `test-react-integration.md`
("Remote import/compile commands"). The importer defaults to
`--std=c++26`, while Microphysics is C++20; a few rules differ between them (e.g. enum/floating arithmetic became ill-formed).
The probes and the campaigns used the default.

Only code reachable from `--root` is lowered, so host-only code in the same
TU only has to parse. Namespace-scope global *definitions* are an exception:
they are rejected even when unreachable (LOWERING/059).

Under `__loom__`, `#include <hip/hip_runtime.h>` gives `__global__`,
`__device__`, `blockIdx`/`threadIdx`/`blockDim`, but no math and no atomics.
`#include <loomcxx/atomic.h>` gives atomics. There are no C++ std headers.

## Fails → rewrite

| Construct | Status | Rewrite |
|---|---|---|
| `<cmath>`, `<limits>`, `<utility>`, any std header | fails (PARSE/036) | Own helpers; `__DBL_MAX__` etc.; hand-written numeric limits |
| `std::exp/log/pow`, `__builtin_*` math | fails | `support.H` `mpl::math` bindings → `loom::scalar::*` or the OCML recipes |
| f64 `exp, log, pow, tanh, sin, cos` (LOWERING/036, #1989), `cbrt, log10` (TARGET/001); also `expm1, log1p, atan2, exp2` | fail (m_f64_*) | OCML-transcribed recipes: `f64_math.H` (exp/log/cbrt), `f64_math_extra.H` (pow/log10/tanh); transcribe the OCML routine for the rest, statement for statement; never `exp(y*log(x))` for pow. A native-lowering branch exists upstream; re-probe |
| f64 `copysign` | works (m_f64_copysignf) | — |
| `int(double)` conversion | works (`v_cvt_i32_f64`) | — |
| f64 rint/floor/trunc | work natively (`roundevenf` → `v_rndne_f64`) | — (`f64_math.H`'s rint emulation and its comment predate this; swapping it is numerically identical but cost 1.2% once through register allocation) |
| f64 ldexp/frexp | no op | bit-level helpers, as in `f64_math.H` |
| Scalar reference parameters (`double&`, `int&`) and scalar member slices (`f(s.T)`) | import, `static` or not (r01, r02, r04, r11, r13; #1993) | Keep them |
| Reference to a whole struct (`burn_t&`, `const Burn&`, `template<class T> f(T&)`) | fails even `static`/`const` (r03, r05: "require an automatic SSA object") | Pointers (`burn_t*`) |
| References to buffer elements (`f(out[i])`), array refs (`const double (&)[N]`), overlapping mutable refs (`f(s, s.T)`), local reference bindings (`double& r = s.x`), reference returns (incl. `amrex::max`'s `const T&`) | fail (r06-r10, r12) | Pointers, by-value returns, explicit copies |
| Member function calls, `operator()`, `operator[]` | fail | Free functions taking a pointer; flat arrays with explicit index math |
| `Array1D/Array2D` (ref-returning accessors, base classes) | fail | Plain `Real a[N]`, `Real m[N*N]` row-major, 0-based |
| Lambdas (`constexpr_for`, `ParallelFor`) | fail | Indexed loops; template-struct recursion for compile-time unrolling |
| `switch` | fails at runtime (folds if `constexpr`) | if/else chain |
| range-for, `for (;;)` | fail | Indexed loops; `while (true)` + `break` |
| `return` in a loop | fails (LOWERING/059) | `done`/result code + `break`, return after the loop |
| Loop with several exits (`break`s/`continue` inside branches) under divergent control | may import and then fail the AMDGPU backend: TARGET/034 "branch lowering constraint 'masked_region_exits_by_cfg' is not satisfied" (VODE's dvode step loop, 2026-10-18) | Single-exit form: `bool done = false; while (!done) { ... }`, each former exit sets `done` after the same side effects, later code guarded by `if (!done)`, a `continue` becomes a guard on the rest of the body (keep its exact predicate, negated with `!(...)` so NaNs behave the same) |
| Uninitialized locals | fail | Initialize every local; `T x{}` for structs |
| Base classes, operator overloading (autodiff `dual`) | fail | POD dual structs with free functions (see translation guide, "Autodiff") |
| Constructors with arguments (`D x(1.0)`) | fail (o07: records must be complete aggregates) | Aggregate init `D x{1.0, 0.0}` |
| Structured bindings (`auto [f, r] = rhs_term(...)`) | fail (o08) | Return a 2-field struct; read `.first`/`.second` |
| `constexpr` local of record type (`constexpr rhs_t data = rhs_data(rate);`) | fails (f02: "unsupported C++ type"), though the `constexpr` function with a `switch` is fine | One `constexpr` scalar per field used: `constexpr int expA = rhs_data(rate).exponent_A;` |
| Enum in arithmetic with `double` | fails to parse (d10, PARSE/036): correct under the default `--std=c++26`; the pinned parser also rejects it in C++20/23 modes (reported fixed upstream as `0d27041f`, a commit not found in public hrx-system history as of 2026-10-18; re-probe d10) | `static_cast<int>(e)`; aprox13 doesn't use it |
| Floating-point contraction | none: `a * b + c` → `v_mul_f64` + `v_add_f64` under `--fp-model=precise` (m_fp_contract) | — |
| Namespace-scope non-constexpr globals (runtime params, tables) | fail | Kernel arguments: a params POD pointer and table pointers |
| Pointer stored in a struct, read through a pointer to that struct (`const Tables* t; t->f[i]`), including a kernel argument pointing at a struct of pointers | fails (c03, c06: "stored pointers and references require an object storage representation") | Pass the tables struct **by value** (c04) or one flat buffer with compile-time offsets (c05); kernel takes the table pointers as separate arguments |
| Namespace-scope `constexpr` *arrays*, `static constexpr` member arrays | fail | Local `const Real t[N] = {...}` inside the function, or an if-chain/`constexpr` function; scalar `constexpr` is fine |
| Ternaries nested > 32 deep | fail (region nesting limit) | Outline into helpers preserving lazy evaluation |
| Stage loops with break/continue | a simple retry loop imports without help (f07); complex stage loops hit SUBRANGE/026 in the campaigns | `MPL_ASSUME(0 <= i && i < N)` on the counter (next row) |
| Private-array index with data-dependent loop bounds (triangular LU loops, VODE order loops), an index read from a struct field (`vs->el[vs->L]`), or a loop with `break` (its compound condition hides the bound, even with constant limits) | imports, then **compile** fails SUBRANGE/026 "access byte end may reach ..." even though probe f06 passed in isolation (network LU, 2026-10-18) | `MPL_ASSUME(lo <= i && i <= hi)` (support.H: `__builtin_assume` under Loom, no-op on host) at the top of each loop body / after each derived index, using the array's true static bounds; loops with compile-time bounds don't need it. Operands must be local scalars and constants: no struct fields (`vs->L` fails: "predicate member paths ..."; copy it into a `const int` local and index with that), no arithmetic on the variable (`j + 1 <= N` fails: "predicate operands must be retained scalar values ..."; move the constant to the other side). Rely on the smoke compile, not the isolated probe |

This table is the authoritative place for `MPL_ASSUME`/SUBRANGE/026; other
files point here.

## Using `f64_math.H` / `f64_math_extra.H`

`f64_math.H` (kept verbatim from the chemistry kernel) expects a `DEVICE`
macro defined before inclusion; `support.H` defines it (`__device__` under
Loom) next to its own `MPL_DEVICE`/`MPL_KERNEL`. Both headers put their
functions in `chemistry::f64_math`. Keep the files and their license headers
as is; bind them from the kernel's own math namespace (`mpl::math`), as
`support.H` does.

## Lint

`scripts/loom-lint/` (ast-grep, syntax-only) mirrors this file; update both
when Phase 1 probes change what imports. Errors must be fixed;
`no-reference` warnings are allowed only for forms the probes showed import.
Read the output, not just the exit code: a rule file with a YAML error makes
ast-grep print `Error: Cannot parse rule` and skip the scan. After adding or
editing a rule, run `scripts/loom-lint/tests/selftest.sh` (every rule must
parse and fire on `tests/should_fire.cpp`; add a line there for a new rule).
The `no-stock-namespace` rule matters most for agents: validators include the
stock headers, so a leftover `autodiff::`, `amrex::`, `RHS::` or `*_rp::` in
kernel code still compiles and validates on the host, and only the lint or the
Loom import catches it. Write native host fallbacks as the `#else` branch of
`#if defined(__loom__)` (or `#ifdef __loom__`); the lint skips those. (The
lint is clean on the accepted chemistry kernel, adv-mi210-12's seed.)

## Works

`if constexpr`; recursive `powi<N>` templates; static member function
templates; `enum class`; anonymous namespaces; structs with arrays by value;
large automatic structs with default member initializers; local 2-D arrays
and pivoted LU; pointers to local structs passed to helpers; struct return by
value; `break`/`continue`; `while (true)`; GCC `__atomic_*` builtins;
`__builtin_bit_cast`; f64 add/sub/mul/div/sqrt/fma.

## Kernel root

`assets/root/burn_kernel.cpp` is the template: `MPL_KERNEL void
mpl::burn_zones(const BurnParams* params, <14 helmholtz table pointers>,
Real* state, <strides, box corner and extents, component indices>, int*
n_rhs, <strides>, Real dt, int* zone_success)`, one lane per zone
(`idx = blockIdx.x * blockDim.x + threadIdx.x`, guarded by `idx < ncells`).
It rebuilds `EosTables` by value from the table pointers. On the host the
same function takes the zone index as a trailing `int`. The Loom-side root
name uses `.` where C++ uses `::` (`mpl.burn_zones`).
