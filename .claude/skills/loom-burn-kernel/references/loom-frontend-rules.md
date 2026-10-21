# Loom C++ frontend: construct → rewrite

Status as of public hrx-system main `fa8eba4554` (2026-10-20; probed
2026-10-21 on gfx90a), measured
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
TU only has to parse. Namespace-scope *mutable* global definitions are an
exception: they are rejected even when unreachable (LOWERING/059, d21).
Unreachable `constexpr` tables are fine (d20).

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
| `const` reference to a whole struct (`const burn_t&`, `const T&` template parameters) | imports (r05), but calling a member's accessor through it does not, even the `const` overload (`vs.y(n)` with `const dvode_t& vs`, r23: "mutation requires an owned automatic source binding"); plain fields and raw array members read fine (r24) | Keep it for structs without accessor members; a struct with `Array1D`/`Array2D` members goes by (const) pointer (`vs->y(n)`, o19) |
| Mutable reference to a whole struct (`burn_t&`, `dvode_t&`, `template<class T> f(T&)`) | imports when every write through it is to a plain array member's elements (`s.xn[n] = ...`, r16), a whole-struct assignment (`r = rate_t{...}`, r21), or inside the struct's own member functions (`s.bump(dt)`, o13, o14); fails for scalar-field writes (`s.T = ...`, `+=`, `++`; r03, r15) and for writes through a member's accessor (`vs.y(n) += ...`, r14, r19): "mutation requires an owned automatic source binding" | Pointer (`burn_t*`, `dvode_t*`) for any function that writes through it; member accessors work through the pointer (`vs->y(n) += ...`, o15). Keep `&` for read-only use |
| References to buffer elements (`f(out[i])`), array reference parameters (`const double (&)[N]`), reference-returning free accessors (`double& at(A&, int)`), local references to array elements (`double& r = a[2]`) | import (r06, r07, r09, r17; hrx-system #2077) | Keep them |
| Overlapping mutable refs (`f(s, s.T)`), local reference to a scalar member (`double& r = s.x`), a `const T&` return that selects between two locals (`amrex::max`/`min`) | fail (r08, r10 at import; r12 imports, then AMDGPU/023 `memory_access.private_root`) | Pointers or explicit copies; min/max return by value (r18) |
| Member function calls (`a.f(x)`, `p->f(x)`), `operator()` returning `T&`/`const T&` | import (o01, o02, o10; hrx-system #2040) | Keep them |
| `Array1D/Array2D`-style structs (array member + 1-based `operator()`; no base class) | import by value, `&` and `const&`, with writes through the accessor (o11, o12) | Keep the 1-based accessor and the source's index expressions; write a local aggregate of the same shape in the kernel namespace (the kernel cannot include AMReX) |
| Lambdas (`constexpr_for`, `ParallelFor`) | fail | Indexed loops; template-struct recursion for compile-time unrolling |
| `switch` | fails at runtime (folds if `constexpr`) | if/else chain |
| range-for | fails | Indexed loops |
| `for (;;)`, `for (init;; step)` | imports (f04, f12; hrx-system #2105) | — |
| `return` in a loop | fails (LOWERING/059) | `done`/result code + `break`, return after the loop |
| Early `return` in an `if` nested inside another `if`, with code after both (stock `actual_log_screen`) | fails import (f24: "conditional returns require at most one fallthrough arm"); a single-level early return imports (f10, f25) | `bool done` flag, later code under `if (!done)`, one return at the end (`screen.H`) |
| Loop with several exits (`break`s/`continue` inside branches) under divergent control | compiles, single (f21) or nested (f22, the dvstep corrector shape: fails TARGET/034 `masked_region_merge_vgpr_values` at `fada6788a1` and before, compiles at `fa8eba4554`) | Keep the stock loop shape; only a `return` inside the loop needs rewriting (previous row) |
| Uninitialized locals | scalar and array locals import (`double r;`, `bool converged;`, `double a[14];`: d04, d13, d15; hrx-system #2105); struct and pointer locals fail (d14: "uninitialized locals require scalar or native vector object storage") | `T x{}` for structs; scalars and arrays may stay uninitialized as in the source |
| Base classes; overloaded arithmetic operators, free or member (autodiff `dual`) | fail (o06; o03, o09: "overloaded binary operation requires a Loom operation") | Flatten bases; POD dual structs with free functions (see translation guide, "Autodiff") |
| Constructors with arguments (`D x(1.0)`) | fail (o07: records must be complete aggregates) | Aggregate init `D x{1.0, 0.0}` |
| Structured bindings (`auto [f, r] = rhs_term(...)`) | fail (o08) | Return a 2-field struct; read `.first`/`.second` |
| `constexpr` local of record type (`constexpr rhs_t data = rhs_data(rate);`) | fails (f02: "unsupported C++ type"), though the `constexpr` function with a `switch` is fine | One `constexpr` scalar per field used: `constexpr int expA = rhs_data(rate).exponent_A;` |
| Enum in arithmetic with `double` | fails to parse (d10, PARSE/036): correct under the default `--std=c++26`; the pinned parser also rejects it in C++20/23 modes (reported fixed upstream as `0d27041f`, a commit not found in public hrx-system history as of 2026-10-18; re-probe d10) | `static_cast<int>(e)`; aprox13 doesn't use it |
| Floating-point contraction | none: `a * b + c` → `v_mul_f64` + `v_add_f64` under `--fp-model=precise` (m_fp_contract) | — |
| Namespace-scope non-constexpr globals (runtime params, tables) | fail | Kernel arguments: a params POD pointer and table pointers |
| Pointer stored in a struct, read through a pointer to that struct (`const Tables* t; t->f[i]`), including a kernel argument pointing at a struct of pointers | fails (c03, c06: "stored pointers and references require an object storage representation") | Pass the tables struct **by value** (c04) or one flat buffer with compile-time offsets (c05); kernel takes the table pointers as separate arguments |
| Namespace-scope `constexpr` arrays (`aion[]`, `zion[]`, rate tables; also `inline constexpr`, 2-D, arrays of aggregates, `extern const` with its definition in the TU) | import with runtime indices since hrx-system #2158/#2163 (d01, d16-d20, d23); stored as read-only data, read with global loads | Keep them at namespace scope as the stock HIP build does (`network_data.H` emits `aion`/`aion_inv`/`zion` this way) |
| `static constexpr` member arrays; block-scope `static constexpr` arrays | fail (d09: "member-dependent initialization ..."; d22: "local storage duration must be automatic or __shared__") | Move the table to namespace scope (row above), or a local `const Real t[N] = {...}` |
| Ternaries nested > 32 deep | import and compile (f11, 40 deep; hrx-system #2079) | — |
| Stage loops with break/continue | a simple retry loop imports without help (f07), and so does a constant-bound loop with a `break` over a private local array (f17); the same loop over a struct member read through a pointer (`vs->y(i)` in VODE's species checks) still fails SUBRANGE/026 at `fa8eba4554` | `MPL_ASSUME(0 <= i && i < N)` on the counter (next row) |
| Private-array index with data-dependent loop bounds (triangular LU loops, VODE order loops), an index read from a struct field (`vs->el[vs->L]`), or a loop with `break` over a member array reached through a pointer (its compound condition hides the bound, even with constant limits) | imports, then **compile** fails SUBRANGE/026 "access byte end may reach ..." even though probe f06 passed in isolation (network LU, 2026-10-18) | First choice for counted loops with compile-time limits: write the bound as a strict `<` (`for (int j = k + 1; j < neqs + 1; ++j)`, `for (int j = 1; j < k; ++j)` with `k` from another such loop); Loom proves those indices in range and no assumption is needed (network LU, `network_rhs.H`; with `<=`, dgesl's back substitution fails SUBRANGE/026; also at `31bf125137`). Otherwise (runtime bounds from struct fields, as in VODE): `MPL_ASSUME(lo <= i && i <= hi)` (support.H: `__builtin_assume` under Loom, no-op on host) at the top of each loop body / after each derived index, using the array's true static bounds. Operands must be local scalars, constants and side-effect-free arithmetic on them (`j + 1 <= n` imports since hrx-system #2105, f13); no struct fields (`vs->L` fails: "predicate member paths ..."; copy it into a `const int` local and index with that, f14). Rely on the smoke compile, not the isolated probe |

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
`#if defined(__loom__)` (or `#ifdef __loom__`), or as an `#elif` after it
(e.g. `#elif defined(MPL_HIP_STANDALONE)` in `support.H`); the lint skips
every branch of that chain except the first. (The
lint is clean on the accepted chemistry kernel, adv-mi210-12's seed.)

## Works

`if constexpr`; recursive `powi<N>` templates; static member function
templates; member functions called with `.` or `->` (non-virtual);
`operator()` accessors returning `T&`/`const T&` (`Array1D`/`Array2D` shape);
`enum class`; anonymous namespaces; structs with arrays by value;
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
