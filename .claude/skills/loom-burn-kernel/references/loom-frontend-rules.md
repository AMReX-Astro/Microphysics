# Loom C++ frontend: construct → rewrite

Status as of the pin, hrx-system main `920eb7f3ff` (probed 2026-11-07 on
gfx90a; the previous pins were main `cfd80e56de` and, before it, Loom 0.1.0,
tag `v0.1.0` = `d01d1c18b1`), measured
with this skill's `probes/` kit (results and probe names in
`probes/README.md`). **Re-probe at your commit** (SKILL.md Phase 1):
rows marked "changed recently" or anything older than the commit you use
may be stale.

Since `cfd80e56de`, ten constructs that failed before import, compile and
give results bit-identical to hipcc at `920eb7f3ff` (d14, d22, f02, l01,
o08, o09, r08, r10, r12, r25 and stricter variants; `probes/gpucheck/`), so
their rows below say "Keep it"; checked as single-construct kernels, not
yet inside a full burn. Base classes and pointers stored in structs still
fail. Since Loom 0.1.0, constructors with arguments (o07) and
assumptions that read a struct field (f14) import. Since
`fa8eba4554`, mutable struct references, a runtime `switch`, `return`
inside a loop, nested early returns and free operator overloads import; the
shipped modules still use the older rewrites, which remain valid. Tools from `1468c06034` through `970d39ac1e` miscompile burn
kernels (`zone_success` lost, ROCm/hrx-system#2337, fixed at `5e6af0ec1a`,
in 0.1.0 and at the pin); don't use them.

## Import invocation

Command lines for import and compile are in `test-react-integration.md`
("Import/compile commands"). The importer defaults to
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
| f64 rint/floor/trunc | work natively (`roundevenf` → `v_rndne_f64`) | — (`f64_math.H` uses `roundevenf` for OCML's rint) |
| f64 ldexp/frexp | native since hrx-system #2131: `loom::scalar::ldexpf`/`frexpf` or `__builtin_ldexp`/`__builtin_frexp` → `v_ldexp_f64`, `v_frexp_mant_f64`, `v_frexp_exp_i32_f64` (m_f64_ldexpf, m_f64_frexpf, m_f64_builtin_frexp_ldexp) | — (the OCML recipes use them, as OCML does) |
| `__builtin_rint`, `__builtin_rintf` | fail import (m_f64_builtin_rint) | `loom::scalar::roundevenf` (f64 and f32: `v_rndne_f64`/`v_rndne_f32`, m_f64_roundevenf, m_f32_roundevenf) |
| Raw f64 reciprocal (`__builtin_amdgcn_rcp`, OCML's `MATH_FAST_RCP`/`MATH_FAST_DIV` seed) | unknown builtin (m_f64_amdgcn_rcp); `approximate::divf` (afn) and `loom::divide_reciprocal` (arcp) are full IEEE divisions with `v_div_scale`/`v_div_fixup` (m_f64_div_afn, m_f64_divide_reciprocal) | IEEE division; the GPU bit check (in the `loom-burn-kernel-tools` follow-up) found no input where it differs from OCML's refined reciprocal |
| Scalar reference parameters (`double&`, `int&`) and scalar member slices (`f(s.T)`) | import, `static` or not (r01, r02, r04, r11, r13; #1993) | Keep them |
| `const` reference to a whole struct (`const burn_t&`, `const T&` template parameters) | imports (r05), including member accessor calls (`vs.y(n)` with `const dvode_t& vs`, r23), plain fields and raw array members (r24) | Keep it; a `const` pointer also works (`vs->y(n)`, o19) |
| Mutable reference to a whole struct (`burn_t&`, `dvode_t&`, `template<class T> f(T&)`) | imports, with writes to scalar fields (`s.T = ...`, r15), array members (r16), member accessors (`vs.y(n) += ...`, r14, r19), whole-struct assignment (r21) and member functions (o13, o14), through template parameters too (r03); failed for scalar-field and accessor writes at `fa8eba4554` | Keep it. Pointers (`burn_t*`, `dvode_t*`, as in the shipped modules) also work (o15) |
| References to buffer elements (`f(out[i])`), array reference parameters (`const double (&)[N]`), reference-returning free accessors (`double& at(A&, int)`), local references to array elements (`double& r = a[2]`) | import (r06, r07, r09, r17; hrx-system #2077) | Keep them |
| Overlapping mutable refs (`f(s, s.T)`), local reference to a scalar member (`double& r = s.x`), a `const T&` return that selects between two locals (`amrex::max`/`min`) | work since `920eb7f3ff` (r08, r08b, r10, r10b, r12, r12b; GPU check: bit-identical to hipcc, `probes/gpucheck/`); failed at `cfd80e56de` (r08, r10 at import; r12 imported, then AMDGPU/023 `memory_access.private_root`) | Keep them. Older ports use pointers or explicit copies and by-value min/max (r18), which still work |
| Member function calls (`a.f(x)`, `p->f(x)`), `operator()` returning `T&`/`const T&` | import (o01, o02, o10; hrx-system #2040) | Keep them |
| `Array1D/Array2D`-style structs (array member + 1-based `operator()`; no base class) | import by value, `&` and `const&`, with writes through the accessor (o11, o12) | Keep the 1-based accessor and the source's index expressions; write a local aggregate of the same shape in the kernel namespace (the kernel cannot include AMReX) |
| Lambdas (`constexpr_for`-style: a lambda passed to a template and called there) | work since `920eb7f3ff`, capturing by reference and by value (l01, l01b; GPU check: bit-identical to hipcc, `probes/gpucheck/`); failed at `cfd80e56de`. Not probed: lambdas stored or returned, `ParallelFor` (host-side; the kernel root replaces it) | Keep `constexpr_for` loops. Older ports use indexed loops or template-struct recursion, which still work |
| `switch` | imports and compiles at runtime (f01); failed at `fa8eba4554` | Keep it |
| range-for | fails | Indexed loops |
| `for (;;)`, `for (init;; step)` | imports (f04, f12; hrx-system #2105) | — |
| `return` in a loop | imports and compiles (f03); failed at `fa8eba4554` | Keep it |
| Early `return` in an `if` nested inside another `if`, with code after both (stock `actual_log_screen`) | imports and compiles (f24), as does a single-level early return (f10, f25); f24 failed at `fa8eba4554` | Keep it (`screen.H` keeps the `done` flag it needed at `fa8eba4554`) |
| Loop with several exits (`break`s/`continue` inside branches) under divergent control | compiles, single (f21) or nested (f22, the dvstep corrector shape) | Keep the stock loop shape; only a `return` inside the loop needs rewriting (previous row) |
| Uninitialized locals | scalar, array and struct locals import (`double r;`, `bool converged;`, `double a[14];`, `Burn s;`: d04, d13, d15, d14), and so does an uninitialized scalar passed as a `T&` out-argument (r25); d14 and r25 since `920eb7f3ff` (GPU check: bit-identical to hipcc, `probes/gpucheck/`; at `cfd80e56de` they failed with "uninitialized locals require scalar or native vector object storage" and "runtime value-reference arguments require an automatic SSA object"). Uninitialized pointer locals are not probed | Keep them uninitialized as in the source; initialize pointer locals |
| Base classes; member arithmetic operators (autodiff `dual`) | base classes fail (o06: "record values require a complete baseless non-union type without virtual members"); member operators work since `920eb7f3ff`, including `operator+=` returning `*this` (o09, o09b; GPU check: bit-identical to hipcc, `probes/gpucheck/`; at `cfd80e56de`: "addressed operator must return a reference"); free `operator*`/`operator+` on an aggregate import (o03) | Flatten bases. Duals: member or free operators both work; the shipped `dual.H` uses free functions (translation guide, "Autodiff") |
| Constructors with arguments (`D x(1.0)`) | import since `cfd80e56de` (o07; user-provided constructors on trivially copyable records, hrx-system #2523); failed in 0.1.0 | Keep them; the shipped modules use aggregate init `D x{1.0, 0.0}`, which also imports |
| Structured bindings (`auto [f, r] = rhs_term(...)`, `auto& [a, b] = p`) | work since `920eb7f3ff`, by value and by reference (o08, o08b; GPU check: bit-identical to hipcc, `probes/gpucheck/`); failed at `cfd80e56de` | Keep them |
| `constexpr` local of record type (`constexpr rhs_t data = rhs_data(rate);`) | works since `920eb7f3ff` (f02; GPU check: bit-identical to hipcc, `probes/gpucheck/`); failed at `cfd80e56de` ("unsupported C++ type") | Keep it. Older ports bind one `constexpr` scalar per field (`constexpr int expA = rhs_data(rate).exponent_A;`), which still works |
| Enum in arithmetic with `double` | fails to parse (d10, PARSE/036): correct under the default `--std=c++26`; the pinned parser also rejects it in C++20/23 modes (reported fixed upstream as `0d27041f`, a commit not found in public hrx-system history as of 2026-10-18; re-probe d10) | `static_cast<int>(e)`; aprox13 doesn't use it |
| Floating-point contraction | none: `a * b + c` → `v_mul_f64` + `v_add_f64` under `--fp-model=precise` (m_fp_contract) | — |
| Namespace-scope non-constexpr globals (runtime params, tables) | fail | Kernel arguments: a params POD pointer and table pointers |
| Pointer stored in a struct, read through a pointer to that struct (`const Tables* t; t->f[i]`), including a kernel argument pointing at a struct of pointers | fails (c03, c06: "stored pointers and references require an object storage representation") | Pass the tables struct **by value** (c04) or one flat buffer with compile-time offsets (c05); kernel takes the table pointers as separate arguments |
| Namespace-scope `constexpr` arrays (`aion[]`, `zion[]`, rate tables; also `inline constexpr`, 2-D, arrays of aggregates, `extern const` with its definition in the TU) | import with runtime indices since hrx-system #2158/#2163 (d01, d16-d20, d23); stored as read-only data, read with global loads | Keep them at namespace scope as the stock HIP build does (`network_data.H` emits `aion`/`aion_inv`/`zion` this way) |
| `static constexpr` member arrays; block-scope `static constexpr` arrays | member arrays fail (d09: "member-dependent initialization ..."); block-scope arrays work since `920eb7f3ff` (d22; GPU check: bit-identical to hipcc, `probes/gpucheck/`; at `cfd80e56de`: "local storage duration must be automatic or __shared__") | Member arrays: move the table to namespace scope (row above), or a local `const Real t[N] = {...}`. Block-scope arrays: keep them |
| Ternaries nested > 32 deep | import and compile (f11, 40 deep; hrx-system #2079) | — |
| Stage loops with break/continue | a simple retry loop imports without help (f07), and so does a constant-bound loop with a `break` over a private local array (f17); the same loop over a struct member read through a pointer (`vs->y(i)` in VODE's species checks) still fails SUBRANGE/026 at `920eb7f3ff` (full burn without the two assumptions, 2026-11-07; also at `cfd80e56de` and in 0.1.0) | `MPL_ASSUME(0 <= i && i < N)` on the counter (next row) |
| Private-array index with data-dependent loop bounds (triangular LU loops, VODE order loops), an index read from a struct field (`vs->el[vs->L]`), or a loop with `break` over a member array reached through a pointer (its compound condition hides the bound, even with constant limits) | imports, then **compile** fails SUBRANGE/026 "access byte end may reach ..." even though probe f06 passed in isolation (network LU, 2026-10-18) | First choice for counted loops with compile-time limits: write the bound as a strict `<` (`for (int j = k + 1; j < neqs + 1; ++j)`, `for (int j = 1; j < k; ++j)` with `k` from another such loop); Loom proves those indices in range and no assumption is needed (network LU, `network_rhs.H`; with `<=`, dgesl's back substitution fails SUBRANGE/026). Otherwise (runtime bounds from struct fields, as in VODE): `MPL_ASSUME(lo <= i && i <= hi)` (support.H: `__builtin_assume` under Loom, no-op on host) at the top of each loop body / after each derived index, using the array's true static bounds. Operands must be local scalars, constants and side-effect-free arithmetic on them (`j + 1 <= n` imports since hrx-system #2105, f13); since `cfd80e56de` also direct struct fields read through a pointer (`__builtin_assume(1 <= vs->L && vs->L <= 6)`, f14, hrx-system #2518; in 0.1.0 it failed "predicate member paths ...", so the shipped `vode.H` copies the field into a `const int` local first, which still works). Rely on the smoke compile, not the isolated probe |

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
when Phase 1 probes change what imports (relax a rule only after the
construct passes `probes/gpucheck/`, and add it to
`tests/should_not_fire.cpp`). Errors must be fixed; the `no-uninit-local`
warning (an uninitialized pointer local) means initialize the pointer.
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
