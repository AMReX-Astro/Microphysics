# Loom frontend probe kit

`src/` holds one-construct probes (root `probe::k`) for the constructs a
Microphysics burn uses; `probe.sh` imports and compiles each with one frozen
tool set on the GPU host and writes `probe-summary.txt`. See the header of
`probe.sh` for how to ship and run it. Add probes for anything a port's
inventory finds that is not covered here, and keep new probes in this kit.

## Results at hrx-system fada6788a1 (gfx90a, run 2026-10-20)

Public `ROCm/hrx-system` main `fada6788a1` (2026-10-20), tools built with
`scripts/build_loom_tools.sh` and both `patches/` applied (`loom-import-cxx`
md5 `e6d89e5b…`, `loom-compile` md5 `6179c3a7…`; md5s differ with the build
environment, so record yours). All 106 probes run in a few seconds. Rerun the
kit whenever the pinned commit changes and update this table.

What changed since the previous pin (`31bf125137`, 2026-10-19): hrx-system
#2105 (uninitialized scalar locals, conditionless `for`, signed counted
loops, assumptions over derived expressions). d04, f04 went from import-fail
to ok; d13, d15, f12, f13 are new and ok, d14 and f14 new and failing. f15
(no assumptions) still fails SUBRANGE/026 and f16 (with them) compiles.
f20-f23 are new: f20 guards a regression from #2074 that drops an assumed
range when the index was loaded from a `short` field (VODE's `L`, `NQ`), which
`patches/assume-cast-chain.patch` fixes (the unpatched tools fail f20 and the
aprox13 full burn); f21-f23 pin down which multi-exit loops still need the
single-exit form. Everything else is unchanged apart from a few VGPR counts.

What these results decide:
- **References:** scalar `double&`/`int&` parameters and scalar member slices
  (`f(s.T)`) import (r01, r02, r04, r11, r13), as do `const` references to
  whole structs (r05), references to buffer elements (r06), array reference
  parameters (r09), reference-returning accessors (r07), and local references
  to array elements (r17). Through a **mutable** reference to a whole struct,
  only three kinds of write import: elements of a plain array member
  (`s.xn[n] = ...`, r16), whole-struct assignment (`r = rate_t{...}`, r21),
  and calls to the struct's own member functions, which may write any field
  (`s.bump(dt)`, o13, o14). Everything else fails
  with "mutation requires an owned automatic source binding": scalar fields
  (`s.T = ...`, `+=`, `++s.n_rhs`; r03, r15), a member's accessor
  (`vs.y(n) += ...`, r14, r19), or passing a member on to a mutating
  reference; also with `static`, template, nested-call, buffer-element and
  `*p` arguments. Through a pointer all of these work, including member
  accessors (`vs->y(n) += ...`, o15). So `burn_t&`, `dvode_t&` and
  `rhs_state_t&` that are written still become pointers; read-only
  `const burn_t&` can stay. Calling a member's accessor through any struct
  reference fails, the `const` overload included (`vs.y(n)` with
  `const dvode_t& vs`, r23), while raw array members read fine (r24) and a
  `const` pointer works (`vs->y(n)`, o19): structs with `Array1D`/`Array2D`
  members always go by pointer. Still failing:
  overlapping mutable arguments `f(s, s.T)` (r08), a local reference to a
  scalar member `double& r = s.x` (r10), and a `const T&` return selecting
  between two locals (`amrex::max`, r12: imports, then AMDGPU/023); make
  min/max return by value (r18).
- **Member functions and accessors:** member calls through `.` and `->` (o01,
  o10) and `operator()` returning `T&`/`const T&` (o02) import. AMReX-style
  `Array1D<T,1,N>`/`Array2D` structs (an array member plus 1-based
  `operator()`) work as locals, by `&` and by `const&`, including writes
  through the accessor such as an in-place LU (o11, o12), and as members of
  a struct reached through a pointer (o15). Integrator and `rhs` state can
  keep its 1-based accessors and source index expressions.
- **Operator overloading** (`operator*`, `operator+`, free or member) still
  fails (o03, o09): autodiff duals stay POD structs with free functions (o04,
  o05, `assets/kernel/dual.H`).
- **Nested ternaries:** 40 deep compiles (f11); the old 32-level verifier
  limit is gone.
- **Math:** f64 exp, log, cbrt, pow, log10, tanh, sin, cos all fail and need
  OCML recipes (`assets/kernel/f64_math.H`: exp/log/cbrt;
  `f64_math_extra.H`: pow/log10/tanh; sin/cos not transcribed yet). sqrt, fma, floor,
  trunc, rint and copysign are native.
- **No contraction:** `a * b + c` compiles to `v_mul_f64` + `v_add_f64`.
- **`constexpr RhsT d = rhs_data(rate)` fails** (f02): a constexpr local of
  record type is unsupported, though the constexpr function with a `switch`
  is fine. Bind each needed field to its own scalar instead:
  `constexpr int expA = rhs_data(rate).exponent_A;`.
- **Enum in arithmetic with `double`** fails to parse (d10). That is correct
  at the importer's default `--std=c++26` (P2864R2 made it ill-formed); if a
  translation uses it, cast to `int`. Microphysics' aprox13 configuration
  doesn't use the construct.
- **Pointers stored in structs** can't be read through a pointer to the
  struct (c03, c06): pass a tables struct **by value** (c04), or use one flat
  buffer with compile-time offsets (c05). Kernel arguments must be the table
  pointers themselves, not a pointer to a struct holding them.
- **Records** must be plain aggregates: no constructors with arguments (o07),
  no bases (o06). Structured bindings fail (o08). Lambdas fail (l01).
- **Uninitialized locals and `for (;;)`:** scalar and array locals may stay
  uninitialized (d04, d13, d15) and conditionless `for` imports (f04, f12);
  struct locals still need `T x{}` (d14).
- **Bounds assumptions:** the simple retry loop imports without
  `__builtin_assume` (f07), and the isolated triangular loop compiles (f06),
  but an index derived from a struct field with runtime loop bounds still
  fails SUBRANGE/026 without assumptions (f15) and compiles with them (f16).
  Assumptions may now use arithmetic (`j + 1 <= n`, f13) but still not struct
  fields (f14). With unpatched `fada6788a1` tools, an assumption on a local
  loaded from a `short` field is dropped (f20): build with
  `patches/assume-cast-chain.patch`.
- **Loop exits:** `return` inside a loop still fails import (f03). A single
  loop with several `break`s and a `continue` compiles (f21, the dvode step
  loop's shape), but an inner loop with several `break`s inside an outer loop
  fails the AMDGPU backend with TARGET/034 (f22); the single-exit form of the
  inner loop compiles (f23).

| Probe | Construct | Result | Detail |
|---|---|---|---|
| c00_control | baseline | ok | vgpr=8 scratch=0B |
| c01_params_pod_ptr | `const Params*` runtime-param POD | ok | vgpr=7 scratch=0B |
| c02_table_ptr_3d | helmholtz `f[j][i][9]` via pointer, clamp, local `fi[36]` | ok | vgpr=16 scratch=288B |
| c03_struct_of_ptrs_by_ptr | struct of table pointers passed by pointer (`t->a[i]`) | import-fail | error [LOWERING/059]: stored pointers and references require an object storage representation |
| c04_struct_of_ptrs_by_value | struct of table pointers passed by value (`t.a[i]`) | ok | vgpr=6 scratch=0B |
| c05_flat_table_offsets | one table buffer, compile-time offsets | ok | vgpr=6 scratch=0B |
| c06_ptr_member_of_kernel_arg_struct | kernel argument: pointer to a struct of pointers in device memory | import-fail | error [LOWERING/059]: stored pointers and references require an object storage representation |
| d01_ns_constexpr_array | `aion[]`, `zion[]`, `aion_inv[]` | import-fail | error [LOWERING/059]: constants require a numeric scalar type |
| d02_local_const_array | local `const double t[N]` | ok | vgpr=16 scratch=104B |
| d03_ns_mutable_global | runtime param global | import-fail | error [LOWERING/059]: global storage definitions require a global-storage projection |
| d04_uninit_local | `bool converged;` | ok | vgpr=8 scratch=0B |
| d05_struct_default_init_braces | `T x{}` | ok | vgpr=8 scratch=0B |
| d06_struct_copy_assign | `old{s}; s = old;` | ok | vgpr=4 scratch=0B |
| d07_limits_and_nan | `__DBL_MAX__`, NaN via bit_cast | ok | vgpr=9 scratch=0B |
| d08_large_struct_runtime_index | dvode_t-sized local struct | ok | vgpr=10 scratch=2872B |
| d09_static_constexpr_member_array | `static constexpr` member array | import-fail | error [LOWERING/059]: member-dependent initialization requires source object initialization semantics |
| d10_enum_u8_base | `enum : unsigned char` | import-fail | error [PARSE/036]: invalid operands of types 'double' and '::probe::Rates' to binary operator '+' |
| d11_short_and_bool_arrays | `short pivot[14]`, `bool` arrays | ok | vgpr=10 scratch=42B |
| d12_f64_to_int_conversion | `int(double)` | ok | vgpr=10 scratch=0B |
| d13_uninit_local_double | `double r;` assigned on both branches | ok | vgpr=8 scratch=0B |
| d14_uninit_local_struct | `Burn s;` without initializer | import-fail | error [LOWERING/059]: uninitialized locals require scalar or native vector object storage |
| d15_uninit_local_array | `double a[14];` without initializer | ok | vgpr=6 scratch=112B |
| f01_switch_runtime | runtime `switch` | import-fail | error [LOWERING/059]: unsupported statement: switch-statement |
| f02_switch_in_constexpr_fn | `rhs_data(rate)` folded | import-fail | error [LOWERING/059]: unsupported C++ type: ::probe::RhsT |
| f03_return_in_loop | `return` in `while(true)` | import-fail | error [LOWERING/059]: returns inside loops require a structured loop-exit projection |
| f04_for_ever | `for(;;)` | ok | vgpr=10 scratch=0B |
| f05_vode_retry_shape | nested `while(true)`, result code | ok | vgpr=20 scratch=0B |
| f06_runtime_triangular_loops | `for k=NQ..1, j=k..NQ` | ok | vgpr=12 scratch=672B |
| f07_retry_loop_no_assume | retry loop with break/continue, no `__builtin_assume` | ok | vgpr=12 scratch=0B |
| f08_retry_loop_assume | same, with `__builtin_assume` | ok | vgpr=12 scratch=0B |
| f09_newton_for_break | EOS Newton `for` + `break` | ok | vgpr=26 scratch=0B |
| f10_early_return_outside_loop | sneut5 early return | ok | vgpr=8 scratch=0B |
| f11_ternary_nest_40 | 40-deep nested ternary | ok | vgpr=86 scratch=240B |
| f12_for_no_cond_continue | `for (m = 0;; ++m)` with `continue`/`break` | ok | vgpr=10 scratch=0B |
| f13_assume_derived_expr | `__builtin_assume(j + 1 <= n)` | ok | vgpr=10 scratch=120B |
| f14_assume_struct_field | `__builtin_assume(1 <= vs->L && vs->L <= 6)` | import-fail | error [LOWERING/059]: predicate member paths require direct non-static record fields |
| f15_runtime_bound_index_no_assume | dvjust-shape runtime-bound loops and `vs->el[L]`, no assumptions | compile-fail | error [SUBRANGE/026]: view.load access byte end may reach 17179869192, exceeding storage root byte extent 848 |
| f16_runtime_bound_index_assume | same as f15 with `__builtin_assume` on every index | ok | vgpr=10 scratch=848B |
| f17_break_loop_private_index_no_assume | constant-bound loop with data-dependent `break`, no assumptions | ok | vgpr=12 scratch=112B |
| f20_assume_short_field_member_index | `const int La = vs->L` (`short` field) assumed, indexes an `Array1D` member through a pointer | ok | vgpr=8 scratch=176B (with `patches/assume-cast-chain.patch`; unpatched: compile-fail SUBRANGE/026 "access byte end may reach 262152") |
| f21_multi_exit_loop | dvode step-loop shape: `while (true)`, 4 `break`s + `continue` in data-dependent branches | ok | vgpr=20 scratch=64B |
| f22_nested_multi_exit_loops | inner `while (true)` with 3 `break`s nested in a multi-exit outer loop (dvstep corrector shape) | compile-fail | error [TARGET/034]: rejected 'cfg.cond_br': branch lowering constraint 'masked_region_merge_vgpr_values' is not satisfied |
| f23_nested_single_exit_form | f22 with both loops in the single-exit `done` form | ok | vgpr=18 scratch=64B |
| l01_lambda_constexpr_for | lambdas | import-fail | error [LOWERING/059]: rvalue-reference parameters require consuming object identity |
| l02_template_struct_unroll | template-struct recursion | ok | vgpr=30 scratch=0B |
| m_f64_absf | f64 fabs | ok | vgpr=4 scratch=0B \| |
| m_f64_cbrtf | f64 cbrt | compile-fail | error [TARGET/001]: has no target-low contract |
| m_f64_copysignf | f64 copysign | ok | vgpr=6 scratch=0B \| |
| m_f64_cosf | f64 cos | compile-fail | error [LOWERING/036]: legalize-math policy amdgpu-math rejected scalar.cosf for cosf in scalar lanes of f64 under math.element.f32 |
| m_f64_expf | f64 exp | compile-fail | error [LOWERING/036]: legalize-math policy amdgpu-math rejected scalar.expf for expf in scalar lanes of f64 under math.element.f32 |
| m_f64_floorf | f64 floor | ok | vgpr=4 scratch=0B \| v_floor_f64:1 |
| m_f64_fmaf | f64 fma | ok | vgpr=8 scratch=0B \| v_fma_f64:1 |
| m_f64_log10f | f64 log10 | compile-fail | error [TARGET/001]: has no target-low contrac |
| m_f64_logf | f64 log | compile-fail | error [LOWERING/036]: legalize-math policy amdgpu-math rejected scalar.logf for logf in scalar lanes of f64 under math.element.f32 |
| m_f64_math_h_recipes | existing exp/log/cbrt recipes | ok | vgpr=56 scratch=0B \| v_add_f64:49 v_fma_f64:36 v_mul_f64:19 v_div_scale_f64:4 v_cmp_lt_f64:3 v_rcp_f64:2 |
| m_f64_powf | f64 pow (non-integer exponent) | compile-fail | error [LOWERING/036]: legalize-math policy amdgpu-math rejected scalar.powf for powf in scalar lanes of f64 under math.element.f32 |
| m_f64_roundevenf | f64 rint (roundeven) | ok | vgpr=4 scratch=0B \| v_rndne_f64:1 |
| m_f64_sinf | f64 sin | compile-fail | error [LOWERING/036]: legalize-math policy amdgpu-math rejected scalar.sinf for sinf in scalar lanes of f64 under math.element.f32 |
| m_f64_sqrtf | f64 sqrt | ok | vgpr=10 scratch=0B \| v_fma_f64:7 v_mul_f64:2 v_ldexp_f64:2 v_rsq_f64:1 v_cmp_lt_f64:1 v_cmp_class_f64:1 |
| m_f64_tanhf | f64 tanh | compile-fail | error [LOWERING/036]: legalize-math policy amdgpu-math rejected scalar.tanhf for tanhf in scalar lanes of f64 under math.element.f32 |
| m_f64_truncf | f64 trunc | ok | vgpr=4 scratch=0B \| v_trunc_f64:1 |
| m_fp_contract | `a * b + c` under `--fp-model=precise` | ok | vgpr=8 scratch=0B \| v_mul_f64:1 v_add_f64:1 |
| m_powi_template | `powi<N>` recursion | ok | vgpr=8 scratch=0B \| v_mul_f64:5 v_add_f64:1 |
| o01_member_fn_call | `a.add(x)` | ok | vgpr=6 scratch=0B |
| o02_operator_call_index | `Array1D::operator()` | ok | vgpr=6 scratch=0B |
| o03_dual_operator_overload | `operator*`/`operator+` on Dual | import-fail | error [LOWERING/059]: overloaded binary operation requires a Loom operation |
| o04_dual_free_functions | POD Dual / Dual3 + free fns, by value | ok | vgpr=10 scratch=0B |
| o05_dual3_array_free_functions | POD Dual / Dual3 + free fns, by value | ok | vgpr=12 scratch=96B |
| o06_base_class | `Grad : Base` | import-fail | error [LOWERING/059]: record values require a complete aggregate without unions, bases or nontrivial lifecycle operations |
| o07_ctor_with_args | `D x(1.0)` | import-fail | error [LOWERING/059]: record values require a complete aggregate without unions, bases or nontrivial lifecycle operations |
| o08_structured_binding | `auto [f, r] = term(x)` | import-fail | error [LOWERING/059]: unsupported local declaration |
| o09_dual_member_operators | member `operator*`/`operator+`/`operator+=` on Dual | import-fail | error [LOWERING/059]: overloaded binary operation requires a Loom operation |
| o10_arrow_member_call | `a->add(x)`, `a.mean()` through a pointer to a local | ok | vgpr=14 scratch=0B |
| o11_array1d_by_ref | `Array1D<T,1,N>` accessor, passed by `&` and `const auto&` | ok | vgpr=10 scratch=232B |
| o12_array2d_lu | `Array2D` accessor, written in place through `&` (dgefa shape) | ok | vgpr=18 scratch=144B |
| o13_member_fn_mutates_via_struct_ref | `s.bump(dt)` mutating scalar fields, `s` a `Burn&` | ok | vgpr=8 scratch=0B |
| o14_struct_ref_member_fns | `vs.add_y(n, v); vs.tick(h)` through `Dvode&` | ok | vgpr=10 scratch=128B |
| o15_ptr_member_accessor_write | `vs->y(n) +=`, `vs->t +=` through a pointer | ok | vgpr=10 scratch=128B |
| o16_ptr_member_to_accessor_ref | `upd(vs->y)`: Array1D member through a pointer into `Array1D&`, written | ok | vgpr=10 scratch=1936B |
| o17_array2d_member_lu_runtime_bounds | `dgefa(vs->jac, nn)` with runtime bounds + `__builtin_assume` | ok | vgpr=18 scratch=1936B |
| o18_array1d_of_rate_dual | `Array1D<rate_t<Dual>, 1, N>` element assign/copy | ok | vgpr=20 scratch=0B |
| o19_const_ptr_member_accessor_read | `vs->y(n)` through `const Dvode*` | ok | vgpr=8 scratch=248B |
| r01_static_scalar_ref_out | `static f(double&)` on locals | ok | vgpr=8 scratch=0B |
| r02_inline_external_ref | `inline` (AMREX_INLINE) non-static fn with scalar `double&` | ok | vgpr=6 scratch=0B |
| r03_template_struct_ref | `template<class BurnT> inline f(BurnT&)` | import-fail | error [LOWERING/059]: mutation requires an owned automatic source binding |
| r04_anon_ns_ref | anon-namespace fn with scalar `double&` | ok | vgpr=6 scratch=0B |
| r05_const_struct_ref_read | `static f(const Burn&)` | ok | vgpr=8 scratch=112B |
| r06_ref_to_buffer_element | `f(out[i])` | ok | vgpr=6 scratch=0B |
| r07_ref_return_accessor | `double& at(A&, int)` | ok | vgpr=4 scratch=0B |
| r08_overlapping_mutable_refs | `f(s, s.T)` | import-fail | error [LOWERING/059]: mutable value-reference arguments must select nonoverlapping automatic objects or member slices |
| r09_array_ref_param | `const double (&)[13]` | ok | vgpr=6 scratch=0B |
| r10_local_ref_binding | `double& r = s.x;` | import-fail | error [LOWERING/059]: address-of requires an existing storage-backed element |
| r11_ref_chain | ref forwarded through 2 static calls | ok | vgpr=10 scratch=0B |
| r12_const_ref_return_max | `amrex::max` returning `const T&` | compile-fail | error [AMDGPU/023]: rejected 'view.load': source-to-low constraint 'memory_access.private_root' is not satisfied |
| r13_struct_member_slice_ref | `bump(s.T)` | ok | vgpr=8 scratch=0B |
| r14_struct_ref_member_update | `dvode_t&` through 2 template layers: `vs.y(n) +=` and `vs.t +=` | import-fail | error [LOWERING/059]: mutation requires an owned automatic source binding |
| r15_struct_ref_scalar_field_write | `s.T = ...` through `Burn&` | import-fail | error [LOWERING/059]: mutation requires an owned automatic source binding |
| r16_struct_ref_array_member_write | `s.xn[n] = ...` through `Burn&` | ok | vgpr=8 scratch=128B |
| r17_local_ref_array_elem | `double& r = a[2];` | ok | vgpr=4 scratch=0B |
| r18_value_return_max | `max(const T&, const T&)` returning `T` by value | ok | vgpr=6 scratch=0B |
| r19_struct_ref_member_accessor_write | `vs.y(n) += ...` through `Dvode&` | import-fail | error [LOWERING/059]: mutation requires an owned automatic source binding |
| r20_const_ref_from_deref | `const Burn&` bound to `*p` | ok | vgpr=10 scratch=0B |
| r21_whole_struct_assign_through_ref | `r = rate_t{...}` through `rate_t&` | ok | vgpr=6 scratch=0B |
| r22_const_ref_to_kernel_arg_struct | `const BurnParams&` bound to `*params` (kernel argument) | ok | vgpr=12 scratch=0B |
| r23_const_struct_ref_member_accessor_read | `vs.y(n)` through `const Dvode&` | import-fail | error [LOWERING/059]: mutation requires an owned automatic source binding |
| r24_const_struct_ref_member_raw_array | `vs.y.arr[n]` through `const Dvode&` | ok | vgpr=8 scratch=248B |
