# Loom frontend probe kit

`src/` holds one-construct probes (root `probe::k`) for the constructs a
Microphysics burn uses; `probe.sh` imports and compiles each with one frozen
tool set on the GPU host and writes `probe-summary.txt`. See the header of
`probe.sh` for how to ship and run it. Add probes for anything a port's
inventory finds that is not covered here, and keep new probes in this kit.

## Results at hrx-system 5e6af0ec1a (gfx90a, run 2026-10-27)

Public `ROCm/hrx-system` main `5e6af0ec1a` (2026-10-27), stock tools built
with `scripts/build_loom_tools.sh` and no patches (`loom-import-cxx` md5
`986c76d4…`, `loom-compile` md5 `a0f2186a…`; md5s differ with the build
environment, so record yours). All 124 probes run in a few seconds. Rerun the
kit whenever the pinned commit changes and update this table; its last column
holds the result at the previous pin, `fa8eba4554`.

The full aprox13 + helmholtz + VODE burn (neutrinos off, 4096 zones) built
with these tools is bit-exact against stock HIP on the same GPU: state,
`n_rhs`/`n_step` and success in every zone, for the default inputs, the prad
limiter and `unit_test.tmax=1.e-15`. That includes the `zone_success`
miscompile from `1468c06034` through `970d39ac1e` (ROCm/hrx-system#2337),
which is fixed here; don't use tools from that range.

Changes from `fa8eba4554` (the shipped modules were written for that pin and
still import and validate unchanged; the rewrites below are no longer
required for new code):
- Writes through a **mutable reference to a whole struct** import: scalar
  fields (r15), member accessors (r14, r19) and template `BurnT&` parameters
  (r03). A member accessor through a `const` struct reference imports too
  (r23).
- Free operator overloads on a POD dual import (o03); member operators still
  fail (o09, with a new message).
- A runtime `switch` (f01), `return` inside a loop (f03) and an early
  `return` nested two `if`s deep (f24) import and compile.
- Two VGPR counts moved by one or two (d08, o17); nothing else changed.

What these results decide (at `5e6af0ec1a`):
- **References:** scalar `double&`/`int&` parameters and scalar member slices
  (`f(s.T)`) import (r01, r02, r04, r11, r13), as do `const` references to
  whole structs (r05), references to buffer elements (r06), array reference
  parameters (r09), reference-returning accessors (r07), and local references
  to array elements (r17). Through a **mutable** reference to a whole struct,
  writes to scalar fields (`s.T = ...`, r15), array members (r16), member
  accessors (`vs.y(n) += ...`, r14, r19), whole-struct assignment (r21) and
  member-function calls (o13, o14) all import, through plain and template
  parameters (r03), and so does a member's accessor read through a `const`
  struct reference (r23). Pointers work as well (o15, o19), so the shipped
  modules' pointer parameters for written structs can stay. Still failing:
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
- **Operator overloading:** free `operator*`/`operator+` on a POD dual
  import (o03); member operators still fail (o09). The shipped autodiff
  duals keep POD structs with free functions (o04, o05,
  `assets/kernel/dual.H`).
- **Nested ternaries:** 40 deep compiles (f11).
- **Math:** f64 exp, log, cbrt, pow, log10, tanh, sin, cos all fail and need
  OCML recipes (`assets/kernel/f64_math.H`: exp/log/cbrt;
  `f64_math_extra.H`: pow/log10/tanh; sin/cos not transcribed yet). sqrt, fma, floor,
  trunc, rint (`roundevenf`), copysign, ldexp and frexp are native, so the
  recipes call the same instructions as OCML; only OCML's fast reciprocal
  (`v_rcp_f64` + Newton steps) is replaced by IEEE division.
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
- **Tables:** namespace-scope `constexpr` arrays import with runtime indices
  (d01, d16-d20, d23), so `aion[]`, `zion[]`, `aion_inv[]` and rate tables keep
  the stock namespace-scope form; `static` arrays (member or block scope)
  still fail (d09, d22) and move to namespace scope.
- **Uninitialized locals and `for (;;)`:** scalar and array locals may stay
  uninitialized (d04, d13, d15) and conditionless `for` imports (f04, f12);
  struct locals still need `T x{}` (d14).
- **Bounds assumptions:** the simple retry loop imports without
  `__builtin_assume` (f07), and the isolated triangular loop compiles (f06),
  but an index derived from a struct field with runtime loop bounds still
  fails SUBRANGE/026 without assumptions (f15) and compiles with them (f16).
  Assumptions may now use arithmetic (`j + 1 <= n`, f13) but still not struct
  fields (f14). f17 (a `break` loop over a private local array) passes
  without assumptions, but the same loop over `vs->y(i)` in the full VODE
  burn failed SUBRANGE/026 at `fa8eba4554` (not re-checked since), so VODE
  keeps those `MPL_ASSUME`s. An assumption on a local loaded from a `short` field holds
  (f20).
- **Loop exits:** `return` inside a loop imports (f03), and so does a
  runtime `switch` (f01). Loops
  with several `break`s and a `continue` compile, alone (f21, the dvode step
  loop's shape) or nested (f22, the dvstep corrector shape), so they keep the
  stock shape.
- **Early returns:** a `return` outside loops imports in a single `if` (f10,
  f25) and inside an `if` nested in another `if` with code after both (f24).
  `screen.H`'s `actual_log_screen` keeps the `done` flag it needed at
  `fa8eba4554`.

| Probe | Construct | Result | Detail | At previous pin `fa8eba4554` |
|---|---|---|---|---|
| c00_control | baseline | ok | vgpr=8 scratch=0B | same |
| c01_params_pod_ptr | `const Params*` runtime-param POD | ok | vgpr=7 scratch=0B | same |
| c02_table_ptr_3d | helmholtz `f[j][i][9]` via pointer, clamp, local `fi[36]` | ok | vgpr=16 scratch=288B | same |
| c03_struct_of_ptrs_by_ptr | struct of table pointers passed by pointer (`t->a[i]`) | import-fail | error [LOWERING/059]: stored pointers and references require an object storage representation | same |
| c04_struct_of_ptrs_by_value | struct of table pointers passed by value (`t.a[i]`) | ok | vgpr=6 scratch=0B | same |
| c05_flat_table_offsets | one table buffer, compile-time offsets | ok | vgpr=6 scratch=0B | same |
| c06_ptr_member_of_kernel_arg_struct | kernel argument: pointer to a struct of pointers in device memory | import-fail | error [LOWERING/059]: stored pointers and references require an object storage representation | same |
| d01_ns_constexpr_array | `aion[]`, `zion[]`, `aion_inv[]` | ok | vgpr=14 scratch=0B | same |
| d02_local_const_array | local `const double t[N]` | ok | vgpr=28 scratch=104B | same |
| d03_ns_mutable_global | runtime param global | import-fail | error [LOWERING/059]: global storage definitions require a global-storage projection | same |
| d04_uninit_local | `bool converged;` | ok | vgpr=8 scratch=0B | same |
| d05_struct_default_init_braces | `T x{}` | ok | vgpr=8 scratch=0B | same |
| d06_struct_copy_assign | `old{s}; s = old;` | ok | vgpr=4 scratch=0B | same |
| d07_limits_and_nan | `__DBL_MAX__`, NaN via bit_cast | ok | vgpr=9 scratch=0B | same |
| d08_large_struct_runtime_index | dvode_t-sized local struct | ok | vgpr=12 scratch=2872B | vgpr=10 scratch=2872B |
| d09_static_constexpr_member_array | `static constexpr` member array | import-fail | error [LOWERING/059]: member-dependent initialization requires source object initialization semantics | same |
| d10_enum_u8_base | `enum : unsigned char` | import-fail | error [PARSE/036]: invalid operands of types 'double' and '::probe::Rates' to binary operator '+' | same |
| d11_short_and_bool_arrays | `short pivot[14]`, `bool` arrays | ok | vgpr=10 scratch=42B | same |
| d12_f64_to_int_conversion | `int(double)` | ok | vgpr=10 scratch=0B | same |
| d13_uninit_local_double | `double r;` assigned on both branches | ok | vgpr=8 scratch=0B | same |
| d14_uninit_local_struct | `Burn s;` without initializer | import-fail | error [LOWERING/059]: uninitialized locals require scalar or native vector object storage | same |
| d15_uninit_local_array | `double a[14];` without initializer | ok | vgpr=6 scratch=112B | same |
| d16_ns_inline_constexpr_array | `inline constexpr double aion[N]` at namespace scope | ok | vgpr=14 scratch=0B | same |
| d17_ns_constexpr_2d_array | namespace-scope `constexpr` 2-D table, two runtime indices | ok | vgpr=8 scratch=0B | same |
| d18_ns_constexpr_struct_array | namespace-scope `constexpr` array of aggregates, runtime index | ok | vgpr=10 scratch=0B | same |
| d19_ns_constexpr_array_loop_index | `aion[n]`, `zion[n]` read in a constant-bound loop | ok | vgpr=20 scratch=0B | same |
| d20_ns_constexpr_array_unreachable | namespace-scope `constexpr` array the root never reads | ok | vgpr=6 scratch=0B | same |
| d21_ns_mutable_global_unreachable | mutable namespace-scope global the root never reads | import-fail | error [LOWERING/059]: global storage definitions require a global-storage projection | same |
| d22_fn_static_constexpr_array | block-scope `static constexpr` array in a device function | import-fail | error [LOWERING/059]: local storage duration must be automatic or __shared__ | same |
| d23_extern_const_array_defined | `extern const` array declared, then defined in the TU | ok | vgpr=14 scratch=0B | same |
| f01_switch_runtime | runtime `switch` | **ok** | vgpr=8 scratch=0B | import-fail: error [LOWERING/059]: unsupported statement: switch-statement |
| f02_switch_in_constexpr_fn | `rhs_data(rate)` folded | import-fail | error [LOWERING/059]: unsupported C++ type: ::probe::RhsT | same |
| f03_return_in_loop | `return` in `while(true)` | **ok** | vgpr=9 scratch=8B | import-fail: error [LOWERING/059]: returns inside loops require a structured loop-exit projection |
| f04_for_ever | `for(;;)` | ok | vgpr=10 scratch=0B | same |
| f05_vode_retry_shape | nested `while(true)`, result code | ok | vgpr=22 scratch=0B | same |
| f06_runtime_triangular_loops | `for k=NQ..1, j=k..NQ` | ok | vgpr=12 scratch=672B | same |
| f07_retry_loop_no_assume | retry loop with break/continue, no `__builtin_assume` | ok | vgpr=12 scratch=0B | same |
| f08_retry_loop_assume | same, with `__builtin_assume` | ok | vgpr=12 scratch=0B | same |
| f09_newton_for_break | EOS Newton `for` + `break` | ok | vgpr=26 scratch=0B | same |
| f10_early_return_outside_loop | sneut5 early return | ok | vgpr=8 scratch=0B | same |
| f11_ternary_nest_40 | 40-deep nested ternary | ok | vgpr=86 scratch=240B | same |
| f12_for_no_cond_continue | `for (m = 0;; ++m)` with `continue`/`break` | ok | vgpr=10 scratch=0B | same |
| f13_assume_derived_expr | `__builtin_assume(j + 1 <= n)` | ok | vgpr=10 scratch=120B | same |
| f14_assume_struct_field | `__builtin_assume(1 <= vs->L && vs->L <= 6)` | import-fail | error [LOWERING/059]: predicate member paths require direct non-static record fields | same |
| f15_runtime_bound_index_no_assume | dvjust-shape runtime-bound loops and `vs->el[L]`, no assumptions | compile-fail | error [SUBRANGE/026]: view.load access byte end may reach 17179869192, exceeding storage root byte extent 848 | same |
| f16_runtime_bound_index_assume | same as f15 with `__builtin_assume` on every index | ok | vgpr=12 scratch=848B | same |
| f17_break_loop_private_index_no_assume | constant-bound loop with data-dependent `break`, no assumptions | ok | vgpr=12 scratch=112B | same |
| f20_assume_short_field_member_index | `const int La = vs->L` (`short` field) assumed, indexes an `Array1D` member through a pointer | ok | vgpr=8 scratch=176B | same |
| f21_multi_exit_loop | dvode step-loop shape: `while (true)`, 4 `break`s + `continue` in data-dependent branches | ok | vgpr=22 scratch=64B | same |
| f22_nested_multi_exit_loops | inner `while (true)` with 3 `break`s nested in a multi-exit outer loop (dvstep corrector shape) | ok | vgpr=20 scratch=64B | same |
| f23_nested_single_exit_form | f22 with both loops in the single-exit `done` form | ok | vgpr=18 scratch=64B | same |
| f24_nested_conditional_return | `return` in an `if` nested in a runtime `if`, code after both (`actual_log_screen`) | **ok** | vgpr=8 scratch=0B | import-fail: error [LOWERING/059]: conditional returns require at most one fallthrough arm; shared continuations need a scoped exit projectio |
| f25_single_level_conditional_return | f24 with the conditions merged into one `if` | ok | vgpr=8 scratch=0B | same |
| l01_lambda_constexpr_for | lambdas | import-fail | error [LOWERING/059]: rvalue-reference parameters require consuming object identity | same |
| l02_template_struct_unroll | template-struct recursion | ok | vgpr=30 scratch=0B | same |
| m_f32_roundevenf | f32 `loom::scalar::roundevenf` (OCML `BUILTIN_RINT_F32`, cbrt) | ok | vgpr=2 scratch=0B \| v_rndne_f32:1 | same |
| m_f64_absf | f64 fabs | ok | vgpr=4 scratch=0B | same |
| m_f64_amdgcn_rcp | `__builtin_amdgcn_rcp` on f64 (OCML `MATH_FAST_RCP` seed) | import-fail | error [PARSE/036]: unknown builtin function '__builtin_amdgcn_rcp' | same |
| m_f64_builtin_frexp_ldexp | `__builtin_frexp` with an `int*` out-parameter, then `__builtin_ldexp` (OCML spelling) | ok | vgpr=5 scratch=0B \| v_ldexp_f64:1 v_frexp_mant_f64:1 v_frexp_exp_i32_f64:1 | same |
| m_f64_builtin_rint | `__builtin_rint` on f64 (OCML `BUILTIN_RINT_F64`) | import-fail | error [LOWERING/059]: call requires a concrete namespace-scope function declaration | same |
| m_f64_cbrtf | f64 cbrt | compile-fail | error [TARGET/001]: has no target-low contract | same |
| m_f64_copysignf | f64 copysign | ok | vgpr=6 scratch=0B | same |
| m_f64_cosf | f64 cos | compile-fail | error [LOWERING/036]: legalize-math policy amdgpu-math rejected scalar.cosf for cosf in scalar lanes of f64 under math.element.f32 | same |
| m_f64_div_afn | f64 `loom::scalar::approximate::divf` (afn) | ok | vgpr=14 scratch=0B \| v_fma_f64:5 v_div_scale_f64:2 v_rcp_f64:1 v_mul_f64:1 v_div_fmas_f64:1 v_div_fixup_f64:1 | same |
| m_f64_divide_reciprocal | f64 `loom::divide_reciprocal` (arcp, hrx-system #2131) | ok | vgpr=14 scratch=0B \| v_fma_f64:3 v_div_scale_f64:2 v_rcp_f64:1 v_mul_f64:1 v_div_fmas_f64:1 v_div_fixup_f64:1 | same |
| m_f64_expf | f64 exp | compile-fail | error [LOWERING/036]: legalize-math policy amdgpu-math rejected scalar.expf for expf in scalar lanes of f64 under math.element.f32 | same |
| m_f64_floorf | f64 floor | ok | vgpr=4 scratch=0B \| v_floor_f64:1 | same |
| m_f64_fmaf | f64 fma | ok | vgpr=8 scratch=0B \| v_fma_f64:1 | same |
| m_f64_frexpf | f64 `loom::scalar::frexpf` (hrx-system #2131) | ok | vgpr=6 scratch=0B \| v_frexp_mant_f64:1 v_frexp_exp_i32_f64:1 | same |
| m_f64_ldexpf | f64 `loom::scalar::ldexpf` (hrx-system #2131) | ok | vgpr=4 scratch=0B \| v_ldexp_f64:1 | same |
| m_f64_log10f | f64 log10 | compile-fail | error [TARGET/001]: has no target-low contrac | same |
| m_f64_logf | f64 log | compile-fail | error [LOWERING/036]: legalize-math policy amdgpu-math rejected scalar.logf for logf in scalar lanes of f64 under math.element.f32 | same |
| m_f64_math_h_recipes | `f64_math.H` exp/log/cbrt recipes | ok | vgpr=50 scratch=0B \| v_add_f64:47 v_fma_f64:36 v_mul_f64:11 v_ldexp_f64:6 v_div_scale_f64:4 v_cmp_lt_f64:3 | same |
| m_f64_powf | f64 pow (non-integer exponent) | compile-fail | error [LOWERING/036]: legalize-math policy amdgpu-math rejected scalar.powf for powf in scalar lanes of f64 under math.element.f32 | same |
| m_f64_roundevenf | f64 rint (roundeven) | ok | vgpr=4 scratch=0B \| v_rndne_f64:1 | same |
| m_f64_sinf | f64 sin | compile-fail | error [LOWERING/036]: legalize-math policy amdgpu-math rejected scalar.sinf for sinf in scalar lanes of f64 under math.element.f32 | same |
| m_f64_sqrtf | f64 sqrt | ok | vgpr=10 scratch=0B \| v_fma_f64:7 v_mul_f64:2 v_ldexp_f64:2 v_rsq_f64:1 v_cmp_lt_f64:1 v_cmp_class_f64:1 | same |
| m_f64_tanhf | f64 tanh | compile-fail | error [LOWERING/036]: legalize-math policy amdgpu-math rejected scalar.tanhf for tanhf in scalar lanes of f64 under math.element.f32 | same |
| m_f64_truncf | f64 trunc | ok | vgpr=4 scratch=0B \| v_trunc_f64:1 | same |
| m_fp_contract | `a * b + c` under `--fp-model=precise` | ok | vgpr=8 scratch=0B \| v_mul_f64:1 v_add_f64:1 | same |
| m_powi_template | `powi<N>` recursion | ok | vgpr=8 scratch=0B \| v_mul_f64:5 v_add_f64:1 | same |
| o01_member_fn_call | `a.add(x)` | ok | vgpr=6 scratch=0B | same |
| o02_operator_call_index | `Array1D::operator()` | ok | vgpr=6 scratch=0B | same |
| o03_dual_operator_overload | `operator*`/`operator+` on Dual | **ok** | vgpr=10 scratch=0B | import-fail: error [LOWERING/059]: overloaded binary operation requires a Loom operation |
| o04_dual_free_functions | POD Dual / Dual3 + free fns, by value | ok | vgpr=10 scratch=0B | same |
| o05_dual3_array_free_functions | POD Dual / Dual3 + free fns, by value | ok | vgpr=12 scratch=96B | same |
| o06_base_class | `Grad : Base` | import-fail | error [LOWERING/059]: record values require a complete aggregate without unions, bases or nontrivial lifecycle operations | same |
| o07_ctor_with_args | `D x(1.0)` | import-fail | error [LOWERING/059]: record values require a complete aggregate without unions, bases or nontrivial lifecycle operations | same |
| o08_structured_binding | `auto [f, r] = term(x)` | import-fail | error [LOWERING/059]: unsupported local declaration | same |
| o09_dual_member_operators | member `operator*`/`operator+`/`operator+=` on Dual | import-fail | error [LOWERING/059]: addressed operator must return an lvalue reference | import-fail: error [LOWERING/059]: overloaded binary operation requires a Loom operation |
| o10_arrow_member_call | `a->add(x)`, `a.mean()` through a pointer to a local | ok | vgpr=14 scratch=0B | same |
| o11_array1d_by_ref | `Array1D<T,1,N>` accessor, passed by `&` and `const auto&` | ok | vgpr=10 scratch=232B | same |
| o12_array2d_lu | `Array2D` accessor, written in place through `&` (dgefa shape) | ok | vgpr=14 scratch=144B | same |
| o13_member_fn_mutates_via_struct_ref | `s.bump(dt)` mutating scalar fields, `s` a `Burn&` | ok | vgpr=8 scratch=0B | same |
| o14_struct_ref_member_fns | `vs.add_y(n, v); vs.tick(h)` through `Dvode&` | ok | vgpr=10 scratch=128B | same |
| o15_ptr_member_accessor_write | `vs->y(n) +=`, `vs->t +=` through a pointer | ok | vgpr=10 scratch=128B | same |
| o16_ptr_member_to_accessor_ref | `upd(vs->y)`: Array1D member through a pointer into `Array1D&`, written | ok | vgpr=10 scratch=1936B | same |
| o17_array2d_member_lu_runtime_bounds | `dgefa(vs->jac, nn)` with runtime bounds + `__builtin_assume` | ok | vgpr=19 scratch=1936B | vgpr=18 scratch=1936B |
| o18_array1d_of_rate_dual | `Array1D<rate_t<Dual>, 1, N>` element assign/copy | ok | vgpr=20 scratch=0B | same |
| o19_const_ptr_member_accessor_read | `vs->y(n)` through `const Dvode*` | ok | vgpr=8 scratch=248B | same |
| r01_static_scalar_ref_out | `static f(double&)` on locals | ok | vgpr=8 scratch=0B | same |
| r02_inline_external_ref | `inline` (AMREX_INLINE) non-static fn with scalar `double&` | ok | vgpr=6 scratch=0B | same |
| r03_template_struct_ref | `template<class BurnT> inline f(BurnT&)` | **ok** | vgpr=12 scratch=128B | import-fail: error [LOWERING/059]: mutation requires an owned automatic source binding |
| r04_anon_ns_ref | anon-namespace fn with scalar `double&` | ok | vgpr=6 scratch=0B | same |
| r05_const_struct_ref_read | `static f(const Burn&)` | ok | vgpr=8 scratch=112B | same |
| r06_ref_to_buffer_element | `f(out[i])` | ok | vgpr=6 scratch=0B | same |
| r07_ref_return_accessor | `double& at(A&, int)` | ok | vgpr=4 scratch=0B | same |
| r08_overlapping_mutable_refs | `f(s, s.T)` | import-fail | error [LOWERING/059]: mutable value-reference arguments must select nonoverlapping automatic objects or member slices | same |
| r09_array_ref_param | `const double (&)[13]` | ok | vgpr=6 scratch=0B | same |
| r10_local_ref_binding | `double& r = s.x;` | import-fail | error [LOWERING/059]: address-of requires an existing storage-backed element | same |
| r11_ref_chain | ref forwarded through 2 static calls | ok | vgpr=10 scratch=0B | same |
| r12_const_ref_return_max | `amrex::max` returning `const T&` | compile-fail | error [AMDGPU/023]: reject | same |
| r13_struct_member_slice_ref | `bump(s.T)` | ok | vgpr=8 scratch=0B | same |
| r14_struct_ref_member_update | `dvode_t&` through 2 template layers: `vs.y(n) +=` and `vs.t +=` | **ok** | vgpr=10 scratch=128B | import-fail: error [LOWERING/059]: mutation requires an owned automatic source binding |
| r15_struct_ref_scalar_field_write | `s.T = ...` through `Burn&` | **ok** | vgpr=10 scratch=0B | import-fail: error [LOWERING/059]: mutation requires an owned automatic source binding |
| r16_struct_ref_array_member_write | `s.xn[n] = ...` through `Burn&` | ok | vgpr=10 scratch=128B | same |
| r17_local_ref_array_elem | `double& r = a[2];` | ok | vgpr=4 scratch=0B | same |
| r18_value_return_max | `max(const T&, const T&)` returning `T` by value | ok | vgpr=6 scratch=0B | same |
| r19_struct_ref_member_accessor_write | `vs.y(n) += ...` through `Dvode&` | **ok** | vgpr=10 scratch=128B | import-fail: error [LOWERING/059]: mutation requires an owned automatic source binding |
| r20_const_ref_from_deref | `const Burn&` bound to `*p` | ok | vgpr=10 scratch=0B | same |
| r21_whole_struct_assign_through_ref | `r = rate_t{...}` through `rate_t&` | ok | vgpr=6 scratch=0B | same |
| r22_const_ref_to_kernel_arg_struct | `const BurnParams&` bound to `*params` (kernel argument) | ok | vgpr=12 scratch=0B | same |
| r23_const_struct_ref_member_accessor_read | `vs.y(n)` through `const Dvode&` | **ok** | vgpr=8 scratch=248B | import-fail: error [LOWERING/059]: mutation requires an owned automatic source binding |
| r24_const_struct_ref_member_raw_array | `vs.y.arr[n]` through `const Dvode&` | ok | vgpr=8 scratch=248B | same |
