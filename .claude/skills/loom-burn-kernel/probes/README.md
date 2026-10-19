# Loom frontend probe kit

`src/` holds one-construct probes (root `probe::k`) for the constructs a
Microphysics burn uses; `probe.sh` imports and compiles each with one frozen
tool set on the GPU host and writes `probe-summary.txt`. See the header of
`probe.sh` for how to ship and run it. Add probes for anything a port's
inventory finds that is not covered here, and keep new probes in this kit.

## Results at hrx-system 7cf21edd85 (gfx90a, run 2026-10-18)

`7cf21edd85` is a merge commit on the author's fork branch `react-kernels`; it changes nothing outside `experimental/`, so its compiler sources equal public `ROCm/hrx-system` main `4031286f0b` (the commit to check out in a fresh clone). Tools built from either should give the same probe results; their md5s differ with the build environment, so record yours.

Tools: the adv-mi210-12 campaign's frozen build (`loom-import-cxx` md5
`5f6207ac…`, `loom-compile` md5 `4aafb6f4…`). Raw summary and logs are kept
in the first port's (untracked, local) agent-notes directory
(`probes/results-7cf21edd85/`); they are not shipped.
All 70 probes run in a few seconds. Rerun the kit whenever the pinned commit changes and
update this table.

What these results decide:
- **References:** scalar `double&`/`int&` parameters and scalar member slices
  (`f(s.T)`) import, in `static` *and* non-static `inline` functions (r01,
  r02, r11, r13). A reference to a whole struct does not, even `static` and
  `const` (r03, r05: "require an automatic SSA object"), and neither do
  references to buffer elements, arrays, overlapping arguments, local
  reference bindings or reference returns. So `burn_t&`, `dvode_t&` and
  `rhs_state_t&` become pointers; scalar out-parameters can stay.
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
  at the importer's default `--std=c++26` (P2864R2 made it ill-formed). The
  pinned parser also rejects it under `--std=c++20`/`c++23`, where it is
  valid but deprecated; reported fixed upstream as `0d27041f`, a commit not found in public
  hrx-system history as of 2026-10-18 (the draft issue lives in the original
  author's local notes); re-probe d10. Microphysics'
  aprox13 configuration doesn't use the construct (no
  `-Wdeprecated-enum-float-conversion` warnings); if a translation does, cast
  to `int`.
- **Pointers stored in structs** can't be read through a pointer to the
  struct (c03, c06): pass a tables struct **by value** (c04), or use one flat
  buffer with compile-time offsets (c05). Kernel arguments must be the table
  pointers themselves, not a pointer to a struct holding them. (c03-c06 were
  added 2026-10-18 after the EOS smoke import failed on `const EosTables*`.)
- **Records** must be plain aggregates: no constructors with arguments (o07),
  no bases (o06). Structured bindings fail (o08). POD duals with free
  functions work (o04, o05).
- The simple retry loop imports without `__builtin_assume` (f07), and the
  isolated triangular loop compiles (f06), but real stage and LU loops still
  needed `MPL_ASSUME` (`references/loom-frontend-rules.md`, SUBRANGE row).

| Probe | Construct | Result | Detail |
|---|---|---|---|
| c00_control | baseline | ok | vgpr=8 scratch=0B |
| c01_params_pod_ptr | `const Params*` runtime-param POD | ok | vgpr=10 scratch=0B |
| c02_table_ptr_3d | helmholtz `f[j][i][9]` via pointer, clamp, local `fi[36]` | ok | vgpr=15 scratch=288B |
| d01_ns_constexpr_array | `aion[]`, `zion[]`, `aion_inv[]` | import-fail | error [LOWERING/059]: constants require a numeric scalar type |
| d02_local_const_array | local `const double t[N]` | ok | vgpr=16 scratch=104B |
| d03_ns_mutable_global | runtime param global | import-fail | error [LOWERING/059]: global storage definitions require a global-storage projection |
| d04_uninit_local | `bool converged;` | import-fail | error [LOWERING/059]: locals require initializers |
| d05_struct_default_init_braces | `T x{}` | ok | vgpr=8 scratch=0B |
| d06_struct_copy_assign | `old{s}; s = old;` | ok | vgpr=4 scratch=0B |
| d07_limits_and_nan | `__DBL_MAX__`, NaN via bit_cast | ok | vgpr=9 scratch=0B |
| d08_large_struct_runtime_index | dvode_t-sized local struct | ok | vgpr=10 scratch=2872B |
| d09_static_constexpr_member_array | `static constexpr` member array | import-fail | error [LOWERING/059]: member-dependent initialization requires source object initialization semantics |
| d10_enum_u8_base | `enum : unsigned char` | import-fail | error [PARSE/036]: invalid operands of types 'double' and '::probe::Rates' to binary operator '+' |
| d11_short_and_bool_arrays | `short pivot[14]`, `bool` arrays | ok | vgpr=10 scratch=42B |
| d12_f64_to_int_conversion | `int(double)` | ok | vgpr=10 scratch=0B |
| f01_switch_runtime | runtime `switch` | import-fail | error [LOWERING/059]: unsupported statement: switch-statement |
| f02_switch_in_constexpr_fn | `rhs_data(rate)` folded | import-fail | error [LOWERING/059]: unsupported C++ type: ::probe::RhsT |
| f03_return_in_loop | `return` in `while(true)` | import-fail | error [LOWERING/059]: returns inside loops require a structured loop-exit projection |
| f04_for_ever | `for(;;)` | import-fail | error [LOWERING/059]: for loops require a condition |
| f05_vode_retry_shape | nested `while(true)`, result code | ok | vgpr=20 scratch=0B |
| f06_runtime_triangular_loops | `for k=NQ..1, j=k..NQ` | ok | vgpr=11 scratch=672B |
| f07_retry_loop_no_assume | retry loop with break/continue, no `__builtin_assume` | ok | vgpr=12 scratch=0B |
| f08_retry_loop_assume | same, with `__builtin_assume` | ok | vgpr=12 scratch=0B |
| f09_newton_for_break | EOS Newton `for` + `break` | ok | vgpr=26 scratch=0B |
| f10_early_return_outside_loop | sneut5 early return | ok | vgpr=8 scratch=0B |
| l01_lambda_constexpr_for | lambdas | import-fail | error [LOWERING/059]: rvalue-reference parameters require consuming object identity |
| l02_template_struct_unroll | template-struct recursion | ok | vgpr=30 scratch=0B |
| m_f64_absf | f64 fabs | ok | vgpr=4 scratch=0B \| |
| m_f64_cbrtf | f64 cbrt | compile-fail | error [TARGET/001]: has no target-low contract |
| m_f64_copysignf | f64 copysign | ok | vgpr=6 scratch=0B \| |
| m_f64_cosf | f64 cos | compile-fail | error [LOWERING/036]: legalize-math policy amdgpu-math rejected scalar.cosf for cosf in scalar lanes of f64 under math.element.f32 |
| m_f64_expf | f64 exp | compile-fail | error [LOWERING/036]: legalize-math policy amdgpu-math rejected scalar.expf for expf in scalar lanes of f64 under math.element.f32 |
| m_f64_floorf | f64 floor | ok | vgpr=4 scratch=0B \| v_floor_f64:1 |
| m_f64_fmaf | f64 fma | ok | vgpr=8 scratch=0B \| v_fma_f64:1 |
| m_f64_log10f | f64 log10 | compile-fail | error [TARGET/001]: has no target-low contract |
| m_f64_logf | f64 log | compile-fail | error [LOWERING/036]: legalize-math policy amdgpu-math rejected scalar.logf for logf in scalar lanes of f64 under math.element.f32 |
| m_f64_math_h_recipes | existing exp/log/cbrt recipes | ok | vgpr=52 scratch=0B \| v_add_f64:49 v_fma_f64:36 v_mul_f64:19 v_div_scale_f64:4 v_cmp_lt_f64:3 v_rcp_f64:2 |
| m_f64_powf | f64 pow (non-integer exponent) | compile-fail | error [LOWERING/036]: legalize-math policy amdgpu-math rejected scalar.powf for powf in scalar lanes of f64 under math.element.f32 |
| m_f64_roundevenf | f64 rint (roundeven) | ok | vgpr=4 scratch=0B \| v_rndne_f64:1 |
| m_f64_sinf | f64 sin | compile-fail | error [LOWERING/036]: legalize-math policy amdgpu-math rejected scalar.sinf for sinf in scalar lanes of f64 under math.element.f32 |
| m_f64_sqrtf | f64 sqrt | ok | vgpr=12 scratch=0B \| v_fma_f64:7 v_mul_f64:2 v_ldexp_f64:2 v_rsq_f64:1 v_cmp_lt_f64:1 v_cmp_class_f64:1 |
| m_f64_tanhf | f64 tanh | compile-fail | error [LOWERING/036]: legalize-math policy amdgpu-math rejected scalar.tanhf for tanhf in scalar lanes of f64 under math.element.f32 |
| m_f64_truncf | f64 trunc | ok | vgpr=4 scratch=0B \| v_trunc_f64:1 |
| m_fp_contract | `a * b + c` under `--fp-model=precise` | ok | vgpr=8 scratch=0B \| v_mul_f64:1 v_add_f64:1 |
| m_powi_template | `powi<N>` recursion | ok | vgpr=8 scratch=0B \| v_mul_f64:5 v_add_f64:1 |
| o01_member_fn_call | `a.add(x)` | import-fail | error [LOWERING/059]: void call must resolve to a function symbol |
| o02_operator_call_index | `Array1D::operator()` | import-fail | error [LOWERING/059]: address-of requires an existing storage-backed element |
| o03_dual_operator_overload | `operator*`/`operator+` on Dual | import-fail | error [LOWERING/059]: overloaded binary operation requires a Loom operation |
| o04_dual_free_functions | POD Dual / Dual3 + free fns, by value | ok | vgpr=10 scratch=0B |
| o05_dual3_array_free_functions | POD Dual / Dual3 + free fns, by value | ok | vgpr=12 scratch=96B |
| o06_base_class | `Grad : Base` | import-fail | error [LOWERING/059]: record values require a complete aggregate without unions, bases or nontrivial lifecycle operations |
| o07_ctor_with_args | `D x(1.0)` | import-fail | error [LOWERING/059]: record values require a complete aggregate without unions, bases or nontrivial lifecycle operations |
| o08_structured_binding | `auto [f, r] = term(x)` | import-fail | error [LOWERING/059]: unsupported local declaration |
| r01_static_scalar_ref_out | `static f(double&)` on locals | ok | vgpr=8 scratch=0B |
| r02_inline_external_ref | `inline` (AMREX_INLINE) non-static fn with scalar `double&` | ok | vgpr=6 scratch=0B |
| r03_template_struct_ref | `template<class BurnT> inline f(BurnT&)` | import-fail | error [LOWERING/059]: runtime value-reference arguments require an automatic SSA object |
| r04_anon_ns_ref | anon-namespace fn with scalar `double&` | ok | vgpr=6 scratch=0B |
| r05_const_struct_ref_read | `static f(const Burn&)` | import-fail | error [LOWERING/059]: runtime value-reference arguments require an automatic SSA object |
| r06_ref_to_buffer_element | `f(out[i])` | import-fail | error [LOWERING/059]: runtime value-reference arguments require an automatic SSA object |
| r07_ref_return_accessor | `double& at(A&, int)` | import-fail | error [LOWERING/059]: address-of requires an existing storage-backed element |
| r08_overlapping_mutable_refs | `f(s, s.T)` | import-fail | error [LOWERING/059]: mutable value-reference arguments must select nonoverlapping automatic objects or member slices |
| r09_array_ref_param | `const double (&)[13]` | import-fail | error [LOWERING/059]: array reference parameters require addressable storage; pass a view or pointer instead |
| r10_local_ref_binding | `double& r = s.x;` | import-fail | error [LOWERING/059]: unsupported C++ type: double& |
| r11_ref_chain | ref forwarded through 2 static calls | ok | vgpr=10 scratch=0B |
| r12_const_ref_return_max | `amrex::max` returning `const T&` | import-fail | error [LOWERING/059]: runtime reference returns require object identity; return an admitted value instead |
| r13_struct_member_slice_ref | `bump(s.T)` | ok | vgpr=8 scratch=0B |
| c03_struct_of_ptrs_by_ptr | struct of table pointers passed by pointer (`t->a[i]`) | import-fail | error [LOWERING/059]: stored pointers and references require an object storage representation |
| c04_struct_of_ptrs_by_value | struct of table pointers passed by value (`t.a[i]`) | ok | |
| c05_flat_table_offsets | one table buffer, compile-time offsets | ok | |
| c06_ptr_member_of_kernel_arg_struct | kernel argument: pointer to a struct of pointers in device memory | import-fail | error [LOWERING/059]: stored pointers and references require an object storage representation |
