# Common preamble for every translation task

Paste this at the top of every module prompt (fill in the `{...}` fields,
and replace `<skill>` with the skill's directory, e.g.
`.claude/skills/loom-burn-kernel` or `.agents/skills/loom-burn-kernel`).
It collects the rules and the corrections that had to be sent to agents
mid-flight during the first helmholtz + VODE port (2026-10-18); sending
them up front is what saves the time. With sub-agents, give each module
prompt to its own agent; working alone (or in a tool without sub-agents),
treat each prompt as a task spec, do them in order, and skip the "report"
parts that only matter for hand-off.

---

You are translating part of the AMReX-Astro Microphysics burn
({NETWORK} + {EOS} + {INTEGRATOR}, {OPTIONS}) into the restricted C++ subset
accepted by the Loom compiler. Repo: {REPO}, branch `{BRANCH}` (checked out).
Do not commit or switch branches. Kernel directory: `{KERNEL_DIR}`.

Read first:
- `{KERNEL_DIR}/README.md`: translation rules, file ownership, build and
  validation commands. Edit only the files you own; report anything you need
  changed elsewhere.
- `<skill>/probes/README.md`: what Loom accepts.
- `<skill>/references/loom-frontend-rules.md` and
  `references/microphysics-translation.md`.
- `{KERNEL_DIR}/dual.H` header comment: the only allowed way to write
  autodiff arithmetic.
- Reachability and specialized source:
  `{INVENTORY_DIR}/` (`summary.md`, `hits.tsv`, `specialized/*.cpp` from
  C++ Insights). `constexpr_for` loops appear fully unrolled there: take loop
  structure from the original source and use the specialized form for
  resolved types, branches and constants.

Hard rules (each one cost a failed Loom import or compile at least once):
1. Same computation, same order, statement for statement. Follow the HIP/GPU
   branch of every `#ifdef` (`AMREX_USE_GPU` defined, no host printing, no
   `ALLOW_JACOBIAN_CACHING`).
2. A struct that a function writes through is passed by pointer (`burn_t*`,
   `dvode_t*`): writing a scalar field or a member's accessor through a
   mutable struct reference fails import. A struct with `Array1D`/`Array2D`
   members always goes by pointer (`const dvode_t*` when read-only): calling
   a member's accessor through any struct reference fails, `const` included.
   Other `const` struct references, scalar
   `Real&`/`int&` parameters, array references, references to array
   elements, member function calls and `Array1D`/`Array2D`-style 1-based
   `operator()` accessors all import: keep them and the source's index
   expressions. Small structs are returned by value; `min`/`max` return by
   value, not `const T&`.
3. **A struct that holds pointers (tables) is passed BY VALUE, never through a
   pointer** (`const Tables* t; t->f[i]` fails import). Kernel arguments are
   the table pointers themselves.
4. No lambdas/`constexpr_for` (template-struct recursion or plain loops), no
   `constexpr` locals of record type (bind each field to a `constexpr`
   scalar), no namespace-scope arrays, no `switch`, no `return`
   inside loops (flags + `break`), no range-for, no structured bindings, no
   constructors with arguments, no uninitialized struct locals, no
   `std::`, no stock namespaces (`amrex::`, `autodiff::`, `admath::`, `gcem::`,
   `RHS::`, `*_rp::` ...): validators include the stock headers, so a
   leftover stock call still compiles and validates on the host and only the
   lint or Loom catches it.
5. **Bounds assumptions**: write counted loops with compile-time limits with
   a strict `<` bound (`j < neqs + 1`), which Loom proves in range; every
   other private-array index whose range the compiler can't see (runtime
   bounds from struct fields) needs an `MPL_ASSUME` (otherwise `loom-compile`
   fails with SUBRANGE/026). Placement and form rules:
   `references/loom-frontend-rules.md`.
6. **Loops with several exits** (multiple `break`s, a `continue`, the
   `break`s that replace stock `return`s) can fail the AMDGPU backend with
   TARGET/034 (`masked_region_exits_by_cfg`, `masked_region_merge_vgpr_values`) when
   one sits inside another (probes f21-f23). Use the single-exit form for
   every inner loop with more than one exit (the outermost loop may keep the
   stock `break`s, as `vode.H`'s `dvode` does): `bool done = false; while (!done) {...}`,
   each exit sets `done` after its side effects, later code in the body is
   guarded by `if (!done)`, and a `continue` becomes a guard (keep its exact
   predicate under `!(...)` so NaN behaviour is unchanged).
7. Math only through `mpl::math`; never substitute identities
   (`pow = exp(y*log x)`); integer powers through `math::powi<N>`.
8. Runtime parameters through a params struct of scalars (`const P*`), fields
   named `<namespace>_<name>`. The host `fill_*_params()` functions live in
   `unit_test/test_react/loom_params.H` (from
   `<skill>/assets/harness/loom_params.H`), except the
   EOS ones, which live in `unit_test/test_react/loom_eos_params.H`; never in a
   validator; report new fields so they are added there.

Validation (zero mismatches, bit for bit):
- `{KERNEL_DIR}/validate_<module>.H` with `long validate_<module>()`, using the
  helpers in `unit_test/test_react/loom_validate.H`. Start from the matching
  template in `<skill>/assets/validators/`.
- Compare against the stock functions over the full test_react range
  (rho 1e6-1e9, T 1e6-1e15, and beyond), the test_react compositions
  (`unit_test/react_util.H`) plus random ones, and every runtime switch the
  code branches on (toggle the stock global, restore it).
- Include a sensitivity check: perturb one operation, confirm mismatches,
  revert, rerun clean. Do it in a copy of the kernel dir, never in the real
  port (recipe: `references/test-react-integration.md`, "Sensitivity check").
- Put each per-sample comparison in an `__attribute__((noinline))` function
  (the templates do). Otherwise g++ -O3 can constant-fold a stock libm call
  (exp, cbrt, pow) for a loop-peeled constant input with MPFR, correctly
  rounded, while the Loom side calls libm at run time: false 1-ulp mismatches
  at one sample (`references/footguns.md`, "GCC constant-folds libm calls").
- Build in your OWN directories (other modules may build at the same time):
  `make -j4 {MAKE_VARS} XTRA_CXXFLAGS=-ffp-contract=off
   USE_LOOM_VALIDATE=TRUE LOOM_KERNEL_DIR=loom/<port> XTRA_CPPFLAGS=-DLOOM_VALIDATE_<MODULE>
   TMP_BUILD_DIR=tmp_loom_<module> EBASE=loomval_<module>`
  (`LOOM_KERNEL_DIR` is `{KERNEL_DIR}` relative to `unit_test/test_react`;
  `-ffp-contract=off` is required: SKILL.md Phase 0; rerun make once if a fresh TMP_BUILD_DIR fails with "No rule to make target
  extern_parameters.cpp"), then
  `./loomval_<module>3d.gnu.ex {INPUTS} unit_test.loom_validate=1 integrator.use_jacobian_caching=0`
  (default `COMP=gnu`, the same compiler as the stock CPU reference and CI;
  every native validation build uses it, `references/test-react-integration.md`).
- Lint your kernel headers (not validators):
  `ast-grep scan -c <skill>/scripts/loom-lint/sgconfig.yml <files>`.
  Read the output: "Error: Cannot parse rule" means the scan didn't run.
  0 errors required.

As a sub-agent: never run `loom-import-cxx`, `loom-compile` or any Loom tool,
and never ssh; whoever coordinates the port smoke-compiles your module on
the GPU host with `scripts/loom_smoke.sh` as soon as it validates (one Loom
compile at a time). Working alone, do that smoke compile yourself, on the
GPU host only, with `loom_smoke.sh` (`wait` before the next one). Never `make clean` or
delete build directories other than your own.

Final report (under ~500 words): what was translated, the API other modules
call (exact signatures, array layouts, index conventions), params/tables
fields, line counts, validation lines verbatim with sample counts, lint
result, any construct you used that `probes/README.md` doesn't cover, open
questions.
