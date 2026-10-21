# Seven dwarfs: burn configurations for exercising Loom

Seven `unit_test/test_react` configurations chosen together to push the Loom
compiler through as many different code patterns as Microphysics offers:
network size and style, integrator control flow, linear algebra, table
lookups, screening and neutrino math, EOS, and unit system. Use them as the
test set for a new Loom pin or a compiler change, roughly in the bring-up
order at the end.

They don't have to be physically standard configurations, because every
dwarf is checked against stock Microphysics built with the same settings.
They do have to succeed in stock test_react, which aborts with "Integration
failed" if any zone fails. All seven were verified that way on 2026-10-22
(CPU, OpenMP, gcc 14, 16³ = 4096 zones, moth).

## The seven

| # | Network | Integrator | Screening | Neutrinos | EOS | Inputs | Max ρ(J)·tmax | Max RHS calls per zone |
|---|---|---|---|---|---|---|---|---|
| 1 | aprox13 | VODE | screen5 | sneut5 | helmholtz | `inputs_aprox13` | 1.8×10⁵ | 1,163 |
| 2 | triple_alpha_plus_cago | VODE | null | off | gamma_law | `inputs_3alpha` | 18 | 2,183 |
| 3 | iso7 | QSS | chabrier1998 | kipp | helmholtz | `inputs_iso7` | 121 | 95,843 |
| 4 | sn160 | VODE | screen5 | off | helmholtz | `inputs_sn160` + `unit_test.tmax=5.e-7` | ~10¹⁰ | 4,132 |
| 5 | ase | Rosenbrock | chugunov2007 | sneut5 | helmholtz | custom (below) | 6.4×10⁹ | 4,364 |
| 6 | aprox19 | RKC | screen5 | sneut5 | helmholtz | `inputs_aprox19` | 3.1×10³ | 11,407 |
| 7 | metal_chem | Rosenbrock | — | — | metal_chem | custom (below) | 3×10⁷ | 90,529 |

ρ(J)·tmax is the largest Jacobian spectral radius over the grid times the
burn time (see "Checking stiffness"). Build each from `unit_test/test_react`
with, for example:

```bash
make -j USE_OMP=TRUE NETWORK_DIR=aprox13 EOS_DIR=helmholtz INTEGRATOR_DIR=VODE \
     SCREEN_METHOD=screen5 USE_NEUTRINOS=TRUE NEUTRINO_METHOD=sneut5
```

(For correctness reference builds add `XTRA_CXXFLAGS=-ffp-contract=off`, as
for any port. Use one source tree per configuration: test_react's
`tmp_build_dir` isn't keyed by network, so builds in the same tree collide.)

## Why each one

**1. aprox13 / VODE / screen5 / sneut5: the reference.** This is the
production default in Castro and the configuration this skill was first
ported with (neutrinos were off for that port). It covers the
`networks/rhs.H` framework: rates driven by `constexpr` data (`rhs_data`,
intermediate rates), compile-time unrolling, and an LU with a compile-time
sparsity pattern and no pivoting. VODE adds variable order, multi-exit loops
and runtime-bounded indices (the SUBRANGE/026 and branch-lowering trouble so
far), screen5 adds autodiff dual numbers, and helmholtz adds 3D table
lookups. Turning on sneut5 adds its piecewise branches and `sin`/`cos`,
which have no f64 recipe yet, so this dwarf is the one that forces those.

**2. triple_alpha_plus_cago / VODE / null / gamma_law: the canary.** The
smallest state (5 equations), no screening and no EOS tables, so a failure
points at the frontend or integrator, not at physics volume. Its compiles take
minutes, so it is the first thing to run after any compiler bump. It shares
VODE with dwarf 1 but differs in every physics module, which separates
integrator problems from network problems.

**3. iso7 / QSS / chabrier1998 / kipp: explicit, no linear algebra.** QSS
splits each derivative into production and loss terms and adapts the step
with predictor-corrector iterations, with no Jacobian or linear solve. It
needs an `rhs.H` network (it calls `RHS::rhs`). chabrier1998 brings a
different mix of transcendental functions than screen5, and kipp is the
other neutrino model, written as autodiff templates.

**4. sn160 / VODE / screen5: scale.** 160 species and 1,626 rates (1,458
ReacLib, 168 tabulated weak rates), with partition functions and a
67,947-line generated `actual_rhs.H`: a huge straight-line program for
instruction scheduling, register allocation and compile time. Its 161×161
dense Jacobian alone is 207 KB per lane, over gfx90a's 131,056 B private
limit, so this dwarf can't run at one GPU lane per zone. It is the test case
for a team-per-zone or shared-memory layout.

**5. ase / Rosenbrock / chugunov2007 / sneut5: pivoted LU and table rates.**
Rosenbrock adds a pivoted dense LU, where the pivot row comes from the data
and becomes a private-array index, and a Butcher tableau of constant arrays.
ase adds tabulated weak rates (search a table index, then 2D interpolation in
ρY_e and T, all data-dependent loads), derived reverse rates and partition
functions. chugunov2007 adds more `pow`/`exp`/`sqrt` work in screening.

**6. aprox19 / RKC / screen5 / sneut5: explicit Chebyshev stepping.** RKC
needs no Jacobian or linear solve, but it estimates the spectral radius by
power iteration (a data-dependent loop) and chooses a runtime stage count.
aprox19 is an `rhs.H` network with tabulated weak rates (Langanke), and it is
the one candidate that is mildly stiff on its standard grid.

**7. metal_chem / Rosenbrock: astrochemistry.** Gas-phase chemistry instead
of nuclear burning: 35 species in number-density mode, its own EOS, cosmic-ray
terms, a dust opacity table (`Semenov_PlanckOpacity.dat`, which must be in the
run directory), rates dense in `exp(a*log T)` forms, and a very wide dynamic
range of values. It is the dwarf closest to GPUAstroChem/JAFF kernels.
Rosenbrock rather than VODE because VODE was much slower here (131 s for
Rosenbrock; VODE had not finished after 12 minutes).

**Coverage.** Integrators: VODE, Rosenbrock, QSS, RKC. Screening: screen5,
chabrier1998, chugunov2007, null. Neutrinos: sneut5, kipp, off. EOS:
helmholtz, gamma_law, metal_chem. Units: mass fractions and number
densities. Networks: `rhs.H` (1, 2, 3, 6), pynucastro (4, 5), chemistry (7).

## Constraints behind the choices

- **NSE can't be a test_react dwarf.** `integration/Make.package` stops with
  "NSE with Strang integration is not supported" for `USE_NSE_TABLE` and
  `USE_NSE_NET`, and test_react only supports Strang splitting (`do_react`
  fills a Strang `burn_t`; SDC has its own tests). `USE_NSE_NET` would also
  require null or chabrier1998 screening. Simplified SDC is ruled out for the
  same reason.
- **BackwardEuler is not used.** It is not used in production, and on
  triple_alpha it ran out of steps (`IERR_TOO_MANY_STEPS`) in the hot corner
  of `inputs_3alpha` even with a 10× step cap.
- **RKC only where the problem is mildly stiff.** Check ρ(J)·tmax first (next
  section). sn160 is never mildly stiff on its grid (2.8×10⁷ even in the
  coldest zone); an sn160/RKC run had not finished after 20 minutes.
- **QSS never near NSE.** iso7/QSS passes on `inputs_iso7` (T ≤ 5×10⁸) but
  fails on the same grid extended to T = 5×10⁹.
- **sn160 needs the shorter burn.** With `tmax` = 10⁻⁶ s, two zones near
  ρ ≈ 1.8×10⁷ (T₀ = 2.2×10⁷ and 5.4×10⁷ K) stall at t ≈ 7.6–8.2×10⁻⁷ s with
  either the default 150,000-step cap or 10⁶. With 5×10⁻⁷ s every zone passes
  (192 s on 24 cores). Rosenbrock on sn160 had not finished after 15 minutes.

## Checking stiffness

For a test_react grid, evaluate the Jacobian at each zone's initial state
(after `eos(eos_input_rt, ...)`) and take its eigenvalues λ. The deciding
number for an explicit or Chebyshev method is ρ(J)·tmax, where ρ(J) = max|λ|:
RKC's stage count for a step h is about √(ρ(J)·h / 0.653) + 1. Treat
ρ(J)·tmax ≲ 10⁴ (about 125 stages over the whole burn) as mildly stiff. The
ratio max|Re λ| / min|Re λ| isn't useful here: near-zero modes put it near
10¹¹–10¹² for every network.

Eigenvalues don't depend on whether species are written as X or Y, because
that change is a diagonal similarity transform, so the analytic `actual_jac`
(in Y) works directly. Networks without an analytic Jacobian (metal_chem's is
empty; its inputs use `integrator.jacobian = 2`) need finite differences of
`actual_rhs`; `numerical_jac` doesn't compile with metal_chem because its
`actual_rhs` takes a non-const `burn_t`.

Measured on each network's test_react grid (8³ zones):

| Network (inputs, tmax) | Median ρ(J)·tmax | Max ρ(J)·tmax | Mildly stiff? |
|---|---|---|---|
| triple_alpha (`inputs_3alpha`, 10⁻³) | 3×10⁻⁸ | 18 | yes |
| iso7 (`inputs_iso7`, 10⁻⁶) | 5×10⁻¹¹ | 121 | yes |
| aprox19 (`inputs_aprox19`, 10⁻⁹) | 0.07 | 3.1×10³ | yes |
| aprox13 (`inputs_aprox13`, 10⁻⁵) | 1.4×10⁻⁴ | 1.8×10⁵ | no (hot corner) |
| nova (`inputs_nova`, 10⁻⁴) | 18 | 2.0×10⁵ | no |
| subch_simple (aprox13 grid, 10⁻⁵) | 54 | 2.4×10⁸ | no |
| ase (custom, 10⁻⁵) | 5.6×10⁴ | 6.4×10⁹ | no |
| CNO_He_burn (aprox13 grid, H/He/C, 10⁻⁵) | 2.2×10⁵ | 7.1×10¹⁰ | no |
| sn160 (`inputs_sn160`, 10⁻⁶) | 1.3×10⁹ | 3.9×10¹⁰ | no, anywhere |
| metal_chem (custom, 3.15×10¹³) | 4.4×10⁶ | 3.1×10⁷ | no |

## Custom inputs

`inputs_ase` (crosses ase's NSE temperature range, but NSE stays off):

```
n_cell = 16
prefix = react_ase_
unit_test.small_dens = 1.0e0
unit_test.dens_min = 1.e6
unit_test.dens_max = 1.e9
unit_test.temp_min = 1.e8
unit_test.temp_max = 7.e9
unit_test.tmax = 1.e-5
unit_test.primary_species_1 = "helium-4"
unit_test.primary_species_2 = "carbon-12"
unit_test.primary_species_3 = "oxygen-16"
```

`inputs_metal_chem` (ISM-like densities and temperatures, a 1 Myr burn,
integrator settings from `unit_test/burn_cell_metal_chem`; species names must
be the lowercase `spec_names`, e.g. `h`, not `H`):

```
n_cell = 16
prefix = react_metal_chem_
unit_test.small_dens = 1.e-60
unit_test.small_temp = 1.e0
unit_test.dens_min = 1.e-25
unit_test.dens_max = 1.e-19
unit_test.temp_min = 1.e1
unit_test.temp_max = 1.e4
unit_test.tmax = 3.15e13
unit_test.primary_species_1 = "h"
unit_test.primary_species_2 = "he"
integrator.use_number_densities = 1
integrator.subtract_internal_energy = 0
integrator.do_species_clip = 0
integrator.SMALL_X_SAFE = 1e-60
integrator.use_jacobian_caching = 1
integrator.X_reject_buffer = 1e100
integrator.jacobian = 2
integrator.renormalize_abundances = 0
integrator.rtol_spec = 1.0e-4
integrator.atol_spec = 1.0e-10
integrator.rtol_enuc = 1.0e-4
integrator.atol_enuc = 1.0e-10
network.redshift = 0.0
network.metallicity = 1
network.dust2gas_ratio = 1
network.small_x = 1e-60
```

## Bring-up order

2 → 1 (turning on sneut5) → 3 → 7 → 5 → 6 → 4. sn160 goes last because it
likely needs a different kernel layout. The Loom features each one forces
first: sin/cos f64 recipes (1), pivoted LU and the weak-rate table search (5),
very large straight-line code or a team-per-zone layout (4).
