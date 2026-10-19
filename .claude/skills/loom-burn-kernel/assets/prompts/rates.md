# Prompt: per-network rate functions (rates.H)

Prepend `_common.md`. Module name for builds: `rates`.

Scope: the rate functions this network's `evaluate_analytical_rate` reaches
(often in the shared `rates/` library; check
`networks/{NETWORK}/actual_network.H`), plus `interfaces/tfactors.H`
(`tf_t`, `get_tfactors`) if not already present, and `RatesParams` for every
runtime switch the rates read (e.g. a `network.use_<variant>` rate-choice flag).

Declare exactly the interface in `assets/NEW_NETWORK.md` ("The generic
modules include ..."): `RatesParams` (with the two `screening_*` fields
`screen.H` reads, plus one field per runtime switch the rates read),
`rate_t<T>`, `tf_t<T>`, `get_tfactors_{r,d}`,
`rate_<name>_{r,d}(const tf_t<T>*, const Real den[, const RatesParams*])`
returning `rate_t<T>` by value, and the 2-argument
`evaluate_analytical_rate_{r,d}<rate>(const tf_t<T>*, const RatesParams*)`
with den = 1 inside. Indices and constants come from the generated
`network_data.H` (`Rates::`, `NumRates`); don't hand-write them.

Every function exists twice: `_r` (number_t = Real, used by rhs()) mirroring
the plain stock arithmetic exactly (never reshape it the dual way: `a/b` stays
`a/b`), and `_d` (number_t = autodiff::dual, used by jac()) written only with
`dual.H` calls in the order its guide prescribes. Use separate names, not
overloads. Screen factors are generated, not computed: if stock folds
`calculate_screen_factor` with gcem at compile time, use the generated values;
if anything gcem-based runs at runtime, transcribe gcem statement for
statement (HIP runs gcem too).

Validate with the shipped `validators/validate_rates.H` (network-neutral:
tfactors, the `evaluate_analytical_rate` dispatch for every rate index in
both instantiations, every screened pair, over T 1e-3..1e3 GK and beyond).
Fill its hooks in `validators/validate_port.H`: every branch temperature of
the rates you translated (T9 cutoffs; the validator brackets each by ±64
ulps), each runtime switch as a parameter variant, and direct calls of each
`rate_<name>_{r,d}` at several densities.
