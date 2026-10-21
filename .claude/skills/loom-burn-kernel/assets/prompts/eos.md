# Prompt: equation of state (eos.H) — only for an EOS other than helmholtz

Prepend `_common.md`. Module name for builds: `eos`. (helmholtz is done:
`assets/modules/eos.H`, network data from the generated `network_data.H`.)

Scope: everything the burn reaches in `EOS/{EOS}/`, `EOS/eos_composition.H`,
`interfaces/eos.H` (eos(), reset_inputs and reset_* helpers: the runtime
switch becomes an if/else chain), `interfaces/eos_override.H`, and the state
types/input modes the integrator actually instantiates (check
`integration/integrator_setup_*.H`, `integrator_type*.H`). Tables become a
`struct EosTables` of `const Real*` passed BY VALUE with index macros matching
the stock C layout; runtime parameters an `EosParams` of scalars. Keep every
runtime-gated branch; list them in the report. Validate every eos_t field the
burn reads, all input modes the burn uses, beyond the table edges, and each
step you split out. Start from `assets/validators/validate_eos.H`.
