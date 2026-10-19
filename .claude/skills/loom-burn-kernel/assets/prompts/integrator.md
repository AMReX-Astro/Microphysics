# Prompt: ODE integrator + burner — only for an integrator other than VODE

Prepend `_common.md`. Module name for builds: `integrator`. (VODE + Strang
burner are done: `assets/modules/vode.H`, `burner.H`.)

Scope: `interfaces/burner.H`, `integration/integrator.H` (incl. burn retry),
`integration/{INTEGRATOR}/`, the Strang setup/cleanup/type/rhs glue in
`integration/`, `normalize_abundances_burn`, and anything in
`integration/utils/` the path really reaches (check: VODE computes its
numerical Jacobian inline, so `numerical_jacobian.H` was dead code for it).
Integrator state becomes a POD struct of flat, source-indexed arrays passed by
pointer; every `return` inside a loop becomes a result code + `break` with each
exit's side effects preserved; `MPL_ASSUME` on every loop over the runtime
method order and every loop with a `break`; read order/size fields into
`const int` locals before indexing. Entry point:
`burner(burn_t*, Real dt, const BurnParams*, EosTables)`.

Validate with `assets/validators/validate_integrator.H`: whole burns on the
test_react grids (n_cell 16 and 32) reproducing `do_react` exactly, every
`burn_t` field, plus hard zones and variants that force retries, failures and
each integrator switch; run with `integrator.use_jacobian_caching=0`.
