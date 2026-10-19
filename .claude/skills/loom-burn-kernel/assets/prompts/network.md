# Prompt: per-network glue (network.H)

Prepend `_common.md`. Module name for builds: `network`.

The generic `rhs.H` machinery is already translated:
`<skill>/assets/modules/network_rhs.H` (copy it into
the kernel directory unchanged). Its header comment lists the hooks the
per-network `network.H` must define BEFORE including it: `struct NetParams`,
`RatesParams net_rates_params(const NetParams*)`,
`net_analytical_rate_{r,d}<rate>` (the network's `evaluate_analytical_rate`),
`postprocess_rate_{r,d}<rate>` and `ener_gener_rate<spec>`, all taking
`(const burn_t*, const Real* y, ...)` because several networks read rho, y,
eta, y_e there.

Translate `networks/{NETWORK}/actual_network.H`'s `postprocess_rate` and
`ener_gener_rate` (both instantiations where templated) and the NetParams
fields for every runtime parameter the network path reads (`network.use_tables`
must be 0: not ported; the launcher rejects it).

Validate with the shipped `validators/validate_network.H` (network-neutral:
stock `RHS::rhs` and `RHS::jac` vs `actual_rhs`/`actual_jac` on every
element, plus `dgefa`/`dgesl`, over states built with the stock EOS,
T 1e4..1e16, rho 1e2..1e12, compositions built from NumSpec, random states).
Its per-network inputs are the `PORT:` hooks in `validators/validate_port.H`:
every network switch as a parameter variant (shared with the rates
validator). Add checks of your own only for paths the template doesn't
reach.
