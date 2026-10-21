# Prompt: extend dual.H — only for autodiff shapes it doesn't cover

Prepend `_common.md`. Module name for builds: `dual`.

`assets/kernel/dual.H` reproduces `autodiff::dual` and `dual_array<1,3>`
bit for bit for the expression shapes listed in its header (validated by
`assets/harness/validate_dual.H`). A new network may use shapes it doesn't
document (number^dual, dual^dual in use, NumberDualDualMulOp, new math
functions, a wider dual_array). For each: read the evaluation path in
`util/autodiff/autodiff/forward/dual/dual.hpp` under
`AUTODIFF_STRICT_ASSOCIATIVITY` (util/microphysics_autodiff.H), add POD free
functions that copy its formulas and order exactly, extend the header guide
with worked examples from the stock code, and extend `validate_dual.H` with
primitive tests and verbatim stock expressions. Never derive the calculus
yourself. Rerun the whole dual validator: 0 mismatches.
