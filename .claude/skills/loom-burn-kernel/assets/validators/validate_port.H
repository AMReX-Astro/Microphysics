// Per-port hooks for validate_rates.H and validate_network.H (host only).
// Template: copy into the kernel directory next to them and fill in the
// lines marked "PORT:". Everything else in those two validators is
// network-neutral (it loops over Rates::NumRates, NumSpec and the screened
// pairs from network_data.H).
#ifndef VALIDATE_PORT_H
#define VALIDATE_PORT_H

#include <string>
#include <vector>

#include <extern_parameters.H>
#include <microphysics_autodiff.H>
#include <loom_validate.H>

#include "rates.H"

namespace loomval_port {

using mpl::Real;

// stock autodiff dual vs mpl::Dual, value and derivative bit for bit
inline void cmp_dual(LoomCheck& c, const char* what, const autodiff::dual& s, const mpl::Dual& l)
{
    const std::string w(what);
    loom_check_bits(c, (w + ".val").c_str(), autodiff::val(s), l.val);
    loom_check_bits(c, (w + ".grad").c_str(), autodiff::derivative(s), l.grad);
}

// Runtime parameters that select rate or network branches (network_rp::...).
// The validators run every check once per variant: apply_param_variant(v)
// for v = 0 .. num_param_variants - 1, then restore_params(). Screening's
// enable_debye_huckel_skip is toggled by the validators themselves.
// PORT: e.g. for network.use_c12ag_deboer17:
//   constexpr int num_param_variants = 2;
//   inline bool saved_deboer = false;
//   inline void save_params() { saved_deboer = network_rp::use_c12ag_deboer17; }
//   inline void apply_param_variant(int v) { network_rp::use_c12ag_deboer17 = (v == 1); }
//   inline void restore_params() { network_rp::use_c12ag_deboer17 = saved_deboer; }
constexpr int num_param_variants = 1;
inline void save_params() {}
inline void apply_param_variant(int v) { (void) v; }
inline void restore_params() {}

// Temperatures where a rate switches branch (each is bracketed by +-64 ulps
// by validate_rates.H). PORT: one entry per branch threshold of the
// translated rates, e.g. rate_triplealf's t9 = 0.08: t.push_back(0.08e9);
inline void branch_temperatures(std::vector<Real>& t) { (void) t; }

// Direct calls of the translated rate functions at several densities, for
// rates whose stock function takes a density argument that
// evaluate_analytical_rate fixes to 1. Called once per temperature and
// parameter variant with the stock and Loom tfactors (Real and dual).
// PORT: for each rate_<name>(tf, den, fr, rr) the network uses, e.g.
//   for (Real den : {1.0, 1.0e-3, 3.7e4, 1.0e8, 2.5e11}) {
//       Real fr = 0.0, rr = 0.0;
//       rate_triplealf(st, den, fr, rr);
//       const mpl::rate_t<Real> l = mpl::rate_triplealf_r(lt, den);
//       loom_check_bits(c, "triplealf_r.fr", fr, l.fr);
//       loom_check_bits(c, "triplealf_r.rr", rr, l.rr);
//       ... the same for the dual versions with cmp_dual() ...
//   }
template <typename StockTf, typename StockTfD>
inline void direct_rate_checks(LoomCheck& c, const StockTf& st, const StockTfD& std_,
                               const mpl::tf_t<Real>& lt, const mpl::tf_t<mpl::Dual>& ltd,
                               const mpl::RatesParams& rp)
{
    (void) c; (void) st; (void) std_; (void) lt; (void) ltd; (void) rp;
}

} // namespace loomval_port

#endif
