// POD dual numbers for the Loom burn kernel (see README.md, rule 10).
//
// Mirrors util/autodiff/autodiff/forward/dual/dual.hpp as Microphysics builds
// it (AUTODIFF_STRICT_ASSOCIATIVITY, util/microphysics_autodiff.H:18), bit for
// bit. Validated natively by validate_dual.H.
//
// ===========================================================================
// 1. Types on the reachable burn path (rhs.H network + helmholtz + VODE, no nu)
// ===========================================================================
//
//   autodiff::dual  (= Dual<double,double>)              -> mpl::Dual
//       get_tfactors, every rate_* the network reaches, postprocess_rate,
//       apply_density_scaling / apply_screening / apply_branching, and ALL of
//       screening: networks/rhs.H jac() builds rhs_state_t<autodiff::dual>,
//       so fill_plasma_state, actual_screen5, debye_huckel, actual_log_screen
//       and actual_screen run on autodiff::dual (d/dT only).
//       rhs() uses rhs_state_t<amrex::Real>: no duals there at all.
//
//   autodiff::dual_array<1,3> (= Dual<double, GradArray<double,1,3>>,
//   screening_dual_t)                                     -> mpl::Dual3
//       NOT reachable in this configuration. It only appears in the
//       non-template actual_screen/actual_log_screen overloads at
//       screening/screen.H:749,763 (called by unit_test/screen_derivs_check,
//       not by rhs.H) and in neutrinos/sneut5.H (USE_NEUTRINOS=FALSE).
//       Provided anyway (dual3_*) and validated, for configurations that
//       need it. grad keeps the source's 1-based indices (README rule 4):
//       grad[1..3] = (T, M1, M2); grad[0] is unused and always 0.
//
//   Expression nodes instantiated (C++ Insights, inventory-nonu/specialized):
//   AddOp, MulOp, NumberDualMulOp, NegOp, InvOp (only inside a/b), ExpOp,
//   LogOp, SqrtOp, CbrtOp, PowOp (dual^number only), plus amrex::Math::powi
//   on a dual (debye_huckel), admath::min/max, comparisons, operator+=/*=.
//   Not instantiated: NumberDualDualMulOp (never built under strict
//   associativity), AbsOp, Log10Op, number^dual, dual^dual. dual_abs,
//   dual_log10 and dual_pow are provided and validated anyway.
//
// ===========================================================================
// 2. API (namespace mpl; every function takes and returns by value)
// ===========================================================================
//
//   struct Dual  { Real val; Real grad; };
//   struct Dual3 { Real val; Real grad[4]; };       // grad[1..3]
//
//   Each dual3_* function below exists with the same meaning; it differs only
//   where GradArray's own operators order operands differently (handled
//   inside; translators need not care).
//
//   Construction / extraction
//     dual_const(n)         number n as a dual: val = n, grad = 0
//                           (`number_t x = n;`, `x = n;`, `return n;`,
//                           `rates.rr = 0.0;`, and every number operand that
//                           the rules below turn into a full dual)
//     dual_seed(x)          val = x, grad = 1   (`dual t = x; autodiff::seed(t);`)
//     dual3_seed(d, k)      d with grad[k] = 1  (autodiff::seed_array(t, m1, m2)
//                           is dual3_seed(t,1), dual3_seed(m1,2), dual3_seed(m2,3))
//     d.val                 autodiff::val(d), `static_cast<Real>(d)`
//     d.grad                autodiff::derivative(d) (Dual3: d.grad[k] = derivative(d)(k))
//     Dual d{};             default-constructed autodiff dual (0, 0)
//
//   Accumulator primitives: "s OP= <leaf>" exactly as autodiff's assign*:
//     dual_add(s, d)        s += dual d
//     dual_sub(s, d)        s -= dual d
//     dual_add_num(s, n)    s += number n  (grad untouched)
//     dual_add_ndm(s, n, d) s += n*d       (NumberDualMulExpr leaf, no temporary)
//     dual_mul(s, d)        s *= dual d    (grad = s.grad*d.val; grad += s.val*d.grad;
//                                           val = s.val*d.val)
//     dual_mul_num(s, n)    s *= number n  (val *= n; grad *= n) -- also n*d itself
//     dual_div(s, d)        s /= dual d    (val /= d.val; grad -= val*d.grad;
//                                           grad /= d.val)
//     dual_div_num(s, n)    s /= number n
//     dual_pow_num(s, n)    s = pow(s, number n)
//     dual_pow(s, d)        s = pow(s, dual d)  (also number^dual, see rules)
//   Unary (autodiff apply<Op>):
//     dual_neg, dual_inv (1/x as a lone InvExpr, NOT a/b), dual_exp, dual_log,
//     dual_log10, dual_sqrt, dual_cbrt, dual_abs
//   Others:
//     dual_powi<N>(x)       admath::powi<N>(x) == amrex::Math::powi on a dual
//     dual_min(x, y)        admath::min: (x.val <= y.val) ? x : y
//     dual_max(x, y)        admath::max: (x.val >= y.val) ? x : y
//                           (NOT mpl::math::min/max: other tie and NaN rules.)
//                           A number argument becomes dual_const(n) in place.
//                           This is the overload taken whenever an argument is
//                           a number, an expression, a temporary or a non-const
//                           dual -- every admath::min/max on the burn path.
//     dual_min_std/dual_max_std  only if BOTH arguments are `const` dual
//                           lvalues of the same type: C++ then picks amrex::min
//                           /max (= std::min/max: (y<x)?y:x, (x<y)?y:x).
//   Comparisons (`d > 1.6`, `a <= b`): compare .val with the plain operator.
//
// ===========================================================================
// 3. Translation guide: how autodiff evaluates an expression
// ===========================================================================
//
// Translate every dual-valued expression in two steps: (A) work out the
// expression TREE that autodiff's operator overloads build (it is not always
// the tree C++ precedence suggests), then (B) evaluate it with EVAL below.
// Only grouping and operand order affect the result: independent
// sub-expressions are pure values, so the order in which you compute them is
// free. Keep each number sub-expression (e.g. `rc * 1.35e-07_rt`,
// `5.0_rt/6.0_rt`, `zion[n-1]*y(n)`) as the same plain-double expression.
//
// (A) Tree building. C++ groups operators as usual (left to right for
//     + - * /). Each operator then rewrites its operands (number = double or
//     int, not a dual; dual = an existing Dual variable or a function's
//     returned Dual; expr = anything else that involves a dual):
//
//   -X     X number            -> the number -X (plain double negation)
//          X == -Y             -> Y
//          X == n*d (NDM)      -> NDM(-n, d)
//          otherwise           -> Neg(X)
//   L + R  L == -A and R == -B -> -(A + B)   (apply the -X rule to A + B)
//          otherwise           -> Add(L, R)  (no reordering)
//   L - R                      -> L + (-R)   (apply -R first, then +)
//   L * R  L == -A and R == -B -> A * B      (apply * again)
//          L number, R == -X   -> (-L) * X   (apply * again)
//          L number, R dual    -> NDM(L, R)  (NumberDualMulExpr)
//          otherwise           -> Mul(L, R)  (includes dual*number, number*expr)
//   L / R  R == Inv(X)         -> L * X
//          otherwise           -> L * Inv(R) (also for R a number!)
//   admath::f(X)               -> f(X)       f in exp, log, sqrt, cbrt, ...
//   admath::pow(L, R)          -> Pow(L, R)
//   autodiff::eval(X)          -> X evaluated to a dual now (EVAL(X))
//
//   Note what these rules imply:
//     a - 2.0          -> Add(a, -2.0)            (a number leaf, not "sub")
//     c - n*d          -> Add(c, NDM(-n, d))
//     -32.12*x         -> NDM(-32.12, x)          (-32.12 is a literal)
//     0.8*exp(..)      -> Mul(0.8, Exp(..))       (NOT NDM: exp(..) is no dual)
//     1.0/x            -> Mul(1.0, Inv(x))
//     x/0.025          -> Mul(x, Inv(0.025))
//     n*d1*d2          -> Mul(NDM(n, d1), d2)
//
// (B) Evaluation. EVAL(E) yields a Dual; ADD/SUB/MUL/DIV/POW(s, R) update the
//     accumulator s with operand R:
//
//   EVAL(number n)     = dual_const(n)
//   EVAL(dual d)       = d
//   EVAL(NDM(n, d))    = dual_mul_num(d, n)
//   EVAL(Add(L, R))    = ADD(EVAL(L), R)
//   EVAL(Mul(L, R))    = MUL(EVAL(L), R)
//   EVAL(Pow(L, R))    = POW(EVAL(L), R)
//   EVAL(Neg(X))       = dual_neg(EVAL(X))
//   EVAL(Inv(X))       = dual_inv(EVAL(X))        (a lone Inv; rare)
//   EVAL(f(X))         = dual_f(EVAL(X))
//
//   ADD(s, R): R number -> dual_add_num(s, R)
//              R dual   -> dual_add(s, R)
//              R == -X  -> SUB(s, X)
//              R == NDM(n, d) -> dual_add_ndm(s, n, d)
//              else     -> dual_add(s, EVAL(R))
//   SUB(s, R): R dual   -> dual_sub(s, R)
//              R == -X  -> ADD(s, X)
//              else     -> dual_sub(s, EVAL(R))
//              (SUB only arises from ADD(s, -X); X is never a number or an
//               NDM there because the -X rule already rewrote those.)
//   MUL(s, R): R number -> dual_mul_num(s, R)
//              R dual   -> dual_mul(s, R)
//              R == -X  -> dual_neg(MUL(s, X))
//              R == Inv(X) -> DIV(s, X)
//              else (NDM too!) -> dual_mul(s, EVAL(R))
//   DIV(s, R): R number -> dual_div_num(s, R)
//              R dual   -> dual_div(s, R)
//              R == -X  -> dual_neg(DIV(s, X))
//              R == Inv(X) -> MUL(s, X)
//              else     -> dual_div(s, EVAL(R))
//   POW(s, R): R number -> dual_pow_num(s, R)
//              else     -> dual_pow(s, EVAL(R))
//
// (C) Statements.
//   number_t x = E;  x = E;  return E;   -> x = EVAL(E)
//   x += E;   -> x = dual_add(x, EVAL(E))     EVAL(number n) = dual_const(n):
//   x -= E;   -> x = dual_sub(x, EVAL(E))     compound assignment ALWAYS uses a
//   x *= E;   -> x = dual_mul(x, EVAL(E))     full dual temporary, so
//   x /= E;   -> x = dual_div(x, EVAL(E))     `x *= 1.7` is dual_mul(x,
//                                             dual_const(1.7)), NOT dual_mul_num,
//                                             and `x += n` is NOT dual_add_num.
//                                             (They differ in signed zeros,
//                                             inf and NaN of the gradient.)
//   admath::min/max(A, B)  -> dual_min/dual_max(EVAL(A), EVAL(B))
//   admath::powi<N>(E)     -> dual_powi<N>(EVAL(E))
//   if (d > c)             -> if (d.val > c)
//
// (D) Worked examples from the stock code (Dual unless stated; `tf.x` is the
//     Dual field x of the translated tf struct).
//
//  interfaces/tfactors.H get_tfactors
//    tf.t9  = temp * 1.0e-9_rt;        t9  = dual_mul_num(temp, 1.0e-9);
//    tf.t92 = tf.t9*tf.t9;             t92 = dual_mul(t9, t9);
//    tf.t912 = admath::sqrt(tf.t9);    t912 = dual_sqrt(t9);
//    tf.t913 = admath::cbrt(tf.t9);    t913 = dual_cbrt(t9);
//    tf.t9i = 1.0e0_rt/tf.t9;          t9i = dual_div(dual_const(1.0), t9);
//    tf.lnt9 = admath::log(tf.t9);     lnt9 = dual_log(t9);
//
//  Analytic rate functions (shapes common in hand-written rate libraries;
//  c1, c2, ... stand for the source's numeric literals, which stay exactly
//  as written; Real locals such as q1, rc, den stay plain doubles; rev is
//  a Dual where it multiplies as one).
//
//  A fit with a polynomial prefactor, a 1/x and an exponential:
//    aa = 1.0e0_rt + c1*tf.t9i23;
//      -> Add(1.0, NDM(c1, t9i23))
//      aa = dual_add_ndm(dual_const(1.0), c1, t9i23);
//    bb = tf.t92*aa*aa;   bb = dual_mul(dual_mul(t92, aa), aa);   (left fold)
//    cc = admath::exp(-c2*tf.t9i13 - tf.t92*q1);        (q1 Real)
//      -> Exp(Add(NDM(-c2, t9i13), Neg(Mul(t92, q1))))
//      cc = dual_exp(dual_sub(dual_mul_num(t9i13, -c2),
//                             dual_mul_num(t92, q1)));
//    gg = c3 * tf.t9i32 * admath::exp(-c4*tf.t9i);
//      -> Mul(NDM(c3, t9i32), Exp(NDM(-c4, t9i)))
//      gg = dual_mul(dual_mul_num(t9i32, c3),
//                    dual_exp(dual_mul_num(t9i, -c4)));
//    zz = 1.0e0_rt/bb;    zz = dual_div(dual_const(1.0), bb);
//    term = c5*f1 + c6*f2 + gg + hh;
//      term = dual_add(dual_add(dual_add_ndm(dual_mul_num(f1, c5),
//                                            c6, f2), gg), hh);
//    term *= c7;          term = dual_mul(term, dual_const(c7));  (compound:
//                                                  full dual, not dual_mul_num)
//    fr = term * den;     fr = dual_mul_num(term, den);
//    rr = rev * term;     rr = dual_mul(rev, term);   (rev a Dual here)
//
//  A long product whose first factor is Real (left fold, Real first):
//    term = t0 * t1 * ... * t6;   (t0 Real, t1..t6 Dual)
//      s = dual_mul_num(t1, t0); s = dual_mul(s, t2); ... s = dual_mul(s, t6);
//    t6 = admath::pow(tf.t9, a6);  dual_pow_num(t9, a6)
//
//  Number sub-expressions, eval(), dual/number and nested exp(-pow(...)):
//    xx = rc * c1 * tf.t9i32 * admath::exp(-c2*tf.t9i);
//      xx = dual_mul(dual_mul_num(t9i32, rc * c1),
//                    dual_exp(dual_mul_num(t9i, -c2)));
//      (`rc * c1` is a number sub-expression: keep it as one double product)
//    r2 = aa + bb * cc;     r2 = dual_add(aa, dual_mul(bb, cc));
//    term = c3 * r2 * rb + xx;
//      term = dual_add(dual_mul(dual_mul_num(r2, c3), rb), xx);
//    uu = c4*admath::exp(-admath::pow(c5*tf.t9i, c6));
//      -> Mul(c4, Exp(Neg(Pow(NDM(c5, t9i), c6))))
//      uu = dual_mul(dual_const(c4),
//             dual_exp(dual_neg(dual_pow_num(dual_mul_num(t9i, c5), c6))));
//    vv = c7*admath::exp(-admath::pow(tf.t9/c8, c9));
//      tf.t9/c8 -> Mul(t9, Inv(c8)) -> dual_div_num(t9, c8)
//    zz = 1.0e0_rt + vv;       zz = dual_add(dual_const(1.0), vv);
//    f1 = c10 + yy * zzinv;    f1 = dual_add(dual_const(c10), dual_mul(yy, zzinv));
//    fr = autodiff::eval(term * den) * den;
//      fr = dual_mul_num(dual_mul_num(term, den), den);
//    if (tf.t9 > c11)  ->  if (t9.val > c11)
//
//  A number divided by a dual inside an exponential, and a rate with no
//  reverse term:
//    term = c1 * t9a56 * tf.t9i32 * admath::exp(-c2/t9a13 - c3*tf.t93);
//      arg = dual_add_ndm(dual_div(dual_const(-c2), t9a13), -c3, t93);
//      term = dual_mul(dual_mul(dual_mul_num(t9a56, c1), t9i32), dual_exp(arg));
//    t9a56 = admath::pow(t9a,5.0_rt/6.0_rt);  dual_pow_num(t9a, 5.0/6.0)
//    fr = den * term;   fr = dual_mul_num(term, den);
//    rr = 0.0e0_rt;     rr = dual_const(0.0);
//
//  A sum of n*d terms with a trailing dual*dual product:
//    admath::exp(-c1 * tf.t9i13 - c2*tf.t923 - c3*tf.t943 + c4*tf.t9*tf.t9)
//      s = dual_mul_num(t9i13, -c1); s = dual_add_ndm(s, -c2, t923);
//      s = dual_add_ndm(s, -c3, t943);
//      s = dual_add(s, dual_mul(dual_mul_num(t9, c4), t9));  -> dual_exp(s)
//
//  Real-times-sum, a number fraction times exp, a sum over a dual, and a
//  number*dual*dual product:
//    term3 = rc * (aa + bb);  term3 = dual_mul(dual_const(rc), dual_add(aa, bb));
//    ee = 1.0_rt/3.0_rt*admath::exp(-c1*tf.t9i);
//      ee = dual_mul(dual_const(1.0/3.0), dual_exp(dual_mul_num(t9i, -c1)));
//    term = (term1 + term2)/gg;  term = dual_div(dual_add(term1, term2), gg);
//    fr = den * rev * term;      fr = dual_mul(dual_mul_num(rev, den), term);
//                                (den Real, rev Dual: Mul(NDM(den, rev), term))
//
//  eval() then a trailing number factor, and a number first factor:
//    bb = autodiff::eval(c1 * tf.t9i32 * admath::exp(-c2*tf.t9i))*c3;
//      bb = dual_mul_num(dual_mul(dual_mul_num(t9i32, c1),
//                                 dual_exp(dual_mul_num(t9i, -c2))), c3);
//    gg = c4 * admath::pow(tf.t9,c5) * admath::exp(-c6*tf.t9i);
//      gg = dual_mul(dual_mul(dual_const(c4), dual_pow_num(t9, c5)),
//                    dual_exp(dual_mul_num(t9i, -c6)));
//      (pow(..) is not a dual, so c4*pow(..) is Mul, not NDM)
//
//  A clipped-temperature polynomial (z-polynomial rates):
//    z  = admath::min(c1, tf.t9);  z = dual_min(dual_const(c1), t9);
//    aa = 1.0e0_rt + c2*z + c3*z2 - c4*z3;
//      aa = dual_add_ndm(dual_add_ndm(dual_add_ndm(dual_const(1.0), c2, z),
//                                     c3, z2), -c4, z3);
//    admath::exp(-c5 * tf.t9i13 * aa)
//      dual_exp(dual_mul(dual_mul_num(t9i13, -c5), aa))
//
//  A network's postprocess_rate combining two sequential rates (r1, r2) into
//  one effective rate:
//    number_t rate_ir = 0.0e0_rt;          rate_ir = dual_const(0.0);
//    number_t denom = rates1.rr + rates2.fr;  denom = dual_add(r1.rr, r2.fr);
//    if (denom > 1.0e-30_rt)               if (denom.val > 1.0e-30)
//    number_t zz = 1.0_rt / denom;         zz = dual_div(dual_const(1.0), denom);
//    rate_ir = rates1.rr * zz;             rate_ir = dual_mul(r1.rr, zz);
//    rates.fr += rates1.fr * (1.0_rt - rate_ir);
//      fr = dual_add(fr, dual_mul(r1.fr, dual_sub(dual_const(1.0), rate_ir)));
//    rates.rr += rates2.rr * rate_ir;      rr = dual_add(rr, dual_mul(r2.rr, rate_ir));
//
//  networks/rhs.H
//    rates.fr *= density_term_forward;  (Real; also when it is 1.0!)
//      fr = dual_mul(fr, dual_const(density_term_forward));
//    rates.fr *= sc;  (sc dual)          fr = dual_mul(fr, sc);
//    const number_t sc = sc1 * sc2;      sc = dual_mul(sc1, sc2);
//    rates.fr *= data.forward_branching_ratio;  fr = dual_mul(fr, dual_const(ratio));
//    autodiff::derivative(rates.fr)      fr.grad
//    autodiff::dual dual_temp = burn_state.T; autodiff::seed(dual_temp);
//      dual_temp = dual_seed(T);
//
//  screening/screen.H fill_plasma_state (number_t = temp_t = Dual)
//    number_t m1 = 0.0_rt;  m1 += zion[n-1]*y(n);
//      m1 = dual_const(0.0); m1 = dual_add(m1, dual_const(zion[n-1]*y(n)));
//    number_t temp_ad = temp;              temp_ad = temp;
//    zbar = m1 * abar;                     zbar = dual_mul_num(m1, abar);
//    tempi = 1.0_rt / temp_ad;             tempi = dual_div(dual_const(1.0), temp_ad);
//    pp = admath::sqrt(rr*tempi*(z2bar + zbar));
//      pp = dual_sqrt(dual_mul(dual_mul_num(tempi, rr), dual_add(z2bar, zbar)));
//    state.qlam0z = 1.88e8_rt * tempi * pp;  dual_mul(dual_mul_num(tempi, 1.88e8), pp)
//    state.taufac = co2 * admath::cbrt(tempi); dual_mul(dual_const(co2), dual_cbrt(tempi))
//    xni = admath::cbrt(rr * zbar);        xni = dual_cbrt(dual_mul_num(zbar, rr));
//    state.n_e = zbar * rr * static_cast<amrex::Real>(C::n_A);
//      n_e = dual_mul_num(dual_mul_num(zbar, rr), n_A);
//    (Dual3: m1, m2 and temp_ad are seeded with dual3_seed(.., 1/2/3) right
//     after `number_t temp_ad = temp;`, as seed_array does.)
//
//  screening/screen.H actual_screen5
//    gamef = qq * gamp;           gamef = dual_mul_num(gamp, qq);
//    tau12 = state.taufac * scn_fac.aznut;  tau12 = dual_mul_num(taufac, aznut);
//    alph12 = gamef / tau12;      alph12 = dual_div(gamef, tau12);
//    if (alph12 > 1.6_rt) { alph12 = 1.6e0_rt; gamef = 1.6e0_rt * tau12; ...
//      if (alph12.val > 1.6) { alph12 = dual_const(1.6);
//                              gamef = dual_mul_num(tau12, 1.6); ...
//    gamp14 = admath::sqrt(admath::sqrt(gamp));  dual_sqrt(dual_sqrt(gamp))
//    cc = 0.896434e0_rt * gamp * scn_fac.zhat - 3.44740e0_rt * gamp14 * scn_fac.zhat2
//         - 0.5551e0_rt * (admath::log(gamp) + scn_fac.lzav) - 2.996e0_rt;
//      -> Add(Add(Add(Mul(NDM(.896434,gamp),zhat), Neg(Mul(NDM(3.4474,gamp14),zhat2))),
//                 Neg(Mul(0.5551, Add(Log(gamp), lzav)))), -2.996)
//      s  = dual_mul_num(dual_mul_num(gamp, 0.896434), zhat);
//      s  = dual_sub(s, dual_mul_num(dual_mul_num(gamp14, 3.44740), zhat2));
//      s  = dual_sub(s, dual_mul(dual_const(0.5551), dual_add_num(dual_log(gamp), lzav)));
//      cc = dual_add_num(s, -2.996);
//    a3 = alph12 * alph12 * alph12;   a3 = dual_mul(dual_mul(alph12, alph12), alph12);
//    rr = (5.0_rt/32.0_rt) - alph12*(0.014e0_rt + 0.0128e0_rt*alph12);
//      rr = dual_sub(dual_const(5.0/32.0),
//             dual_mul(alph12, dual_add_ndm(dual_const(0.014), 0.0128, alph12)));
//    tt = -0.0098e0_rt + 0.0048e0_rt*alph12;  dual_add_ndm(dual_const(-0.0098), 0.0048, alph12)
//    uu = 0.0055e0_rt + alph12*tt;    dual_add(dual_const(0.0055), dual_mul(alph12, tt))
//    vv = gamef * alph12 * uu;        dual_mul(dual_mul(gamef, alph12), uu)
//    h12 = cc - a3 * (ss + vv);       h12 = dual_sub(cc, dual_mul(a3, dual_add(ss, vv)));
//    xlgfac = admath::max(1.0_rt - 0.0562e0_rt*a3, 0.77_rt);
//      xlgfac = dual_max(dual_add_ndm(dual_const(1.0), -0.0562, a3), dual_const(0.77));
//    h12 += admath::log(xlgfac);      h12 = dual_add(h12, dual_log(xlgfac));
//    weight1 = inv_dg * (gamefs - gamef);
//      weight1 = dual_mul(dual_const(inv_dg), dual_sub(dual_const(gamefs), gamef));
//    weight2 = inv_dg * (gamef - gamefx);
//      weight2 = dual_mul(dual_const(inv_dg), dual_add_num(gamef, -gamefx));
//    h12 = h12w * weight1 + h12 * weight2;
//      h12 = dual_add(dual_mul(h12w, weight1), dual_mul(h12, weight2));
//    return admath::max(h12, 0.0_rt); return dual_max(h12, dual_const(0.0));
//
//  screening/screen.H debye_huckel / actual_screen
//    h_DH = z1z2 * admath::sqrt(3 * admath::powi<3>(Gamma_e) * z_ratio);
//      (3 is an int: NDM(3, powi-result))
//      h_DH = dual_mul(dual_const(z1z2),
//               dual_sqrt(dual_mul(dual_mul_num(dual_powi<3>(Gamma_e), 3.0), z_ratio)));
//    z_ratio = state.z2bar / state.zbar;  dual_div(z2bar, zbar)
//    if (log_scor <= 0.0) { return 1.0_rt; }
//      if (log_scor.val <= 0.0) -> result dual_const(1.0)
//    log_scor = admath::min(log_scor, log_scor_max);
//      log_scor = dual_min(log_scor, dual_const(300.0));
//    return admath::exp(log_scor);   dual_exp(log_scor)
//
// Warnings for translators:
//  * number*expr (expr not a plain dual) is a full dual multiply with a
//    zero-gradient constant: dual_mul(dual_const(n), EVAL(expr)), not
//    dual_mul_num. Only number*dual (an existing dual) is the cheap NDM.
//  * dual*number and number*dual are both dual_mul_num(d, n).
//  * a/b is never "multiply by inverse": use DIV (dual_div / dual_div_num).
//    dual_inv is only for a lone Inv node, which the burn path never makes.
//  * Compound assignment wraps the right side in a full dual (section C).
//  * dual_min/dual_max use autodiff's <= / >= rules, not std::min/max.
//  * Every multiply by a constant that the source executes must be executed,
//    even if the constant is 1.0 (density scaling with exponent 0): it can
//    turn a -0 gradient into +0 or val=inf into a NaN gradient.
// ===========================================================================

#ifndef MPL_DUAL_H
#define MPL_DUAL_H

#include "support.H"

namespace mpl {

struct Dual {
    Real val;
    Real grad;
};

// autodiff::dual_array<1,3>; grad[0] unused (always 0).
struct Dual3 {
    Real val;
    Real grad[4];
};

// ---------------------------------------------------------------------------
// autodiff::dual
// ---------------------------------------------------------------------------

// from util/autodiff/autodiff/forward/dual/dual.hpp:989 assign (number)
MPL_DEVICE Dual dual_const(Real n) {
    Dual r{n, 0.0};
    return r;
}

// from util/autodiff/autodiff/forward/utils/derivative.hpp:146 seed
MPL_DEVICE Dual dual_seed(Real x) {
    Dual r{x, 1.0};
    return r;
}

// from dual.hpp:1109 assignAdd (number)
MPL_DEVICE Dual dual_add_num(Dual s, Real n) {
    const Real v = s.val + n;
    Dual r{v, s.grad};
    return r;
}

// from dual.hpp:1113 assignAdd (dual)
MPL_DEVICE Dual dual_add(Dual s, Dual d) {
    const Real v = s.val + d.val;
    const Real g = s.grad + d.grad;
    Dual r{v, g};
    return r;
}

// from dual.hpp:1122 assignAdd (number * dual)
MPL_DEVICE Dual dual_add_ndm(Dual s, Real n, Dual d) {
    const Real v = s.val + n * d.val;
    const Real g = s.grad + n * d.grad;
    Dual r{v, g};
    return r;
}

// from dual.hpp:1179 assignSub (dual)
MPL_DEVICE Dual dual_sub(Dual s, Dual d) {
    const Real v = s.val - d.val;
    const Real g = s.grad - d.grad;
    Dual r{v, g};
    return r;
}

// from dual.hpp:1241 assignMul (number); also assign(number * dual) + scale
MPL_DEVICE Dual dual_mul_num(Dual s, Real n) {
    const Real v = s.val * n;
    const Real g = s.grad * n;
    Dual r{v, g};
    return r;
}

// from dual.hpp:1246 assignMul (dual)
MPL_DEVICE Dual dual_mul(Dual s, Dual d) {
    const Real g1 = s.grad * d.val;
    const Real g2 = g1 + s.val * d.grad;
    const Real v = s.val * d.val;
    Dual r{v, g2};
    return r;
}

// from dual.hpp:1316 assignDiv (number)
MPL_DEVICE Dual dual_div_num(Dual s, Real n) {
    const Real v = s.val / n;
    const Real g = s.grad / n;
    Dual r{v, g};
    return r;
}

// from dual.hpp:1321 assignDiv (dual)
MPL_DEVICE Dual dual_div(Dual s, Dual d) {
    const Real v = s.val / d.val;
    const Real g1 = s.grad - v * d.grad;
    const Real g2 = g1 / d.val;
    Dual r{v, g2};
    return r;
}

// from dual.hpp:1605 apply NegOp (and negate)
MPL_DEVICE Dual dual_neg(Dual s) {
    Dual r{-s.val, -s.grad};
    return r;
}

// from dual.hpp:1612 apply InvOp
MPL_DEVICE Dual dual_inv(Dual s) {
    const Real v = 1.0 / s.val;
    const Real g = s.grad * (-v * v);
    Dual r{v, g};
    return r;
}

// from dual.hpp:1687 apply ExpOp
MPL_DEVICE Dual dual_exp(Dual s) {
    const Real v = math::exp(s.val);
    const Real g = s.grad * v;
    Dual r{v, g};
    return r;
}

// from dual.hpp:1694 apply LogOp
MPL_DEVICE Dual dual_log(Dual s) {
    const Real aux = 1.0 / s.val;
    const Real v = math::log(s.val);
    const Real g = s.grad * aux;
    Dual r{v, g};
    return r;
}

// from dual.hpp:1702 apply Log10Op
MPL_DEVICE Dual dual_log10(Dual s) {
    constexpr Real ln10 = 2.3025850929940456840179914546843;
    const Real aux = 1.0 / (ln10 * s.val);
    const Real v = math::log10(s.val);
    const Real g = s.grad * aux;
    Dual r{v, g};
    return r;
}

// from dual.hpp:1711 apply SqrtOp
MPL_DEVICE Dual dual_sqrt(Dual s) {
    const Real v = math::sqrt(s.val);
    const Real g = s.grad * (0.5 / v);
    Dual r{v, g};
    return r;
}

// from util/microphysics_autodiff.H:345 apply CbrtOp
MPL_DEVICE Dual dual_cbrt(Dual s) {
    const Real v = math::cbrt(s.val);
    const Real g = s.grad * (1.0 / (3.0 * v * v));
    Dual r{v, g};
    return r;
}

// from dual.hpp:1718 apply AbsOp
MPL_DEVICE Dual dual_abs(Dual s) {
    const Real sgn = s.val < 0.0 ? -1.0 : (s.val > 0.0 ? 1.0 : 0.0);
    const Real g = s.grad * sgn;
    const Real v = math::abs(s.val);
    Dual r{v, g};
    return r;
}

// from dual.hpp:1393 assignPow (number)
MPL_DEVICE Dual dual_pow_num(Dual s, Real n) {
    const Real aux = math::pow(s.val, n);
    const Real g = s.grad * (n * aux / s.val);
    Dual r{aux, g};
    return r;
}

// from dual.hpp:1399 assignPow (dual)
MPL_DEVICE Dual dual_pow(Dual s, Dual d) {
    const Real aux1 = math::pow(s.val, d.val);
    const Real aux2 = math::log(s.val);
    const Real g1 = s.grad * (d.val / s.val);
    const Real g2 = g1 + aux2 * d.grad;
    const Real g3 = g2 * aux1;
    Dual r{aux1, g3};
    return r;
}

// from dual.hpp:921 min: (x <= y) ? x : y
MPL_DEVICE Dual dual_min(Dual x, Dual y) {
    Dual r = y;
    if (x.val <= y.val) {
        r = x;
    }
    return r;
}

// from dual.hpp:929 max: (x >= y) ? x : y
MPL_DEVICE Dual dual_max(Dual x, Dual y) {
    Dual r = y;
    if (x.val >= y.val) {
        r = x;
    }
    return r;
}

// admath::min/max with two `const` lvalue duals of the same type resolve to
// amrex::min/max (AMReX_Algorithm.H) = std::min/max, not to autodiff's:
// min = (y < x) ? y : x, max = (x < y) ? y : x. Not on the burn path.
MPL_DEVICE Dual dual_min_std(Dual x, Dual y) {
    Dual r = x;
    if (y.val < x.val) {
        r = y;
    }
    return r;
}

MPL_DEVICE Dual dual_max_std(Dual x, Dual y) {
    Dual r = x;
    if (x.val < y.val) {
        r = y;
    }
    return r;
}

// from util/microphysics_autodiff.H:379 powi -> amrex::Math::powi (AMReX_Math.H)
// on a dual: x*x is Mul(x, x), x*powi<N-1>(x) is Mul(x, dual), T(1)/p is
// Mul(dual_const(1), Inv(p)).
template <int N>
MPL_DEVICE Dual dual_powi(Dual x) {
    if constexpr (N < 0) {
        return dual_div(dual_const(1.0), dual_powi<-N>(x));
    } else if constexpr (N == 0) {
        return dual_const(1.0);
    } else if constexpr (N == 1) {
        return x;
    } else if constexpr (N == 2) {
        return dual_mul(x, x);
    } else if constexpr (N % 2 == 0) {
        return dual_powi<2>(dual_powi<N / 2>(x));
    } else {
        return dual_mul(x, dual_powi<N - 1>(x));
    }
}

// ---------------------------------------------------------------------------
// autodiff::dual_array<1,3> (GradArray<double,1,3> gradient). GradArray's
// operators: g *= s is g[i]*s; s * g is g[i]*s; -g is g[i]*-1.0.
// ---------------------------------------------------------------------------

MPL_DEVICE Dual3 dual3_const(Real n) {
    Dual3 r{};
    r.val = n;
    return r;
}

// from util/microphysics_autodiff.H:434 seed_array: grad(k) = 1
MPL_DEVICE Dual3 dual3_seed(Dual3 d, int k) {
    Dual3 r = d;
    r.grad[k] = 1.0;
    return r;
}

MPL_DEVICE Dual3 dual3_add_num(Dual3 s, Real n) {
    Dual3 r = s;
    r.val = s.val + n;
    return r;
}

MPL_DEVICE Dual3 dual3_add(Dual3 s, Dual3 d) {
    Dual3 r{};
    r.val = s.val + d.val;
    for (int i = 1; i <= 3; ++i) {
        r.grad[i] = s.grad[i] + d.grad[i];
    }
    return r;
}

MPL_DEVICE Dual3 dual3_add_ndm(Dual3 s, Real n, Dual3 d) {
    Dual3 r{};
    r.val = s.val + n * d.val;
    for (int i = 1; i <= 3; ++i) {
        r.grad[i] = s.grad[i] + d.grad[i] * n;
    }
    return r;
}

MPL_DEVICE Dual3 dual3_sub(Dual3 s, Dual3 d) {
    Dual3 r{};
    r.val = s.val - d.val;
    for (int i = 1; i <= 3; ++i) {
        r.grad[i] = s.grad[i] - d.grad[i];
    }
    return r;
}

MPL_DEVICE Dual3 dual3_mul_num(Dual3 s, Real n) {
    Dual3 r{};
    r.val = s.val * n;
    for (int i = 1; i <= 3; ++i) {
        r.grad[i] = s.grad[i] * n;
    }
    return r;
}

MPL_DEVICE Dual3 dual3_mul(Dual3 s, Dual3 d) {
    Dual3 r{};
    for (int i = 1; i <= 3; ++i) {
        const Real g1 = s.grad[i] * d.val;
        r.grad[i] = g1 + d.grad[i] * s.val;
    }
    r.val = s.val * d.val;
    return r;
}

MPL_DEVICE Dual3 dual3_div_num(Dual3 s, Real n) {
    Dual3 r{};
    r.val = s.val / n;
    for (int i = 1; i <= 3; ++i) {
        r.grad[i] = s.grad[i] / n;
    }
    return r;
}

MPL_DEVICE Dual3 dual3_div(Dual3 s, Dual3 d) {
    Dual3 r{};
    const Real v = s.val / d.val;
    r.val = v;
    for (int i = 1; i <= 3; ++i) {
        const Real g1 = s.grad[i] - d.grad[i] * v;
        r.grad[i] = g1 / d.val;
    }
    return r;
}

MPL_DEVICE Dual3 dual3_neg(Dual3 s) {
    Dual3 r{};
    r.val = -s.val;
    for (int i = 1; i <= 3; ++i) {
        r.grad[i] = s.grad[i] * -1.0;
    }
    return r;
}

// scale every gradient component by f (GradArray operator*=)
MPL_DEVICE Dual3 dual3_scale_grad(Dual3 s, Real v, Real f) {
    Dual3 r{};
    r.val = v;
    for (int i = 1; i <= 3; ++i) {
        r.grad[i] = s.grad[i] * f;
    }
    return r;
}

MPL_DEVICE Dual3 dual3_inv(Dual3 s) {
    const Real v = 1.0 / s.val;
    return dual3_scale_grad(s, v, -v * v);
}

MPL_DEVICE Dual3 dual3_exp(Dual3 s) {
    const Real v = math::exp(s.val);
    return dual3_scale_grad(s, v, v);
}

MPL_DEVICE Dual3 dual3_log(Dual3 s) {
    const Real aux = 1.0 / s.val;
    const Real v = math::log(s.val);
    return dual3_scale_grad(s, v, aux);
}

MPL_DEVICE Dual3 dual3_log10(Dual3 s) {
    constexpr Real ln10 = 2.3025850929940456840179914546843;
    const Real aux = 1.0 / (ln10 * s.val);
    const Real v = math::log10(s.val);
    return dual3_scale_grad(s, v, aux);
}

MPL_DEVICE Dual3 dual3_sqrt(Dual3 s) {
    const Real v = math::sqrt(s.val);
    return dual3_scale_grad(s, v, 0.5 / v);
}

MPL_DEVICE Dual3 dual3_cbrt(Dual3 s) {
    const Real v = math::cbrt(s.val);
    return dual3_scale_grad(s, v, 1.0 / (3.0 * v * v));
}

MPL_DEVICE Dual3 dual3_abs(Dual3 s) {
    const Real sgn = s.val < 0.0 ? -1.0 : (s.val > 0.0 ? 1.0 : 0.0);
    const Real v = math::abs(s.val);
    return dual3_scale_grad(s, v, sgn);
}

MPL_DEVICE Dual3 dual3_pow_num(Dual3 s, Real n) {
    const Real aux = math::pow(s.val, n);
    return dual3_scale_grad(s, aux, n * aux / s.val);
}

MPL_DEVICE Dual3 dual3_pow(Dual3 s, Dual3 d) {
    const Real aux1 = math::pow(s.val, d.val);
    const Real aux2 = math::log(s.val);
    const Real f = d.val / s.val;
    Dual3 r{};
    r.val = aux1;
    for (int i = 1; i <= 3; ++i) {
        const Real g1 = s.grad[i] * f;
        const Real g2 = g1 + d.grad[i] * aux2;
        r.grad[i] = g2 * aux1;
    }
    return r;
}

MPL_DEVICE Dual3 dual3_min(Dual3 x, Dual3 y) {
    Dual3 r = y;
    if (x.val <= y.val) {
        r = x;
    }
    return r;
}

MPL_DEVICE Dual3 dual3_max(Dual3 x, Dual3 y) {
    Dual3 r = y;
    if (x.val >= y.val) {
        r = x;
    }
    return r;
}

MPL_DEVICE Dual3 dual3_min_std(Dual3 x, Dual3 y) {
    Dual3 r = x;
    if (y.val < x.val) {
        r = y;
    }
    return r;
}

MPL_DEVICE Dual3 dual3_max_std(Dual3 x, Dual3 y) {
    Dual3 r = x;
    if (x.val < y.val) {
        r = y;
    }
    return r;
}

template <int N>
MPL_DEVICE Dual3 dual3_powi(Dual3 x) {
    if constexpr (N < 0) {
        return dual3_div(dual3_const(1.0), dual3_powi<-N>(x));
    } else if constexpr (N == 0) {
        return dual3_const(1.0);
    } else if constexpr (N == 1) {
        return x;
    } else if constexpr (N == 2) {
        return dual3_mul(x, x);
    } else if constexpr (N % 2 == 0) {
        return dual3_powi<2>(dual3_powi<N / 2>(x));
    } else {
        return dual3_mul(x, dual3_powi<N - 1>(x));
    }
}

} // namespace mpl

#endif
