// Helmholtz EOS for the Loom burn kernel (helmholtz; any network via the generated network_data.H;
// neutrinos off). See README.md for the translation rules.
//
// Translated from (Microphysics; first used in the first helmholtz port):
//   EOS/helmholtz/actual_eos.H, EOS/helmholtz/actual_eos_data.H,
//   EOS/eos_composition.H, interfaces/eos.H, interfaces/eos_override.H,
//   interfaces/eos_type.H, interfaces/eos_data.H.
//
// State types. The burn instantiates the EOS for exactly two state types,
// burn_t (eos(eos_input_rt) in integrator_setup, eos(eos_input_re) in
// update_thermodynamics) and eos_re_extra_t (eos(eos_input_re) plus
// composition_derivatives in the Strang jac when
// integrator.correct_jacobian_for_const_e is set). Both have the same EOS
// trait set (interfaces/eos_type.H concepts):
//   has_xn, has_energy, has_eta, has_dedA, has_dedZ        -> true
//   has_pressure, has_entropy, has_enthalpy, has_dpdA, has_dpdZ,
//   has_xne_xnp, has_pele_ppos                              -> false
// The functions below are templates over the state type T (used through a
// T*) with every `if constexpr (has_*<T>)` resolved for that trait set:
// dead branches are dropped, live ones are kept statement for statement.
// T must provide rho, T, e, xn[eos_nspec], dedr, dedT, dedA, dedZ, cv, mu,
// mu_e, y_e, eta, abar, zbar (mpl::eos_re_extra_t below, or the kernel's
// burn_t mirror). validate_eos.H static_asserts the trait set on the stock
// types. Do not instantiate these for a state type with a different set.
//
// Input mode. eos_input_t stays a runtime value (an int with the stock
// enumerator values); the stock runtime switch in reset_inputs is an if/else
// chain and every mode keeps its (possibly empty, for this trait set) branch.
// The burn uses eos_input_rt and eos_input_re.
//
// Runtime parameters come in through EosParams, the helmholtz tables
// through EosTables (device pointers with the stock C layout). GPU branches
// are followed (AMREX_USE_GPU defined): no input checks, no printing.
#ifndef MPL_EOS_H
#define MPL_EOS_H

#include "support.H"
#include "network_data.H"

namespace mpl {

// network_properties.H NumSpec (network_data.H, generated)
constexpr int eos_nspec = NumSpec;

// interfaces/eos_type.H: enum eos_input_t : std::uint8_t
constexpr int eos_input_rt = 0;
constexpr int eos_input_rh = 1;
constexpr int eos_input_tp = 2;
constexpr int eos_input_rp = 3;
constexpr int eos_input_re = 4;
constexpr int eos_input_ps = 5;
constexpr int eos_input_ph = 6;
constexpr int eos_input_th = 7;

// interfaces/eos_type.H: namespace EOS, enum eos_iter : std::uint8_t
namespace EOS {
constexpr int itemp = 1;
constexpr int idens = 2;
constexpr int iener = 3;
constexpr int ienth = 4;
constexpr int ientr = 5;
constexpr int ipres = 6;
} // namespace EOS

// EOS/helmholtz/actual_eos_data.H: table geometry and constants
namespace helmholtz {
constexpr int imax = 541;
constexpr int jmax = 201;

constexpr Real tlo   = 3.0e0;
constexpr Real thi   = 13.0e0;
constexpr Real tstp  = (thi - tlo) / (static_cast<Real>(jmax-1));
constexpr Real tstpi = 1.0e0 / tstp;
constexpr Real dlo   = -12.0e0;
constexpr Real dhi   = 15.0e0;
constexpr Real dstp  = (dhi - dlo) / (static_cast<Real>(imax-1));
constexpr Real dstpi = 1.0e0 / dstp;

// 2006 CODATA physical constants
constexpr Real h = 6.6260689633e-27;
constexpr Real avo_eos = 6.0221417930e23;
constexpr Real kerg = 1.380650424e-16;
constexpr Real amu = 1.66053878283e-24;
} // namespace helmholtz

// std::numbers::pi_v<double>
constexpr Real eos_pi = 3.141592653589793238462643383279502884;

// Flat indices into the helmholtz tables, matching the stock C arrays
//   f[jmax][imax][9], dpdf/ef/xf[jmax][imax][4]   (row-major).
#define MPL_HELM_F(j, i, m) ((j) * (::mpl::helmholtz::imax * 9) + (i) * 9 + (m))
#define MPL_HELM_F4(j, i, m) ((j) * (::mpl::helmholtz::imax * 4) + (i) * 4 + (m))

// Runtime parameters the EOS reads (filled by fill_eos_params in
// validate_eos.H from eos_rp::*, helmholtz::* and EOSData::*).
struct EosParams {
    // eos_rp (extern_parameters)
    Real eos_prad_limiter_rho_c;
    Real eos_prad_limiter_delta_rho;
    // helmholtz:: (actual_eos_data.H, set in actual_eos_init)
    int helmholtz_do_coulomb;
    int helmholtz_input_is_constant;
    Real helmholtz_ttol;
    Real helmholtz_dtol;
    // EOSData:: (eos_data.H, set in eos_init)
    Real EOSData_mintemp;
    Real EOSData_maxtemp;
    Real EOSData_mindens;
    Real EOSData_maxdens;
    Real EOSData_mine;
    Real EOSData_maxe;
};

// helmholtz tables (device pointers; layouts as the stock arrays).
// dpdf and xf are only read by state types with pressure / xne (not on the
// burn path); they are carried for completeness.
struct EosTables {
    const Real* f;        // [jmax][imax][9], MPL_HELM_F
    const Real* dpdf;     // [jmax][imax][4], MPL_HELM_F4 (unused here)
    const Real* ef;       // [jmax][imax][4], MPL_HELM_F4
    const Real* xf;       // [jmax][imax][4], MPL_HELM_F4 (unused here)
    const Real* d;        // [imax]
    const Real* t;        // [jmax]
    const Real* dt_sav;   // [jmax]
    const Real* dt2_sav;  // [jmax]
    const Real* dti_sav;  // [jmax]
    const Real* dt2i_sav; // [jmax]
    const Real* dd_sav;   // [imax]
    const Real* dd2_sav;  // [imax]
    const Real* ddi_sav;  // [imax]
    const Real* dd2i_sav; // [imax]
};

// from interfaces/eos_type.H:107 eos_re_extra_t (eos_re_t base flattened in
// declaration order; NAUX_NET = 0, no NSE_NET)
struct eos_re_extra_t {
    Real rho{};
    Real T{};
    Real e{};
    Real xn[eos_nspec]{};
    Real dedT{};
    Real dedr{};

    Real cv{};

    Real mu{};
    Real mu_e{};
    Real y_e{};
    Real eta{};

    Real abar{};
    Real zbar{};

    Real dedA{};
    Real dedZ{};
};

// from EOS/eos_composition.H:18 eos_xderivs_t
struct eos_xderivs_t {
    Real dedX[eos_nspec];
    Real dpdX[eos_nspec];
    Real dhdX[eos_nspec];
};

// ---------------------------------------------------------------------------
// EOS/eos_composition.H
// ---------------------------------------------------------------------------

// from EOS/eos_composition.H:50 composition (no AUX_THERMO)
template <typename T>
MPL_DEVICE void composition(T* state)
{
    // network_properties.H zion[], aion_inv[] (network_data.H lists)
    const Real zion[eos_nspec] = {MPL_NET_ZION_LIST};
    const Real aion_inv[eos_nspec] = {MPL_NET_AION_INV_LIST};

    state->y_e = 0.0;
    state->abar = 0.0;
    for (int n = 0; n < eos_nspec; n++) {
        state->y_e += state->xn[n] * zion[n] * aion_inv[n];
        state->abar += state->xn[n] * aion_inv[n];
    }
    state->abar = 1.0 / state->abar;

    state->mu_e = 1.0 / state->y_e;
    state->zbar = state->abar * state->y_e;
}

// from EOS/eos_composition.H:82 composition_derivatives
// (has_dpdA/has_dpdZ false -> dpdX = 0; has_dedA && has_dedZ true;
//  has_pressure false -> dhdX = 0). Stock returns eos_xderivs_t by value.
template <typename T>
MPL_DEVICE void composition_derivatives(const T* state, eos_xderivs_t* state_xderivs)
{
    // network_properties.H aion[], zion[], aion_inv[] (network_data.H lists)
    const Real aion[eos_nspec] = {MPL_NET_AION_LIST};
    const Real zion[eos_nspec] = {MPL_NET_ZION_LIST};
    const Real aion_inv[eos_nspec] = {MPL_NET_AION_INV_LIST};

    for (int n = 0; n < eos_nspec; n++) {
        state_xderivs->dpdX[n] = 0.0;

        state_xderivs->dedX[n] =
            state->dedA * (state->abar * aion_inv[n]) * (aion[n] - state->abar) +
            state->dedZ * (state->abar * aion_inv[n]) * (zion[n] - state->zbar);

        state_xderivs->dhdX[n] = 0.0;
    }
}

// ---------------------------------------------------------------------------
// EOS/helmholtz/actual_eos.H: interpolation polynomials
// ---------------------------------------------------------------------------

// from EOS/helmholtz/actual_eos.H:51 psi0
MPL_DEVICE Real psi0(Real z)
{
    return z * z * z * (z * (-6.0e0 * z + 15.0e0) -10.0e0) + 1.0e0;
}

// from EOS/helmholtz/actual_eos.H:57 dpsi0
MPL_DEVICE Real dpsi0(Real z)
{
    return z * z * (z * (-30.0e0*z + 60.0e0) - 30.0e0);
}

// from EOS/helmholtz/actual_eos.H:63 ddpsi0
MPL_DEVICE Real ddpsi0(Real z)
{
    return z * (z * (-120.0e0 * z + 180.0e0) - 60.0e0);
}

// from EOS/helmholtz/actual_eos.H:70 psi1
MPL_DEVICE Real psi1(Real z)
{
    return z * ( z * z * (z * (-3.0e0 * z + 8.0e0) - 6.0e0) + 1.0e0);
}

// from EOS/helmholtz/actual_eos.H:76 dpsi1
MPL_DEVICE Real dpsi1(Real z)
{
    return z * z * ( z * (-15.0e0 * z + 32.0e0) - 18.0e0) + 1.0e0;
}

// from EOS/helmholtz/actual_eos.H:82 ddpsi1
MPL_DEVICE Real ddpsi1(Real z)
{
    return z * (z * (-60.0e0 * z + 96.0e0) - 36.0e0);
}

// from EOS/helmholtz/actual_eos.H:89 psi2
MPL_DEVICE Real psi2(Real z)
{
    return 0.5e0 * z * z * ( z * ( z * (-z + 3.0e0) - 3.0e0) + 1.0e0);
}

// from EOS/helmholtz/actual_eos.H:95 dpsi2
MPL_DEVICE Real dpsi2(Real z)
{
    return 0.5e0 * z * (z * (z * (-5.0e0 * z + 12.0e0) - 9.0e0) + 2.0e0);
}

// from EOS/helmholtz/actual_eos.H:101 ddpsi2
MPL_DEVICE Real ddpsi2(Real z)
{
    return 0.5e0 * (z * ( z * (-20.0e0 * z + 36.0e0) - 18.0e0) + 2.0e0);
}

// from EOS/helmholtz/actual_eos.H:107 fwt
MPL_DEVICE void fwt(const Real* fi, const Real* wt, Real* fwtr)
{
    fwtr[0] = fi[ 0]*wt[0] + fi[ 1]*wt[1] + fi[ 2]*wt[2] + fi[18]*wt[3] + fi[19]*wt[4] + fi[20]*wt[5];
    fwtr[1] = fi[ 3]*wt[0] + fi[ 5]*wt[1] + fi[ 7]*wt[2] + fi[21]*wt[3] + fi[23]*wt[4] + fi[25]*wt[5];
    fwtr[2] = fi[ 4]*wt[0] + fi[ 6]*wt[1] + fi[ 8]*wt[2] + fi[22]*wt[3] + fi[24]*wt[4] + fi[26]*wt[5];
    fwtr[3] = fi[ 9]*wt[0] + fi[10]*wt[1] + fi[11]*wt[2] + fi[27]*wt[3] + fi[28]*wt[4] + fi[29]*wt[5];
    fwtr[4] = fi[12]*wt[0] + fi[14]*wt[1] + fi[16]*wt[2] + fi[30]*wt[3] + fi[32]*wt[4] + fi[34]*wt[5];
    fwtr[5] = fi[13]*wt[0] + fi[15]*wt[1] + fi[17]*wt[2] + fi[31]*wt[3] + fi[33]*wt[4] + fi[35]*wt[5];
}

// from EOS/helmholtz/actual_eos.H:123 xpsi0
MPL_DEVICE Real xpsi0(Real z)
{
    return z * z * (2.0e0 * z - 3.0e0) + 1.0;
}

// from EOS/helmholtz/actual_eos.H:130 xpsi1
MPL_DEVICE Real xpsi1(Real z)
{
    return z * (z * (z - 2.0e0) + 1.0e0);
}

// ---------------------------------------------------------------------------
// EOS/helmholtz/actual_eos.H: physics contributions
// ---------------------------------------------------------------------------

// One term of the constexpr_for<0, 16> sums in apply_electrons:
// electron_table_indexing::map[k] = {dj, di, n} (actual_eos.H:144).
#define MPL_HELM_EF_TERM(k, dj, di, n) \
    etaele = etaele + (tab.ef[MPL_HELM_F4(jat + (dj), iat + (di), (n))] * wdt[(k)])

// from EOS/helmholtz/actual_eos.H:155 apply_electrons
// (has_pressure, has_dpdA, has_dpdZ, has_xne_xnp, has_pele_ppos, has_entropy
//  false: the dpdf / xf sums and the pressure/entropy updates are dropped)
template <typename T>
MPL_DEVICE void apply_electrons(T* state, EosTables tab)
{
    // assume complete ionization
    Real ytot1 = 1.0e0 / state->abar;

    // enter the table with ye*den
    Real din = state->y_e * state->rho;

    // hash locate this temperature and density
    int jat = int((math::log10(state->T) - helmholtz::tlo) * helmholtz::tstpi) + 1;
    jat = math::clamp(jat, 1, helmholtz::jmax-1) - 1;
    int iat = int((math::log10(din) - helmholtz::dlo) * helmholtz::dstpi) + 1;
    iat = math::clamp(iat, 1, helmholtz::imax-1) - 1;

    Real fi[36];

    // access the table locations only once
    for (int i = 0; i < 9; ++i) {
        fi[i     ] = tab.f[MPL_HELM_F(jat  , iat  , i)]; // f, ft, ftt, fd, fdd, fdt, fddt, fdtt, fddtt
        fi[i +  9] = tab.f[MPL_HELM_F(jat  , iat+1, i)];
        fi[i + 18] = tab.f[MPL_HELM_F(jat+1, iat  , i)];
        fi[i + 27] = tab.f[MPL_HELM_F(jat+1, iat+1, i)];
    }

    // various differences
    Real xt  = math::max((state->T - tab.t[jat]) * tab.dti_sav[jat], 0.0e0);
    Real xd  = math::max((din - tab.d[iat]) * tab.ddi_sav[iat], 0.0e0);
    Real mxt = 1.0e0 - xt;
    Real mxd = 1.0e0 - xd;

    // the six density basis functions
    Real sid[6];

    sid[0] =  psi0(xd);
    sid[1] =  psi1(xd) * tab.dd_sav[iat];
    sid[2] =  psi2(xd) * tab.dd2_sav[iat];

    sid[3] =  psi0(mxd);
    sid[4] = -psi1(mxd) * tab.dd_sav[iat];
    sid[5] =  psi2(mxd) * tab.dd2_sav[iat];

    // derivatives of the weight functions
    Real dsid[6];

    dsid[0] =  dpsi0(xd) * tab.ddi_sav[iat];
    dsid[1] =  dpsi1(xd);
    dsid[2] =  dpsi2(xd) * tab.dd_sav[iat];

    dsid[3] = -dpsi0(mxd) * tab.ddi_sav[iat];
    dsid[4] =  dpsi1(mxd);
    dsid[5] = -dpsi2(mxd) * tab.dd_sav[iat];

    Real fwtr[6];
    Real sit[6];

    {
        // temperature basis functions

        sit[0] = psi0(xt);
        sit[1] = psi1(xt) * tab.dt_sav[jat];
        sit[2] = psi2(xt) * tab.dt2_sav[jat];

        sit[3] =  psi0(mxt);
        sit[4] = -psi1(mxt) * tab.dt_sav[jat];
        sit[5] =  psi2(mxt) * tab.dt2_sav[jat];

        fwt(fi, sit, fwtr);
    }

    Real free = 0.e0;
    Real df_d = 0.e0;

    for (int i = 0; i <= 5; ++i) {
        // the free energy
        free = free + fwtr[i] * sid[i];

        // derivative with respect to density
        df_d = df_d + fwtr[i] * dsid[i];
    }

    {
        // temperature derivative of weight functions
        // this was originally called dsit
        sit[0] =  dpsi0(xt) * tab.dti_sav[jat];
        sit[1] =  dpsi1(xt);
        sit[2] =  dpsi2(xt) * tab.dt_sav[jat];

        sit[3] = -dpsi0(mxt) * tab.dti_sav[jat];
        sit[4] =  dpsi1(mxt);
        sit[5] = -dpsi2(mxt) * tab.dt_sav[jat];

        fwt(fi, sit, fwtr);
    }

    Real df_t = 0.e0;
    Real df_dt = 0.e0;

    for (int i = 0; i <= 5; ++i) {
        // derivative with respect to temperature
        df_t += fwtr[i] * sid[i];

        // derivative with respect to temperature and density
        df_dt += fwtr[i] * dsid[i];
    }

    {
        // second derivatives of the weight functions
        // this was originally called ddsit

        sit[0] =  ddpsi0(xt) * tab.dt2i_sav[jat];
        sit[1] =  ddpsi1(xt) * tab.dti_sav[jat];
        sit[2] =  ddpsi2(xt);

        sit[3] =  ddpsi0(mxt) * tab.dt2i_sav[jat];
        sit[4] = -ddpsi1(mxt) * tab.dti_sav[jat];
        sit[5] =  ddpsi2(mxt);

        fwt(fi, sit, fwtr);
    }

    Real df_tt = 0.e0;
    for (int i = 0; i <= 5; ++i) {
        // derivative with respect to temperature**2
        df_tt = df_tt + fwtr[i] * sid[i];
    }

    // now get the pressure derivative with density, chemical potential, and
    // electron positron number densities
    // get the interpolation weight functions
    sit[0] = xpsi0(xt);
    sit[1] = xpsi1(xt) * tab.dt_sav[jat];

    sit[2] = xpsi0(mxt);
    sit[3] = -xpsi1(mxt) * tab.dt_sav[jat];

    sid[0] = xpsi0(xd);
    sid[1] = xpsi1(xd) * tab.dd_sav[iat];

    sid[2] = xpsi0(mxd);
    sid[3] = -xpsi1(mxd) * tab.dd_sav[iat];

    // Reuse subexpressions that would go into computing the
    // cubic interpolation.
    Real wdt[16];

    for (int i = 0; i <= 3; ++i) {
        wdt[i     ] = sid[0] * sit[i];
        wdt[i +  4] = sid[1] * sit[i];
        wdt[i +  8] = sid[2] * sit[i];
        wdt[i + 12] = sid[3] * sit[i];
    }

    // pressure derivative with density: has_pressure false, sum dropped
    // (dpepdd is only read under has_pressure)

    // electron chemical potential etaele (has_eta true): the
    // constexpr_for<0, 16> sum, unrolled in k order
    Real etaele{};
    MPL_HELM_EF_TERM( 0, 0, 0, 0);
    MPL_HELM_EF_TERM( 1, 0, 0, 1);
    MPL_HELM_EF_TERM( 2, 1, 0, 0);
    MPL_HELM_EF_TERM( 3, 1, 0, 1);
    MPL_HELM_EF_TERM( 4, 0, 0, 2);
    MPL_HELM_EF_TERM( 5, 0, 0, 3);
    MPL_HELM_EF_TERM( 6, 1, 0, 2);
    MPL_HELM_EF_TERM( 7, 1, 0, 3);
    MPL_HELM_EF_TERM( 8, 0, 1, 0);
    MPL_HELM_EF_TERM( 9, 0, 1, 1);
    MPL_HELM_EF_TERM(10, 1, 1, 0);
    MPL_HELM_EF_TERM(11, 1, 1, 1);
    MPL_HELM_EF_TERM(12, 0, 1, 2);
    MPL_HELM_EF_TERM(13, 0, 1, 3);
    MPL_HELM_EF_TERM(14, 1, 1, 2);
    MPL_HELM_EF_TERM(15, 1, 1, 3);

    // electron + positron number densities: has_xne_xnp false, dropped

    // the desired electron-positron thermodynamic quantities
    Real x       = din * din;
    Real pele    = x * df_d;
    Real dpepdt  = x * df_dt;
    (void) pele;
    (void) dpepdt;
    // dpepda / dpepdz: has_pressure false, dropped

    x            = state->y_e * state->y_e;
    Real sele    = -df_t * state->y_e;
    Real dsepdt  = -df_tt * state->y_e;
    Real dsepdd  = -df_dt * x;

    Real eele    = state->y_e*free + state->T * sele;
    Real deepdt  = state->T * dsepdt;
    Real deepdd  = x * df_d + state->T * dsepdd;

    // has_dedA
    Real deepda{};
    {
        Real dsepda = ytot1 * (state->y_e * df_dt * din - sele);
        deepda = -state->y_e * ytot1 * (free +  df_d * din) + state->T * dsepda;
    }

    // has_dedZ
    Real deepdz{};
    {
        Real dsepdz = -ytot1 * (state->y_e * df_dt * state->rho  + df_t);
        deepdz = ytot1* (free + state->y_e * df_d * state->rho) + state->T * dsepdz;
    }

    // has_energy (pressure / entropy updates dropped)
    state->e    = state->e + eele;
    state->dedT = state->dedT + deepdt;
    state->dedr = state->dedr + deepdd;
    state->dedA = state->dedA + deepda;
    state->dedZ = state->dedZ + deepdz;

    // has_eta
    state->eta = etaele;
}

#undef MPL_HELM_EF_TERM

// from EOS/helmholtz/actual_eos.H:452 apply_ions
// (has_entropy, has_pressure false: sion and the pressure update dropped)
template <typename T>
MPL_DEVICE void apply_ions(T* state)
{
    Real deni = 1.0e0 / state->rho;
    Real tempi = 1.0e0 / state->T;
    (void) tempi;

    Real ytot1   = 1.0e0 / state->abar;
    Real xni     = helmholtz::avo_eos * ytot1 * state->rho;
    Real dxnidd  = helmholtz::avo_eos * ytot1;
    Real dxnida  = -xni * ytot1;

    Real kt = helmholtz::kerg * state->T;

    Real pion    = xni * kt;
    Real dpiondd = dxnidd * kt;
    Real dpiondt = xni * helmholtz::kerg;
    Real dpionda = dxnida * kt;
    Real dpiondz = 0.0e0;
    (void) dpiondz;

    Real eion    = 1.5e0 * pion * deni;
    Real deiondd = (1.5e0 * dpiondd - eion) * deni;
    Real deiondt = 1.5e0 * dpiondt * deni;
    Real deionda = 1.5e0 * dpionda * deni;
    Real deiondz = 0.0e0;

    // has_energy
    state->e    = state->e + eion;
    state->dedT = state->dedT + deiondt;
    state->dedr = state->dedr + deiondd;
    state->dedA = state->dedA + deionda;
    state->dedZ = state->dedZ + deiondz;
}

// from EOS/helmholtz/actual_eos.H:533 apply_radiation (no RADIATION)
// (has_entropy, has_pressure false: srad and the pressure set dropped)
template <typename T>
MPL_DEVICE void apply_radiation(T* state, const EosParams* params)
{
    constexpr Real clight  = 2.99792458e10;
    constexpr Real ssol    = 5.67051e-5;
    constexpr Real asol    = 4.0e0 * ssol / clight;
    constexpr Real asoli3  = asol/3.0e0;

    Real deni = 1.0e0 / state->rho;
    Real tempi = 1.0e0 / state->T;

    Real prad = asoli3 * state->T * state->T * state->T * state->T;

    // optional smoothing of the radiation pressure at low density
    // (eos.prad_limiter_rho_c / eos.prad_limiter_delta_rho, default off)
    if (params->eos_prad_limiter_rho_c > 0.0e0 && params->eos_prad_limiter_delta_rho > 0.0e0) {
        prad = prad * 0.5e0 * (1.0e0 + math::tanh((state->rho - params->eos_prad_limiter_rho_c) / params->eos_prad_limiter_delta_rho));
    }

    Real dpraddd = 0.0e0;
    (void) dpraddd;
    Real dpraddt = 4.0e0 * prad * tempi;

    Real erad    = 3.0e0 * prad*deni;
    Real deraddd = -erad * deni;
    Real deraddt = 3.0e0 * dpraddt * deni;

    // Note that unlike the other terms, radiation
    // sets these terms instead of adding to them,
    // since it comes first.

    // has_energy
    state->e    = erad;
    state->dedr = deraddd;
    state->dedT = deraddt;
    state->dedA = 0.0;
    state->dedZ = 0.0;
}

// from EOS/helmholtz/actual_eos.H:618 apply_coulomb_corrections
// (has_entropy, has_pressure, has_dpdA, has_dpdZ false; has_energy,
//  has_dedA, has_dedZ true)
template <typename T>
MPL_DEVICE void apply_coulomb_corrections(T* state)
{
    // Constants used for the Coulomb corrections
    constexpr Real a1 = -0.898004e0;
    constexpr Real b1 =  0.96786e0;
    constexpr Real c1 =  0.220703e0;
    constexpr Real d1 = -0.86097e0;
    constexpr Real a2 =  0.29561e0;
    constexpr Real b2 =  1.9885e0;
    constexpr Real c2 =  0.288675e0;
    constexpr Real qe   = 4.8032042712e-10;
    constexpr Real esqu = qe * qe;
    constexpr Real onethird = 1.0e0/3.0e0;
    constexpr Real forth = 4.0e0/3.0e0;
    constexpr Real pi = eos_pi;

    Real pcoul    = 0.e0;
    Real dpcouldd = 0.e0;
    Real dpcouldt = 0.e0;
    Real ecoul    = 0.e0;
    Real decouldd = 0.e0;
    Real decouldt = 0.e0;

    Real dpcoulda = 0.e0;
    Real dpcouldz = 0.e0;
    Real decoulda = 0.e0;
    Real decouldz = 0.e0;

    Real s = 0.0; // temporary variables
    Real z = 0.0;

    // uniform background corrections only
    // from yakovlev & shalybkov 1989
    // lami is the average ion separation
    // plasg is the plasma coupling parameter

    Real ytot1 = 1.0e0 / state->abar;
    Real xni     = helmholtz::avo_eos * ytot1 * state->rho;
    Real dxnidd  = helmholtz::avo_eos * ytot1;
    Real dxnida  = -xni * ytot1;

    Real kt      = helmholtz::kerg * state->T;
    Real ktinv   = 1.0e0 / kt;

    z             = forth * pi;
    s             = z * xni;
    Real dsdd     = z * dxnidd;
    Real dsda     = z * dxnida;

    Real inv_lami = math::cbrt(s);
    Real lami     = 1.0e0 / inv_lami;
    z             = -onethird * lami;
    Real lamidd   = z * dsdd / s;
    Real lamida   = z * dsda / s;

    Real plasg    = state->zbar * state->zbar * esqu * ktinv * inv_lami;
    z             = -plasg * inv_lami;
    Real plasgdd  = z * lamidd;
    Real plasgdt  = -plasg*ktinv * helmholtz::kerg;
    Real plasgda  = z * lamida;
    Real plasgdz  = 2.0e0 * plasg/state->zbar;

    // .yakovlev & shalybkov 1989 equations 82, 85, 86, 87
    if (plasg >= 1.0e0)
    {
        Real x = math::pow(plasg, 0.25e0);
        Real y = helmholtz::avo_eos * ytot1 * helmholtz::kerg;
        ecoul    = y * state->T * (a1 * plasg + b1 * x + c1 / x + d1);
        pcoul    = onethird * state->rho * ecoul;

        y        = helmholtz::avo_eos*ytot1*kt*(a1 + 0.25e0/plasg*(b1*x - c1/x));
        decouldd = y * plasgdd;
        decouldt = y * plasgdt + ecoul/state->T;

        decoulda = y * plasgda - ecoul * ytot1;
        decouldz = y * plasgdz;

        y        = onethird * state->rho;
        dpcouldd = onethird * ecoul + y * decouldd;
        dpcouldt = y * decouldt;

        dpcoulda = y * decoulda;
        dpcouldz = y * decouldz;

        // yakovlev & shalybkov 1989 equations 102, 103, 104
    }
    else if (plasg < 1.0e0)
    {

        Real pion    = xni * kt;
        Real dpiondd = dxnidd * kt;
        Real dpiondt = xni * helmholtz::kerg;
        Real dpionda = dxnida * kt;
        Real dpiondz = 0.0e0;

        Real x = plasg * math::sqrt(plasg);
        Real y = math::pow(plasg, b2);
        z        = c2 * x - onethird * a2 * y;
        pcoul    = -pion * z;
        ecoul    = 3.0e0 * pcoul / state->rho;

        s        = 1.5e0 * c2 * x / plasg - onethird * a2 * b2 * y / plasg;
        dpcouldd = -dpiondd * z - pion * s * plasgdd;
        dpcouldt = -dpiondt * z - pion * s * plasgdt;
        dpcoulda = -dpionda * z - pion * s * plasgda;
        dpcouldz = -dpiondz * z - pion * s * plasgdz;

        s        = 3.0e0 / state->rho;
        decouldd = s * dpcouldd - ecoul / state->rho;
        decouldt = s * dpcouldt;
        decoulda = s * dpcoulda;
        decouldz = s * dpcouldz;
    }

    // Disable Coulomb corrections if they cause
    // the energy or pressure to go negative.

    // std::numeric_limits<amrex::Real>::max(); has_pressure false
    Real p_temp = __DBL_MAX__;
    Real e_temp = __DBL_MAX__;

    // has_energy
    e_temp = state->e + ecoul;

    if (p_temp <= 0.0e0 || e_temp <= 0.0e0)
    {
        pcoul    = 0.0e0;
        dpcouldd = 0.0e0;
        dpcouldt = 0.0e0;
        ecoul    = 0.0e0;
        decouldd = 0.0e0;
        decouldt = 0.0e0;

        dpcoulda = 0.0e0;
        dpcouldz = 0.0e0;
        decoulda = 0.0e0;
        decouldz = 0.0e0;
    }
    (void) pcoul;
    (void) dpcouldd;
    (void) dpcouldt;
    (void) dpcoulda;
    (void) dpcouldz;

    // has_energy
    state->e    = state->e + ecoul;
    state->dedr = state->dedr + decouldd;
    state->dedT = state->dedT + decouldt;
    state->dedA = state->dedA + decoulda;
    state->dedZ = state->dedZ + decouldz;
}

// ---------------------------------------------------------------------------
// EOS/helmholtz/actual_eos.H: Newton iteration driver
// ---------------------------------------------------------------------------

// from EOS/helmholtz/actual_eos.H:827 prepare_for_iterations
// (only the eos_input_re branch has a live body for this trait set; the
//  #ifndef AMREX_USE_GPU else-abort is dropped)
template <typename T>
MPL_DEVICE void prepare_for_iterations(int input, const T* state,
                                       bool& single_iter, Real& v_want, Real& v1_want, Real& v2_want,
                                       int& var, int& dvar, int& var1, int& var2)
{
    (void) v1_want;
    (void) v2_want;
    (void) var1;
    (void) var2;

    single_iter = true;

    if (input == eos_input_rt) {
        // Nothing to do here.
    }
    else if (input == eos_input_rh) {
        // has_enthalpy false
    }
    else if (input == eos_input_tp) {
        // has_pressure false
    }
    else if (input == eos_input_rp) {
        // has_pressure false
    }
    else if (input == eos_input_re) {
        // has_energy
        v_want = state->e;
        var  = EOS::iener;
        dvar = EOS::itemp;
    }
    else if (input == eos_input_ps) {
        // has_pressure && has_entropy false
    }
    else if (input == eos_input_ph) {
        // has_pressure && has_enthalpy false
    }
    else if (input == eos_input_th) {
        // has_enthalpy false
    }
}

// from EOS/helmholtz/actual_eos.H:921 single_iter_update
// (only the iener branches have live bodies for this trait set)
template <typename T>
MPL_DEVICE void single_iter_update(T* state, int var, int dvar,
                                   Real v_want, bool& converged, const EosParams* params)
{
    Real x = 0.0;
    Real v = 0.0;
    Real dvdx = 0.0;
    Real xtol = 0.0;
    Real smallx = 0.0;

    if (dvar == EOS::itemp) {

        x = state->T;
        smallx = params->EOSData_mintemp;
        xtol = params->helmholtz_ttol;

        if (var == EOS::ipres) {
            // has_pressure false
        }
        else if (var == EOS::iener) {
            v    = state->e;
            dvdx = state->dedT;
        }
        else if (var == EOS::ientr) {
            // has_entropy false
        }
        else if (var == EOS::ienth) {
            // has_enthalpy false
        }

    }
    else // dvar == density
    {

        x = state->rho;
        smallx = params->EOSData_mindens;
        xtol = params->helmholtz_dtol;

        if (var == EOS::ipres) {
            // has_pressure false
        }
        else if (var == EOS::iener) {
            v    = state->e;
            dvdx = state->dedr;
        }
        else if (var == EOS::ientr) {
            // has_entropy false
        }
        else if (var == EOS::ienth) {
            // has_enthalpy false
        }

    }

    // Now do the calculation for the next guess for T/rho
    Real xnew = x - (v - v_want) / dvdx;

    // Don't let the temperature/density change by more than a factor of two
    xnew = math::clamp(xnew, 0.5 * x, 2.0 * x);

    // Don't let us freeze/evacuate
    xnew = math::max(smallx, xnew);

    // Store the new temperature/density

    if (dvar == EOS::itemp) {
        state->T = xnew;
    }
    else
    {
        state->rho = xnew;
    }

    // Compute the error from the last iteration

    Real error = math::abs(xnew - x);

    if (error < xtol * x) converged = true;
}

// from EOS/helmholtz/actual_eos.H:1026 double_iter_update
// (only the iener branches have live bodies for this trait set; with this
//  trait set prepare_for_iterations never clears single_iter, so this is
//  statically reachable but not called)
template <typename T>
MPL_DEVICE void double_iter_update(T* state, int var1, int var2,
                                   Real v1_want, Real v2_want, bool& converged,
                                   const EosParams* params)
{
    Real v1 = 0.0, dv1dt = 0.0, dv1dr = 0.0, v2 = 0.0, dv2dt = 0.0, dv2dr = 0.0;

    // Figure out which variables we're using

    Real told = state->T;
    Real rold = state->rho;

    if (var1 == EOS::ipres) {
        // has_pressure false
    }
    else if (var1 == EOS::iener) {
        v1    = state->e;
        dv1dt = state->dedT;
        dv1dr = state->dedr;
    }
    else if (var1 == EOS::ientr) {
        // has_entropy false
    }
    else if (var1 == EOS::ienth) {
        // has_enthalpy false
    }

    if (var2 == EOS::ipres) {
        // has_pressure false
    }
    else if (var2 == EOS::iener) {
        v2    = state->e;
        dv2dt = state->dedT;
        dv2dr = state->dedr;
    }
    else if (var2 == EOS::ientr) {
        // has_entropy false
    }
    else if (var2 == EOS::ienth) {
        // has_enthalpy false
    }

    // Two functions, f and g, to iterate over
    Real v1i = v1_want - v1;
    Real v2i = v2_want - v2;

    // note that dfi/dT = - df/dT
    Real delr = (-v1i * dv2dt + v2i * dv1dt) / (dv2dr * dv1dt - dv2dt * dv1dr);

    Real rnew = rold + delr;

    Real tnew = told + (v1i - dv1dr * delr) / dv1dt;

    // Don't let the temperature or density change by more
    // than a factor of two
    tnew = math::clamp(tnew, 0.5e0 * told, 2.0e0 * told);
    rnew = math::clamp(rnew, 0.5e0 * rold, 2.0e0 * rold);

    // Don't let us freeze or evacuate
    tnew = math::max(params->EOSData_mintemp, tnew);
    rnew = math::max(params->EOSData_mindens, rnew);

    // Store the new temperature and density
    state->rho = rnew;
    state->T = tnew;

    // Compute the errors
    Real error1 = math::abs(rnew - rold);
    Real error2 = math::abs(tnew - told);

    if (error1 < params->helmholtz_dtol * rold && error2 < params->helmholtz_ttol * told) converged = true;
}

// from EOS/helmholtz/actual_eos.H:1137 finalize_state
// (has_pressure, has_enthalpy, has_entropy false)
template <typename T>
MPL_DEVICE void finalize_state(int input, T* state,
                               Real v_want, Real v1_want, Real v2_want,
                               const EosParams* params)
{
    (void) v1_want;
    (void) v2_want;

    // define mu -- the total mean molecular weight including both electrons and ions
    state->mu = 1.0e0 / (1.0e0 / state->abar + 1.0e0 / state->mu_e);

    // Specific heats (has_energy; Gamma_1 needs has_pressure)
    state->cv = state->dedT;

    if (params->helmholtz_input_is_constant) {

        if (input == eos_input_rh) {
            // has_enthalpy false
        }
        else if (input == eos_input_tp) {
            // has_pressure false
        }
        else if (input == eos_input_rp) {
            // has_pressure false
        }
        else if (input == eos_input_re) {
            // has_energy
            state->e = v_want;
        }
        else if (input == eos_input_ps) {
            // has_pressure, has_entropy false
        }
        else if (input == eos_input_ph) {
            // has_pressure, has_enthalpy false
        }
        else if (input == eos_input_th) {
            // has_enthalpy false
        }

    }
}

// from EOS/helmholtz/actual_eos.H:1235 actual_eos
template <typename T>
MPL_DEVICE void actual_eos(int input, T* state, const EosParams* params, EosTables tab)
{
    constexpr int max_newton = 100;

    // stock: `bool single_iter, converged; int var1, var2;` uninitialized;
    // single_iter and converged are assigned before use, var1/var2 are only
    // read by double_iter_update (not reached with this trait set).
    bool single_iter = false;
    bool converged = false;
    int var{}, dvar{}, var1{}, var2{};
    Real v_want{}, v1_want{}, v2_want{};

    prepare_for_iterations(input, state, single_iter, v_want, v1_want, v2_want, var, dvar, var1, var2);

    converged = false;

    // Only take a single step if we're coming in with both rho and T;
    // in this call we just want the EOS to fill the other thermodynamic quantities.

    if (input == eos_input_rt) converged = true;

    // Iterate until converged.

    for (int iter = 1; iter <= max_newton; ++iter) {

        // Radiation must come first since it initializes the
        // state instead of adding to it.

        apply_radiation(state, params);

        apply_ions(state);

        apply_electrons(state, tab);

        if (params->helmholtz_do_coulomb) {
            apply_coulomb_corrections(state);
        }

        // has_enthalpy false: no h update

        if (converged) {
            break;
        }
        else if (single_iter) {
            single_iter_update(state, var, dvar, v_want, converged, params);
        }
        else {
            double_iter_update(state, var1, var2, v1_want, v2_want, converged, params);
        }

    }

    finalize_state(input, state, v_want, v1_want, v2_want, params);
}

// ---------------------------------------------------------------------------
// interfaces/eos_override.H, interfaces/eos.H
// ---------------------------------------------------------------------------

// from interfaces/eos_override.H:7 eos_override
template <typename T>
MPL_DEVICE void eos_override(T* state)
{
    (void) state;
}

// from interfaces/eos.H:70 eos_reset
template <typename T>
MPL_DEVICE void eos_reset(T* state, bool& has_been_reset, const EosParams* params, EosTables tab)
{
    state->T = math::min(params->EOSData_maxtemp, math::max(params->EOSData_mintemp, state->T));
    state->rho = math::min(params->EOSData_maxdens, math::max(params->EOSData_mindens, state->rho));

    actual_eos(eos_input_rt, state, params, tab);

    has_been_reset = true;
}

// from interfaces/eos.H:83 reset_rho
template <typename T>
MPL_DEVICE void reset_rho(T* state, const EosParams* params)
{
    state->rho = math::min(params->EOSData_maxdens, math::max(params->EOSData_mindens, state->rho));
}

// from interfaces/eos.H:91 reset_T
template <typename T>
MPL_DEVICE void reset_T(T* state, const EosParams* params)
{
    state->T = math::min(params->EOSData_maxtemp, math::max(params->EOSData_mintemp, state->T));
}

// from interfaces/eos.H:98 reset_e (has_energy)
template <typename T>
MPL_DEVICE void reset_e(T* state, bool& has_been_reset, const EosParams* params, EosTables tab)
{
    if (state->e < params->EOSData_mine || state->e > params->EOSData_maxe) {
        eos_reset(state, has_been_reset, params, tab);
    }
}

// interfaces/eos.H:109 reset_h, :120 reset_s, :131 reset_p have empty
// bodies for this trait set (has_enthalpy, has_entropy, has_pressure false)
// and are not translated; their call sites in reset_inputs are marked.

// from interfaces/eos.H:142 reset_inputs (runtime switch -> if/else chain)
template <typename T>
MPL_DEVICE void reset_inputs(int input, T* state, bool& has_been_reset,
                             const EosParams* params, EosTables tab)
{
    // Reset the input quantities to valid values. For inputs other than rho and T,
    // this will evolve an EOS call, which will negate the need to do the main EOS call.

    if (input == eos_input_rt) {
        reset_rho(state, params);
        reset_T(state, params);
    }
    else if (input == eos_input_rh) {
        reset_rho(state, params);
        // reset_h: no-op (has_enthalpy false)
    }
    else if (input == eos_input_tp) {
        reset_T(state, params);
        // reset_p: no-op (has_pressure false)
    }
    else if (input == eos_input_rp) {
        reset_rho(state, params);
        // reset_p: no-op (has_pressure false)
    }
    else if (input == eos_input_re) {
        reset_rho(state, params);
        reset_e(state, has_been_reset, params, tab);
    }
    else if (input == eos_input_ps) {
        // reset_p, reset_s: no-ops (has_pressure, has_entropy false)
    }
    else if (input == eos_input_ph) {
        // reset_p, reset_h: no-ops (has_pressure, has_enthalpy false)
    }
    else if (input == eos_input_th) {
        reset_T(state, params);
        // reset_h: no-op (has_enthalpy false)
    }
}

// from interfaces/eos.H:396 eos
// (the #ifndef AMREX_USE_GPU initialized check is dropped)
template <typename T>
MPL_DEVICE void eos(int input, T* state, const EosParams* params, EosTables tab)
{
    bool has_been_reset = false;

    // has_xn: Get abar, zbar, etc.
    composition(state);

    // Force the inputs to be valid.
    reset_inputs(input, state, has_been_reset, params, tab);

    // Allow the user to override any details of the
    // EOS state. This should generally occur right
    // before the actual_eos call.
    eos_override(state);

    // Call the EOS.

    if (!has_been_reset) {
        actual_eos(input, state, params, tab);
    }
}

} // namespace mpl

#endif
