#ifndef ACTUAL_EOS_H
#define ACTUAL_EOS_H

#include <string>
#include <string_view>
#include <iostream>
#include <fstream>
#include <sstream>
#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Loop.H>

#include <extern_parameters.H>
#include <fundamental_constants.H>
#include <eos_type.H>
#include <eos_data.H>
#include <actual_eos_data.H>
#include <cmath>
#include <vector>
#if defined(USE_NEW_HELMHOLTZ_TABLE)
#include <mp_math.H>
#endif

// Frank Timmes Helmholtz based Equation of State
// https://cococubed.com/

// given a temperature temp [K], density den [g/cm**3], and a composition
// characterized by abar and zbar, this routine returns most of the other
// thermodynamic quantities. of prime interest is the pressure [erg/cm**3],
// specific thermal energy [erg/gr], the entropy [erg/g/K], along with
// their derivatives with respect to temperature, density, abar, and zbar.
// other quantities such the normalized chemical potential eta (plus its
// derivatives), number density of electrons and positron pair (along
// with their derivatives), adiabatic indices, specific heats, and
// relativistically correct sound speed are also returned.
//
// this routine assumes planckian photons, an ideal gas of ions,
// and an electron-positron gas with an arbitrary degree of relativity
// and degeneracy. interpolation in a table of the helmholtz free energy
// is used to return the electron-positron thermodynamic quantities.
// all other derivatives are analytic.
//
// references: cox & giuli chapter 24 ; timmes & swesty apj 1999

constexpr std::string_view eos_name = "helmholtz";

// quintic hermite polynomial functions
// psi0 and its derivatives
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real psi0 (amrex::Real z) noexcept
{
    return z * z * z * (z * (-6.0e0_rt * z + 15.0e0_rt) -10.0e0_rt) + 1.0e0_rt;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real dpsi0 (amrex::Real z) noexcept
{
    return z * z * (z * (-30.0e0_rt*z + 60.0e0_rt) - 30.0e0_rt);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real ddpsi0 (amrex::Real z) noexcept
{
    return z * (z * (-120.0e0_rt * z + 180.0e0_rt) - 60.0e0_rt);
}

// psi1 and its derivatives
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real psi1 (amrex::Real z) noexcept
{
    return z * ( z * z * (z * (-3.0e0_rt * z + 8.0e0_rt) - 6.0e0_rt) + 1.0e0_rt);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real dpsi1 (amrex::Real z) noexcept
{
    return z * z * ( z * (-15.0e0_rt * z + 32.0e0_rt) - 18.0e0_rt) + 1.0e0_rt;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real ddpsi1 (amrex::Real z) noexcept
{
    return z * (z * (-60.0e0_rt * z + 96.0e0_rt) - 36.0e0_rt);
}

// psi2 and its derivatives
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real psi2 (amrex::Real z) noexcept
{
    return 0.5e0_rt * z * z * ( z * ( z * (-z + 3.0e0_rt) - 3.0e0_rt) + 1.0e0_rt);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real dpsi2 (amrex::Real z) noexcept
{
    return 0.5e0_rt * z * (z * (z * (-5.0e0_rt * z + 12.0e0_rt) - 9.0e0_rt) + 2.0e0_rt);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real ddpsi2 (amrex::Real z) noexcept
{
    return 0.5e0_rt * (z * ( z * (-20.0e0_rt * z + 36.0e0_rt) - 18.0e0_rt) + 2.0e0_rt);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void fwt (const amrex::Real* AMREX_RESTRICT fi,
          const amrex::Real* AMREX_RESTRICT wt,
          amrex::Real* AMREX_RESTRICT fwtr) noexcept
{

    fwtr[0] = fi[ 0]*wt[0] + fi[ 1]*wt[1] + fi[ 2]*wt[2] + fi[18]*wt[3] + fi[19]*wt[4] + fi[20]*wt[5];
    fwtr[1] = fi[ 3]*wt[0] + fi[ 5]*wt[1] + fi[ 7]*wt[2] + fi[21]*wt[3] + fi[23]*wt[4] + fi[25]*wt[5];
    fwtr[2] = fi[ 4]*wt[0] + fi[ 6]*wt[1] + fi[ 8]*wt[2] + fi[22]*wt[3] + fi[24]*wt[4] + fi[26]*wt[5];
    fwtr[3] = fi[ 9]*wt[0] + fi[10]*wt[1] + fi[11]*wt[2] + fi[27]*wt[3] + fi[28]*wt[4] + fi[29]*wt[5];
    fwtr[4] = fi[12]*wt[0] + fi[14]*wt[1] + fi[16]*wt[2] + fi[30]*wt[3] + fi[32]*wt[4] + fi[34]*wt[5];
    fwtr[5] = fi[13]*wt[0] + fi[15]*wt[1] + fi[17]*wt[2] + fi[31]*wt[3] + fi[33]*wt[4] + fi[35]*wt[5];
}

// cubic hermite polynomial functions
// psi0 & derivatives
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real xpsi0(amrex::Real z) noexcept
{
    return z * z * (2.0e0_rt * z - 3.0e0_rt) + 1.0_rt;
}

// psi1 & derivatives
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real xpsi1(amrex::Real z) noexcept
{
    return z * (z * (z - 2.0e0_rt) + 1.0e0_rt);
}


namespace electron_table_indexing {

    struct Term {
        int dj;   // 0 or 1
        int di;   // 0 or 1
        int n;    // 0..3
    };

    constexpr std::array<Term, 16> map = {{
            {0,0,0}, {0,0,1}, {1,0,0}, {1,0,1},
            {0,0,2}, {0,0,3}, {1,0,2}, {1,0,3},
            {0,1,0}, {0,1,1}, {1,1,0}, {1,1,1},
            {0,1,2}, {0,1,3}, {1,1,2}, {1,1,3}}};

}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void apply_electrons (T& state)
{
    using namespace helmholtz;

    // assume complete ionization
    [[maybe_unused]] amrex::Real ytot1 = 1.0e0_rt / state.abar;

    // enter the table with ye*den
    amrex::Real din = state.y_e * state.rho;

    // hash locate this temperature and density
#if defined(USE_NEW_HELMHOLTZ_TABLE)
    // we use the fastmath variant of log
    int jat = int((mp::fastlg2(state.T) - tlo) * tstpi) + 1;
    int iat = int((mp::fastlg2(din) - dlo) * dstpi) + 1;
#else
    int jat = int((std::log10(state.T) - tlo) * tstpi) + 1;
    int iat = int((std::log10(din) - dlo) * dstpi) + 1;
#endif
    jat = amrex::Clamp(jat, 1, jmax-1) - 1;
    iat = amrex::Clamp(iat, 1, imax-1) - 1;

    amrex::Real fi[36];

    // access the table locations only once
    for (int i = 0; i < 9; ++i) {
        fi[i     ] = f[jat  ][iat  ][i]; // f, ft, ftt, fd, fdd, fdt, fddt, fdtt, fddtt
        fi[i +  9] = f[jat  ][iat+1][i];
        fi[i + 18] = f[jat+1][iat  ][i];
        fi[i + 27] = f[jat+1][iat+1][i];
    }

    // various differences
    amrex::Real xt  = amrex::max((state.T - t[jat]) * dti_sav[jat], 0.0e0_rt);
    amrex::Real xd  = amrex::max((din - d[iat]) * ddi_sav[iat], 0.0e0_rt);
    amrex::Real mxt = 1.0e0_rt - xt;
    amrex::Real mxd = 1.0e0_rt - xd;

    // the six density basis functions
    amrex::Real sid[6];

    sid[0] =  psi0(xd);
    sid[1] =  psi1(xd) * dd_sav[iat];
    sid[2] =  psi2(xd) * dd2_sav[iat];

    sid[3] =  psi0(mxd);
    sid[4] = -psi1(mxd) * dd_sav[iat];
    sid[5] =  psi2(mxd) * dd2_sav[iat];

    // derivatives of the weight functions
    amrex::Real dsid[6];

    dsid[0] =  dpsi0(xd) * ddi_sav[iat];
    dsid[1] =  dpsi1(xd);
    dsid[2] =  dpsi2(xd) * dd_sav[iat];

    dsid[3] = -dpsi0(mxd) * ddi_sav[iat];
    dsid[4] =  dpsi1(mxd);
    dsid[5] = -dpsi2(mxd) * dd_sav[iat];

    // This array saves some subexpressions that go into
    // computing the biquintic polynomial. Instead of explicitly
    // constructing it in full, we'll use these subexpressions
    // and then compute the result as
    // (table data * temperature terms) * density terms.

    amrex::Real fwtr[6];
    amrex::Real sit[6];

    {
        // temperature basis functions

        sit[0] = psi0(xt);
        sit[1] = psi1(xt) * dt_sav[jat];
        sit[2] = psi2(xt) * dt2_sav[jat];

        sit[3] =  psi0(mxt);
        sit[4] = -psi1(mxt) * dt_sav[jat];
        sit[5] =  psi2(mxt) * dt2_sav[jat];

        fwt(fi, sit, fwtr);
    }

    amrex::Real free = 0.e0_rt;
    amrex::Real df_d = 0.e0_rt;

    for (int i = 0; i <= 5; ++i) {
        // the free energy
        free = free + fwtr[i] * sid[i];

        // derivative with respect to density
        df_d = df_d + fwtr[i] * dsid[i];
    }

    {
        // temperature derivative of weight functions
        // this was originally called dsit
        sit[0] =  dpsi0(xt) * dti_sav[jat];
        sit[1] =  dpsi1(xt);
        sit[2] =  dpsi2(xt) * dt_sav[jat];

        sit[3] = -dpsi0(mxt) * dti_sav[jat];
        sit[4] =  dpsi1(mxt);
        sit[5] = -dpsi2(mxt) * dt_sav[jat];

        fwt(fi, sit, fwtr);
    }

    amrex::Real df_t = 0.e0_rt;
    amrex::Real df_dt = 0.e0_rt;

    for (int i = 0; i <= 5; ++i) {
        // derivative with respect to temperature
        df_t += fwtr[i] * sid[i];

        // derivative with respect to temperature and density
        df_dt += fwtr[i] * dsid[i];
    }

    {
        // second derivatives of the weight functions
        // this was originally called ddsit

        sit[0] =  ddpsi0(xt) * dt2i_sav[jat];
        sit[1] =  ddpsi1(xt) * dti_sav[jat];
        sit[2] =  ddpsi2(xt);

        sit[3] =  ddpsi0(mxt) * dt2i_sav[jat];
        sit[4] = -ddpsi1(mxt) * dti_sav[jat];
        sit[5] =  ddpsi2(mxt);

        fwt(fi, sit, fwtr);
    }

    amrex::Real df_tt = 0.e0_rt;
    for (int i = 0; i <= 5; ++i) {
        // derivative with respect to temperature**2
        df_tt = df_tt + fwtr[i] * sid[i];
    }

    // now get the pressure derivative with density, chemical potential, and
    // electron positron number densities
    // get the interpolation weight functions
    sit[0] = xpsi0(xt);
    sit[1] = xpsi1(xt) * dt_sav[jat];

    sit[2] = xpsi0(mxt);
    sit[3] = -xpsi1(mxt) * dt_sav[jat];

    sid[0] = xpsi0(xd);
    sid[1] = xpsi1(xd) * dd_sav[iat];

    sid[2] = xpsi0(mxd);
    sid[3] = -xpsi1(mxd) * dd_sav[iat];


    // Reuse subexpressions that would go into computing the
    // cubic interpolation.
    amrex::Real wdt[16];

    for (int i = 0; i <= 3; ++i) {
        wdt[i     ] = sid[0] * sit[i];
        wdt[i +  4] = sid[1] * sit[i];
        wdt[i +  8] = sid[2] * sit[i];
        wdt[i + 12] = sid[3] * sit[i];
    }

    // Read in the tabular data for the pressure derivatives.
    // We have some freedom in how we store it in the local
    // array. We choose here to index it such that we can
    // immediately evaluate the cubic interpolant below as
    // fi * wdt, which ensures that we have the right combination
    // of grid points and derivatives at grid points to evaluate
    // the interpolation correctly. This is handled by map,
    // giving the j offset, i offset, and component

    // pressure derivative with density
    [[maybe_unused]] amrex::Real dpepdd = 0.0e0_rt;
    if constexpr (has_pressure<T>) {
        amrex::constexpr_for<0, 16>([&] (auto K) {
            constexpr int k = K;
            constexpr auto m = electron_table_indexing::map[k];

            dpepdd += dpdf[jat+m.dj][iat+m.di][m.n] * wdt[k];
        });
        dpepdd = amrex::max(state.y_e * dpepdd, 0.0e0_rt);
    }

    // electron chemical potential etaele
    [[maybe_unused]] amrex::Real etaele{};
    if constexpr (has_eta<T>) {
        amrex::constexpr_for<0, 16>([&] (auto K) {
            constexpr int k = K;
            constexpr auto m = electron_table_indexing::map[k];

            etaele += ef[jat+m.dj][iat+m.di][m.n] * wdt[k];
        });
    }

    // electron + positron number densities
    [[maybe_unused]] amrex::Real xnefer{};
    if constexpr (has_xne_xnp<T>) {
        amrex::constexpr_for<0, 16>([&] (auto K) {
            constexpr int k = K;
            constexpr auto m = electron_table_indexing::map[k];

            xnefer += xf[jat+m.dj][iat+m.di][m.n] * wdt[k];
        });
    }

    // the desired electron-positron thermodynamic quantities

    // dpepdd at high temperatures and low densities is below the
    // floating point limit of the subtraction of two large terms.
    // since state.dpdr doesn't enter the maxwell relations at all, use the
    // bicubic interpolation done above instead of this one
    amrex::Real x       = din * din;
    amrex::Real pele    = x * df_d;
    amrex::Real dpepdt  = x * df_dt;
    [[maybe_unused]] amrex::Real dpepda{};
    [[maybe_unused]] amrex::Real dpepdz{};
    if constexpr (has_pressure<T>) {
        if constexpr (has_dpdA<T> || has_dpdZ<T>) {
            amrex::Real s       = dpepdd/state.y_e - 2.0e0_rt * din * df_d;
            if constexpr (has_dpdA<T>) {
                dpepda  = -ytot1 * (2.0e0_rt * pele + s * din);
            }
            if constexpr (has_dpdZ<T>) {
                dpepdz  = state.rho*ytot1*(2.0e0_rt * din * df_d  +  s);
            }
        }
    }

    x            = state.y_e * state.y_e;
    amrex::Real sele    = -df_t * state.y_e;
    amrex::Real dsepdt  = -df_tt * state.y_e;
    amrex::Real dsepdd  = -df_dt * x;

    amrex::Real eele    = state.y_e*free + state.T * sele;
    amrex::Real deepdt  = state.T * dsepdt;
    amrex::Real deepdd  = x * df_d + state.T * dsepdd;

    [[maybe_unused]] amrex::Real deepda{};
    if constexpr (has_dedA<T>) {
        amrex::Real dsepda = ytot1 * (state.y_e * df_dt * din - sele);
        deepda = -state.y_e * ytot1 * (free +  df_d * din) + state.T * dsepda;
    }

    [[maybe_unused]] amrex::Real deepdz{};
    if constexpr (has_dedZ<T>) {
        amrex::Real dsepdz = -ytot1 * (state.y_e * df_dt * state.rho  + df_t);
        deepdz = ytot1* (free + state.y_e * df_d * state.rho) + state.T * dsepdz;
    }

    if constexpr (has_pressure<T>) {
        state.p    = state.p + pele;
        state.dpdT = state.dpdT + dpepdt;
        state.dpdr = state.dpdr + dpepdd;
        if constexpr (has_dpdA<T>) {
            state.dpdA = state.dpdA + dpepda;
        }
        if constexpr (has_dpdZ<T>) {
            state.dpdZ = state.dpdZ + dpepdz;
        }
    }

    if constexpr (has_entropy<T>) {
        state.s    = state.s + sele;
        state.dsdT = state.dsdT + dsepdt;
        state.dsdr = state.dsdr + dsepdd;
    }

    if constexpr (has_energy<T>) {
        state.e    = state.e + eele;
        state.dedT = state.dedT + deepdt;
        state.dedr = state.dedr + deepdd;
        if constexpr (has_dedA<T>) {
            state.dedA = state.dedA + deepda;
        }
        if constexpr (has_dedZ<T>) {
            state.dedZ = state.dedZ + deepdz;
        }
    }

    if constexpr (has_eta<T>) {
        state.eta = etaele;
    }

    if constexpr (has_xne_xnp<T>) {
        state.xne = xnefer;
        state.xnp = 0.0e0_rt;
    }

    if constexpr (has_pele_ppos<T>) {
        state.pele = pele;
        state.ppos = 0.0e0_rt;
    }
}



template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void apply_ions (T& state)
{
    using namespace helmholtz;

    constexpr amrex::Real pi = std::numbers::pi_v<amrex::Real>;
    constexpr amrex::Real sioncon = (2.0e0_rt * pi * amu * kerg)/(h*h);
    constexpr amrex::Real kergavo = kerg * avo_eos;

    amrex::Real deni = 1.0e0_rt / state.rho;
    amrex::Real tempi = 1.0e0_rt / state.T;

    amrex::Real ytot1   = 1.0e0_rt / state.abar;
    amrex::Real xni     = avo_eos * ytot1 * state.rho;
    amrex::Real dxnidd  = avo_eos * ytot1;
    [[maybe_unused]] amrex::Real dxnida  = -xni * ytot1;

    amrex::Real kt = kerg * state.T;

    amrex::Real pion    = xni * kt;
    amrex::Real dpiondd = dxnidd * kt;
    amrex::Real dpiondt = xni * kerg;
    [[maybe_unused]] amrex::Real dpionda = dxnida * kt;
    [[maybe_unused]] amrex::Real dpiondz = 0.0e0_rt;

    amrex::Real eion    = 1.5e0_rt * pion * deni;
    amrex::Real deiondd = (1.5e0_rt * dpiondd - eion) * deni;
    amrex::Real deiondt = 1.5e0_rt * dpiondt * deni;
    [[maybe_unused]] amrex::Real deionda = 1.5e0_rt * dpionda * deni;
    [[maybe_unused]] amrex::Real deiondz = 0.0e0_rt;

    [[maybe_unused]] amrex::Real sion{};
    [[maybe_unused]] amrex::Real dsiondd{};
    [[maybe_unused]] amrex::Real dsiondt{};
    if constexpr (has_entropy<T>) {
        amrex::Real x       = state.abar * state.abar * std::sqrt(state.abar) * deni / avo_eos;
        amrex::Real s       = sioncon * state.T;
        amrex::Real z       = x * s * std::sqrt(s);
        amrex::Real y       = std::log(z);
        sion    = (pion * deni + eion) * tempi + kergavo * ytot1 * y;
        dsiondd = (dpiondd * deni - pion * deni * deni + deiondd) * tempi -
            kergavo * deni * ytot1;
        dsiondt = (dpiondt * deni + deiondt) * tempi -
            (pion * deni + eion) * tempi * tempi +
            1.5e0_rt * kergavo * tempi * ytot1;
    }

    if constexpr (has_pressure<T>) {
        state.p    = state.p + pion;
        state.dpdT = state.dpdT + dpiondt;
        state.dpdr = state.dpdr + dpiondd;
        if constexpr (has_dpdA<T>) {
            state.dpdA = state.dpdA + dpionda;
        }
        if constexpr (has_dpdZ<T>) {
            state.dpdZ = state.dpdZ + dpiondz;
        }
    }

    if constexpr (has_energy<T>) {
        state.e    = state.e + eion;
        state.dedT = state.dedT + deiondt;
        state.dedr = state.dedr + deiondd;
        if constexpr (has_dedA<T>) {
            state.dedA = state.dedA + deionda;
        }
        if constexpr (has_dedZ<T>) {
            state.dedZ = state.dedZ + deiondz;
        }
    }

    if constexpr (has_entropy<T>) {
        state.s    = state.s + sion;
        state.dsdT = state.dsdT + dsiondt;
        state.dsdr = state.dsdr + dsiondd;
    }
}



template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void apply_radiation (T& state)
{
    using namespace helmholtz;

    constexpr amrex::Real clight  = 2.99792458e10_rt;
#ifdef RADIATION
    constexpr amrex::Real ssol    = 0.0e0_rt;
#else
    constexpr amrex::Real ssol    = 5.67051e-5_rt;
#endif
    constexpr amrex::Real asol    = 4.0e0_rt * ssol / clight;
    constexpr amrex::Real asoli3  = asol/3.0e0_rt;

    amrex::Real deni = 1.0e0_rt / state.rho;
    amrex::Real tempi = 1.0e0_rt / state.T;

    amrex::Real prad = asoli3 * state.T * state.T * state.T * state.T;

    // In low density material, this radiation pressure becomes unphysically high.
    // For rho ~ 1 g/cc and T ~ 1e9, then this radiation pressure will lead to a
    // sound speed ~ 0.1c, and this is a problem because the codes using this EOS
    // are primarily non-relativistic. Thus we can optionally smooth out the radiation
    // pressure such that it disappears at low densities. Since all terms below depend
    // on prad, this will result in the radiation term effectively vanishing below the
    // cutoff density. For simplicity we ignore the effect this has on the derivatives.

    if (eos_rp::prad_limiter_rho_c > 0.0e0_rt && eos_rp::prad_limiter_delta_rho > 0.0e0_rt) {
        prad = prad * 0.5e0_rt * (1.0e0_rt + std::tanh((state.rho - eos_rp::prad_limiter_rho_c) / eos_rp::prad_limiter_delta_rho));
    }

    amrex::Real dpraddd = 0.0e0_rt;
    amrex::Real dpraddt = 4.0e0_rt * prad * tempi;

    amrex::Real erad    = 3.0e0_rt * prad*deni;
    amrex::Real deraddd = -erad * deni;
    amrex::Real deraddt = 3.0e0_rt * dpraddt * deni;

    [[maybe_unused]] amrex::Real srad{};
    [[maybe_unused]] amrex::Real dsraddd{};
    [[maybe_unused]] amrex::Real dsraddt{};
    if constexpr (has_entropy<T>) {
        srad    = (prad * deni + erad) * tempi;
        dsraddd = (dpraddd * deni - prad * deni * deni + deraddd) * tempi;
        dsraddt = (dpraddt * deni + deraddt - srad) * tempi;
    }

    // Note that unlike the other terms, radiation
    // sets these terms instead of adding to them,
    // since it comes first.

    if constexpr (has_pressure<T>) {
        state.p    = prad;
        state.dpdr = dpraddd;
        state.dpdT = dpraddt;
        if constexpr (has_dpdA<T>) {
            state.dpdA = 0.0_rt;
        }
        if constexpr (has_dpdZ<T>) {
            state.dpdZ = 0.0_rt;
        }
    }

    if constexpr (has_energy<T>) {
        state.e    = erad;
        state.dedr = deraddd;
        state.dedT = deraddt;
        if constexpr (has_dedA<T>) {
            state.dedA = 0.0_rt;
        }
        if constexpr (has_dedZ<T>) {
            state.dedZ = 0.0_rt;
        }
    }

    if constexpr (has_entropy<T>) {
        state.s    = srad;
        state.dsdr = dsraddd;
        state.dsdT = dsraddt;
    }
}



template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void apply_coulomb_corrections (T& state)
{
    using namespace helmholtz;

    // Constants used for the Coulomb corrections
    constexpr amrex::Real a1 = -0.898004e0_rt;
    constexpr amrex::Real b1 =  0.96786e0_rt;
    constexpr amrex::Real c1 =  0.220703e0_rt;
    constexpr amrex::Real d1 = -0.86097e0_rt;
    constexpr amrex::Real e1 =  2.5269e0_rt;
    constexpr amrex::Real a2 =  0.29561e0_rt;
    constexpr amrex::Real b2 =  1.9885e0_rt;
    constexpr amrex::Real c2 =  0.288675e0_rt;
    constexpr amrex::Real qe   = 4.8032042712e-10_rt;
    constexpr amrex::Real esqu = qe * qe;
    constexpr amrex::Real onethird = 1.0e0_rt/3.0e0_rt;
    constexpr amrex::Real forth = 4.0e0_rt/3.0e0_rt;
    constexpr amrex::Real pi = std::numbers::pi_v<amrex::Real>;

    [[maybe_unused]] amrex::Real pcoul    = 0.e0_rt;
    [[maybe_unused]] amrex::Real dpcouldd = 0.e0_rt;
    [[maybe_unused]] amrex::Real dpcouldt = 0.e0_rt;
    [[maybe_unused]] amrex::Real ecoul    = 0.e0_rt;
    [[maybe_unused]] amrex::Real decouldd = 0.e0_rt;
    [[maybe_unused]] amrex::Real decouldt = 0.e0_rt;
    [[maybe_unused]] amrex::Real scoul    = 0.e0_rt;
    [[maybe_unused]] amrex::Real dscouldd = 0.e0_rt;
    [[maybe_unused]] amrex::Real dscouldt = 0.e0_rt;

    [[maybe_unused]] amrex::Real dpcoulda = 0.e0_rt;
    [[maybe_unused]] amrex::Real dpcouldz = 0.e0_rt;
    [[maybe_unused]] amrex::Real decoulda = 0.e0_rt;
    [[maybe_unused]] amrex::Real decouldz = 0.e0_rt;

    amrex::Real s, z; // temporary variables

    // uniform background corrections only
    // from yakovlev & shalybkov 1989
    // lami is the average ion separation
    // plasg is the plasma coupling parameter

    amrex::Real ytot1 = 1.0e0_rt / state.abar;
    amrex::Real xni     = avo_eos * ytot1 * state.rho;
    amrex::Real dxnidd  = avo_eos * ytot1;
    [[maybe_unused]] amrex::Real dxnida  = -xni * ytot1;

    amrex::Real kt      = kerg * state.T;
    amrex::Real ktinv   = 1.0e0_rt / kt;

    z             = forth * pi;
    s             = z * xni;
    amrex::Real dsdd     = z * dxnidd;
    [[maybe_unused]] amrex::Real dsda     = z * dxnida;

    amrex::Real inv_lami = std::cbrt(s);
    amrex::Real lami     = 1.0e0_rt / inv_lami;
    z             = -onethird * lami;
    amrex::Real lamidd   = z * dsdd / s;
    [[maybe_unused]] amrex::Real lamida   = z * dsda / s;

    amrex::Real plasg    = state.zbar * state.zbar * esqu * ktinv * inv_lami;
    z             = -plasg * inv_lami;
    amrex::Real plasgdd  = z * lamidd;
    amrex::Real plasgdt  = -plasg*ktinv * kerg;
    [[maybe_unused]] amrex::Real plasgda  = z * lamida;
    [[maybe_unused]] amrex::Real plasgdz  = 2.0e0_rt * plasg/state.zbar;

    // .yakovlev & shalybkov 1989 equations 82, 85, 86, 87
    if (plasg >= 1.0e0_rt)
    {
        amrex::Real x = std::pow(plasg, 0.25e0_rt);
        amrex::Real y = avo_eos * ytot1 * kerg;
        ecoul    = y * state.T * (a1 * plasg + b1 * x + c1 / x + d1);
        pcoul    = onethird * state.rho * ecoul;

        if constexpr (has_entropy<T>) {
            scoul  = -y * (3.0e0_rt * b1 * x - 5.0e0_rt*c1 / x +
                           d1 * (std::log(plasg) - 1.0e0_rt) - e1);
        }

        y        = avo_eos*ytot1*kt*(a1 + 0.25e0_rt/plasg*(b1*x - c1/x));
        decouldd = y * plasgdd;
        decouldt = y * plasgdt + ecoul/state.T;

        decoulda = y * plasgda - ecoul * ytot1;
        decouldz = y * plasgdz;

        y        = onethird * state.rho;
        dpcouldd = onethird * ecoul + y * decouldd;
        dpcouldt = y * decouldt;

        dpcoulda = y * decoulda;
        dpcouldz = y * decouldz;

       if constexpr (has_entropy<T>) {
           y = -avo_eos * kerg / (state.abar * plasg) *
               (0.75e0_rt * b1 * x + 1.25e0_rt * c1 / x + d1);
           dscouldd = y * plasgdd;
           dscouldt = y * plasgdt;
       }

       // yakovlev & shalybkov 1989 equations 102, 103, 104
    }
    else if (plasg < 1.0e0_rt)
    {

        amrex::Real pion    = xni * kt;
        amrex::Real dpiondd = dxnidd * kt;
        amrex::Real dpiondt = xni * kerg;
        amrex::Real dpionda = dxnida * kt;
        amrex::Real dpiondz = 0.0e0_rt;

        amrex::Real x = plasg * std::sqrt(plasg);
        amrex::Real y = std::pow(plasg, b2);
        z        = c2 * x - onethird * a2 * y;
        pcoul    = -pion * z;
        ecoul    = 3.0e0_rt * pcoul / state.rho;
        if constexpr (has_entropy<T>) {
            scoul = -avo_eos * ytot1 * kerg * (c2 * x -a2 * (b2 - 1.0e0_rt) / b2 * y);
        }

        s        = 1.5e0_rt * c2 * x / plasg - onethird * a2 * b2 * y / plasg;
        dpcouldd = -dpiondd * z - pion * s * plasgdd;
        dpcouldt = -dpiondt * z - pion * s * plasgdt;
        dpcoulda = -dpionda * z - pion * s * plasgda;
        dpcouldz = -dpiondz * z - pion * s * plasgdz;

        s        = 3.0e0_rt / state.rho;
        decouldd = s * dpcouldd - ecoul / state.rho;
        decouldt = s * dpcouldt;
        decoulda = s * dpcoulda;
        decouldz = s * dpcouldz;

        if constexpr (has_entropy<T>) {
            s = -avo_eos * kerg / (state.abar * plasg) *
                (1.5e0_rt * c2 * x - a2 * (b2 - 1.0e0_rt) * y);
            dscouldd = s * plasgdd;
            dscouldt = s * plasgdt;
        }
    }

    // Disable Coulomb corrections if they cause
    // the energy or pressure to go negative.

    amrex::Real p_temp = std::numeric_limits<amrex::Real>::max();
    amrex::Real e_temp = std::numeric_limits<amrex::Real>::max();

    if constexpr (has_pressure<T>) {
        p_temp = state.p + pcoul;
    }
    if constexpr (has_energy<T>) {
        e_temp = state.e + ecoul;
    }

    if (p_temp <= 0.0e0_rt || e_temp <= 0.0e0_rt)
    {
        pcoul    = 0.0e0_rt;
        dpcouldd = 0.0e0_rt;
        dpcouldt = 0.0e0_rt;
        ecoul    = 0.0e0_rt;
        decouldd = 0.0e0_rt;
        decouldt = 0.0e0_rt;
        if constexpr (has_entropy<T>) {
            scoul    = 0.0e0_rt;
            dscouldd = 0.0e0_rt;
            dscouldt = 0.0e0_rt;
        }

        dpcoulda = 0.0e0_rt;
        dpcouldz = 0.0e0_rt;
        decoulda = 0.0e0_rt;
        decouldz = 0.0e0_rt;
    }

    if constexpr (has_pressure<T>) {
        state.p    = state.p + pcoul;
        state.dpdr = state.dpdr + dpcouldd;
        state.dpdT = state.dpdT + dpcouldt;
        if constexpr (has_dpdA<T>) {
            state.dpdA = state.dpdA + dpcoulda;
        }
        if constexpr (has_dpdZ<T>) {
            state.dpdZ = state.dpdZ + dpcouldz;
        }
    }

    if constexpr (has_energy<T>) {
        state.e    = state.e + ecoul;
        state.dedr = state.dedr + decouldd;
        state.dedT = state.dedT + decouldt;
        if constexpr (has_dedA<T>) {
            state.dedA = state.dedA + decoulda;
        }
        if constexpr (has_dedZ<T>) {
            state.dedZ = state.dedZ + decouldz;
        }
    }

    if constexpr (has_entropy<T>) {
        state.s    = state.s + scoul;
        state.dsdr = state.dsdr + dscouldd;
        state.dsdT = state.dsdT + dscouldt;
    }
}



template <typename I, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void prepare_for_iterations (I input, const T& state,
                             bool& single_iter, amrex::Real& v_want, amrex::Real& v1_want, amrex::Real& v2_want,
                             int& var, int& dvar, int& var1, int& var2)
{
    using namespace helmholtz;
    using namespace EOS;

    single_iter = true;

    if (input == eos_input_rt) {

        // Nothing to do here.

    }
    else if (input == eos_input_rh) {

        if constexpr (has_enthalpy<T>) {
            v_want = state.h;
            var  = ienth;
            dvar = itemp;
        }

    }
    else if (input == eos_input_tp) {

        if constexpr (has_pressure<T>) {
            v_want = state.p;
            var  = ipres;
            dvar = idens;
        }

    }
    else if (input == eos_input_rp) {

        if constexpr (has_pressure<T>) {
            v_want = state.p;
            var  = ipres;
            dvar = itemp;
        }

    }
    else if (input == eos_input_re) {

        if constexpr (has_energy<T>) {
            v_want = state.e;
            var  = iener;
            dvar = itemp;
        }

    }
    else if (input == eos_input_ps) {

        if constexpr (has_pressure<T> && has_entropy<T>) {
            single_iter = false;
            v1_want = state.p;
            v2_want = state.s;
            var1 = ipres;
            var2 = ientr;
        }

    }
    else if (input == eos_input_ph) {

        if constexpr (has_pressure<T> && has_enthalpy<T>) {
            single_iter = false;
            v1_want = state.p;
            v2_want = state.h;
            var1 = ipres;
            var2 = ienth;
        }

    }
    else if (input == eos_input_th) {

        if constexpr (has_enthalpy<T>) {
            v_want = state.h;
            var  = ienth;
            dvar = idens;
        }

    }
#ifndef AMREX_USE_GPU
    else {

        amrex::Error("Unknown EOS input");

    }
#endif
}



template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void single_iter_update (T& state, int var, int dvar,
                         amrex::Real v_want, bool& converged)
{

    using namespace helmholtz;
    using namespace EOS;

    amrex::Real x, v = 0.0_rt, dvdx = 0.0_rt, xtol, smallx;

    if (dvar == itemp) {

        x = state.T;
        smallx = EOSData::mintemp;
        xtol = ttol;

        if (var == ipres) {
            if constexpr (has_pressure<T>) {
                v    = state.p;
                dvdx = state.dpdT;
            }
        }
        else if (var == iener) {
            if constexpr (has_energy<T>) {
                v    = state.e;
                dvdx = state.dedT;
            }
        }
        else if (var == ientr) {
            if constexpr (has_entropy<T>) {
                v    = state.s;
                dvdx = state.dsdT;
            }
        }
        else if (var == ienth) {
            if constexpr (has_enthalpy<T>) {
                v    = state.h;
                dvdx = state.dhdT;
            }
        }

    }
    else // dvar == density
    {

        x = state.rho;
        smallx = EOSData::mindens;
        xtol = dtol;

        if (var == ipres) {
            if constexpr (has_pressure<T>) {
                v    = state.p;
                dvdx = state.dpdr;
            }
        }
        else if (var == iener) {
            if constexpr (has_energy<T>) {
                v    = state.e;
                dvdx = state.dedr;
            }
        }
        else if (var == ientr) {
            if constexpr (has_entropy<T>) {
                v    = state.s;
                dvdx = state.dsdr;
            }
        }
        else if (var == ienth) {
            if constexpr (has_enthalpy<T>) {
                v    = state.h;
                dvdx = state.dhdr;
            }
        }

    }

    // Now do the calculation for the next guess for T/rho
    amrex::Real xnew = x - (v - v_want) / dvdx;

    // Don't let the temperature/density change by more than a factor of two
    xnew = amrex::Clamp(xnew, 0.5_rt * x, 2.0_rt * x);

    // Don't let us freeze/evacuate
    xnew = amrex::max(smallx, xnew);

    // Store the new temperature/density

    if (dvar == itemp) {
        state.T = xnew;
    }
    else
    {
        state.rho = xnew;
    }

    // Compute the error from the last iteration

    amrex::Real error = std::abs(xnew - x);

    if (error < xtol * x) converged = true;
}



template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void double_iter_update (T& state, int var1, int var2,
                         amrex::Real v1_want, amrex::Real v2_want, bool& converged)
{
    using namespace helmholtz;
    using namespace EOS;

    amrex::Real v1 = 0.0_rt, dv1dt = 0.0_rt, dv1dr = 0.0_rt, v2 = 0.0_rt, dv2dt = 0.0_rt, dv2dr = 0.0_rt;

    // Figure out which variables we're using

    amrex::Real told = state.T;
    amrex::Real rold = state.rho;

    if (var1 == ipres) {
        if constexpr (has_pressure<T>) {
            v1    = state.p;
            dv1dt = state.dpdT;
            dv1dr = state.dpdr;
        }
    }
    else if (var1 == iener) {
        if constexpr (has_energy<T>) {
            v1    = state.e;
            dv1dt = state.dedT;
            dv1dr = state.dedr;
        }
    }
    else if (var1 == ientr) {
        if constexpr (has_entropy<T>) {
            v1    = state.s;
            dv1dt = state.dsdT;
            dv1dr = state.dsdr;
        }
    }
    else if (var1 == ienth) {
        if constexpr (has_enthalpy<T>) {
            v1    = state.h;
            dv1dt = state.dhdT;
            dv1dr = state.dhdr;
        }
    }

    if (var2 == ipres) {
        if constexpr (has_pressure<T>) {
            v2    = state.p;
            dv2dt = state.dpdT;
            dv2dr = state.dpdr;
        }
    }
    else if (var2 == iener) {
        if constexpr (has_energy<T>) {
            v2    = state.e;
            dv2dt = state.dedT;
            dv2dr = state.dedr;
        }
    }
    else if (var2 == ientr) {
        if constexpr (has_entropy<T>) {
            v2    = state.s;
            dv2dt = state.dsdT;
            dv2dr = state.dsdr;
        }
    }
    else if (var2 == ienth) {
        if constexpr (has_enthalpy<T>) {
            v2    = state.h;
            dv2dt = state.dhdT;
            dv2dr = state.dhdr;
        }
    }

    // Two functions, f and g, to iterate over
    amrex::Real v1i = v1_want - v1;
    amrex::Real v2i = v2_want - v2;

    //
    // 0 = f + dfdr * delr + dfdt * delt
    // 0 = g + dgdr * delr + dgdt * delt
    //

    // note that dfi/dT = - df/dT
    amrex::Real delr = (-v1i * dv2dt + v2i * dv1dt) / (dv2dr * dv1dt - dv2dt * dv1dr);

    amrex::Real rnew = rold + delr;

    amrex::Real tnew = told + (v1i - dv1dr * delr) / dv1dt;

    // Don't let the temperature or density change by more
    // than a factor of two
    tnew = amrex::Clamp(tnew, 0.5e0_rt * told, 2.0e0_rt * told);
    rnew = amrex::Clamp(rnew, 0.5e0_rt * rold, 2.0e0_rt * rold);

    // Don't let us freeze or evacuate
    tnew = amrex::max(EOSData::mintemp, tnew);
    rnew = amrex::max(EOSData::mindens, rnew);

    // Store the new temperature and density
    state.rho = rnew;
    state.T = tnew;

    // Compute the errors
    amrex::Real error1 = std::abs(rnew - rold);
    amrex::Real error2 = std::abs(tnew - told);

    if (error1 < dtol * rold && error2 < ttol * told) converged = true;
}



template <typename I, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void finalize_state (I input, T& state,
                     amrex::Real v_want, amrex::Real v1_want, amrex::Real v2_want)
{
    using namespace helmholtz;

    // define mu -- the total mean molecular weight including both electrons and ions
    state.mu = 1.0e0_rt / (1.0e0_rt / state.abar + 1.0e0_rt / state.mu_e);

    // Calculate some remaining derivatives
    if constexpr (has_pressure<T>) {
        state.dpde = state.dpdT / state.dedT;
        state.dpdr_e = state.dpdr - state.dpdT * state.dedr / state.dedT;
    }

    // Specific heats and Gamma_1
    if constexpr (has_energy<T>) {
        state.cv = state.dedT;

        if constexpr (has_pressure<T>) {
            amrex::Real chit = state.T / state.p * state.dpdT;
            amrex::Real chid = state.dpdr * state.rho / state.p;

            state.gam1 = (chit * (state.p / state.rho)) * (chit / (state.T * state.cv)) + chid;
            state.cp = state.cv * state.gam1 / chid;
        }
    }

    // Use the non-relativistic version of the sound speed, cs = sqrt(gam_1 * P / rho).
    // This replaces the relativistic version that comes out of helmeos.
    if constexpr (has_pressure<T>) {
        state.cs = std::sqrt(state.gam1 * state.p / state.rho);
    }

    if (input_is_constant) {

       if (input == eos_input_rh) {  // NOLINT(bugprone-branch-clone)

           if constexpr (has_enthalpy<T>) {
               state.h = v_want;
           }

       }
       else if (input == eos_input_tp) {  // NOLINT(bugprone-branch-clone)

           if constexpr (has_pressure<T>) {
               state.p = v_want;
           }

       }
       else if (input == eos_input_rp) {

           if constexpr (has_pressure<T>) {
               state.p = v_want;
           }

       }
       else if (input == eos_input_re) {

           if constexpr (has_energy<T>) {
               state.e = v_want;
           }

       }
       else if (input == eos_input_ps) {

           if constexpr (has_pressure<T>) {
               state.p = v1_want;
           }
           if constexpr (has_entropy<T>) {
               state.s = v2_want;
           }

       }
       else if (input == eos_input_ph) {

           if constexpr (has_pressure<T>) {
               state.p = v1_want;
           }
           if constexpr (has_enthalpy<T>) {
               state.h = v2_want;
           }

       }
       else if (input == eos_input_th) {

           if constexpr (has_enthalpy<T>) {
               state.h = v_want;
           }

       }

    }
}



template <typename I, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_eos (I input, T& state)
{
    static_assert(std::is_same_v<I, eos_input_t>, "input must be an eos_input_t");

    using namespace helmholtz;

    constexpr int max_newton = 100;

    bool single_iter, converged;
    int var{}, dvar{}, var1, var2;
    amrex::Real v_want{}, v1_want{}, v2_want{};

    prepare_for_iterations(input, state, single_iter, v_want, v1_want, v2_want, var, dvar, var1, var2);

    converged = false;

    // Only take a single step if we're coming in with both rho and T;
    // in this call we just want the EOS to fill the other thermodynamic quantities.

    if (input == eos_input_rt) converged = true;

    // Iterate until converged.

    for (int iter = 1; iter <= max_newton; ++iter) {

        // Radiation must come first since it initializes the
        // state instead of adding to it.

        apply_radiation(state);

        apply_ions(state);

        apply_electrons(state);

        if (do_coulomb) {
            apply_coulomb_corrections(state);
        }

        // Calculate enthalpy the usual way, h = e + p / rho.

        if constexpr (has_enthalpy<T>) {
            state.h = state.e + state.p / state.rho;
            state.dhdr = state.dedr + state.dpdr / state.rho - state.p / (state.rho * state.rho);
            state.dhdT = state.dedT + state.dpdT / state.rho;
        }

        if (converged) {
            break;
        }
        else if (single_iter) {
            single_iter_update(state, var, dvar, v_want, converged);
        }
        else {
            double_iter_update(state, var1, var2, v1_want, v2_want, converged);
        }

    }

    finalize_state(input, state, v_want, v1_want, v2_want);
}



AMREX_INLINE
void actual_eos_init ()
{
    using namespace helmholtz;

    amrex::Real dth, dt2, dti, dt2i;
    amrex::Real dd, dd2, ddi, dd2i;

    // Read in the runtime parameters

    input_is_constant = eos_rp::eos_input_is_constant;
    do_coulomb = eos_rp::use_eos_coulomb;
    ttol = eos_rp::eos_ttol;
    dtol = eos_rp::eos_dtol;

    //    read the helmholtz free energy table

    for (int j = 0; j < jmax; ++j) {
        amrex::Real tsav = tlo + j * tstp;
#if defined(USE_NEW_HELMHOLTZ_TABLE)
        t[j] = mp::fastpow2(tsav);
#else
        t[j] = std::pow(10.0e0_rt, tsav);
#endif

        for (int i = 0; i < imax; ++i) {
            amrex::Real dsav = dlo + i * dstp;
#if defined(USE_NEW_HELMHOLTZ_TABLE)
            d[i] = mp::fastpow2(dsav);
#else
            d[i] = std::pow(10.0e0_rt, dsav);
#endif
        }
    }

    // it does not work on all machines (for GPUs) broadcast to other
    // procs from managed memory.  So instead we'll read into a local
    // buffer, broadcast that, and then copy that into the managed
    // memory.

    amrex::Vector<amrex::Real> f_local(static_cast<size_t>(9) * imax * jmax);
    amrex::Vector<amrex::Real> dpdf_local(static_cast<size_t>(4) * imax * jmax);
    amrex::Vector<amrex::Real> ef_local(static_cast<size_t>(4) * imax * jmax);
    amrex::Vector<amrex::Real> xf_local(static_cast<size_t>(4) * imax * jmax);

    if (amrex::ParallelDescriptor::IOProcessor()) {

        // open the table
        std::ifstream table;

#if defined(USE_NEW_HELMHOLTZ_TABLE)
        table.open("helm_table_fastmath_p128_q800.dat");
#else
        table.open("helm_table.dat");
#endif

        if (!table.is_open()) {
            // the table was not present or we could not open it; abort
            amrex::Error("helm_table.dat could not be opened");
        }

        std::string line;

        // read in the free energy table
        int idx = 0;
        for (int j = 0; j < jmax; ++j) {
            for (int i = 0; i < imax; ++i) {
                std::getline(table, line);
                if (line.empty()) {
                    amrex::Error("Error reading free energy from helm_table.dat");
                }
                std::istringstream data(line);
                data >> f_local[idx] >> f_local[idx+3] >> f_local[idx+1]
                     >> f_local[idx+4] >> f_local[idx+2] >> f_local[idx+5]
                     >> f_local[idx+6] >> f_local[idx+7] >> f_local[idx+8];
                idx += 9;
            }
        }

        // read the pressure derivative with density table
        idx = 0;
        for (int j = 0; j < jmax; ++j) {
            for (int i = 0; i < imax; ++i) {
                std::getline(table, line);
                if (line.empty()) {
                    amrex::Error("Error reading pressure derivative from helm_table.dat");
                }
                std::istringstream data(line);
                data >> dpdf_local[idx] >> dpdf_local[idx+2]
                     >> dpdf_local[idx+1] >> dpdf_local[idx+3];
                idx += 4;
            }
        }

        // read the electron chemical potential table
        idx = 0;
        for (int j = 0; j < jmax; ++j) {
            for (int i = 0; i < imax; ++i) {
                std::getline(table, line);
                if (line.empty()) {
                    amrex::Error("Error reading electron chemical potential from helm_table.dat");
                }
                std::istringstream data(line);
                data >> ef_local[idx] >> ef_local[idx+2]
                     >> ef_local[idx+1] >> ef_local[idx+3];
                idx += 4;
            }
        }

        // read the number density table
        idx = 0;
        for (int j = 0; j < jmax; ++j) {
            for (int i = 0; i < imax; ++i) {
                std::getline(table, line);
                if (line.empty()) {
                    amrex::Error("Error reading number density from helm_table.dat");
                }
                std::istringstream data(line);
                data >> xf_local[idx] >> xf_local[idx+2]
                     >> xf_local[idx+1] >> xf_local[idx+3];
                idx += 4;
            }
        }

        table.close();

    }

    amrex::ParallelDescriptor::Bcast(f_local.data(),    static_cast<size_t>(9) * imax * jmax);
    amrex::ParallelDescriptor::Bcast(dpdf_local.data(), static_cast<size_t>(4) * imax * jmax);
    amrex::ParallelDescriptor::Bcast(ef_local.data(),   static_cast<size_t>(4) * imax * jmax);
    amrex::ParallelDescriptor::Bcast(xf_local.data(),   static_cast<size_t>(4) * imax * jmax);

    // now copy into managed memory
    int idx = 0;
    for (int j = 0; j < jmax; ++j) {  // NOLINT(modernize-loop-convert)
        for (int i = 0; i < imax; ++i) {
            for (int m = 0; m < 9; ++m) {
                f[j][i][m] = f_local[idx];
                idx++;
            }
        }
    }

    idx = 0;
    for (int j = 0; j < jmax; ++j) {
        for (int i = 0; i < imax; ++i) {
            for (int m = 0; m < 4; ++m) {
                dpdf[j][i][m] = dpdf_local[idx];
                ef[j][i][m] = ef_local[idx];
                xf[j][i][m] = xf_local[idx];
                idx++;
            }
        }
    }

    // construct the temperature and density deltas and their inverses
    for (int j = 0; j < jmax-1; ++j)
    {
        dth         = t[j+1] - t[j];
        dt2         = dth * dth;
        dti         = 1.0e0_rt / dth;
        dt2i        = 1.0e0_rt / dt2;
        dt_sav[j]   = dth;
        dt2_sav[j]  = dt2;
        dti_sav[j]  = dti;
        dt2i_sav[j] = dt2i;
    }

    for (int i = 0; i < imax-1; ++i)
    {
        dd          = d[i+1] - d[i];
        dd2         = dd * dd;
        ddi         = 1.0e0_rt / dd;
        dd2i        = 1.0e0_rt / dd2;
        dd_sav[i]   = dd;
        dd2_sav[i]  = dd2;
        ddi_sav[i]  = ddi;
        dd2i_sav[i] = dd2i;
    }

    // Set up the minimum and maximum possible densities.

#if defined(USE_NEW_HELMHOLTZ_TABLE)
    EOSData::mintemp = mp::fastpow2(tlo);
    EOSData::maxtemp = mp::fastpow2(thi);
    EOSData::mindens = mp::fastpow2(dlo);
    EOSData::maxdens = mp::fastpow2(dhi);
#else
    EOSData::mintemp = std::pow(10.e0_rt, tlo);
    EOSData::maxtemp = std::pow(10.e0_rt, thi);
    EOSData::mindens = std::pow(10.e0_rt, dlo);
    EOSData::maxdens = std::pow(10.e0_rt, dhi);
#endif

}



AMREX_INLINE
void actual_eos_finalize ()
{
}



template <typename I>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool is_input_valid (I input)
{
  amrex::ignore_unused(input);

  static_assert(std::is_same_v<I, eos_input_t>, "input must be an eos_input_t");

  bool valid = true;

  return valid;
}

#endif
