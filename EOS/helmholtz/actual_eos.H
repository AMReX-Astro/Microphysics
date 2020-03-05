#ifndef _actual_eos_H_
#define _actual_eos_H_

#include <string>
#include <iostream>
#include <sstream>
#include <extern_parameters.H>
#include <fundamental_constants.H>
#include <actual_eos_data.H>
#include <eos_type.H>
#include <eos_data.H>
#include <cmath>

// Frank Timmes Helmholtz based Equation of State
// http://cococubed.asu.edu/
  
// given a temperature temp [K], density den [g/cm**3], and a composition
// characterized by abar and zbar, this routine returns most of the other
// thermodynamic quantities. of prime interest is the pressure [erg/cm**3],
// specific thermal energy [erg/gr], the entropy [erg/g/K], along with
// their derivatives with respect to temperature, density, abar, and zbar.
// other quantites such the normalized chemical potential eta (plus its
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

const std::string eos_name = "helmholtz";

// quintic hermite polynomial functions
// psi0 and its derivatives
AMREX_GPU_HOST_DEVICE inline
Real psi0 (Real z)
{
    return z * z * z * (z * (-6.0e0_rt * z + 15.0e0_rt) -10.0e0_rt) + 1.0e0_rt;
}

AMREX_GPU_HOST_DEVICE inline
Real dpsi0 (Real z)
{
    return z * z * (z * (-30.0e0_rt*z + 60.0e0_rt) - 30.0e0_rt);
}

AMREX_GPU_HOST_DEVICE inline
Real ddpsi0 (Real z)
{
    return z * (z * (-120.0e0_rt * z + 180.0e0_rt) - 60.0e0_rt);
}

// psi1 and its derivatives
AMREX_GPU_HOST_DEVICE inline
Real psi1 (Real z)
{
    return z * ( z * z * (z * (-3.0e0_rt * z + 8.0e0_rt) - 6.0e0_rt) + 1.0e0_rt);
}

AMREX_GPU_HOST_DEVICE inline
Real dpsi1 (Real z)
{
    return z * z * ( z * (-15.0e0_rt * z + 32.0e0_rt) - 18.0e0_rt) + 1.0e0_rt;
}

AMREX_GPU_HOST_DEVICE inline
Real ddpsi1 (Real z)
{
    return z * (z * (-60.0e0_rt * z + 96.0e0_rt) - 36.0e0_rt);
}

// psi2 and its derivatives
AMREX_GPU_HOST_DEVICE inline
Real psi2 (Real z)
{
    return 0.5e0_rt * z * z * ( z * ( z * (-z + 3.0e0_rt) - 3.0e0_rt) + 1.0e0_rt);
}

AMREX_GPU_HOST_DEVICE inline
Real dpsi2 (Real z)
{
    return 0.5e0_rt * z * (z * (z * (-5.0e0_rt * z + 12.0e0_rt) - 9.0e0_rt) + 2.0e0_rt);
}

AMREX_GPU_HOST_DEVICE inline
Real ddpsi2 (Real z)
{
    return 0.5e0_rt * (z * ( z * (-20.0e0_rt * z + 36.0e0_rt) - 18.0e0_rt) + 2.0e0_rt);
}

AMREX_GPU_HOST_DEVICE inline
void fwt (const Real* fi, const Real* wt, Real* fwtr)
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
AMREX_GPU_HOST_DEVICE inline
Real xpsi0(Real z)
{
    return z * z * (2.0e0_rt * z - 3.0e0_rt) + 1.0_rt;
}

AMREX_GPU_HOST_DEVICE inline
Real xdpsi0(Real z)
{
    return z * (6.0e0_rt * z - 6.0e0_rt);
}

// psi1 & derivatives
AMREX_GPU_HOST_DEVICE inline
Real xpsi1(Real z)
{
    return z * (z * (z - 2.0e0_rt) + 1.0e0_rt);
}

AMREX_GPU_HOST_DEVICE inline
Real xdpsi1(Real z)
{
    return z * (3.0e0_rt * z - 4.0e0_rt) + 1.0e0_rt;
}



AMREX_GPU_HOST_DEVICE inline
void apply_electrons(eos_t& state)
{

    // assume complete ionization
    Real ytot1 = 1.0e0_rt / state.abar;

    // define mu -- the total mean molecular weight including both electrons and ions
    state.mu = 1.0e0_rt / (1.0e0_rt / state.abar + 1.0e0_rt / state.mu_e);

    // enter the table with ye*den
    Real din = state.y_e * state.rho;

    // hash locate this temperature and density
    int jat = int((std::log10(state.T) - tlo) * tstpi) + 1;
    jat = amrex::max(1, amrex::min(jat, jtmax-1)) - 1;
    int iat = int((std::log10(din) - dlo) * dstpi) + 1;
    iat = amrex::max(1, amrex::min(iat, itmax-1)) - 1;

    Real fi[36];

    // access the table locations only once
    for (int i = 0; i < 9; ++i) {
        fi[i     ] = _f[jat  ][iat  ][i]; // f, ft, ftt, fd, fdd, fdt, fddt, fdtt, fddtt
        fi[i +  9] = _f[jat  ][iat+1][i];
        fi[i + 18] = _f[jat+1][iat  ][i];
        fi[i + 27] = _f[jat+1][iat+1][i];
    }

    // various differences
    Real xt  = amrex::max((state.T - _t[jat]) * dti_sav[jat], 0.0e0_rt);
    Real xd  = amrex::max((din - _d[iat]) * ddi_sav[iat], 0.0e0_rt);
    Real mxt = 1.0e0_rt - xt;
    Real mxd = 1.0e0_rt - xd;

    // the six density and six temperature basis functions
    Real sit[6];

    sit[0] = psi0(xt);
    sit[1] = psi1(xt) * dt_sav[jat];
    sit[2] = psi2(xt) * dt2_sav[jat];

    sit[3] =  psi0(mxt);
    sit[4] = -psi1(mxt) * dt_sav[jat];
    sit[5] =  psi2(mxt) * dt2_sav[jat];

    Real sid[6];

    sid[0] =  psi0(xd);
    sid[1] =  psi1(xd) * dd_sav[iat];
    sid[2] =  psi2(xd) * dd2_sav[iat];

    sid[3] =  psi0(mxd);
    sid[4] = -psi1(mxd) * dd_sav[iat];
    sid[5] =  psi2(mxd) * dd2_sav[iat];

    // derivatives of the weight functions
    Real dsit[6];

    dsit[0] =  dpsi0(xt) * dti_sav[jat];
    dsit[1] =  dpsi1(xt);
    dsit[2] =  dpsi2(xt) * dt_sav[jat];

    dsit[3] = -dpsi0(mxt) * dti_sav[jat];
    dsit[4] =  dpsi1(mxt);
    dsit[5] = -dpsi2(mxt) * dt_sav[jat];

    Real dsid[6];

    dsid[0] =  dpsi0(xd) * ddi_sav[iat];
    dsid[1] =  dpsi1(xd);
    dsid[2] =  dpsi2(xd) * dd_sav[iat];

    dsid[3] = -dpsi0(mxd) * ddi_sav[iat];
    dsid[4] =  dpsi1(mxd);
    dsid[5] = -dpsi2(mxd) * dd_sav[iat];

    // second derivatives of the weight functions
    Real ddsit[6];

    ddsit[0] =  ddpsi0(xt) * dt2i_sav[jat];
    ddsit[1] =  ddpsi1(xt) * dti_sav[jat];
    ddsit[2] =  ddpsi2(xt);

    ddsit[3] =  ddpsi0(mxt) * dt2i_sav[jat];
    ddsit[4] = -ddpsi1(mxt) * dti_sav[jat];
    ddsit[5] =  ddpsi2(mxt);

    // This array saves some subexpressions that go into
    // computing the biquintic polynomial. Instead of explicitly
    // constructing it in full, we'll use these subexpressions
    // and then compute the result as
    // (table data * temperature terms) * density terms.

    Real fwtr[6];

    fwt(fi, sit, fwtr);

    Real free = 0.e0_rt;
    Real df_d = 0.e0_rt;

    for (int i = 0; i <= 5; ++i) {
        // the free energy
        free = free + fwtr[i] * sid[i];

        // derivative with respect to density
        df_d = df_d + fwtr[i] * dsid[i];
    }

    fwt(fi, dsit, fwtr);

    Real df_t = 0.e0_rt;
    Real df_dt = 0.e0_rt;

    for (int i = 0; i <= 5; ++i) {
        // derivative with respect to temperature
        df_t += fwtr[i] * sid[i];

        // derivative with respect to temperature and density
        df_dt += fwtr[i] * dsid[i];
    }

    fwt(fi, ddsit, fwtr);

    Real df_tt = 0.e0_rt;
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

    // derivatives of weight functions
    dsit[0] = xdpsi0(xt) * dti_sav[jat];
    dsit[1] = xdpsi1(xt);

    dsit[2] = -xdpsi0(mxt) * dti_sav[jat];
    dsit[3] = xdpsi1(mxt);

    dsid[0] = xdpsi0(xd) * ddi_sav[iat];
    dsid[1] = xdpsi1(xd);

    dsid[2] = -xdpsi0(mxd) * ddi_sav[iat];
    dsid[3] = xdpsi1(mxd);

    // Reuse subexpressions that would go into computing the
    // cubic interpolation.
    Real wdt[16];

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
    // the interpolation correctly. Alternate indexing schemes are
    // possible if we were to reorder wdt.
    fi[ 0] = _dpdf[jat  ][iat  ][0];
    fi[ 1] = _dpdf[jat  ][iat  ][1];
    fi[ 4] = _dpdf[jat  ][iat  ][2];
    fi[ 5] = _dpdf[jat  ][iat  ][3];

    fi[ 8] = _dpdf[jat  ][iat+1][0];
    fi[ 9] = _dpdf[jat  ][iat+1][1];
    fi[12] = _dpdf[jat  ][iat+1][2];
    fi[13] = _dpdf[jat  ][iat+1][3];

    fi[ 2] = _dpdf[jat+1][iat  ][0];
    fi[ 3] = _dpdf[jat+1][iat  ][1];
    fi[ 6] = _dpdf[jat+1][iat  ][2];
    fi[ 7] = _dpdf[jat+1][iat  ][3];

    fi[10] = _dpdf[jat+1][iat+1][0];
    fi[11] = _dpdf[jat+1][iat+1][1];
    fi[14] = _dpdf[jat+1][iat+1][2];
    fi[15] = _dpdf[jat+1][iat+1][3];

    // pressure derivative with density
    Real dpepdd = 0.0e0_rt;
    for (int i = 0; i <= 15; ++i) {
        dpepdd = dpepdd + fi[i] * wdt[i];
    }
    dpepdd = amrex::max(state.y_e * dpepdd, 0.0e0_rt);

    // Read in the tabular data for the electron chemical potential.
    fi[ 0] = _ef[jat  ][iat  ][0];
    fi[ 1] = _ef[jat  ][iat  ][1];
    fi[ 4] = _ef[jat  ][iat  ][2];
    fi[ 5] = _ef[jat  ][iat  ][3];

    fi[ 8] = _ef[jat  ][iat+1][0];
    fi[ 9] = _ef[jat  ][iat+1][1];
    fi[12] = _ef[jat  ][iat+1][2];
    fi[13] = _ef[jat  ][iat+1][3];

    fi[ 2] = _ef[jat+1][iat  ][0];
    fi[ 3] = _ef[jat+1][iat  ][1];
    fi[ 6] = _ef[jat+1][iat  ][2];
    fi[ 7] = _ef[jat+1][iat  ][3];

    fi[10] = _ef[jat+1][iat+1][0];
    fi[11] = _ef[jat+1][iat+1][1];
    fi[14] = _ef[jat+1][iat+1][2];
    fi[15] = _ef[jat+1][iat+1][3];

    // electron chemical potential etaele
    Real etaele = 0.0e0_rt;
    for (int i = 0; i <= 15; ++i) {
        etaele = etaele + fi[i] * wdt[i];
    }

    // Read in the tabular data for the number density.
    fi[ 0] = _xf[jat  ][iat  ][0];
    fi[ 1] = _xf[jat  ][iat  ][1];
    fi[ 4] = _xf[jat  ][iat  ][2];
    fi[ 5] = _xf[jat  ][iat  ][3];

    fi[ 8] = _xf[jat  ][iat+1][0];
    fi[ 9] = _xf[jat  ][iat+1][1];
    fi[12] = _xf[jat  ][iat+1][2];
    fi[13] = _xf[jat  ][iat+1][3];

    fi[ 2] = _xf[jat+1][iat  ][0];
    fi[ 3] = _xf[jat+1][iat  ][1];
    fi[ 6] = _xf[jat+1][iat  ][2];
    fi[ 7] = _xf[jat+1][iat  ][3];

    fi[10] = _xf[jat+1][iat+1][0];
    fi[11] = _xf[jat+1][iat+1][1];
    fi[14] = _xf[jat+1][iat+1][2];
    fi[15] = _xf[jat+1][iat+1][3];

    // electron + positron number densities
    Real xnefer = 0.0e0_rt;
    for (int i = 0; i <= 15; ++i) {
        xnefer = xnefer + fi[i] * wdt[i];
    }

    // the desired electron-positron thermodynamic quantities

    // dpepdd at high temperatures and low densities is below the
    // floating point limit of the subtraction of two large terms.
    // since state.dpdr doesn't enter the maxwell relations at all, use the
    // bicubic interpolation done above instead of this one
    Real x       = din * din;
    Real pele    = x * df_d;
    Real dpepdt  = x * df_dt;
    Real s       = dpepdd/state.y_e - 2.0e0_rt * din * df_d;
#ifdef EXTRA_THERMO
    Real dpepda  = -ytot1 * (2.0e0_rt * pele + s * din);
    Real dpepdz  = state.rho*ytot1*(2.0e0_rt * din * df_d  +  s);
#endif

    x            = state.y_e * state.y_e;
    Real sele    = -df_t * state.y_e;
    Real dsepdt  = -df_tt * state.y_e;
    Real dsepdd  = -df_dt * x;
#ifdef EXTRA_THERMO
    Real dsepda  = ytot1 * (state.y_e * df_dt * din - sele);
    Real dsepdz  = -ytot1 * (state.y_e * df_dt * state.rho  + df_t);
#endif

    Real eele    = state.y_e*free + state.T * sele;
    Real deepdt  = state.T * dsepdt;
    Real deepdd  = x * df_d + state.T * dsepdd;
#ifdef EXTRA_THERMO
    Real deepda  = -state.y_e * ytot1 * (free +  df_d * din) + state.T * dsepda;
    Real deepdz  = ytot1* (free + state.y_e * df_d * state.rho) + state.T * dsepdz;
#endif

    state.p    = state.p + pele;
    state.dpdT = state.dpdT + dpepdt;
    state.dpdr = state.dpdr + dpepdd;
#ifdef EXTRA_THERMO
    state.dpdA = state.dpdA + dpepda;
    state.dpdZ = state.dpdZ + dpepdz;
#endif

    state.s    = state.s + sele;
    state.dsdT = state.dsdT + dsepdt;
    state.dsdr = state.dsdr + dsepdd;

    state.e    = state.e + eele;
    state.dedT = state.dedT + deepdt;
    state.dedr = state.dedr + deepdd;
#ifdef EXTRA_THERMO
    state.dedA = state.dedA + deepda;
    state.dedZ = state.dedZ + deepdz;
#endif

    state.eta = etaele;
    state.xne = xnefer;
    state.xnp = 0.0e0_rt;

    state.pele = pele;
    state.ppos = 0.0e0_rt;

}



AMREX_GPU_HOST_DEVICE inline
void apply_ions(eos_t& state)
{

    const Real pi      = 3.1415926535897932384e0_rt;
    const Real sioncon = (2.0e0_rt * pi * amu * kerg)/(h*h);
    const Real kergavo = kerg * avo_eos;

    Real deni = 1.0e0_rt / state.rho;
    Real tempi = 1.0e0_rt / state.T;

    Real ytot1   = 1.0e0_rt / state.abar;
    Real xni     = avo_eos * ytot1 * state.rho;
    Real dxnidd  = avo_eos * ytot1;
    Real dxnida  = -xni * ytot1;

    Real kt = kerg * state.T;

    Real pion    = xni * kt;
    Real dpiondd = dxnidd * kt;
    Real dpiondt = xni * kerg;
#ifdef EXTRA_THERMO
    Real dpionda = dxnida * kt;
    Real dpiondz = 0.0e0_rt;
#endif

    Real eion    = 1.5e0_rt * pion * deni;
    Real deiondd = (1.5e0_rt * dpiondd - eion) * deni;
    Real deiondt = 1.5e0_rt * dpiondt * deni;
#ifdef EXTRA_THERMO
    Real deionda = 1.5e0_rt * dpionda * deni;
    Real deiondz = 0.0e0_rt;
#endif

    Real x       = state.abar * state.abar * std::sqrt(state.abar) * deni / avo_eos;
    Real s       = sioncon * state.T;
    Real z       = x * s * std::sqrt(s);
    Real y       = std::log(z);
    Real sion    = (pion * deni + eion) * tempi + kergavo * ytot1 * y;
    Real dsiondd = (dpiondd * deni - pion * deni * deni + deiondd) * tempi -
                   kergavo * deni * ytot1;
    Real dsiondt = (dpiondt * deni + deiondt) * tempi -  
                   (pion * deni + eion) * tempi * tempi +
                   1.5e0_rt * kergavo * tempi * ytot1;

    state.p    = state.p + pion;
    state.dpdT = state.dpdT + dpiondt;
    state.dpdr = state.dpdr + dpiondd;
#ifdef EXTRA_THERMO
    state.dpdA = state.dpdA + dpionda;
    state.dpdZ = state.dpdZ + dpiondz;
#endif

    state.e    = state.e + eion;
    state.dedT = state.dedT + deiondt;
    state.dedr = state.dedr + deiondd;
#ifdef EXTRA_THERMO
    state.dedA = state.dedA + deionda;
    state.dedZ = state.dedZ + deiondz;
#endif

    state.s    = state.s + sion;
    state.dsdT = state.dsdT + dsiondt;
    state.dsdr = state.dsdr + dsiondd;

}



AMREX_GPU_HOST_DEVICE inline
void apply_radiation(eos_t& state)
{

    const Real clight  = 2.99792458e10_rt;
#ifdef RADIATION
    const Real ssol    = 0.0e0_rt;
#else
    const Real ssol    = 5.67051e-5_rt;
#endif
    const Real asol    = 4.0e0_rt * ssol / clight;
    const Real asoli3  = asol/3.0e0_rt;

    Real deni = 1.0e0_rt / state.rho;
    Real tempi = 1.0e0_rt / state.T;

    Real prad = asoli3 * state.T * state.T * state.T * state.T;

    // In low density material, this radiation pressure becomes unphysically high.
    // For rho ~ 1 g/cc and T ~ 1e9, then this radiation pressure will lead to a
    // sound speed ~ 0.1c, and this is a problem because the codes using this EOS
    // are primarily non-relativistic. Thus we can optionally smooth out the radiation
    // pressure such that it disappears at low densities. Since all terms below depend
    // on prad, this will result in the radiation term effectively vanishing below the
    // cutoff density. For simplicity we ignore the effect this has on the derivatives.

    if (prad_limiter_rho_c > 0.0e0_rt && prad_limiter_delta_rho > 0.0e0_rt) {
        prad = prad * 0.5e0_rt * (1.0e0_rt + std::tanh((state.rho - prad_limiter_rho_c) / prad_limiter_delta_rho));
    }

    Real dpraddd = 0.0e0_rt;
    Real dpraddt = 4.0e0_rt * prad * tempi;
#ifdef EXTRA_THERMO
    Real dpradda = 0.0e0_rt;
    Real dpraddz = 0.0e0_rt;
#endif

    Real erad    = 3.0e0_rt * prad*deni;
    Real deraddd = -erad * deni;
    Real deraddt = 3.0e0_rt * dpraddt * deni;
#ifdef EXTRA_THERMO
    Real deradda = 0.0e0_rt;
    Real deraddz = 0.0e0_rt;
#endif

    Real srad    = (prad * deni + erad) * tempi;
    Real dsraddd = (dpraddd * deni - prad * deni * deni + deraddd) * tempi;
    Real dsraddt = (dpraddt * deni + deraddt - srad) * tempi;

    // Note that unlike the other terms, radiation
    // sets these terms instead of adding to them,
    // since it comes first.

    state.p    = prad;
    state.dpdr = dpraddd;
    state.dpdT = dpraddt;
#ifdef EXTRA_THERMO
    state.dpdA = dpradda;
    state.dpdZ = dpraddz;
#endif

    state.e    = erad;
    state.dedr = deraddd;
    state.dedT = deraddt;
#ifdef EXTRA_THERMO
    state.dedA = deradda;
    state.dedZ = deraddz;
#endif

    state.s    = srad;
    state.dsdr = dsraddd;
    state.dsdT = dsraddt;

}



AMREX_GPU_HOST_DEVICE inline
void apply_coulomb_corrections(eos_t& state)
{

    // Constants used for the Coulomb corrections
    const Real a1 = -0.898004e0_rt;
    const Real b1 =  0.96786e0_rt;
    const Real c1 =  0.220703e0_rt;
    const Real d1 = -0.86097e0_rt;
    const Real e1 =  2.5269e0_rt;
    const Real a2 =  0.29561e0_rt;
    const Real b2 =  1.9885e0_rt;
    const Real c2 =  0.288675e0_rt;
    const Real qe   = 4.8032042712e-10_rt;
    const Real esqu = qe * qe;
    const Real onethird = 1.0e0_rt/3.0e0_rt;
    const Real forth = 4.0e0_rt/3.0e0_rt;
    const Real pi    = 3.1415926535897932384e0_rt;

    Real pcoul    = 0.e0_rt;
    Real dpcouldd = 0.e0_rt;
    Real dpcouldt = 0.e0_rt;
    Real ecoul    = 0.e0_rt;
    Real decouldd = 0.e0_rt;
    Real decouldt = 0.e0_rt;
    Real scoul    = 0.e0_rt;
    Real dscouldd = 0.e0_rt;
    Real dscouldt = 0.e0_rt;

#ifdef EXTRA_THERMO
    Real dpcoulda = 0.e0_rt;
    Real dpcouldz = 0.e0_rt;
    Real decoulda = 0.e0_rt;
    Real decouldz = 0.e0_rt;
#endif

    Real s, x, y, z; // temporary variables

    // uniform background corrections only
    // from yakovlev & shalybkov 1989
    // lami is the average ion seperation
    // plasg is the plasma coupling parameter

    Real ytot1 = 1.0e0_rt / state.abar;
    Real xni     = avo_eos * ytot1 * state.rho;
    Real dxnidd  = avo_eos * ytot1;
    Real dxnida  = -xni * ytot1;

    Real kt      = kerg * state.T;
    Real ktinv   = 1.0e0_rt / kt;

    z             = forth * pi;
    s             = z * xni;
    Real dsdd     = z * dxnidd;
    Real dsda     = z * dxnida;

    Real lami     = 1.0e0_rt / std::pow(s, onethird);
    Real inv_lami = 1.0e0_rt / lami;
    z             = -onethird * lami;
    Real lamidd   = z * dsdd / s;
    Real lamida   = z * dsda / s;

    Real plasg    = state.zbar * state.zbar * esqu * ktinv * inv_lami;
    z             = -plasg * inv_lami;
    Real plasgdd  = z * lamidd;
    Real plasgda  = z * lamida;
    Real plasgdt  = -plasg*ktinv * kerg;
    Real plasgdz  = 2.0e0_rt * plasg/state.zbar;

    // .yakovlev & shalybkov 1989 equations 82, 85, 86, 87
    if (plasg >= 1.0e0_rt)
    {
        x        = std::pow(plasg, 0.25e0_rt);
        y        = avo_eos * ytot1 * kerg;
        ecoul    = y * state.T * (a1 * plasg + b1 * x + c1 / x + d1);
        pcoul    = onethird * state.rho * ecoul;
        scoul    = -y * (3.0e0_rt * b1 * x - 5.0e0_rt*c1 / x +
                    d1 * (std::log(plasg) - 1.0e0_rt) - e1);

        y        = avo_eos*ytot1*kt*(a1 + 0.25e0_rt/plasg*(b1*x - c1/x));
        decouldd = y * plasgdd;
        decouldt = y * plasgdt + ecoul/state.T;

#ifdef EXTRA_THERMO
        decoulda = y * plasgda - ecoul/state.abar;
        decouldz = y * plasgdz;
#endif

        y        = onethird * state.rho;
        dpcouldd = onethird * ecoul + y * decouldd;
        dpcouldt = y * decouldt;

#ifdef EXTRA_THERMO
        dpcoulda = y * decoulda;
        dpcouldz = y * decouldz;
#endif

       y        = -avo_eos * kerg / (state.abar * plasg) *
                   (0.75e0_rt * b1 * x + 1.25e0_rt * c1 / x + d1);
       dscouldd = y * plasgdd;
       dscouldt = y * plasgdt;

       // yakovlev & shalybkov 1989 equations 102, 103, 104
    }
    else if (plasg < 1.0e0_rt)
    {

        Real pion    = xni * kt;
        Real dpiondd = dxnidd * kt;
        Real dpiondt = xni * kerg;
#ifdef EXTRA_THERMO
        Real dpionda = dxnida * kt;
        Real dpiondz = 0.0e0_rt;
#endif

        x        = plasg * std::sqrt(plasg);
        y        = std::pow(plasg, b2);
        z        = c2 * x - onethird * a2 * y;
        pcoul    = -pion * z;
        ecoul    = 3.0e0_rt * pcoul / state.rho;
        scoul    = -avo_eos / state.abar * kerg * (c2 * x -a2 * (b2 - 1.0e0_rt) / b2 * y);

        s        = 1.5e0_rt * c2 * x / plasg - onethird * a2 * b2 * y / plasg;
        dpcouldd = -dpiondd * z - pion * s * plasgdd;
        dpcouldt = -dpiondt * z - pion * s * plasgdt;
#ifdef EXTRA_THERMO
        dpcoulda = -dpionda * z - pion * s * plasgda;
        dpcouldz = -dpiondz * z - pion * s * plasgdz;
#endif

        s        = 3.0e0_rt / state.rho;
        decouldd = s * dpcouldd - ecoul / state.rho;
        decouldt = s * dpcouldt;
#ifdef EXTRA_THERMO
        decoulda = s * dpcoulda;
        decouldz = s * dpcouldz;
#endif

        s        = -avo_eos * kerg / (state.abar * plasg) * 
                    (1.5e0_rt * c2 * x - a2 * (b2 - 1.0e0_rt) * y);
        dscouldd = s * plasgdd;
        dscouldt = s * plasgdt;
    }

    // Disable Coulomb corrections if they cause
    // the energy or pressure to go negative.

    Real p_temp = state.p + pcoul;
    Real e_temp = state.e + ecoul;

    if (p_temp <= 0.0e0_rt || e_temp <= 0.0e0_rt)
    {

        pcoul    = 0.0e0_rt;
        dpcouldd = 0.0e0_rt;
        dpcouldt = 0.0e0_rt;
        ecoul    = 0.0e0_rt;
        decouldd = 0.0e0_rt;
        decouldt = 0.0e0_rt;
        scoul    = 0.0e0_rt;
        dscouldd = 0.0e0_rt;
        dscouldt = 0.0e0_rt;

#ifdef EXTRA_THERMO
        dpcoulda = 0.0e0_rt;
        dpcouldz = 0.0e0_rt;
        decoulda = 0.0e0_rt;
        decouldz = 0.0e0_rt;
#endif

    }

    state.p    = state.p + pcoul;
    state.dpdr = state.dpdr + dpcouldd;
    state.dpdT = state.dpdT + dpcouldt;
#ifdef EXTRA_THERMO
    state.dpdA = state.dpdA + dpcoulda;
    state.dpdZ = state.dpdZ + dpcouldz;
#endif

    state.e    = state.e + ecoul;
    state.dedr = state.dedr + decouldd;
    state.dedT = state.dedT + decouldt;
#ifdef EXTRA_THERMO
    state.dedA = state.dedA + decoulda;
    state.dedZ = state.dedZ + decouldz;
#endif

    state.s    = state.s + scoul;
    state.dsdr = state.dsdr + dscouldd;
    state.dsdT = state.dsdT + dscouldt;

}



AMREX_GPU_HOST_DEVICE inline
void prepare_for_iterations(eos_input_t input, eos_t& state,
                            bool& single_iter, Real& v_want, Real& v1_want, Real& v2_want,
                            int& var, int& dvar, int& var1, int& var2)
{

    single_iter = true;

    if (input == eos_input_rt) {

        // Nothing to do here.

    }
    else if (input == eos_input_rh) {

        v_want = state.h;
        var  = ienth;
        dvar = itemp;

    }
    else if (input == eos_input_tp) {

        v_want = state.p;
        var  = ipres;
        dvar = idens;

    }
    else if (input == eos_input_rp) {

        v_want = state.p;
        var  = ipres;
        dvar = itemp;

    }
    else if (input == eos_input_re) {

        v_want = state.e;
        var  = iener;
        dvar = itemp;

    }
    else if (input == eos_input_ps) {

        single_iter = false;
        v1_want = state.p;
        v2_want = state.s;
        var1 = ipres;
        var2 = ientr;

    }
    else if (input == eos_input_ph) {

        single_iter = false;
        v1_want = state.p;
        v2_want = state.h;
        var1 = ipres;
        var2 = ienth;

    }
    else if (input == eos_input_th) {

        v_want = state.h;
        var  = ienth;
        dvar = idens;

    }

}




AMREX_GPU_HOST_DEVICE inline
void single_iter_update(eos_t& state, int var, int dvar,
                        Real v_want, bool& converged)
{

    Real x, v, dvdx, xtol, smallx;

    if (dvar == itemp) {

        x = state.T;
        smallx = EOSData::mintemp;
        xtol = ttol;

        if (var == ipres) {
            v    = state.p;
            dvdx = state.dpdT;
        }
        else if (var == iener) {
            v    = state.e;
            dvdx = state.dedT;
        }
        else if (var == ientr) {
            v    = state.s;
            dvdx = state.dsdT;
        }
        else if (var == ienth) {
            v    = state.h;
            dvdx = state.dhdT;
        }

    }
    else // dvar == density
    {

        x = state.rho;
        smallx = EOSData::mindens;
        xtol = dtol;

        if (var == ipres) {
            v    = state.p;
            dvdx = state.dpdr;
        }
        else if (var == iener) {
            v    = state.e;
            dvdx = state.dedr;
        }
        else if (var == ientr) {
            v    = state.s;
            dvdx = state.dsdr;
        }
        else if (var == ienth) {
            v    = state.h;
            dvdx = state.dhdr;
        }

    }

    // Now do the calculation for the next guess for T/rho
    Real xnew = x - (v - v_want) / dvdx;

    // Don't let the temperature/density change by more than a factor of two
    xnew = amrex::max(0.5_rt * x, amrex::min(xnew, 2.0_rt * x));

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

    Real error = std::abs((xnew - x) / x);

    if (error < xtol) converged = true;

}



AMREX_GPU_HOST_DEVICE inline
void double_iter_update(eos_t& state, int var1, int var2,
                        Real v1_want, Real v2_want, bool& converged)
{

    Real v1, dv1dt, dv1dr, v2, dv2dt, dv2dr;

    // Figure out which variables we're using

    Real told = state.T;
    Real rold = state.rho;

    if (var1 == ipres) {
        v1    = state.p;
        dv1dt = state.dpdT;
        dv1dr = state.dpdr;
    }
    else if (var1 == iener) {
        v1    = state.e;
        dv1dt = state.dedT;
        dv1dr = state.dedr;
    }
    else if (var1 == ientr) {
        v1    = state.s;
        dv1dt = state.dsdT;
        dv1dr = state.dsdr;
    }
    else if (var1 == ienth) {
        v1    = state.h;
        dv1dt = state.dhdT;
        dv1dr = state.dhdr;
    }

    if (var2 == ipres) {
        v2    = state.p;
        dv2dt = state.dpdT;
        dv2dr = state.dpdr;
    }
    else if (var2 == iener) {
        v2    = state.e;
        dv2dt = state.dedT;
        dv2dr = state.dedr;
    }
    else if (var2 == ientr) {
        v2    = state.s;
        dv2dt = state.dsdT;
        dv2dr = state.dsdr;
    }
    else if (var2 == ienth) {
        v2    = state.h;
        dv2dt = state.dhdT;
        dv2dr = state.dhdr;
    }

    // Two functions, f and g, to iterate over
    Real v1i = v1_want - v1;
    Real v2i = v2_want - v2;

    //
    // 0 = f + dfdr * delr + dfdt * delt
    // 0 = g + dgdr * delr + dgdt * delt
    //

    // note that dfi/dT = - df/dT
    Real delr = (-v1i * dv2dt + v2i * dv1dt) / (dv2dr * dv1dt - dv2dt * dv1dr);

    Real rnew = rold + delr;

    Real tnew = told + (v1i - dv1dr * delr) / dv1dt;

    // Don't let the temperature or density change by more
    // than a factor of two
    tnew = amrex::max(0.5e0_rt * told, amrex::min(tnew, 2.0e0_rt * told));
    rnew = amrex::max(0.5e0_rt * rold, amrex::min(rnew, 2.0e0_rt * rold));

    // Don't let us freeze or evacuate
    tnew = amrex::max(EOSData::mintemp, tnew);
    rnew = amrex::max(EOSData::mindens, rnew);

    // Store the new temperature and density
    state.rho = rnew;
    state.T = tnew;

    // Compute the errors
    Real error1 = std::abs((rnew - rold) / rold);
    Real error2 = std::abs((tnew - told) / told);

    if (error1 < dtol && error2 < ttol) converged = true;

}



AMREX_GPU_HOST_DEVICE inline
void finalize_state (eos_input_t input, eos_t& state,
                     Real v_want, Real v1_want, Real v2_want)
{

    // Calculate some remaining derivatives
    state.dpde = state.dpdT / state.dedT;
    state.dpdr_e = state.dpdr - state.dpdT * state.dedr / state.dedT;

    // Specific heats and Gamma_1
    Real chit = state.T / state.p * state.dpdT;
    Real chid = state.dpdr * state.rho / state.p;

    state.cv = state.dedT;
    state.gam1 = (chit * (state.p / state.rho)) * (chit / (state.T * state.cv)) + chid;
    state.cp = state.cv * state.gam1 / chid;

    // Use the non-relativistic version of the sound speed, cs = sqrt(gam_1 * P / rho).
    // This replaces the relativistic version that comes out of helmeos.
    state.cs = std::sqrt(state.gam1 * state.p / state.rho);

    if (input_is_constant) {

       if (input == eos_input_rh) {

           state.h = v_want;

       }
       else if (input == eos_input_tp) {

           state.p = v_want;

       }
       else if (input == eos_input_rp) {

           state.p = v_want;

       }
       else if (input == eos_input_re) {

           state.e = v_want;

       }
       else if (input == eos_input_ps) {

           state.p = v1_want;
           state.s = v2_want;

       }
       else if (input == eos_input_ph) {

           state.p = v1_want;
           state.h = v2_want;

       }
       else if (input == eos_input_th) {

           state.h = v_want;

       }

    }

}



AMREX_GPU_HOST_DEVICE inline
void actual_eos(eos_input_t input, eos_t& state)
{

    const int max_newton = 100;

    bool single_iter, converged;
    int var, dvar, var1, var2;
    Real v_want, v1_want, v2_want;

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

        state.h = state.e + state.p / state.rho;
        state.dhdr = state.dedr + state.dpdr / state.rho - state.p / (state.rho * state.rho);
        state.dhdT = state.dedT + state.dpdT / state.rho;

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



inline
void actual_eos_init ()
{

    Real dth, dt2, dti, dt2i;
    Real dd, dd2, ddi, dd2i;
    Real tsav, dsav;

    // Read in the runtime parameters

    input_is_constant = eos_input_is_constant;
    do_coulomb = use_eos_coulomb;
    ttol = eos_ttol;
    dtol = eos_dtol;

    //    read the helmholtz free energy table
    itmax = imax;
    jtmax = jmax;
    tlo   = 3.0e0_rt;
    thi   = 13.0e0_rt;
    tstp  = (thi - tlo) / ((Real) (jmax-1));
    tstpi = 1.0e0_rt / tstp;
    dlo   = -12.0e0_rt;
    dhi   = 15.0e0_rt;
    dstp  = (dhi - dlo) / ((Real) (imax-1));
    dstpi = 1.0e0_rt / dstp;

    for (int j = 0; j < jmax; ++j) {
        tsav = tlo + j * tstp;
        _t[j] = std::pow(10.0e0_rt, tsav);
        for (int i = 0; i < imax; ++i) {
            dsav = dlo + i * dstp;
            _d[i] = std::pow(10.0e0_rt, dsav);
        }
    }

    if (amrex::ParallelDescriptor::IOProcessor()) {

        // open the table
        std::ifstream table;
        table.open("helm_table.dat");

        std::string line;

        // read in the free energy table
        for (int j = 0; j < jmax; ++j) {
            for (int i = 0; i < imax; ++i) {
                std::getline(table, line);
                std::istringstream data(line);
                data >> _f[j][i][0] >> _f[j][i][3] >> _f[j][i][1] >> _f[j][i][4] >> _f[j][i][2]
                     >> _f[j][i][5] >> _f[j][i][6] >> _f[j][i][7] >> _f[j][i][8];
            }
        }

        // read the pressure derivative with density table
        for (int j = 0; j < jmax; ++j) {
            for (int i = 0; i < imax; ++i) {
                std::getline(table, line);
                std::istringstream data(line);
                data >> _dpdf[j][i][0] >> _dpdf[j][i][2] >> _dpdf[j][i][1] >> _dpdf[j][i][3];
            }
        }

        // read the electron chemical potential table
        for (int j = 0; j < jmax; ++j) {
            for (int i = 0; i < imax; ++i) {
                std::getline(table, line);
                std::istringstream data(line);
                data >> _ef[j][i][0] >> _ef[j][i][2] >> _ef[j][i][1] >> _ef[j][i][3];
            }
        }

        // read the number density table
        for (int j = 0; j < jmax; ++j) {
            for (int i = 0; i < imax; ++i) {
                std::getline(table, line);
                std::istringstream data(line);
                data >> _xf[j][i][0] >> _xf[j][i][2] >> _xf[j][i][1] >> _xf[j][i][3];
            }
        }

        table.close();

    }

    amrex::ParallelDescriptor::Bcast(&_f[0][0][0],    9 * imax * jmax);
    amrex::ParallelDescriptor::Bcast(&_dpdf[0][0][0], 4 * imax * jmax);
    amrex::ParallelDescriptor::Bcast(&_ef[0][0][0],   4 * imax * jmax);
    amrex::ParallelDescriptor::Bcast(&_xf[0][0][0],   4 * imax * jmax);

    // construct the temperature and density deltas and their inverses
    for (int j = 0; j < jmax-1; ++j)
    {
        dth         = _t[j+1] - _t[j];
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
        dd          = _d[i+1] - _d[i];
        dd2         = dd * dd;
        ddi         = 1.0e0_rt / dd;
        dd2i        = 1.0e0_rt / dd2;
        dd_sav[i]   = dd;
        dd2_sav[i]  = dd2;
        ddi_sav[i]  = ddi;
        dd2i_sav[i] = dd2i;
    }

    // Set up the minimum and maximum possible densities.

    EOSData::mintemp = std::pow(10.e0_rt, tlo);
    EOSData::maxtemp = std::pow(10.e0_rt, thi);
    EOSData::mindens = std::pow(10.e0_rt, dlo);
    EOSData::maxdens = std::pow(10.e0_rt, dhi);

}



inline
void actual_eos_finalize ()
{
}

#endif

