#ifndef SNEUT5_H
#define SNEUT5_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Math.H>

#include <microphysics_autodiff.H>
#include <extern_parameters.H>

using namespace amrex::literals;

namespace nu_constants {

    // numerical constants
    constexpr amrex::Real fac1   = 5.0e0_rt * M_PI / 3.0e0_rt;
    constexpr amrex::Real fac2   = 10.0e0_rt * M_PI;
    constexpr amrex::Real fac3   = M_PI / 5.0e0_rt;
    constexpr amrex::Real oneth  = 1.0e0_rt/3.0e0_rt;
    constexpr amrex::Real twoth  = 2.0e0_rt/3.0e0_rt;
    constexpr amrex::Real sixth  = 1.0e0_rt/6.0e0_rt;

    // theta is sin**2(theta_weinberg) = 0.2319 plus/minus 0.00005 (1996)
    // xnufam is the number of neutrino flavors = 3.02 plus/minus 0.005 (1998)
    // change theta and xnufam if need be, and the changes will automatically
    // propagate through the routine. cv and ca are the vector and axial currents.

    constexpr amrex::Real theta  = 0.2319e0_rt;
    constexpr amrex::Real xnufam = 3.0e0_rt;
    constexpr amrex::Real cv     = 0.5e0_rt + 2.0e0_rt * theta;
    constexpr amrex::Real cvp    = 1.0e0_rt - cv;
    constexpr amrex::Real ca     = 0.5e0_rt;
    constexpr amrex::Real cap    = 1.0e0_rt - ca;
    constexpr amrex::Real tfac1  = cv*cv + ca*ca + (xnufam-1.0e0_rt) * (cvp*cvp+cap*cap);
    constexpr amrex::Real tfac2  = cv*cv - ca*ca + (xnufam-1.0e0_rt) * (cvp*cvp - cap*cap);
    constexpr amrex::Real tfac3  = tfac2/tfac1;
    constexpr amrex::Real tfac4  = 0.5e0_rt * tfac1;
    constexpr amrex::Real tfac5  = 0.5e0_rt * tfac2;
    constexpr amrex::Real tfac6  = cv*cv + 1.5e0_rt*ca*ca + (xnufam - 1.0e0_rt)*(cvp*cvp + 1.5e0_rt*cap*cap);

    // N.B. if these are passed to an autodiff function (either admath::* or
    // overloaded operators), they need to be cast to a prvalue for CUDA
    // compilation to work

}


template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
number_t ifermi12(const number_t& f)
{

    // this routine applies a rational function expansion to get the inverse
    // fermi-dirac integral of order 1/2 when it is equal to f.
    // maximum error is 4.19e-9_rt.   reference: antia apjs 84,101 1993

    // declare work variables
    number_t rn,den,ff;

    // the return value
    number_t ifermi12r;

    // load the coefficients of the expansion from Table 8 of Antia

    constexpr amrex::Real an{0.5e0_rt};
    constexpr int m1{4};
    constexpr int k1{3};
    constexpr int m2{6};
    constexpr int k2{5};

    const amrex::Array1D<amrex::Real, 1, 5> a1 = {
        1.999266880833e4_rt,
        5.702479099336e3_rt,
        6.610132843877e2_rt,
        3.818838129486e1_rt,
        1.0e0_rt};

    const amrex::Array1D<amrex::Real, 1, 4> b1 = {
        1.771804140488e4_rt,
        -2.014785161019e3_rt,
        9.130355392717e1_rt,
        -1.670718177489e0_rt};

    const amrex::Array1D<amrex::Real, 1, 7> a2 = {
        -1.277060388085e-2_rt,
        7.187946804945e-2_rt,
        -4.262314235106e-1_rt,
        4.997559426872e-1_rt,
        -1.285579118012e0_rt,
        -3.930805454272e-1_rt,
        1.0e0_rt};

    const amrex::Array1D<amrex::Real, 1, 6> b2 = {
        -9.745794806288e-3_rt,
        5.485432756838e-2_rt,
        -3.299466243260e-1_rt,
        4.077841975923e-1_rt,
        -1.145531476975e0_rt,
        -6.067091689181e-2_rt};

    if (f < 4.0e0_rt) {

        // build sum_{i=0, m1} a_i x**i.  This is the numerator
        // in Eq. 4 of Antia.
        //
        // with our 1-based indexing, this expression is
        // a[1] + a[2] * f + a[3] * f**2 + ... a[m1+1] * f**m1
        //
        // we do the sum starting with the largest term and working
        // on a single power of f each iteration.
        //
        // in the starting rn here, the leading f is actually
        // a1(m1+1) * f, but that element is 1
        rn  = f + a1(m1);

        for (int i = m1 - 1; i >= 1; --i) {
            rn  = rn*f + a1(i);
        }

        // now we do the denominator in Eq. 4.  None of the coefficients
        // are 1, so we loop over all
        den = b1(k1+1);

        for (int i = k1; i >= 1; --i) {
            den = den*f + b1(i);
        }

        // Eq. 6 of Antia

        ifermi12r = admath::log(f * rn/den);

    } else {

        ff = 1.0e0_rt/admath::pow(f, 1.0e0_rt/(1.0e0_rt + an));

        // this construction is the same as above, but using the
        // second set of coefficients

        rn = ff + a2(m2);

        for (int i = m2 - 1; i >= 1; --i) {
            rn = rn*ff + a2(i);
        }

        den = b2(k2+1);

        for (int i = k2; i >= 1; --i) {
            den = den*ff + b2(i);
        }

        ifermi12r = rn/(den*ff);

    }

    return ifermi12r;
}

/// type for values that only depend on temperature
template <typename number_t>
using dt_t = std::conditional_t<
    std::is_same_v<number_t, autodiff::dual_array<1, 3>>,
    autodiff::dual_array<1, 1>,
    number_t>;

/// type for values that only depend on composition
template <typename number_t>
using dcomp_t = std::conditional_t<
    std::is_same_v<number_t, autodiff::dual_array<1, 3>>,
    autodiff::dual_array<2, 3>,
    number_t>;

template <typename number_t>
struct sneutf_t {

    amrex::Real den;
    dt_t<number_t> temp;
    dcomp_t<number_t> zbar;

    amrex::Real deni;
    dcomp_t<number_t> abari;

    dcomp_t<number_t> ye;

    dt_t<number_t> t9;
    dt_t<number_t> xl;
    dt_t<number_t> xl2;
    dt_t<number_t> xl3;
    dt_t<number_t> xl4;
    dt_t<number_t> xl5;
    dt_t<number_t> xl6;
    dt_t<number_t> xl8;
    dt_t<number_t> xl9;
    dt_t<number_t> xlm1;
    dt_t<number_t> xlm2;
    dt_t<number_t> xlm3;
    dcomp_t<number_t> rm;

    number_t zeta;
    number_t zeta2;
    number_t zeta3;
};


template <typename number_t>
AMREX_GPU_HOST_DEVICE inline
sneutf_t<number_t> get_sneut_factors(const amrex::Real den, const number_t& temp,
                                     const number_t& abar, const number_t& zbar) {

    using autodiff::eval;

    constexpr amrex::Real con1   = 1.0e0_rt/5.9302e0_rt;

    sneutf_t<number_t> sf{};

    sf.den = den;
    sf.temp = microphysics_autodiff::narrow_array<dt_t<number_t>>(temp);
    sf.zbar = microphysics_autodiff::narrow_array<dcomp_t<number_t>>(zbar);

    // to avoid lots of divisions
    sf.deni  = 1.0e0_rt / den;
    sf.abari = 1.0e0_rt / microphysics_autodiff::narrow_array<dcomp_t<number_t>>(abar);

    // some composition variables
    sf.ye    = sf.zbar * sf.abari;

    // some frequent factors
    sf.t9     = sf.temp * 1.0e-9_rt;
    sf.xl     = sf.t9 * con1;
    sf.xl2    = sf.xl * sf.xl;
    sf.xl3    = sf.xl2 * sf.xl;
    sf.xl4    = sf.xl3 * sf.xl;
    sf.xl5    = sf.xl4 * sf.xl;
    sf.xl6    = sf.xl5 * sf.xl;
    const auto xl7 = eval(sf.xl6 * sf.xl);
    sf.xl8    = xl7 * sf.xl;
    sf.xl9    = sf.xl8 * sf.xl;
    sf.xlm1   = 1.0e0_rt / sf.xl;
    sf.xlm2   = sf.xlm1 * sf.xlm1;
    sf.xlm3   = sf.xlm1 * sf.xlm2;

    sf.rm     = den * sf.ye;

    const auto a0 = eval(sf.rm * 1.0e-9_rt);
    const auto a1 = eval(admath::pow(a0, amrex::Real(nu_constants::oneth)));
    sf.zeta = a1 * sf.xlm1;

    sf.zeta2 = sf.zeta * sf.zeta;
    sf.zeta3 = sf.zeta2 * sf.zeta;

    return sf;

}



template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
number_t nu_pair(const sneutf_t<number_t>& sf) {

    // pair neutrino section
    // for reactions like e+ + e- => nu_e + nubar_e

    using autodiff::eval;

    // equation 2.8
    const auto gl = eval(1.0e0_rt - 13.04e0_rt * sf.xl2 +133.5e0_rt * sf.xl4 +1534.0e0_rt * sf.xl6 + 918.6e0_rt * sf.xl8);

    // equation 2.7

    const auto tmp_1 = eval(6.002e19_rt + 2.084e20_rt * sf.zeta + 1.872e21_rt * sf.zeta2);
    number_t tmp_2;
    if (sf.t9 < 10.0_rt) {
        tmp_2 = admath::exp(-5.5924e0_rt * sf.zeta);
    } else {
        tmp_2 = admath::exp(-4.9924e0_rt * sf.zeta);
    }

    const auto xnum = eval(tmp_1 * tmp_2);

    dt_t<number_t> tmp_3;
    if (sf.t9 < 10.0_rt) {
        tmp_3 = 9.383e-1_rt * sf.xlm1 - 4.141e-1_rt * sf.xlm2 + 5.829e-2_rt * sf.xlm3;
    } else {
        tmp_3 = 1.2383e0_rt * sf.xlm1 - 8.141e-1_rt * sf.xlm2;
    }

    const auto xden = eval(sf.zeta3 + tmp_3);

    const auto tmp_4 = eval(1.0e0_rt / xden);
    const auto fpair = eval(xnum * tmp_4);

    // equation 2.6
    const auto xlp5 = eval(admath::sqrt(sf.xl));
    const auto tmp_5 = eval(10.7480e0_rt * sf.xl2 + 0.3967e0_rt * xlp5 + 1.005e0_rt);

    const auto c = eval(1.0e0_rt / (7.692e7_rt * sf.xl3 + 9.715e6_rt * xlp5));
    const auto tmp_7 = eval(1.0e0_rt + sf.rm * c);

    const auto qpair = eval(1.0e0_rt / tmp_5 * admath::pow(tmp_7, -0.3e0_rt));

    // equation 2.5
    number_t spair = admath::exp(-2.0e0_rt * sf.xlm1) * fpair * gl;

    spair *= amrex::Real(nu_constants::tfac4) * (1.0e0_rt + amrex::Real(nu_constants::tfac3) * qpair);

    return spair;

}

template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
number_t nu_plasma(const sneutf_t<number_t>& sf) {

    // plasma neutrino section
    // for collective reactions like gamma_plasmon => nu_e + nubar_e
    // equation 4.6

    using autodiff::eval;

    const auto tmp_1 = eval(admath::pow(1.019e-6_rt * sf.rm, amrex::Real(nu_constants::twoth)));

    const auto tmp_2 = eval(1.0e0_rt / (sf.temp * sf.temp * admath::sqrt(1.0e0_rt + tmp_1)));

    const auto gl2   = eval(1.1095e11_rt * sf.rm * tmp_2);

    const auto gl = eval(admath::sqrt(gl2));
    const auto gl12 = eval(admath::sqrt(gl));
    const auto gl32 = eval(gl * gl12);
    const auto gl72 = eval(gl2 * gl32);
    const auto gl6 = eval(gl2 * gl2 * gl2);

    // equation 4.7
    const auto ft = eval(2.4e0_rt + 0.6e0_rt * gl12 + 0.51e0_rt * gl + 1.25e0_rt * gl32);

    // equation 4.8
    const auto tmp_3 = eval(8.6e0_rt * gl2 + 1.35e0_rt * gl72);
    const auto tmp_4 = eval(225.0e0_rt - 17.0e0_rt * gl + gl2);

    const auto fl = eval(tmp_3 * (1.0e0_rt / tmp_4));

    // equation 4.9 and 4.10
    const auto cc = eval(admath::log10(2.0e0_rt * sf.rm));
    const auto xlnt = eval(admath::log10(sf.temp));

    const auto xnum = eval(amrex::Real(nu_constants::sixth) * (17.5e0_rt + cc - 3.0e0_rt * xlnt));
    const auto xden = eval(amrex::Real(nu_constants::sixth) * (-24.5e0_rt + cc + 3.0e0_rt * xlnt));

    // equation 4.11
    number_t fxy;
    if (admath::abs(xnum) > 0.7e0_rt || xden < 0.0e0_rt) {
        fxy   = 1.0e0_rt;
    } else {

        const auto c = eval(admath::min(0.0e0_rt, xden - 1.6e0_rt + 1.25e0_rt * xnum));
        const auto a3 = eval(c / (0.57e0_rt - 0.25e0_rt * xnum));

        const auto a1 = eval(0.39e0_rt - 1.25e0_rt * xnum - 0.35e0_rt * admath::sin(4.5e0_rt * xnum));
        const auto b1 = eval(0.3e0_rt * admath::exp(-admath::powi<2>(4.5e0_rt*xnum + 0.9e0_rt)));
        fxy = 1.05e0_rt + (a1 - b1) * admath::exp(-a3 * a3);
    }

    // equation 4.1 and 4.5
    number_t splas = eval((ft + fl) * fxy * admath::exp(-gl) * gl6);
    splas *= eval(0.93153e0_rt * 3.0e21_rt * sf.xl9);

    return splas;

}

template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
number_t nu_photo(const sneutf_t<number_t>& sf) {

    // photoneutrino process section
    // for reactions like e- + gamma => e- + nu_e + nubar_e
    //                    e+ + gamma => e+ + nu_e + nubar_e
    // equation 3.8 for tau, equation 3.6 for cc,
    // and table 2 written out for speed

    using autodiff::eval;

    dt_t<number_t> tau;
    dt_t<number_t> cc;
    amrex::Real c00, c01, c02, c03, c04, c05, c06;
    amrex::Real c10, c11, c12, c13, c14, c15, c16;
    amrex::Real c20, c21, c22, c23, c24, c25, c26;
    amrex::Real dd01, dd02, dd03, dd04, dd05;
    amrex::Real dd11, dd12, dd13, dd14, dd15;
    amrex::Real dd21, dd22, dd23, dd24, dd25;

    if (sf.temp < 1.0e8_rt) {

        // note: we already bailed above for T < 1.e7, so this is really 1.e7 <= T < 1.e8

        tau  =  admath::log10(sf.temp * 1.0e-7_rt);
        cc   =  0.5654e0_rt + tau;
        c00  =  1.008e11_rt;
        c01  =  0.0e0_rt;
        c02  =  0.0e0_rt;
        c03  =  0.0e0_rt;
        c04  =  0.0e0_rt;
        c05  =  0.0e0_rt;
        c06  =  0.0e0_rt;
        c10  =  8.156e10_rt;
        c11  =  9.728e8_rt;
        c12  = -3.806e9_rt;
        c13  = -4.384e9_rt;
        c14  = -5.774e9_rt;
        c15  = -5.249e9_rt;
        c16  = -5.153e9_rt;
        c20  =  1.067e11_rt;
        c21  = -9.782e9_rt;
        c22  = -7.193e9_rt;
        c23  = -6.936e9_rt;
        c24  = -6.893e9_rt;
        c25  = -7.041e9_rt;
        c26  = -7.193e9_rt;
        dd01 =  0.0e0_rt;
        dd02 =  0.0e0_rt;
        dd03 =  0.0e0_rt;
        dd04 =  0.0e0_rt;
        dd05 =  0.0e0_rt;
        dd11 = -1.879e10_rt;
        dd12 = -9.667e9_rt;
        dd13 = -5.602e9_rt;
        dd14 = -3.370e9_rt;
        dd15 = -1.825e9_rt;
        dd21 = -2.919e10_rt;
        dd22 = -1.185e10_rt;
        dd23 = -7.270e9_rt;
        dd24 = -4.222e9_rt;
        dd25 = -1.560e9_rt;

    } else if (sf.temp >= 1.0e8_rt && sf.temp < 1.0e9_rt) {

        tau  =  admath::log10(sf.temp * 1.0e-8_rt);
        cc   =  1.5654e0_rt;
        c00  =  9.889e10_rt;
        c01  = -4.524e8_rt;
        c02  = -6.088e6_rt;
        c03  =  4.269e7_rt;
        c04  =  5.172e7_rt;
        c05  =  4.910e7_rt;
        c06  =  4.388e7_rt;
        c10  =  1.813e11_rt;
        c11  = -7.556e9_rt;
        c12  = -3.304e9_rt;
        c13  = -1.031e9_rt;
        c14  = -1.764e9_rt;
        c15  = -1.851e9_rt;
        c16  = -1.928e9_rt;
        c20  =  9.750e10_rt;
        c21  =  3.484e10_rt;
        c22  =  5.199e9_rt;
        c23  = -1.695e9_rt;
        c24  = -2.865e9_rt;
        c25  = -3.395e9_rt;
        c26  = -3.418e9_rt;
        dd01 = -1.135e8_rt;
        dd02 =  1.256e8_rt;
        dd03 =  5.149e7_rt;
        dd04 =  3.436e7_rt;
        dd05 =  1.005e7_rt;
        dd11 =  1.652e9_rt;
        dd12 = -3.119e9_rt;
        dd13 = -1.839e9_rt;
        dd14 = -1.458e9_rt;
        dd15 = -8.956e8_rt;
        dd21 = -1.548e10_rt;
        dd22 = -9.338e9_rt;
        dd23 = -5.899e9_rt;
        dd24 = -3.035e9_rt;
        dd25 = -1.598e9_rt;

    } else {

        // T > 1.e9

        tau  =  admath::log10(sf.t9);
        cc   =  1.5654e0_rt;
        c00  =  9.581e10_rt;
        c01  =  4.107e8_rt;
        c02  =  2.305e8_rt;
        c03  =  2.236e8_rt;
        c04  =  1.580e8_rt;
        c05  =  2.165e8_rt;
        c06  =  1.721e8_rt;
        c10  =  1.459e12_rt;
        c11  =  1.314e11_rt;
        c12  = -1.169e11_rt;
        c13  = -1.765e11_rt;
        c14  = -1.867e11_rt;
        c15  = -1.983e11_rt;
        c16  = -1.896e11_rt;
        c20  =  2.424e11_rt;
        c21  = -3.669e9_rt;
        c22  = -8.691e9_rt;
        c23  = -7.967e9_rt;
        c24  = -7.932e9_rt;
        c25  = -7.987e9_rt;
        c26  = -8.333e9_rt;
        dd01 =  4.724e8_rt;
        dd02 =  2.976e8_rt;
        dd03 =  2.242e8_rt;
        dd04 =  7.937e7_rt;
        dd05 =  4.859e7_rt;
        dd11 = -7.094e11_rt;
        dd12 = -3.697e11_rt;
        dd13 = -2.189e11_rt;
        dd14 = -1.273e11_rt;
        dd15 = -5.705e10_rt;
        dd21 = -2.254e10_rt;
        dd22 = -1.551e10_rt;
        dd23 = -7.793e9_rt;
        dd24 = -4.489e9_rt;
        dd25 = -2.185e9_rt;

    }

    // equation 3.7

    const auto [sin1, cos1] = admath::sincos(amrex::Real(nu_constants::fac1) * tau);

    const auto cos1_sq = eval(cos1 * cos1);
    const auto sin1_sq = eval(sin1 * sin1);

    // double, triple, etc. angle formulas
    // sin/cos (2 fac1 tau)
    const auto sin2 = eval(2.0_rt * sin1 * cos1);
    const auto cos2 = eval(2.0_rt * cos1_sq - 1.0_rt);

    // sin/cos (3 fac1 tau)
    const auto sin3 = eval(sin1 * (3.0_rt - 4.0_rt * sin1 * sin1));
    const auto cos3 = eval(cos1 * (4.0_rt * cos1_sq - 3.0_rt));

    // sin/cos (4 fac1 tau) -- use double angle on sin2/cos2
    const auto sin4 = eval(2.0_rt * sin2 * cos2);
    const auto cos4 = eval(2.0_rt * cos2 * cos2 - 1.0_rt);

    // sin/cos (5 fac1 tau)
    const auto sin5 = eval(sin1 * (5.0_rt - sin1_sq * (20.0_rt - 16.0_rt * sin1_sq)));
    const auto cos5 = eval(cos1 * (cos1_sq * (16.0_rt * cos1_sq - 20.0_rt) + 5.0_rt));

    // cos (10 pi tau)
    const auto last = eval(admath::cos(amrex::Real(nu_constants::fac2) * tau));

    const auto a0 = eval(0.5e0_rt * c00
        + c01 * cos1 + dd01 * sin1 + c02 * cos2 + dd02 * sin2
        + c03 * cos3 + dd03 * sin3 + c04 * cos4 + dd04 * sin4
        + c05 * cos5 + dd05 * sin5 + 0.5e0_rt * c06 * last);
    const auto a1 = eval(0.5e0_rt * c10
        + c11 * cos1 + dd11 * sin1 + c12 * cos2 + dd12 * sin2
        + c13 * cos3 + dd13 * sin3 + c14 * cos4 + dd14 * sin4
        + c15 * cos5 + dd15 * sin5 + 0.5e0_rt * c16 * last);
    const auto a2 = eval(0.5e0_rt * c20
        + c21 * cos1 + dd21 * sin1 + c22 * cos2 + dd22 * sin2
        + c23 * cos3 + dd23 * sin3 + c24 * cos4 + dd24 * sin4
        + c25 * cos5 + dd25 * sin5 + 0.5e0_rt * c26 * last);

    // equation 3.4
    const auto z = eval(admath::exp(-cc * sf.zeta));

    const auto fnum = eval((a0 + a1 * sf.zeta + a2 * sf.zeta2) * z);

    const auto fden = eval(sf.zeta3 + 6.290e-3_rt * sf.xlm1 + 7.483e-3_rt * sf.xlm2 + 3.061e-4_rt * sf.xlm3);

    const auto fphot = eval(fnum * (1.0e0_rt / fden));

    // equation 3.3
    const auto qnum = eval(0.666e0_rt * admath::pow(1.0e0_rt + 2.045e0_rt * sf.xl, -2.066e0_rt));

    const auto dum = eval(1.875e8_rt * sf.xl + 1.653e8_rt * sf.xl2 + 8.499e8_rt * sf.xl3 - 1.604e8_rt * sf.xl4);

    const auto qden = eval(1.0e0_rt + sf.rm * (1.0e0_rt / dum));

    const auto qphot = eval(qnum * (1.0e0_rt / qden));

    // equation 3.2
    number_t sphot = sf.xl5 * fphot * sf.rm;
    sphot *= amrex::Real(nu_constants::tfac4) * (1.0e0_rt - amrex::Real(nu_constants::tfac3) * qphot);

    if (sphot <= 0.0_rt) {
       sphot   = 0.0e0_rt;
    }

    return sphot;
}

template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
number_t nu_brem(const sneutf_t<number_t>& sf) {

    // bremsstrahlung neutrino section
    // for reactions like e- + (z,a) => e- + (z,a) + nu + nubar
    //                    n  + n     => n + n + nu + nubar
    //                    n  + p     => n + p + nu + nubar
    // equation 4.3

    using autodiff::eval;

    const amrex::Real den6 = sf.den * 1.0e-6_rt;
    const auto t8      = eval(sf.temp * 1.0e-8_rt);
    const auto t812    = eval(admath::sqrt(t8));
    const auto t832    = eval(t8 * t812);
    const auto t82     = eval(t8*t8);
    const auto t83     = eval(t82*t8);
    const auto t85     = eval(t82*t83);
    const auto t86     = eval(t85*t8);
    const auto t8m1    = eval(1.0e0_rt/t8);
    const auto t8m2    = eval(t8m1*t8m1);
    const auto t8m3    = eval(t8m2*t8m1);
    const auto t8m5    = eval(t8m3*t8m2);

    const amrex::Real tfermi = 5.9302e9_rt * (std::sqrt(1.0e0_rt + 1.018e0_rt * std::pow(den6 * autodiff::val(sf.ye), nu_constants::twoth)) - 1.0e0_rt);

    number_t sbrem{};

    // "weak" degenerate electrons only
    if (sf.temp > 0.3e0_rt * tfermi) {

        // equation 5.3
        const auto tmp_1 = eval(7.05e6_rt * t832 + 5.12e4_rt * t83);

        const auto tmp_1a = eval(1.0e0_rt / tmp_1);
        const auto eta = eval(sf.rm * tmp_1a);

        const auto etam1 = eval(1.0e0_rt/eta);
        const auto etam2 = eval(etam1 * etam1);

        // equation 5.2
        const auto tmp_2 = eval(23.5e0_rt + 6.83e4_rt * t8m2 + 7.81e8_rt * t8m5);
        const auto tmp_3 = eval(1.0e0_rt / tmp_2);

        const auto tmp_4a = eval(1.0e0_rt + 1.47e0_rt * etam1 + 3.29e-2_rt * etam2);

        const auto tmp_5a = eval(1.26e0_rt * (1.0e0_rt + etam1));

        const auto tmp_4 = eval(1.0e0_rt/tmp_4a);
        const auto tmp_5 = eval(tmp_5a * tmp_4);

        const auto fbrem = eval(tmp_3 + tmp_5);

        // equation 5.9
        const auto tmp_6 = eval(230.0e0_rt + 6.7e5_rt * t8m2 + 7.66e9_rt * t8m5);

        const auto tmp_7 = eval(1.0e0_rt + sf.rm * 1.0e-9_rt);
        const auto tmp_8a = eval(tmp_6 * tmp_7);

        const auto tmp_8 = eval(1.0e0_rt / tmp_8a);

        const auto c00 = eval(7.75e5_rt * t832 + 247.0e0_rt * admath::pow(t8, 3.85e0_rt));
        const auto c01 = eval(4.07e0_rt + 0.0240e0_rt * admath::pow(t8, 1.4e0_rt));
        const auto c02 = eval(4.59e-5_rt * admath::pow(t8, -0.110e0_rt));

        const auto rmi = eval(1.0e0_rt / sf.rm);
        const auto tmp_9 = eval(c00 * rmi + c01 + c02 * std::pow(sf.den, 0.656e0_rt));

        const auto tmp_9a = eval(1.0e0_rt / tmp_9);

        const auto gbrem = eval(tmp_8 + tmp_9a);

        // equation 5.1
        const auto tmp_10 = eval(0.5738e0_rt * sf.zbar * sf.ye * t86 * sf.den);

        const auto tmp_11 = eval(amrex::Real(nu_constants::tfac4) * fbrem - amrex::Real(nu_constants::tfac5) * gbrem);
        sbrem = tmp_10 * tmp_11;

    } else {

        // liquid metal with c12 parameters (not too different for other elements)
        // equation 5.18 and 5.16

        const amrex::Real u = amrex::Real(nu_constants::fac3) * (std::log10(sf.den) - 3.0e0_rt);

        // compute the expensive trig functions of equation 5.21 only once

        const auto [sin1, cos1] = amrex::Math::sincos(u);

        // double, triple, etc. angle formulas
        // sin/cos (2 u)
        const amrex::Real sin2 = 2.0_rt * sin1 * cos1;
        const amrex::Real cos2 = 2.0_rt * cos1 * cos1 - 1.0_rt;

        // sin/cos (3 u)
        const amrex::Real sin3 = sin1 * (3.0_rt - 4.0_rt * sin1 * sin1);
        const amrex::Real cos3 = cos1 * (4.0_rt * cos1 * cos1 - 3.0_rt);

        // sin/cos (4 u) -- use double angle on sin2/cos2
        const amrex::Real sin4 = 2.0_rt * sin2 * cos2;
        const amrex::Real cos4 = 2.0_rt * cos2 * cos2 - 1.0_rt;

        // sin/cos (5 u)
        const amrex::Real cos5 = cos1 * (cos1 * cos1 * (16.0_rt * cos1 * cos1 - 20.0_rt) + 5.0_rt);

        // equation 5.21
        const amrex::Real fb = 0.5e0_rt * 0.17946e0_rt + 0.00945e0_rt * u + 0.34529e0_rt
            - 0.05821e0_rt * cos1 - 0.04969e0_rt * sin1
            - 0.01089e0_rt * cos2 - 0.01584e0_rt * sin2
            - 0.01147e0_rt * cos3 - 0.00504e0_rt * sin3
            - 0.00656e0_rt * cos4 - 0.00281e0_rt * sin4
            - 0.00519e0_rt * cos5;

        // equation 5.22
        const amrex::Real ft = 0.5e0_rt * 0.06781e0_rt - 0.02342e0_rt * u + 0.24819e0_rt
            - 0.00944e0_rt * cos1 - 0.02213e0_rt * sin1
            - 0.01289e0_rt * cos2 - 0.01136e0_rt * sin2
            - 0.00589e0_rt * cos3 - 0.00467e0_rt * sin3
            - 0.00404e0_rt * cos4 - 0.00131e0_rt * sin4
            - 0.00330e0_rt * cos5;

        // equation 5.23
        const amrex::Real gb = 0.5e0_rt * 0.00766e0_rt - 0.01259e0_rt * u + 0.07917e0_rt
            - 0.00710e0_rt * cos1 + 0.02300e0_rt * sin1
            - 0.00028e0_rt * cos2 - 0.01078e0_rt * sin2
            + 0.00232e0_rt * cos3 + 0.00118e0_rt * sin3
            + 0.00044e0_rt * cos4 - 0.00089e0_rt * sin4
            + 0.00158e0_rt * cos5;

        // equation 5.24
        const amrex::Real gt = -0.5e0_rt * 0.00769e0_rt  - 0.00829e0_rt * u + 0.05211e0_rt
            + 0.00356e0_rt * cos1 + 0.01052e0_rt * sin1
            - 0.00184e0_rt * cos2 - 0.00354e0_rt * sin2
            + 0.00146e0_rt * cos3 - 0.00014e0_rt * sin3
            + 0.00031e0_rt * cos4 - 0.00018e0_rt * sin4
            + 0.00069e0_rt * cos5;

        const auto dum = eval(2.275e-1_rt * sf.zbar * sf.zbar * t8m1 * admath::pow(den6*sf.abari, amrex::Real(nu_constants::oneth)));

        const auto gm1 = eval(1.0e0_rt / dum);
        const auto gm13 = eval(admath::pow(gm1, amrex::Real(nu_constants::oneth)));
        const auto gm23 = eval(gm13 * gm13);

        // equation 5.25 and 5.26
        const auto v = eval(-0.05483e0_rt - 0.01946e0_rt * gm13 + 1.86310e0_rt * gm23 - 0.78873e0_rt * gm1);

        const auto w = eval(-0.06711e0_rt + 0.06859e0_rt * gm13 + 1.74360e0_rt * gm23 - 0.74498e0_rt * gm1);

        // equation 5.19 and 5.20
        const auto fliq = eval(v * fb + (1.0e0_rt - v) * ft);

        const auto gliq = eval(w * gb + (1.0e0_rt - w) * gt);

        // equation 5.17
        const auto dum_2 = eval(0.5738e0_rt * sf.zbar * sf.ye * t86 * sf.den);

        const auto z  = eval(amrex::Real(nu_constants::tfac4)*fliq - amrex::Real(nu_constants::tfac5)*gliq);
        sbrem = dum_2 * z;
    }
    return sbrem;

}

template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
number_t nu_recomb(const sneutf_t<number_t>& sf) {

    // recombination neutrino section
    // for reactions like e- (continuum) => e- (bound) + nu_e + nubar_e
    // equation 6.11 solved for nu

    using autodiff::eval;

    const auto xnum   = eval(1.10520e8_rt * sf.den * sf.ye / (sf.temp * admath::sqrt(sf.temp)));

    // the chemical potential
    const auto nu   = ifermi12(xnum);
    const auto nu2  = eval(nu * nu);
    const auto nu3  = eval(nu2 * nu);

    amrex::Real a1{}, a2{}, a3{}, b{}, c{}, d{}, f1{}, f2{}, f3{};

    // table 12
    if (nu >= -20.0_rt && nu < 0.0_rt) {
       a1 = 1.51e-2_rt;
       a2 = 2.42e-1_rt;
       a3 = 1.21e0_rt;
       b  = 3.71e-2_rt;
       c  = 9.06e-1_rt;
       d  = 9.28e-1_rt;
       f1 = 0.0e0_rt;
       f2 = 0.0e0_rt;
       f3 = 0.0e0_rt;
    } else if (nu >= 0.0_rt && nu <= 10.0_rt) {
       a1 = 1.23e-2_rt;
       a2 = 2.66e-1_rt;
       a3 = 1.30e0_rt;
       b  = 1.17e-1_rt;
       c  = 8.97e-1_rt;
       d  = 1.77e-1_rt;
       f1 = -1.20e-2_rt;
       f2 = 2.29e-2_rt;
       f3 = -1.04e-3_rt;
    }

    number_t sreco{};
    // equation 6.7, 6.13 and 6.12
    if (nu >= -20.0_rt &&  nu <= 10.0_rt) {

        // equation 6.7
        const auto tempi  = eval(1.0e0_rt / sf.temp);
        const auto zeta   = eval(1.579e5_rt * sf.zbar * sf.zbar * tempi);

        // equation 6.13
        const auto tmp_1 = eval(1.0e0_rt / (1.0e0_rt + f1 * nu + f2 * nu2 + f3 * nu3));
        const auto z = eval(zeta * tmp_1);

        // equation 6.12
        const auto tmp_2 = eval(1.0e0_rt / z);
        const auto tmp_3 = eval(admath::pow(z, -2.25_rt));
        const auto tmp_4 = eval(admath::pow(z, -4.55_rt));
        const auto tmp_5 = eval(a1 * tmp_2 + a2 * tmp_3 + a3 * tmp_4);

        const auto tmp_6 = eval(admath::exp(c * nu));
        const auto tmp_7 = eval(b * tmp_6 * (1.0e0_rt + d * z));
        const auto gum = eval(1.0e0_rt + tmp_7);

        const auto tmp_8 = eval(admath::exp(nu));
        const auto tmp_9 = eval(1.0e0_rt / gum);
        const auto bigj = eval(tmp_5 * tmp_8 * tmp_9);

        // equation 6.5
        const auto tmp_10 = eval(admath::exp(zeta + nu));
        const auto tmp_11 = eval(1.0e0_rt + tmp_10);
        const auto tmp_12 = eval(1.0e0_rt/tmp_11);

        sreco = amrex::Real(nu_constants::tfac6) * 2.649e-18_rt * sf.ye * admath::powi<13>(sf.zbar) * sf.den * bigj * tmp_12;
    }
    return sreco;

}

template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
number_t sneut5(const number_t& temp, const amrex::Real den,
                const number_t& abar, const number_t& zbar)
{
    /*
    this routine computes thermal neutrino losses from the analytic fits of
    itoh et al. apjs 102, 411, 1996, and also returns their derivatives.

    input:
    temp = temperature
    den  = density
    abar = mean atomic weight
    zbar = mean charge

    output:
    snu  = total neutrino loss rate in erg/g/sec
    */

    if (temp < 1.0e7_rt) {
        return 0.0_rt;
    }

    auto sf = get_sneut_factors(den, temp, abar, zbar);

    number_t spair = nu_pair(sf);
    number_t splas = nu_plasma(sf);
    number_t sphot = nu_photo(sf);
    number_t sbrem = nu_brem(sf);
    number_t sreco = 0.0_rt;
    if (neutrino_cooling_rp::include_recomb) {
        sreco = nu_recomb(sf);
    }

    // convert from erg/cm^3/s to erg/g/s
    // comment these out to duplicate the itoh et al plots

    spair   *= sf.deni;
    splas   *= sf.deni;
    sphot   *= sf.deni;
    sbrem   *= sf.deni;
    sreco   *= sf.deni;

    // the total neutrino loss rate
    number_t snu = splas + spair + sphot + sbrem;
    if (neutrino_cooling_rp::include_recomb) {
        snu += sreco;
    }
    return snu;
}

template <int do_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void sneut5(const amrex::Real temp, const amrex::Real den,
            const amrex::Real abar, const amrex::Real zbar,
            amrex::Real& snu, amrex::Real& dsnudt, amrex::Real& dsnudd,
            amrex::Real& dsnuda, amrex::Real& dsnudz)
{
    // autodiff wrapper
    if constexpr (do_derivatives) {
        using dual_t = autodiff::dual_array<1, 3>;
        dual_t temp_dual = temp;
        dual_t abar_dual = abar;
        dual_t zbar_dual = zbar;
        autodiff::seed_array(temp_dual, abar_dual, zbar_dual);
        dual_t snu_dual = sneut5(temp_dual, den, abar_dual, zbar_dual);
        snu = autodiff::val(snu_dual);
        const auto& grad = autodiff::derivative(snu_dual);
        dsnudt = grad(1);
        dsnudd = 0.0e0_rt;
        dsnuda = grad(2);
        dsnudz = grad(3);
    } else {
        snu = sneut5(temp, den, abar, zbar);
        dsnudt = 0.0e0_rt;
        dsnudd = 0.0e0_rt;
        dsnuda = 0.0e0_rt;
        dsnudz = 0.0e0_rt;
    }

}

#endif
