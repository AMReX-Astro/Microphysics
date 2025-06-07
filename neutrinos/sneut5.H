#ifndef SNEUT5_H
#define SNEUT5_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Math.H>

#include <microphysics_autodiff.H>

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

template <typename number_t>
struct sneutf_t {

    amrex::Real den;
    number_t temp;
    number_t zbar;

    amrex::Real deni;
    number_t abari;

    number_t ye;

    number_t t9;
    number_t xl;
    number_t xl2;
    number_t xl3;
    number_t xl4;
    number_t xl5;
    number_t xl6;
    number_t xl8;
    number_t xl9;
    number_t xlm1;
    number_t xlm2;
    number_t xlm3;
    number_t rm;

    number_t zeta;
    number_t zeta2;
    number_t zeta3;
};


template <typename number_t>
AMREX_GPU_HOST_DEVICE inline
sneutf_t<number_t> get_sneut_factors(const amrex::Real den, const number_t& temp,
                                     const number_t& abar, const number_t& zbar) {

    constexpr amrex::Real con1   = 1.0e0_rt/5.9302e0_rt;

    sneutf_t<number_t> sf{};

    sf.den = den;
    sf.temp = temp;
    sf.zbar = zbar;

    // to avoid lots of divisions
    sf.deni  = 1.0e0_rt / den;
    sf.abari = 1.0e0_rt / abar;

    // some composition variables
    sf.ye    = zbar * sf.abari;

    // some frequent factors
    sf.t9     = temp * 1.0e-9_rt;
    sf.xl     = sf.t9 * con1;
    sf.xl2    = sf.xl * sf.xl;
    sf.xl3    = sf.xl2 * sf.xl;
    sf.xl4    = sf.xl3 * sf.xl;
    sf.xl5    = sf.xl4 * sf.xl;
    sf.xl6    = sf.xl5 * sf.xl;
    number_t xl7 = sf.xl6 * sf.xl;
    sf.xl8    = xl7 * sf.xl;
    sf.xl9    = sf.xl8 * sf.xl;
    sf.xlm1   = 1.0e0_rt / sf.xl;
    sf.xlm2   = sf.xlm1 * sf.xlm1;
    sf.xlm3   = sf.xlm1 * sf.xlm2;

    sf.rm     = den * sf.ye;

    number_t a0 = sf.rm * 1.0e-9_rt;
    number_t a1 = admath::pow(a0, amrex::Real(nu_constants::oneth));
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

    // equation 2.8
    number_t gl = 1.0e0_rt - 13.04e0_rt * sf.xl2 +133.5e0_rt * sf.xl4 +1534.0e0_rt * sf.xl6 + 918.6e0_rt * sf.xl8;

    // equation 2.7

    number_t a1 = 6.002e19_rt + 2.084e20_rt * sf.zeta + 1.872e21_rt * sf.zeta2;
    number_t b1;
    if (sf.t9 < 10.0_rt) {
        b1 = admath::exp(-5.5924e0_rt * sf.zeta);
    } else {
        b1 = admath::exp(-4.9924e0_rt * sf.zeta);
    }

    number_t xnum = a1 * b1;

    if (sf.t9 < 10.0_rt) {
        a1 = 9.383e-1_rt * sf.xlm1 - 4.141e-1_rt * sf.xlm2 + 5.829e-2_rt * sf.xlm3;
    } else {
        a1 = 1.2383e0_rt * sf.xlm1 - 8.141e-1_rt * sf.xlm2;
    }

    number_t xden = sf.zeta3 + a1;

    a1 = 1.0e0_rt / xden;
    number_t fpair = xnum * a1;

    // equation 2.6
    number_t xlp5 = admath::sqrt(sf.xl);
    a1 = 10.7480e0_rt * sf.xl2 + 0.3967e0_rt * xlp5 + 1.005e0_rt;
    xnum = 1.0e0_rt / a1;

    a1 = 7.692e7_rt * sf.xl3 + 9.715e6_rt * xlp5;
    number_t c = 1.0e0_rt / a1;
    b1 = 1.0e0_rt + sf.rm * c;

    xden = admath::pow(b1, -0.3e0_rt);

    number_t qpair = xnum * xden;

    // equation 2.5
    a1 = admath::exp(-2.0e0_rt * sf.xlm1);

    number_t spair = a1 * fpair;
    spair *= gl;

    a1 = amrex::Real(nu_constants::tfac4) * (1.0e0_rt + amrex::Real(nu_constants::tfac3) * qpair);
    spair *= a1;

    return spair;

}

template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
number_t nu_plasma(const sneutf_t<number_t>& sf) {

    // plasma neutrino section
    // for collective reactions like gamma_plasmon => nu_e + nubar_e
    // equation 4.6

    number_t a1 = 1.019e-6_rt * sf.rm;
    number_t a2 = admath::pow(a1, amrex::Real(nu_constants::twoth));

    number_t b1 = admath::sqrt(1.0e0_rt + a2);

    number_t c00 = 1.0e0_rt / (sf.temp * sf.temp * b1);

    number_t gl2   = 1.1095e11_rt * sf.rm * c00;

    number_t gl = admath::sqrt(gl2);
    number_t gl12 = admath::sqrt(gl);
    number_t gl32 = gl * gl12;
    number_t gl72 = gl2 * gl32;
    number_t gl6 = gl2 * gl2 * gl2;

    // equation 4.7
    number_t ft = 2.4e0_rt + 0.6e0_rt * gl12 + 0.51e0_rt * gl + 1.25e0_rt * gl32;

    // equation 4.8
    a1 = 8.6e0_rt * gl2 + 1.35e0_rt * gl72;
    b1 = 225.0e0_rt - 17.0e0_rt * gl + gl2;

    number_t c = 1.0e0_rt / b1;
    number_t fl = a1 * c;

    // equation 4.9 and 4.10
    number_t cc = admath::log10(2.0e0_rt * sf.rm);
    number_t xlnt = admath::log10(sf.temp);

    number_t xnum = amrex::Real(nu_constants::sixth) * (17.5e0_rt + cc - 3.0e0_rt * xlnt);
    number_t xden = amrex::Real(nu_constants::sixth) * (-24.5e0_rt + cc + 3.0e0_rt * xlnt);

    // equation 4.11
    number_t fxy;
    if (admath::abs(xnum) > 0.7e0_rt || xden < 0.0e0_rt) {
        fxy   = 1.0e0_rt;
    } else {

        a1 = 0.39e0_rt - 1.25e0_rt * xnum - 0.35e0_rt * admath::sin(4.5e0_rt * xnum);

        b1 = 0.3e0_rt * admath::exp(-amrex::Math::powi<2>(autodiff::eval(4.5e0_rt*xnum + 0.9e0_rt)));

        c = admath::min(0.0e0_rt, xden - 1.6e0_rt + 1.25e0_rt * xnum);

        number_t d = 0.57e0_rt - 0.25e0_rt * xnum;
        number_t a3 = c / d;
        c00 = admath::exp(-a3 * a3);

        fxy = 1.05e0_rt + (a1 - b1) * c00;
    }

    // equation 4.1 and 4.5
    number_t splas = (ft + fl) * fxy;

    a2 = admath::exp(-gl);
    splas *= a2;
    splas *= gl6;

    a2 = 0.93153e0_rt * 3.0e21_rt * sf.xl9;
    splas *= a2;

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

    number_t tau;
    number_t cc;
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

    // double, triple, etc. angle formulas
    // sin/cos (2 fac1 tau)
    const number_t sin2 = 2.0_rt * sin1 * cos1;
    const number_t cos2 = 2.0_rt * cos1 * cos1 - 1.0_rt;

    // sin/cos (3 fac1 tau)
    const number_t sin3 = sin1 * (3.0_rt - 4.0_rt * sin1 * sin1);
    const number_t cos3 = cos1 * (4.0_rt * cos1 * cos1 - 3.0_rt);

    // sin/cos (4 fac1 tau) -- use double angle on sin2/cos2
    const number_t sin4 = 2.0_rt * sin2 * cos2;
    const number_t cos4 = 2.0_rt * cos2 * cos2 - 1.0_rt;

    // sin/cos (5 fac1 tau)
    const number_t sin5 = sin1 * (5.0_rt - sin1 * sin1 * (20.0_rt - 16.0_rt * sin1 * sin1));
    const number_t cos5 = cos1 * (cos1 * cos1 * (16.0_rt * cos1 * cos1 - 20.0_rt) + 5.0_rt);

    // cos (10 pi tau)
    const number_t last = admath::cos(amrex::Real(nu_constants::fac2) * tau);

    number_t a0 = 0.5e0_rt * c00
        + c01 * cos1 + dd01 * sin1 + c02 * cos2 + dd02 * sin2
        + c03 * cos3 + dd03 * sin3 + c04 * cos4 + dd04 * sin4
        + c05 * cos5 + dd05 * sin5 + 0.5e0_rt * c06 * last;
    number_t a1 = 0.5e0_rt * c10
        + c11 * cos1 + dd11 * sin1 + c12 * cos2 + dd12 * sin2
        + c13 * cos3 + dd13 * sin3 + c14 * cos4 + dd14 * sin4
        + c15 * cos5 + dd15 * sin5 + 0.5e0_rt * c16 * last;
    number_t a2 = 0.5e0_rt * c20
        + c21 * cos1 + dd21 * sin1 + c22 * cos2 + dd22 * sin2
        + c23 * cos3 + dd23 * sin3 + c24 * cos4 + dd24 * sin4
        + c25 * cos5 + dd25 * sin5 + 0.5e0_rt * c26 * last;

    // equation 3.4
    number_t dum = a0 + a1 * sf.zeta + a2 * sf.zeta2;

    number_t z = admath::exp(-cc * sf.zeta);

    number_t xnum = dum * z;

    number_t xden = sf.zeta3 + 6.290e-3_rt * sf.xlm1 + 7.483e-3_rt * sf.xlm2 + 3.061e-4_rt * sf.xlm3;

    dum = 1.0e0_rt / xden;
    number_t fphot = xnum * dum;

    // equation 3.3
    a0 = 1.0e0_rt + 2.045e0_rt * sf.xl;
    xnum = 0.666e0_rt * admath::pow(a0, -2.066e0_rt);

    dum = 1.875e8_rt * sf.xl + 1.653e8_rt * sf.xl2 + 8.499e8_rt * sf.xl3 - 1.604e8_rt * sf.xl4;

    z  = 1.0e0_rt / dum;
    xden = 1.0e0_rt + sf.rm * z;

    z  = 1.0e0_rt / xden;
    number_t qphot = xnum * z;

    // equation 3.2
    number_t sphot = sf.xl5 * fphot;

    sphot *= sf.rm;

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

    amrex::Real den6 = sf.den * 1.0e-6_rt;
    number_t t8      = sf.temp * 1.0e-8_rt;
    number_t t812    = admath::sqrt(t8);
    number_t t832    = t8 * t812;
    number_t t82     = t8*t8;
    number_t t83     = t82*t8;
    number_t t85     = t82*t83;
    number_t t86     = t85*t8;
    number_t t8m1    = 1.0e0_rt/t8;
    number_t t8m2    = t8m1*t8m1;
    number_t t8m3    = t8m2*t8m1;
    number_t t8m5    = t8m3*t8m2;

    number_t tfermi = 5.9302e9_rt * (admath::sqrt(1.0e0_rt + 1.018e0_rt * admath::pow(den6 * sf.ye, amrex::Real(nu_constants::twoth))) - 1.0e0_rt);

    number_t sbrem{};

    // "weak" degenerate electrons only
    if (sf.temp > 0.3e0_rt * tfermi) {

        // equation 5.3
        number_t dum = 7.05e6_rt * t832 + 5.12e4_rt * t83;

        number_t z = 1.0e0_rt / dum;
        number_t eta = sf.rm * z;

        number_t etam1 = 1.0e0_rt/eta;
        number_t etam2 = etam1 * etam1;

        // equation 5.2
        number_t a0 = 23.5e0_rt + 6.83e4_rt * t8m2 + 7.81e8_rt * t8m5;
        number_t xnum = 1.0e0_rt / a0;

        dum = 1.0e0_rt + 1.47e0_rt * etam1 + 3.29e-2_rt * etam2;

        number_t c00 = 1.26e0_rt * (1.0e0_rt + etam1);

        z = 1.0e0_rt/dum;
        number_t xden = c00 * z;

        number_t fbrem = xnum + xden;

        // equation 5.9
        a0 = 230.0e0_rt + 6.7e5_rt * t8m2 + 7.66e9_rt * t8m5;

        z = 1.0e0_rt + sf.rm * 1.0e-9_rt;
        dum = a0 * z;

        xnum = 1.0e0_rt / dum;

        c00 = 7.75e5_rt * t832 + 247.0e0_rt * admath::pow(t8, 3.85e0_rt);
        number_t c01 = 4.07e0_rt + 0.0240e0_rt * admath::pow(t8, 1.4e0_rt);
        number_t c02 = 4.59e-5_rt * admath::pow(t8, -0.110e0_rt);

        number_t rmi = 1.0e0_rt / sf.rm;
        dum = c00 * rmi + c01 + c02 * std::pow(sf.den, 0.656e0_rt);

        xden  = 1.0e0_rt / dum;

        number_t gbrem = xnum + xden;

        // equation 5.1
        dum = 0.5738e0_rt * sf.zbar * sf.ye * t86 * sf.den;

        z = amrex::Real(nu_constants::tfac4) * fbrem - amrex::Real(nu_constants::tfac5) * gbrem;
        sbrem = dum * z;

    } else {

        // liquid metal with c12 parameters (not too different for other elements)
        // equation 5.18 and 5.16

        amrex::Real u = amrex::Real(nu_constants::fac3) * (std::log10(sf.den) - 3.0e0_rt);

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
        amrex::Real fb = 0.5e0_rt * 0.17946e0_rt + 0.00945e0_rt * u + 0.34529e0_rt
            - 0.05821e0_rt * cos1 - 0.04969e0_rt * sin1
            - 0.01089e0_rt * cos2 - 0.01584e0_rt * sin2
            - 0.01147e0_rt * cos3 - 0.00504e0_rt * sin3
            - 0.00656e0_rt * cos4 - 0.00281e0_rt * sin4
            - 0.00519e0_rt * cos5;

        // equation 5.22
        amrex::Real ft = 0.5e0_rt * 0.06781e0_rt - 0.02342e0_rt * u + 0.24819e0_rt
            - 0.00944e0_rt * cos1 - 0.02213e0_rt * sin1
            - 0.01289e0_rt * cos2 - 0.01136e0_rt * sin2
            - 0.00589e0_rt * cos3 - 0.00467e0_rt * sin3
            - 0.00404e0_rt * cos4 - 0.00131e0_rt * sin4
            - 0.00330e0_rt * cos5;

        // equation 5.23
        amrex::Real gb = 0.5e0_rt * 0.00766e0_rt - 0.01259e0_rt * u + 0.07917e0_rt
            - 0.00710e0_rt * cos1 + 0.02300e0_rt * sin1
            - 0.00028e0_rt * cos2 - 0.01078e0_rt * sin2
            + 0.00232e0_rt * cos3 + 0.00118e0_rt * sin3
            + 0.00044e0_rt * cos4 - 0.00089e0_rt * sin4
            + 0.00158e0_rt * cos5;

        // equation 5.24
        amrex::Real gt = -0.5e0_rt * 0.00769e0_rt  - 0.00829e0_rt * u + 0.05211e0_rt
            + 0.00356e0_rt * cos1 + 0.01052e0_rt * sin1
            - 0.00184e0_rt * cos2 - 0.00354e0_rt * sin2
            + 0.00146e0_rt * cos3 - 0.00014e0_rt * sin3
            + 0.00031e0_rt * cos4 - 0.00018e0_rt * sin4
            + 0.00069e0_rt * cos5;

        number_t dum = 2.275e-1_rt * sf.zbar * sf.zbar * t8m1 * admath::pow(den6*sf.abari, amrex::Real(nu_constants::oneth));

        number_t gm1 = 1.0e0_rt / dum;
        number_t gm13 = admath::pow(gm1, amrex::Real(nu_constants::oneth));
        number_t gm23 = gm13 * gm13;

        // equation 5.25 and 5.26
        number_t v = -0.05483e0_rt - 0.01946e0_rt * gm13 + 1.86310e0_rt * gm23 - 0.78873e0_rt * gm1;

        number_t w = -0.06711e0_rt + 0.06859e0_rt * gm13 + 1.74360e0_rt * gm23 - 0.74498e0_rt * gm1;

        // equation 5.19 and 5.20
        number_t fliq = v * fb + (1.0e0_rt - v) * ft;

        number_t gliq = w * gb + (1.0e0_rt - w) * gt;

        // equation 5.17
        dum = 0.5738e0_rt * sf.zbar * sf.ye * t86 * sf.den;

        number_t z  = amrex::Real(nu_constants::tfac4)*fliq - amrex::Real(nu_constants::tfac5)*gliq;
        sbrem = dum * z;
    }
    return sbrem;

}

template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
number_t nu_recomb(const sneutf_t<number_t>& sf) {

    // recombination neutrino section
    // for reactions like e- (continuum) => e- (bound) + nu_e + nubar_e
    // equation 6.11 solved for nu

    number_t xnum   = 1.10520e8_rt * sf.den * sf.ye / (sf.temp * admath::sqrt(sf.temp));

    // the chemical potential
    number_t nu   = ifermi12(xnum);
    number_t nu2  = nu * nu;
    number_t nu3  = nu2 * nu;

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
    // equation 6.7, 6.13 and 6.14
    if (nu >= -20.0_rt &&  nu <= 10.0_rt) {

        number_t tempi  = 1.0e0_rt / sf.temp;
        number_t zeta   = 1.579e5_rt * sf.zbar * sf.zbar * tempi;

        number_t c00 = 1.0e0_rt / (1.0e0_rt + f1 * nu + f2 * nu2 + f3 * nu3);
        number_t dum = zeta * c00;

        number_t z = 1.0e0_rt / dum;
        number_t dd00 = admath::pow(dum, -2.25_rt);
        number_t dd01 = admath::pow(dum, -4.55_rt);
        c00  = a1 * z + a2 * dd00 + a3 * dd01;

        z = admath::exp(c * nu);
        dd00 = b * z * (1.0e0_rt + d * dum);
        number_t gum = 1.0e0_rt + dd00;

        z   = admath::exp(nu);
        number_t tmp  = 1.0e0_rt / gum;

        number_t bigj = c00 * z * tmp;

        // equation 6.5
        z     = admath::exp(zeta + nu);
        dum   = 1.0e0_rt + z;
        tmp   = 1.0e0_rt/dum;

        sreco   = amrex::Real(nu_constants::tfac6) * 2.649e-18_rt * sf.ye * amrex::Math::powi<13>(sf.zbar) * sf.den * bigj * tmp;
    }
    return sreco;

}

template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
number_t sneut5(const number_t& temp, const amrex::Real den,
                const number_t& abar, const number_t& zbar, amrex::Real& pair,
                amrex::Real& phot, amrex::Real& plas, amrex::Real& brem)
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

    pair = autodiff::val(spair);
    phot = autodiff::val(sphot);
    plas = autodiff::val(splas);
    brem = autodiff::val(sbrem);

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
            amrex::Real& dsnuda, amrex::Real& dsnudz, amrex::Real& pair,
            amrex::Real& phot, amrex::Real& plas, amrex::Real& brem)
{
    // autodiff wrapper
    if constexpr (do_derivatives) {
        using dual_t = autodiff::dual_array<1, 3>;
        dual_t temp_dual = temp;
        dual_t abar_dual = abar;
        dual_t zbar_dual = zbar;
        autodiff::seed_array(temp_dual, abar_dual, zbar_dual);
        dual_t snu_dual = sneut5(temp_dual, den, abar_dual, zbar_dual,
                                 pair, phot, plas, brem);
        snu = autodiff::val(snu_dual);
        const auto& grad = autodiff::derivative(snu_dual);
        dsnudt = grad(1);
        dsnudd = 0.0e0_rt;
        dsnuda = grad(2);
        dsnudz = grad(3);
    } else {
        snu = sneut5(temp, den, abar, zbar, pair, phot, plas, brem);
        dsnudt = 0.0e0_rt;
        dsnudd = 0.0e0_rt;
        dsnuda = 0.0e0_rt;
        dsnudz = 0.0e0_rt;
    }

}

template <int do_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void sneut5(const amrex::Real temp, const amrex::Real den,
            const amrex::Real abar, const amrex::Real zbar,
            amrex::Real& snu, amrex::Real& dsnudt, amrex::Real& dsnudd,
            amrex::Real& dsnuda, amrex::Real& dsnudz)
            {
                amrex::Real pair, phot, plas, brem;
                sneut5<do_derivatives>(temp, den, abar, zbar,
                                       snu, dsnudt, dsnudd, dsnuda, dsnudz,
                                       pair, phot, plas, brem);
            }

#endif
