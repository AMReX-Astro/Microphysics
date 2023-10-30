#ifndef SNEUT5_H
#define SNEUT5_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Math.H>

using namespace amrex;

namespace nu_constants {

    // numerical constants
    constexpr Real fac1   = 5.0e0_rt * M_PI / 3.0e0_rt;
    constexpr Real fac2   = 10.0e0_rt * M_PI;
    constexpr Real fac3   = M_PI / 5.0e0_rt;
    constexpr Real oneth  = 1.0e0_rt/3.0e0_rt;
    constexpr Real twoth  = 2.0e0_rt/3.0e0_rt;
    constexpr Real sixth  = 1.0e0_rt/6.0e0_rt;
    constexpr Real iln10  = 4.342944819032518e-1_rt;

    // theta is sin**2(theta_weinberg) = 0.2319 plus/minus 0.00005 (1996)
    // xnufam is the number of neutrino flavors = 3.02 plus/minus 0.005 (1998)
    // change theta and xnufam if need be, and the changes will automatically
    // propagate through the routine. cv and ca are the vector and axial currents.

    constexpr Real theta  = 0.2319e0_rt;
    constexpr Real xnufam = 3.0e0_rt;
    constexpr Real cv     = 0.5e0_rt + 2.0e0_rt * theta;
    constexpr Real cvp    = 1.0e0_rt - cv;
    constexpr Real ca     = 0.5e0_rt;
    constexpr Real cap    = 1.0e0_rt - ca;
    constexpr Real tfac1  = cv*cv + ca*ca + (xnufam-1.0e0_rt) * (cvp*cvp+cap*cap);
    constexpr Real tfac2  = cv*cv - ca*ca + (xnufam-1.0e0_rt) * (cvp*cvp - cap*cap);
    constexpr Real tfac3  = tfac2/tfac1;
    constexpr Real tfac4  = 0.5e0_rt * tfac1;
    constexpr Real tfac5  = 0.5e0_rt * tfac2;
    constexpr Real tfac6  = cv*cv + 1.5e0_rt*ca*ca + (xnufam - 1.0e0_rt)*(cvp*cvp + 1.5e0_rt*cap*cap);

}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real ifermi12(const Real f)
{

    // this routine applies a rational function expansion to get the inverse
    // fermi-dirac integral of order 1/2 when it is equal to f.
    // maximum error is 4.19e-9_rt.   reference: antia apjs 84,101 1993

    // declare work variables
    Real rn,den,ff;

    // the return value
    Real ifermi12r;

    // load the coefficients of the expansion from Table 8 of Antia

    constexpr Real an{0.5e0_rt};
    constexpr int m1{4};
    constexpr int k1{3};
    constexpr int m2{6};
    constexpr int k2{5};

    const Array1D<Real, 1, 5> a1 = {
        1.999266880833e4_rt,
        5.702479099336e3_rt,
        6.610132843877e2_rt,
        3.818838129486e1_rt,
        1.0e0_rt};

    const Array1D<Real, 1, 4> b1 = {
        1.771804140488e4_rt,
        -2.014785161019e3_rt,
        9.130355392717e1_rt,
        -1.670718177489e0_rt};

    const Array1D<Real, 1, 7> a2 = {
        -1.277060388085e-2_rt,
        7.187946804945e-2_rt,
        -4.262314235106e-1_rt,
        4.997559426872e-1_rt,
        -1.285579118012e0_rt,
        -3.930805454272e-1_rt,
        1.0e0_rt};

    const Array1D<Real, 1, 6> b2 = {
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

        ifermi12r = std::log(f * rn/den);

    } else {

        ff = 1.0e0_rt/std::pow(f, 1.0e0_rt/(1.0e0_rt + an));

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

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real zfermim12(const Real x)
{

    // this routine applies a rational function expansion to get the fermi-dirac
    // integral of order -1/2 evaluated at x. maximum error is 1.23e-12_rt.
    // reference: antia apjs 84,101 1993

    // declare work variables
    Real rn,den,xx;

    // return value
    Real zfermim12r;

    // load the coefficients of the expansion from Table 2 of Antia

    constexpr int m1{7};
    constexpr int k1{7};
    constexpr int m2{11};
    constexpr int k2{11};

    const Array1D<Real, 1, 8> a1 = {
        1.71446374704454e7_rt,
        3.88148302324068e7_rt,
        3.16743385304962e7_rt,
        1.14587609192151e7_rt,
        1.83696370756153e6_rt,
        1.14980998186874e5_rt,
        1.98276889924768e3_rt,
        1.0e0_rt};

    const Array1D<Real, 1, 8> b1 = {
        9.67282587452899e6_rt,
        2.87386436731785e7_rt,
        3.26070130734158e7_rt,
        1.77657027846367e7_rt,
        4.81648022267831e6_rt,
        6.13709569333207e5_rt,
        3.13595854332114e4_rt,
        4.35061725080755e2_rt};

    const Array1D<Real, 1, 12> a2 = {
        -4.46620341924942e-15_rt,
        -1.58654991146236e-12_rt,
        -4.44467627042232e-10_rt,
        -6.84738791621745e-8_rt,
        -6.64932238528105e-6_rt,
        -3.69976170193942e-4_rt,
        -1.12295393687006e-2_rt,
        -1.60926102124442e-1_rt,
        -8.52408612877447e-1_rt,
        -7.45519953763928e-1_rt,
        2.98435207466372e0_rt,
        1.0e0_rt};

    const Array1D<Real, 1, 12> b2 = {
        -2.23310170962369e-15_rt,
        -7.94193282071464e-13_rt,
        -2.22564376956228e-10_rt,
        -3.43299431079845e-8_rt,
        -3.33919612678907e-6_rt,
        -1.86432212187088e-4_rt,
        -5.69764436880529e-3_rt,
        -8.34904593067194e-2_rt,
        -4.78770844009440e-1_rt,
        -4.99759250374148e-1_rt,
        1.86795964993052e0_rt,
        4.16485970495288e-1_rt};

    if (x < 2.0e0_rt) {

       xx = std::exp(x);
       rn = xx + a1(m1);

       for (int i = m1 - 1; i >= 1; --i) {
          rn = rn*xx + a1(i);
       }

       den = b1(k1+1);

       for (int i = k1; i >= 1; --i) {
          den = den*xx + b1(i);
       }

       zfermim12r = xx * rn/den;

    } else {

       xx = 1.0e0_rt/(x*x);
       rn = xx + a2(m2);

       for (int i = m2 - 1; i >= 1; --i) {
          rn = rn*xx + a2(i);
       }

       den = b2(k2+1);

       for (int i = k2; i >= 1; --i) {
          den = den*xx + b2(i);
       }

       zfermim12r = std::sqrt(x)*rn/den;

    }

    return zfermim12r;
}

struct sneutf_t {

    Real den;
    Real temp;
    Real abar;
    Real zbar;

    Real deni;
    Real tempi;
    Real abari;
    Real zbari;

    Real ye;

    Real t9;
    Real xl;
    Real xldt;
    Real xlp5;
    Real xl2;
    Real xl3;
    Real xl4;
    Real xl5;
    Real xl6;
    Real xl7;
    Real xl8;
    Real xl9;
    Real xlmp5;
    Real xlm1;
    Real xlm2;
    Real xlm3;
    Real xlm4;
    Real rm;
    Real rmda;
    Real rmdz;
    Real rmi;

    Real zeta;
    Real zeta2;
    Real zeta3;
    Real zetadt;
    Real zetada;
    Real zetadz;
};


template <int do_derivatives>
AMREX_GPU_HOST_DEVICE inline
sneutf_t get_sneut_factors(Real den, Real temp, Real abar, Real zbar) {

    constexpr Real con1   = 1.0e0_rt/5.9302e0_rt;

    sneutf_t sf;

    sf.den = den;
    sf.temp = temp;
    sf.abar = abar;
    sf.zbar = zbar;

    // to avoid lots of divisions
    sf.deni  = 1.0e0_rt / den;
    sf.tempi = 1.0e0_rt / temp;
    sf.abari = 1.0e0_rt / abar;
    sf.zbari = 1.0e0_rt / zbar;

    // some composition variables
    sf.ye    = zbar * sf.abari;

    // some frequent factors
    sf.t9     = temp * 1.0e-9_rt;
    sf.xl     = sf.t9 * con1;
    sf.xldt   = 1.0e-9_rt * con1;
    sf.xlp5   = std::sqrt(sf.xl);
    sf.xl2    = sf.xl * sf.xl;
    sf.xl3    = sf.xl2 * sf.xl;
    sf.xl4    = sf.xl3 * sf.xl;
    sf.xl5    = sf.xl4 * sf.xl;
    sf.xl6    = sf.xl5 * sf.xl;
    sf.xl7    = sf.xl6 * sf.xl;
    sf.xl8    = sf.xl7 * sf.xl;
    sf.xl9    = sf.xl8 * sf.xl;
    sf.xlmp5  = 1.0e0_rt / sf.xlp5;
    sf.xlm1   = 1.0e0_rt / sf.xl;
    sf.xlm2   = sf.xlm1 * sf.xlm1;
    sf.xlm3   = sf.xlm1 * sf.xlm2;
    sf.xlm4   = sf.xlm1 * sf.xlm3;

    sf.rm     = den * sf.ye;
    sf.rmda   = -sf.rm * sf.abari;
    sf.rmdz   = den * sf.abari;
    sf.rmi    = 1.0e0_rt / sf.rm;

    Real a0 = sf.rm * 1.0e-9_rt;
    Real a1 = std::pow(a0, nu_constants::oneth);
    sf.zeta = a1 * sf.xlm1;

    if constexpr (do_derivatives) {
        Real a2 = nu_constants::oneth * a1 * sf.rmi * sf.xlm1;
        sf.zetadt = -a1 * sf.xlm2 * sf.xldt;
        sf.zetada = a2 * sf.rmda;
        sf.zetadz = a2 * sf.rmdz;
    }

    sf.zeta2 = sf.zeta * sf.zeta;
    sf.zeta3 = sf.zeta2 * sf.zeta;

    return sf;

}



template <int do_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nu_pair(const sneutf_t& sf,
             Real& spair, Real& spairdt, Real& spairda, Real& spairdz) {

    // pair neutrino section
    // for reactions like e+ + e- => nu_e + nubar_e

    // equation 2.8
    Real gl = 1.0e0_rt - 13.04e0_rt * sf.xl2 +133.5e0_rt * sf.xl4 +1534.0e0_rt * sf.xl6 + 918.6e0_rt * sf.xl8;
    Real gldt;
    if constexpr (do_derivatives) {
        gldt = sf.xldt*(-26.08e0_rt * sf.xl + 534.0e0_rt * sf.xl3 + 9204.0e0_rt * sf.xl5 + 7348.8e0_rt * sf.xl7);
    }

    // equation 2.7

    Real a1 = 6.002e19_rt + 2.084e20_rt * sf.zeta + 1.872e21_rt * sf.zeta2;
    Real b1, b2;
    if (sf.t9 < 10.0_rt) {
        b1 = std::exp(-5.5924e0_rt * sf.zeta);
        if constexpr (do_derivatives) {
            b2 = -b1 * 5.5924e0_rt;
        }
    } else {
        b1 = std::exp(-4.9924e0_rt * sf.zeta);
        if constexpr (do_derivatives) {
            b2 = -b1 * 4.9924e0_rt;
        }
    }

    Real a2, c;

    Real xnum = a1 * b1;
    Real xnumdt, xnumda, xnumdz;
    if constexpr (do_derivatives) {
        a2 = 2.084e20_rt + 2.0e0_rt * 1.872e21_rt * sf.zeta;
        c = a2 * b1 + a1 * b2;
        xnumdt = c * sf.zetadt;
        xnumda = c * sf.zetada;
        xnumdz = c * sf.zetadz;
    }

    if (sf.t9 < 10.0_rt) {
        a1 = 9.383e-1_rt * sf.xlm1 - 4.141e-1_rt * sf.xlm2 + 5.829e-2_rt * sf.xlm3;
        if constexpr (do_derivatives) {
            a2 = -9.383e-1_rt * sf.xlm2 + 2.0e0_rt * 4.141e-1_rt * sf.xlm3 - 3.0e0_rt * 5.829e-2_rt * sf.xlm4;
        }
    } else {
        a1 = 1.2383e0_rt * sf.xlm1 - 8.141e-1_rt * sf.xlm2;
        if constexpr (do_derivatives) {
            a2 = -1.2383e0_rt * sf.xlm2 + 2.0e0_rt * 8.141e-1_rt * sf.xlm3;
        }
    }

    Real xden = sf.zeta3 + a1;
    Real xdendt, xdenda, xdendz;
    if constexpr (do_derivatives) {
        b1 = 3.0e0_rt * sf.zeta2;
        xdendt = b1 * sf.zetadt + a2 * sf.xldt;
        xdenda = b1 * sf.zetada;
        xdendz = b1 * sf.zetadz;
    }

    a1 = 1.0e0_rt / xden;
    Real fpair = xnum * a1;
    Real fpairdt, fpairda, fpairdz;
    if constexpr (do_derivatives) {
        fpairdt = (xnumdt - fpair * xdendt) * a1;
        fpairda = (xnumda - fpair * xdenda) * a1;
        fpairdz = (xnumdz - fpair * xdendz) * a1;
    }

    // equation 2.6
    a1 = 10.7480e0_rt * sf.xl2 + 0.3967e0_rt * sf.xlp5 + 1.005e0_rt;
    xnum = 1.0e0_rt / a1;
    if constexpr (do_derivatives) {
        a2 = sf.xldt * (2.0e0_rt * 10.7480e0_rt * sf.xl + 0.5e0_rt * 0.3967e0_rt * sf.xlmp5);
        xnumdt = -xnum * xnum * a2;
    }

    a1 = 7.692e7_rt * sf.xl3 + 9.715e6_rt * sf.xlp5;
    c = 1.0e0_rt / a1;
    b1 = 1.0e0_rt + sf.rm * c;

    xden = std::pow(b1, -0.3e0_rt);
    if constexpr (do_derivatives) {
        a2 = sf.xldt * (3.0e0_rt * 7.692e7_rt * sf.xl2 + 0.5e0_rt * 9.715e6_rt * sf.xlmp5);
        Real d = -0.3e0_rt * xden / b1;
        xdendt = -d * sf.rm * c * c * a2;
        xdenda = d * sf.rmda * c;
        xdendz = d * sf.rmdz * c;
    }

    Real qpair = xnum * xden;
    Real qpairdt, qpairda, qpairdz;
    if constexpr (do_derivatives) {
        qpairdt = xnumdt * xden + xnum * xdendt;
        qpairda = xnum * xdenda;
        qpairdz = xnum * xdendz;
    }

    // equation 2.5
    a1 = std::exp(-2.0e0_rt * sf.xlm1);

    spair = a1 * fpair;
    if constexpr (do_derivatives) {
        a2 = a1 * 2.0e0_rt * sf.xlm2 * sf.xldt;
        spairdt = a2 * fpair + a1 * fpairdt;
        spairda = a1 * fpairda;
        spairdz = a1 * fpairdz;
    }

    a1 = spair;
    spair = gl*a1;
    if constexpr (do_derivatives) {
        spairdt = gl * spairdt + gldt * a1;
        spairda = gl * spairda;
        spairdz = gl * spairdz;
    }

    a1 = nu_constants::tfac4 * (1.0e0_rt + nu_constants::tfac3 * qpair);

    Real a3 = spair;
    spair = a1 * a3;
    if constexpr (do_derivatives) {
        a2 = nu_constants::tfac4 * nu_constants::tfac3;
        spairdt = a1 * spairdt + a2 * qpairdt * a3;
        spairda = a1 * spairda + a2 * qpairda * a3;
        spairdz = a1 * spairdz + a2 * qpairdz * a3;
    }

}

template <int do_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nu_plasma(const sneutf_t& sf,
               Real& splas, Real& splasdt, Real& splasda, Real& splasdz) {

    // plasma neutrino section
    // for collective reactions like gamma_plasmon => nu_e + nubar_e
    // equation 4.6

    Real a1 = 1.019e-6_rt * sf.rm;
    Real a2 = std::pow(a1, nu_constants::twoth);

    Real b1 = std::sqrt(1.0e0_rt + a2);

    Real c00 = 1.0e0_rt / (sf.temp * sf.temp * b1);

    Real gl2   = 1.1095e11_rt * sf.rm * c00;
    Real gl2dt, gl2da, gl2dz;
    if constexpr (do_derivatives) {
        Real a3 = nu_constants::twoth * a2 / a1;
        Real b2 = 1.0e0_rt / b1;
        gl2dt = -2.0e0_rt * gl2 * sf.tempi;
        Real d = sf.rm * c00 * b2 * 0.5e0_rt * b2 * a3 * 1.019e-6_rt;
        gl2da = 1.1095e11_rt * (sf.rmda * c00  - d * sf.rmda);
        gl2dz = 1.1095e11_rt * (sf.rmdz * c00  - d * sf.rmdz);
    }

    Real gl = std::sqrt(gl2);
    Real gl12 = std::sqrt(gl);
    Real gl32 = gl * gl12;
    Real gl72 = gl2 * gl32;
    Real gl6 = gl2 * gl2 * gl2;

    // equation 4.7
    Real ft = 2.4e0_rt + 0.6e0_rt * gl12 + 0.51e0_rt * gl + 1.25e0_rt * gl32;
    Real gum = 1.0e0_rt / gl2;
    Real ftdt, ftda, ftdz;
    if constexpr (do_derivatives) {
        a1 = (0.25e0_rt * 0.6e0_rt * gl12 +
              0.5e0_rt * 0.51e0_rt * gl +
              0.75e0_rt * 1.25e0_rt * gl32) * gum;
        ftdt = a1 * gl2dt;
        ftda = a1 * gl2da;
        ftdz = a1 * gl2dz;
    }

    // equation 4.8
    a1 = 8.6e0_rt * gl2 + 1.35e0_rt * gl72;

    b1 = 225.0e0_rt - 17.0e0_rt * gl + gl2;

    Real c = 1.0e0_rt / b1;
    Real fl = a1 * c;
    Real fldt, flda, fldz;
    if constexpr (do_derivatives) {
        a2 = 8.6e0_rt + 1.75e0_rt * 1.35e0_rt * gl72 * gum;
        Real b2 = -0.5e0_rt * 17.0e0_rt * gl * gum + 1.0e0_rt;
        Real d = (a2 - fl * b2) * c;
        fldt = d * gl2dt;
        flda = d * gl2da;
        fldz = d * gl2dz;
    }

    // equation 4.9 and 4.10
    Real cc = std::log10(2.0e0_rt * sf.rm);
    Real xlnt = std::log10(sf.temp);

    Real xnum = nu_constants::sixth * (17.5e0_rt + cc - 3.0e0_rt * xlnt);
    Real xden = nu_constants::sixth * (-24.5e0_rt + cc + 3.0e0_rt * xlnt);
    Real xnumdt, xnumda, xnumdz;
    Real xdendt, xdenda, xdendz;
    if constexpr (do_derivatives) {
        xnumdt = -nu_constants::iln10 * 0.5e0_rt * sf.tempi;
        a2 = nu_constants::iln10 * nu_constants::sixth * sf.rmi;
        xnumda = a2 * sf.rmda;
        xnumdz = a2 * sf.rmdz;

        xdendt = nu_constants::iln10 * 0.5e0_rt * sf.tempi;
        xdenda = a2 * sf.rmda;
        xdendz = a2 * sf.rmdz;
    }

    // equation 4.11
    Real fxy, fxydt, fxyda, fxydz;
    Real dumdt, dumda, dumdz;
    if (std::abs(xnum) > 0.7e0_rt || xden < 0.0e0_rt) {
        fxy   = 1.0e0_rt;
        fxydt = 0.0e0_rt;
        fxydz = 0.0e0_rt;
        fxyda = 0.0e0_rt;

    } else {

        a1 = 0.39e0_rt - 1.25e0_rt * xnum - 0.35e0_rt * std::sin(4.5e0_rt * xnum);

        b1 = 0.3e0_rt * std::exp(-std::pow(4.5e0_rt*xnum + 0.9e0_rt, 2.0_rt));

        c = amrex::min(0.0e0_rt, xden - 1.6e0_rt + 1.25e0_rt * xnum);

        Real b2;
        if constexpr (do_derivatives) {
            a2 = -1.25e0_rt - 4.5e0_rt * 0.35e0_rt * std::cos(4.5e0_rt * xnum);
            b2 = -b1 * 2.0e0_rt * (4.5e0_rt * xnum + 0.9e0_rt) * 4.5e0_rt;
            if (c == 0.0_rt) {
                dumdt = 0.0e0_rt;
                dumda = 0.0e0_rt;
                dumdz = 0.0e0_rt;
            } else {
                dumdt = xdendt + 1.25e0_rt * xnumdt;
                dumda = xdenda + 1.25e0_rt * xnumda;
                dumdz = xdendz + 1.25e0_rt * xnumdz;
            }
        }

        Real d = 0.57e0_rt - 0.25e0_rt * xnum;
        Real a3 = c / d;
        c00 = std::exp(-a3 * a3);

        fxy = 1.05e0_rt + (a1 - b1) * c00;
        if constexpr (do_derivatives) {
            Real f1 = -c00 * 2.0e0_rt * a3 / d;
            Real c01 = f1 * (dumdt + a3 * 0.25e0_rt * xnumdt);
            Real c03 = f1 * (dumda + a3 * 0.25e0_rt * xnumda);
            Real c04 = f1 * (dumdz + a3 * 0.25e0_rt * xnumdz);

            fxydt = (a2 * xnumdt -  b2 * xnumdt) * c00 + (a1 - b1) * c01;
            fxyda = (a2 * xnumda -  b2 * xnumda) * c00 + (a1 - b1) * c03;
            fxydz = (a2 * xnumdz -  b2 * xnumdz) * c00 + (a1 - b1) * c04;
        }
    }

    // equation 4.1 and 4.5
    splas = (ft + fl) * fxy;
    if constexpr (do_derivatives) {
        splasdt = (ftdt + fldt) * fxy + (ft + fl) * fxydt;
        splasda = (ftda + flda) * fxy + (ft + fl) * fxyda;
        splasdz = (ftdz + fldz) * fxy + (ft + fl) * fxydz;
    }

    a2 = std::exp(-gl);

    a1 = splas;
    splas = a2*a1;
    if constexpr (do_derivatives) {
        Real a3 = -0.5e0_rt * a2 * gl * gum;
        splasdt = a2 * splasdt + a3 * gl2dt * a1;
        splasda = a2 * splasda + a3 * gl2da * a1;
        splasdz = a2 * splasdz + a3 * gl2dz * a1;
    }

    a2      = gl6;

    a1      = splas;
    splas   = a2 * a1;
    if constexpr (do_derivatives) {
        Real a3 = 3.0e0_rt * gl6 * gum;
        splasdt = a2 * splasdt + a3 * gl2dt * a1;
        splasda = a2 * splasda + a3 * gl2da * a1;
        splasdz = a2 * splasdz + a3 * gl2dz * a1;
    }

    a2 = 0.93153e0_rt * 3.0e21_rt * sf.xl9;

    a1 = splas;
    splas = a2 * a1;
    if constexpr (do_derivatives) {
        Real a3 = 0.93153e0_rt * 3.0e21_rt * 9.0e0_rt * sf.xl8 * sf.xldt;
        splasdt = a2 * splasdt + a3 * a1;
        splasda = a2 * splasda;
        splasdz = a2 * splasdz;
    }
}

template <int do_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nu_photo(const sneutf_t& sf,
              Real& sphot, Real& sphotdt, Real& sphotda, Real& sphotdz) {

    // photoneutrino process section
    // for reactions like e- + gamma => e- + nu_e + nubar_e
    //                    e+ + gamma => e+ + nu_e + nubar_e
    // equation 3.8 for tau, equation 3.6 for cc,
    // and table 2 written out for speed

    Real tau;
    Real cc;
    Real c00, c01, c02, c03, c04, c05, c06;
    Real c10, c11, c12, c13, c14, c15, c16;
    Real c20, c21, c22, c23, c24, c25, c26;
    Real dd01, dd02, dd03, dd04, dd05;
    Real dd11, dd12, dd13, dd14, dd15;
    Real dd21, dd22, dd23, dd24, dd25;

    if (sf.temp < 1.0e8_rt) {

        // note: we already bailed above for T < 1.e7, so this is really 1.e7 <= T < 1.e8

        tau  =  std::log10(sf.temp * 1.0e-7_rt);
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

        tau  =  std::log10(sf.temp * 1.0e-8_rt);
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

        tau  =  std::log10(sf.t9);
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

    Real taudt = nu_constants::iln10 * sf.tempi;

    // equation 3.7

    const auto [cos1, sin1] = amrex::Math::sincos(nu_constants::fac1 * tau);

    // double, triple, etc. angle formulas
    // sin/cos (2 fac1 tau)
    const Real sin2 = 2.0_rt * sin1 * cos1;
    const Real cos2 = 2.0_rt * cos1 * cos1 - 1.0_rt;

    // sin/cos (3 fac1 tau)
    const Real sin3 = sin1 * (3.0_rt - 4.0_rt * sin1 * sin1);
    const Real cos3 = cos1 * (4.0_rt * cos1 * cos1 - 3.0_rt);

    // sin/cos (4 fac1 tau) -- use double angle on sin2/cos2
    const Real sin4 = 2.0_rt * sin2 * cos2;
    const Real cos4 = 2.0_rt * cos2 * cos2 - 1.0_rt;

    // sin/cos (5 fac1 tau)
    const Real sin5 = sin1 * (5.0_rt - sin1 * sin1 * (20.0_rt - 16.0_rt * sin1 * sin1));
    const Real cos5 = cos1 * (cos1 * cos1 * (16.0_rt * cos1 * cos1 - 20.0_rt) + 5.0_rt);

    const auto [xast, last] = amrex::Math::sincos(nu_constants::fac2 * tau);

    Real a0 = 0.5e0_rt * c00
        + c01 * cos1 + dd01 * sin1 + c02 * cos2 + dd02 * sin2
        + c03 * cos3 + dd03 * sin3 + c04 * cos4 + dd04 * sin4
        + c05 * cos5 + dd05 * sin5 + 0.5e0_rt * c06 * last;
    Real a1 = 0.5e0_rt * c10
        + c11 * cos1 + dd11 * sin1 + c12 * cos2 + dd12 * sin2
        + c13 * cos3 + dd13 * sin3 + c14 * cos4 + dd14 * sin4
        + c15 * cos5 + dd15 * sin5 + 0.5e0_rt * c16 * last;
    Real a2 = 0.5e0_rt * c20
        + c21 * cos1 + dd21 * sin1 + c22 * cos2 + dd22 * sin2
        + c23 * cos3 + dd23 * sin3 + c24 * cos4 + dd24 * sin4
        + c25 * cos5 + dd25 * sin5 + 0.5e0_rt * c26 * last;

    Real f0, f1, f2;
    if constexpr (do_derivatives) {
        f0 =  taudt * nu_constants::fac1 *
            (-c01 * sin1 + dd01 * cos1
             - c02 * sin2 * 2.0e0_rt + dd02 * cos2 * 2.0e0_rt
             - c03 * sin3 * 3.0e0_rt + dd03 * cos3 * 3.0e0_rt
             - c04 * sin4 * 4.0e0_rt + dd04 * cos4 * 4.0e0_rt
             - c05 * sin5 * 5.0e0_rt + dd05 * cos5 * 5.0e0_rt)
            - 0.5e0_rt * c06 * xast * nu_constants::fac2 * taudt;

        f1 = taudt * nu_constants::fac1 *
            (-c11 * sin1 + dd11 * cos1
             - c12 * sin2 * 2.0e0_rt + dd12 * cos2 * 2.0e0_rt
             - c13 * sin3 * 3.0e0_rt + dd13 * cos3 * 3.0e0_rt
             - c14 * sin4 * 4.0e0_rt + dd14 * cos4 * 4.0e0_rt
             - c15 * sin5 * 5.0e0_rt + dd15*cos5*5.0e0_rt)
            - 0.5e0_rt*c16*xast*nu_constants::fac2*taudt;

        f2 = taudt * nu_constants::fac1 *
            (-c21 * sin1 + dd21 * cos1
             - c22 * sin2 * 2.0e0_rt + dd22 * cos2 * 2.0e0_rt
             - c23 * sin3 * 3.0e0_rt + dd23 * cos3 * 3.0e0_rt
             - c24 * sin4 * 4.0e0_rt + dd24 * cos4 * 4.0e0_rt
             - c25 * sin5 * 5.0e0_rt + dd25 * cos5 * 5.0e0_rt)
            - 0.5e0_rt * c26 * xast * nu_constants::fac2 * taudt;
    }

    // equation 3.4
    Real dum = a0 + a1 * sf.zeta + a2 * sf.zeta2;
    Real dumdt, dumda, dumdz;
    if constexpr (do_derivatives) {
        dumdt = f0 + f1 * sf.zeta + a1 * sf.zetadt + f2 * sf.zeta2 + 2.0e0_rt * a2 * sf.zeta * sf.zetadt;
        dumda = a1 * sf.zetada + 2.0e0_rt * a2 * sf.zeta * sf.zetada;
        dumdz = a1 * sf.zetadz + 2.0e0_rt * a2 * sf.zeta * sf.zetadz;
    }

    Real z = std::exp(-cc * sf.zeta);

    Real xnum = dum * z;
    Real xnumdt, xnumda, xnumdz;
    if constexpr (do_derivatives) {
        xnumdt = dumdt * z - dum * z * cc * sf.zetadt;
        xnumda = dumda * z - dum * z * cc * sf.zetada;
        xnumdz = dumdz * z - dum * z * cc * sf.zetadz;
    }

    Real xden = sf.zeta3 + 6.290e-3_rt * sf.xlm1 + 7.483e-3_rt * sf.xlm2 + 3.061e-4_rt * sf.xlm3;

    dum  = 3.0e0_rt * sf.zeta2;
    Real xdendt, xdenda, xdendz;
    if constexpr (do_derivatives) {
        xdendt = dum * sf.zetadt - sf.xldt * (6.290e-3_rt * sf.xlm2
                                           + 2.0e0_rt * 7.483e-3_rt * sf.xlm3
                                           + 3.0e0_rt * 3.061e-4_rt * sf.xlm4);
        xdenda = dum * sf.zetada;
        xdendz = dum * sf.zetadz;
    }
    dum = 1.0e0_rt / xden;
    Real fphot = xnum * dum;
    Real fphotdt, fphotda, fphotdz;
    if constexpr (do_derivatives) {
        fphotdt = (xnumdt - fphot * xdendt) * dum;
        fphotda = (xnumda - fphot * xdenda) * dum;
        fphotdz = (xnumdz - fphot * xdendz) * dum;
    }

    // equation 3.3
    a0 = 1.0e0_rt + 2.045e0_rt * sf.xl;
    xnum = 0.666e0_rt * std::pow(a0, -2.066e0_rt);
    if constexpr (do_derivatives) {
        xnumdt = -2.066e0_rt * xnum / a0 * 2.045e0_rt * sf.xldt;
    }

    dum = 1.875e8_rt * sf.xl + 1.653e8_rt * sf.xl2 + 8.499e8_rt * sf.xl3 - 1.604e8_rt * sf.xl4;
    if constexpr (do_derivatives) {
        dumdt = sf.xldt * (1.875e8_rt
                           + 2.0e0_rt * 1.653e8_rt * sf.xl
                           + 3.0e0_rt * 8.499e8_rt * sf.xl2
                           - 4.0e0_rt * 1.604e8_rt * sf.xl3);
    }

    z  = 1.0e0_rt / dum;
    xden = 1.0e0_rt + sf.rm * z;
    if constexpr (do_derivatives) {
        xdendt = -sf.rm * z * z * dumdt;
        xdenda = sf.rmda * z;
        xdendz = sf.rmdz * z;
    }

    z  = 1.0e0_rt / xden;
    Real qphot = xnum * z;
    Real qphotdt;
    if constexpr (do_derivatives) {
        qphotdt = (xnumdt - qphot * xdendt) * z;
    }
    dum = -qphot * z;
    Real qphotda, qphotdz;
    if constexpr (do_derivatives) {
        qphotda = dum * xdenda;
        qphotdz = dum * xdendz;
    }

    // equation 3.2
    sphot = sf.xl5 * fphot;
    if constexpr (do_derivatives) {
        sphotdt = 5.0e0_rt * sf.xl4 * sf.xldt * fphot + sf.xl5 * fphotdt;
        sphotda = sf.xl5 * fphotda;
        sphotdz = sf.xl5 * fphotdz;
    }

    a1  = sphot;
    sphot = sf.rm * a1;
    if constexpr (do_derivatives) {
        sphotdt = sf.rm * sphotdt;
        sphotda = sf.rm * sphotda + sf.rmda * a1;
        sphotdz = sf.rm * sphotdz + sf.rmdz * a1;
    }

    a1 = nu_constants::tfac4 * (1.0e0_rt - nu_constants::tfac3 * qphot);
    Real a3 = sphot;
    sphot = a1 * a3;
    if constexpr (do_derivatives) {
        a2 = -nu_constants::tfac4 * nu_constants::tfac3;
        sphotdt = a1 * sphotdt + a2 * qphotdt * a3;
        sphotda = a1 * sphotda + a2 * qphotda * a3;
        sphotdz = a1 * sphotdz + a2 * qphotdz * a3;
    }

    if (sphot <= 0.0_rt) {
       sphot   = 0.0e0_rt;
       if constexpr (do_derivatives) {
           sphotdt = 0.0e0_rt;
           sphotda = 0.0e0_rt;
           sphotdz = 0.0e0_rt;
       }
    }
}

template <int do_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nu_brem(const sneutf_t& sf,
             Real& sbrem, Real& sbremdt, Real& sbremda, Real& sbremdz) {

    // bremsstrahlung neutrino section
    // for reactions like e- + (z,a) => e- + (z,a) + nu + nubar
    //                    n  + n     => n + n + nu + nubar
    //                    n  + p     => n + p + nu + nubar
    // equation 4.3

    Real den6   = sf.den * 1.0e-6_rt;
    Real t8     = sf.temp * 1.0e-8_rt;
    Real t812   = std::sqrt(t8);
    Real t832   = t8 * t812;
    Real t82    = t8*t8;
    Real t83    = t82*t8;
    Real t85    = t82*t83;
    Real t86    = t85*t8;
    Real t8m1   = 1.0e0_rt/t8;
    Real t8m2   = t8m1*t8m1;
    Real t8m3   = t8m2*t8m1;
    Real t8m5   = t8m3*t8m2;
    Real t8m6   = t8m5*t8m1;

    Real tfermi = 5.9302e9_rt * (std::sqrt(1.0e0_rt + 1.018e0_rt * std::pow(den6 * sf.ye, nu_constants::twoth)) - 1.0e0_rt);

    // "weak" degenerate electrons only
    if (sf.temp > 0.3e0_rt * tfermi) {

        // equation 5.3
        Real dum = 7.05e6_rt * t832 + 5.12e4_rt * t83;
        Real dumdt;
        if constexpr (do_derivatives) {
            dumdt = (1.5e0_rt * 7.05e6_rt * t812 + 3.0e0_rt * 5.12e4_rt * t82) * 1.0e-8_rt;
        }

        Real z = 1.0e0_rt / dum;
        Real eta = sf.rm * z;
        Real etadt, etada, etadz;
        if constexpr (do_derivatives) {
           etadt = -sf.rm*z*z*dumdt;
           etada = sf.rmda*z;
           etadz = sf.rmdz*z;
       }

        Real etam1 = 1.0e0_rt/eta;
        Real etam2 = etam1 * etam1;
        Real etam3 = etam2 * etam1;

        // equation 5.2
        Real a0 = 23.5e0_rt + 6.83e4_rt * t8m2 + 7.81e8_rt * t8m5;
        Real f0 = (-2.0e0_rt * 6.83e4_rt * t8m3 - 5.0e0_rt * 7.81e8_rt * t8m6) * 1.0e-8_rt;
        Real xnum = 1.0e0_rt / a0;

        dum = 1.0e0_rt + 1.47e0_rt * etam1 + 3.29e-2_rt * etam2;
        Real dumda, dumdz;
        if constexpr (do_derivatives) {
            z = -1.47e0_rt * etam2 - 2.0e0_rt * 3.29e-2_rt * etam3;
            dumdt = z*etadt;
            dumda = z*etada;
            dumdz = z*etadz;
        }

        Real c00 = 1.26e0_rt * (1.0e0_rt + etam1);
        Real c01, c02, c03, c04;
        if constexpr (do_derivatives) {
            z     = -1.26e0_rt*etam2;
            c01   = z*etadt;
            c03   = z*etada;
            c04   = z*etadz;
        }

        z = 1.0e0_rt/dum;
        Real xden = c00 * z;
        Real xdendt, xdenda, xdendz;
        if constexpr (do_derivatives) {
            xdendt = (c01 - xden * dumdt) * z;
            xdenda = (c03 - xden * dumda) * z;
            xdendz = (c04 - xden * dumdz) * z;
        }

        Real fbrem = xnum + xden;
        Real fbremdt, fbremda, fbremdz;
        if constexpr (do_derivatives) {
            fbremdt = -xnum*xnum*f0 + xdendt;
            fbremda = xdenda;
            fbremdz = xdendz;
        }

        // equation 5.9
        a0 = 230.0e0_rt + 6.7e5_rt * t8m2 + 7.66e9_rt * t8m5;
        f0 = (-2.0e0_rt * 6.7e5_rt * t8m3 - 5.0e0_rt * 7.66e9_rt * t8m6) * 1.0e-8_rt;

        z = 1.0e0_rt + sf.rm * 1.0e-9_rt;
        dum = a0 * z;
        if constexpr (do_derivatives) {
            dumdt = f0 * z;
            z = a0 * 1.0e-9_rt;
            dumda = z * sf.rmda;
            dumdz = z * sf.rmdz;
        }

        xnum = 1.0e0_rt / dum;
        Real xnumdt, xnumda, xnumdz;
        if constexpr (do_derivatives) {
            z = -xnum * xnum;
            xnumdt = z * dumdt;
            xnumda = z * dumda;
            xnumdz = z * dumdz;
        }

        c00 = 7.75e5_rt * t832 + 247.0e0_rt * std::pow(t8, 3.85e0_rt);
        c01 = 4.07e0_rt + 0.0240e0_rt * std::pow(t8, 1.4e0_rt);
        c02 = 4.59e-5_rt * std::pow(t8, -0.110e0_rt);

        Real dd00, dd01, dd02;
        if constexpr (do_derivatives) {
            dd00  = (1.5e0_rt * 7.75e5_rt * t812 + 3.85e0_rt * 247.0e0_rt *
                     std::pow(t8, 2.85e0_rt)) * 1.0e-8_rt;

            dd01  = 1.4e0_rt * 0.0240e0_rt * std::pow(t8, 0.4e0_rt) * 1.0e-8_rt;
            dd02  = -0.11e0_rt * 4.59e-5_rt * std::pow(t8, -1.11e0_rt) * 1.0e-8_rt;
        }

        z = std::pow(sf.den, 0.656e0_rt);
        dum = c00 * sf.rmi + c01 + c02 * z;
        if constexpr (do_derivatives) {
            dumdt = dd00 * sf.rmi + dd01 + dd02 * z;
            z     = -c00 * sf.rmi * sf.rmi;
            dumda = z * sf.rmda;
            dumdz = z * sf.rmdz;
        }

        xden  = 1.0e0_rt / dum;
        if constexpr (do_derivatives) {
            z = -xden * xden;
            xdendt = z * dumdt;
            xdenda = z * dumda;
            xdendz = z * dumdz;
        }

        Real gbrem = xnum + xden;
        Real gbremdt, gbremda, gbremdz;
        if constexpr (do_derivatives) {
            gbremdt = xnumdt + xdendt;
            gbremda = xnumda + xdenda;
            gbremdz = xnumdz + xdendz;
        }

        // equation 5.1
        dum = 0.5738e0_rt * sf.zbar * sf.ye * t86 * sf.den;
        if constexpr (do_derivatives) {
            dumdt = 0.5738e0_rt * sf.zbar * sf.ye * 6.0e0_rt * t85 * sf.den * 1.0e-8_rt;
            dumda = -dum * sf.abari;
            dumdz = 0.5738e0_rt * 2.0e0_rt * sf.ye * t86 * sf.den;
        }

        z = nu_constants::tfac4 * fbrem - nu_constants::tfac5 * gbrem;
        sbrem = dum * z;
        if constexpr (do_derivatives) {
            sbremdt = dumdt * z + dum * (nu_constants::tfac4 * fbremdt - nu_constants::tfac5 * gbremdt);
            sbremda = dumda * z + dum * (nu_constants::tfac4 * fbremda - nu_constants::tfac5 * gbremda);
            sbremdz = dumdz * z + dum * (nu_constants::tfac4 * fbremdz - nu_constants::tfac5 * gbremdz);
        }

        // liquid metal with c12 parameters (not too different for other elements)
        // equation 5.18 and 5.16

    } else {

        Real u = nu_constants::fac3 * (std::log10(sf.den) - 3.0e0_rt);
        //a0    = iln10*fac3*sf.deni;

        // compute the expensive trig functions of equation 5.21 only once
        Real cos1 = std::cos(u);
        Real cos2 = std::cos(2.0e0_rt*u);
        Real cos3 = std::cos(3.0e0_rt*u);
        Real cos4 = std::cos(4.0e0_rt*u);
        Real cos5 = std::cos(5.0e0_rt*u);

        Real sin1 = std::sin(u);
        Real sin2 = std::sin(2.0e0_rt*u);
        Real sin3 = std::sin(3.0e0_rt*u);
        Real sin4 = std::sin(4.0e0_rt*u);
        //sin5 = std::sin(5.0e0_rt*u);

        // equation 5.21
        Real fb = 0.5e0_rt * 0.17946e0_rt + 0.00945e0_rt * u + 0.34529e0_rt
            - 0.05821e0_rt * cos1 - 0.04969e0_rt * sin1
            - 0.01089e0_rt * cos2 - 0.01584e0_rt * sin2
            - 0.01147e0_rt * cos3 - 0.00504e0_rt * sin3
            - 0.00656e0_rt * cos4 - 0.00281e0_rt * sin4
            - 0.00519e0_rt * cos5;

        // c00 =  a0*(0.00945e0_rt
        //      + 0.05821e0_rt*sin1       - 0.04969e0_rt*cos1
        //      + 0.01089e0_rt*sin2*2.0e0_rt - 0.01584e0_rt*cos2*2.0e0_rt
        //      + 0.01147e0_rt*sin3*3.0e0_rt - 0.00504e0_rt*cos3*3.0e0_rt
        //      + 0.00656e0_rt*sin4*4.0e0_rt - 0.00281e0_rt*cos4*4.0e0_rt
        //      + 0.00519e0_rt*sin5*5.0e0_rt);

        // equation 5.22
        Real ft = 0.5e0_rt * 0.06781e0_rt - 0.02342e0_rt * u + 0.24819e0_rt
            - 0.00944e0_rt * cos1 - 0.02213e0_rt * sin1
            - 0.01289e0_rt * cos2 - 0.01136e0_rt * sin2
            - 0.00589e0_rt * cos3 - 0.00467e0_rt * sin3
            - 0.00404e0_rt * cos4 - 0.00131e0_rt * sin4
            - 0.00330e0_rt * cos5;

        // c01 = a0*(-0.02342e0_rt
        //      + 0.00944e0_rt*sin1       - 0.02213e0_rt*cos1
        //      + 0.01289e0_rt*sin2*2.0e0_rt - 0.01136e0_rt*cos2*2.0e0_rt
        //      + 0.00589e0_rt*sin3*3.0e0_rt - 0.00467e0_rt*cos3*3.0e0_rt
        //      + 0.00404e0_rt*sin4*4.0e0_rt - 0.00131e0_rt*cos4*4.0e0_rt
        //      + 0.00330e0_rt*sin5*5.0e0_rt);

        // equation 5.23
        Real gb = 0.5e0_rt * 0.00766e0_rt - 0.01259e0_rt * u + 0.07917e0_rt
            - 0.00710e0_rt * cos1 + 0.02300e0_rt * sin1
            - 0.00028e0_rt * cos2 - 0.01078e0_rt * sin2
            + 0.00232e0_rt * cos3 + 0.00118e0_rt * sin3
            + 0.00044e0_rt * cos4 - 0.00089e0_rt * sin4
            + 0.00158e0_rt * cos5;

        // c02 = a0*(-0.01259e0_rt
        //      + 0.00710e0_rt*sin1       + 0.02300e0_rt*cos1
        //      + 0.00028e0_rt*sin2*2.0e0_rt - 0.01078e0_rt*cos2*2.0e0_rt
        //      - 0.00232e0_rt*sin3*3.0e0_rt + 0.00118e0_rt*cos3*3.0e0_rt
        //      - 0.00044e0_rt*sin4*4.0e0_rt - 0.00089e0_rt*cos4*4.0e0_rt
        //      - 0.00158e0_rt*sin5*5.0e0_rt);

        // equation 5.24
        Real gt = -0.5e0_rt * 0.00769e0_rt  - 0.00829e0_rt * u + 0.05211e0_rt
            + 0.00356e0_rt * cos1 + 0.01052e0_rt * sin1
            - 0.00184e0_rt * cos2 - 0.00354e0_rt * sin2
            + 0.00146e0_rt * cos3 - 0.00014e0_rt * sin3
            + 0.00031e0_rt * cos4 - 0.00018e0_rt * sin4
            + 0.00069e0_rt * cos5;

        // c03 = a0*(-0.00829e0_rt
        //      - 0.00356e0_rt*sin1       + 0.01052e0_rt*cos1
        //      + 0.00184e0_rt*sin2*2.0e0_rt - 0.00354e0_rt*cos2*2.0e0_rt
        //      - 0.00146e0_rt*sin3*3.0e0_rt - 0.00014e0_rt*cos3*3.0e0_rt
        //      - 0.00031e0_rt*sin4*4.0e0_rt - 0.00018e0_rt*cos4*4.0e0_rt
        //      - 0.00069e0_rt*sin5*5.0e0_rt);

        Real dum = 2.275e-1_rt * sf.zbar * sf.zbar * t8m1 * std::pow(den6*sf.abari, nu_constants::oneth);
        Real dumdt, dumda, dumdz;
        if constexpr (do_derivatives) {
            dumdt = -dum*sf.tempi;
            dumda = -nu_constants::oneth*dum*sf.abari;
            dumdz = 2.0e0_rt*dum*sf.zbari;
        }

        Real gm1 = 1.0e0_rt / dum;
        Real gm2 = gm1 * gm1;
        Real gm13 = std::pow(gm1, nu_constants::oneth);
        Real gm23 = gm13 * gm13;
        Real gm43 = gm13 * gm1;
        Real gm53 = gm23 * gm1;

        // equation 5.25 and 5.26
        Real v = -0.05483e0_rt - 0.01946e0_rt * gm13 + 1.86310e0_rt * gm23 - 0.78873e0_rt * gm1;
        Real a0 = nu_constants::oneth * 0.01946e0_rt * gm43 - nu_constants::twoth * 1.86310e0_rt * gm53 + 0.78873e0_rt * gm2;

        Real w = -0.06711e0_rt + 0.06859e0_rt * gm13 + 1.74360e0_rt * gm23 - 0.74498e0_rt * gm1;
        Real a1 = -nu_constants::oneth*0.06859e0_rt * gm43 - nu_constants::twoth * 1.74360e0_rt * gm53 + 0.74498e0_rt * gm2;

        // equation 5.19 and 5.20
        Real fliq = v * fb + (1.0e0_rt - v) * ft;
        Real fliqdt, fliqda, fliqdz;
        if constexpr (do_derivatives) {
           fliqdt = a0 * dumdt * (fb - ft);
           fliqda = a0 * dumda * (fb - ft);
           fliqdz = a0 * dumdz * (fb - ft);
        }

        Real gliq = w * gb + (1.0e0_rt - w) * gt;
        Real gliqdt, gliqda, gliqdz;
        if constexpr (do_derivatives) {
            gliqdt = a1 * dumdt*(gb - gt);
            gliqda = a1 * dumda*(gb - gt);
            gliqdz = a1 * dumdz*(gb - gt);
        }

        // equation 5.17
        dum = 0.5738e0_rt * sf.zbar * sf.ye * t86 * sf.den;
        if constexpr (do_derivatives) {
            dumdt = 0.5738e0_rt * sf.zbar * sf.ye * 6.0e0_rt * t85 * sf.den * 1.0e-8_rt;
            dumda = -dum * sf.abari;
            dumdz = 0.5738e0_rt * 2.0e0_rt * sf.ye * t86 * sf.den;
        }

        Real z  = nu_constants::tfac4*fliq - nu_constants::tfac5*gliq;
        sbrem = dum * z;
        if constexpr (do_derivatives) {
            sbremdt = dumdt*z + dum*(nu_constants::tfac4*fliqdt - nu_constants::tfac5*gliqdt);
            sbremda = dumda*z + dum*(nu_constants::tfac4*fliqda - nu_constants::tfac5*gliqda);
            sbremdz = dumdz*z + dum*(nu_constants::tfac4*fliqdz - nu_constants::tfac5*gliqdz);
        }
    }
}

template <int do_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nu_recomb(const sneutf_t& sf,
               Real& sreco, Real& srecodt, Real& srecoda, Real& srecodz) {


    // recombination neutrino section
    // for reactions like e- (continuum) => e- (bound) + nu_e + nubar_e
    // equation 6.11 solved for nu

    Real xnum   = 1.10520e8_rt * sf.den * sf.ye / (sf.temp * std::sqrt(sf.temp));
    Real xnumdt{0.0};
    Real xnumda{0.0};
    Real xnumdz{0.0};

    if constexpr (do_derivatives) {
        xnumdt = -1.50e0_rt*xnum*sf.tempi;
        xnumda = -xnum*sf.abari;
        xnumdz = xnum*sf.zbari;
    }

    // the chemical potential
    Real nu   = ifermi12(xnum);

    // a0 is d(nu)/d(xnum)
    Real a0 = 1.0e0_rt/(0.5e0_rt*zfermim12(nu));
    Real nudt, nuda, nudz;
    if constexpr (do_derivatives) {
        nudt = a0*xnumdt;
        nuda = a0*xnumda;
        nudz = a0*xnumdz;
    }
    Real nu2  = nu * nu;
    Real nu3  = nu2 * nu;

    Real a1, a2, a3, b, c, d, f1, f2, f3;

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

    // equation 6.7, 6.13 and 6.14
    if (nu >= -20.0_rt &&  nu <= 10.0_rt) {

        Real zeta   = 1.579e5_rt * sf.zbar * sf.zbar * sf.tempi;
        Real zetadt, zetada, zetadz;
        if constexpr (do_derivatives) {
            zetadt = -zeta * sf.tempi;
            zetada = 0.0e0_rt;
            zetadz = 2.0e0_rt * zeta * sf.zbari;
        }

        Real c00 = 1.0e0_rt / (1.0e0_rt + f1 * nu + f2 * nu2 + f3 * nu3);
        Real c01 = f1 + f2 * 2.0e0_rt * nu + f3 * 3.0e0_rt * nu2;
        Real dum = zeta * c00;
        Real dumdt, dumda, dumdz;
        if constexpr (do_derivatives) {
            dumdt = zetadt * c00 + zeta * c01 * nudt;
            dumda = zeta * c01 * nuda;
            dumdz = zetadz * c00 + zeta * c01 * nudz;
        }

        Real z = 1.0e0_rt / dum;
        Real dd00 = std::pow(dum, -2.25_rt);
        Real dd01 = std::pow(dum, -4.55_rt);
        c00  = a1 * z + a2 * dd00 + a3 * dd01;
        c01 = -(a1 * z + 2.25_rt * a2 * dd00 + 4.55_rt * a3 * dd01) * z;

        z = std::exp(c * nu);
        dd00 = b * z * (1.0e0_rt + d * dum);
        Real gum = 1.0e0_rt + dd00;
        Real gumdt, gumda, gumdz;
        if constexpr (do_derivatives) {
            gumdt  = dd00 * c * nudt + b * z * d * dumdt;
            gumda  = dd00 * c * nuda + b * z * d * dumda;
            gumdz  = dd00 * c * nudz + b * z * d * dumdz;
        }

        z   = std::exp(nu);
        a1  = 1.0e0_rt / gum;

        Real bigj = c00 * z * a1;
        Real bigjdt, bigjda, bigjdz;
        if constexpr (do_derivatives) {
            bigjdt = c01 * dumdt * z * a1 + c00 * z * nudt * a1 - c00 * z * a1 * a1 * gumdt;
            bigjda = c01 * dumda * z * a1 + c00 * z * nuda * a1 - c00 * z * a1 * a1 * gumda;
            bigjdz = c01 * dumdz * z * a1 + c00 * z * nudz * a1 - c00 * z * a1 * a1 * gumdz;
        }

        // equation 6.5
        z     = std::exp(zeta + nu);
        dum   = 1.0e0_rt + z;
        a1    = 1.0e0_rt/dum;
        a2    = 1.0e0_rt/bigj;

        sreco   = nu_constants::tfac6 * 2.649e-18_rt * sf.ye * std::pow(sf.zbar, 13.0_rt) * sf.den * bigj * a1;
        if constexpr (do_derivatives) {
            srecodt = sreco * (bigjdt * a2 - z * (zetadt + nudt) * a1);
            srecoda = sreco * (-1.0e0_rt * sf.abari + bigjda * a2 - z * (zetada + nuda) * a1);
            srecodz = sreco * (14.0e0_rt * sf.zbari + bigjdz * a2 - z * (zetadz + nudz) * a1);
        }
    }

}

template <int do_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void sneut5(const Real temp, const Real den,
            const Real abar, const Real zbar,
            Real& snu, Real& dsnudt, Real& dsnudd,
            Real& dsnuda, Real& dsnudz)
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
    snu    = total neutrino loss rate in erg/g/sec
    dsnudt = derivative of snu with temperature
    dsnudd = derivative of snu with density
    dsnuda = derivative of snu with abar
    dsnudz = derivative of snu with zbar
    */

    // initialize
    Real spair{0.0e0_rt};
    Real spairdt{0.0e0_rt};
    Real spairda{0.0e0_rt};
    Real spairdz{0.0e0_rt};

    Real splas{0.0e0_rt};
    Real splasdt{0.0e0_rt};
    Real splasda{0.0e0_rt};
    Real splasdz{0.0e0_rt};

    Real sphot{0.0e0_rt};
    Real sphotdt{0.0e0_rt};
    Real sphotda{0.0e0_rt};
    Real sphotdz{0.0e0_rt};

    Real sbrem{0.0e0_rt};
    Real sbremdt{0.0e0_rt};
    Real sbremda{0.0e0_rt};
    Real sbremdz{0.0e0_rt};

    Real sreco{0.0e0_rt};
    Real srecodt{0.0e0_rt};
    Real srecoda{0.0e0_rt};
    Real srecodz{0.0e0_rt};

    snu     = 0.0e0_rt;
    dsnudt  = 0.0e0_rt;
    dsnudd  = 0.0e0_rt;
    dsnuda  = 0.0e0_rt;
    dsnudz  = 0.0e0_rt;

    if (temp < 1.0e7_rt) return;

    auto sf = get_sneut_factors<do_derivatives>(den, temp, abar, zbar);

    nu_pair<do_derivatives>(sf, spair, spairdt, spairda, spairdz);

    nu_plasma<do_derivatives>(sf, splas, splasdt, splasda, splasdz);

    nu_photo<do_derivatives>(sf, sphot, sphotdt, sphotda, sphotdz);

    nu_brem<do_derivatives>(sf, sbrem, sbremdt, sbremda, sbremdz);

    nu_recomb<do_derivatives>(sf, sreco, srecodt, srecoda, srecodz);

    // convert from erg/cm^3/s to erg/g/s
    // comment these out to duplicate the itoh et al plots

    spair   = spair*sf.deni;
    if constexpr (do_derivatives) {
        spairdt = spairdt*sf.deni;
        spairda = spairda*sf.deni;
        spairdz = spairdz*sf.deni;
    }

    splas   = splas*sf.deni;
    if constexpr (do_derivatives) {
        splasdt = splasdt*sf.deni;
        splasda = splasda*sf.deni;
        splasdz = splasdz*sf.deni;
    }

    sphot   = sphot*sf.deni;
    if constexpr (do_derivatives) {
        sphotdt = sphotdt*sf.deni;
        sphotda = sphotda*sf.deni;
        sphotdz = sphotdz*sf.deni;
    }

    sbrem   = sbrem*sf.deni;
    if constexpr (do_derivatives) {
        sbremdt = sbremdt*sf.deni;
        sbremda = sbremda*sf.deni;
        sbremdz = sbremdz*sf.deni;
    }

    sreco   = sreco*sf.deni;
    if constexpr (do_derivatives) {
        srecodt = srecodt*sf.deni;
        srecoda = srecoda*sf.deni;
        srecodz = srecodz*sf.deni;
    }

    // the total neutrino loss rate
    snu    =  splas + spair + sphot + sbrem + sreco;
    if constexpr (do_derivatives) {
        dsnudt =  splasdt + spairdt + sphotdt + sbremdt + srecodt;
        dsnuda =  splasda + spairda + sphotda + sbremda + srecoda;
        dsnudz =  splasdz + spairdz + sphotdz + sbremdz + srecodz;
    }

}

#endif
