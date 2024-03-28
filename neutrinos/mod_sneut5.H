#ifndef SNEUT5_H
#define SNEUT5_H

#define _USE_MATH_DEFINES

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

using namespace amrex;

template <int do_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void sneut5(const Real temp, const Real dens, const Real abar,
            const Real zbar, Real& snu, Real& dsnudt, Real& dsnudd,
            Real& dsnuda, Real& dsnudz)
{

Real t9, t8, rho, ep_1, ep_2, w0, gamma, lambda,
     w1, w2, rho_bar, n_e, mu_e;

t9 = temp * 1.0e-9_rt;
t8 = temp * 1.0e-8_rt;
rho = dens;
snu = 0.0e0_rt;
mu_e = abar/zbar;
dsnudt  = 0.0e0_rt;
dsnudd  = 0.0e0_rt;
dsnuda  = 0.0e0_rt;
dsnudz  = 0.0e0_rt;

// constants 
constexpr Real q_e = 4.8032e-10_rt; // Electron charge cm^3/2 g^1/2 s^-1
constexpr Real m_e = 9.1094e-28_rt; // Mass of electron in g
//constexpr Real n_e = rho/(m_e * 1.66054e-24_rt);
//constexpr Real mu_e = abar/zbar;
constexpr Real hbar = 1.0546e-27_rt; // h/2pi in cm^2 g s^-1
constexpr Real c = 3.0e10_rt; // speed of light in cm s^-1
constexpr Real quart = 1.0e0_rt/4.0e0_rt;
constexpr Real twothd = 2.0e0_rt/ 3.0e0_rt;
constexpr Real k = 1.3807e-16_rt; // Boltzmann constant cm^2 g s^-2 K^-1

n_e = rho/(m_e * 1.66054e-24_rt);
// initialize
// Epsilon (ep_) has units (erg g^-1 s^-1)

Real ep_pair{0.0e0_rt};
Real ep_pairdt{0.0e0_rt};
Real ep_pairda{0.0e0_rt};
Real ep_pairdz{0.0e0_rt};

Real ep_photo{0.0e0_rt};
Real ep_photdt{0.0e0_rt};
Real ep_photda{0.0e0_rt};
Real ep_photdz{0.0e0_rt};

Real ep_plasma{0.0e0_rt};
Real ep_plasmadt{0.0e0_rt};
Real ep_plasmada{0.0e0_rt};
Real ep_plasmadz{0.0e0_rt};

Real ep_brem{0.0e0_rt};
Real ep_bremdt{0.0e0_rt};
Real ep_bremda{0.0e0_rt};
Real ep_bremdz{0.0e0_rt};

if (temp < 1.0e7_rt) return;

// Pair annhilation neutrinos
// for reactions like e+ + e- => nu + nubar
// eq 18.81
if(t9 < 1){
    ep_pair = (4.9e8_rt/rho) * std::pow(t9,3.0e0_rt) * std::exp(-11.86e0_rt*t9);
} else if(t9 > 3.0e0_rt){
    ep_pair = (4.45e15_rt/rho) * std::pow(t9, 9.0e0_rt);
}

// Photoneutrinos
// for reactions like 
// eq 18.82
ep_1 = 1.103e13_rt * std::pow(rho, -1.0e0_rt) * std::pow(t9, 9.0e0_rt) * std::exp(-5.93e0_rt/t9);
ep_2 = 0.976_rt * 1e8_rt * std::pow(t9, 8.0e0_rt) / (1.0e0_rt + 4.2e0_rt * t9);
rho_bar = 6.446e-6_rt * rho * std::pow(t9, -1.0e0_rt) / (1.0e0_rt + 4.2e0_rt * t9);
ep_photo = ep_1 + ep_2 * 1.0e0_rt/(mu_e + rho_bar);

// Plasmaneutrinos
// for reactions like gamma_plasm -> nu + nubar
// eq. 18.83
w1 = (4.0e0_rt * M_PI * std::pow(q_e, 2.0e0_rt) * n_e)/m_e;
w2 = 1 + (std::pow(hbar, 2.0e0_rt)/std::pow(m_e * c, 2.0e0_rt)) * std::pow((3.0e0_rt * std::pow(M_PI, 2.0e0_rt) * n_e), twothd);
w0 = std::sqrt(w1) * std::pow(w2, quart);

gamma = (hbar * w0)/(k * temp);
lambda = (k * temp)/(m_e * std::pow(c,2.0e0_rt));

// eq 18.85
if(gamma < 1){
    ep_plasma = 3.356e19_rt * std::pow(rho, -1.0e0_rt) * std::pow(lambda, 6.0e0_rt) * (1.0e0_rt + 0.0158_rt*std::pow(gamma,2.0e0_rt)) * std::pow(t9,3.0e0_rt);
} else if(gamma > 1){
    ep_plasma = 5.252e20_rt * std::pow(rho, -1.0e0_rt) * std::pow(lambda, 7.5e0_rt) * std::pow(t9, 1.5e0_rt) * std::exp(-lambda);
}
// Bremsstrahlung neutrinos
// for large rho
// eq 18.86
ep_brem = 0.76e0_rt * (std::pow(zbar, 2.0e0_rt)/abar) * std::pow(t8, 6.0e0_rt);

// the total neutrino loss rate

snu = ep_pair + ep_photo + ep_plasma + ep_brem;
if constexpr (do_derivatives) {
    dsnudt = 0;
    dsnuda = 0;
    dsnudz = 0;
    dsnudd = 0;
} 

}

#endif