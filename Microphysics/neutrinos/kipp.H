// This routine calculates the neutrino cooling rates using the aproximations from Kippenhahn et al. 1989
// using autodiff

#ifndef KIP_NEUT_H
#define KIP_NEUT_H

#define _USE_MATH_DEFINES

#include<iostream>
#include<AMReX_REAL.H>
#include<fundamental_constants.H>
#include<microphysics_autodiff.H>
// #include<tuple>

using namespace amrex::literals;

template<typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
// pair annhilation
number_t kip_pair(const number_t& rho, const number_t& T9, const number_t& abar, const number_t& zbar){
    // eq. 18.81
    if (T9 < 1.0){
        return (4.9e8_rt/rho) * amrex::Math::powi<3>(T9) * admath::exp(-11.86*T9);
    }
    if (T9 > 3.0){
        return (4.45e15_rt/rho) * amrex::Math::powi<9>(T9);
    }
    return 0.0e0_rt; // Need to figire out what to do in the range 1.0 < T9 < 3.0
}

template<typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
// photoneutrinos
number_t kip_phot(const number_t& rho, const number_t& T9, const number_t& abar, const number_t& zbar){
    // eq. 18.82
    number_t ep_1 = (1.103e13_rt/rho) * amrex::Math::powi<9>(T9) * admath::exp(-5.93e0_rt / T9);
    number_t ep_2 = 0.976e8_rt * amrex::Math::powi<8>(T9) / (1.0e0_rt + 4.2e0_rt * T9);
    number_t rho_bar = 6.446e-6_rt * rho / (T9 + 4.2e0_rt * amrex::Math::powi<2>(T9));
    number_t mu_e = abar / zbar;

    return ep_1 + ep_2 / (mu_e + rho_bar);
}

template<typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
// plasmaneutrinos
number_t kip_plas(const number_t& rho, const number_t& T9, const number_t& abar, const number_t& zbar){
    number_t n_e = rho / (C::m_u);
    // assuming non-degeneracy
    // 18.83
    number_t w0 = admath::sqrt((4.0e0_rt * M_PI * gcem::pow(C::q_e, 2.0e0_rt) * n_e) / C::m_e);
    number_t gamma = (C::hbar * w0) / (C::k_B * T9);
    number_t lambda = (C::k_B * T9) / (C::m_e * gcem::pow(C::c_light, 2.0e0_rt));
    // 18.85
    if (gamma < 1) {
        return (3.356e19_rt/rho) * amrex::Math::powi<6>(lambda) * (1.0e0_rt + 0.0158_rt * amrex::Math::powi<2>(gamma)) * amrex::Math::powi<3>(T9);
    } else {
        return (5.252e20_rt/rho) * admath::pow(lambda, 7.5e0_rt) * admath::pow(T9, 1.5e0_rt) * admath::exp(-gamma);
    }
}

template<typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
// bremsstrahlung neutrinos
number_t kip_brem(const number_t& rho, const number_t& T9, const number_t& abar, const number_t& zbar){
    // eq. 18.86
    number_t T8 = T9 * 10.0e0_rt;
    return 0.76e0_rt * (amrex::Math::powi<2>(zbar)/abar) * amrex::Math::powi<6>(T8);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void kipp(const amrex::Real& temp, const amrex::Real& rho, const amrex::Real& abar,
    const amrex::Real& zbar, amrex::Real& snu, amrex::Real& dsnudt, amrex::Real& dsnudrho, amrex::Real& dsnudz, amrex::Real& dsnuda){

    amrex::Real T9 = temp * 1.0e-9_rt;

    using number_t = autodiff::dual_array<1,4>;

    number_t rho_dual = rho;
    number_t T9_dual = T9;
    number_t abar_dual = abar;
    number_t zbar_dual = zbar;

    autodiff::seed_array(rho_dual, T9_dual, abar_dual, zbar_dual);

    number_t kip_pair_dual = kip_pair(rho_dual, T9_dual, abar_dual, zbar_dual);
    number_t kip_phot_dual = kip_phot(rho_dual, T9_dual, abar_dual, zbar_dual);
    number_t kip_plas_dual = kip_plas(rho_dual, T9_dual, abar_dual, zbar_dual);
    number_t kip_brem_dual = kip_brem(rho_dual, T9_dual, abar_dual, zbar_dual);

    // sum of all neutrino cooling rates

    number_t snu_dual = kip_pair_dual + kip_phot_dual + kip_plas_dual + kip_brem_dual;

    // std::tie(dsnudrho, dsnudt, dsnuda, dsnudz) = autodiff::derivative(snu_dual);

    auto derivatives = autodiff::derivative(snu_dual);
    dsnudrho = derivatives(1);
    dsnudt = derivatives(2);
    dsnuda = derivatives(3);
    dsnudz = derivatives(4);

    snu = autodiff::val(snu_dual);
}
#endif // KIP_NEUT_H