// This routine calculates the neutrino cooling rates using the approximations from Kippenhahn et al. 1989
// using autodiff

#ifndef KIP_NEUT_H
#define KIP_NEUT_H

#include<iostream>
#include<AMReX_REAL.H>
#include<fundamental_constants.H>
#include<microphysics_autodiff.H>

using namespace amrex::literals;

template<typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
// pair annihilation
number_t kip_pair(const number_t& temp, const amrex::Real& rho,
                  const number_t& abar, const number_t& zbar){
    number_t T9 = temp * 1.0e-9_rt;
    // eq. 18.81 - fixed using woosley lecture notes
    if (T9 < 2.0e0_rt){
        return (4.9e18_rt/rho) * amrex::Math::powi<3>(T9) * admath::exp(-11.86/T9);
    }
    else{
        return (3.2e15_rt/rho) * amrex::Math::powi<9>(T9);
    }
}

template<typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
// photoneutrinos
number_t kip_phot(const number_t& temp, const amrex::Real& rho,
                  const number_t& abar, const number_t& zbar){
    number_t T9 = temp * 1.0e-9_rt;
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
number_t kip_plas(const number_t& temp, const amrex::Real& rho,
                  const number_t& abar, const number_t& zbar){
    number_t n_e = rho / (C::m_u * abar/zbar);
    // assuming non-degeneracy
    // 18.83 replaced by aprox in woosley lecture notes
    number_t w0 = admath::sqrt((4.0e0_rt * M_PI * amrex::Math::powi<2>(C::q_e) * n_e) / C::m_e);
    number_t gamma = (C::hbar * w0) / (C::k_B * temp);
    number_t lambda = (C::k_B * temp) / (C::m_e * amrex::Math::powi<2>(C::c_light));
    // 18.85
   if (gamma <= 1.0e0_rt) {
       // return (7.4e21_rt * amrex::Math::powi<6>(C::hbar * w0/(C::m_e * amrex::Math::powi<2>(C::c_light)))
       //         * amrex::Math::powi<3>(lambda))/rho;
        number_t epsilon_p = (C::hbar * w0) / (C::m_e * C::c_light * C::c_light);
        return (7.4e21_rt * amrex::Math::powi<6>(epsilon_p) * amrex::Math::powi<3>(lambda)) / rho;

   } else{
        return (3.3e21_rt * admath::pow((C::hbar * w0/(C::m_e * amrex::Math::powi<2>(C::c_light))),7.5e0_rt)
                * admath::pow(lambda, 3.0e0_rt/2.0e0_rt) * admath::exp(-gamma))/rho;
    }

}

template<typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
// bremsstrahlung neutrinos
number_t kip_brem(const number_t& temp, const amrex::Real& rho,
                  const number_t& abar, const number_t& zbar){
    // eq. 18.86
    number_t T8 = temp * 1.0e-8_rt;
    return 0.76e0_rt * (amrex::Math::powi<2>(zbar)/abar) * amrex::Math::powi<6>(T8);
}

template<typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
number_t kipp(const number_t& temp, const amrex::Real rho,
              const number_t& abar, const number_t& zbar, amrex::Real& pair,
              amrex::Real& phot, amrex::Real& plas, amrex::Real& brem){

    if (temp < 1.0e7_rt){
        return 0.0e0_rt;
    }

    // calculating individual contributions of all processes
    number_t kip_pair_dual = kip_pair(temp, rho, abar, zbar);
    number_t kip_phot_dual = kip_phot(temp, rho, abar, zbar);
    number_t kip_plas_dual = kip_plas(temp, rho, abar, zbar);
    number_t kip_brem_dual = kip_brem(temp, rho, abar, zbar);

    // extracting contributions to store in plotfile
    pair = autodiff::val(kip_pair_dual);
    phot = autodiff::val(kip_phot_dual);
    plas = autodiff::val(kip_plas_dual);
    brem = autodiff::val(kip_brem_dual);

    // total neutrino cooling rate
    number_t snu = kip_pair_dual + kip_phot_dual + kip_plas_dual + kip_brem_dual;

    return snu;

}

template <int do_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void kipp(const amrex::Real& temp, const amrex::Real& rho, const amrex::Real& abar,
          const amrex::Real& zbar, amrex::Real& snu, amrex::Real& dsnudt,
          amrex::Real& dsnudrho, amrex::Real& dsnuda, amrex::Real& dsnudz,
          amrex::Real& pair, amrex::Real& phot, amrex::Real& plas, amrex::Real& brem) {

    /*
        input:
        temp = temperature
        rho = density
        abar = mean atomic weight
        zbar = mean charge

        output:
        snu = total neutrino cooling rate in erg/g/s
        dsnudt = derivative of snu with respect to temp
        dsnudz = derivative of snu with respect to zbar
        dsnuda = derivative of snu with respect to abar
        pair = pair annihilation contribution in erg/g/s
        phot = photoneutrino contribution in erg/g/s
        plas = plasma neutrino contribution in erg/g/s
        brem = bremsstrahlung contribution rate in erg/g/s
    */

    // autodiff wrapper
    if constexpr (do_derivatives){

    using dual_t = autodiff::dual_array<1,3>;
    dual_t temp_dual = temp;
    dual_t abar_dual = abar;
    dual_t zbar_dual = zbar;

    autodiff::seed_array(temp_dual, abar_dual, zbar_dual);
    dual_t snu_dual = kipp(temp_dual, rho, abar_dual, zbar_dual, pair, phot, plas, brem);

    snu = autodiff::val(snu_dual);
    const auto& grad = autodiff::derivative(snu_dual);
    dsnudt = grad(1);
    dsnudrho = 0.0e0_rt;
    dsnuda = grad(2);
    dsnudz = grad(3);
    } else {
        snu = kipp(temp, rho, abar, zbar, pair, phot, plas, brem);
        dsnudt = 0.0e0_rt;
        dsnudrho = 0.0e0_rt;
        dsnuda = 0.0e0_rt;
        dsnudz = 0.0e0_rt;
    }
}

//overloading to pass arguments in all cases
template <int do_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void kipp(const amrex::Real temp, const amrex::Real den,
          const amrex::Real abar, const amrex::Real zbar,
          amrex::Real& snu, amrex::Real& dsnudt, amrex::Real& dsnudd,
          amrex::Real& dsnuda, amrex::Real& dsnudz){

    amrex::Real pair, phot, plas, brem;
    kipp<do_derivatives>(temp, den, abar, zbar, snu, dsnudt, dsnudd, dsnuda,
                         dsnudz, pair, phot, plas, brem);
}

#endif // KIP_NEUT_H