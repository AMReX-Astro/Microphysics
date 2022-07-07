#ifndef REACLIB_RATES_H
#define REACLIB_RATES_H

#include <AMReX.H>
#include <AMReX_Print.H>

#include <actual_network.H>

using namespace Rates;

struct rate_eval_t {
    Array1D<Real, 1, NumRates>  screened_rates;
    Array1D<Real, 1, NumRates>  dscreened_rates_dT;
    Array1D<Real, NrateReaclib+1, NrateReaclib+NrateTabular> add_energy_rate;
};

struct tf_t {
    Real T9;
    Real T9i;
    Real T943i;
    Real T923i;
    Real T913i;
    Real T913;
    Real T923;
    Real T953;
    Real lnT9;
};

AMREX_GPU_HOST_DEVICE AMREX_INLINE
tf_t evaluate_tfactors(const Real T)
{

    tf_t tf;
    tf.T9 = T / 1.e9_rt;
    tf.T9i = 1.0_rt / tf.T9;
    tf.T913 = std::cbrt(tf.T9);
    tf.T913i = 1.0_rt / tf.T913;
    tf.T923i = tf.T913i * tf.T913i;
    tf.T943i = tf.T9i * tf.T913i;
    tf.T923 = tf.T913 * tf.T913;
    tf.T953 = tf.T9 * tf.T923;
    tf.lnT9 = std::log(tf.T9);

    return tf;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_he4_he4__c12(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // he4 + he4 + he4 --> c12

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // fy05r
    ln_set_rate =  -24.3505 + -4.12656 * tfactors.T9i + -13.49 * tfactors.T913i + 21.4259 * tfactors.T913
                         + -1.34769 * tfactors.T9 + 0.0879816 * tfactors.T953 + -13.1653 * tfactors.lnT9;

    dln_set_rate_dT9 =  4.12656 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -13.49 * tfactors.T943i + (1.0/3.0) * 21.4259 * tfactors.T923i
                              + -1.34769 + (5.0/3.0) * 0.0879816 * tfactors.T923 + -13.1653 * tfactors.T9i;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

    // fy05r
    ln_set_rate =  -11.7884 + -1.02446 * tfactors.T9i + -23.57 * tfactors.T913i + 20.4886 * tfactors.T913
                         + -12.9882 * tfactors.T9 + -20.0 * tfactors.T953 + -2.16667 * tfactors.lnT9;

    dln_set_rate_dT9 =  1.02446 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -23.57 * tfactors.T943i + (1.0/3.0) * 20.4886 * tfactors.T923i
                              + -12.9882 + (5.0/3.0) * -20.0 * tfactors.T923 + -2.16667 * tfactors.T9i;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

    // fy05n
    ln_set_rate =  -0.971052 + -37.06 * tfactors.T913i + 29.3493 * tfactors.T913
                         + -115.507 * tfactors.T9 + -10.0 * tfactors.T953 + -1.33333 * tfactors.lnT9;

    dln_set_rate_dT9 =  + -(1.0/3.0) * -37.06 * tfactors.T943i + (1.0/3.0) * 29.3493 * tfactors.T923i
                              + -115.507 + (5.0/3.0) * -10.0 * tfactors.T923 + -1.33333 * tfactors.T9i;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_c12__o16(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // c12 + he4 --> o16

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // nac2 
    ln_set_rate =  254.634 + -1.84097 * tfactors.T9i + 103.411 * tfactors.T913i + -420.567 * tfactors.T913
                         + 64.0874 * tfactors.T9 + -12.4624 * tfactors.T953 + 137.303 * tfactors.lnT9;

    dln_set_rate_dT9 =  1.84097 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 103.411 * tfactors.T943i + (1.0/3.0) * -420.567 * tfactors.T923i
                              + 64.0874 + (5.0/3.0) * -12.4624 * tfactors.T923 + 137.303 * tfactors.T9i;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

    // nac2 
    ln_set_rate =  69.6526 + -1.39254 * tfactors.T9i + 58.9128 * tfactors.T913i + -148.273 * tfactors.T913
                         + 9.08324 * tfactors.T9 + -0.541041 * tfactors.T953 + 70.3554 * tfactors.lnT9;

    dln_set_rate_dT9 =  1.39254 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 58.9128 * tfactors.T943i + (1.0/3.0) * -148.273 * tfactors.T923i
                              + 9.08324 + (5.0/3.0) * -0.541041 * tfactors.T923 + 70.3554 * tfactors.T9i;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_n14__f18(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // n14 + he4 --> f18

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10n
    ln_set_rate =  21.5339 + -36.2504 * tfactors.T913i
                         + -5.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    dln_set_rate_dT9 =  + -(1.0/3.0) * -36.2504 * tfactors.T943i
                              + (5.0/3.0) * -5.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

    // il10r
    ln_set_rate =  13.8995 + -10.9656 * tfactors.T9i + -5.6227 * tfactors.T913i
                         + -1.5 * tfactors.lnT9;

    dln_set_rate_dT9 =  10.9656 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -5.6227 * tfactors.T943i
                              + -1.5 * tfactors.T9i;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

    // il10r
    ln_set_rate =  0.196838 + -5.16034 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    dln_set_rate_dT9 =  5.16034 * tfactors.T9i * tfactors.T9i
                              + -1.5 * tfactors.T9i;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_f18__p_ne21(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // f18 + he4 --> p + ne21

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // rpsmr
    ln_set_rate =  49.7863 + -1.84559 * tfactors.T9i + 21.4461 * tfactors.T913i + -73.252 * tfactors.T913
                         + 2.42329 * tfactors.T9 + -0.077278 * tfactors.T953 + 40.7604 * tfactors.lnT9;

    dln_set_rate_dT9 =  1.84559 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 21.4461 * tfactors.T943i + (1.0/3.0) * -73.252 * tfactors.T923i
                              + 2.42329 + (5.0/3.0) * -0.077278 * tfactors.T923 + 40.7604 * tfactors.T9i;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_c12__n13(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // c12 + p --> n13

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ls09n
    ln_set_rate =  17.1482 + -13.692 * tfactors.T913i + -0.230881 * tfactors.T913
                         + 4.44362 * tfactors.T9 + -3.15898 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    dln_set_rate_dT9 =  + -(1.0/3.0) * -13.692 * tfactors.T943i + (1.0/3.0) * -0.230881 * tfactors.T923i
                              + 4.44362 + (5.0/3.0) * -3.15898 * tfactors.T923 + -0.666667 * tfactors.T9i;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

    // ls09r
    ln_set_rate =  17.5428 + -3.77849 * tfactors.T9i + -5.10735 * tfactors.T913i + -2.24111 * tfactors.T913
                         + 0.148883 * tfactors.T9 + -1.5 * tfactors.lnT9;

    dln_set_rate_dT9 =  3.77849 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -5.10735 * tfactors.T943i + (1.0/3.0) * -2.24111 * tfactors.T923i
                              + 0.148883 + -1.5 * tfactors.T9i;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_n13__p_o16(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // n13 + he4 --> p + o16

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cf88n
    ln_set_rate =  40.4644 + -35.829 * tfactors.T913i + -0.530275 * tfactors.T913
                         + -0.982462 * tfactors.T9 + 0.0808059 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    dln_set_rate_dT9 =  + -(1.0/3.0) * -35.829 * tfactors.T943i + (1.0/3.0) * -0.530275 * tfactors.T923i
                              + -0.982462 + (5.0/3.0) * 0.0808059 * tfactors.T923 + -0.666667 * tfactors.T9i;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_o16__ne20(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // o16 + he4 --> ne20

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // co10r
    ln_set_rate =  9.50848 + -12.7643 * tfactors.T9i + -3.65925 * tfactors.T913
                         + 0.714224 * tfactors.T9 + -0.00107508 * tfactors.T953 + -1.5 * tfactors.lnT9;

    dln_set_rate_dT9 =  12.7643 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -3.65925 * tfactors.T923i
                              + 0.714224 + (5.0/3.0) * -0.00107508 * tfactors.T923 + -1.5 * tfactors.T9i;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

    // co10r
    ln_set_rate =  3.88571 + -10.3585 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    dln_set_rate_dT9 =  10.3585 * tfactors.T9i * tfactors.T9i
                              + -1.5 * tfactors.T9i;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

    // co10n
    ln_set_rate =  23.903 + -39.7262 * tfactors.T913i + -0.210799 * tfactors.T913
                         + 0.442879 * tfactors.T9 + -0.0797753 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    dln_set_rate_dT9 =  + -(1.0/3.0) * -39.7262 * tfactors.T943i + (1.0/3.0) * -0.210799 * tfactors.T923i
                              + 0.442879 + (5.0/3.0) * -0.0797753 * tfactors.T923 + -0.666667 * tfactors.T9i;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_c14__o18(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // c14 + he4 --> o18

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10r
    ln_set_rate =  11.8309 + -10.3983 * tfactors.T9i + -3.83188 * tfactors.T913
                         + 1.64358 * tfactors.T9 + -0.177785 * tfactors.T953 + -1.5 * tfactors.lnT9;

    dln_set_rate_dT9 =  10.3983 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -3.83188 * tfactors.T923i
                              + 1.64358 + (5.0/3.0) * -0.177785 * tfactors.T923 + -1.5 * tfactors.T9i;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

    // il10r
    ln_set_rate =  -23.805 + -2.06876 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    dln_set_rate_dT9 =  2.06876 * tfactors.T9i * tfactors.T9i
                              + -1.5 * tfactors.T9i;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

    // il10n
    ln_set_rate =  18.4877 + -31.7222 * tfactors.T913i + 11.3923 * tfactors.T913
                         + -9.92249 * tfactors.T9 + -2.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    dln_set_rate_dT9 =  + -(1.0/3.0) * -31.7222 * tfactors.T943i + (1.0/3.0) * 11.3923 * tfactors.T923i
                              + -9.92249 + (5.0/3.0) * -2.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;

}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
fill_reaclib_rates(const tf_t& tfactors, rate_eval_t& rate_eval)
{

    Real rate;
    Real drate_dT;

    rate_he4_he4_he4__c12(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_he4_he4__c12) = rate;
    rate_eval.dscreened_rates_dT(k_he4_he4_he4__c12) = drate_dT;

    rate_he4_c12__o16(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_c12__o16) = rate;
    rate_eval.dscreened_rates_dT(k_he4_c12__o16) = drate_dT;

    rate_he4_n14__f18(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_n14__f18) = rate;
    rate_eval.dscreened_rates_dT(k_he4_n14__f18) = drate_dT;

    rate_he4_f18__p_ne21(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_f18__p_ne21) = rate;
    rate_eval.dscreened_rates_dT(k_he4_f18__p_ne21) = drate_dT;

    rate_p_c12__n13(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_c12__n13) = rate;
    rate_eval.dscreened_rates_dT(k_p_c12__n13) = drate_dT;

    rate_he4_n13__p_o16(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_n13__p_o16) = rate;
    rate_eval.dscreened_rates_dT(k_he4_n13__p_o16) = drate_dT;

    rate_he4_o16__ne20(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_o16__ne20) = rate;
    rate_eval.dscreened_rates_dT(k_he4_o16__ne20) = drate_dT;

    rate_he4_c14__o18(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_c14__o18) = rate;
    rate_eval.dscreened_rates_dT(k_he4_c14__o18) = drate_dT;


}

#endif
