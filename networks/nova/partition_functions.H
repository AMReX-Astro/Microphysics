#ifndef PARTITION_FUNCTIONS_H
#define PARTITION_FUNCTIONS_H

#include <numbers>

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Math.H>
#include <AMReX_ConstexprFor.H>

#include <tfactors.H>
#include <fundamental_constants.H>
#include <network_properties.H>
#include <interp_tools.H>

using namespace amrex::literals;
using namespace Species;

namespace part_fun {

    constexpr int npts_1 = 72;

    // this is T9

    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> temp_array_1 = {
        0.01, 0.15, 0.2, 0.3, 0.4,
        0.5, 0.6, 0.7, 0.8, 0.9,
        1.0, 1.5, 2.0, 2.5, 3.0,
        3.5, 4.0, 4.5, 5.0, 6.0,
        7.0, 8.0, 9.0, 10.0, 12.0,
        14.0, 16.0, 18.0, 20.0, 22.0,
        24.0, 26.0, 28.0, 30.0, 35.0,
        40.0, 45.0, 50.0, 55.0, 60.0,
        65.0, 70.0, 75.0, 80.0, 85.0,
        90.0, 95.0, 100.0, 105.0, 110.0,
        115.0, 120.0, 125.0, 130.0, 135.0,
        140.0, 145.0, 150.0, 155.0, 160.0,
        165.0, 170.0, 175.0, 180.0, 190.0,
        200.0, 210.0, 220.0, 230.0, 240.0,
        250.0, 275.0,
    };


    // O14

    // this is log(partition function)

    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> O14_pf_array = {
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.04879016416943205,
        0.10436001532424286, 0.1823215567939546, 0.27763173659827955, 0.3852624007906449, 0.5007752879124892,
        0.6097655716208943, 0.712949807856125, 0.8153648132841944, 0.9082585601768908, 1.1085626195212777,
        1.269760544863939, 1.3937663759585917, 1.4906543764441336, 1.5644405465033646, 1.6193882432872684,
        1.6582280766035324, 1.6845453849209058, 1.6974487897568136, 1.7029282555214393, 1.6974487897568136,
        1.6845453849209058, 1.665818245870208, 1.6409365794934714, 1.6094379124341003, 1.5748464676644813,
        1.534714366238164, 1.4883995840570443, 1.4422019930581866, 1.3887912413184778, 1.33500106673234,
        1.2781522025001875, 1.2178757094949273, 1.1568811967920856, 1.091923300517313, 1.0260415958332743,
        0.9555114450274363, 0.8878912573524571, 0.8153648132841944, 0.7419373447293773, 0.5933268452777344,
        0.4382549309311553, 0.28517894223366247, 0.12221763272424911, -0.040821994520255166, -0.2045671657412744,
        -0.3696154552144672, -0.7874578600311866,
    };

    constexpr amrex::Real O14_pf_threshold_T9 = 10.0;

    // O15

    // this is log(partition function)

    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> O15_pf_array = {
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.04879016416943205,
        0.10436001532424286, 0.1906203596086497, 0.3148107398400336, 0.47000362924573563, 0.6418538861723947,
        0.8329091229351039, 1.0296194171811581, 1.2325602611778486, 1.43746264769429, 1.9572739077056285,
        2.4849066497880004, 3.0155349008501706, 3.55820113047182, 4.109233174715851, 4.672828834461906,
        5.241747015059643, 5.814130531825066, 6.391917113392602, 6.966024187106113, 7.549609165154532,
        8.125630988477065, 8.702842538302868, 9.277999020449997, 9.852194258148577, 10.42821608147111,
        11.00043178410354, 11.571194373094205, 12.144197241812087, 12.7098736543592, 13.279367126213993,
        13.845069360205818, 14.414347059052979, 14.975531474761027, 15.541842221869874, 16.105010411409665,
        16.666217059468007, 17.22995316637645, 17.793321303930423, 18.35560874720865, 19.47847103810022,
        20.602227508569356, 21.723897717254317, 22.846724264043022, 23.969756828847586, 25.09397871172002,
        26.22043406223961, 29.03956608598326,
    };

    constexpr amrex::Real O15_pf_threshold_T9 = 10.0;

    // O16

    // this is log(partition function)

    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> O16_pf_array = {
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.02955880224154443,
        0.08617769624105241, 0.16551443847757333, 0.27763173659827955, 0.4252677354043441, 0.6043159668533296,
        0.7975071958841882, 1.007957920399979, 1.2296405510745139, 1.4609379041156563, 2.0681277817795625,
        2.70805020110221, 3.3775875160230218, 4.062165663857866, 4.7535901911063645, 5.4510384535657,
        6.142037405587356, 6.828712071641684, 7.512071245835466, 8.185907481482324, 8.85366542803745,
        9.517825071724143, 10.177324218165856, 10.831706855275558, 11.48246625748552, 12.128111104060462,
        12.774223335915433, 13.41652441595382, 14.054527458434775, 14.690979295318174, 15.324022551808413,
        15.956752500549557, 16.588099280204055, 17.216707939626428, 17.842646370492925, 18.469470908121796,
        19.09362521719479, 19.715407911546766, 20.337603356134426, 20.954377557909798, 22.19802884605391,
        23.43796058076729, 24.670655986211848, 25.905049387238495, 27.136724794113768, 28.368185181905268,
        29.602320498988682, 32.68726692127686,
    };

    constexpr amrex::Real O16_pf_threshold_T9 = 10.0;

    // O17

    // this is log(partition function)

    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> O17_pf_array = {
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.009950330853168092, 0.009950330853168092,
        0.01980262729617973, 0.02955880224154443, 0.03922071315328133, 0.03922071315328133, 0.058268908123975824,
        0.0769610411361284, 0.09531017980432493, 0.12221763272424911, 0.13976194237515863, 0.2151113796169455,
        0.30010459245033816, 0.4054651081081644, 0.5364933705145685, 0.6830968447064438, 0.8458682675776092,
        1.0260415958332743, 1.2178757094949273, 1.4182774069729414, 1.631199404215613, 2.1972245773362196,
        2.803360380906535, 3.4339872044851463, 4.07753744390572, 4.727387818712341, 5.3706380281276624,
        6.018593214496234, 6.660575149839686, 7.299797366758161, 7.937374696163295, 8.573573524852344,
        9.20833836930551, 9.84161214881804, 10.474467099121865, 11.107960231903714, 11.744037185933616,
        12.375815420117268, 13.008074231002201, 13.641157170819497, 14.27293540500315, 14.907433858481587,
        15.541842221869874, 16.176364559082295, 16.811242831518264, 17.445170652418238, 18.081003376382203,
        18.713350357915186, 19.35284482498281, 19.98721115537519, 20.6234455016642, 21.895747974180978,
        23.17427093505873, 24.448959264183063, 25.72721214289187, 27.01298140785541, 28.298850488504204,
        29.585466167433697, 32.8184069217693,
    };

    constexpr amrex::Real O17_pf_threshold_T9 = 2.0;

    // F17

    // this is log(partition function)

    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> F17_pf_array = {
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.009950330853168092, 0.01980262729617973, 0.02955880224154443, 0.04879016416943205,
        0.058268908123975824, 0.0769610411361284, 0.08617769624105241, 0.09531017980432493, 0.12221763272424911,
        0.13976194237515863, 0.1570037488096647, 0.17395330712343798, 0.19885085874516517, 0.26236426446749106,
        0.3293037471426003, 0.41871033485818504, 0.5247285289349821, 0.6471032420585384, 0.7884573603642703,
        0.9400072584914712, 1.1019400787607843, 1.2781522025001875, 1.4678743481123135, 1.9810014688665833,
        2.5416019934645457, 3.1441522786722644, 3.7612001156935624, 4.388257184424518, 5.0238805208462765,
        5.655991810819852, 6.2878585601617845, 6.917705609835305, 7.549609165154532, 8.183118079394745,
        8.814330422638774, 9.441452092939569, 10.077440859659566, 10.707728780601661, 11.33857207782545,
        11.970350312009105, 12.601487417784837, 13.233904752137237, 13.864300722133706, 14.498607402670718,
        15.131918791619999, 15.765697263786848, 16.3957273875566, 17.034386382832476, 17.667783558986645,
        18.301897207962398, 18.939474537367534, 19.571252771551187, 20.207427671356875, 21.479387816667746,
        22.755353682242777, 24.033808850340435, 25.311289864099532, 26.592562750080187, 27.877881193860073,
        29.165735482166713, 32.393195050726305,
    };

    constexpr amrex::Real F17_pf_threshold_T9 = 1.0;

    // F18

    // this is log(partition function)

    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, npts_1-1> F18_pf_array = {
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.01980262729617973, 0.058268908123975824, 0.11332868530700327,
        0.19885085874516517, 0.28517894223366247, 0.3852624007906449, 0.4762341789963717, 0.6523251860396901,
        0.8064758658669485, 0.9360933591703348, 1.0438040521731147, 1.1378330018213911, 1.5303947050936475,
        1.742219023667919, 1.9387416595767009, 2.125847914493992, 2.312535423847214, 2.4932054526026954,
        2.681021528714291, 2.8791984572980396, 3.0819099697950434, 3.292126286607793, 3.856510295497887,
        4.468204330914934, 5.111987788356544, 5.768320995793772, 6.431331081933479, 7.098375638590786,
        7.757906208351747, 8.420682291035394, 9.080231686629162, 9.740968623038354, 10.395130356885344,
        11.052476048529304, 11.711776323715394, 12.367340793126296, 13.025852477023484, 13.68426227135432,
        14.346138809026444, 15.003353980360327, 15.665538935316304, 16.325109820342647, 16.98938901690174,
        17.654962870557586, 18.319754825362406, 18.985994553002424, 19.650321295026878, 20.316800228504665,
        20.985630101413904, 21.65148513968584, 22.32265341352701, 22.991259485170836, 24.334183749590636,
        25.679092894547672, 27.02771463936839, 28.377709063416525, 29.735155270198756, 31.093627125719348,
        32.45933485323085, 35.88424876151032,
    };

    constexpr amrex::Real F18_pf_threshold_T9 = 1.5;


    // interpolation routine

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void interpolate_pf(const amrex::Real t9, int idx,
                        const T& temp_array, const T& pf_array,
                        amrex::Real& logpf, amrex::Real& dlogpf_dT9) {

        // if the index is -1, we are either outside of the bounds or
        // we never computed it.  try recomputing here
        if (idx == -1) {
            idx = interp_net::find_index(t9, temp_array);
        }

        // if the index is still -1, then we are out-of-bounds, so
        // just skip the interpolation

        if (idx != -1) {

            // we have idx such that:
            // temp_array[idx] <= t9 < temp_array[idx+1]
            AMREX_ASSERT(t9 >= temp_array(idx) && t9 < temp_array(idx+1));

            // find the derivative: dlog(pf)/dT9 -- this is
            // (log(pf_{i+1}) - log(pf_i)) / (T_{i+1} - T_i)

            dlogpf_dT9 = (pf_array(idx+1) - pf_array(idx)) /
                           (temp_array(idx+1) - temp_array(idx));

            // find log(PF)

            logpf = pf_array(idx) + dlogpf_dT9 * (t9 - temp_array(idx));

        } else {

            // T < the smallest T or >= the largest T in the partition function table
            logpf = 0.0;
            dlogpf_dT9 = 0.0;

        }

    }

    struct pf_cache_t {
        // Store the coefficient and derivative adjacent in memory, as they're
        // always accessed at the same time.
        amrex::Array2D<amrex::Real, 1, NumSpecTotal, 1, 2, amrex::Order::C> data;
        int index_temp_array_1{-1};

        AMREX_GPU_HOST_DEVICE AMREX_INLINE
        pf_cache_t() {
            // The entries will be default-initialized to -900,
            // log(partition_function) should never be -900.
            for (int i = 1; i <= NumSpecTotal; ++i) {
                    data(i, 1) = -900.0_rt;
                    data(i, 2) = -900.0_rt;
            }
        }
    };

}

// get the partition function for nucleus inuc.  Here pf_cache
// should already contain the index into the temp / pf arrays

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_partition_function(const int inuc,
                            [[maybe_unused]] const amrex::Real T9,
                            const part_fun::pf_cache_t& pf_cache,
                            amrex::Real& logpf, amrex::Real& dlogpf_dT9) {

    logpf = 0.0_rt;
    dlogpf_dT9 = 0.0_rt;

    // inuc is the 1-based index for the species

    switch (inuc) {

    case O14:
        if (T9 > part_fun::O14_pf_threshold_T9) {
            part_fun::interpolate_pf(T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::O14_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case O15:
        if (T9 > part_fun::O15_pf_threshold_T9) {
            part_fun::interpolate_pf(T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::O15_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case O16:
        if (T9 > part_fun::O16_pf_threshold_T9) {
            part_fun::interpolate_pf(T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::O16_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case O17:
        if (T9 > part_fun::O17_pf_threshold_T9) {
            part_fun::interpolate_pf(T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::O17_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case F17:
        if (T9 > part_fun::F17_pf_threshold_T9) {
            part_fun::interpolate_pf(T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::F17_pf_array, logpf, dlogpf_dT9);
        }
        break;

    case F18:
        if (T9 > part_fun::F18_pf_threshold_T9) {
            part_fun::interpolate_pf(T9, pf_cache.index_temp_array_1, part_fun::temp_array_1, part_fun::F18_pf_array, logpf, dlogpf_dT9);
        }
        break;


    default:

        logpf = 0.0_rt;
        dlogpf_dT9 = 0.0_rt;

    }

}


// a version of get_partition_function where we have not cached the
// index into the arrays. Uses tfactors instead of T9 as argument,
// and returns pf and dpf_dT instead of logpf and dlogpf_dT9

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_partition_function(const int inuc,
                            [[maybe_unused]] const tf_t& tfactors,
                            amrex::Real& pf, amrex::Real& dpf_dT) {

    // create a dummy cache.  This will initialize the indices
    // to -1 and we will force a search for the correct index
    // when we do the interpolation.

    part_fun::pf_cache_t pf_cache;

    amrex::Real logpf, dlogpf_dT9;
    get_partition_function(inuc, tfactors.T9, pf_cache, logpf, dlogpf_dT9);

    // Convert logpf to pf and dlogpf_dT9 to dpf_dT
    pf = std::exp(logpf);
    dpf_dT = pf * dlogpf_dT9 * 1.e-9_rt;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_partition_function_cached(const int inuc, const amrex::Real T9,
                                   part_fun::pf_cache_t& pf_cache,
                                   amrex::Real& logpf, amrex::Real& dlogpf_dT9) {
    if (pf_cache.data(inuc, 1) != -900.0_rt) {
        // present in cache
        amrex::ignore_unused(T9);
        logpf = pf_cache.data(inuc, 1);
        dlogpf_dT9 = pf_cache.data(inuc, 2);
    } else {
        get_partition_function(inuc, T9, pf_cache, logpf, dlogpf_dT9);
        pf_cache.data(inuc, 1) = logpf;
        pf_cache.data(inuc, 2) = dlogpf_dT9;
    }
}


template <int spec>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr amrex::Real get_spin_state() {

    static_assert(spec >= 1 && spec <= NumSpec);

    // Set the spin states (2J + 1) of the element

    if constexpr (
        spec == He4 || spec == C12 || spec == O14 ||
        spec == O16
                  )
    {
        return 1.0_rt;
    }
    else if constexpr (
        spec == H1 || spec == He3 || spec == C13 ||
        spec == N13 || spec == N15 || spec == O15
                       )
    {
        return 2.0_rt;
    }
    else if constexpr (
        spec == H2 || spec == N14 || spec == F18
                       )
    {
        return 3.0_rt;
    }
    else if constexpr (
        spec == Be7
                       )
    {
        return 4.0_rt;
    }
    else if constexpr (
        spec == B8
                       )
    {
        return 5.0_rt;
    }
    else if constexpr (
        spec == O17 || spec == F17
                       )
    {
        return 6.0_rt;
    }

    // Return -1 if we don't recognize the species.
    return -1.0_rt;
}


// Legacy (non-templated) interface

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real get_spin_state (int spec) {

    amrex::Real s = -1.0_rt;

    // Set the spin states (2J + 1) of the element

    amrex::constexpr_for<1, NumSpec+1>([&] (auto n) {
        if (n == spec) {
            s = get_spin_state<n>();
        }
    });

    return s;
}


#endif
