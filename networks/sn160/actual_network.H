#ifndef ACTUAL_NETWORK_H
#define ACTUAL_NETWORK_H

#include <random>
#include <string_view>

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Loop.H>

#include <extern_parameters.H>
#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex::literals;


inline
void actual_network_init() {


}

constexpr std::string_view network_name = "pynucastro-cxx";

namespace network
{

    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr amrex::Real bion () {
        using namespace Species;

        static_assert(spec >= 1 && spec <= NumSpec);

        // Set the binding energy of the element

        if constexpr (spec == N) {
            return 0.0_rt;
        }
        else if constexpr (spec == H1) {
            return 0.0_rt;
        }
        else if constexpr (spec == H2) {
            return 2.2245662689997516_rt;
        }
        else if constexpr (spec == He3) {
            return 7.718041347999588_rt;
        }
        else if constexpr (spec == He4) {
            return 28.295662457999697_rt;
        }
        else if constexpr (spec == Li6) {
            return 31.99398709199886_rt;
        }
        else if constexpr (spec == Li7) {
            return 39.245080591999795_rt;
        }
        else if constexpr (spec == Be7) {
            return 37.60083855599987_rt;
        }
        else if constexpr (spec == Be9) {
            return 58.1640247560008_rt;
        }
        else if constexpr (spec == B8) {
            return 37.7372096199997_rt;
        }
        else if constexpr (spec == B10) {
            return 64.750834819999_rt;
        }
        else if constexpr (spec == B11) {
            return 76.20505591999972_rt;
        }
        else if constexpr (spec == C12) {
            return 92.16173498399803_rt;
        }
        else if constexpr (spec == C13) {
            return 97.10804378399916_rt;
        }
        else if constexpr (spec == C14) {
            return 105.28447818399945_rt;
        }
        else if constexpr (spec == N13) {
            return 94.10522604799917_rt;
        }
        else if constexpr (spec == N14) {
            return 104.65860734799753_rt;
        }
        else if constexpr (spec == N15) {
            return 115.49190414799887_rt;
        }
        else if constexpr (spec == O14) {
            return 98.73189611199996_rt;
        }
        else if constexpr (spec == O15) {
            return 111.95539521199862_rt;
        }
        else if constexpr (spec == O16) {
            return 127.6193154119992_rt;
        }
        else if constexpr (spec == O17) {
            return 131.76239561199873_rt;
        }
        else if constexpr (spec == O18) {
            return 139.8077658120019_rt;
        }
        else if constexpr (spec == F17) {
            return 128.21958437599824_rt;
        }
        else if constexpr (spec == F18) {
            return 137.36950247599816_rt;
        }
        else if constexpr (spec == F19) {
            return 147.80136567599766_rt;
        }
        else if constexpr (spec == Ne18) {
            return 132.14265544000227_rt;
        }
        else if constexpr (spec == Ne19) {
            return 143.7795235400008_rt;
        }
        else if constexpr (spec == Ne20) {
            return 160.64482384000075_rt;
        }
        else if constexpr (spec == Ne21) {
            return 167.40598973999658_rt;
        }
        else if constexpr (spec == Ne22) {
            return 177.77024384000106_rt;
        }
        else if constexpr (spec == Na21) {
            return 163.07672270400144_rt;
        }
        else if constexpr (spec == Na22) {
            return 174.14457080400098_rt;
        }
        else if constexpr (spec == Na23) {
            return 186.56435240400242_rt;
        }
        else if constexpr (spec == Mg23) {
            return 181.7258218679999_rt;
        }
        else if constexpr (spec == Mg24) {
            return 198.2570479679962_rt;
        }
        else if constexpr (spec == Mg25) {
            return 205.5875680680001_rt;
        }
        else if constexpr (spec == Mg26) {
            return 216.68065016799665_rt;
        }
        else if constexpr (spec == Al25) {
            return 200.52841103199896_rt;
        }
        else if constexpr (spec == Al26) {
            return 211.89389913199557_rt;
        }
        else if constexpr (spec == Al27) {
            return 224.95193723199915_rt;
        }
        else if constexpr (spec == Si28) {
            return 236.53684539599638_rt;
        }
        else if constexpr (spec == Si29) {
            return 245.01044789599834_rt;
        }
        else if constexpr (spec == Si30) {
            return 255.61964649599756_rt;
        }
        else if constexpr (spec == Si31) {
            return 262.207042596001_rt;
        }
        else if constexpr (spec == Si32) {
            return 271.4070106959989_rt;
        }
        else if constexpr (spec == P29) {
            return 239.28581936000046_rt;
        }
        else if constexpr (spec == P30) {
            return 250.60519745999775_rt;
        }
        else if constexpr (spec == P31) {
            return 262.9161999600037_rt;
        }
        else if constexpr (spec == P32) {
            return 270.8518536599986_rt;
        }
        else if constexpr (spec == P33) {
            return 280.95569175999844_rt;
        }
        else if constexpr (spec == S32) {
            return 271.78016372399725_rt;
        }
        else if constexpr (spec == S33) {
            return 280.42180302400084_rt;
        }
        else if constexpr (spec == S34) {
            return 291.83895282400044_rt;
        }
        else if constexpr (spec == S35) {
            return 298.8247909240017_rt;
        }
        else if constexpr (spec == S36) {
            return 308.714039024002_rt;
        }
        else if constexpr (spec == Cl33) {
            return 274.05689768799857_rt;
        }
        else if constexpr (spec == Cl34) {
            return 285.5650057879975_rt;
        }
        else if constexpr (spec == Cl35) {
            return 298.20976388800045_rt;
        }
        else if constexpr (spec == Cl36) {
            return 306.78956198800006_rt;
        }
        else if constexpr (spec == Cl37) {
            return 317.10042008799786_rt;
        }
        else if constexpr (spec == Ar36) {
            return 306.7167469519991_rt;
        }
        else if constexpr (spec == Ar37) {
            return 315.50420305199077_rt;
        }
        else if constexpr (spec == Ar38) {
            return 327.3426711519933_rt;
        }
        else if constexpr (spec == Ar39) {
            return 333.9411592519973_rt;
        }
        else if constexpr (spec == Ar40) {
            return 343.81037735200516_rt;
        }
        else if constexpr (spec == K37) {
            return 308.5743760159967_rt;
        }
        else if constexpr (spec == K38) {
            return 320.64625411599985_rt;
        }
        else if constexpr (spec == K39) {
            return 333.7240072160057_rt;
        }
        else if constexpr (spec == K40) {
            return 341.52363031599816_rt;
        }
        else if constexpr (spec == K41) {
            return 351.6189974159933_rt;
        }
        else if constexpr (spec == Ca40) {
            return 342.05218528000114_rt;
        }
        else if constexpr (spec == Ca41) {
            return 350.4150113799915_rt;
        }
        else if constexpr (spec == Ca42) {
            return 361.8957094799989_rt;
        }
        else if constexpr (spec == Ca43) {
            return 369.82860757999885_rt;
        }
        else if constexpr (spec == Ca44) {
            return 380.95975567999994_rt;
        }
        else if constexpr (spec == Ca45) {
            return 388.3745737799982_rt;
        }
        else if constexpr (spec == Ca46) {
            return 398.77329187999567_rt;
        }
        else if constexpr (spec == Ca47) {
            return 406.0497099799977_rt;
        }
        else if constexpr (spec == Ca48) {
            return 416.0011960799966_rt;
        }
        else if constexpr (spec == Sc43) {
            return 366.82549054399715_rt;
        }
        else if constexpr (spec == Sc44) {
            return 376.5247086439995_rt;
        }
        else if constexpr (spec == Sc45) {
            return 387.8523267439959_rt;
        }
        else if constexpr (spec == Sc46) {
            return 396.61294484399696_rt;
        }
        else if constexpr (spec == Sc47) {
            return 407.25946294399176_rt;
        }
        else if constexpr (spec == Sc48) {
            return 415.49798104399815_rt;
        }
        else if constexpr (spec == Sc49) {
            return 425.6276991440027_rt;
        }
        else if constexpr (spec == Ti44) {
            return 375.47496160800074_rt;
        }
        else if constexpr (spec == Ti45) {
            return 385.0079797080034_rt;
        }
        else if constexpr (spec == Ti46) {
            return 398.19726780800556_rt;
        }
        else if constexpr (spec == Ti47) {
            return 407.0779259080009_rt;
        }
        else if constexpr (spec == Ti48) {
            return 418.70458400799544_rt;
        }
        else if constexpr (spec == Ti49) {
            return 426.84696210800035_rt;
        }
        else if constexpr (spec == Ti50) {
            return 437.7861402080016_rt;
        }
        else if constexpr (spec == Ti51) {
            return 444.1585883080042_rt;
        }
        else if constexpr (spec == V46) {
            return 390.3625507720062_rt;
        }
        else if constexpr (spec == V47) {
            return 403.3650388719979_rt;
        }
        else if constexpr (spec == V48) {
            return 413.90728697199665_rt;
        }
        else if constexpr (spec == V49) {
            return 425.4628050719984_rt;
        }
        else if constexpr (spec == V50) {
            return 434.79516317199887_rt;
        }
        else if constexpr (spec == V51) {
            return 445.84635127200454_rt;
        }
        else if constexpr (spec == V52) {
            return 453.1575893720001_rt;
        }
        else if constexpr (spec == Cr48) {
            return 411.4679399359957_rt;
        }
        else if constexpr (spec == Cr49) {
            return 422.05065803600155_rt;
        }
        else if constexpr (spec == Cr50) {
            return 435.0509361359945_rt;
        }
        else if constexpr (spec == Cr51) {
            return 444.31160423599795_rt;
        }
        else if constexpr (spec == Cr52) {
            return 456.35172233599593_rt;
        }
        else if constexpr (spec == Cr53) {
            return 464.2911504360018_rt;
        }
        else if constexpr (spec == Cr54) {
            return 474.01022853599716_rt;
        }
        else if constexpr (spec == Mn50) {
            return 426.63411910000286_rt;
        }
        else if constexpr (spec == Mn51) {
            return 440.321747199996_rt;
        }
        else if constexpr (spec == Mn52) {
            return 450.8612552999912_rt;
        }
        else if constexpr (spec == Mn53) {
            return 462.9114833999993_rt;
        }
        else if constexpr (spec == Mn54) {
            return 471.85070149999956_rt;
        }
        else if constexpr (spec == Mn55) {
            return 482.07635959999607_rt;
        }
        else if constexpr (spec == Fe52) {
            return 447.6996182639923_rt;
        }
        else if constexpr (spec == Fe53) {
            return 458.38633636399754_rt;
        }
        else if constexpr (spec == Fe54) {
            return 471.76475446399854_rt;
        }
        else if constexpr (spec == Fe55) {
            return 481.06287256399816_rt;
        }
        else if constexpr (spec == Fe56) {
            return 492.2599506639962_rt;
        }
        else if constexpr (spec == Fe57) {
            return 499.9061287639925_rt;
        }
        else if constexpr (spec == Fe58) {
            return 509.9507268639936_rt;
        }
        else if constexpr (spec == Co53) {
            return 449.3158893279979_rt;
        }
        else if constexpr (spec == Co54) {
            return 462.7379074280034_rt;
        }
        else if constexpr (spec == Co55) {
            return 476.82912552799826_rt;
        }
        else if constexpr (spec == Co56) {
            return 486.91094362799777_rt;
        }
        else if constexpr (spec == Co57) {
            return 498.28746172798856_rt;
        }
        else if constexpr (spec == Co58) {
            return 506.8603798280019_rt;
        }
        else if constexpr (spec == Co59) {
            return 517.3141979280044_rt;
        }
        else if constexpr (spec == Ni56) {
            return 483.9956965919919_rt;
        }
        else if constexpr (spec == Ni57) {
            return 494.2434146919986_rt;
        }
        else if constexpr (spec == Ni58) {
            return 506.4596327920008_rt;
        }
        else if constexpr (spec == Ni59) {
            return 515.4588508919987_rt;
        }
        else if constexpr (spec == Ni60) {
            return 526.8465689919976_rt;
        }
        else if constexpr (spec == Ni61) {
            return 534.6666870919944_rt;
        }
        else if constexpr (spec == Ni62) {
            return 545.2624051920066_rt;
        }
        else if constexpr (spec == Ni63) {
            return 552.1002232920073_rt;
        }
        else if constexpr (spec == Ni64) {
            return 561.7576413920033_rt;
        }
        else if constexpr (spec == Cu57) {
            return 484.68606765599543_rt;
        }
        else if constexpr (spec == Cu58) {
            return 497.11628575600486_rt;
        }
        else if constexpr (spec == Cu59) {
            return 509.8782038560021_rt;
        }
        else if constexpr (spec == Cu60) {
            return 519.9363219559964_rt;
        }
        else if constexpr (spec == Cu61) {
            return 531.6464400559926_rt;
        }
        else if constexpr (spec == Cu62) {
            return 540.5211581559997_rt;
        }
        else if constexpr (spec == Cu63) {
            return 551.384876256001_rt;
        }
        else if constexpr (spec == Cu64) {
            return 559.3006943560031_rt;
        }
        else if constexpr (spec == Cu65) {
            return 569.2113124559983_rt;
        }
        else if constexpr (spec == Zn59) {
            return 499.95305681999895_rt;
        }
        else if constexpr (spec == Zn60) {
            return 514.9831749199948_rt;
        }
        else if constexpr (spec == Zn61) {
            return 525.2289930199986_rt;
        }
        else if constexpr (spec == Zn62) {
            return 538.1194111200093_rt;
        }
        else if constexpr (spec == Zn63) {
            return 547.2360292200028_rt;
        }
        else if constexpr (spec == Zn64) {
            return 559.0979473199986_rt;
        }
        else if constexpr (spec == Zn65) {
            return 567.0772654199973_rt;
        }
        else if constexpr (spec == Zn66) {
            return 578.135783519996_rt;
        }
        else if constexpr (spec == Ga62) {
            return 528.155964084006_rt;
        }
        else if constexpr (spec == Ga63) {
            return 540.787382184004_rt;
        }
        else if constexpr (spec == Ga64) {
            return 551.1444002839926_rt;
        }
        else if constexpr (spec == Ge63) {
            return 530.377935147997_rt;
        }
        else if constexpr (spec == Ge64) {
            return 545.8452532479932_rt;
        }


        // Return zero if we don't recognize the species.
        return 0.0_rt;
    }

    template<int spec>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr amrex::Real mion () {
        using namespace Species;

        static_assert(spec >= 1 && spec <= NumSpec);

        // Set the mass of the element

        if constexpr (spec == N) {
            return 1.674927498034172e-24_rt;
        }
        else if constexpr (spec == H1) {
            return 1.6735328377636005e-24_rt;
        }
        else if constexpr (spec == H2) {
            return 3.3444946862179e-24_rt;
        }
        else if constexpr (spec == He3) {
            return 5.008234515140786e-24_rt;
        }
        else if constexpr (spec == He4) {
            return 6.646479071584587e-24_rt;
        }
        else if constexpr (spec == Li6) {
            return 9.988346544883529e-24_rt;
        }
        else if constexpr (spec == Li7) {
            return 1.1650347794645127e-23_rt;
        }
        else if constexpr (spec == Be7) {
            return 1.1651884262042064e-23_rt;
        }
        else if constexpr (spec == Be9) {
            return 1.4965082049084434e-23_rt;
        }
        else if constexpr (spec == B8) {
            return 1.332517399630266e-23_rt;
        }
        else if constexpr (spec == B10) {
            return 1.6626872831361982e-23_rt;
        }
        else if constexpr (spec == B11) {
            return 1.828138132559934e-23_rt;
        }
        else if constexpr (spec == C12) {
            return 1.99264687992e-23_rt;
        }
        else if constexpr (spec == C13) {
            return 2.1592578700883826e-23_rt;
        }
        else if constexpr (spec == C14) {
            return 2.3252930380658474e-23_rt;
        }
        else if constexpr (spec == N13) {
            return 2.1596537049448796e-23_rt;
        }
        else if constexpr (spec == N14) {
            return 2.3252651436495096e-23_rt;
        }
        else if constexpr (spec == N15) {
            return 2.490826682883827e-23_rt;
        }
        else if constexpr (spec == O14) {
            return 2.3261822098665393e-23_rt;
        }
        else if constexpr (spec == O15) {
            return 2.4913176568383383e-23_rt;
        }
        else if constexpr (spec == O16) {
            return 2.6560180592333686e-23_rt;
        }
        else if constexpr (spec == O17) {
            return 2.822772237905709e-23_rt;
        }
        else if constexpr (spec == O18) {
            return 2.9888307701990364e-23_rt;
        }
        else if constexpr (spec == F17) {
            return 2.8232643353472436e-23_rt;
        }
        else if constexpr (spec == F18) {
            return 2.989125964092377e-23_rt;
        }
        else if constexpr (spec == F19) {
            return 3.1547590653659727e-23_rt;
        }
        else if constexpr (spec == Ne18) {
            return 2.9899182681834417e-23_rt;
        }
        else if constexpr (spec == Ne19) {
            return 3.155336557821978e-23_rt;
        }
        else if constexpr (spec == Ne20) {
            return 3.3198227947612416e-23_rt;
        }
        else if constexpr (spec == Ne21) {
            return 3.4861102572650884e-23_rt;
        }
        else if constexpr (spec == Ne22) {
            return 3.651755410955497e-23_rt;
        }
        else if constexpr (spec == Na21) {
            return 3.486742553187392e-23_rt;
        }
        else if constexpr (spec == Na22) {
            return 3.652262279854593e-23_rt;
        }
        else if constexpr (spec == Na23) {
            return 3.817541002484691e-23_rt;
        }
        else if constexpr (spec == Mg23) {
            return 3.8182640828719474e-23_rt;
        }
        else if constexpr (spec == Mg24) {
            return 3.9828098739467446e-23_rt;
        }
        else if constexpr (spec == Mg25) {
            return 4.148995839845366e-23_rt;
        }
        else if constexpr (spec == Mg26) {
            return 4.314511068143473e-23_rt;
        }
        else if constexpr (spec == Al25) {
            return 4.149758250478667e-23_rt;
        }
        else if constexpr (spec == Al26) {
            return 4.315224917996432e-23_rt;
        }
        else if constexpr (spec == Al27) {
            return 4.480389861070653e-23_rt;
        }
        else if constexpr (spec == Si28) {
            return 4.6456779473820677e-23_rt;
        }
        else if constexpr (spec == Si29) {
            return 4.811660140333913e-23_rt;
        }
        else if constexpr (spec == Si30) {
            return 4.977261628701015e-23_rt;
        }
        else if constexpr (spec == Si31) {
            return 5.143580068485421e-23_rt;
        }
        else if constexpr (spec == Si32) {
            return 5.309432775007637e-23_rt;
        }
        else if constexpr (spec == P29) {
            return 4.812541182037512e-23_rt;
        }
        else if constexpr (spec == P30) {
            return 4.978016069409398e-23_rt;
        }
        else if constexpr (spec == P31) {
            return 5.14331418367544e-23_rt;
        }
        else if constexpr (spec == P32) {
            return 5.30939227471144e-23_rt;
        }
        else if constexpr (spec == P33) {
            return 5.475083851770546e-23_rt;
        }
        else if constexpr (spec == S32) {
            return 5.309087322384128e-23_rt;
        }
        else if constexpr (spec == S33) {
            return 5.475039560055494e-23_rt;
        }
        else if constexpr (spec == S34) {
            return 5.64049701803872e-23_rt;
        }
        else if constexpr (spec == S35) {
            return 5.806744429084988e-23_rt;
        }
        else if constexpr (spec == S36) {
            return 5.97247426028627e-23_rt;
        }
        else if constexpr (spec == Cl33) {
            return 5.476034741466158e-23_rt;
        }
        else if constexpr (spec == Cl34) {
            return 5.641475984659597e-23_rt;
        }
        else if constexpr (spec == Cl35) {
            return 5.806714601585715e-23_rt;
        }
        else if constexpr (spec == Cl36) {
            return 5.972677863452325e-23_rt;
        }
        else if constexpr (spec == Cl37) {
            return 6.13833253584433e-23_rt;
        }
        else if constexpr (spec == Ar36) {
            return 5.972551377884467e-23_rt;
        }
        else if constexpr (spec == Ar37) {
            return 6.138477621350144e-23_rt;
        }
        else if constexpr (spec == Ar38) {
            return 6.303859972524341e-23_rt;
        }
        else if constexpr (spec == Ar39) {
            return 6.470176434980143e-23_rt;
        }
        else if constexpr (spec == Ar40) {
            return 6.635909836853254e-23_rt;
        }
        else if constexpr (spec == K37) {
            return 6.139573509201139e-23_rt;
        }
        else if constexpr (spec == K38) {
            return 6.304914251263421e-23_rt;
        }
        else if constexpr (spec == K39) {
            return 6.470075679819664e-23_rt;
        }
        else if constexpr (spec == K40) {
            return 6.636178020512743e-23_rt;
        }
        else if constexpr (spec == K41) {
            return 6.801871107664764e-23_rt;
        }
        else if constexpr (spec == Ca40) {
            return 6.635944331004904e-23_rt;
        }
        else if constexpr (spec == Ca41) {
            return 6.801946271643761e-23_rt;
        }
        else if constexpr (spec == Ca42) {
            return 6.967392401113526e-23_rt;
        }
        else if constexpr (spec == Ca43) {
            return 7.133470983379846e-23_rt;
        }
        else if constexpr (spec == Ca44) {
            return 7.298979425797081e-23_rt;
        }
        else if constexpr (spec == Ca45) {
            return 7.465150364212236e-23_rt;
        }
        else if constexpr (spec == Ca46) {
            return 7.630789374136599e-23_rt;
        }
        else if constexpr (spec == Ca47) {
            return 7.796984984592748e-23_rt;
        }
        else if constexpr (spec == Ca48) {
            return 7.962703720862762e-23_rt;
        }
        else if constexpr (spec == Sc43) {
            return 7.133866871591414e-23_rt;
        }
        else if constexpr (spec == Sc44) {
            return 7.299630578717192e-23_rt;
        }
        else if constexpr (spec == Sc45) {
            return 7.465103997175656e-23_rt;
        }
        else if constexpr (spec == Sc46) {
            return 7.631035024949397e-23_rt;
        }
        else if constexpr (spec == Sc47) {
            return 7.796629860511342e-23_rt;
        }
        else if constexpr (spec == Sc48) {
            return 7.962653961064011e-23_rt;
        }
        else if constexpr (spec == Sc49) {
            return 8.128340924594064e-23_rt;
        }
        else if constexpr (spec == Ti44) {
            return 7.299678247096977e-23_rt;
        }
        else if constexpr (spec == Ti45) {
            return 7.465471582063894e-23_rt;
        }
        else if constexpr (spec == Ti46) {
            return 7.630613127700392e-23_rt;
        }
        else if constexpr (spec == Ti47) {
            return 7.796522756400429e-23_rt;
        }
        else if constexpr (spec == Ti48) {
            return 7.961942866136786e-23_rt;
        }
        else if constexpr (spec == Ti49) {
            return 8.12798410520117e-23_rt;
        }
        else if constexpr (spec == Ti50) {
            return 8.293526769379317e-23_rt;
        }
        else if constexpr (spec == Ti51) {
            return 8.459883527125195e-23_rt;
        }
        else if constexpr (spec == V46) {
            return 7.631870326846012e-23_rt;
        }
        else if constexpr (spec == V47) {
            return 7.797045172607207e-23_rt;
        }
        else if constexpr (spec == V48) {
            return 7.962658595985007e-23_rt;
        }
        else if constexpr (spec == V49) {
            return 8.128091387578276e-23_rt;
        }
        else if constexpr (spec == V50) {
            return 8.293920493439312e-23_rt;
        }
        else if constexpr (spec == V51) {
            return 8.459443190021272e-23_rt;
        }
        else if constexpr (spec == V52) {
            return 8.625632593248613e-23_rt;
        }
        else if constexpr (spec == Cr48) {
            return 7.962953983065421e-23_rt;
        }
        else if constexpr (spec == Cr49) {
            return 8.128560192010424e-23_rt;
        }
        else if constexpr (spec == Cr50) {
            return 8.293735431739904e-23_rt;
        }
        else if constexpr (spec == Cr51) {
            return 8.459577317504256e-23_rt;
        }
        else if constexpr (spec == Cr52) {
            return 8.624923721300802e-23_rt;
        }
        else if constexpr (spec == Cr53) {
            return 8.791001139488887e-23_rt;
        }
        else if constexpr (spec == Cr54) {
            return 8.956761306248089e-23_rt;
        }
        else if constexpr (spec == Mn50) {
            return 8.295096399635981e-23_rt;
        }
        else if constexpr (spec == Mn51) {
            return 8.460149108098279e-23_rt;
        }
        else if constexpr (spec == Mn52) {
            return 8.625763019925445e-23_rt;
        }
        else if constexpr (spec == Mn53) {
            return 8.79110762145079e-23_rt;
        }
        else if constexpr (spec == Mn54) {
            return 8.957006810882611e-23_rt;
        }
        else if constexpr (spec == Mn55) {
            return 9.122676671554189e-23_rt;
        }
        else if constexpr (spec == Fe52) {
            return 8.626187166893794e-23_rt;
        }
        else if constexpr (spec == Fe53) {
            return 8.791774836154814e-23_rt;
        }
        else if constexpr (spec == Fe54) {
            return 8.95688266630639e-23_rt;
        }
        else if constexpr (spec == Fe55) {
            return 9.122717876001846e-23_rt;
        }
        else if constexpr (spec == Fe56) {
            return 9.288214565329033e-23_rt;
        }
        else if constexpr (spec == Fe57) {
            return 9.454344260077968e-23_rt;
        }
        else if constexpr (spec == Fe58) {
            return 9.620046397626297e-23_rt;
        }
        else if constexpr (spec == Co53) {
            return 8.793252324182073e-23_rt;
        }
        else if constexpr (spec == Co54) {
            return 8.958352381927671e-23_rt;
        }
        else if constexpr (spec == Co55) {
            return 9.123333143937473e-23_rt;
        }
        else if constexpr (spec == Co56) {
            return 9.289028646418131e-23_rt;
        }
        else if constexpr (spec == Co57) {
            return 9.454493347659797e-23_rt;
        }
        else if constexpr (spec == Co58) {
            return 9.620457835997809e-23_rt;
        }
        else if constexpr (spec == Co59) {
            return 9.786087023454982e-23_rt;
        }
        else if constexpr (spec == Ni56) {
            return 9.289408870379396e-23_rt;
        }
        else if constexpr (spec == Ni57) {
            return 9.455074798498772e-23_rt;
        }
        else if constexpr (spec == Ni58) {
            return 9.620389809618879e-23_rt;
        }
        else if constexpr (spec == Ni59) {
            return 9.786278303079172e-23_rt;
        }
        else if constexpr (spec == Ni60) {
            return 9.951741007739486e-23_rt;
        }
        else if constexpr (spec == Ni61) {
            return 1.0117839694866957e-22_rt;
        }
        else if constexpr (spec == Ni62) {
            return 1.0283443586351462e-22_rt;
        }
        else if constexpr (spec == Ni63) {
            return 1.0449717384359494e-22_rt;
        }
        else if constexpr (spec == Ni64) {
            return 1.0615488543012106e-22_rt;
        }
        else if constexpr (spec == Cu57) {
            return 9.456639084334996e-23_rt;
        }
        else if constexpr (spec == Cu58) {
            return 9.62191594648998e-23_rt;
        }
        else if constexpr (spec == Cu59) {
            return 9.787133677749025e-23_rt;
        }
        else if constexpr (spec == Cu60) {
            return 9.952833405138436e-23_rt;
        }
        else if constexpr (spec == Cu61) {
            return 1.0118238636778396e-22_rt;
        }
        else if constexpr (spec == Cu62) {
            return 1.0284149324379613e-22_rt;
        }
        else if constexpr (spec == Cu63) {
            return 1.044970544052462e-22_rt;
        }
        else if constexpr (spec == Cu64) {
            return 1.0615787067577502e-22_rt;
        }
        else if constexpr (spec == Cu65) {
            return 1.0781513089230255e-22_rt;
        }
        else if constexpr (spec == Zn59) {
            return 9.788763529890726e-23_rt;
        }
        else if constexpr (spec == Zn60) {
            return 9.953576917772707e-23_rt;
        }
        else if constexpr (spec == Zn61) {
            return 1.0119243184597849e-22_rt;
        }
        else if constexpr (spec == Zn62) {
            return 1.0284438008651201e-22_rt;
        }
        else if constexpr (spec == Zn63) {
            return 1.0450305573660534e-22_rt;
        }
        else if constexpr (spec == Zn64) {
            return 1.0615683744492522e-22_rt;
        }
        else if constexpr (spec == Zn65) {
            return 1.07817540516422e-22_rt;
        }
        else if constexpr (spec == Zn66) {
            return 1.0947275441532975e-22_rt;
        }
        else if constexpr (spec == Ga62) {
            return 1.0286074688388061e-22_rt;
        }
        else if constexpr (spec == Ga63) {
            return 1.0451315683385183e-22_rt;
        }
        else if constexpr (spec == Ga64) {
            return 1.0616962127009758e-22_rt;
        }
        else if constexpr (spec == Ge63) {
            return 1.0453031869843748e-22_rt;
        }
        else if constexpr (spec == Ge64) {
            return 1.0617767319746517e-22_rt;
        }


        // Return zero if we don't recognize the species.
        return 0.0_rt;
    }

    // Legacy (non-templated) interfaces

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real bion (int spec) {
        using namespace Species;

        amrex::Real b = 0.0_rt;

        // Set the binding energy of the element
        amrex::constexpr_for<1, NumSpec+1>([&] (auto n) {
            if (n == spec) {
                b = bion<n>();
            }
        });

        return b;
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real mion (int spec) {
        using namespace Species;

        amrex::Real m = 0.0_rt;

        amrex::constexpr_for<1, NumSpec+1>([&] (auto n) {
            if (n == spec) {
                m = mion<n>();
            }
        });

        return m;
    }
}

namespace Rates
{

    enum NetworkRates :
    std::uint16_t
    {
        k_Be7_to_Li7_electron_capture_reaclib = 1,
        k_C14_to_N14_beta_neg_reaclib = 2,
        k_N13_to_C13_beta_pos_reaclib = 3,
        k_O14_to_N14_beta_pos_reaclib = 4,
        k_O15_to_N15_beta_pos_reaclib = 5,
        k_F17_to_O17_beta_pos_reaclib = 6,
        k_F18_to_O18_beta_pos_reaclib = 7,
        k_Ne18_to_F18_beta_pos_reaclib = 8,
        k_Ne19_to_F19_beta_pos_reaclib = 9,
        k_Zn59_to_Cu59_beta_pos_reaclib = 10,
        k_Ge63_to_Ga63_beta_pos_reaclib = 11,
        k_B8_to_He4_He4_beta_pos_reaclib = 12,
        k_Zn59_to_p_Ni58_beta_pos_reaclib = 13,
        k_p_n_to_d_reaclib = 14,
        k_p_p_to_d_beta_pos_reaclib = 15,
        k_p_p_to_d_electron_capture_reaclib = 16,
        k_d_p_to_He3_reaclib = 17,
        k_d_d_to_He4_reaclib = 18,
        k_d_He4_to_Li6_reaclib = 19,
        k_He3_n_to_He4_reaclib = 20,
        k_He3_p_to_He4_beta_pos_reaclib = 21,
        k_He3_He4_to_Be7_reaclib = 22,
        k_Li6_n_to_Li7_reaclib = 23,
        k_Li6_p_to_Be7_reaclib = 24,
        k_Li6_He4_to_B10_reaclib = 25,
        k_Li7_He4_to_B11_reaclib = 26,
        k_Be7_p_to_B8_reaclib = 27,
        k_Be9_p_to_B10_reaclib = 28,
        k_B10_n_to_B11_reaclib = 29,
        k_B11_p_to_C12_reaclib = 30,
        k_C12_n_to_C13_reaclib = 31,
        k_C12_p_to_N13_reaclib = 32,
        k_C12_He4_to_O16_reaclib = 33,
        k_C13_n_to_C14_reaclib = 34,
        k_C13_p_to_N14_reaclib = 35,
        k_C14_p_to_N15_reaclib = 36,
        k_C14_He4_to_O18_reaclib = 37,
        k_N13_n_to_N14_reaclib = 38,
        k_N13_p_to_O14_reaclib = 39,
        k_N14_n_to_N15_reaclib = 40,
        k_N14_p_to_O15_reaclib = 41,
        k_N14_He4_to_F18_reaclib = 42,
        k_N15_p_to_O16_reaclib = 43,
        k_N15_He4_to_F19_reaclib = 44,
        k_O14_n_to_O15_reaclib = 45,
        k_O14_He4_to_Ne18_reaclib = 46,
        k_O15_n_to_O16_reaclib = 47,
        k_O15_He4_to_Ne19_reaclib = 48,
        k_O16_n_to_O17_reaclib = 49,
        k_O16_p_to_F17_reaclib = 50,
        k_O16_He4_to_Ne20_reaclib = 51,
        k_O17_n_to_O18_reaclib = 52,
        k_O17_p_to_F18_reaclib = 53,
        k_O17_He4_to_Ne21_reaclib = 54,
        k_O18_p_to_F19_reaclib = 55,
        k_O18_He4_to_Ne22_reaclib = 56,
        k_F17_n_to_F18_reaclib = 57,
        k_F17_p_to_Ne18_reaclib = 58,
        k_F17_He4_to_Na21_reaclib = 59,
        k_F18_n_to_F19_reaclib = 60,
        k_F18_p_to_Ne19_reaclib = 61,
        k_F18_He4_to_Na22_reaclib = 62,
        k_F19_p_to_Ne20_reaclib = 63,
        k_F19_He4_to_Na23_reaclib = 64,
        k_Ne18_n_to_Ne19_reaclib = 65,
        k_Ne19_n_to_Ne20_reaclib = 66,
        k_Ne19_He4_to_Mg23_reaclib = 67,
        k_Ne20_n_to_Ne21_reaclib = 68,
        k_Ne20_p_to_Na21_reaclib = 69,
        k_Ne20_He4_to_Mg24_reaclib = 70,
        k_Ne21_n_to_Ne22_reaclib = 71,
        k_Ne21_p_to_Na22_reaclib = 72,
        k_Ne21_He4_to_Mg25_reaclib = 73,
        k_Ne22_p_to_Na23_reaclib = 74,
        k_Ne22_He4_to_Mg26_reaclib = 75,
        k_Na21_n_to_Na22_reaclib = 76,
        k_Na21_He4_to_Al25_reaclib = 77,
        k_Na22_n_to_Na23_reaclib = 78,
        k_Na22_p_to_Mg23_reaclib = 79,
        k_Na22_He4_to_Al26_reaclib = 80,
        k_Na23_p_to_Mg24_reaclib = 81,
        k_Na23_He4_to_Al27_reaclib = 82,
        k_Mg23_n_to_Mg24_reaclib = 83,
        k_Mg24_n_to_Mg25_reaclib = 84,
        k_Mg24_p_to_Al25_reaclib = 85,
        k_Mg24_He4_to_Si28_reaclib = 86,
        k_Mg25_n_to_Mg26_reaclib = 87,
        k_Mg25_p_to_Al26_reaclib = 88,
        k_Mg25_He4_to_Si29_reaclib = 89,
        k_Mg26_p_to_Al27_reaclib = 90,
        k_Mg26_He4_to_Si30_reaclib = 91,
        k_Al25_n_to_Al26_reaclib = 92,
        k_Al25_He4_to_P29_reaclib = 93,
        k_Al26_n_to_Al27_reaclib = 94,
        k_Al26_He4_to_P30_reaclib = 95,
        k_Al27_p_to_Si28_reaclib = 96,
        k_Al27_He4_to_P31_reaclib = 97,
        k_Si28_n_to_Si29_reaclib = 98,
        k_Si28_p_to_P29_reaclib = 99,
        k_Si28_He4_to_S32_reaclib = 100,
        k_Si29_n_to_Si30_reaclib = 101,
        k_Si29_p_to_P30_reaclib = 102,
        k_Si29_He4_to_S33_reaclib = 103,
        k_Si30_n_to_Si31_reaclib = 104,
        k_Si30_p_to_P31_reaclib = 105,
        k_Si30_He4_to_S34_reaclib = 106,
        k_Si31_n_to_Si32_reaclib = 107,
        k_Si31_p_to_P32_reaclib = 108,
        k_Si31_He4_to_S35_reaclib = 109,
        k_Si32_p_to_P33_reaclib = 110,
        k_Si32_He4_to_S36_reaclib = 111,
        k_P29_n_to_P30_reaclib = 112,
        k_P29_He4_to_Cl33_reaclib = 113,
        k_P30_n_to_P31_reaclib = 114,
        k_P30_He4_to_Cl34_reaclib = 115,
        k_P31_n_to_P32_reaclib = 116,
        k_P31_p_to_S32_reaclib = 117,
        k_P31_He4_to_Cl35_reaclib = 118,
        k_P32_n_to_P33_reaclib = 119,
        k_P32_p_to_S33_reaclib = 120,
        k_P32_He4_to_Cl36_reaclib = 121,
        k_P33_p_to_S34_reaclib = 122,
        k_P33_He4_to_Cl37_reaclib = 123,
        k_S32_n_to_S33_reaclib = 124,
        k_S32_p_to_Cl33_reaclib = 125,
        k_S32_He4_to_Ar36_reaclib = 126,
        k_S33_n_to_S34_reaclib = 127,
        k_S33_p_to_Cl34_reaclib = 128,
        k_S33_He4_to_Ar37_reaclib = 129,
        k_S34_n_to_S35_reaclib = 130,
        k_S34_p_to_Cl35_reaclib = 131,
        k_S34_He4_to_Ar38_reaclib = 132,
        k_S35_n_to_S36_reaclib = 133,
        k_S35_p_to_Cl36_reaclib = 134,
        k_S35_He4_to_Ar39_reaclib = 135,
        k_S36_p_to_Cl37_reaclib = 136,
        k_S36_He4_to_Ar40_reaclib = 137,
        k_Cl33_n_to_Cl34_reaclib = 138,
        k_Cl33_He4_to_K37_reaclib = 139,
        k_Cl34_n_to_Cl35_reaclib = 140,
        k_Cl34_He4_to_K38_reaclib = 141,
        k_Cl35_n_to_Cl36_reaclib = 142,
        k_Cl35_p_to_Ar36_reaclib = 143,
        k_Cl35_He4_to_K39_reaclib = 144,
        k_Cl36_n_to_Cl37_reaclib = 145,
        k_Cl36_p_to_Ar37_reaclib = 146,
        k_Cl36_He4_to_K40_reaclib = 147,
        k_Cl37_p_to_Ar38_reaclib = 148,
        k_Cl37_He4_to_K41_reaclib = 149,
        k_Ar36_n_to_Ar37_reaclib = 150,
        k_Ar36_p_to_K37_reaclib = 151,
        k_Ar36_He4_to_Ca40_reaclib = 152,
        k_Ar37_n_to_Ar38_reaclib = 153,
        k_Ar37_p_to_K38_reaclib = 154,
        k_Ar37_He4_to_Ca41_reaclib = 155,
        k_Ar38_n_to_Ar39_reaclib = 156,
        k_Ar38_p_to_K39_reaclib = 157,
        k_Ar38_He4_to_Ca42_reaclib = 158,
        k_Ar39_n_to_Ar40_reaclib = 159,
        k_Ar39_p_to_K40_reaclib = 160,
        k_Ar39_He4_to_Ca43_reaclib = 161,
        k_Ar40_p_to_K41_reaclib = 162,
        k_Ar40_He4_to_Ca44_reaclib = 163,
        k_K37_n_to_K38_reaclib = 164,
        k_K38_n_to_K39_reaclib = 165,
        k_K39_n_to_K40_reaclib = 166,
        k_K39_p_to_Ca40_reaclib = 167,
        k_K39_He4_to_Sc43_reaclib = 168,
        k_K40_n_to_K41_reaclib = 169,
        k_K40_p_to_Ca41_reaclib = 170,
        k_K40_He4_to_Sc44_reaclib = 171,
        k_K41_p_to_Ca42_reaclib = 172,
        k_K41_He4_to_Sc45_reaclib = 173,
        k_Ca40_n_to_Ca41_reaclib = 174,
        k_Ca40_He4_to_Ti44_reaclib = 175,
        k_Ca41_n_to_Ca42_reaclib = 176,
        k_Ca41_He4_to_Ti45_reaclib = 177,
        k_Ca42_n_to_Ca43_reaclib = 178,
        k_Ca42_p_to_Sc43_reaclib = 179,
        k_Ca42_He4_to_Ti46_reaclib = 180,
        k_Ca43_n_to_Ca44_reaclib = 181,
        k_Ca43_p_to_Sc44_reaclib = 182,
        k_Ca43_He4_to_Ti47_reaclib = 183,
        k_Ca44_n_to_Ca45_reaclib = 184,
        k_Ca44_p_to_Sc45_reaclib = 185,
        k_Ca44_He4_to_Ti48_reaclib = 186,
        k_Ca45_n_to_Ca46_reaclib = 187,
        k_Ca45_p_to_Sc46_reaclib = 188,
        k_Ca45_He4_to_Ti49_reaclib = 189,
        k_Ca46_n_to_Ca47_reaclib = 190,
        k_Ca46_p_to_Sc47_reaclib = 191,
        k_Ca46_He4_to_Ti50_reaclib = 192,
        k_Ca47_n_to_Ca48_reaclib = 193,
        k_Ca47_p_to_Sc48_reaclib = 194,
        k_Ca47_He4_to_Ti51_reaclib = 195,
        k_Ca48_p_to_Sc49_reaclib = 196,
        k_Sc43_n_to_Sc44_reaclib = 197,
        k_Sc43_p_to_Ti44_reaclib = 198,
        k_Sc43_He4_to_V47_reaclib = 199,
        k_Sc44_n_to_Sc45_reaclib = 200,
        k_Sc44_p_to_Ti45_reaclib = 201,
        k_Sc44_He4_to_V48_reaclib = 202,
        k_Sc45_n_to_Sc46_reaclib = 203,
        k_Sc45_p_to_Ti46_reaclib = 204,
        k_Sc45_He4_to_V49_reaclib = 205,
        k_Sc46_n_to_Sc47_reaclib = 206,
        k_Sc46_p_to_Ti47_reaclib = 207,
        k_Sc46_He4_to_V50_reaclib = 208,
        k_Sc47_n_to_Sc48_reaclib = 209,
        k_Sc47_p_to_Ti48_reaclib = 210,
        k_Sc47_He4_to_V51_reaclib = 211,
        k_Sc48_n_to_Sc49_reaclib = 212,
        k_Sc48_p_to_Ti49_reaclib = 213,
        k_Sc48_He4_to_V52_reaclib = 214,
        k_Sc49_p_to_Ti50_reaclib = 215,
        k_Ti44_n_to_Ti45_reaclib = 216,
        k_Ti44_He4_to_Cr48_reaclib = 217,
        k_Ti45_n_to_Ti46_reaclib = 218,
        k_Ti45_p_to_V46_reaclib = 219,
        k_Ti45_He4_to_Cr49_reaclib = 220,
        k_Ti46_n_to_Ti47_reaclib = 221,
        k_Ti46_p_to_V47_reaclib = 222,
        k_Ti46_He4_to_Cr50_reaclib = 223,
        k_Ti47_n_to_Ti48_reaclib = 224,
        k_Ti47_p_to_V48_reaclib = 225,
        k_Ti47_He4_to_Cr51_reaclib = 226,
        k_Ti48_n_to_Ti49_reaclib = 227,
        k_Ti48_p_to_V49_reaclib = 228,
        k_Ti48_He4_to_Cr52_reaclib = 229,
        k_Ti49_n_to_Ti50_reaclib = 230,
        k_Ti49_p_to_V50_reaclib = 231,
        k_Ti49_He4_to_Cr53_reaclib = 232,
        k_Ti50_n_to_Ti51_reaclib = 233,
        k_Ti50_p_to_V51_reaclib = 234,
        k_Ti50_He4_to_Cr54_reaclib = 235,
        k_Ti51_p_to_V52_reaclib = 236,
        k_V46_n_to_V47_reaclib = 237,
        k_V46_He4_to_Mn50_reaclib = 238,
        k_V47_n_to_V48_reaclib = 239,
        k_V47_p_to_Cr48_reaclib = 240,
        k_V47_He4_to_Mn51_reaclib = 241,
        k_V48_n_to_V49_reaclib = 242,
        k_V48_p_to_Cr49_reaclib = 243,
        k_V48_He4_to_Mn52_reaclib = 244,
        k_V49_n_to_V50_reaclib = 245,
        k_V49_p_to_Cr50_reaclib = 246,
        k_V49_He4_to_Mn53_reaclib = 247,
        k_V50_n_to_V51_reaclib = 248,
        k_V50_p_to_Cr51_reaclib = 249,
        k_V50_He4_to_Mn54_reaclib = 250,
        k_V51_n_to_V52_reaclib = 251,
        k_V51_p_to_Cr52_reaclib = 252,
        k_V51_He4_to_Mn55_reaclib = 253,
        k_V52_p_to_Cr53_reaclib = 254,
        k_Cr48_n_to_Cr49_reaclib = 255,
        k_Cr48_He4_to_Fe52_reaclib = 256,
        k_Cr49_n_to_Cr50_reaclib = 257,
        k_Cr49_p_to_Mn50_reaclib = 258,
        k_Cr49_He4_to_Fe53_reaclib = 259,
        k_Cr50_n_to_Cr51_reaclib = 260,
        k_Cr50_p_to_Mn51_reaclib = 261,
        k_Cr50_He4_to_Fe54_reaclib = 262,
        k_Cr51_n_to_Cr52_reaclib = 263,
        k_Cr51_p_to_Mn52_reaclib = 264,
        k_Cr51_He4_to_Fe55_reaclib = 265,
        k_Cr52_n_to_Cr53_reaclib = 266,
        k_Cr52_p_to_Mn53_reaclib = 267,
        k_Cr52_He4_to_Fe56_reaclib = 268,
        k_Cr53_n_to_Cr54_reaclib = 269,
        k_Cr53_p_to_Mn54_reaclib = 270,
        k_Cr53_He4_to_Fe57_reaclib = 271,
        k_Cr54_p_to_Mn55_reaclib = 272,
        k_Cr54_He4_to_Fe58_reaclib = 273,
        k_Mn50_n_to_Mn51_reaclib = 274,
        k_Mn50_He4_to_Co54_reaclib = 275,
        k_Mn51_n_to_Mn52_reaclib = 276,
        k_Mn51_p_to_Fe52_reaclib = 277,
        k_Mn51_He4_to_Co55_reaclib = 278,
        k_Mn52_n_to_Mn53_reaclib = 279,
        k_Mn52_p_to_Fe53_reaclib = 280,
        k_Mn52_He4_to_Co56_reaclib = 281,
        k_Mn53_n_to_Mn54_reaclib = 282,
        k_Mn53_p_to_Fe54_reaclib = 283,
        k_Mn53_He4_to_Co57_reaclib = 284,
        k_Mn54_n_to_Mn55_reaclib = 285,
        k_Mn54_p_to_Fe55_reaclib = 286,
        k_Mn54_He4_to_Co58_reaclib = 287,
        k_Mn55_p_to_Fe56_reaclib = 288,
        k_Mn55_He4_to_Co59_reaclib = 289,
        k_Fe52_n_to_Fe53_reaclib = 290,
        k_Fe52_p_to_Co53_reaclib = 291,
        k_Fe52_He4_to_Ni56_reaclib = 292,
        k_Fe53_n_to_Fe54_reaclib = 293,
        k_Fe53_p_to_Co54_reaclib = 294,
        k_Fe53_He4_to_Ni57_reaclib = 295,
        k_Fe54_n_to_Fe55_reaclib = 296,
        k_Fe54_p_to_Co55_reaclib = 297,
        k_Fe54_He4_to_Ni58_reaclib = 298,
        k_Fe55_n_to_Fe56_reaclib = 299,
        k_Fe55_p_to_Co56_reaclib = 300,
        k_Fe55_He4_to_Ni59_reaclib = 301,
        k_Fe56_n_to_Fe57_reaclib = 302,
        k_Fe56_p_to_Co57_reaclib = 303,
        k_Fe56_He4_to_Ni60_reaclib = 304,
        k_Fe57_n_to_Fe58_reaclib = 305,
        k_Fe57_p_to_Co58_reaclib = 306,
        k_Fe57_He4_to_Ni61_reaclib = 307,
        k_Fe58_p_to_Co59_reaclib = 308,
        k_Fe58_He4_to_Ni62_reaclib = 309,
        k_Co53_n_to_Co54_reaclib = 310,
        k_Co53_He4_to_Cu57_reaclib = 311,
        k_Co54_n_to_Co55_reaclib = 312,
        k_Co54_He4_to_Cu58_reaclib = 313,
        k_Co55_n_to_Co56_reaclib = 314,
        k_Co55_p_to_Ni56_reaclib = 315,
        k_Co55_He4_to_Cu59_reaclib = 316,
        k_Co56_n_to_Co57_reaclib = 317,
        k_Co56_p_to_Ni57_reaclib = 318,
        k_Co56_He4_to_Cu60_reaclib = 319,
        k_Co57_n_to_Co58_reaclib = 320,
        k_Co57_p_to_Ni58_reaclib = 321,
        k_Co57_He4_to_Cu61_reaclib = 322,
        k_Co58_n_to_Co59_reaclib = 323,
        k_Co58_p_to_Ni59_reaclib = 324,
        k_Co58_He4_to_Cu62_reaclib = 325,
        k_Co59_p_to_Ni60_reaclib = 326,
        k_Co59_He4_to_Cu63_reaclib = 327,
        k_Ni56_n_to_Ni57_reaclib = 328,
        k_Ni56_p_to_Cu57_reaclib = 329,
        k_Ni56_He4_to_Zn60_reaclib = 330,
        k_Ni57_n_to_Ni58_reaclib = 331,
        k_Ni57_p_to_Cu58_reaclib = 332,
        k_Ni57_He4_to_Zn61_reaclib = 333,
        k_Ni58_n_to_Ni59_reaclib = 334,
        k_Ni58_p_to_Cu59_reaclib = 335,
        k_Ni58_He4_to_Zn62_reaclib = 336,
        k_Ni59_n_to_Ni60_reaclib = 337,
        k_Ni59_p_to_Cu60_reaclib = 338,
        k_Ni59_He4_to_Zn63_reaclib = 339,
        k_Ni60_n_to_Ni61_reaclib = 340,
        k_Ni60_p_to_Cu61_reaclib = 341,
        k_Ni60_He4_to_Zn64_reaclib = 342,
        k_Ni61_n_to_Ni62_reaclib = 343,
        k_Ni61_p_to_Cu62_reaclib = 344,
        k_Ni61_He4_to_Zn65_reaclib = 345,
        k_Ni62_n_to_Ni63_reaclib = 346,
        k_Ni62_p_to_Cu63_reaclib = 347,
        k_Ni62_He4_to_Zn66_reaclib = 348,
        k_Ni63_n_to_Ni64_reaclib = 349,
        k_Ni63_p_to_Cu64_reaclib = 350,
        k_Ni64_p_to_Cu65_reaclib = 351,
        k_Cu57_n_to_Cu58_reaclib = 352,
        k_Cu58_n_to_Cu59_reaclib = 353,
        k_Cu58_p_to_Zn59_reaclib = 354,
        k_Cu58_He4_to_Ga62_reaclib = 355,
        k_Cu59_n_to_Cu60_reaclib = 356,
        k_Cu59_p_to_Zn60_reaclib = 357,
        k_Cu59_He4_to_Ga63_reaclib = 358,
        k_Cu60_n_to_Cu61_reaclib = 359,
        k_Cu60_p_to_Zn61_reaclib = 360,
        k_Cu60_He4_to_Ga64_reaclib = 361,
        k_Cu61_n_to_Cu62_reaclib = 362,
        k_Cu61_p_to_Zn62_reaclib = 363,
        k_Cu62_n_to_Cu63_reaclib = 364,
        k_Cu62_p_to_Zn63_reaclib = 365,
        k_Cu63_n_to_Cu64_reaclib = 366,
        k_Cu63_p_to_Zn64_reaclib = 367,
        k_Cu64_n_to_Cu65_reaclib = 368,
        k_Cu64_p_to_Zn65_reaclib = 369,
        k_Cu65_p_to_Zn66_reaclib = 370,
        k_Zn59_n_to_Zn60_reaclib = 371,
        k_Zn59_He4_to_Ge63_reaclib = 372,
        k_Zn60_n_to_Zn61_reaclib = 373,
        k_Zn60_He4_to_Ge64_reaclib = 374,
        k_Zn61_n_to_Zn62_reaclib = 375,
        k_Zn61_p_to_Ga62_reaclib = 376,
        k_Zn62_n_to_Zn63_reaclib = 377,
        k_Zn62_p_to_Ga63_reaclib = 378,
        k_Zn63_n_to_Zn64_reaclib = 379,
        k_Zn63_p_to_Ga64_reaclib = 380,
        k_Zn64_n_to_Zn65_reaclib = 381,
        k_Zn65_n_to_Zn66_reaclib = 382,
        k_Ga62_n_to_Ga63_reaclib = 383,
        k_Ga62_p_to_Ge63_reaclib = 384,
        k_Ga63_n_to_Ga64_reaclib = 385,
        k_Ga63_p_to_Ge64_reaclib = 386,
        k_Ge63_n_to_Ge64_reaclib = 387,
        k_d_d_to_n_He3_reaclib = 388,
        k_He3_d_to_p_He4_reaclib = 389,
        k_Li6_p_to_He4_He3_reaclib = 390,
        k_Li6_d_to_n_Be7_reaclib = 391,
        k_Li6_d_to_p_Li7_reaclib = 392,
        k_Li7_p_to_He4_He4_reaclib = 393,
        k_Be7_n_to_p_Li7_reaclib = 394,
        k_Be7_n_to_He4_He4_reaclib = 395,
        k_Be9_p_to_He4_Li6_reaclib = 396,
        k_Be9_He4_to_n_C12_reaclib = 397,
        k_B10_n_to_He4_Li7_reaclib = 398,
        k_B10_p_to_He4_Be7_reaclib = 399,
        k_B10_He4_to_n_N13_reaclib = 400,
        k_B10_He4_to_p_C13_reaclib = 401,
        k_B11_He4_to_n_N14_reaclib = 402,
        k_B11_He4_to_p_C14_reaclib = 403,
        k_C12_C12_to_p_Na23_reaclib = 404,
        k_C12_C12_to_He4_Ne20_reaclib = 405,
        k_C13_p_to_n_N13_reaclib = 406,
        k_C13_d_to_n_N14_reaclib = 407,
        k_C13_He4_to_n_O16_reaclib = 408,
        k_C14_d_to_n_N15_reaclib = 409,
        k_N13_He4_to_p_O16_reaclib = 410,
        k_N14_n_to_p_C14_reaclib = 411,
        k_N14_p_to_n_O14_reaclib = 412,
        k_N14_He4_to_n_F17_reaclib = 413,
        k_N15_p_to_n_O15_reaclib = 414,
        k_N15_p_to_He4_C12_reaclib = 415,
        k_O14_He4_to_p_F17_reaclib = 416,
        k_O15_n_to_He4_C12_reaclib = 417,
        k_O16_C12_to_p_Al27_reaclib = 418,
        k_O16_C12_to_He4_Mg24_reaclib = 419,
        k_O16_O16_to_p_P31_reaclib = 420,
        k_O16_O16_to_He4_Si28_reaclib = 421,
        k_O17_n_to_He4_C14_reaclib = 422,
        k_O17_p_to_He4_N14_reaclib = 423,
        k_O17_He4_to_n_Ne20_reaclib = 424,
        k_O18_p_to_He4_N15_reaclib = 425,
        k_O18_He4_to_n_Ne21_reaclib = 426,
        k_F17_n_to_p_O17_reaclib = 427,
        k_F18_n_to_p_O18_reaclib = 428,
        k_F18_n_to_He4_N15_reaclib = 429,
        k_F18_p_to_He4_O15_reaclib = 430,
        k_F18_He4_to_p_Ne21_reaclib = 431,
        k_F19_p_to_n_Ne19_reaclib = 432,
        k_F19_p_to_He4_O16_reaclib = 433,
        k_F19_He4_to_p_Ne22_reaclib = 434,
        k_Ne18_n_to_p_F18_reaclib = 435,
        k_Ne18_n_to_He4_O15_reaclib = 436,
        k_Ne18_He4_to_p_Na21_reaclib = 437,
        k_Ne19_n_to_He4_O16_reaclib = 438,
        k_Ne19_He4_to_p_Na22_reaclib = 439,
        k_Ne20_p_to_He4_F17_reaclib = 440,
        k_Ne20_C12_to_p_P31_reaclib = 441,
        k_Ne20_C12_to_He4_Si28_reaclib = 442,
        k_Ne21_He4_to_n_Mg24_reaclib = 443,
        k_Ne22_He4_to_n_Mg25_reaclib = 444,
        k_Na21_n_to_p_Ne21_reaclib = 445,
        k_Na21_n_to_He4_F18_reaclib = 446,
        k_Na22_n_to_p_Ne22_reaclib = 447,
        k_Na22_n_to_He4_F19_reaclib = 448,
        k_Na22_He4_to_p_Mg25_reaclib = 449,
        k_Na23_p_to_n_Mg23_reaclib = 450,
        k_Na23_p_to_He4_Ne20_reaclib = 451,
        k_Na23_He4_to_p_Mg26_reaclib = 452,
        k_Mg23_n_to_He4_Ne20_reaclib = 453,
        k_Mg23_n_to_C12_C12_reaclib = 454,
        k_Mg23_He4_to_p_Al26_reaclib = 455,
        k_Mg24_p_to_He4_Na21_reaclib = 456,
        k_Mg25_He4_to_n_Si28_reaclib = 457,
        k_Mg26_He4_to_n_Si29_reaclib = 458,
        k_Al25_n_to_p_Mg25_reaclib = 459,
        k_Al25_n_to_He4_Na22_reaclib = 460,
        k_Al25_He4_to_p_Si28_reaclib = 461,
        k_Al26_n_to_p_Mg26_reaclib = 462,
        k_Al26_n_to_He4_Na23_reaclib = 463,
        k_Al26_He4_to_p_Si29_reaclib = 464,
        k_Al27_p_to_He4_Mg24_reaclib = 465,
        k_Al27_He4_to_n_P30_reaclib = 466,
        k_Al27_He4_to_p_Si30_reaclib = 467,
        k_Si31_p_to_n_P31_reaclib = 468,
        k_Si31_He4_to_n_S34_reaclib = 469,
        k_P29_n_to_p_Si29_reaclib = 470,
        k_P29_n_to_He4_Al26_reaclib = 471,
        k_P29_He4_to_p_S32_reaclib = 472,
        k_P30_n_to_p_Si30_reaclib = 473,
        k_P30_He4_to_p_S33_reaclib = 474,
        k_P31_p_to_He4_Si28_reaclib = 475,
        k_P31_He4_to_p_S34_reaclib = 476,
        k_P32_n_to_p_Si32_reaclib = 477,
        k_P32_p_to_n_S32_reaclib = 478,
        k_P32_p_to_He4_Si29_reaclib = 479,
        k_P33_p_to_He4_Si30_reaclib = 480,
        k_S32_n_to_He4_Si29_reaclib = 481,
        k_S33_n_to_p_P33_reaclib = 482,
        k_S33_n_to_He4_Si30_reaclib = 483,
        k_S34_He4_to_n_Ar37_reaclib = 484,
        k_S34_He4_to_p_Cl37_reaclib = 485,
        k_S35_n_to_He4_Si32_reaclib = 486,
        k_S35_p_to_He4_P32_reaclib = 487,
        k_S35_He4_to_n_Ar38_reaclib = 488,
        k_S36_p_to_He4_P33_reaclib = 489,
        k_Cl33_n_to_p_S33_reaclib = 490,
        k_Cl33_n_to_He4_P30_reaclib = 491,
        k_Cl33_He4_to_p_Ar36_reaclib = 492,
        k_Cl34_n_to_p_S34_reaclib = 493,
        k_Cl34_n_to_He4_P31_reaclib = 494,
        k_Cl34_He4_to_p_Ar37_reaclib = 495,
        k_Cl35_n_to_p_S35_reaclib = 496,
        k_Cl35_n_to_He4_P32_reaclib = 497,
        k_Cl35_p_to_He4_S32_reaclib = 498,
        k_Cl35_He4_to_p_Ar38_reaclib = 499,
        k_Cl36_n_to_p_S36_reaclib = 500,
        k_Cl36_n_to_He4_P33_reaclib = 501,
        k_Cl36_p_to_He4_S33_reaclib = 502,
        k_Ar36_n_to_p_Cl36_reaclib = 503,
        k_Ar36_n_to_He4_S33_reaclib = 504,
        k_Ar37_n_to_p_Cl37_reaclib = 505,
        k_Ar39_n_to_He4_S36_reaclib = 506,
        k_Ar39_p_to_He4_Cl36_reaclib = 507,
        k_Ar40_p_to_He4_Cl37_reaclib = 508,
        k_K37_n_to_p_Ar37_reaclib = 509,
        k_K37_n_to_He4_Cl34_reaclib = 510,
        k_K37_He4_to_p_Ca40_reaclib = 511,
        k_K38_n_to_p_Ar38_reaclib = 512,
        k_K38_n_to_He4_Cl35_reaclib = 513,
        k_K38_He4_to_p_Ca41_reaclib = 514,
        k_K39_n_to_p_Ar39_reaclib = 515,
        k_K39_n_to_He4_Cl36_reaclib = 516,
        k_K39_p_to_He4_Ar36_reaclib = 517,
        k_K40_n_to_p_Ar40_reaclib = 518,
        k_K40_n_to_He4_Cl37_reaclib = 519,
        k_K40_p_to_n_Ca40_reaclib = 520,
        k_K40_p_to_He4_Ar37_reaclib = 521,
        k_K40_He4_to_p_Ca43_reaclib = 522,
        k_K41_p_to_n_Ca41_reaclib = 523,
        k_K41_p_to_He4_Ar38_reaclib = 524,
        k_K41_He4_to_n_Sc44_reaclib = 525,
        k_K41_He4_to_p_Ca44_reaclib = 526,
        k_Ca40_n_to_He4_Ar37_reaclib = 527,
        k_Ca41_n_to_He4_Ar38_reaclib = 528,
        k_Ca42_n_to_He4_Ar39_reaclib = 529,
        k_Ca42_p_to_He4_K39_reaclib = 530,
        k_Ca43_n_to_He4_Ar40_reaclib = 531,
        k_Ca43_He4_to_n_Ti46_reaclib = 532,
        k_Ca45_He4_to_n_Ti48_reaclib = 533,
        k_Ca47_p_to_n_Sc47_reaclib = 534,
        k_Ca47_He4_to_n_Ti50_reaclib = 535,
        k_Sc43_n_to_p_Ca43_reaclib = 536,
        k_Sc43_n_to_He4_K40_reaclib = 537,
        k_Sc43_p_to_He4_Ca40_reaclib = 538,
        k_Sc43_He4_to_p_Ti46_reaclib = 539,
        k_Sc44_n_to_p_Ca44_reaclib = 540,
        k_Sc44_p_to_He4_Ca41_reaclib = 541,
        k_Sc44_He4_to_p_Ti47_reaclib = 542,
        k_Sc45_n_to_p_Ca45_reaclib = 543,
        k_Sc45_p_to_He4_Ca42_reaclib = 544,
        k_Sc45_He4_to_p_Ti48_reaclib = 545,
        k_Sc46_n_to_p_Ca46_reaclib = 546,
        k_Sc46_p_to_n_Ti46_reaclib = 547,
        k_Sc46_p_to_He4_Ca43_reaclib = 548,
        k_Sc46_He4_to_n_V49_reaclib = 549,
        k_Sc46_He4_to_p_Ti49_reaclib = 550,
        k_Sc47_p_to_He4_Ca44_reaclib = 551,
        k_Sc47_He4_to_p_Ti50_reaclib = 552,
        k_Sc48_n_to_p_Ca48_reaclib = 553,
        k_Sc48_p_to_n_Ti48_reaclib = 554,
        k_Sc48_p_to_He4_Ca45_reaclib = 555,
        k_Sc48_He4_to_n_V51_reaclib = 556,
        k_Sc48_He4_to_p_Ti51_reaclib = 557,
        k_Sc49_p_to_n_Ti49_reaclib = 558,
        k_Sc49_p_to_He4_Ca46_reaclib = 559,
        k_Ti44_n_to_p_Sc44_reaclib = 560,
        k_Ti44_n_to_He4_Ca41_reaclib = 561,
        k_Ti44_He4_to_p_V47_reaclib = 562,
        k_Ti45_n_to_p_Sc45_reaclib = 563,
        k_Ti45_n_to_He4_Ca42_reaclib = 564,
        k_Ti45_He4_to_p_V48_reaclib = 565,
        k_Ti47_n_to_p_Sc47_reaclib = 566,
        k_Ti47_n_to_He4_Ca44_reaclib = 567,
        k_Ti49_n_to_He4_Ca46_reaclib = 568,
        k_Ti49_He4_to_n_Cr52_reaclib = 569,
        k_Ti51_n_to_He4_Ca48_reaclib = 570,
        k_Ti51_p_to_n_V51_reaclib = 571,
        k_Ti51_He4_to_n_Cr54_reaclib = 572,
        k_V46_n_to_p_Ti46_reaclib = 573,
        k_V46_n_to_He4_Sc43_reaclib = 574,
        k_V46_He4_to_p_Cr49_reaclib = 575,
        k_V47_n_to_p_Ti47_reaclib = 576,
        k_V47_n_to_He4_Sc44_reaclib = 577,
        k_V47_He4_to_p_Cr50_reaclib = 578,
        k_V48_n_to_p_Ti48_reaclib = 579,
        k_V48_n_to_He4_Sc45_reaclib = 580,
        k_V48_He4_to_p_Cr51_reaclib = 581,
        k_V49_n_to_p_Ti49_reaclib = 582,
        k_V49_p_to_He4_Ti46_reaclib = 583,
        k_V49_He4_to_p_Cr52_reaclib = 584,
        k_V50_n_to_p_Ti50_reaclib = 585,
        k_V50_n_to_He4_Sc47_reaclib = 586,
        k_V50_p_to_n_Cr50_reaclib = 587,
        k_V50_p_to_He4_Ti47_reaclib = 588,
        k_V50_He4_to_p_Cr53_reaclib = 589,
        k_V51_p_to_He4_Ti48_reaclib = 590,
        k_V52_n_to_He4_Sc49_reaclib = 591,
        k_V52_p_to_n_Cr52_reaclib = 592,
        k_V52_p_to_He4_Ti49_reaclib = 593,
        k_V52_He4_to_n_Mn55_reaclib = 594,
        k_Cr48_n_to_p_V48_reaclib = 595,
        k_Cr48_n_to_He4_Ti45_reaclib = 596,
        k_Cr48_He4_to_p_Mn51_reaclib = 597,
        k_Cr49_n_to_p_V49_reaclib = 598,
        k_Cr49_n_to_He4_Ti46_reaclib = 599,
        k_Cr49_He4_to_p_Mn52_reaclib = 600,
        k_Cr50_n_to_He4_Ti47_reaclib = 601,
        k_Cr51_n_to_p_V51_reaclib = 602,
        k_Cr51_n_to_He4_Ti48_reaclib = 603,
        k_Cr53_n_to_He4_Ti50_reaclib = 604,
        k_Cr54_p_to_He4_V51_reaclib = 605,
        k_Mn50_n_to_p_Cr50_reaclib = 606,
        k_Mn50_n_to_He4_V47_reaclib = 607,
        k_Mn50_He4_to_p_Fe53_reaclib = 608,
        k_Mn51_n_to_p_Cr51_reaclib = 609,
        k_Mn51_n_to_He4_V48_reaclib = 610,
        k_Mn51_He4_to_p_Fe54_reaclib = 611,
        k_Mn52_n_to_p_Cr52_reaclib = 612,
        k_Mn52_n_to_He4_V49_reaclib = 613,
        k_Mn52_He4_to_p_Fe55_reaclib = 614,
        k_Mn53_n_to_p_Cr53_reaclib = 615,
        k_Mn53_n_to_He4_V50_reaclib = 616,
        k_Mn53_p_to_He4_Cr50_reaclib = 617,
        k_Mn53_He4_to_p_Fe56_reaclib = 618,
        k_Mn54_n_to_p_Cr54_reaclib = 619,
        k_Mn54_n_to_He4_V51_reaclib = 620,
        k_Mn54_p_to_He4_Cr51_reaclib = 621,
        k_Mn55_p_to_He4_Cr52_reaclib = 622,
        k_Fe52_n_to_p_Mn52_reaclib = 623,
        k_Fe52_n_to_He4_Cr49_reaclib = 624,
        k_Fe52_He4_to_p_Co55_reaclib = 625,
        k_Fe53_n_to_p_Mn53_reaclib = 626,
        k_Fe53_n_to_He4_Cr50_reaclib = 627,
        k_Fe53_He4_to_p_Co56_reaclib = 628,
        k_Fe54_n_to_p_Mn54_reaclib = 629,
        k_Fe54_n_to_He4_Cr51_reaclib = 630,
        k_Fe55_n_to_p_Mn55_reaclib = 631,
        k_Fe55_n_to_He4_Cr52_reaclib = 632,
        k_Fe56_n_to_He4_Cr53_reaclib = 633,
        k_Fe57_n_to_He4_Cr54_reaclib = 634,
        k_Fe57_p_to_He4_Mn54_reaclib = 635,
        k_Fe58_p_to_He4_Mn55_reaclib = 636,
        k_Co53_n_to_p_Fe53_reaclib = 637,
        k_Co53_n_to_He4_Mn50_reaclib = 638,
        k_Co53_He4_to_p_Ni56_reaclib = 639,
        k_Co54_n_to_p_Fe54_reaclib = 640,
        k_Co54_n_to_He4_Mn51_reaclib = 641,
        k_Co54_He4_to_p_Ni57_reaclib = 642,
        k_Co55_n_to_p_Fe55_reaclib = 643,
        k_Co55_n_to_He4_Mn52_reaclib = 644,
        k_Co55_He4_to_p_Ni58_reaclib = 645,
        k_Co56_n_to_p_Fe56_reaclib = 646,
        k_Co56_n_to_He4_Mn53_reaclib = 647,
        k_Co56_He4_to_p_Ni59_reaclib = 648,
        k_Co57_n_to_p_Fe57_reaclib = 649,
        k_Co57_n_to_He4_Mn54_reaclib = 650,
        k_Co57_p_to_He4_Fe54_reaclib = 651,
        k_Co57_He4_to_p_Ni60_reaclib = 652,
        k_Co58_n_to_p_Fe58_reaclib = 653,
        k_Co58_n_to_He4_Mn55_reaclib = 654,
        k_Co58_p_to_He4_Fe55_reaclib = 655,
        k_Co59_p_to_He4_Fe56_reaclib = 656,
        k_Ni56_n_to_p_Co56_reaclib = 657,
        k_Ni56_n_to_He4_Fe53_reaclib = 658,
        k_Ni57_n_to_p_Co57_reaclib = 659,
        k_Ni57_n_to_He4_Fe54_reaclib = 660,
        k_Ni58_n_to_p_Co58_reaclib = 661,
        k_Ni58_n_to_He4_Fe55_reaclib = 662,
        k_Ni59_n_to_p_Co59_reaclib = 663,
        k_Ni59_n_to_He4_Fe56_reaclib = 664,
        k_Ni60_n_to_He4_Fe57_reaclib = 665,
        k_Ni61_n_to_He4_Fe58_reaclib = 666,
        k_Ni61_p_to_He4_Co58_reaclib = 667,
        k_Ni62_p_to_He4_Co59_reaclib = 668,
        k_Cu57_n_to_p_Ni57_reaclib = 669,
        k_Cu57_n_to_He4_Co54_reaclib = 670,
        k_Cu57_He4_to_p_Zn60_reaclib = 671,
        k_Cu58_n_to_p_Ni58_reaclib = 672,
        k_Cu58_n_to_He4_Co55_reaclib = 673,
        k_Cu59_n_to_p_Ni59_reaclib = 674,
        k_Cu59_n_to_He4_Co56_reaclib = 675,
        k_Cu59_p_to_He4_Ni56_reaclib = 676,
        k_Cu60_n_to_p_Ni60_reaclib = 677,
        k_Cu60_n_to_He4_Co57_reaclib = 678,
        k_Cu60_p_to_He4_Ni57_reaclib = 679,
        k_Cu61_n_to_p_Ni61_reaclib = 680,
        k_Cu61_n_to_He4_Co58_reaclib = 681,
        k_Cu61_p_to_He4_Ni58_reaclib = 682,
        k_Cu62_n_to_p_Ni62_reaclib = 683,
        k_Cu62_n_to_He4_Co59_reaclib = 684,
        k_Cu62_p_to_He4_Ni59_reaclib = 685,
        k_Cu63_n_to_p_Ni63_reaclib = 686,
        k_Cu63_p_to_He4_Ni60_reaclib = 687,
        k_Cu64_n_to_p_Ni64_reaclib = 688,
        k_Cu64_p_to_He4_Ni61_reaclib = 689,
        k_Cu65_p_to_He4_Ni62_reaclib = 690,
        k_Zn59_n_to_p_Cu59_reaclib = 691,
        k_Zn59_n_to_He4_Ni56_reaclib = 692,
        k_Zn60_n_to_p_Cu60_reaclib = 693,
        k_Zn60_n_to_He4_Ni57_reaclib = 694,
        k_Zn61_n_to_p_Cu61_reaclib = 695,
        k_Zn61_n_to_He4_Ni58_reaclib = 696,
        k_Zn61_p_to_He4_Cu58_reaclib = 697,
        k_Zn62_n_to_p_Cu62_reaclib = 698,
        k_Zn62_n_to_He4_Ni59_reaclib = 699,
        k_Zn62_p_to_He4_Cu59_reaclib = 700,
        k_Zn63_n_to_p_Cu63_reaclib = 701,
        k_Zn63_n_to_He4_Ni60_reaclib = 702,
        k_Zn63_p_to_He4_Cu60_reaclib = 703,
        k_Zn64_n_to_p_Cu64_reaclib = 704,
        k_Zn64_n_to_He4_Ni61_reaclib = 705,
        k_Zn64_p_to_He4_Cu61_reaclib = 706,
        k_Zn65_n_to_p_Cu65_reaclib = 707,
        k_Zn65_n_to_He4_Ni62_reaclib = 708,
        k_Zn65_p_to_He4_Cu62_reaclib = 709,
        k_Zn66_n_to_He4_Ni63_reaclib = 710,
        k_Zn66_p_to_He4_Cu63_reaclib = 711,
        k_Ga62_n_to_p_Zn62_reaclib = 712,
        k_Ga62_n_to_He4_Cu59_reaclib = 713,
        k_Ga62_p_to_He4_Zn59_reaclib = 714,
        k_Ga63_n_to_p_Zn63_reaclib = 715,
        k_Ga63_n_to_He4_Cu60_reaclib = 716,
        k_Ga63_p_to_He4_Zn60_reaclib = 717,
        k_Ga64_n_to_p_Zn64_reaclib = 718,
        k_Ga64_n_to_He4_Cu61_reaclib = 719,
        k_Ga64_p_to_He4_Zn61_reaclib = 720,
        k_Ge63_n_to_p_Ga63_reaclib = 721,
        k_Ge63_n_to_He4_Zn60_reaclib = 722,
        k_Ge64_n_to_p_Ga64_reaclib = 723,
        k_Ge64_n_to_He4_Zn61_reaclib = 724,
        k_He3_He3_to_p_p_He4_reaclib = 725,
        k_Li7_d_to_n_He4_He4_reaclib = 726,
        k_Be7_d_to_p_He4_He4_reaclib = 727,
        k_Be9_p_to_d_He4_He4_reaclib = 728,
        k_B8_n_to_p_He4_He4_reaclib = 729,
        k_B11_p_to_3He4_reaclib = 730,
        k_Li7_He3_to_n_p_He4_He4_reaclib = 731,
        k_Be7_He3_to_p_p_He4_He4_reaclib = 732,
        k_Be9_p_to_n_p_He4_He4_reaclib = 733,
        k_He4_p_n_to_Li6_reaclib = 734,
        k_He4_He4_n_to_Be9_reaclib = 735,
        k_He4_He4_He4_to_C12_reaclib = 736,
        k_p_p_n_to_p_d_reaclib = 737,
        k_Na21_to_Ne21_electron_capture_weaktab = 738,
        k_Ne21_to_Na21_beta_neg_weaktab = 739,
        k_Na22_to_Ne22_electron_capture_weaktab = 740,
        k_Ne22_to_Na22_beta_neg_weaktab = 741,
        k_Mg23_to_Na23_electron_capture_weaktab = 742,
        k_Na23_to_Mg23_beta_neg_weaktab = 743,
        k_Al25_to_Mg25_electron_capture_weaktab = 744,
        k_Mg25_to_Al25_beta_neg_weaktab = 745,
        k_Al26_to_Mg26_electron_capture_weaktab = 746,
        k_Mg26_to_Al26_beta_neg_weaktab = 747,
        k_P29_to_Si29_electron_capture_weaktab = 748,
        k_Si29_to_P29_beta_neg_weaktab = 749,
        k_P30_to_Si30_electron_capture_weaktab = 750,
        k_Si30_to_P30_beta_neg_weaktab = 751,
        k_P31_to_Si31_electron_capture_weaktab = 752,
        k_Si31_to_P31_beta_neg_weaktab = 753,
        k_P32_to_S32_beta_neg_weaktab = 754,
        k_P32_to_Si32_electron_capture_weaktab = 755,
        k_S32_to_P32_electron_capture_weaktab = 756,
        k_Si32_to_P32_beta_neg_weaktab = 757,
        k_Cl33_to_S33_electron_capture_weaktab = 758,
        k_P33_to_S33_beta_neg_weaktab = 759,
        k_S33_to_Cl33_beta_neg_weaktab = 760,
        k_S33_to_P33_electron_capture_weaktab = 761,
        k_Cl34_to_S34_electron_capture_weaktab = 762,
        k_S34_to_Cl34_beta_neg_weaktab = 763,
        k_Cl35_to_S35_electron_capture_weaktab = 764,
        k_S35_to_Cl35_beta_neg_weaktab = 765,
        k_Ar36_to_Cl36_electron_capture_weaktab = 766,
        k_Cl36_to_Ar36_beta_neg_weaktab = 767,
        k_Cl36_to_S36_electron_capture_weaktab = 768,
        k_S36_to_Cl36_beta_neg_weaktab = 769,
        k_Ar37_to_Cl37_electron_capture_weaktab = 770,
        k_Ar37_to_K37_beta_neg_weaktab = 771,
        k_Cl37_to_Ar37_beta_neg_weaktab = 772,
        k_K37_to_Ar37_electron_capture_weaktab = 773,
        k_Ar38_to_K38_beta_neg_weaktab = 774,
        k_K38_to_Ar38_electron_capture_weaktab = 775,
        k_Ar39_to_K39_beta_neg_weaktab = 776,
        k_K39_to_Ar39_electron_capture_weaktab = 777,
        k_Ar40_to_K40_beta_neg_weaktab = 778,
        k_Ca40_to_K40_electron_capture_weaktab = 779,
        k_K40_to_Ar40_electron_capture_weaktab = 780,
        k_K40_to_Ca40_beta_neg_weaktab = 781,
        k_Ca41_to_K41_electron_capture_weaktab = 782,
        k_K41_to_Ca41_beta_neg_weaktab = 783,
        k_Ca43_to_Sc43_beta_neg_weaktab = 784,
        k_Sc43_to_Ca43_electron_capture_weaktab = 785,
        k_Ca44_to_Sc44_beta_neg_weaktab = 786,
        k_Sc44_to_Ca44_electron_capture_weaktab = 787,
        k_Sc44_to_Ti44_beta_neg_weaktab = 788,
        k_Ti44_to_Sc44_electron_capture_weaktab = 789,
        k_Co53_to_Fe53_electron_capture_weaktab = 790,
        k_Fe53_to_Co53_beta_neg_weaktab = 791,
        k_Cu57_to_Ni57_electron_capture_weaktab = 792,
        k_Ni57_to_Cu57_beta_neg_weaktab = 793,
        k_Ca45_to_Sc45_beta_neg_weaktab = 794,
        k_Sc45_to_Ca45_electron_capture_weaktab = 795,
        k_Sc45_to_Ti45_beta_neg_weaktab = 796,
        k_Ti45_to_Sc45_electron_capture_weaktab = 797,
        k_Ca46_to_Sc46_beta_neg_weaktab = 798,
        k_Sc46_to_Ca46_electron_capture_weaktab = 799,
        k_Sc46_to_Ti46_beta_neg_weaktab = 800,
        k_Ti46_to_Sc46_electron_capture_weaktab = 801,
        k_Ti46_to_V46_beta_neg_weaktab = 802,
        k_V46_to_Ti46_electron_capture_weaktab = 803,
        k_Ca47_to_Sc47_beta_neg_weaktab = 804,
        k_Sc47_to_Ca47_electron_capture_weaktab = 805,
        k_Sc47_to_Ti47_beta_neg_weaktab = 806,
        k_Ti47_to_Sc47_electron_capture_weaktab = 807,
        k_Ti47_to_V47_beta_neg_weaktab = 808,
        k_V47_to_Ti47_electron_capture_weaktab = 809,
        k_Ca48_to_Sc48_beta_neg_weaktab = 810,
        k_Cr48_to_V48_electron_capture_weaktab = 811,
        k_Sc48_to_Ca48_electron_capture_weaktab = 812,
        k_Sc48_to_Ti48_beta_neg_weaktab = 813,
        k_Ti48_to_Sc48_electron_capture_weaktab = 814,
        k_Ti48_to_V48_beta_neg_weaktab = 815,
        k_V48_to_Cr48_beta_neg_weaktab = 816,
        k_V48_to_Ti48_electron_capture_weaktab = 817,
        k_Cr49_to_V49_electron_capture_weaktab = 818,
        k_Sc49_to_Ti49_beta_neg_weaktab = 819,
        k_Ti49_to_Sc49_electron_capture_weaktab = 820,
        k_Ti49_to_V49_beta_neg_weaktab = 821,
        k_V49_to_Cr49_beta_neg_weaktab = 822,
        k_V49_to_Ti49_electron_capture_weaktab = 823,
        k_Cr50_to_Mn50_beta_neg_weaktab = 824,
        k_Cr50_to_V50_electron_capture_weaktab = 825,
        k_Mn50_to_Cr50_electron_capture_weaktab = 826,
        k_Ti50_to_V50_beta_neg_weaktab = 827,
        k_V50_to_Cr50_beta_neg_weaktab = 828,
        k_V50_to_Ti50_electron_capture_weaktab = 829,
        k_Cr51_to_Mn51_beta_neg_weaktab = 830,
        k_Cr51_to_V51_electron_capture_weaktab = 831,
        k_Mn51_to_Cr51_electron_capture_weaktab = 832,
        k_Ti51_to_V51_beta_neg_weaktab = 833,
        k_V51_to_Cr51_beta_neg_weaktab = 834,
        k_V51_to_Ti51_electron_capture_weaktab = 835,
        k_Cr52_to_Mn52_beta_neg_weaktab = 836,
        k_Cr52_to_V52_electron_capture_weaktab = 837,
        k_Fe52_to_Mn52_electron_capture_weaktab = 838,
        k_Mn52_to_Cr52_electron_capture_weaktab = 839,
        k_Mn52_to_Fe52_beta_neg_weaktab = 840,
        k_V52_to_Cr52_beta_neg_weaktab = 841,
        k_Cr53_to_Mn53_beta_neg_weaktab = 842,
        k_Fe53_to_Mn53_electron_capture_weaktab = 843,
        k_Mn53_to_Cr53_electron_capture_weaktab = 844,
        k_Mn53_to_Fe53_beta_neg_weaktab = 845,
        k_Co54_to_Fe54_electron_capture_weaktab = 846,
        k_Cr54_to_Mn54_beta_neg_weaktab = 847,
        k_Fe54_to_Co54_beta_neg_weaktab = 848,
        k_Fe54_to_Mn54_electron_capture_weaktab = 849,
        k_Mn54_to_Cr54_electron_capture_weaktab = 850,
        k_Mn54_to_Fe54_beta_neg_weaktab = 851,
        k_Co55_to_Fe55_electron_capture_weaktab = 852,
        k_Fe55_to_Co55_beta_neg_weaktab = 853,
        k_Fe55_to_Mn55_electron_capture_weaktab = 854,
        k_Mn55_to_Fe55_beta_neg_weaktab = 855,
        k_Co56_to_Fe56_electron_capture_weaktab = 856,
        k_Co56_to_Ni56_beta_neg_weaktab = 857,
        k_Fe56_to_Co56_beta_neg_weaktab = 858,
        k_Ni56_to_Co56_electron_capture_weaktab = 859,
        k_Co57_to_Fe57_electron_capture_weaktab = 860,
        k_Co57_to_Ni57_beta_neg_weaktab = 861,
        k_Fe57_to_Co57_beta_neg_weaktab = 862,
        k_Ni57_to_Co57_electron_capture_weaktab = 863,
        k_Co58_to_Fe58_electron_capture_weaktab = 864,
        k_Co58_to_Ni58_beta_neg_weaktab = 865,
        k_Cu58_to_Ni58_electron_capture_weaktab = 866,
        k_Fe58_to_Co58_beta_neg_weaktab = 867,
        k_Ni58_to_Co58_electron_capture_weaktab = 868,
        k_Ni58_to_Cu58_beta_neg_weaktab = 869,
        k_Co59_to_Ni59_beta_neg_weaktab = 870,
        k_Cu59_to_Ni59_electron_capture_weaktab = 871,
        k_Ni59_to_Co59_electron_capture_weaktab = 872,
        k_Ni59_to_Cu59_beta_neg_weaktab = 873,
        k_Cu60_to_Ni60_electron_capture_weaktab = 874,
        k_Cu60_to_Zn60_beta_neg_weaktab = 875,
        k_Ni60_to_Cu60_beta_neg_weaktab = 876,
        k_Zn60_to_Cu60_electron_capture_weaktab = 877,
        k_Cu61_to_Ni61_electron_capture_weaktab = 878,
        k_Cu61_to_Zn61_beta_neg_weaktab = 879,
        k_Ni61_to_Cu61_beta_neg_weaktab = 880,
        k_Zn61_to_Cu61_electron_capture_weaktab = 881,
        k_Cu62_to_Ni62_electron_capture_weaktab = 882,
        k_Cu62_to_Zn62_beta_neg_weaktab = 883,
        k_Ga62_to_Zn62_electron_capture_weaktab = 884,
        k_Ni62_to_Cu62_beta_neg_weaktab = 885,
        k_Zn62_to_Cu62_electron_capture_weaktab = 886,
        k_Zn62_to_Ga62_beta_neg_weaktab = 887,
        k_Cu63_to_Ni63_electron_capture_weaktab = 888,
        k_Cu63_to_Zn63_beta_neg_weaktab = 889,
        k_Ga63_to_Zn63_electron_capture_weaktab = 890,
        k_Ni63_to_Cu63_beta_neg_weaktab = 891,
        k_Zn63_to_Cu63_electron_capture_weaktab = 892,
        k_Zn63_to_Ga63_beta_neg_weaktab = 893,
        k_Cu64_to_Ni64_electron_capture_weaktab = 894,
        k_Cu64_to_Zn64_beta_neg_weaktab = 895,
        k_Ga64_to_Ge64_beta_neg_weaktab = 896,
        k_Ga64_to_Zn64_electron_capture_weaktab = 897,
        k_Ge64_to_Ga64_electron_capture_weaktab = 898,
        k_Ni64_to_Cu64_beta_neg_weaktab = 899,
        k_Zn64_to_Cu64_electron_capture_weaktab = 900,
        k_Zn64_to_Ga64_beta_neg_weaktab = 901,
        k_Cu65_to_Zn65_beta_neg_weaktab = 902,
        k_Zn65_to_Cu65_electron_capture_weaktab = 903,
        k_n_to_p_beta_neg_weaktab = 904,
        k_p_to_n_electron_capture_weaktab = 905,
        k_d_to_n_p_derived = 906,
        k_He3_to_p_d_derived = 907,
        k_He4_to_n_He3_derived = 908,
        k_He4_to_d_d_derived = 909,
        k_Li6_to_He4_d_derived = 910,
        k_Li7_to_n_Li6_derived = 911,
        k_Be7_to_p_Li6_derived = 912,
        k_Be7_to_He4_He3_derived = 913,
        k_B8_to_p_Be7_derived = 914,
        k_B10_to_p_Be9_derived = 915,
        k_B10_to_He4_Li6_derived = 916,
        k_B11_to_n_B10_derived = 917,
        k_B11_to_He4_Li7_derived = 918,
        k_C12_to_p_B11_derived = 919,
        k_C13_to_n_C12_derived = 920,
        k_C14_to_n_C13_derived = 921,
        k_N13_to_p_C12_derived = 922,
        k_N14_to_n_N13_derived = 923,
        k_N14_to_p_C13_derived = 924,
        k_N15_to_n_N14_derived = 925,
        k_N15_to_p_C14_derived = 926,
        k_O14_to_p_N13_derived = 927,
        k_O15_to_n_O14_derived = 928,
        k_O15_to_p_N14_derived = 929,
        k_O16_to_n_O15_derived = 930,
        k_O16_to_p_N15_derived = 931,
        k_O16_to_He4_C12_derived = 932,
        k_O17_to_n_O16_derived = 933,
        k_O18_to_n_O17_derived = 934,
        k_O18_to_He4_C14_derived = 935,
        k_F17_to_p_O16_derived = 936,
        k_F18_to_n_F17_derived = 937,
        k_F18_to_p_O17_derived = 938,
        k_F18_to_He4_N14_derived = 939,
        k_F19_to_n_F18_derived = 940,
        k_F19_to_p_O18_derived = 941,
        k_F19_to_He4_N15_derived = 942,
        k_Ne18_to_p_F17_derived = 943,
        k_Ne18_to_He4_O14_derived = 944,
        k_Ne19_to_n_Ne18_derived = 945,
        k_Ne19_to_p_F18_derived = 946,
        k_Ne19_to_He4_O15_derived = 947,
        k_Ne20_to_n_Ne19_derived = 948,
        k_Ne20_to_p_F19_derived = 949,
        k_Ne20_to_He4_O16_derived = 950,
        k_Ne21_to_n_Ne20_derived = 951,
        k_Ne21_to_He4_O17_derived = 952,
        k_Ne22_to_n_Ne21_derived = 953,
        k_Ne22_to_He4_O18_derived = 954,
        k_Na21_to_p_Ne20_derived = 955,
        k_Na21_to_He4_F17_derived = 956,
        k_Na22_to_n_Na21_derived = 957,
        k_Na22_to_p_Ne21_derived = 958,
        k_Na22_to_He4_F18_derived = 959,
        k_Na23_to_n_Na22_derived = 960,
        k_Na23_to_p_Ne22_derived = 961,
        k_Na23_to_He4_F19_derived = 962,
        k_Mg23_to_p_Na22_derived = 963,
        k_Mg23_to_He4_Ne19_derived = 964,
        k_Mg24_to_n_Mg23_derived = 965,
        k_Mg24_to_p_Na23_derived = 966,
        k_Mg24_to_He4_Ne20_derived = 967,
        k_Mg25_to_n_Mg24_derived = 968,
        k_Mg25_to_He4_Ne21_derived = 969,
        k_Mg26_to_n_Mg25_derived = 970,
        k_Mg26_to_He4_Ne22_derived = 971,
        k_Al25_to_p_Mg24_derived = 972,
        k_Al25_to_He4_Na21_derived = 973,
        k_Al26_to_n_Al25_derived = 974,
        k_Al26_to_p_Mg25_derived = 975,
        k_Al26_to_He4_Na22_derived = 976,
        k_Al27_to_n_Al26_derived = 977,
        k_Al27_to_p_Mg26_derived = 978,
        k_Al27_to_He4_Na23_derived = 979,
        k_Si28_to_p_Al27_derived = 980,
        k_Si28_to_He4_Mg24_derived = 981,
        k_Si29_to_n_Si28_derived = 982,
        k_Si29_to_He4_Mg25_derived = 983,
        k_Si30_to_n_Si29_derived = 984,
        k_Si30_to_He4_Mg26_derived = 985,
        k_Si31_to_n_Si30_derived = 986,
        k_Si32_to_n_Si31_derived = 987,
        k_P29_to_p_Si28_derived = 988,
        k_P29_to_He4_Al25_derived = 989,
        k_P30_to_n_P29_derived = 990,
        k_P30_to_p_Si29_derived = 991,
        k_P30_to_He4_Al26_derived = 992,
        k_P31_to_n_P30_derived = 993,
        k_P31_to_p_Si30_derived = 994,
        k_P31_to_He4_Al27_derived = 995,
        k_P32_to_n_P31_derived = 996,
        k_P32_to_p_Si31_derived = 997,
        k_P33_to_n_P32_derived = 998,
        k_P33_to_p_Si32_derived = 999,
        k_S32_to_p_P31_derived = 1000,
        k_S32_to_He4_Si28_derived = 1001,
        k_S33_to_n_S32_derived = 1002,
        k_S33_to_p_P32_derived = 1003,
        k_S33_to_He4_Si29_derived = 1004,
        k_S34_to_n_S33_derived = 1005,
        k_S34_to_p_P33_derived = 1006,
        k_S34_to_He4_Si30_derived = 1007,
        k_S35_to_n_S34_derived = 1008,
        k_S35_to_He4_Si31_derived = 1009,
        k_S36_to_n_S35_derived = 1010,
        k_S36_to_He4_Si32_derived = 1011,
        k_Cl33_to_p_S32_derived = 1012,
        k_Cl33_to_He4_P29_derived = 1013,
        k_Cl34_to_n_Cl33_derived = 1014,
        k_Cl34_to_p_S33_derived = 1015,
        k_Cl34_to_He4_P30_derived = 1016,
        k_Cl35_to_n_Cl34_derived = 1017,
        k_Cl35_to_p_S34_derived = 1018,
        k_Cl35_to_He4_P31_derived = 1019,
        k_Cl36_to_n_Cl35_derived = 1020,
        k_Cl36_to_p_S35_derived = 1021,
        k_Cl36_to_He4_P32_derived = 1022,
        k_Cl37_to_n_Cl36_derived = 1023,
        k_Cl37_to_p_S36_derived = 1024,
        k_Cl37_to_He4_P33_derived = 1025,
        k_Ar36_to_p_Cl35_derived = 1026,
        k_Ar36_to_He4_S32_derived = 1027,
        k_Ar37_to_n_Ar36_derived = 1028,
        k_Ar37_to_p_Cl36_derived = 1029,
        k_Ar37_to_He4_S33_derived = 1030,
        k_Ar38_to_n_Ar37_derived = 1031,
        k_Ar38_to_p_Cl37_derived = 1032,
        k_Ar38_to_He4_S34_derived = 1033,
        k_Ar39_to_n_Ar38_derived = 1034,
        k_Ar39_to_He4_S35_derived = 1035,
        k_Ar40_to_n_Ar39_derived = 1036,
        k_Ar40_to_He4_S36_derived = 1037,
        k_K37_to_p_Ar36_derived = 1038,
        k_K37_to_He4_Cl33_derived = 1039,
        k_K38_to_n_K37_derived = 1040,
        k_K38_to_p_Ar37_derived = 1041,
        k_K38_to_He4_Cl34_derived = 1042,
        k_K39_to_n_K38_derived = 1043,
        k_K39_to_p_Ar38_derived = 1044,
        k_K39_to_He4_Cl35_derived = 1045,
        k_K40_to_n_K39_derived = 1046,
        k_K40_to_p_Ar39_derived = 1047,
        k_K40_to_He4_Cl36_derived = 1048,
        k_K41_to_n_K40_derived = 1049,
        k_K41_to_p_Ar40_derived = 1050,
        k_K41_to_He4_Cl37_derived = 1051,
        k_Ca40_to_p_K39_derived = 1052,
        k_Ca40_to_He4_Ar36_derived = 1053,
        k_Ca41_to_n_Ca40_derived = 1054,
        k_Ca41_to_p_K40_derived = 1055,
        k_Ca41_to_He4_Ar37_derived = 1056,
        k_Ca42_to_n_Ca41_derived = 1057,
        k_Ca42_to_p_K41_derived = 1058,
        k_Ca42_to_He4_Ar38_derived = 1059,
        k_Ca43_to_n_Ca42_derived = 1060,
        k_Ca43_to_He4_Ar39_derived = 1061,
        k_Ca44_to_n_Ca43_derived = 1062,
        k_Ca44_to_He4_Ar40_derived = 1063,
        k_Ca45_to_n_Ca44_derived = 1064,
        k_Ca46_to_n_Ca45_derived = 1065,
        k_Ca47_to_n_Ca46_derived = 1066,
        k_Ca48_to_n_Ca47_derived = 1067,
        k_Sc43_to_p_Ca42_derived = 1068,
        k_Sc43_to_He4_K39_derived = 1069,
        k_Sc44_to_n_Sc43_derived = 1070,
        k_Sc44_to_p_Ca43_derived = 1071,
        k_Sc44_to_He4_K40_derived = 1072,
        k_Sc45_to_n_Sc44_derived = 1073,
        k_Sc45_to_p_Ca44_derived = 1074,
        k_Sc45_to_He4_K41_derived = 1075,
        k_Sc46_to_n_Sc45_derived = 1076,
        k_Sc46_to_p_Ca45_derived = 1077,
        k_Sc47_to_n_Sc46_derived = 1078,
        k_Sc47_to_p_Ca46_derived = 1079,
        k_Sc48_to_n_Sc47_derived = 1080,
        k_Sc48_to_p_Ca47_derived = 1081,
        k_Sc49_to_n_Sc48_derived = 1082,
        k_Sc49_to_p_Ca48_derived = 1083,
        k_Ti44_to_p_Sc43_derived = 1084,
        k_Ti44_to_He4_Ca40_derived = 1085,
        k_Ti45_to_n_Ti44_derived = 1086,
        k_Ti45_to_p_Sc44_derived = 1087,
        k_Ti45_to_He4_Ca41_derived = 1088,
        k_Ti46_to_n_Ti45_derived = 1089,
        k_Ti46_to_p_Sc45_derived = 1090,
        k_Ti46_to_He4_Ca42_derived = 1091,
        k_Ti47_to_n_Ti46_derived = 1092,
        k_Ti47_to_p_Sc46_derived = 1093,
        k_Ti47_to_He4_Ca43_derived = 1094,
        k_Ti48_to_n_Ti47_derived = 1095,
        k_Ti48_to_p_Sc47_derived = 1096,
        k_Ti48_to_He4_Ca44_derived = 1097,
        k_Ti49_to_n_Ti48_derived = 1098,
        k_Ti49_to_p_Sc48_derived = 1099,
        k_Ti49_to_He4_Ca45_derived = 1100,
        k_Ti50_to_n_Ti49_derived = 1101,
        k_Ti50_to_p_Sc49_derived = 1102,
        k_Ti50_to_He4_Ca46_derived = 1103,
        k_Ti51_to_n_Ti50_derived = 1104,
        k_Ti51_to_He4_Ca47_derived = 1105,
        k_V46_to_p_Ti45_derived = 1106,
        k_V47_to_n_V46_derived = 1107,
        k_V47_to_p_Ti46_derived = 1108,
        k_V47_to_He4_Sc43_derived = 1109,
        k_V48_to_n_V47_derived = 1110,
        k_V48_to_p_Ti47_derived = 1111,
        k_V48_to_He4_Sc44_derived = 1112,
        k_V49_to_n_V48_derived = 1113,
        k_V49_to_p_Ti48_derived = 1114,
        k_V49_to_He4_Sc45_derived = 1115,
        k_V50_to_n_V49_derived = 1116,
        k_V50_to_p_Ti49_derived = 1117,
        k_V50_to_He4_Sc46_derived = 1118,
        k_V51_to_n_V50_derived = 1119,
        k_V51_to_p_Ti50_derived = 1120,
        k_V51_to_He4_Sc47_derived = 1121,
        k_V52_to_n_V51_derived = 1122,
        k_V52_to_p_Ti51_derived = 1123,
        k_V52_to_He4_Sc48_derived = 1124,
        k_Cr48_to_p_V47_derived = 1125,
        k_Cr48_to_He4_Ti44_derived = 1126,
        k_Cr49_to_n_Cr48_derived = 1127,
        k_Cr49_to_p_V48_derived = 1128,
        k_Cr49_to_He4_Ti45_derived = 1129,
        k_Cr50_to_n_Cr49_derived = 1130,
        k_Cr50_to_p_V49_derived = 1131,
        k_Cr50_to_He4_Ti46_derived = 1132,
        k_Cr51_to_n_Cr50_derived = 1133,
        k_Cr51_to_p_V50_derived = 1134,
        k_Cr51_to_He4_Ti47_derived = 1135,
        k_Cr52_to_n_Cr51_derived = 1136,
        k_Cr52_to_p_V51_derived = 1137,
        k_Cr52_to_He4_Ti48_derived = 1138,
        k_Cr53_to_n_Cr52_derived = 1139,
        k_Cr53_to_p_V52_derived = 1140,
        k_Cr53_to_He4_Ti49_derived = 1141,
        k_Cr54_to_n_Cr53_derived = 1142,
        k_Cr54_to_He4_Ti50_derived = 1143,
        k_Mn50_to_p_Cr49_derived = 1144,
        k_Mn50_to_He4_V46_derived = 1145,
        k_Mn51_to_n_Mn50_derived = 1146,
        k_Mn51_to_p_Cr50_derived = 1147,
        k_Mn51_to_He4_V47_derived = 1148,
        k_Mn52_to_n_Mn51_derived = 1149,
        k_Mn52_to_p_Cr51_derived = 1150,
        k_Mn52_to_He4_V48_derived = 1151,
        k_Mn53_to_n_Mn52_derived = 1152,
        k_Mn53_to_p_Cr52_derived = 1153,
        k_Mn53_to_He4_V49_derived = 1154,
        k_Mn54_to_n_Mn53_derived = 1155,
        k_Mn54_to_p_Cr53_derived = 1156,
        k_Mn54_to_He4_V50_derived = 1157,
        k_Mn55_to_n_Mn54_derived = 1158,
        k_Mn55_to_p_Cr54_derived = 1159,
        k_Mn55_to_He4_V51_derived = 1160,
        k_Fe52_to_p_Mn51_derived = 1161,
        k_Fe52_to_He4_Cr48_derived = 1162,
        k_Fe53_to_n_Fe52_derived = 1163,
        k_Fe53_to_p_Mn52_derived = 1164,
        k_Fe53_to_He4_Cr49_derived = 1165,
        k_Fe54_to_n_Fe53_derived = 1166,
        k_Fe54_to_p_Mn53_derived = 1167,
        k_Fe54_to_He4_Cr50_derived = 1168,
        k_Fe55_to_n_Fe54_derived = 1169,
        k_Fe55_to_p_Mn54_derived = 1170,
        k_Fe55_to_He4_Cr51_derived = 1171,
        k_Fe56_to_n_Fe55_derived = 1172,
        k_Fe56_to_p_Mn55_derived = 1173,
        k_Fe56_to_He4_Cr52_derived = 1174,
        k_Fe57_to_n_Fe56_derived = 1175,
        k_Fe57_to_He4_Cr53_derived = 1176,
        k_Fe58_to_n_Fe57_derived = 1177,
        k_Fe58_to_He4_Cr54_derived = 1178,
        k_Co53_to_p_Fe52_derived = 1179,
        k_Co54_to_n_Co53_derived = 1180,
        k_Co54_to_p_Fe53_derived = 1181,
        k_Co54_to_He4_Mn50_derived = 1182,
        k_Co55_to_n_Co54_derived = 1183,
        k_Co55_to_p_Fe54_derived = 1184,
        k_Co55_to_He4_Mn51_derived = 1185,
        k_Co56_to_n_Co55_derived = 1186,
        k_Co56_to_p_Fe55_derived = 1187,
        k_Co56_to_He4_Mn52_derived = 1188,
        k_Co57_to_n_Co56_derived = 1189,
        k_Co57_to_p_Fe56_derived = 1190,
        k_Co57_to_He4_Mn53_derived = 1191,
        k_Co58_to_n_Co57_derived = 1192,
        k_Co58_to_p_Fe57_derived = 1193,
        k_Co58_to_He4_Mn54_derived = 1194,
        k_Co59_to_n_Co58_derived = 1195,
        k_Co59_to_p_Fe58_derived = 1196,
        k_Co59_to_He4_Mn55_derived = 1197,
        k_Ni56_to_p_Co55_derived = 1198,
        k_Ni56_to_He4_Fe52_derived = 1199,
        k_Ni57_to_n_Ni56_derived = 1200,
        k_Ni57_to_p_Co56_derived = 1201,
        k_Ni57_to_He4_Fe53_derived = 1202,
        k_Ni58_to_n_Ni57_derived = 1203,
        k_Ni58_to_p_Co57_derived = 1204,
        k_Ni58_to_He4_Fe54_derived = 1205,
        k_Ni59_to_n_Ni58_derived = 1206,
        k_Ni59_to_p_Co58_derived = 1207,
        k_Ni59_to_He4_Fe55_derived = 1208,
        k_Ni60_to_n_Ni59_derived = 1209,
        k_Ni60_to_p_Co59_derived = 1210,
        k_Ni60_to_He4_Fe56_derived = 1211,
        k_Ni61_to_n_Ni60_derived = 1212,
        k_Ni61_to_He4_Fe57_derived = 1213,
        k_Ni62_to_n_Ni61_derived = 1214,
        k_Ni62_to_He4_Fe58_derived = 1215,
        k_Ni63_to_n_Ni62_derived = 1216,
        k_Ni64_to_n_Ni63_derived = 1217,
        k_Cu57_to_p_Ni56_derived = 1218,
        k_Cu57_to_He4_Co53_derived = 1219,
        k_Cu58_to_n_Cu57_derived = 1220,
        k_Cu58_to_p_Ni57_derived = 1221,
        k_Cu58_to_He4_Co54_derived = 1222,
        k_Cu59_to_n_Cu58_derived = 1223,
        k_Cu59_to_p_Ni58_derived = 1224,
        k_Cu59_to_He4_Co55_derived = 1225,
        k_Cu60_to_n_Cu59_derived = 1226,
        k_Cu60_to_p_Ni59_derived = 1227,
        k_Cu60_to_He4_Co56_derived = 1228,
        k_Cu61_to_n_Cu60_derived = 1229,
        k_Cu61_to_p_Ni60_derived = 1230,
        k_Cu61_to_He4_Co57_derived = 1231,
        k_Cu62_to_n_Cu61_derived = 1232,
        k_Cu62_to_p_Ni61_derived = 1233,
        k_Cu62_to_He4_Co58_derived = 1234,
        k_Cu63_to_n_Cu62_derived = 1235,
        k_Cu63_to_p_Ni62_derived = 1236,
        k_Cu63_to_He4_Co59_derived = 1237,
        k_Cu64_to_n_Cu63_derived = 1238,
        k_Cu64_to_p_Ni63_derived = 1239,
        k_Cu65_to_n_Cu64_derived = 1240,
        k_Cu65_to_p_Ni64_derived = 1241,
        k_Zn59_to_p_Cu58_derived = 1242,
        k_Zn60_to_n_Zn59_derived = 1243,
        k_Zn60_to_p_Cu59_derived = 1244,
        k_Zn60_to_He4_Ni56_derived = 1245,
        k_Zn61_to_n_Zn60_derived = 1246,
        k_Zn61_to_p_Cu60_derived = 1247,
        k_Zn61_to_He4_Ni57_derived = 1248,
        k_Zn62_to_n_Zn61_derived = 1249,
        k_Zn62_to_p_Cu61_derived = 1250,
        k_Zn62_to_He4_Ni58_derived = 1251,
        k_Zn63_to_n_Zn62_derived = 1252,
        k_Zn63_to_p_Cu62_derived = 1253,
        k_Zn63_to_He4_Ni59_derived = 1254,
        k_Zn64_to_n_Zn63_derived = 1255,
        k_Zn64_to_p_Cu63_derived = 1256,
        k_Zn64_to_He4_Ni60_derived = 1257,
        k_Zn65_to_n_Zn64_derived = 1258,
        k_Zn65_to_p_Cu64_derived = 1259,
        k_Zn65_to_He4_Ni61_derived = 1260,
        k_Zn66_to_n_Zn65_derived = 1261,
        k_Zn66_to_p_Cu65_derived = 1262,
        k_Zn66_to_He4_Ni62_derived = 1263,
        k_Ga62_to_p_Zn61_derived = 1264,
        k_Ga62_to_He4_Cu58_derived = 1265,
        k_Ga63_to_n_Ga62_derived = 1266,
        k_Ga63_to_p_Zn62_derived = 1267,
        k_Ga63_to_He4_Cu59_derived = 1268,
        k_Ga64_to_n_Ga63_derived = 1269,
        k_Ga64_to_p_Zn63_derived = 1270,
        k_Ga64_to_He4_Cu60_derived = 1271,
        k_Ge63_to_p_Ga62_derived = 1272,
        k_Ge63_to_He4_Zn59_derived = 1273,
        k_Ge64_to_n_Ge63_derived = 1274,
        k_Ge64_to_p_Ga63_derived = 1275,
        k_Ge64_to_He4_Zn60_derived = 1276,
        k_Li6_to_n_p_He4_derived = 1277,
        k_Be9_to_n_He4_He4_derived = 1278,
        k_C12_to_3He4_derived = 1279,
        k_He3_n_to_d_d_derived = 1280,
        k_He3_He4_to_p_Li6_derived = 1281,
        k_He4_p_to_d_He3_derived = 1282,
        k_He4_He4_to_n_Be7_derived = 1283,
        k_He4_He4_to_p_Li7_derived = 1284,
        k_Li6_He4_to_p_Be9_derived = 1285,
        k_Li7_p_to_n_Be7_derived = 1286,
        k_Li7_p_to_d_Li6_derived = 1287,
        k_Li7_He4_to_n_B10_derived = 1288,
        k_Be7_n_to_d_Li6_derived = 1289,
        k_Be7_He4_to_p_B10_derived = 1290,
        k_C12_n_to_He4_Be9_derived = 1291,
        k_C12_He4_to_n_O15_derived = 1292,
        k_C12_He4_to_p_N15_derived = 1293,
        k_C12_C12_to_n_Mg23_derived = 1294,
        k_C13_p_to_He4_B10_derived = 1295,
        k_C14_p_to_n_N14_derived = 1296,
        k_C14_p_to_He4_B11_derived = 1297,
        k_C14_He4_to_n_O17_derived = 1298,
        k_N13_n_to_p_C13_derived = 1299,
        k_N13_n_to_He4_B10_derived = 1300,
        k_N14_n_to_d_C13_derived = 1301,
        k_N14_n_to_He4_B11_derived = 1302,
        k_N14_He4_to_p_O17_derived = 1303,
        k_N15_n_to_d_C14_derived = 1304,
        k_N15_He4_to_n_F18_derived = 1305,
        k_N15_He4_to_p_O18_derived = 1306,
        k_O14_n_to_p_N14_derived = 1307,
        k_O15_n_to_p_N15_derived = 1308,
        k_O15_He4_to_n_Ne18_derived = 1309,
        k_O15_He4_to_p_F18_derived = 1310,
        k_O16_n_to_He4_C13_derived = 1311,
        k_O16_p_to_He4_N13_derived = 1312,
        k_O16_He4_to_n_Ne19_derived = 1313,
        k_O16_He4_to_p_F19_derived = 1314,
        k_O17_p_to_n_F17_derived = 1315,
        k_O18_p_to_n_F18_derived = 1316,
        k_F17_n_to_He4_N14_derived = 1317,
        k_F17_p_to_He4_O14_derived = 1318,
        k_F17_He4_to_p_Ne20_derived = 1319,
        k_F18_p_to_n_Ne18_derived = 1320,
        k_F18_He4_to_n_Na21_derived = 1321,
        k_F19_He4_to_n_Na22_derived = 1322,
        k_Ne19_n_to_p_F19_derived = 1323,
        k_Ne20_n_to_He4_O17_derived = 1324,
        k_Ne20_He4_to_n_Mg23_derived = 1325,
        k_Ne20_He4_to_p_Na23_derived = 1326,
        k_Ne20_He4_to_C12_C12_derived = 1327,
        k_Ne21_n_to_He4_O18_derived = 1328,
        k_Ne21_p_to_n_Na21_derived = 1329,
        k_Ne21_p_to_He4_F18_derived = 1330,
        k_Ne22_p_to_n_Na22_derived = 1331,
        k_Ne22_p_to_He4_F19_derived = 1332,
        k_Na21_p_to_He4_Ne18_derived = 1333,
        k_Na21_He4_to_p_Mg24_derived = 1334,
        k_Na22_p_to_He4_Ne19_derived = 1335,
        k_Na22_He4_to_n_Al25_derived = 1336,
        k_Na23_p_to_C12_C12_derived = 1337,
        k_Na23_He4_to_n_Al26_derived = 1338,
        k_Mg23_n_to_p_Na23_derived = 1339,
        k_Mg24_n_to_He4_Ne21_derived = 1340,
        k_Mg24_He4_to_p_Al27_derived = 1341,
        k_Mg24_He4_to_C12_O16_derived = 1342,
        k_Mg25_n_to_He4_Ne22_derived = 1343,
        k_Mg25_p_to_n_Al25_derived = 1344,
        k_Mg25_p_to_He4_Na22_derived = 1345,
        k_Mg26_p_to_n_Al26_derived = 1346,
        k_Mg26_p_to_He4_Na23_derived = 1347,
        k_Al26_p_to_He4_Mg23_derived = 1348,
        k_Al26_He4_to_n_P29_derived = 1349,
        k_Al27_p_to_C12_O16_derived = 1350,
        k_Si28_n_to_He4_Mg25_derived = 1351,
        k_Si28_p_to_He4_Al25_derived = 1352,
        k_Si28_He4_to_p_P31_derived = 1353,
        k_Si28_He4_to_C12_Ne20_derived = 1354,
        k_Si28_He4_to_O16_O16_derived = 1355,
        k_Si29_n_to_He4_Mg26_derived = 1356,
        k_Si29_p_to_n_P29_derived = 1357,
        k_Si29_p_to_He4_Al26_derived = 1358,
        k_Si29_He4_to_n_S32_derived = 1359,
        k_Si29_He4_to_p_P32_derived = 1360,
        k_Si30_p_to_n_P30_derived = 1361,
        k_Si30_p_to_He4_Al27_derived = 1362,
        k_Si30_He4_to_n_S33_derived = 1363,
        k_Si30_He4_to_p_P33_derived = 1364,
        k_Si32_p_to_n_P32_derived = 1365,
        k_Si32_He4_to_n_S35_derived = 1366,
        k_P30_n_to_He4_Al27_derived = 1367,
        k_P30_He4_to_n_Cl33_derived = 1368,
        k_P31_n_to_p_Si31_derived = 1369,
        k_P31_p_to_C12_Ne20_derived = 1370,
        k_P31_p_to_O16_O16_derived = 1371,
        k_P31_He4_to_n_Cl34_derived = 1372,
        k_P32_He4_to_n_Cl35_derived = 1373,
        k_P32_He4_to_p_S35_derived = 1374,
        k_P33_p_to_n_S33_derived = 1375,
        k_P33_He4_to_n_Cl36_derived = 1376,
        k_P33_He4_to_p_S36_derived = 1377,
        k_S32_n_to_p_P32_derived = 1378,
        k_S32_p_to_He4_P29_derived = 1379,
        k_S32_He4_to_p_Cl35_derived = 1380,
        k_S33_p_to_n_Cl33_derived = 1381,
        k_S33_p_to_He4_P30_derived = 1382,
        k_S33_He4_to_n_Ar36_derived = 1383,
        k_S33_He4_to_p_Cl36_derived = 1384,
        k_S34_n_to_He4_Si31_derived = 1385,
        k_S34_p_to_n_Cl34_derived = 1386,
        k_S34_p_to_He4_P31_derived = 1387,
        k_S35_p_to_n_Cl35_derived = 1388,
        k_S36_p_to_n_Cl36_derived = 1389,
        k_S36_He4_to_n_Ar39_derived = 1390,
        k_Cl34_He4_to_n_K37_derived = 1391,
        k_Cl35_He4_to_n_K38_derived = 1392,
        k_Cl36_p_to_n_Ar36_derived = 1393,
        k_Cl36_He4_to_n_K39_derived = 1394,
        k_Cl36_He4_to_p_Ar39_derived = 1395,
        k_Cl37_p_to_n_Ar37_derived = 1396,
        k_Cl37_p_to_He4_S34_derived = 1397,
        k_Cl37_He4_to_n_K40_derived = 1398,
        k_Cl37_He4_to_p_Ar40_derived = 1399,
        k_Ar36_p_to_He4_Cl33_derived = 1400,
        k_Ar36_He4_to_p_K39_derived = 1401,
        k_Ar37_n_to_He4_S34_derived = 1402,
        k_Ar37_p_to_n_K37_derived = 1403,
        k_Ar37_p_to_He4_Cl34_derived = 1404,
        k_Ar37_He4_to_n_Ca40_derived = 1405,
        k_Ar37_He4_to_p_K40_derived = 1406,
        k_Ar38_n_to_He4_S35_derived = 1407,
        k_Ar38_p_to_n_K38_derived = 1408,
        k_Ar38_p_to_He4_Cl35_derived = 1409,
        k_Ar38_He4_to_n_Ca41_derived = 1410,
        k_Ar38_He4_to_p_K41_derived = 1411,
        k_Ar39_p_to_n_K39_derived = 1412,
        k_Ar39_He4_to_n_Ca42_derived = 1413,
        k_Ar40_p_to_n_K40_derived = 1414,
        k_Ar40_He4_to_n_Ca43_derived = 1415,
        k_K39_He4_to_p_Ca42_derived = 1416,
        k_K40_He4_to_n_Sc43_derived = 1417,
        k_Ca40_n_to_p_K40_derived = 1418,
        k_Ca40_p_to_He4_K37_derived = 1419,
        k_Ca40_He4_to_p_Sc43_derived = 1420,
        k_Ca41_n_to_p_K41_derived = 1421,
        k_Ca41_p_to_He4_K38_derived = 1422,
        k_Ca41_He4_to_n_Ti44_derived = 1423,
        k_Ca41_He4_to_p_Sc44_derived = 1424,
        k_Ca42_He4_to_n_Ti45_derived = 1425,
        k_Ca42_He4_to_p_Sc45_derived = 1426,
        k_Ca43_p_to_n_Sc43_derived = 1427,
        k_Ca43_p_to_He4_K40_derived = 1428,
        k_Ca43_He4_to_p_Sc46_derived = 1429,
        k_Ca44_p_to_n_Sc44_derived = 1430,
        k_Ca44_p_to_He4_K41_derived = 1431,
        k_Ca44_He4_to_n_Ti47_derived = 1432,
        k_Ca44_He4_to_p_Sc47_derived = 1433,
        k_Ca45_p_to_n_Sc45_derived = 1434,
        k_Ca45_He4_to_p_Sc48_derived = 1435,
        k_Ca46_p_to_n_Sc46_derived = 1436,
        k_Ca46_He4_to_n_Ti49_derived = 1437,
        k_Ca46_He4_to_p_Sc49_derived = 1438,
        k_Ca48_p_to_n_Sc48_derived = 1439,
        k_Ca48_He4_to_n_Ti51_derived = 1440,
        k_Sc43_He4_to_n_V46_derived = 1441,
        k_Sc44_n_to_He4_K41_derived = 1442,
        k_Sc44_p_to_n_Ti44_derived = 1443,
        k_Sc44_He4_to_n_V47_derived = 1444,
        k_Sc45_p_to_n_Ti45_derived = 1445,
        k_Sc45_He4_to_n_V48_derived = 1446,
        k_Sc47_n_to_p_Ca47_derived = 1447,
        k_Sc47_p_to_n_Ti47_derived = 1448,
        k_Sc47_He4_to_n_V50_derived = 1449,
        k_Sc49_He4_to_n_V52_derived = 1450,
        k_Ti45_He4_to_n_Cr48_derived = 1451,
        k_Ti46_n_to_p_Sc46_derived = 1452,
        k_Ti46_n_to_He4_Ca43_derived = 1453,
        k_Ti46_p_to_n_V46_derived = 1454,
        k_Ti46_p_to_He4_Sc43_derived = 1455,
        k_Ti46_He4_to_n_Cr49_derived = 1456,
        k_Ti46_He4_to_p_V49_derived = 1457,
        k_Ti47_p_to_n_V47_derived = 1458,
        k_Ti47_p_to_He4_Sc44_derived = 1459,
        k_Ti47_He4_to_n_Cr50_derived = 1460,
        k_Ti47_He4_to_p_V50_derived = 1461,
        k_Ti48_n_to_p_Sc48_derived = 1462,
        k_Ti48_n_to_He4_Ca45_derived = 1463,
        k_Ti48_p_to_n_V48_derived = 1464,
        k_Ti48_p_to_He4_Sc45_derived = 1465,
        k_Ti48_He4_to_n_Cr51_derived = 1466,
        k_Ti48_He4_to_p_V51_derived = 1467,
        k_Ti49_n_to_p_Sc49_derived = 1468,
        k_Ti49_p_to_n_V49_derived = 1469,
        k_Ti49_p_to_He4_Sc46_derived = 1470,
        k_Ti49_He4_to_p_V52_derived = 1471,
        k_Ti50_n_to_He4_Ca47_derived = 1472,
        k_Ti50_p_to_n_V50_derived = 1473,
        k_Ti50_p_to_He4_Sc47_derived = 1474,
        k_Ti50_He4_to_n_Cr53_derived = 1475,
        k_Ti51_p_to_He4_Sc48_derived = 1476,
        k_V47_p_to_He4_Ti44_derived = 1477,
        k_V47_He4_to_n_Mn50_derived = 1478,
        k_V48_p_to_n_Cr48_derived = 1479,
        k_V48_p_to_He4_Ti45_derived = 1480,
        k_V48_He4_to_n_Mn51_derived = 1481,
        k_V49_n_to_He4_Sc46_derived = 1482,
        k_V49_p_to_n_Cr49_derived = 1483,
        k_V49_He4_to_n_Mn52_derived = 1484,
        k_V50_He4_to_n_Mn53_derived = 1485,
        k_V51_n_to_p_Ti51_derived = 1486,
        k_V51_n_to_He4_Sc48_derived = 1487,
        k_V51_p_to_n_Cr51_derived = 1488,
        k_V51_He4_to_n_Mn54_derived = 1489,
        k_V51_He4_to_p_Cr54_derived = 1490,
        k_Cr49_p_to_He4_V46_derived = 1491,
        k_Cr49_He4_to_n_Fe52_derived = 1492,
        k_Cr50_n_to_p_V50_derived = 1493,
        k_Cr50_p_to_n_Mn50_derived = 1494,
        k_Cr50_p_to_He4_V47_derived = 1495,
        k_Cr50_He4_to_n_Fe53_derived = 1496,
        k_Cr50_He4_to_p_Mn53_derived = 1497,
        k_Cr51_p_to_n_Mn51_derived = 1498,
        k_Cr51_p_to_He4_V48_derived = 1499,
        k_Cr51_He4_to_n_Fe54_derived = 1500,
        k_Cr51_He4_to_p_Mn54_derived = 1501,
        k_Cr52_n_to_p_V52_derived = 1502,
        k_Cr52_n_to_He4_Ti49_derived = 1503,
        k_Cr52_p_to_n_Mn52_derived = 1504,
        k_Cr52_p_to_He4_V49_derived = 1505,
        k_Cr52_He4_to_n_Fe55_derived = 1506,
        k_Cr52_He4_to_p_Mn55_derived = 1507,
        k_Cr53_p_to_n_Mn53_derived = 1508,
        k_Cr53_p_to_He4_V50_derived = 1509,
        k_Cr53_He4_to_n_Fe56_derived = 1510,
        k_Cr54_n_to_He4_Ti51_derived = 1511,
        k_Cr54_p_to_n_Mn54_derived = 1512,
        k_Cr54_He4_to_n_Fe57_derived = 1513,
        k_Mn50_He4_to_n_Co53_derived = 1514,
        k_Mn51_p_to_He4_Cr48_derived = 1515,
        k_Mn51_He4_to_n_Co54_derived = 1516,
        k_Mn52_p_to_n_Fe52_derived = 1517,
        k_Mn52_p_to_He4_Cr49_derived = 1518,
        k_Mn52_He4_to_n_Co55_derived = 1519,
        k_Mn53_p_to_n_Fe53_derived = 1520,
        k_Mn53_He4_to_n_Co56_derived = 1521,
        k_Mn54_p_to_n_Fe54_derived = 1522,
        k_Mn54_He4_to_n_Co57_derived = 1523,
        k_Mn54_He4_to_p_Fe57_derived = 1524,
        k_Mn55_n_to_He4_V52_derived = 1525,
        k_Mn55_p_to_n_Fe55_derived = 1526,
        k_Mn55_He4_to_n_Co58_derived = 1527,
        k_Mn55_He4_to_p_Fe58_derived = 1528,
        k_Fe53_p_to_n_Co53_derived = 1529,
        k_Fe53_p_to_He4_Mn50_derived = 1530,
        k_Fe53_He4_to_n_Ni56_derived = 1531,
        k_Fe54_p_to_n_Co54_derived = 1532,
        k_Fe54_p_to_He4_Mn51_derived = 1533,
        k_Fe54_He4_to_n_Ni57_derived = 1534,
        k_Fe54_He4_to_p_Co57_derived = 1535,
        k_Fe55_p_to_n_Co55_derived = 1536,
        k_Fe55_p_to_He4_Mn52_derived = 1537,
        k_Fe55_He4_to_n_Ni58_derived = 1538,
        k_Fe55_He4_to_p_Co58_derived = 1539,
        k_Fe56_p_to_n_Co56_derived = 1540,
        k_Fe56_p_to_He4_Mn53_derived = 1541,
        k_Fe56_He4_to_n_Ni59_derived = 1542,
        k_Fe56_He4_to_p_Co59_derived = 1543,
        k_Fe57_p_to_n_Co57_derived = 1544,
        k_Fe57_He4_to_n_Ni60_derived = 1545,
        k_Fe58_p_to_n_Co58_derived = 1546,
        k_Fe58_He4_to_n_Ni61_derived = 1547,
        k_Co54_He4_to_n_Cu57_derived = 1548,
        k_Co55_p_to_He4_Fe52_derived = 1549,
        k_Co55_He4_to_n_Cu58_derived = 1550,
        k_Co56_p_to_n_Ni56_derived = 1551,
        k_Co56_p_to_He4_Fe53_derived = 1552,
        k_Co56_He4_to_n_Cu59_derived = 1553,
        k_Co57_p_to_n_Ni57_derived = 1554,
        k_Co57_He4_to_n_Cu60_derived = 1555,
        k_Co58_p_to_n_Ni58_derived = 1556,
        k_Co58_He4_to_n_Cu61_derived = 1557,
        k_Co58_He4_to_p_Ni61_derived = 1558,
        k_Co59_p_to_n_Ni59_derived = 1559,
        k_Co59_He4_to_n_Cu62_derived = 1560,
        k_Co59_He4_to_p_Ni62_derived = 1561,
        k_Ni56_p_to_He4_Co53_derived = 1562,
        k_Ni56_He4_to_n_Zn59_derived = 1563,
        k_Ni56_He4_to_p_Cu59_derived = 1564,
        k_Ni57_p_to_n_Cu57_derived = 1565,
        k_Ni57_p_to_He4_Co54_derived = 1566,
        k_Ni57_He4_to_n_Zn60_derived = 1567,
        k_Ni57_He4_to_p_Cu60_derived = 1568,
        k_Ni58_p_to_n_Cu58_derived = 1569,
        k_Ni58_p_to_He4_Co55_derived = 1570,
        k_Ni58_He4_to_n_Zn61_derived = 1571,
        k_Ni58_He4_to_p_Cu61_derived = 1572,
        k_Ni59_p_to_n_Cu59_derived = 1573,
        k_Ni59_p_to_He4_Co56_derived = 1574,
        k_Ni59_He4_to_n_Zn62_derived = 1575,
        k_Ni59_He4_to_p_Cu62_derived = 1576,
        k_Ni60_p_to_n_Cu60_derived = 1577,
        k_Ni60_p_to_He4_Co57_derived = 1578,
        k_Ni60_He4_to_n_Zn63_derived = 1579,
        k_Ni60_He4_to_p_Cu63_derived = 1580,
        k_Ni61_p_to_n_Cu61_derived = 1581,
        k_Ni61_He4_to_n_Zn64_derived = 1582,
        k_Ni61_He4_to_p_Cu64_derived = 1583,
        k_Ni62_p_to_n_Cu62_derived = 1584,
        k_Ni62_He4_to_n_Zn65_derived = 1585,
        k_Ni62_He4_to_p_Cu65_derived = 1586,
        k_Ni63_p_to_n_Cu63_derived = 1587,
        k_Ni63_He4_to_n_Zn66_derived = 1588,
        k_Ni64_p_to_n_Cu64_derived = 1589,
        k_Cu58_He4_to_p_Zn61_derived = 1590,
        k_Cu59_p_to_n_Zn59_derived = 1591,
        k_Cu59_He4_to_n_Ga62_derived = 1592,
        k_Cu59_He4_to_p_Zn62_derived = 1593,
        k_Cu60_p_to_n_Zn60_derived = 1594,
        k_Cu60_He4_to_n_Ga63_derived = 1595,
        k_Cu60_He4_to_p_Zn63_derived = 1596,
        k_Cu61_p_to_n_Zn61_derived = 1597,
        k_Cu61_He4_to_n_Ga64_derived = 1598,
        k_Cu61_He4_to_p_Zn64_derived = 1599,
        k_Cu62_p_to_n_Zn62_derived = 1600,
        k_Cu62_He4_to_p_Zn65_derived = 1601,
        k_Cu63_p_to_n_Zn63_derived = 1602,
        k_Cu63_He4_to_p_Zn66_derived = 1603,
        k_Cu64_p_to_n_Zn64_derived = 1604,
        k_Cu65_p_to_n_Zn65_derived = 1605,
        k_Zn59_He4_to_p_Ga62_derived = 1606,
        k_Zn60_p_to_He4_Cu57_derived = 1607,
        k_Zn60_He4_to_n_Ge63_derived = 1608,
        k_Zn60_He4_to_p_Ga63_derived = 1609,
        k_Zn61_He4_to_n_Ge64_derived = 1610,
        k_Zn61_He4_to_p_Ga64_derived = 1611,
        k_Zn62_p_to_n_Ga62_derived = 1612,
        k_Zn63_p_to_n_Ga63_derived = 1613,
        k_Zn64_p_to_n_Ga64_derived = 1614,
        k_Ga63_p_to_n_Ge63_derived = 1615,
        k_Ga64_p_to_n_Ge64_derived = 1616,
        k_d_p_to_n_p_p_derived = 1617,
        k_He4_p_p_to_He3_He3_derived = 1618,
        k_He4_He4_n_to_d_Li7_derived = 1619,
        k_He4_He4_p_to_n_B8_derived = 1620,
        k_He4_He4_p_to_d_Be7_derived = 1621,
        k_d_He4_He4_to_p_Be9_derived = 1622,
        k_He4_He4_He4_to_p_B11_derived = 1623,
        k_He4_He4_p_n_to_He3_Li7_derived = 1624,
        k_He4_He4_p_n_to_p_Be9_derived = 1625,
        k_He4_He4_p_p_to_He3_Be7_derived = 1626,
        NumRates = k_He4_He4_p_p_to_He3_Be7_derived
    };

    enum ScreenPairs :
    std::uint16_t
    {
        k_p_Cl36 = 1,
        k_p_F18 = 2,
        k_p_P33 = 3,
        k_p_Ti49 = 4,
        k_p_Si31 = 5,
        k_He4_Mn53 = 6,
        k_p_Cu65 = 7,
        k_p_Li7 = 8,
        k_He4_Be9 = 9,
        k_He4_Mn54 = 10,
        k_p_Ca47 = 11,
        k_He4_Ne19 = 12,
        k_C12_O16 = 13,
        k_He4_Ti48 = 14,
        k_He3_Be7 = 15,
        k_He4_Li7 = 16,
        k_p_Ca42 = 17,
        k_p_Mg24 = 18,
        k_p_Sc43 = 19,
        k_p_K40 = 20,
        k_p_Na22 = 21,
        k_He4_Ca46 = 22,
        k_p_Ar38 = 23,
        k_p_Fe56 = 24,
        k_p_Co58 = 25,
        k_p_Cu61 = 26,
        k_p_Cr54 = 27,
        k_p_S36 = 28,
        k_He4_Ti44 = 29,
        k_p_V52 = 30,
        k_He4_Co57 = 31,
        k_p_Ne20 = 32,
        k_He4_Fe55 = 33,
        k_He4_Ca42 = 34,
        k_He4_Mg24 = 35,
        k_He4_C13 = 36,
        k_He4_K40 = 37,
        k_He4_Na22 = 38,
        k_p_Cr49 = 39,
        k_He4_Ar38 = 40,
        k_p_V47 = 41,
        k_He4_Cr53 = 42,
        k_He4_V51 = 43,
        k_p_Zn63 = 44,
        k_He4_B11 = 45,
        k_p_Cu63 = 46,
        k_p_Sc45 = 47,
        k_He4_Co53 = 48,
        k_p_Ni61 = 49,
        k_p_Co55 = 50,
        k_C12_C12 = 51,
        k_d_d = 52,
        k_p_Fe53 = 53,
        k_He4_Cl33 = 54,
        k_d_He3 = 55,
        k_He4_Cr49 = 56,
        k_p_Al27 = 57,
        k_He4_Ti46 = 58,
        k_p_Li6 = 59,
        k_p_Mg25 = 60,
        k_p_Ne22 = 61,
        k_p_Ca43 = 62,
        k_He4_V47 = 63,
        k_He4_O18 = 64,
        k_He4_Ni60 = 65,
        k_He4_Cu63 = 66,
        k_He4_Sc45 = 67,
        k_He4_Si32 = 68,
        k_p_Mn52 = 69,
        k_He4_Ar40 = 70,
        k_He4_Fe58 = 71,
        k_He4_Ni61 = 72,
        k_p_Ga62 = 73,
        k_p_Zn60 = 74,
        k_He4_Zn61 = 75,
        k_He4_Li6 = 76,
        k_He4_Ni56 = 77,
        k_p_P31 = 78,
        k_He4_Ne22 = 79,
        k_p_Ti47 = 80,
        k_p_Si29 = 81,
        k_p_Ni58 = 82,
        k_p_O16 = 83,
        k_p_Cu60 = 84,
        k_He4_Cr51 = 85,
        k_p_N14 = 86,
        k_He4_Mn52 = 87,
        k_p_Ca45 = 88,
        k_He4_V49 = 89,
        k_d_Be8 = 90,
        k_O16_O16 = 91,
        k_He4_O15 = 92,
        k_He4_Ni63 = 93,
        k_He4_P31 = 94,
        k_He4_N13 = 95,
        k_He4_Ti47 = 96,
        k_He4_Si29 = 97,
        k_He4_B10 = 98,
        k_p_Ga64 = 99,
        k_p_Zn62 = 100,
        k_He4_Ca45 = 101,
        k_p_Na21 = 102,
        k_p_Ar37 = 103,
        k_p_S34 = 104,
        k_p_Cr52 = 105,
        k_p_V50 = 106,
        k_p_Cl37 = 107,
        k_p_F19 = 108,
        k_He4_Ca40 = 109,
        k_p_Sc48 = 110,
        k_He2_Be8 = 111,
        k_p_Ni64 = 112,
        k_He4_K38 = 113,
        k_He4_Ar36 = 114,
        k_He4_Cl36 = 115,
        k_He4_F18 = 116,
        k_He4_P33 = 117,
        k_He4_S34 = 118,
        k_He4_Si31 = 119,
        k_He4_Ti49 = 120,
        k_He4_V50 = 121,
        k_He4_Cr52 = 122,
        k_He4_Ca47 = 123,
        k_He4_Sc48 = 124,
        k_p_Mn55 = 125,
        k_d_C14 = 126,
        k_He4_Cu61 = 127,
        k_p_d = 128,
        k_p_Ne21 = 129,
        k_p_He3 = 130,
        k_He4_Sc43 = 131,
        k_He4_Fe56 = 132,
        k_He4_Co58 = 133,
        k_He4_Cr54 = 134,
        k_He4_S36 = 135,
        k_He4_Mn55 = 136,
        k_p_Ca48 = 137,
        k_He4_V52 = 138,
        k_He4_Ne20 = 139,
        k_p_Cu58 = 140,
        k_p_Co56 = 141,
        k_He4_Mn50 = 142,
        k_p_Si28 = 143,
        k_p_Fe54 = 144,
        k_d_He4 = 145,
        k_He3_Li7 = 146,
        k_p_Ca44 = 147,
        k_p_Mg26 = 148,
        k_p_Zn65 = 149,
        k_p_Be8 = 150,
        k_He4_Cu58 = 151,
        k_He4_Co55 = 152,
        k_He4_Fe53 = 153,
        k_He4_O14 = 154,
        k_He4_Al27 = 155,
        k_He4_Ca43 = 156,
        k_He4_Mg25 = 157,
        k_p_Ga63 = 158,
        k_p_Zn61 = 159,
        k_p_Cl35 = 160,
        k_p_F17 = 161,
        k_p_P32 = 162,
        k_p_S33 = 163,
        k_p_Si30 = 164,
        k_p_Ni59 = 165,
        k_He4_Zn60 = 166,
        k_p_Sc47 = 167,
        k_He4_Ne18 = 168,
        k_He2_He4 = 169,
        k_He4_Cl34 = 170,
        k_p_C12 = 171,
        k_He4_Cu60 = 172,
        k_He4_Ni58 = 173,
        k_He4_O16 = 174,
        k_He4_F17 = 175,
        k_He4_N14 = 176,
        k_He4_P32 = 177,
        k_He4_Si30 = 178,
        k_p_S35 = 179,
        k_p_Ti45 = 180,
        k_d_C13 = 181,
        k_p_Sc49 = 182,
        k_p_p = 183,
        k_He4_Mg23 = 184,
        k_p_O17 = 185,
        k_p_C14 = 186,
        k_He4_Na21 = 187,
        k_p_N15 = 188,
        k_He4_Ar37 = 189,
        k_p_K41 = 190,
        k_p_Co59 = 191,
        k_He4_Cl37 = 192,
        k_He4_F19 = 193,
        k_p_Ar39 = 194,
        k_He4_S35 = 195,
        k_He4_P29 = 196,
        k_p_Fe52 = 197,
        k_He4_Sc49 = 198,
        k_p_Al26 = 199,
        k_He4_C14 = 200,
        k_p_Cr50 = 201,
        k_p_S32 = 202,
        k_p_V48 = 203,
        k_p_Zn64 = 204,
        k_p_He4 = 205,
        k_p_Cu64 = 206,
        k_p_Sc46 = 207,
        k_p_Ti51 = 208,
        k_He4_Co54 = 209,
        k_p_Ni62 = 210,
        k_He4_Al26 = 211,
        k_He4_Ne21 = 212,
        k_He4_S32 = 213,
        k_p_Cu59 = 214,
        k_p_Ni57 = 215,
        k_d_Be7 = 216,
        k_He4_Ti51 = 217,
        k_He4_Ca48 = 218,
        k_p_Mn53 = 219,
        k_p_Zn66 = 220,
        k_He4_Zn59 = 221,
        k_p_Be9 = 222,
        k_He4_P30 = 223,
        k_He4_Co56 = 224,
        k_He4_Cu59 = 225,
        k_He4_Si28 = 226,
        k_He4_Fe54 = 227,
        k_He4_Ni57 = 228,
        k_p_Ti48 = 229,
        k_He4_Ca44 = 230,
        k_He4_Mg26 = 231,
        k_He4_Be8 = 232,
        k_p_Ca46 = 233,
        k_He3_He3 = 234,
        k_p_Ca41 = 235,
        k_p_Co57 = 236,
        k_p_K39 = 237,
        k_He3_He4 = 238,
        k_p_Fe55 = 239,
        k_He4_Cl35 = 240,
        k_p_C13 = 241,
        k_d_Li7 = 242,
        k_He4_S33 = 243,
        k_p_B11 = 244,
        k_He4_Ni59 = 245,
        k_p_Cr53 = 246,
        k_p_Ti50 = 247,
        k_p_V51 = 248,
        k_He4_Sc47 = 249,
        k_He4_Al25 = 250,
        k_He4_C12 = 251,
        k_He4_Ca41 = 252,
        k_He4_K39 = 253,
        k_p_Ti46 = 254,
        k_p_Cu62 = 255,
        k_p_Sc44 = 256,
        k_He4_Ti50 = 257,
        k_p_Ni60 = 258,
        k_p_Fe57 = 259,
        k_p_O18 = 260,
        k_p_Si32 = 261,
        k_p_Fe58 = 262,
        k_p_Ar40 = 263,
        k_He4_Cr48 = 264,
        k_He4_Ti45 = 265,
        k_He4_V46 = 266,
        k_He4_O17 = 267,
        k_p_Ni56 = 268,
        k_He4_N15 = 269,
        k_He4_K41 = 270,
        k_He4_Sc44 = 271,
        k_p_Mn51 = 272,
        k_He4_Ar39 = 273,
        k_He4_Fe57 = 274,
        k_He4_Co59 = 275,
        k_He4_Cu62 = 276,
        k_p_Be7 = 277,
        k_p_Cr51 = 278,
        k_p_Na23 = 279,
        k_p_V49 = 280,
        k_He4_Cu57 = 281,
        k_He4_Fe52 = 282,
        k_d_Li6 = 283,
        k_p_Ni63 = 284,
        k_He4_K37 = 285,
        k_He4_Cr50 = 286,
        k_p_N13 = 287,
        k_He4_Mn51 = 288,
        k_He4_V48 = 289,
        k_p_B10 = 290,
        k_He4_Be7 = 291,
        k_He4_He4 = 292,
        k_He4_Na23 = 293,
        k_He4_Sc46 = 294,
        k_He4_Ni62 = 295,
        k_C12_Ne20 = 296,
        k_p_Mn54 = 297,
        k_p_Ca40 = 298,
        k_p_Ar36 = 299,
        NumScreenPairs = k_p_Ar36
    };

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "Be7_to_Li7_electron_capture_reaclib",  // 1,
        "C14_to_N14_beta_neg_reaclib",  // 2,
        "N13_to_C13_beta_pos_reaclib",  // 3,
        "O14_to_N14_beta_pos_reaclib",  // 4,
        "O15_to_N15_beta_pos_reaclib",  // 5,
        "F17_to_O17_beta_pos_reaclib",  // 6,
        "F18_to_O18_beta_pos_reaclib",  // 7,
        "Ne18_to_F18_beta_pos_reaclib",  // 8,
        "Ne19_to_F19_beta_pos_reaclib",  // 9,
        "Zn59_to_Cu59_beta_pos_reaclib",  // 10,
        "Ge63_to_Ga63_beta_pos_reaclib",  // 11,
        "B8_to_He4_He4_beta_pos_reaclib",  // 12,
        "Zn59_to_p_Ni58_beta_pos_reaclib",  // 13,
        "p_n_to_d_reaclib",  // 14,
        "p_p_to_d_beta_pos_reaclib",  // 15,
        "p_p_to_d_electron_capture_reaclib",  // 16,
        "d_p_to_He3_reaclib",  // 17,
        "d_d_to_He4_reaclib",  // 18,
        "d_He4_to_Li6_reaclib",  // 19,
        "He3_n_to_He4_reaclib",  // 20,
        "He3_p_to_He4_beta_pos_reaclib",  // 21,
        "He3_He4_to_Be7_reaclib",  // 22,
        "Li6_n_to_Li7_reaclib",  // 23,
        "Li6_p_to_Be7_reaclib",  // 24,
        "Li6_He4_to_B10_reaclib",  // 25,
        "Li7_He4_to_B11_reaclib",  // 26,
        "Be7_p_to_B8_reaclib",  // 27,
        "Be9_p_to_B10_reaclib",  // 28,
        "B10_n_to_B11_reaclib",  // 29,
        "B11_p_to_C12_reaclib",  // 30,
        "C12_n_to_C13_reaclib",  // 31,
        "C12_p_to_N13_reaclib",  // 32,
        "C12_He4_to_O16_reaclib",  // 33,
        "C13_n_to_C14_reaclib",  // 34,
        "C13_p_to_N14_reaclib",  // 35,
        "C14_p_to_N15_reaclib",  // 36,
        "C14_He4_to_O18_reaclib",  // 37,
        "N13_n_to_N14_reaclib",  // 38,
        "N13_p_to_O14_reaclib",  // 39,
        "N14_n_to_N15_reaclib",  // 40,
        "N14_p_to_O15_reaclib",  // 41,
        "N14_He4_to_F18_reaclib",  // 42,
        "N15_p_to_O16_reaclib",  // 43,
        "N15_He4_to_F19_reaclib",  // 44,
        "O14_n_to_O15_reaclib",  // 45,
        "O14_He4_to_Ne18_reaclib",  // 46,
        "O15_n_to_O16_reaclib",  // 47,
        "O15_He4_to_Ne19_reaclib",  // 48,
        "O16_n_to_O17_reaclib",  // 49,
        "O16_p_to_F17_reaclib",  // 50,
        "O16_He4_to_Ne20_reaclib",  // 51,
        "O17_n_to_O18_reaclib",  // 52,
        "O17_p_to_F18_reaclib",  // 53,
        "O17_He4_to_Ne21_reaclib",  // 54,
        "O18_p_to_F19_reaclib",  // 55,
        "O18_He4_to_Ne22_reaclib",  // 56,
        "F17_n_to_F18_reaclib",  // 57,
        "F17_p_to_Ne18_reaclib",  // 58,
        "F17_He4_to_Na21_reaclib",  // 59,
        "F18_n_to_F19_reaclib",  // 60,
        "F18_p_to_Ne19_reaclib",  // 61,
        "F18_He4_to_Na22_reaclib",  // 62,
        "F19_p_to_Ne20_reaclib",  // 63,
        "F19_He4_to_Na23_reaclib",  // 64,
        "Ne18_n_to_Ne19_reaclib",  // 65,
        "Ne19_n_to_Ne20_reaclib",  // 66,
        "Ne19_He4_to_Mg23_reaclib",  // 67,
        "Ne20_n_to_Ne21_reaclib",  // 68,
        "Ne20_p_to_Na21_reaclib",  // 69,
        "Ne20_He4_to_Mg24_reaclib",  // 70,
        "Ne21_n_to_Ne22_reaclib",  // 71,
        "Ne21_p_to_Na22_reaclib",  // 72,
        "Ne21_He4_to_Mg25_reaclib",  // 73,
        "Ne22_p_to_Na23_reaclib",  // 74,
        "Ne22_He4_to_Mg26_reaclib",  // 75,
        "Na21_n_to_Na22_reaclib",  // 76,
        "Na21_He4_to_Al25_reaclib",  // 77,
        "Na22_n_to_Na23_reaclib",  // 78,
        "Na22_p_to_Mg23_reaclib",  // 79,
        "Na22_He4_to_Al26_reaclib",  // 80,
        "Na23_p_to_Mg24_reaclib",  // 81,
        "Na23_He4_to_Al27_reaclib",  // 82,
        "Mg23_n_to_Mg24_reaclib",  // 83,
        "Mg24_n_to_Mg25_reaclib",  // 84,
        "Mg24_p_to_Al25_reaclib",  // 85,
        "Mg24_He4_to_Si28_reaclib",  // 86,
        "Mg25_n_to_Mg26_reaclib",  // 87,
        "Mg25_p_to_Al26_reaclib",  // 88,
        "Mg25_He4_to_Si29_reaclib",  // 89,
        "Mg26_p_to_Al27_reaclib",  // 90,
        "Mg26_He4_to_Si30_reaclib",  // 91,
        "Al25_n_to_Al26_reaclib",  // 92,
        "Al25_He4_to_P29_reaclib",  // 93,
        "Al26_n_to_Al27_reaclib",  // 94,
        "Al26_He4_to_P30_reaclib",  // 95,
        "Al27_p_to_Si28_reaclib",  // 96,
        "Al27_He4_to_P31_reaclib",  // 97,
        "Si28_n_to_Si29_reaclib",  // 98,
        "Si28_p_to_P29_reaclib",  // 99,
        "Si28_He4_to_S32_reaclib",  // 100,
        "Si29_n_to_Si30_reaclib",  // 101,
        "Si29_p_to_P30_reaclib",  // 102,
        "Si29_He4_to_S33_reaclib",  // 103,
        "Si30_n_to_Si31_reaclib",  // 104,
        "Si30_p_to_P31_reaclib",  // 105,
        "Si30_He4_to_S34_reaclib",  // 106,
        "Si31_n_to_Si32_reaclib",  // 107,
        "Si31_p_to_P32_reaclib",  // 108,
        "Si31_He4_to_S35_reaclib",  // 109,
        "Si32_p_to_P33_reaclib",  // 110,
        "Si32_He4_to_S36_reaclib",  // 111,
        "P29_n_to_P30_reaclib",  // 112,
        "P29_He4_to_Cl33_reaclib",  // 113,
        "P30_n_to_P31_reaclib",  // 114,
        "P30_He4_to_Cl34_reaclib",  // 115,
        "P31_n_to_P32_reaclib",  // 116,
        "P31_p_to_S32_reaclib",  // 117,
        "P31_He4_to_Cl35_reaclib",  // 118,
        "P32_n_to_P33_reaclib",  // 119,
        "P32_p_to_S33_reaclib",  // 120,
        "P32_He4_to_Cl36_reaclib",  // 121,
        "P33_p_to_S34_reaclib",  // 122,
        "P33_He4_to_Cl37_reaclib",  // 123,
        "S32_n_to_S33_reaclib",  // 124,
        "S32_p_to_Cl33_reaclib",  // 125,
        "S32_He4_to_Ar36_reaclib",  // 126,
        "S33_n_to_S34_reaclib",  // 127,
        "S33_p_to_Cl34_reaclib",  // 128,
        "S33_He4_to_Ar37_reaclib",  // 129,
        "S34_n_to_S35_reaclib",  // 130,
        "S34_p_to_Cl35_reaclib",  // 131,
        "S34_He4_to_Ar38_reaclib",  // 132,
        "S35_n_to_S36_reaclib",  // 133,
        "S35_p_to_Cl36_reaclib",  // 134,
        "S35_He4_to_Ar39_reaclib",  // 135,
        "S36_p_to_Cl37_reaclib",  // 136,
        "S36_He4_to_Ar40_reaclib",  // 137,
        "Cl33_n_to_Cl34_reaclib",  // 138,
        "Cl33_He4_to_K37_reaclib",  // 139,
        "Cl34_n_to_Cl35_reaclib",  // 140,
        "Cl34_He4_to_K38_reaclib",  // 141,
        "Cl35_n_to_Cl36_reaclib",  // 142,
        "Cl35_p_to_Ar36_reaclib",  // 143,
        "Cl35_He4_to_K39_reaclib",  // 144,
        "Cl36_n_to_Cl37_reaclib",  // 145,
        "Cl36_p_to_Ar37_reaclib",  // 146,
        "Cl36_He4_to_K40_reaclib",  // 147,
        "Cl37_p_to_Ar38_reaclib",  // 148,
        "Cl37_He4_to_K41_reaclib",  // 149,
        "Ar36_n_to_Ar37_reaclib",  // 150,
        "Ar36_p_to_K37_reaclib",  // 151,
        "Ar36_He4_to_Ca40_reaclib",  // 152,
        "Ar37_n_to_Ar38_reaclib",  // 153,
        "Ar37_p_to_K38_reaclib",  // 154,
        "Ar37_He4_to_Ca41_reaclib",  // 155,
        "Ar38_n_to_Ar39_reaclib",  // 156,
        "Ar38_p_to_K39_reaclib",  // 157,
        "Ar38_He4_to_Ca42_reaclib",  // 158,
        "Ar39_n_to_Ar40_reaclib",  // 159,
        "Ar39_p_to_K40_reaclib",  // 160,
        "Ar39_He4_to_Ca43_reaclib",  // 161,
        "Ar40_p_to_K41_reaclib",  // 162,
        "Ar40_He4_to_Ca44_reaclib",  // 163,
        "K37_n_to_K38_reaclib",  // 164,
        "K38_n_to_K39_reaclib",  // 165,
        "K39_n_to_K40_reaclib",  // 166,
        "K39_p_to_Ca40_reaclib",  // 167,
        "K39_He4_to_Sc43_reaclib",  // 168,
        "K40_n_to_K41_reaclib",  // 169,
        "K40_p_to_Ca41_reaclib",  // 170,
        "K40_He4_to_Sc44_reaclib",  // 171,
        "K41_p_to_Ca42_reaclib",  // 172,
        "K41_He4_to_Sc45_reaclib",  // 173,
        "Ca40_n_to_Ca41_reaclib",  // 174,
        "Ca40_He4_to_Ti44_reaclib",  // 175,
        "Ca41_n_to_Ca42_reaclib",  // 176,
        "Ca41_He4_to_Ti45_reaclib",  // 177,
        "Ca42_n_to_Ca43_reaclib",  // 178,
        "Ca42_p_to_Sc43_reaclib",  // 179,
        "Ca42_He4_to_Ti46_reaclib",  // 180,
        "Ca43_n_to_Ca44_reaclib",  // 181,
        "Ca43_p_to_Sc44_reaclib",  // 182,
        "Ca43_He4_to_Ti47_reaclib",  // 183,
        "Ca44_n_to_Ca45_reaclib",  // 184,
        "Ca44_p_to_Sc45_reaclib",  // 185,
        "Ca44_He4_to_Ti48_reaclib",  // 186,
        "Ca45_n_to_Ca46_reaclib",  // 187,
        "Ca45_p_to_Sc46_reaclib",  // 188,
        "Ca45_He4_to_Ti49_reaclib",  // 189,
        "Ca46_n_to_Ca47_reaclib",  // 190,
        "Ca46_p_to_Sc47_reaclib",  // 191,
        "Ca46_He4_to_Ti50_reaclib",  // 192,
        "Ca47_n_to_Ca48_reaclib",  // 193,
        "Ca47_p_to_Sc48_reaclib",  // 194,
        "Ca47_He4_to_Ti51_reaclib",  // 195,
        "Ca48_p_to_Sc49_reaclib",  // 196,
        "Sc43_n_to_Sc44_reaclib",  // 197,
        "Sc43_p_to_Ti44_reaclib",  // 198,
        "Sc43_He4_to_V47_reaclib",  // 199,
        "Sc44_n_to_Sc45_reaclib",  // 200,
        "Sc44_p_to_Ti45_reaclib",  // 201,
        "Sc44_He4_to_V48_reaclib",  // 202,
        "Sc45_n_to_Sc46_reaclib",  // 203,
        "Sc45_p_to_Ti46_reaclib",  // 204,
        "Sc45_He4_to_V49_reaclib",  // 205,
        "Sc46_n_to_Sc47_reaclib",  // 206,
        "Sc46_p_to_Ti47_reaclib",  // 207,
        "Sc46_He4_to_V50_reaclib",  // 208,
        "Sc47_n_to_Sc48_reaclib",  // 209,
        "Sc47_p_to_Ti48_reaclib",  // 210,
        "Sc47_He4_to_V51_reaclib",  // 211,
        "Sc48_n_to_Sc49_reaclib",  // 212,
        "Sc48_p_to_Ti49_reaclib",  // 213,
        "Sc48_He4_to_V52_reaclib",  // 214,
        "Sc49_p_to_Ti50_reaclib",  // 215,
        "Ti44_n_to_Ti45_reaclib",  // 216,
        "Ti44_He4_to_Cr48_reaclib",  // 217,
        "Ti45_n_to_Ti46_reaclib",  // 218,
        "Ti45_p_to_V46_reaclib",  // 219,
        "Ti45_He4_to_Cr49_reaclib",  // 220,
        "Ti46_n_to_Ti47_reaclib",  // 221,
        "Ti46_p_to_V47_reaclib",  // 222,
        "Ti46_He4_to_Cr50_reaclib",  // 223,
        "Ti47_n_to_Ti48_reaclib",  // 224,
        "Ti47_p_to_V48_reaclib",  // 225,
        "Ti47_He4_to_Cr51_reaclib",  // 226,
        "Ti48_n_to_Ti49_reaclib",  // 227,
        "Ti48_p_to_V49_reaclib",  // 228,
        "Ti48_He4_to_Cr52_reaclib",  // 229,
        "Ti49_n_to_Ti50_reaclib",  // 230,
        "Ti49_p_to_V50_reaclib",  // 231,
        "Ti49_He4_to_Cr53_reaclib",  // 232,
        "Ti50_n_to_Ti51_reaclib",  // 233,
        "Ti50_p_to_V51_reaclib",  // 234,
        "Ti50_He4_to_Cr54_reaclib",  // 235,
        "Ti51_p_to_V52_reaclib",  // 236,
        "V46_n_to_V47_reaclib",  // 237,
        "V46_He4_to_Mn50_reaclib",  // 238,
        "V47_n_to_V48_reaclib",  // 239,
        "V47_p_to_Cr48_reaclib",  // 240,
        "V47_He4_to_Mn51_reaclib",  // 241,
        "V48_n_to_V49_reaclib",  // 242,
        "V48_p_to_Cr49_reaclib",  // 243,
        "V48_He4_to_Mn52_reaclib",  // 244,
        "V49_n_to_V50_reaclib",  // 245,
        "V49_p_to_Cr50_reaclib",  // 246,
        "V49_He4_to_Mn53_reaclib",  // 247,
        "V50_n_to_V51_reaclib",  // 248,
        "V50_p_to_Cr51_reaclib",  // 249,
        "V50_He4_to_Mn54_reaclib",  // 250,
        "V51_n_to_V52_reaclib",  // 251,
        "V51_p_to_Cr52_reaclib",  // 252,
        "V51_He4_to_Mn55_reaclib",  // 253,
        "V52_p_to_Cr53_reaclib",  // 254,
        "Cr48_n_to_Cr49_reaclib",  // 255,
        "Cr48_He4_to_Fe52_reaclib",  // 256,
        "Cr49_n_to_Cr50_reaclib",  // 257,
        "Cr49_p_to_Mn50_reaclib",  // 258,
        "Cr49_He4_to_Fe53_reaclib",  // 259,
        "Cr50_n_to_Cr51_reaclib",  // 260,
        "Cr50_p_to_Mn51_reaclib",  // 261,
        "Cr50_He4_to_Fe54_reaclib",  // 262,
        "Cr51_n_to_Cr52_reaclib",  // 263,
        "Cr51_p_to_Mn52_reaclib",  // 264,
        "Cr51_He4_to_Fe55_reaclib",  // 265,
        "Cr52_n_to_Cr53_reaclib",  // 266,
        "Cr52_p_to_Mn53_reaclib",  // 267,
        "Cr52_He4_to_Fe56_reaclib",  // 268,
        "Cr53_n_to_Cr54_reaclib",  // 269,
        "Cr53_p_to_Mn54_reaclib",  // 270,
        "Cr53_He4_to_Fe57_reaclib",  // 271,
        "Cr54_p_to_Mn55_reaclib",  // 272,
        "Cr54_He4_to_Fe58_reaclib",  // 273,
        "Mn50_n_to_Mn51_reaclib",  // 274,
        "Mn50_He4_to_Co54_reaclib",  // 275,
        "Mn51_n_to_Mn52_reaclib",  // 276,
        "Mn51_p_to_Fe52_reaclib",  // 277,
        "Mn51_He4_to_Co55_reaclib",  // 278,
        "Mn52_n_to_Mn53_reaclib",  // 279,
        "Mn52_p_to_Fe53_reaclib",  // 280,
        "Mn52_He4_to_Co56_reaclib",  // 281,
        "Mn53_n_to_Mn54_reaclib",  // 282,
        "Mn53_p_to_Fe54_reaclib",  // 283,
        "Mn53_He4_to_Co57_reaclib",  // 284,
        "Mn54_n_to_Mn55_reaclib",  // 285,
        "Mn54_p_to_Fe55_reaclib",  // 286,
        "Mn54_He4_to_Co58_reaclib",  // 287,
        "Mn55_p_to_Fe56_reaclib",  // 288,
        "Mn55_He4_to_Co59_reaclib",  // 289,
        "Fe52_n_to_Fe53_reaclib",  // 290,
        "Fe52_p_to_Co53_reaclib",  // 291,
        "Fe52_He4_to_Ni56_reaclib",  // 292,
        "Fe53_n_to_Fe54_reaclib",  // 293,
        "Fe53_p_to_Co54_reaclib",  // 294,
        "Fe53_He4_to_Ni57_reaclib",  // 295,
        "Fe54_n_to_Fe55_reaclib",  // 296,
        "Fe54_p_to_Co55_reaclib",  // 297,
        "Fe54_He4_to_Ni58_reaclib",  // 298,
        "Fe55_n_to_Fe56_reaclib",  // 299,
        "Fe55_p_to_Co56_reaclib",  // 300,
        "Fe55_He4_to_Ni59_reaclib",  // 301,
        "Fe56_n_to_Fe57_reaclib",  // 302,
        "Fe56_p_to_Co57_reaclib",  // 303,
        "Fe56_He4_to_Ni60_reaclib",  // 304,
        "Fe57_n_to_Fe58_reaclib",  // 305,
        "Fe57_p_to_Co58_reaclib",  // 306,
        "Fe57_He4_to_Ni61_reaclib",  // 307,
        "Fe58_p_to_Co59_reaclib",  // 308,
        "Fe58_He4_to_Ni62_reaclib",  // 309,
        "Co53_n_to_Co54_reaclib",  // 310,
        "Co53_He4_to_Cu57_reaclib",  // 311,
        "Co54_n_to_Co55_reaclib",  // 312,
        "Co54_He4_to_Cu58_reaclib",  // 313,
        "Co55_n_to_Co56_reaclib",  // 314,
        "Co55_p_to_Ni56_reaclib",  // 315,
        "Co55_He4_to_Cu59_reaclib",  // 316,
        "Co56_n_to_Co57_reaclib",  // 317,
        "Co56_p_to_Ni57_reaclib",  // 318,
        "Co56_He4_to_Cu60_reaclib",  // 319,
        "Co57_n_to_Co58_reaclib",  // 320,
        "Co57_p_to_Ni58_reaclib",  // 321,
        "Co57_He4_to_Cu61_reaclib",  // 322,
        "Co58_n_to_Co59_reaclib",  // 323,
        "Co58_p_to_Ni59_reaclib",  // 324,
        "Co58_He4_to_Cu62_reaclib",  // 325,
        "Co59_p_to_Ni60_reaclib",  // 326,
        "Co59_He4_to_Cu63_reaclib",  // 327,
        "Ni56_n_to_Ni57_reaclib",  // 328,
        "Ni56_p_to_Cu57_reaclib",  // 329,
        "Ni56_He4_to_Zn60_reaclib",  // 330,
        "Ni57_n_to_Ni58_reaclib",  // 331,
        "Ni57_p_to_Cu58_reaclib",  // 332,
        "Ni57_He4_to_Zn61_reaclib",  // 333,
        "Ni58_n_to_Ni59_reaclib",  // 334,
        "Ni58_p_to_Cu59_reaclib",  // 335,
        "Ni58_He4_to_Zn62_reaclib",  // 336,
        "Ni59_n_to_Ni60_reaclib",  // 337,
        "Ni59_p_to_Cu60_reaclib",  // 338,
        "Ni59_He4_to_Zn63_reaclib",  // 339,
        "Ni60_n_to_Ni61_reaclib",  // 340,
        "Ni60_p_to_Cu61_reaclib",  // 341,
        "Ni60_He4_to_Zn64_reaclib",  // 342,
        "Ni61_n_to_Ni62_reaclib",  // 343,
        "Ni61_p_to_Cu62_reaclib",  // 344,
        "Ni61_He4_to_Zn65_reaclib",  // 345,
        "Ni62_n_to_Ni63_reaclib",  // 346,
        "Ni62_p_to_Cu63_reaclib",  // 347,
        "Ni62_He4_to_Zn66_reaclib",  // 348,
        "Ni63_n_to_Ni64_reaclib",  // 349,
        "Ni63_p_to_Cu64_reaclib",  // 350,
        "Ni64_p_to_Cu65_reaclib",  // 351,
        "Cu57_n_to_Cu58_reaclib",  // 352,
        "Cu58_n_to_Cu59_reaclib",  // 353,
        "Cu58_p_to_Zn59_reaclib",  // 354,
        "Cu58_He4_to_Ga62_reaclib",  // 355,
        "Cu59_n_to_Cu60_reaclib",  // 356,
        "Cu59_p_to_Zn60_reaclib",  // 357,
        "Cu59_He4_to_Ga63_reaclib",  // 358,
        "Cu60_n_to_Cu61_reaclib",  // 359,
        "Cu60_p_to_Zn61_reaclib",  // 360,
        "Cu60_He4_to_Ga64_reaclib",  // 361,
        "Cu61_n_to_Cu62_reaclib",  // 362,
        "Cu61_p_to_Zn62_reaclib",  // 363,
        "Cu62_n_to_Cu63_reaclib",  // 364,
        "Cu62_p_to_Zn63_reaclib",  // 365,
        "Cu63_n_to_Cu64_reaclib",  // 366,
        "Cu63_p_to_Zn64_reaclib",  // 367,
        "Cu64_n_to_Cu65_reaclib",  // 368,
        "Cu64_p_to_Zn65_reaclib",  // 369,
        "Cu65_p_to_Zn66_reaclib",  // 370,
        "Zn59_n_to_Zn60_reaclib",  // 371,
        "Zn59_He4_to_Ge63_reaclib",  // 372,
        "Zn60_n_to_Zn61_reaclib",  // 373,
        "Zn60_He4_to_Ge64_reaclib",  // 374,
        "Zn61_n_to_Zn62_reaclib",  // 375,
        "Zn61_p_to_Ga62_reaclib",  // 376,
        "Zn62_n_to_Zn63_reaclib",  // 377,
        "Zn62_p_to_Ga63_reaclib",  // 378,
        "Zn63_n_to_Zn64_reaclib",  // 379,
        "Zn63_p_to_Ga64_reaclib",  // 380,
        "Zn64_n_to_Zn65_reaclib",  // 381,
        "Zn65_n_to_Zn66_reaclib",  // 382,
        "Ga62_n_to_Ga63_reaclib",  // 383,
        "Ga62_p_to_Ge63_reaclib",  // 384,
        "Ga63_n_to_Ga64_reaclib",  // 385,
        "Ga63_p_to_Ge64_reaclib",  // 386,
        "Ge63_n_to_Ge64_reaclib",  // 387,
        "d_d_to_n_He3_reaclib",  // 388,
        "He3_d_to_p_He4_reaclib",  // 389,
        "Li6_p_to_He4_He3_reaclib",  // 390,
        "Li6_d_to_n_Be7_reaclib",  // 391,
        "Li6_d_to_p_Li7_reaclib",  // 392,
        "Li7_p_to_He4_He4_reaclib",  // 393,
        "Be7_n_to_p_Li7_reaclib",  // 394,
        "Be7_n_to_He4_He4_reaclib",  // 395,
        "Be9_p_to_He4_Li6_reaclib",  // 396,
        "Be9_He4_to_n_C12_reaclib",  // 397,
        "B10_n_to_He4_Li7_reaclib",  // 398,
        "B10_p_to_He4_Be7_reaclib",  // 399,
        "B10_He4_to_n_N13_reaclib",  // 400,
        "B10_He4_to_p_C13_reaclib",  // 401,
        "B11_He4_to_n_N14_reaclib",  // 402,
        "B11_He4_to_p_C14_reaclib",  // 403,
        "C12_C12_to_p_Na23_reaclib",  // 404,
        "C12_C12_to_He4_Ne20_reaclib",  // 405,
        "C13_p_to_n_N13_reaclib",  // 406,
        "C13_d_to_n_N14_reaclib",  // 407,
        "C13_He4_to_n_O16_reaclib",  // 408,
        "C14_d_to_n_N15_reaclib",  // 409,
        "N13_He4_to_p_O16_reaclib",  // 410,
        "N14_n_to_p_C14_reaclib",  // 411,
        "N14_p_to_n_O14_reaclib",  // 412,
        "N14_He4_to_n_F17_reaclib",  // 413,
        "N15_p_to_n_O15_reaclib",  // 414,
        "N15_p_to_He4_C12_reaclib",  // 415,
        "O14_He4_to_p_F17_reaclib",  // 416,
        "O15_n_to_He4_C12_reaclib",  // 417,
        "O16_C12_to_p_Al27_reaclib",  // 418,
        "O16_C12_to_He4_Mg24_reaclib",  // 419,
        "O16_O16_to_p_P31_reaclib",  // 420,
        "O16_O16_to_He4_Si28_reaclib",  // 421,
        "O17_n_to_He4_C14_reaclib",  // 422,
        "O17_p_to_He4_N14_reaclib",  // 423,
        "O17_He4_to_n_Ne20_reaclib",  // 424,
        "O18_p_to_He4_N15_reaclib",  // 425,
        "O18_He4_to_n_Ne21_reaclib",  // 426,
        "F17_n_to_p_O17_reaclib",  // 427,
        "F18_n_to_p_O18_reaclib",  // 428,
        "F18_n_to_He4_N15_reaclib",  // 429,
        "F18_p_to_He4_O15_reaclib",  // 430,
        "F18_He4_to_p_Ne21_reaclib",  // 431,
        "F19_p_to_n_Ne19_reaclib",  // 432,
        "F19_p_to_He4_O16_reaclib",  // 433,
        "F19_He4_to_p_Ne22_reaclib",  // 434,
        "Ne18_n_to_p_F18_reaclib",  // 435,
        "Ne18_n_to_He4_O15_reaclib",  // 436,
        "Ne18_He4_to_p_Na21_reaclib",  // 437,
        "Ne19_n_to_He4_O16_reaclib",  // 438,
        "Ne19_He4_to_p_Na22_reaclib",  // 439,
        "Ne20_p_to_He4_F17_reaclib",  // 440,
        "Ne20_C12_to_p_P31_reaclib",  // 441,
        "Ne20_C12_to_He4_Si28_reaclib",  // 442,
        "Ne21_He4_to_n_Mg24_reaclib",  // 443,
        "Ne22_He4_to_n_Mg25_reaclib",  // 444,
        "Na21_n_to_p_Ne21_reaclib",  // 445,
        "Na21_n_to_He4_F18_reaclib",  // 446,
        "Na22_n_to_p_Ne22_reaclib",  // 447,
        "Na22_n_to_He4_F19_reaclib",  // 448,
        "Na22_He4_to_p_Mg25_reaclib",  // 449,
        "Na23_p_to_n_Mg23_reaclib",  // 450,
        "Na23_p_to_He4_Ne20_reaclib",  // 451,
        "Na23_He4_to_p_Mg26_reaclib",  // 452,
        "Mg23_n_to_He4_Ne20_reaclib",  // 453,
        "Mg23_n_to_C12_C12_reaclib",  // 454,
        "Mg23_He4_to_p_Al26_reaclib",  // 455,
        "Mg24_p_to_He4_Na21_reaclib",  // 456,
        "Mg25_He4_to_n_Si28_reaclib",  // 457,
        "Mg26_He4_to_n_Si29_reaclib",  // 458,
        "Al25_n_to_p_Mg25_reaclib",  // 459,
        "Al25_n_to_He4_Na22_reaclib",  // 460,
        "Al25_He4_to_p_Si28_reaclib",  // 461,
        "Al26_n_to_p_Mg26_reaclib",  // 462,
        "Al26_n_to_He4_Na23_reaclib",  // 463,
        "Al26_He4_to_p_Si29_reaclib",  // 464,
        "Al27_p_to_He4_Mg24_reaclib",  // 465,
        "Al27_He4_to_n_P30_reaclib",  // 466,
        "Al27_He4_to_p_Si30_reaclib",  // 467,
        "Si31_p_to_n_P31_reaclib",  // 468,
        "Si31_He4_to_n_S34_reaclib",  // 469,
        "P29_n_to_p_Si29_reaclib",  // 470,
        "P29_n_to_He4_Al26_reaclib",  // 471,
        "P29_He4_to_p_S32_reaclib",  // 472,
        "P30_n_to_p_Si30_reaclib",  // 473,
        "P30_He4_to_p_S33_reaclib",  // 474,
        "P31_p_to_He4_Si28_reaclib",  // 475,
        "P31_He4_to_p_S34_reaclib",  // 476,
        "P32_n_to_p_Si32_reaclib",  // 477,
        "P32_p_to_n_S32_reaclib",  // 478,
        "P32_p_to_He4_Si29_reaclib",  // 479,
        "P33_p_to_He4_Si30_reaclib",  // 480,
        "S32_n_to_He4_Si29_reaclib",  // 481,
        "S33_n_to_p_P33_reaclib",  // 482,
        "S33_n_to_He4_Si30_reaclib",  // 483,
        "S34_He4_to_n_Ar37_reaclib",  // 484,
        "S34_He4_to_p_Cl37_reaclib",  // 485,
        "S35_n_to_He4_Si32_reaclib",  // 486,
        "S35_p_to_He4_P32_reaclib",  // 487,
        "S35_He4_to_n_Ar38_reaclib",  // 488,
        "S36_p_to_He4_P33_reaclib",  // 489,
        "Cl33_n_to_p_S33_reaclib",  // 490,
        "Cl33_n_to_He4_P30_reaclib",  // 491,
        "Cl33_He4_to_p_Ar36_reaclib",  // 492,
        "Cl34_n_to_p_S34_reaclib",  // 493,
        "Cl34_n_to_He4_P31_reaclib",  // 494,
        "Cl34_He4_to_p_Ar37_reaclib",  // 495,
        "Cl35_n_to_p_S35_reaclib",  // 496,
        "Cl35_n_to_He4_P32_reaclib",  // 497,
        "Cl35_p_to_He4_S32_reaclib",  // 498,
        "Cl35_He4_to_p_Ar38_reaclib",  // 499,
        "Cl36_n_to_p_S36_reaclib",  // 500,
        "Cl36_n_to_He4_P33_reaclib",  // 501,
        "Cl36_p_to_He4_S33_reaclib",  // 502,
        "Ar36_n_to_p_Cl36_reaclib",  // 503,
        "Ar36_n_to_He4_S33_reaclib",  // 504,
        "Ar37_n_to_p_Cl37_reaclib",  // 505,
        "Ar39_n_to_He4_S36_reaclib",  // 506,
        "Ar39_p_to_He4_Cl36_reaclib",  // 507,
        "Ar40_p_to_He4_Cl37_reaclib",  // 508,
        "K37_n_to_p_Ar37_reaclib",  // 509,
        "K37_n_to_He4_Cl34_reaclib",  // 510,
        "K37_He4_to_p_Ca40_reaclib",  // 511,
        "K38_n_to_p_Ar38_reaclib",  // 512,
        "K38_n_to_He4_Cl35_reaclib",  // 513,
        "K38_He4_to_p_Ca41_reaclib",  // 514,
        "K39_n_to_p_Ar39_reaclib",  // 515,
        "K39_n_to_He4_Cl36_reaclib",  // 516,
        "K39_p_to_He4_Ar36_reaclib",  // 517,
        "K40_n_to_p_Ar40_reaclib",  // 518,
        "K40_n_to_He4_Cl37_reaclib",  // 519,
        "K40_p_to_n_Ca40_reaclib",  // 520,
        "K40_p_to_He4_Ar37_reaclib",  // 521,
        "K40_He4_to_p_Ca43_reaclib",  // 522,
        "K41_p_to_n_Ca41_reaclib",  // 523,
        "K41_p_to_He4_Ar38_reaclib",  // 524,
        "K41_He4_to_n_Sc44_reaclib",  // 525,
        "K41_He4_to_p_Ca44_reaclib",  // 526,
        "Ca40_n_to_He4_Ar37_reaclib",  // 527,
        "Ca41_n_to_He4_Ar38_reaclib",  // 528,
        "Ca42_n_to_He4_Ar39_reaclib",  // 529,
        "Ca42_p_to_He4_K39_reaclib",  // 530,
        "Ca43_n_to_He4_Ar40_reaclib",  // 531,
        "Ca43_He4_to_n_Ti46_reaclib",  // 532,
        "Ca45_He4_to_n_Ti48_reaclib",  // 533,
        "Ca47_p_to_n_Sc47_reaclib",  // 534,
        "Ca47_He4_to_n_Ti50_reaclib",  // 535,
        "Sc43_n_to_p_Ca43_reaclib",  // 536,
        "Sc43_n_to_He4_K40_reaclib",  // 537,
        "Sc43_p_to_He4_Ca40_reaclib",  // 538,
        "Sc43_He4_to_p_Ti46_reaclib",  // 539,
        "Sc44_n_to_p_Ca44_reaclib",  // 540,
        "Sc44_p_to_He4_Ca41_reaclib",  // 541,
        "Sc44_He4_to_p_Ti47_reaclib",  // 542,
        "Sc45_n_to_p_Ca45_reaclib",  // 543,
        "Sc45_p_to_He4_Ca42_reaclib",  // 544,
        "Sc45_He4_to_p_Ti48_reaclib",  // 545,
        "Sc46_n_to_p_Ca46_reaclib",  // 546,
        "Sc46_p_to_n_Ti46_reaclib",  // 547,
        "Sc46_p_to_He4_Ca43_reaclib",  // 548,
        "Sc46_He4_to_n_V49_reaclib",  // 549,
        "Sc46_He4_to_p_Ti49_reaclib",  // 550,
        "Sc47_p_to_He4_Ca44_reaclib",  // 551,
        "Sc47_He4_to_p_Ti50_reaclib",  // 552,
        "Sc48_n_to_p_Ca48_reaclib",  // 553,
        "Sc48_p_to_n_Ti48_reaclib",  // 554,
        "Sc48_p_to_He4_Ca45_reaclib",  // 555,
        "Sc48_He4_to_n_V51_reaclib",  // 556,
        "Sc48_He4_to_p_Ti51_reaclib",  // 557,
        "Sc49_p_to_n_Ti49_reaclib",  // 558,
        "Sc49_p_to_He4_Ca46_reaclib",  // 559,
        "Ti44_n_to_p_Sc44_reaclib",  // 560,
        "Ti44_n_to_He4_Ca41_reaclib",  // 561,
        "Ti44_He4_to_p_V47_reaclib",  // 562,
        "Ti45_n_to_p_Sc45_reaclib",  // 563,
        "Ti45_n_to_He4_Ca42_reaclib",  // 564,
        "Ti45_He4_to_p_V48_reaclib",  // 565,
        "Ti47_n_to_p_Sc47_reaclib",  // 566,
        "Ti47_n_to_He4_Ca44_reaclib",  // 567,
        "Ti49_n_to_He4_Ca46_reaclib",  // 568,
        "Ti49_He4_to_n_Cr52_reaclib",  // 569,
        "Ti51_n_to_He4_Ca48_reaclib",  // 570,
        "Ti51_p_to_n_V51_reaclib",  // 571,
        "Ti51_He4_to_n_Cr54_reaclib",  // 572,
        "V46_n_to_p_Ti46_reaclib",  // 573,
        "V46_n_to_He4_Sc43_reaclib",  // 574,
        "V46_He4_to_p_Cr49_reaclib",  // 575,
        "V47_n_to_p_Ti47_reaclib",  // 576,
        "V47_n_to_He4_Sc44_reaclib",  // 577,
        "V47_He4_to_p_Cr50_reaclib",  // 578,
        "V48_n_to_p_Ti48_reaclib",  // 579,
        "V48_n_to_He4_Sc45_reaclib",  // 580,
        "V48_He4_to_p_Cr51_reaclib",  // 581,
        "V49_n_to_p_Ti49_reaclib",  // 582,
        "V49_p_to_He4_Ti46_reaclib",  // 583,
        "V49_He4_to_p_Cr52_reaclib",  // 584,
        "V50_n_to_p_Ti50_reaclib",  // 585,
        "V50_n_to_He4_Sc47_reaclib",  // 586,
        "V50_p_to_n_Cr50_reaclib",  // 587,
        "V50_p_to_He4_Ti47_reaclib",  // 588,
        "V50_He4_to_p_Cr53_reaclib",  // 589,
        "V51_p_to_He4_Ti48_reaclib",  // 590,
        "V52_n_to_He4_Sc49_reaclib",  // 591,
        "V52_p_to_n_Cr52_reaclib",  // 592,
        "V52_p_to_He4_Ti49_reaclib",  // 593,
        "V52_He4_to_n_Mn55_reaclib",  // 594,
        "Cr48_n_to_p_V48_reaclib",  // 595,
        "Cr48_n_to_He4_Ti45_reaclib",  // 596,
        "Cr48_He4_to_p_Mn51_reaclib",  // 597,
        "Cr49_n_to_p_V49_reaclib",  // 598,
        "Cr49_n_to_He4_Ti46_reaclib",  // 599,
        "Cr49_He4_to_p_Mn52_reaclib",  // 600,
        "Cr50_n_to_He4_Ti47_reaclib",  // 601,
        "Cr51_n_to_p_V51_reaclib",  // 602,
        "Cr51_n_to_He4_Ti48_reaclib",  // 603,
        "Cr53_n_to_He4_Ti50_reaclib",  // 604,
        "Cr54_p_to_He4_V51_reaclib",  // 605,
        "Mn50_n_to_p_Cr50_reaclib",  // 606,
        "Mn50_n_to_He4_V47_reaclib",  // 607,
        "Mn50_He4_to_p_Fe53_reaclib",  // 608,
        "Mn51_n_to_p_Cr51_reaclib",  // 609,
        "Mn51_n_to_He4_V48_reaclib",  // 610,
        "Mn51_He4_to_p_Fe54_reaclib",  // 611,
        "Mn52_n_to_p_Cr52_reaclib",  // 612,
        "Mn52_n_to_He4_V49_reaclib",  // 613,
        "Mn52_He4_to_p_Fe55_reaclib",  // 614,
        "Mn53_n_to_p_Cr53_reaclib",  // 615,
        "Mn53_n_to_He4_V50_reaclib",  // 616,
        "Mn53_p_to_He4_Cr50_reaclib",  // 617,
        "Mn53_He4_to_p_Fe56_reaclib",  // 618,
        "Mn54_n_to_p_Cr54_reaclib",  // 619,
        "Mn54_n_to_He4_V51_reaclib",  // 620,
        "Mn54_p_to_He4_Cr51_reaclib",  // 621,
        "Mn55_p_to_He4_Cr52_reaclib",  // 622,
        "Fe52_n_to_p_Mn52_reaclib",  // 623,
        "Fe52_n_to_He4_Cr49_reaclib",  // 624,
        "Fe52_He4_to_p_Co55_reaclib",  // 625,
        "Fe53_n_to_p_Mn53_reaclib",  // 626,
        "Fe53_n_to_He4_Cr50_reaclib",  // 627,
        "Fe53_He4_to_p_Co56_reaclib",  // 628,
        "Fe54_n_to_p_Mn54_reaclib",  // 629,
        "Fe54_n_to_He4_Cr51_reaclib",  // 630,
        "Fe55_n_to_p_Mn55_reaclib",  // 631,
        "Fe55_n_to_He4_Cr52_reaclib",  // 632,
        "Fe56_n_to_He4_Cr53_reaclib",  // 633,
        "Fe57_n_to_He4_Cr54_reaclib",  // 634,
        "Fe57_p_to_He4_Mn54_reaclib",  // 635,
        "Fe58_p_to_He4_Mn55_reaclib",  // 636,
        "Co53_n_to_p_Fe53_reaclib",  // 637,
        "Co53_n_to_He4_Mn50_reaclib",  // 638,
        "Co53_He4_to_p_Ni56_reaclib",  // 639,
        "Co54_n_to_p_Fe54_reaclib",  // 640,
        "Co54_n_to_He4_Mn51_reaclib",  // 641,
        "Co54_He4_to_p_Ni57_reaclib",  // 642,
        "Co55_n_to_p_Fe55_reaclib",  // 643,
        "Co55_n_to_He4_Mn52_reaclib",  // 644,
        "Co55_He4_to_p_Ni58_reaclib",  // 645,
        "Co56_n_to_p_Fe56_reaclib",  // 646,
        "Co56_n_to_He4_Mn53_reaclib",  // 647,
        "Co56_He4_to_p_Ni59_reaclib",  // 648,
        "Co57_n_to_p_Fe57_reaclib",  // 649,
        "Co57_n_to_He4_Mn54_reaclib",  // 650,
        "Co57_p_to_He4_Fe54_reaclib",  // 651,
        "Co57_He4_to_p_Ni60_reaclib",  // 652,
        "Co58_n_to_p_Fe58_reaclib",  // 653,
        "Co58_n_to_He4_Mn55_reaclib",  // 654,
        "Co58_p_to_He4_Fe55_reaclib",  // 655,
        "Co59_p_to_He4_Fe56_reaclib",  // 656,
        "Ni56_n_to_p_Co56_reaclib",  // 657,
        "Ni56_n_to_He4_Fe53_reaclib",  // 658,
        "Ni57_n_to_p_Co57_reaclib",  // 659,
        "Ni57_n_to_He4_Fe54_reaclib",  // 660,
        "Ni58_n_to_p_Co58_reaclib",  // 661,
        "Ni58_n_to_He4_Fe55_reaclib",  // 662,
        "Ni59_n_to_p_Co59_reaclib",  // 663,
        "Ni59_n_to_He4_Fe56_reaclib",  // 664,
        "Ni60_n_to_He4_Fe57_reaclib",  // 665,
        "Ni61_n_to_He4_Fe58_reaclib",  // 666,
        "Ni61_p_to_He4_Co58_reaclib",  // 667,
        "Ni62_p_to_He4_Co59_reaclib",  // 668,
        "Cu57_n_to_p_Ni57_reaclib",  // 669,
        "Cu57_n_to_He4_Co54_reaclib",  // 670,
        "Cu57_He4_to_p_Zn60_reaclib",  // 671,
        "Cu58_n_to_p_Ni58_reaclib",  // 672,
        "Cu58_n_to_He4_Co55_reaclib",  // 673,
        "Cu59_n_to_p_Ni59_reaclib",  // 674,
        "Cu59_n_to_He4_Co56_reaclib",  // 675,
        "Cu59_p_to_He4_Ni56_reaclib",  // 676,
        "Cu60_n_to_p_Ni60_reaclib",  // 677,
        "Cu60_n_to_He4_Co57_reaclib",  // 678,
        "Cu60_p_to_He4_Ni57_reaclib",  // 679,
        "Cu61_n_to_p_Ni61_reaclib",  // 680,
        "Cu61_n_to_He4_Co58_reaclib",  // 681,
        "Cu61_p_to_He4_Ni58_reaclib",  // 682,
        "Cu62_n_to_p_Ni62_reaclib",  // 683,
        "Cu62_n_to_He4_Co59_reaclib",  // 684,
        "Cu62_p_to_He4_Ni59_reaclib",  // 685,
        "Cu63_n_to_p_Ni63_reaclib",  // 686,
        "Cu63_p_to_He4_Ni60_reaclib",  // 687,
        "Cu64_n_to_p_Ni64_reaclib",  // 688,
        "Cu64_p_to_He4_Ni61_reaclib",  // 689,
        "Cu65_p_to_He4_Ni62_reaclib",  // 690,
        "Zn59_n_to_p_Cu59_reaclib",  // 691,
        "Zn59_n_to_He4_Ni56_reaclib",  // 692,
        "Zn60_n_to_p_Cu60_reaclib",  // 693,
        "Zn60_n_to_He4_Ni57_reaclib",  // 694,
        "Zn61_n_to_p_Cu61_reaclib",  // 695,
        "Zn61_n_to_He4_Ni58_reaclib",  // 696,
        "Zn61_p_to_He4_Cu58_reaclib",  // 697,
        "Zn62_n_to_p_Cu62_reaclib",  // 698,
        "Zn62_n_to_He4_Ni59_reaclib",  // 699,
        "Zn62_p_to_He4_Cu59_reaclib",  // 700,
        "Zn63_n_to_p_Cu63_reaclib",  // 701,
        "Zn63_n_to_He4_Ni60_reaclib",  // 702,
        "Zn63_p_to_He4_Cu60_reaclib",  // 703,
        "Zn64_n_to_p_Cu64_reaclib",  // 704,
        "Zn64_n_to_He4_Ni61_reaclib",  // 705,
        "Zn64_p_to_He4_Cu61_reaclib",  // 706,
        "Zn65_n_to_p_Cu65_reaclib",  // 707,
        "Zn65_n_to_He4_Ni62_reaclib",  // 708,
        "Zn65_p_to_He4_Cu62_reaclib",  // 709,
        "Zn66_n_to_He4_Ni63_reaclib",  // 710,
        "Zn66_p_to_He4_Cu63_reaclib",  // 711,
        "Ga62_n_to_p_Zn62_reaclib",  // 712,
        "Ga62_n_to_He4_Cu59_reaclib",  // 713,
        "Ga62_p_to_He4_Zn59_reaclib",  // 714,
        "Ga63_n_to_p_Zn63_reaclib",  // 715,
        "Ga63_n_to_He4_Cu60_reaclib",  // 716,
        "Ga63_p_to_He4_Zn60_reaclib",  // 717,
        "Ga64_n_to_p_Zn64_reaclib",  // 718,
        "Ga64_n_to_He4_Cu61_reaclib",  // 719,
        "Ga64_p_to_He4_Zn61_reaclib",  // 720,
        "Ge63_n_to_p_Ga63_reaclib",  // 721,
        "Ge63_n_to_He4_Zn60_reaclib",  // 722,
        "Ge64_n_to_p_Ga64_reaclib",  // 723,
        "Ge64_n_to_He4_Zn61_reaclib",  // 724,
        "He3_He3_to_p_p_He4_reaclib",  // 725,
        "Li7_d_to_n_He4_He4_reaclib",  // 726,
        "Be7_d_to_p_He4_He4_reaclib",  // 727,
        "Be9_p_to_d_He4_He4_reaclib",  // 728,
        "B8_n_to_p_He4_He4_reaclib",  // 729,
        "B11_p_to_3He4_reaclib",  // 730,
        "Li7_He3_to_n_p_He4_He4_reaclib",  // 731,
        "Be7_He3_to_p_p_He4_He4_reaclib",  // 732,
        "Be9_p_to_n_p_He4_He4_reaclib",  // 733,
        "He4_p_n_to_Li6_reaclib",  // 734,
        "He4_He4_n_to_Be9_reaclib",  // 735,
        "He4_He4_He4_to_C12_reaclib",  // 736,
        "p_p_n_to_p_d_reaclib",  // 737,
        "Na21_to_Ne21_electron_capture_weaktab",  // 738,
        "Ne21_to_Na21_beta_neg_weaktab",  // 739,
        "Na22_to_Ne22_electron_capture_weaktab",  // 740,
        "Ne22_to_Na22_beta_neg_weaktab",  // 741,
        "Mg23_to_Na23_electron_capture_weaktab",  // 742,
        "Na23_to_Mg23_beta_neg_weaktab",  // 743,
        "Al25_to_Mg25_electron_capture_weaktab",  // 744,
        "Mg25_to_Al25_beta_neg_weaktab",  // 745,
        "Al26_to_Mg26_electron_capture_weaktab",  // 746,
        "Mg26_to_Al26_beta_neg_weaktab",  // 747,
        "P29_to_Si29_electron_capture_weaktab",  // 748,
        "Si29_to_P29_beta_neg_weaktab",  // 749,
        "P30_to_Si30_electron_capture_weaktab",  // 750,
        "Si30_to_P30_beta_neg_weaktab",  // 751,
        "P31_to_Si31_electron_capture_weaktab",  // 752,
        "Si31_to_P31_beta_neg_weaktab",  // 753,
        "P32_to_S32_beta_neg_weaktab",  // 754,
        "P32_to_Si32_electron_capture_weaktab",  // 755,
        "S32_to_P32_electron_capture_weaktab",  // 756,
        "Si32_to_P32_beta_neg_weaktab",  // 757,
        "Cl33_to_S33_electron_capture_weaktab",  // 758,
        "P33_to_S33_beta_neg_weaktab",  // 759,
        "S33_to_Cl33_beta_neg_weaktab",  // 760,
        "S33_to_P33_electron_capture_weaktab",  // 761,
        "Cl34_to_S34_electron_capture_weaktab",  // 762,
        "S34_to_Cl34_beta_neg_weaktab",  // 763,
        "Cl35_to_S35_electron_capture_weaktab",  // 764,
        "S35_to_Cl35_beta_neg_weaktab",  // 765,
        "Ar36_to_Cl36_electron_capture_weaktab",  // 766,
        "Cl36_to_Ar36_beta_neg_weaktab",  // 767,
        "Cl36_to_S36_electron_capture_weaktab",  // 768,
        "S36_to_Cl36_beta_neg_weaktab",  // 769,
        "Ar37_to_Cl37_electron_capture_weaktab",  // 770,
        "Ar37_to_K37_beta_neg_weaktab",  // 771,
        "Cl37_to_Ar37_beta_neg_weaktab",  // 772,
        "K37_to_Ar37_electron_capture_weaktab",  // 773,
        "Ar38_to_K38_beta_neg_weaktab",  // 774,
        "K38_to_Ar38_electron_capture_weaktab",  // 775,
        "Ar39_to_K39_beta_neg_weaktab",  // 776,
        "K39_to_Ar39_electron_capture_weaktab",  // 777,
        "Ar40_to_K40_beta_neg_weaktab",  // 778,
        "Ca40_to_K40_electron_capture_weaktab",  // 779,
        "K40_to_Ar40_electron_capture_weaktab",  // 780,
        "K40_to_Ca40_beta_neg_weaktab",  // 781,
        "Ca41_to_K41_electron_capture_weaktab",  // 782,
        "K41_to_Ca41_beta_neg_weaktab",  // 783,
        "Ca43_to_Sc43_beta_neg_weaktab",  // 784,
        "Sc43_to_Ca43_electron_capture_weaktab",  // 785,
        "Ca44_to_Sc44_beta_neg_weaktab",  // 786,
        "Sc44_to_Ca44_electron_capture_weaktab",  // 787,
        "Sc44_to_Ti44_beta_neg_weaktab",  // 788,
        "Ti44_to_Sc44_electron_capture_weaktab",  // 789,
        "Co53_to_Fe53_electron_capture_weaktab",  // 790,
        "Fe53_to_Co53_beta_neg_weaktab",  // 791,
        "Cu57_to_Ni57_electron_capture_weaktab",  // 792,
        "Ni57_to_Cu57_beta_neg_weaktab",  // 793,
        "Ca45_to_Sc45_beta_neg_weaktab",  // 794,
        "Sc45_to_Ca45_electron_capture_weaktab",  // 795,
        "Sc45_to_Ti45_beta_neg_weaktab",  // 796,
        "Ti45_to_Sc45_electron_capture_weaktab",  // 797,
        "Ca46_to_Sc46_beta_neg_weaktab",  // 798,
        "Sc46_to_Ca46_electron_capture_weaktab",  // 799,
        "Sc46_to_Ti46_beta_neg_weaktab",  // 800,
        "Ti46_to_Sc46_electron_capture_weaktab",  // 801,
        "Ti46_to_V46_beta_neg_weaktab",  // 802,
        "V46_to_Ti46_electron_capture_weaktab",  // 803,
        "Ca47_to_Sc47_beta_neg_weaktab",  // 804,
        "Sc47_to_Ca47_electron_capture_weaktab",  // 805,
        "Sc47_to_Ti47_beta_neg_weaktab",  // 806,
        "Ti47_to_Sc47_electron_capture_weaktab",  // 807,
        "Ti47_to_V47_beta_neg_weaktab",  // 808,
        "V47_to_Ti47_electron_capture_weaktab",  // 809,
        "Ca48_to_Sc48_beta_neg_weaktab",  // 810,
        "Cr48_to_V48_electron_capture_weaktab",  // 811,
        "Sc48_to_Ca48_electron_capture_weaktab",  // 812,
        "Sc48_to_Ti48_beta_neg_weaktab",  // 813,
        "Ti48_to_Sc48_electron_capture_weaktab",  // 814,
        "Ti48_to_V48_beta_neg_weaktab",  // 815,
        "V48_to_Cr48_beta_neg_weaktab",  // 816,
        "V48_to_Ti48_electron_capture_weaktab",  // 817,
        "Cr49_to_V49_electron_capture_weaktab",  // 818,
        "Sc49_to_Ti49_beta_neg_weaktab",  // 819,
        "Ti49_to_Sc49_electron_capture_weaktab",  // 820,
        "Ti49_to_V49_beta_neg_weaktab",  // 821,
        "V49_to_Cr49_beta_neg_weaktab",  // 822,
        "V49_to_Ti49_electron_capture_weaktab",  // 823,
        "Cr50_to_Mn50_beta_neg_weaktab",  // 824,
        "Cr50_to_V50_electron_capture_weaktab",  // 825,
        "Mn50_to_Cr50_electron_capture_weaktab",  // 826,
        "Ti50_to_V50_beta_neg_weaktab",  // 827,
        "V50_to_Cr50_beta_neg_weaktab",  // 828,
        "V50_to_Ti50_electron_capture_weaktab",  // 829,
        "Cr51_to_Mn51_beta_neg_weaktab",  // 830,
        "Cr51_to_V51_electron_capture_weaktab",  // 831,
        "Mn51_to_Cr51_electron_capture_weaktab",  // 832,
        "Ti51_to_V51_beta_neg_weaktab",  // 833,
        "V51_to_Cr51_beta_neg_weaktab",  // 834,
        "V51_to_Ti51_electron_capture_weaktab",  // 835,
        "Cr52_to_Mn52_beta_neg_weaktab",  // 836,
        "Cr52_to_V52_electron_capture_weaktab",  // 837,
        "Fe52_to_Mn52_electron_capture_weaktab",  // 838,
        "Mn52_to_Cr52_electron_capture_weaktab",  // 839,
        "Mn52_to_Fe52_beta_neg_weaktab",  // 840,
        "V52_to_Cr52_beta_neg_weaktab",  // 841,
        "Cr53_to_Mn53_beta_neg_weaktab",  // 842,
        "Fe53_to_Mn53_electron_capture_weaktab",  // 843,
        "Mn53_to_Cr53_electron_capture_weaktab",  // 844,
        "Mn53_to_Fe53_beta_neg_weaktab",  // 845,
        "Co54_to_Fe54_electron_capture_weaktab",  // 846,
        "Cr54_to_Mn54_beta_neg_weaktab",  // 847,
        "Fe54_to_Co54_beta_neg_weaktab",  // 848,
        "Fe54_to_Mn54_electron_capture_weaktab",  // 849,
        "Mn54_to_Cr54_electron_capture_weaktab",  // 850,
        "Mn54_to_Fe54_beta_neg_weaktab",  // 851,
        "Co55_to_Fe55_electron_capture_weaktab",  // 852,
        "Fe55_to_Co55_beta_neg_weaktab",  // 853,
        "Fe55_to_Mn55_electron_capture_weaktab",  // 854,
        "Mn55_to_Fe55_beta_neg_weaktab",  // 855,
        "Co56_to_Fe56_electron_capture_weaktab",  // 856,
        "Co56_to_Ni56_beta_neg_weaktab",  // 857,
        "Fe56_to_Co56_beta_neg_weaktab",  // 858,
        "Ni56_to_Co56_electron_capture_weaktab",  // 859,
        "Co57_to_Fe57_electron_capture_weaktab",  // 860,
        "Co57_to_Ni57_beta_neg_weaktab",  // 861,
        "Fe57_to_Co57_beta_neg_weaktab",  // 862,
        "Ni57_to_Co57_electron_capture_weaktab",  // 863,
        "Co58_to_Fe58_electron_capture_weaktab",  // 864,
        "Co58_to_Ni58_beta_neg_weaktab",  // 865,
        "Cu58_to_Ni58_electron_capture_weaktab",  // 866,
        "Fe58_to_Co58_beta_neg_weaktab",  // 867,
        "Ni58_to_Co58_electron_capture_weaktab",  // 868,
        "Ni58_to_Cu58_beta_neg_weaktab",  // 869,
        "Co59_to_Ni59_beta_neg_weaktab",  // 870,
        "Cu59_to_Ni59_electron_capture_weaktab",  // 871,
        "Ni59_to_Co59_electron_capture_weaktab",  // 872,
        "Ni59_to_Cu59_beta_neg_weaktab",  // 873,
        "Cu60_to_Ni60_electron_capture_weaktab",  // 874,
        "Cu60_to_Zn60_beta_neg_weaktab",  // 875,
        "Ni60_to_Cu60_beta_neg_weaktab",  // 876,
        "Zn60_to_Cu60_electron_capture_weaktab",  // 877,
        "Cu61_to_Ni61_electron_capture_weaktab",  // 878,
        "Cu61_to_Zn61_beta_neg_weaktab",  // 879,
        "Ni61_to_Cu61_beta_neg_weaktab",  // 880,
        "Zn61_to_Cu61_electron_capture_weaktab",  // 881,
        "Cu62_to_Ni62_electron_capture_weaktab",  // 882,
        "Cu62_to_Zn62_beta_neg_weaktab",  // 883,
        "Ga62_to_Zn62_electron_capture_weaktab",  // 884,
        "Ni62_to_Cu62_beta_neg_weaktab",  // 885,
        "Zn62_to_Cu62_electron_capture_weaktab",  // 886,
        "Zn62_to_Ga62_beta_neg_weaktab",  // 887,
        "Cu63_to_Ni63_electron_capture_weaktab",  // 888,
        "Cu63_to_Zn63_beta_neg_weaktab",  // 889,
        "Ga63_to_Zn63_electron_capture_weaktab",  // 890,
        "Ni63_to_Cu63_beta_neg_weaktab",  // 891,
        "Zn63_to_Cu63_electron_capture_weaktab",  // 892,
        "Zn63_to_Ga63_beta_neg_weaktab",  // 893,
        "Cu64_to_Ni64_electron_capture_weaktab",  // 894,
        "Cu64_to_Zn64_beta_neg_weaktab",  // 895,
        "Ga64_to_Ge64_beta_neg_weaktab",  // 896,
        "Ga64_to_Zn64_electron_capture_weaktab",  // 897,
        "Ge64_to_Ga64_electron_capture_weaktab",  // 898,
        "Ni64_to_Cu64_beta_neg_weaktab",  // 899,
        "Zn64_to_Cu64_electron_capture_weaktab",  // 900,
        "Zn64_to_Ga64_beta_neg_weaktab",  // 901,
        "Cu65_to_Zn65_beta_neg_weaktab",  // 902,
        "Zn65_to_Cu65_electron_capture_weaktab",  // 903,
        "n_to_p_beta_neg_weaktab",  // 904,
        "p_to_n_electron_capture_weaktab",  // 905,
        "d_to_n_p_derived",  // 906,
        "He3_to_p_d_derived",  // 907,
        "He4_to_n_He3_derived",  // 908,
        "He4_to_d_d_derived",  // 909,
        "Li6_to_He4_d_derived",  // 910,
        "Li7_to_n_Li6_derived",  // 911,
        "Be7_to_p_Li6_derived",  // 912,
        "Be7_to_He4_He3_derived",  // 913,
        "B8_to_p_Be7_derived",  // 914,
        "B10_to_p_Be9_derived",  // 915,
        "B10_to_He4_Li6_derived",  // 916,
        "B11_to_n_B10_derived",  // 917,
        "B11_to_He4_Li7_derived",  // 918,
        "C12_to_p_B11_derived",  // 919,
        "C13_to_n_C12_derived",  // 920,
        "C14_to_n_C13_derived",  // 921,
        "N13_to_p_C12_derived",  // 922,
        "N14_to_n_N13_derived",  // 923,
        "N14_to_p_C13_derived",  // 924,
        "N15_to_n_N14_derived",  // 925,
        "N15_to_p_C14_derived",  // 926,
        "O14_to_p_N13_derived",  // 927,
        "O15_to_n_O14_derived",  // 928,
        "O15_to_p_N14_derived",  // 929,
        "O16_to_n_O15_derived",  // 930,
        "O16_to_p_N15_derived",  // 931,
        "O16_to_He4_C12_derived",  // 932,
        "O17_to_n_O16_derived",  // 933,
        "O18_to_n_O17_derived",  // 934,
        "O18_to_He4_C14_derived",  // 935,
        "F17_to_p_O16_derived",  // 936,
        "F18_to_n_F17_derived",  // 937,
        "F18_to_p_O17_derived",  // 938,
        "F18_to_He4_N14_derived",  // 939,
        "F19_to_n_F18_derived",  // 940,
        "F19_to_p_O18_derived",  // 941,
        "F19_to_He4_N15_derived",  // 942,
        "Ne18_to_p_F17_derived",  // 943,
        "Ne18_to_He4_O14_derived",  // 944,
        "Ne19_to_n_Ne18_derived",  // 945,
        "Ne19_to_p_F18_derived",  // 946,
        "Ne19_to_He4_O15_derived",  // 947,
        "Ne20_to_n_Ne19_derived",  // 948,
        "Ne20_to_p_F19_derived",  // 949,
        "Ne20_to_He4_O16_derived",  // 950,
        "Ne21_to_n_Ne20_derived",  // 951,
        "Ne21_to_He4_O17_derived",  // 952,
        "Ne22_to_n_Ne21_derived",  // 953,
        "Ne22_to_He4_O18_derived",  // 954,
        "Na21_to_p_Ne20_derived",  // 955,
        "Na21_to_He4_F17_derived",  // 956,
        "Na22_to_n_Na21_derived",  // 957,
        "Na22_to_p_Ne21_derived",  // 958,
        "Na22_to_He4_F18_derived",  // 959,
        "Na23_to_n_Na22_derived",  // 960,
        "Na23_to_p_Ne22_derived",  // 961,
        "Na23_to_He4_F19_derived",  // 962,
        "Mg23_to_p_Na22_derived",  // 963,
        "Mg23_to_He4_Ne19_derived",  // 964,
        "Mg24_to_n_Mg23_derived",  // 965,
        "Mg24_to_p_Na23_derived",  // 966,
        "Mg24_to_He4_Ne20_derived",  // 967,
        "Mg25_to_n_Mg24_derived",  // 968,
        "Mg25_to_He4_Ne21_derived",  // 969,
        "Mg26_to_n_Mg25_derived",  // 970,
        "Mg26_to_He4_Ne22_derived",  // 971,
        "Al25_to_p_Mg24_derived",  // 972,
        "Al25_to_He4_Na21_derived",  // 973,
        "Al26_to_n_Al25_derived",  // 974,
        "Al26_to_p_Mg25_derived",  // 975,
        "Al26_to_He4_Na22_derived",  // 976,
        "Al27_to_n_Al26_derived",  // 977,
        "Al27_to_p_Mg26_derived",  // 978,
        "Al27_to_He4_Na23_derived",  // 979,
        "Si28_to_p_Al27_derived",  // 980,
        "Si28_to_He4_Mg24_derived",  // 981,
        "Si29_to_n_Si28_derived",  // 982,
        "Si29_to_He4_Mg25_derived",  // 983,
        "Si30_to_n_Si29_derived",  // 984,
        "Si30_to_He4_Mg26_derived",  // 985,
        "Si31_to_n_Si30_derived",  // 986,
        "Si32_to_n_Si31_derived",  // 987,
        "P29_to_p_Si28_derived",  // 988,
        "P29_to_He4_Al25_derived",  // 989,
        "P30_to_n_P29_derived",  // 990,
        "P30_to_p_Si29_derived",  // 991,
        "P30_to_He4_Al26_derived",  // 992,
        "P31_to_n_P30_derived",  // 993,
        "P31_to_p_Si30_derived",  // 994,
        "P31_to_He4_Al27_derived",  // 995,
        "P32_to_n_P31_derived",  // 996,
        "P32_to_p_Si31_derived",  // 997,
        "P33_to_n_P32_derived",  // 998,
        "P33_to_p_Si32_derived",  // 999,
        "S32_to_p_P31_derived",  // 1000,
        "S32_to_He4_Si28_derived",  // 1001,
        "S33_to_n_S32_derived",  // 1002,
        "S33_to_p_P32_derived",  // 1003,
        "S33_to_He4_Si29_derived",  // 1004,
        "S34_to_n_S33_derived",  // 1005,
        "S34_to_p_P33_derived",  // 1006,
        "S34_to_He4_Si30_derived",  // 1007,
        "S35_to_n_S34_derived",  // 1008,
        "S35_to_He4_Si31_derived",  // 1009,
        "S36_to_n_S35_derived",  // 1010,
        "S36_to_He4_Si32_derived",  // 1011,
        "Cl33_to_p_S32_derived",  // 1012,
        "Cl33_to_He4_P29_derived",  // 1013,
        "Cl34_to_n_Cl33_derived",  // 1014,
        "Cl34_to_p_S33_derived",  // 1015,
        "Cl34_to_He4_P30_derived",  // 1016,
        "Cl35_to_n_Cl34_derived",  // 1017,
        "Cl35_to_p_S34_derived",  // 1018,
        "Cl35_to_He4_P31_derived",  // 1019,
        "Cl36_to_n_Cl35_derived",  // 1020,
        "Cl36_to_p_S35_derived",  // 1021,
        "Cl36_to_He4_P32_derived",  // 1022,
        "Cl37_to_n_Cl36_derived",  // 1023,
        "Cl37_to_p_S36_derived",  // 1024,
        "Cl37_to_He4_P33_derived",  // 1025,
        "Ar36_to_p_Cl35_derived",  // 1026,
        "Ar36_to_He4_S32_derived",  // 1027,
        "Ar37_to_n_Ar36_derived",  // 1028,
        "Ar37_to_p_Cl36_derived",  // 1029,
        "Ar37_to_He4_S33_derived",  // 1030,
        "Ar38_to_n_Ar37_derived",  // 1031,
        "Ar38_to_p_Cl37_derived",  // 1032,
        "Ar38_to_He4_S34_derived",  // 1033,
        "Ar39_to_n_Ar38_derived",  // 1034,
        "Ar39_to_He4_S35_derived",  // 1035,
        "Ar40_to_n_Ar39_derived",  // 1036,
        "Ar40_to_He4_S36_derived",  // 1037,
        "K37_to_p_Ar36_derived",  // 1038,
        "K37_to_He4_Cl33_derived",  // 1039,
        "K38_to_n_K37_derived",  // 1040,
        "K38_to_p_Ar37_derived",  // 1041,
        "K38_to_He4_Cl34_derived",  // 1042,
        "K39_to_n_K38_derived",  // 1043,
        "K39_to_p_Ar38_derived",  // 1044,
        "K39_to_He4_Cl35_derived",  // 1045,
        "K40_to_n_K39_derived",  // 1046,
        "K40_to_p_Ar39_derived",  // 1047,
        "K40_to_He4_Cl36_derived",  // 1048,
        "K41_to_n_K40_derived",  // 1049,
        "K41_to_p_Ar40_derived",  // 1050,
        "K41_to_He4_Cl37_derived",  // 1051,
        "Ca40_to_p_K39_derived",  // 1052,
        "Ca40_to_He4_Ar36_derived",  // 1053,
        "Ca41_to_n_Ca40_derived",  // 1054,
        "Ca41_to_p_K40_derived",  // 1055,
        "Ca41_to_He4_Ar37_derived",  // 1056,
        "Ca42_to_n_Ca41_derived",  // 1057,
        "Ca42_to_p_K41_derived",  // 1058,
        "Ca42_to_He4_Ar38_derived",  // 1059,
        "Ca43_to_n_Ca42_derived",  // 1060,
        "Ca43_to_He4_Ar39_derived",  // 1061,
        "Ca44_to_n_Ca43_derived",  // 1062,
        "Ca44_to_He4_Ar40_derived",  // 1063,
        "Ca45_to_n_Ca44_derived",  // 1064,
        "Ca46_to_n_Ca45_derived",  // 1065,
        "Ca47_to_n_Ca46_derived",  // 1066,
        "Ca48_to_n_Ca47_derived",  // 1067,
        "Sc43_to_p_Ca42_derived",  // 1068,
        "Sc43_to_He4_K39_derived",  // 1069,
        "Sc44_to_n_Sc43_derived",  // 1070,
        "Sc44_to_p_Ca43_derived",  // 1071,
        "Sc44_to_He4_K40_derived",  // 1072,
        "Sc45_to_n_Sc44_derived",  // 1073,
        "Sc45_to_p_Ca44_derived",  // 1074,
        "Sc45_to_He4_K41_derived",  // 1075,
        "Sc46_to_n_Sc45_derived",  // 1076,
        "Sc46_to_p_Ca45_derived",  // 1077,
        "Sc47_to_n_Sc46_derived",  // 1078,
        "Sc47_to_p_Ca46_derived",  // 1079,
        "Sc48_to_n_Sc47_derived",  // 1080,
        "Sc48_to_p_Ca47_derived",  // 1081,
        "Sc49_to_n_Sc48_derived",  // 1082,
        "Sc49_to_p_Ca48_derived",  // 1083,
        "Ti44_to_p_Sc43_derived",  // 1084,
        "Ti44_to_He4_Ca40_derived",  // 1085,
        "Ti45_to_n_Ti44_derived",  // 1086,
        "Ti45_to_p_Sc44_derived",  // 1087,
        "Ti45_to_He4_Ca41_derived",  // 1088,
        "Ti46_to_n_Ti45_derived",  // 1089,
        "Ti46_to_p_Sc45_derived",  // 1090,
        "Ti46_to_He4_Ca42_derived",  // 1091,
        "Ti47_to_n_Ti46_derived",  // 1092,
        "Ti47_to_p_Sc46_derived",  // 1093,
        "Ti47_to_He4_Ca43_derived",  // 1094,
        "Ti48_to_n_Ti47_derived",  // 1095,
        "Ti48_to_p_Sc47_derived",  // 1096,
        "Ti48_to_He4_Ca44_derived",  // 1097,
        "Ti49_to_n_Ti48_derived",  // 1098,
        "Ti49_to_p_Sc48_derived",  // 1099,
        "Ti49_to_He4_Ca45_derived",  // 1100,
        "Ti50_to_n_Ti49_derived",  // 1101,
        "Ti50_to_p_Sc49_derived",  // 1102,
        "Ti50_to_He4_Ca46_derived",  // 1103,
        "Ti51_to_n_Ti50_derived",  // 1104,
        "Ti51_to_He4_Ca47_derived",  // 1105,
        "V46_to_p_Ti45_derived",  // 1106,
        "V47_to_n_V46_derived",  // 1107,
        "V47_to_p_Ti46_derived",  // 1108,
        "V47_to_He4_Sc43_derived",  // 1109,
        "V48_to_n_V47_derived",  // 1110,
        "V48_to_p_Ti47_derived",  // 1111,
        "V48_to_He4_Sc44_derived",  // 1112,
        "V49_to_n_V48_derived",  // 1113,
        "V49_to_p_Ti48_derived",  // 1114,
        "V49_to_He4_Sc45_derived",  // 1115,
        "V50_to_n_V49_derived",  // 1116,
        "V50_to_p_Ti49_derived",  // 1117,
        "V50_to_He4_Sc46_derived",  // 1118,
        "V51_to_n_V50_derived",  // 1119,
        "V51_to_p_Ti50_derived",  // 1120,
        "V51_to_He4_Sc47_derived",  // 1121,
        "V52_to_n_V51_derived",  // 1122,
        "V52_to_p_Ti51_derived",  // 1123,
        "V52_to_He4_Sc48_derived",  // 1124,
        "Cr48_to_p_V47_derived",  // 1125,
        "Cr48_to_He4_Ti44_derived",  // 1126,
        "Cr49_to_n_Cr48_derived",  // 1127,
        "Cr49_to_p_V48_derived",  // 1128,
        "Cr49_to_He4_Ti45_derived",  // 1129,
        "Cr50_to_n_Cr49_derived",  // 1130,
        "Cr50_to_p_V49_derived",  // 1131,
        "Cr50_to_He4_Ti46_derived",  // 1132,
        "Cr51_to_n_Cr50_derived",  // 1133,
        "Cr51_to_p_V50_derived",  // 1134,
        "Cr51_to_He4_Ti47_derived",  // 1135,
        "Cr52_to_n_Cr51_derived",  // 1136,
        "Cr52_to_p_V51_derived",  // 1137,
        "Cr52_to_He4_Ti48_derived",  // 1138,
        "Cr53_to_n_Cr52_derived",  // 1139,
        "Cr53_to_p_V52_derived",  // 1140,
        "Cr53_to_He4_Ti49_derived",  // 1141,
        "Cr54_to_n_Cr53_derived",  // 1142,
        "Cr54_to_He4_Ti50_derived",  // 1143,
        "Mn50_to_p_Cr49_derived",  // 1144,
        "Mn50_to_He4_V46_derived",  // 1145,
        "Mn51_to_n_Mn50_derived",  // 1146,
        "Mn51_to_p_Cr50_derived",  // 1147,
        "Mn51_to_He4_V47_derived",  // 1148,
        "Mn52_to_n_Mn51_derived",  // 1149,
        "Mn52_to_p_Cr51_derived",  // 1150,
        "Mn52_to_He4_V48_derived",  // 1151,
        "Mn53_to_n_Mn52_derived",  // 1152,
        "Mn53_to_p_Cr52_derived",  // 1153,
        "Mn53_to_He4_V49_derived",  // 1154,
        "Mn54_to_n_Mn53_derived",  // 1155,
        "Mn54_to_p_Cr53_derived",  // 1156,
        "Mn54_to_He4_V50_derived",  // 1157,
        "Mn55_to_n_Mn54_derived",  // 1158,
        "Mn55_to_p_Cr54_derived",  // 1159,
        "Mn55_to_He4_V51_derived",  // 1160,
        "Fe52_to_p_Mn51_derived",  // 1161,
        "Fe52_to_He4_Cr48_derived",  // 1162,
        "Fe53_to_n_Fe52_derived",  // 1163,
        "Fe53_to_p_Mn52_derived",  // 1164,
        "Fe53_to_He4_Cr49_derived",  // 1165,
        "Fe54_to_n_Fe53_derived",  // 1166,
        "Fe54_to_p_Mn53_derived",  // 1167,
        "Fe54_to_He4_Cr50_derived",  // 1168,
        "Fe55_to_n_Fe54_derived",  // 1169,
        "Fe55_to_p_Mn54_derived",  // 1170,
        "Fe55_to_He4_Cr51_derived",  // 1171,
        "Fe56_to_n_Fe55_derived",  // 1172,
        "Fe56_to_p_Mn55_derived",  // 1173,
        "Fe56_to_He4_Cr52_derived",  // 1174,
        "Fe57_to_n_Fe56_derived",  // 1175,
        "Fe57_to_He4_Cr53_derived",  // 1176,
        "Fe58_to_n_Fe57_derived",  // 1177,
        "Fe58_to_He4_Cr54_derived",  // 1178,
        "Co53_to_p_Fe52_derived",  // 1179,
        "Co54_to_n_Co53_derived",  // 1180,
        "Co54_to_p_Fe53_derived",  // 1181,
        "Co54_to_He4_Mn50_derived",  // 1182,
        "Co55_to_n_Co54_derived",  // 1183,
        "Co55_to_p_Fe54_derived",  // 1184,
        "Co55_to_He4_Mn51_derived",  // 1185,
        "Co56_to_n_Co55_derived",  // 1186,
        "Co56_to_p_Fe55_derived",  // 1187,
        "Co56_to_He4_Mn52_derived",  // 1188,
        "Co57_to_n_Co56_derived",  // 1189,
        "Co57_to_p_Fe56_derived",  // 1190,
        "Co57_to_He4_Mn53_derived",  // 1191,
        "Co58_to_n_Co57_derived",  // 1192,
        "Co58_to_p_Fe57_derived",  // 1193,
        "Co58_to_He4_Mn54_derived",  // 1194,
        "Co59_to_n_Co58_derived",  // 1195,
        "Co59_to_p_Fe58_derived",  // 1196,
        "Co59_to_He4_Mn55_derived",  // 1197,
        "Ni56_to_p_Co55_derived",  // 1198,
        "Ni56_to_He4_Fe52_derived",  // 1199,
        "Ni57_to_n_Ni56_derived",  // 1200,
        "Ni57_to_p_Co56_derived",  // 1201,
        "Ni57_to_He4_Fe53_derived",  // 1202,
        "Ni58_to_n_Ni57_derived",  // 1203,
        "Ni58_to_p_Co57_derived",  // 1204,
        "Ni58_to_He4_Fe54_derived",  // 1205,
        "Ni59_to_n_Ni58_derived",  // 1206,
        "Ni59_to_p_Co58_derived",  // 1207,
        "Ni59_to_He4_Fe55_derived",  // 1208,
        "Ni60_to_n_Ni59_derived",  // 1209,
        "Ni60_to_p_Co59_derived",  // 1210,
        "Ni60_to_He4_Fe56_derived",  // 1211,
        "Ni61_to_n_Ni60_derived",  // 1212,
        "Ni61_to_He4_Fe57_derived",  // 1213,
        "Ni62_to_n_Ni61_derived",  // 1214,
        "Ni62_to_He4_Fe58_derived",  // 1215,
        "Ni63_to_n_Ni62_derived",  // 1216,
        "Ni64_to_n_Ni63_derived",  // 1217,
        "Cu57_to_p_Ni56_derived",  // 1218,
        "Cu57_to_He4_Co53_derived",  // 1219,
        "Cu58_to_n_Cu57_derived",  // 1220,
        "Cu58_to_p_Ni57_derived",  // 1221,
        "Cu58_to_He4_Co54_derived",  // 1222,
        "Cu59_to_n_Cu58_derived",  // 1223,
        "Cu59_to_p_Ni58_derived",  // 1224,
        "Cu59_to_He4_Co55_derived",  // 1225,
        "Cu60_to_n_Cu59_derived",  // 1226,
        "Cu60_to_p_Ni59_derived",  // 1227,
        "Cu60_to_He4_Co56_derived",  // 1228,
        "Cu61_to_n_Cu60_derived",  // 1229,
        "Cu61_to_p_Ni60_derived",  // 1230,
        "Cu61_to_He4_Co57_derived",  // 1231,
        "Cu62_to_n_Cu61_derived",  // 1232,
        "Cu62_to_p_Ni61_derived",  // 1233,
        "Cu62_to_He4_Co58_derived",  // 1234,
        "Cu63_to_n_Cu62_derived",  // 1235,
        "Cu63_to_p_Ni62_derived",  // 1236,
        "Cu63_to_He4_Co59_derived",  // 1237,
        "Cu64_to_n_Cu63_derived",  // 1238,
        "Cu64_to_p_Ni63_derived",  // 1239,
        "Cu65_to_n_Cu64_derived",  // 1240,
        "Cu65_to_p_Ni64_derived",  // 1241,
        "Zn59_to_p_Cu58_derived",  // 1242,
        "Zn60_to_n_Zn59_derived",  // 1243,
        "Zn60_to_p_Cu59_derived",  // 1244,
        "Zn60_to_He4_Ni56_derived",  // 1245,
        "Zn61_to_n_Zn60_derived",  // 1246,
        "Zn61_to_p_Cu60_derived",  // 1247,
        "Zn61_to_He4_Ni57_derived",  // 1248,
        "Zn62_to_n_Zn61_derived",  // 1249,
        "Zn62_to_p_Cu61_derived",  // 1250,
        "Zn62_to_He4_Ni58_derived",  // 1251,
        "Zn63_to_n_Zn62_derived",  // 1252,
        "Zn63_to_p_Cu62_derived",  // 1253,
        "Zn63_to_He4_Ni59_derived",  // 1254,
        "Zn64_to_n_Zn63_derived",  // 1255,
        "Zn64_to_p_Cu63_derived",  // 1256,
        "Zn64_to_He4_Ni60_derived",  // 1257,
        "Zn65_to_n_Zn64_derived",  // 1258,
        "Zn65_to_p_Cu64_derived",  // 1259,
        "Zn65_to_He4_Ni61_derived",  // 1260,
        "Zn66_to_n_Zn65_derived",  // 1261,
        "Zn66_to_p_Cu65_derived",  // 1262,
        "Zn66_to_He4_Ni62_derived",  // 1263,
        "Ga62_to_p_Zn61_derived",  // 1264,
        "Ga62_to_He4_Cu58_derived",  // 1265,
        "Ga63_to_n_Ga62_derived",  // 1266,
        "Ga63_to_p_Zn62_derived",  // 1267,
        "Ga63_to_He4_Cu59_derived",  // 1268,
        "Ga64_to_n_Ga63_derived",  // 1269,
        "Ga64_to_p_Zn63_derived",  // 1270,
        "Ga64_to_He4_Cu60_derived",  // 1271,
        "Ge63_to_p_Ga62_derived",  // 1272,
        "Ge63_to_He4_Zn59_derived",  // 1273,
        "Ge64_to_n_Ge63_derived",  // 1274,
        "Ge64_to_p_Ga63_derived",  // 1275,
        "Ge64_to_He4_Zn60_derived",  // 1276,
        "Li6_to_n_p_He4_derived",  // 1277,
        "Be9_to_n_He4_He4_derived",  // 1278,
        "C12_to_3He4_derived",  // 1279,
        "He3_n_to_d_d_derived",  // 1280,
        "He3_He4_to_p_Li6_derived",  // 1281,
        "He4_p_to_d_He3_derived",  // 1282,
        "He4_He4_to_n_Be7_derived",  // 1283,
        "He4_He4_to_p_Li7_derived",  // 1284,
        "Li6_He4_to_p_Be9_derived",  // 1285,
        "Li7_p_to_n_Be7_derived",  // 1286,
        "Li7_p_to_d_Li6_derived",  // 1287,
        "Li7_He4_to_n_B10_derived",  // 1288,
        "Be7_n_to_d_Li6_derived",  // 1289,
        "Be7_He4_to_p_B10_derived",  // 1290,
        "C12_n_to_He4_Be9_derived",  // 1291,
        "C12_He4_to_n_O15_derived",  // 1292,
        "C12_He4_to_p_N15_derived",  // 1293,
        "C12_C12_to_n_Mg23_derived",  // 1294,
        "C13_p_to_He4_B10_derived",  // 1295,
        "C14_p_to_n_N14_derived",  // 1296,
        "C14_p_to_He4_B11_derived",  // 1297,
        "C14_He4_to_n_O17_derived",  // 1298,
        "N13_n_to_p_C13_derived",  // 1299,
        "N13_n_to_He4_B10_derived",  // 1300,
        "N14_n_to_d_C13_derived",  // 1301,
        "N14_n_to_He4_B11_derived",  // 1302,
        "N14_He4_to_p_O17_derived",  // 1303,
        "N15_n_to_d_C14_derived",  // 1304,
        "N15_He4_to_n_F18_derived",  // 1305,
        "N15_He4_to_p_O18_derived",  // 1306,
        "O14_n_to_p_N14_derived",  // 1307,
        "O15_n_to_p_N15_derived",  // 1308,
        "O15_He4_to_n_Ne18_derived",  // 1309,
        "O15_He4_to_p_F18_derived",  // 1310,
        "O16_n_to_He4_C13_derived",  // 1311,
        "O16_p_to_He4_N13_derived",  // 1312,
        "O16_He4_to_n_Ne19_derived",  // 1313,
        "O16_He4_to_p_F19_derived",  // 1314,
        "O17_p_to_n_F17_derived",  // 1315,
        "O18_p_to_n_F18_derived",  // 1316,
        "F17_n_to_He4_N14_derived",  // 1317,
        "F17_p_to_He4_O14_derived",  // 1318,
        "F17_He4_to_p_Ne20_derived",  // 1319,
        "F18_p_to_n_Ne18_derived",  // 1320,
        "F18_He4_to_n_Na21_derived",  // 1321,
        "F19_He4_to_n_Na22_derived",  // 1322,
        "Ne19_n_to_p_F19_derived",  // 1323,
        "Ne20_n_to_He4_O17_derived",  // 1324,
        "Ne20_He4_to_n_Mg23_derived",  // 1325,
        "Ne20_He4_to_p_Na23_derived",  // 1326,
        "Ne20_He4_to_C12_C12_derived",  // 1327,
        "Ne21_n_to_He4_O18_derived",  // 1328,
        "Ne21_p_to_n_Na21_derived",  // 1329,
        "Ne21_p_to_He4_F18_derived",  // 1330,
        "Ne22_p_to_n_Na22_derived",  // 1331,
        "Ne22_p_to_He4_F19_derived",  // 1332,
        "Na21_p_to_He4_Ne18_derived",  // 1333,
        "Na21_He4_to_p_Mg24_derived",  // 1334,
        "Na22_p_to_He4_Ne19_derived",  // 1335,
        "Na22_He4_to_n_Al25_derived",  // 1336,
        "Na23_p_to_C12_C12_derived",  // 1337,
        "Na23_He4_to_n_Al26_derived",  // 1338,
        "Mg23_n_to_p_Na23_derived",  // 1339,
        "Mg24_n_to_He4_Ne21_derived",  // 1340,
        "Mg24_He4_to_p_Al27_derived",  // 1341,
        "Mg24_He4_to_C12_O16_derived",  // 1342,
        "Mg25_n_to_He4_Ne22_derived",  // 1343,
        "Mg25_p_to_n_Al25_derived",  // 1344,
        "Mg25_p_to_He4_Na22_derived",  // 1345,
        "Mg26_p_to_n_Al26_derived",  // 1346,
        "Mg26_p_to_He4_Na23_derived",  // 1347,
        "Al26_p_to_He4_Mg23_derived",  // 1348,
        "Al26_He4_to_n_P29_derived",  // 1349,
        "Al27_p_to_C12_O16_derived",  // 1350,
        "Si28_n_to_He4_Mg25_derived",  // 1351,
        "Si28_p_to_He4_Al25_derived",  // 1352,
        "Si28_He4_to_p_P31_derived",  // 1353,
        "Si28_He4_to_C12_Ne20_derived",  // 1354,
        "Si28_He4_to_O16_O16_derived",  // 1355,
        "Si29_n_to_He4_Mg26_derived",  // 1356,
        "Si29_p_to_n_P29_derived",  // 1357,
        "Si29_p_to_He4_Al26_derived",  // 1358,
        "Si29_He4_to_n_S32_derived",  // 1359,
        "Si29_He4_to_p_P32_derived",  // 1360,
        "Si30_p_to_n_P30_derived",  // 1361,
        "Si30_p_to_He4_Al27_derived",  // 1362,
        "Si30_He4_to_n_S33_derived",  // 1363,
        "Si30_He4_to_p_P33_derived",  // 1364,
        "Si32_p_to_n_P32_derived",  // 1365,
        "Si32_He4_to_n_S35_derived",  // 1366,
        "P30_n_to_He4_Al27_derived",  // 1367,
        "P30_He4_to_n_Cl33_derived",  // 1368,
        "P31_n_to_p_Si31_derived",  // 1369,
        "P31_p_to_C12_Ne20_derived",  // 1370,
        "P31_p_to_O16_O16_derived",  // 1371,
        "P31_He4_to_n_Cl34_derived",  // 1372,
        "P32_He4_to_n_Cl35_derived",  // 1373,
        "P32_He4_to_p_S35_derived",  // 1374,
        "P33_p_to_n_S33_derived",  // 1375,
        "P33_He4_to_n_Cl36_derived",  // 1376,
        "P33_He4_to_p_S36_derived",  // 1377,
        "S32_n_to_p_P32_derived",  // 1378,
        "S32_p_to_He4_P29_derived",  // 1379,
        "S32_He4_to_p_Cl35_derived",  // 1380,
        "S33_p_to_n_Cl33_derived",  // 1381,
        "S33_p_to_He4_P30_derived",  // 1382,
        "S33_He4_to_n_Ar36_derived",  // 1383,
        "S33_He4_to_p_Cl36_derived",  // 1384,
        "S34_n_to_He4_Si31_derived",  // 1385,
        "S34_p_to_n_Cl34_derived",  // 1386,
        "S34_p_to_He4_P31_derived",  // 1387,
        "S35_p_to_n_Cl35_derived",  // 1388,
        "S36_p_to_n_Cl36_derived",  // 1389,
        "S36_He4_to_n_Ar39_derived",  // 1390,
        "Cl34_He4_to_n_K37_derived",  // 1391,
        "Cl35_He4_to_n_K38_derived",  // 1392,
        "Cl36_p_to_n_Ar36_derived",  // 1393,
        "Cl36_He4_to_n_K39_derived",  // 1394,
        "Cl36_He4_to_p_Ar39_derived",  // 1395,
        "Cl37_p_to_n_Ar37_derived",  // 1396,
        "Cl37_p_to_He4_S34_derived",  // 1397,
        "Cl37_He4_to_n_K40_derived",  // 1398,
        "Cl37_He4_to_p_Ar40_derived",  // 1399,
        "Ar36_p_to_He4_Cl33_derived",  // 1400,
        "Ar36_He4_to_p_K39_derived",  // 1401,
        "Ar37_n_to_He4_S34_derived",  // 1402,
        "Ar37_p_to_n_K37_derived",  // 1403,
        "Ar37_p_to_He4_Cl34_derived",  // 1404,
        "Ar37_He4_to_n_Ca40_derived",  // 1405,
        "Ar37_He4_to_p_K40_derived",  // 1406,
        "Ar38_n_to_He4_S35_derived",  // 1407,
        "Ar38_p_to_n_K38_derived",  // 1408,
        "Ar38_p_to_He4_Cl35_derived",  // 1409,
        "Ar38_He4_to_n_Ca41_derived",  // 1410,
        "Ar38_He4_to_p_K41_derived",  // 1411,
        "Ar39_p_to_n_K39_derived",  // 1412,
        "Ar39_He4_to_n_Ca42_derived",  // 1413,
        "Ar40_p_to_n_K40_derived",  // 1414,
        "Ar40_He4_to_n_Ca43_derived",  // 1415,
        "K39_He4_to_p_Ca42_derived",  // 1416,
        "K40_He4_to_n_Sc43_derived",  // 1417,
        "Ca40_n_to_p_K40_derived",  // 1418,
        "Ca40_p_to_He4_K37_derived",  // 1419,
        "Ca40_He4_to_p_Sc43_derived",  // 1420,
        "Ca41_n_to_p_K41_derived",  // 1421,
        "Ca41_p_to_He4_K38_derived",  // 1422,
        "Ca41_He4_to_n_Ti44_derived",  // 1423,
        "Ca41_He4_to_p_Sc44_derived",  // 1424,
        "Ca42_He4_to_n_Ti45_derived",  // 1425,
        "Ca42_He4_to_p_Sc45_derived",  // 1426,
        "Ca43_p_to_n_Sc43_derived",  // 1427,
        "Ca43_p_to_He4_K40_derived",  // 1428,
        "Ca43_He4_to_p_Sc46_derived",  // 1429,
        "Ca44_p_to_n_Sc44_derived",  // 1430,
        "Ca44_p_to_He4_K41_derived",  // 1431,
        "Ca44_He4_to_n_Ti47_derived",  // 1432,
        "Ca44_He4_to_p_Sc47_derived",  // 1433,
        "Ca45_p_to_n_Sc45_derived",  // 1434,
        "Ca45_He4_to_p_Sc48_derived",  // 1435,
        "Ca46_p_to_n_Sc46_derived",  // 1436,
        "Ca46_He4_to_n_Ti49_derived",  // 1437,
        "Ca46_He4_to_p_Sc49_derived",  // 1438,
        "Ca48_p_to_n_Sc48_derived",  // 1439,
        "Ca48_He4_to_n_Ti51_derived",  // 1440,
        "Sc43_He4_to_n_V46_derived",  // 1441,
        "Sc44_n_to_He4_K41_derived",  // 1442,
        "Sc44_p_to_n_Ti44_derived",  // 1443,
        "Sc44_He4_to_n_V47_derived",  // 1444,
        "Sc45_p_to_n_Ti45_derived",  // 1445,
        "Sc45_He4_to_n_V48_derived",  // 1446,
        "Sc47_n_to_p_Ca47_derived",  // 1447,
        "Sc47_p_to_n_Ti47_derived",  // 1448,
        "Sc47_He4_to_n_V50_derived",  // 1449,
        "Sc49_He4_to_n_V52_derived",  // 1450,
        "Ti45_He4_to_n_Cr48_derived",  // 1451,
        "Ti46_n_to_p_Sc46_derived",  // 1452,
        "Ti46_n_to_He4_Ca43_derived",  // 1453,
        "Ti46_p_to_n_V46_derived",  // 1454,
        "Ti46_p_to_He4_Sc43_derived",  // 1455,
        "Ti46_He4_to_n_Cr49_derived",  // 1456,
        "Ti46_He4_to_p_V49_derived",  // 1457,
        "Ti47_p_to_n_V47_derived",  // 1458,
        "Ti47_p_to_He4_Sc44_derived",  // 1459,
        "Ti47_He4_to_n_Cr50_derived",  // 1460,
        "Ti47_He4_to_p_V50_derived",  // 1461,
        "Ti48_n_to_p_Sc48_derived",  // 1462,
        "Ti48_n_to_He4_Ca45_derived",  // 1463,
        "Ti48_p_to_n_V48_derived",  // 1464,
        "Ti48_p_to_He4_Sc45_derived",  // 1465,
        "Ti48_He4_to_n_Cr51_derived",  // 1466,
        "Ti48_He4_to_p_V51_derived",  // 1467,
        "Ti49_n_to_p_Sc49_derived",  // 1468,
        "Ti49_p_to_n_V49_derived",  // 1469,
        "Ti49_p_to_He4_Sc46_derived",  // 1470,
        "Ti49_He4_to_p_V52_derived",  // 1471,
        "Ti50_n_to_He4_Ca47_derived",  // 1472,
        "Ti50_p_to_n_V50_derived",  // 1473,
        "Ti50_p_to_He4_Sc47_derived",  // 1474,
        "Ti50_He4_to_n_Cr53_derived",  // 1475,
        "Ti51_p_to_He4_Sc48_derived",  // 1476,
        "V47_p_to_He4_Ti44_derived",  // 1477,
        "V47_He4_to_n_Mn50_derived",  // 1478,
        "V48_p_to_n_Cr48_derived",  // 1479,
        "V48_p_to_He4_Ti45_derived",  // 1480,
        "V48_He4_to_n_Mn51_derived",  // 1481,
        "V49_n_to_He4_Sc46_derived",  // 1482,
        "V49_p_to_n_Cr49_derived",  // 1483,
        "V49_He4_to_n_Mn52_derived",  // 1484,
        "V50_He4_to_n_Mn53_derived",  // 1485,
        "V51_n_to_p_Ti51_derived",  // 1486,
        "V51_n_to_He4_Sc48_derived",  // 1487,
        "V51_p_to_n_Cr51_derived",  // 1488,
        "V51_He4_to_n_Mn54_derived",  // 1489,
        "V51_He4_to_p_Cr54_derived",  // 1490,
        "Cr49_p_to_He4_V46_derived",  // 1491,
        "Cr49_He4_to_n_Fe52_derived",  // 1492,
        "Cr50_n_to_p_V50_derived",  // 1493,
        "Cr50_p_to_n_Mn50_derived",  // 1494,
        "Cr50_p_to_He4_V47_derived",  // 1495,
        "Cr50_He4_to_n_Fe53_derived",  // 1496,
        "Cr50_He4_to_p_Mn53_derived",  // 1497,
        "Cr51_p_to_n_Mn51_derived",  // 1498,
        "Cr51_p_to_He4_V48_derived",  // 1499,
        "Cr51_He4_to_n_Fe54_derived",  // 1500,
        "Cr51_He4_to_p_Mn54_derived",  // 1501,
        "Cr52_n_to_p_V52_derived",  // 1502,
        "Cr52_n_to_He4_Ti49_derived",  // 1503,
        "Cr52_p_to_n_Mn52_derived",  // 1504,
        "Cr52_p_to_He4_V49_derived",  // 1505,
        "Cr52_He4_to_n_Fe55_derived",  // 1506,
        "Cr52_He4_to_p_Mn55_derived",  // 1507,
        "Cr53_p_to_n_Mn53_derived",  // 1508,
        "Cr53_p_to_He4_V50_derived",  // 1509,
        "Cr53_He4_to_n_Fe56_derived",  // 1510,
        "Cr54_n_to_He4_Ti51_derived",  // 1511,
        "Cr54_p_to_n_Mn54_derived",  // 1512,
        "Cr54_He4_to_n_Fe57_derived",  // 1513,
        "Mn50_He4_to_n_Co53_derived",  // 1514,
        "Mn51_p_to_He4_Cr48_derived",  // 1515,
        "Mn51_He4_to_n_Co54_derived",  // 1516,
        "Mn52_p_to_n_Fe52_derived",  // 1517,
        "Mn52_p_to_He4_Cr49_derived",  // 1518,
        "Mn52_He4_to_n_Co55_derived",  // 1519,
        "Mn53_p_to_n_Fe53_derived",  // 1520,
        "Mn53_He4_to_n_Co56_derived",  // 1521,
        "Mn54_p_to_n_Fe54_derived",  // 1522,
        "Mn54_He4_to_n_Co57_derived",  // 1523,
        "Mn54_He4_to_p_Fe57_derived",  // 1524,
        "Mn55_n_to_He4_V52_derived",  // 1525,
        "Mn55_p_to_n_Fe55_derived",  // 1526,
        "Mn55_He4_to_n_Co58_derived",  // 1527,
        "Mn55_He4_to_p_Fe58_derived",  // 1528,
        "Fe53_p_to_n_Co53_derived",  // 1529,
        "Fe53_p_to_He4_Mn50_derived",  // 1530,
        "Fe53_He4_to_n_Ni56_derived",  // 1531,
        "Fe54_p_to_n_Co54_derived",  // 1532,
        "Fe54_p_to_He4_Mn51_derived",  // 1533,
        "Fe54_He4_to_n_Ni57_derived",  // 1534,
        "Fe54_He4_to_p_Co57_derived",  // 1535,
        "Fe55_p_to_n_Co55_derived",  // 1536,
        "Fe55_p_to_He4_Mn52_derived",  // 1537,
        "Fe55_He4_to_n_Ni58_derived",  // 1538,
        "Fe55_He4_to_p_Co58_derived",  // 1539,
        "Fe56_p_to_n_Co56_derived",  // 1540,
        "Fe56_p_to_He4_Mn53_derived",  // 1541,
        "Fe56_He4_to_n_Ni59_derived",  // 1542,
        "Fe56_He4_to_p_Co59_derived",  // 1543,
        "Fe57_p_to_n_Co57_derived",  // 1544,
        "Fe57_He4_to_n_Ni60_derived",  // 1545,
        "Fe58_p_to_n_Co58_derived",  // 1546,
        "Fe58_He4_to_n_Ni61_derived",  // 1547,
        "Co54_He4_to_n_Cu57_derived",  // 1548,
        "Co55_p_to_He4_Fe52_derived",  // 1549,
        "Co55_He4_to_n_Cu58_derived",  // 1550,
        "Co56_p_to_n_Ni56_derived",  // 1551,
        "Co56_p_to_He4_Fe53_derived",  // 1552,
        "Co56_He4_to_n_Cu59_derived",  // 1553,
        "Co57_p_to_n_Ni57_derived",  // 1554,
        "Co57_He4_to_n_Cu60_derived",  // 1555,
        "Co58_p_to_n_Ni58_derived",  // 1556,
        "Co58_He4_to_n_Cu61_derived",  // 1557,
        "Co58_He4_to_p_Ni61_derived",  // 1558,
        "Co59_p_to_n_Ni59_derived",  // 1559,
        "Co59_He4_to_n_Cu62_derived",  // 1560,
        "Co59_He4_to_p_Ni62_derived",  // 1561,
        "Ni56_p_to_He4_Co53_derived",  // 1562,
        "Ni56_He4_to_n_Zn59_derived",  // 1563,
        "Ni56_He4_to_p_Cu59_derived",  // 1564,
        "Ni57_p_to_n_Cu57_derived",  // 1565,
        "Ni57_p_to_He4_Co54_derived",  // 1566,
        "Ni57_He4_to_n_Zn60_derived",  // 1567,
        "Ni57_He4_to_p_Cu60_derived",  // 1568,
        "Ni58_p_to_n_Cu58_derived",  // 1569,
        "Ni58_p_to_He4_Co55_derived",  // 1570,
        "Ni58_He4_to_n_Zn61_derived",  // 1571,
        "Ni58_He4_to_p_Cu61_derived",  // 1572,
        "Ni59_p_to_n_Cu59_derived",  // 1573,
        "Ni59_p_to_He4_Co56_derived",  // 1574,
        "Ni59_He4_to_n_Zn62_derived",  // 1575,
        "Ni59_He4_to_p_Cu62_derived",  // 1576,
        "Ni60_p_to_n_Cu60_derived",  // 1577,
        "Ni60_p_to_He4_Co57_derived",  // 1578,
        "Ni60_He4_to_n_Zn63_derived",  // 1579,
        "Ni60_He4_to_p_Cu63_derived",  // 1580,
        "Ni61_p_to_n_Cu61_derived",  // 1581,
        "Ni61_He4_to_n_Zn64_derived",  // 1582,
        "Ni61_He4_to_p_Cu64_derived",  // 1583,
        "Ni62_p_to_n_Cu62_derived",  // 1584,
        "Ni62_He4_to_n_Zn65_derived",  // 1585,
        "Ni62_He4_to_p_Cu65_derived",  // 1586,
        "Ni63_p_to_n_Cu63_derived",  // 1587,
        "Ni63_He4_to_n_Zn66_derived",  // 1588,
        "Ni64_p_to_n_Cu64_derived",  // 1589,
        "Cu58_He4_to_p_Zn61_derived",  // 1590,
        "Cu59_p_to_n_Zn59_derived",  // 1591,
        "Cu59_He4_to_n_Ga62_derived",  // 1592,
        "Cu59_He4_to_p_Zn62_derived",  // 1593,
        "Cu60_p_to_n_Zn60_derived",  // 1594,
        "Cu60_He4_to_n_Ga63_derived",  // 1595,
        "Cu60_He4_to_p_Zn63_derived",  // 1596,
        "Cu61_p_to_n_Zn61_derived",  // 1597,
        "Cu61_He4_to_n_Ga64_derived",  // 1598,
        "Cu61_He4_to_p_Zn64_derived",  // 1599,
        "Cu62_p_to_n_Zn62_derived",  // 1600,
        "Cu62_He4_to_p_Zn65_derived",  // 1601,
        "Cu63_p_to_n_Zn63_derived",  // 1602,
        "Cu63_He4_to_p_Zn66_derived",  // 1603,
        "Cu64_p_to_n_Zn64_derived",  // 1604,
        "Cu65_p_to_n_Zn65_derived",  // 1605,
        "Zn59_He4_to_p_Ga62_derived",  // 1606,
        "Zn60_p_to_He4_Cu57_derived",  // 1607,
        "Zn60_He4_to_n_Ge63_derived",  // 1608,
        "Zn60_He4_to_p_Ga63_derived",  // 1609,
        "Zn61_He4_to_n_Ge64_derived",  // 1610,
        "Zn61_He4_to_p_Ga64_derived",  // 1611,
        "Zn62_p_to_n_Ga62_derived",  // 1612,
        "Zn63_p_to_n_Ga63_derived",  // 1613,
        "Zn64_p_to_n_Ga64_derived",  // 1614,
        "Ga63_p_to_n_Ge63_derived",  // 1615,
        "Ga64_p_to_n_Ge64_derived",  // 1616,
        "d_p_to_n_p_p_derived",  // 1617,
        "He4_p_p_to_He3_He3_derived",  // 1618,
        "He4_He4_n_to_d_Li7_derived",  // 1619,
        "He4_He4_p_to_n_B8_derived",  // 1620,
        "He4_He4_p_to_d_Be7_derived",  // 1621,
        "d_He4_He4_to_p_Be9_derived",  // 1622,
        "He4_He4_He4_to_p_B11_derived",  // 1623,
        "He4_He4_p_n_to_He3_Li7_derived",  // 1624,
        "He4_He4_p_n_to_p_Be9_derived",  // 1625,
        "He4_He4_p_p_to_He3_Be7_derived"  // 1626,
    };

    template<std::uint8_t rate>
    AMREX_GPU_DEVICE AMREX_INLINE
    constexpr amrex::Real get_p_random() {
        return 0.0_rt;
    }

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    // A RatePair data table, 2D array of shape (NumNSERatePairs, 10).
    // Each row corresponds to a RatePair returned by _get_nse_rate_pairs().
    // Use enum names following NetworkSpecies and NetworkRates.

    //  - Columns 1–3: 1-based index of reactant nuclei (forward rate)
    //  - Columns 4–6: 1-based index of product nuclei (forward rate)
    //  - Columns 7-8: 1-based index of non (n, H1, He4) nuclei.
    //  - Column  9:   1-based index of the forward rate
    //  - Column  10:  1-based index of the reverse rate

    using namespace Species;
    using namespace Rates;

    constexpr int H1_index = 2;
    constexpr int N_index = 1;
    constexpr int He4_index = 5;
    constexpr int NumNSERatePairs = 695;

    inline AMREX_GPU_MANAGED amrex::Array2D<short, 1, NumNSERatePairs, 1, 10, amrex::Order::C> rate_pair_data {
        N, H1, -1, H2, -1, -1, H2, -1, k_p_n_to_d_reaclib, k_d_to_n_p_derived,
        H1, H2, -1, He3, -1, -1, H2, He3, k_d_p_to_He3_reaclib, k_He3_to_p_d_derived,
        H2, H2, -1, He4, -1, -1, H2, H2, k_d_d_to_He4_reaclib, k_He4_to_d_d_derived,
        He4, H2, -1, Li6, -1, -1, H2, Li6, k_d_He4_to_Li6_reaclib, k_Li6_to_He4_d_derived,
        N, He3, -1, He4, -1, -1, He3, -1, k_He3_n_to_He4_reaclib, k_He4_to_n_He3_derived,
        He4, He3, -1, Be7, -1, -1, He3, Be7, k_He3_He4_to_Be7_reaclib, k_Be7_to_He4_He3_derived,
        N, Li6, -1, Li7, -1, -1, Li6, Li7, k_Li6_n_to_Li7_reaclib, k_Li7_to_n_Li6_derived,
        H1, Li6, -1, Be7, -1, -1, Li6, Be7, k_Li6_p_to_Be7_reaclib, k_Be7_to_p_Li6_derived,
        He4, Li6, -1, B10, -1, -1, Li6, B10, k_Li6_He4_to_B10_reaclib, k_B10_to_He4_Li6_derived,
        He4, Li7, -1, B11, -1, -1, Li7, B11, k_Li7_He4_to_B11_reaclib, k_B11_to_He4_Li7_derived,
        H1, Be7, -1, B8, -1, -1, Be7, B8, k_Be7_p_to_B8_reaclib, k_B8_to_p_Be7_derived,
        H1, Be9, -1, B10, -1, -1, Be9, B10, k_Be9_p_to_B10_reaclib, k_B10_to_p_Be9_derived,
        N, B10, -1, B11, -1, -1, B10, B11, k_B10_n_to_B11_reaclib, k_B11_to_n_B10_derived,
        H1, B11, -1, C12, -1, -1, B11, C12, k_B11_p_to_C12_reaclib, k_C12_to_p_B11_derived,
        N, C12, -1, C13, -1, -1, C12, C13, k_C12_n_to_C13_reaclib, k_C13_to_n_C12_derived,
        H1, C12, -1, N13, -1, -1, C12, N13, k_C12_p_to_N13_reaclib, k_N13_to_p_C12_derived,
        He4, C12, -1, O16, -1, -1, C12, O16, k_C12_He4_to_O16_reaclib, k_O16_to_He4_C12_derived,
        N, C13, -1, C14, -1, -1, C13, C14, k_C13_n_to_C14_reaclib, k_C14_to_n_C13_derived,
        H1, C13, -1, N14, -1, -1, C13, N14, k_C13_p_to_N14_reaclib, k_N14_to_p_C13_derived,
        H1, C14, -1, N15, -1, -1, C14, N15, k_C14_p_to_N15_reaclib, k_N15_to_p_C14_derived,
        He4, C14, -1, O18, -1, -1, C14, O18, k_C14_He4_to_O18_reaclib, k_O18_to_He4_C14_derived,
        N, N13, -1, N14, -1, -1, N13, N14, k_N13_n_to_N14_reaclib, k_N14_to_n_N13_derived,
        H1, N13, -1, O14, -1, -1, N13, O14, k_N13_p_to_O14_reaclib, k_O14_to_p_N13_derived,
        N, N14, -1, N15, -1, -1, N14, N15, k_N14_n_to_N15_reaclib, k_N15_to_n_N14_derived,
        H1, N14, -1, O15, -1, -1, N14, O15, k_N14_p_to_O15_reaclib, k_O15_to_p_N14_derived,
        He4, N14, -1, F18, -1, -1, N14, F18, k_N14_He4_to_F18_reaclib, k_F18_to_He4_N14_derived,
        H1, N15, -1, O16, -1, -1, N15, O16, k_N15_p_to_O16_reaclib, k_O16_to_p_N15_derived,
        He4, N15, -1, F19, -1, -1, N15, F19, k_N15_He4_to_F19_reaclib, k_F19_to_He4_N15_derived,
        N, O14, -1, O15, -1, -1, O14, O15, k_O14_n_to_O15_reaclib, k_O15_to_n_O14_derived,
        He4, O14, -1, Ne18, -1, -1, O14, Ne18, k_O14_He4_to_Ne18_reaclib, k_Ne18_to_He4_O14_derived,
        N, O15, -1, O16, -1, -1, O15, O16, k_O15_n_to_O16_reaclib, k_O16_to_n_O15_derived,
        He4, O15, -1, Ne19, -1, -1, O15, Ne19, k_O15_He4_to_Ne19_reaclib, k_Ne19_to_He4_O15_derived,
        N, O16, -1, O17, -1, -1, O16, O17, k_O16_n_to_O17_reaclib, k_O17_to_n_O16_derived,
        H1, O16, -1, F17, -1, -1, O16, F17, k_O16_p_to_F17_reaclib, k_F17_to_p_O16_derived,
        He4, O16, -1, Ne20, -1, -1, O16, Ne20, k_O16_He4_to_Ne20_reaclib, k_Ne20_to_He4_O16_derived,
        N, O17, -1, O18, -1, -1, O17, O18, k_O17_n_to_O18_reaclib, k_O18_to_n_O17_derived,
        H1, O17, -1, F18, -1, -1, O17, F18, k_O17_p_to_F18_reaclib, k_F18_to_p_O17_derived,
        He4, O17, -1, Ne21, -1, -1, O17, Ne21, k_O17_He4_to_Ne21_reaclib, k_Ne21_to_He4_O17_derived,
        H1, O18, -1, F19, -1, -1, O18, F19, k_O18_p_to_F19_reaclib, k_F19_to_p_O18_derived,
        He4, O18, -1, Ne22, -1, -1, O18, Ne22, k_O18_He4_to_Ne22_reaclib, k_Ne22_to_He4_O18_derived,
        N, F17, -1, F18, -1, -1, F17, F18, k_F17_n_to_F18_reaclib, k_F18_to_n_F17_derived,
        H1, F17, -1, Ne18, -1, -1, F17, Ne18, k_F17_p_to_Ne18_reaclib, k_Ne18_to_p_F17_derived,
        He4, F17, -1, Na21, -1, -1, F17, Na21, k_F17_He4_to_Na21_reaclib, k_Na21_to_He4_F17_derived,
        N, F18, -1, F19, -1, -1, F18, F19, k_F18_n_to_F19_reaclib, k_F19_to_n_F18_derived,
        H1, F18, -1, Ne19, -1, -1, F18, Ne19, k_F18_p_to_Ne19_reaclib, k_Ne19_to_p_F18_derived,
        He4, F18, -1, Na22, -1, -1, F18, Na22, k_F18_He4_to_Na22_reaclib, k_Na22_to_He4_F18_derived,
        H1, F19, -1, Ne20, -1, -1, F19, Ne20, k_F19_p_to_Ne20_reaclib, k_Ne20_to_p_F19_derived,
        He4, F19, -1, Na23, -1, -1, F19, Na23, k_F19_He4_to_Na23_reaclib, k_Na23_to_He4_F19_derived,
        N, Ne18, -1, Ne19, -1, -1, Ne18, Ne19, k_Ne18_n_to_Ne19_reaclib, k_Ne19_to_n_Ne18_derived,
        N, Ne19, -1, Ne20, -1, -1, Ne19, Ne20, k_Ne19_n_to_Ne20_reaclib, k_Ne20_to_n_Ne19_derived,
        He4, Ne19, -1, Mg23, -1, -1, Ne19, Mg23, k_Ne19_He4_to_Mg23_reaclib, k_Mg23_to_He4_Ne19_derived,
        N, Ne20, -1, Ne21, -1, -1, Ne20, Ne21, k_Ne20_n_to_Ne21_reaclib, k_Ne21_to_n_Ne20_derived,
        H1, Ne20, -1, Na21, -1, -1, Ne20, Na21, k_Ne20_p_to_Na21_reaclib, k_Na21_to_p_Ne20_derived,
        He4, Ne20, -1, Mg24, -1, -1, Ne20, Mg24, k_Ne20_He4_to_Mg24_reaclib, k_Mg24_to_He4_Ne20_derived,
        N, Ne21, -1, Ne22, -1, -1, Ne21, Ne22, k_Ne21_n_to_Ne22_reaclib, k_Ne22_to_n_Ne21_derived,
        H1, Ne21, -1, Na22, -1, -1, Ne21, Na22, k_Ne21_p_to_Na22_reaclib, k_Na22_to_p_Ne21_derived,
        He4, Ne21, -1, Mg25, -1, -1, Ne21, Mg25, k_Ne21_He4_to_Mg25_reaclib, k_Mg25_to_He4_Ne21_derived,
        H1, Ne22, -1, Na23, -1, -1, Ne22, Na23, k_Ne22_p_to_Na23_reaclib, k_Na23_to_p_Ne22_derived,
        He4, Ne22, -1, Mg26, -1, -1, Ne22, Mg26, k_Ne22_He4_to_Mg26_reaclib, k_Mg26_to_He4_Ne22_derived,
        N, Na21, -1, Na22, -1, -1, Na21, Na22, k_Na21_n_to_Na22_reaclib, k_Na22_to_n_Na21_derived,
        He4, Na21, -1, Al25, -1, -1, Na21, Al25, k_Na21_He4_to_Al25_reaclib, k_Al25_to_He4_Na21_derived,
        N, Na22, -1, Na23, -1, -1, Na22, Na23, k_Na22_n_to_Na23_reaclib, k_Na23_to_n_Na22_derived,
        H1, Na22, -1, Mg23, -1, -1, Na22, Mg23, k_Na22_p_to_Mg23_reaclib, k_Mg23_to_p_Na22_derived,
        He4, Na22, -1, Al26, -1, -1, Na22, Al26, k_Na22_He4_to_Al26_reaclib, k_Al26_to_He4_Na22_derived,
        H1, Na23, -1, Mg24, -1, -1, Na23, Mg24, k_Na23_p_to_Mg24_reaclib, k_Mg24_to_p_Na23_derived,
        He4, Na23, -1, Al27, -1, -1, Na23, Al27, k_Na23_He4_to_Al27_reaclib, k_Al27_to_He4_Na23_derived,
        N, Mg23, -1, Mg24, -1, -1, Mg23, Mg24, k_Mg23_n_to_Mg24_reaclib, k_Mg24_to_n_Mg23_derived,
        N, Mg24, -1, Mg25, -1, -1, Mg24, Mg25, k_Mg24_n_to_Mg25_reaclib, k_Mg25_to_n_Mg24_derived,
        H1, Mg24, -1, Al25, -1, -1, Mg24, Al25, k_Mg24_p_to_Al25_reaclib, k_Al25_to_p_Mg24_derived,
        He4, Mg24, -1, Si28, -1, -1, Mg24, Si28, k_Mg24_He4_to_Si28_reaclib, k_Si28_to_He4_Mg24_derived,
        N, Mg25, -1, Mg26, -1, -1, Mg25, Mg26, k_Mg25_n_to_Mg26_reaclib, k_Mg26_to_n_Mg25_derived,
        H1, Mg25, -1, Al26, -1, -1, Mg25, Al26, k_Mg25_p_to_Al26_reaclib, k_Al26_to_p_Mg25_derived,
        He4, Mg25, -1, Si29, -1, -1, Mg25, Si29, k_Mg25_He4_to_Si29_reaclib, k_Si29_to_He4_Mg25_derived,
        H1, Mg26, -1, Al27, -1, -1, Mg26, Al27, k_Mg26_p_to_Al27_reaclib, k_Al27_to_p_Mg26_derived,
        He4, Mg26, -1, Si30, -1, -1, Mg26, Si30, k_Mg26_He4_to_Si30_reaclib, k_Si30_to_He4_Mg26_derived,
        N, Al25, -1, Al26, -1, -1, Al25, Al26, k_Al25_n_to_Al26_reaclib, k_Al26_to_n_Al25_derived,
        He4, Al25, -1, P29, -1, -1, Al25, P29, k_Al25_He4_to_P29_reaclib, k_P29_to_He4_Al25_derived,
        N, Al26, -1, Al27, -1, -1, Al26, Al27, k_Al26_n_to_Al27_reaclib, k_Al27_to_n_Al26_derived,
        He4, Al26, -1, P30, -1, -1, Al26, P30, k_Al26_He4_to_P30_reaclib, k_P30_to_He4_Al26_derived,
        H1, Al27, -1, Si28, -1, -1, Al27, Si28, k_Al27_p_to_Si28_reaclib, k_Si28_to_p_Al27_derived,
        He4, Al27, -1, P31, -1, -1, Al27, P31, k_Al27_He4_to_P31_reaclib, k_P31_to_He4_Al27_derived,
        N, Si28, -1, Si29, -1, -1, Si28, Si29, k_Si28_n_to_Si29_reaclib, k_Si29_to_n_Si28_derived,
        H1, Si28, -1, P29, -1, -1, Si28, P29, k_Si28_p_to_P29_reaclib, k_P29_to_p_Si28_derived,
        He4, Si28, -1, S32, -1, -1, Si28, S32, k_Si28_He4_to_S32_reaclib, k_S32_to_He4_Si28_derived,
        N, Si29, -1, Si30, -1, -1, Si29, Si30, k_Si29_n_to_Si30_reaclib, k_Si30_to_n_Si29_derived,
        H1, Si29, -1, P30, -1, -1, Si29, P30, k_Si29_p_to_P30_reaclib, k_P30_to_p_Si29_derived,
        He4, Si29, -1, S33, -1, -1, Si29, S33, k_Si29_He4_to_S33_reaclib, k_S33_to_He4_Si29_derived,
        N, Si30, -1, Si31, -1, -1, Si30, Si31, k_Si30_n_to_Si31_reaclib, k_Si31_to_n_Si30_derived,
        H1, Si30, -1, P31, -1, -1, Si30, P31, k_Si30_p_to_P31_reaclib, k_P31_to_p_Si30_derived,
        He4, Si30, -1, S34, -1, -1, Si30, S34, k_Si30_He4_to_S34_reaclib, k_S34_to_He4_Si30_derived,
        N, Si31, -1, Si32, -1, -1, Si31, Si32, k_Si31_n_to_Si32_reaclib, k_Si32_to_n_Si31_derived,
        H1, Si31, -1, P32, -1, -1, Si31, P32, k_Si31_p_to_P32_reaclib, k_P32_to_p_Si31_derived,
        He4, Si31, -1, S35, -1, -1, Si31, S35, k_Si31_He4_to_S35_reaclib, k_S35_to_He4_Si31_derived,
        H1, Si32, -1, P33, -1, -1, Si32, P33, k_Si32_p_to_P33_reaclib, k_P33_to_p_Si32_derived,
        He4, Si32, -1, S36, -1, -1, Si32, S36, k_Si32_He4_to_S36_reaclib, k_S36_to_He4_Si32_derived,
        N, P29, -1, P30, -1, -1, P29, P30, k_P29_n_to_P30_reaclib, k_P30_to_n_P29_derived,
        He4, P29, -1, Cl33, -1, -1, P29, Cl33, k_P29_He4_to_Cl33_reaclib, k_Cl33_to_He4_P29_derived,
        N, P30, -1, P31, -1, -1, P30, P31, k_P30_n_to_P31_reaclib, k_P31_to_n_P30_derived,
        He4, P30, -1, Cl34, -1, -1, P30, Cl34, k_P30_He4_to_Cl34_reaclib, k_Cl34_to_He4_P30_derived,
        N, P31, -1, P32, -1, -1, P31, P32, k_P31_n_to_P32_reaclib, k_P32_to_n_P31_derived,
        H1, P31, -1, S32, -1, -1, P31, S32, k_P31_p_to_S32_reaclib, k_S32_to_p_P31_derived,
        He4, P31, -1, Cl35, -1, -1, P31, Cl35, k_P31_He4_to_Cl35_reaclib, k_Cl35_to_He4_P31_derived,
        N, P32, -1, P33, -1, -1, P32, P33, k_P32_n_to_P33_reaclib, k_P33_to_n_P32_derived,
        H1, P32, -1, S33, -1, -1, P32, S33, k_P32_p_to_S33_reaclib, k_S33_to_p_P32_derived,
        He4, P32, -1, Cl36, -1, -1, P32, Cl36, k_P32_He4_to_Cl36_reaclib, k_Cl36_to_He4_P32_derived,
        H1, P33, -1, S34, -1, -1, P33, S34, k_P33_p_to_S34_reaclib, k_S34_to_p_P33_derived,
        He4, P33, -1, Cl37, -1, -1, P33, Cl37, k_P33_He4_to_Cl37_reaclib, k_Cl37_to_He4_P33_derived,
        N, S32, -1, S33, -1, -1, S32, S33, k_S32_n_to_S33_reaclib, k_S33_to_n_S32_derived,
        H1, S32, -1, Cl33, -1, -1, S32, Cl33, k_S32_p_to_Cl33_reaclib, k_Cl33_to_p_S32_derived,
        He4, S32, -1, Ar36, -1, -1, S32, Ar36, k_S32_He4_to_Ar36_reaclib, k_Ar36_to_He4_S32_derived,
        N, S33, -1, S34, -1, -1, S33, S34, k_S33_n_to_S34_reaclib, k_S34_to_n_S33_derived,
        H1, S33, -1, Cl34, -1, -1, S33, Cl34, k_S33_p_to_Cl34_reaclib, k_Cl34_to_p_S33_derived,
        He4, S33, -1, Ar37, -1, -1, S33, Ar37, k_S33_He4_to_Ar37_reaclib, k_Ar37_to_He4_S33_derived,
        N, S34, -1, S35, -1, -1, S34, S35, k_S34_n_to_S35_reaclib, k_S35_to_n_S34_derived,
        H1, S34, -1, Cl35, -1, -1, S34, Cl35, k_S34_p_to_Cl35_reaclib, k_Cl35_to_p_S34_derived,
        He4, S34, -1, Ar38, -1, -1, S34, Ar38, k_S34_He4_to_Ar38_reaclib, k_Ar38_to_He4_S34_derived,
        N, S35, -1, S36, -1, -1, S35, S36, k_S35_n_to_S36_reaclib, k_S36_to_n_S35_derived,
        H1, S35, -1, Cl36, -1, -1, S35, Cl36, k_S35_p_to_Cl36_reaclib, k_Cl36_to_p_S35_derived,
        He4, S35, -1, Ar39, -1, -1, S35, Ar39, k_S35_He4_to_Ar39_reaclib, k_Ar39_to_He4_S35_derived,
        H1, S36, -1, Cl37, -1, -1, S36, Cl37, k_S36_p_to_Cl37_reaclib, k_Cl37_to_p_S36_derived,
        He4, S36, -1, Ar40, -1, -1, S36, Ar40, k_S36_He4_to_Ar40_reaclib, k_Ar40_to_He4_S36_derived,
        N, Cl33, -1, Cl34, -1, -1, Cl33, Cl34, k_Cl33_n_to_Cl34_reaclib, k_Cl34_to_n_Cl33_derived,
        He4, Cl33, -1, K37, -1, -1, Cl33, K37, k_Cl33_He4_to_K37_reaclib, k_K37_to_He4_Cl33_derived,
        N, Cl34, -1, Cl35, -1, -1, Cl34, Cl35, k_Cl34_n_to_Cl35_reaclib, k_Cl35_to_n_Cl34_derived,
        He4, Cl34, -1, K38, -1, -1, Cl34, K38, k_Cl34_He4_to_K38_reaclib, k_K38_to_He4_Cl34_derived,
        N, Cl35, -1, Cl36, -1, -1, Cl35, Cl36, k_Cl35_n_to_Cl36_reaclib, k_Cl36_to_n_Cl35_derived,
        H1, Cl35, -1, Ar36, -1, -1, Cl35, Ar36, k_Cl35_p_to_Ar36_reaclib, k_Ar36_to_p_Cl35_derived,
        He4, Cl35, -1, K39, -1, -1, Cl35, K39, k_Cl35_He4_to_K39_reaclib, k_K39_to_He4_Cl35_derived,
        N, Cl36, -1, Cl37, -1, -1, Cl36, Cl37, k_Cl36_n_to_Cl37_reaclib, k_Cl37_to_n_Cl36_derived,
        H1, Cl36, -1, Ar37, -1, -1, Cl36, Ar37, k_Cl36_p_to_Ar37_reaclib, k_Ar37_to_p_Cl36_derived,
        He4, Cl36, -1, K40, -1, -1, Cl36, K40, k_Cl36_He4_to_K40_reaclib, k_K40_to_He4_Cl36_derived,
        H1, Cl37, -1, Ar38, -1, -1, Cl37, Ar38, k_Cl37_p_to_Ar38_reaclib, k_Ar38_to_p_Cl37_derived,
        He4, Cl37, -1, K41, -1, -1, Cl37, K41, k_Cl37_He4_to_K41_reaclib, k_K41_to_He4_Cl37_derived,
        N, Ar36, -1, Ar37, -1, -1, Ar36, Ar37, k_Ar36_n_to_Ar37_reaclib, k_Ar37_to_n_Ar36_derived,
        H1, Ar36, -1, K37, -1, -1, Ar36, K37, k_Ar36_p_to_K37_reaclib, k_K37_to_p_Ar36_derived,
        He4, Ar36, -1, Ca40, -1, -1, Ar36, Ca40, k_Ar36_He4_to_Ca40_reaclib, k_Ca40_to_He4_Ar36_derived,
        N, Ar37, -1, Ar38, -1, -1, Ar37, Ar38, k_Ar37_n_to_Ar38_reaclib, k_Ar38_to_n_Ar37_derived,
        H1, Ar37, -1, K38, -1, -1, Ar37, K38, k_Ar37_p_to_K38_reaclib, k_K38_to_p_Ar37_derived,
        He4, Ar37, -1, Ca41, -1, -1, Ar37, Ca41, k_Ar37_He4_to_Ca41_reaclib, k_Ca41_to_He4_Ar37_derived,
        N, Ar38, -1, Ar39, -1, -1, Ar38, Ar39, k_Ar38_n_to_Ar39_reaclib, k_Ar39_to_n_Ar38_derived,
        H1, Ar38, -1, K39, -1, -1, Ar38, K39, k_Ar38_p_to_K39_reaclib, k_K39_to_p_Ar38_derived,
        He4, Ar38, -1, Ca42, -1, -1, Ar38, Ca42, k_Ar38_He4_to_Ca42_reaclib, k_Ca42_to_He4_Ar38_derived,
        N, Ar39, -1, Ar40, -1, -1, Ar39, Ar40, k_Ar39_n_to_Ar40_reaclib, k_Ar40_to_n_Ar39_derived,
        H1, Ar39, -1, K40, -1, -1, Ar39, K40, k_Ar39_p_to_K40_reaclib, k_K40_to_p_Ar39_derived,
        He4, Ar39, -1, Ca43, -1, -1, Ar39, Ca43, k_Ar39_He4_to_Ca43_reaclib, k_Ca43_to_He4_Ar39_derived,
        H1, Ar40, -1, K41, -1, -1, Ar40, K41, k_Ar40_p_to_K41_reaclib, k_K41_to_p_Ar40_derived,
        He4, Ar40, -1, Ca44, -1, -1, Ar40, Ca44, k_Ar40_He4_to_Ca44_reaclib, k_Ca44_to_He4_Ar40_derived,
        N, K37, -1, K38, -1, -1, K37, K38, k_K37_n_to_K38_reaclib, k_K38_to_n_K37_derived,
        N, K38, -1, K39, -1, -1, K38, K39, k_K38_n_to_K39_reaclib, k_K39_to_n_K38_derived,
        N, K39, -1, K40, -1, -1, K39, K40, k_K39_n_to_K40_reaclib, k_K40_to_n_K39_derived,
        H1, K39, -1, Ca40, -1, -1, K39, Ca40, k_K39_p_to_Ca40_reaclib, k_Ca40_to_p_K39_derived,
        He4, K39, -1, Sc43, -1, -1, K39, Sc43, k_K39_He4_to_Sc43_reaclib, k_Sc43_to_He4_K39_derived,
        N, K40, -1, K41, -1, -1, K40, K41, k_K40_n_to_K41_reaclib, k_K41_to_n_K40_derived,
        H1, K40, -1, Ca41, -1, -1, K40, Ca41, k_K40_p_to_Ca41_reaclib, k_Ca41_to_p_K40_derived,
        He4, K40, -1, Sc44, -1, -1, K40, Sc44, k_K40_He4_to_Sc44_reaclib, k_Sc44_to_He4_K40_derived,
        H1, K41, -1, Ca42, -1, -1, K41, Ca42, k_K41_p_to_Ca42_reaclib, k_Ca42_to_p_K41_derived,
        He4, K41, -1, Sc45, -1, -1, K41, Sc45, k_K41_He4_to_Sc45_reaclib, k_Sc45_to_He4_K41_derived,
        N, Ca40, -1, Ca41, -1, -1, Ca40, Ca41, k_Ca40_n_to_Ca41_reaclib, k_Ca41_to_n_Ca40_derived,
        He4, Ca40, -1, Ti44, -1, -1, Ca40, Ti44, k_Ca40_He4_to_Ti44_reaclib, k_Ti44_to_He4_Ca40_derived,
        N, Ca41, -1, Ca42, -1, -1, Ca41, Ca42, k_Ca41_n_to_Ca42_reaclib, k_Ca42_to_n_Ca41_derived,
        He4, Ca41, -1, Ti45, -1, -1, Ca41, Ti45, k_Ca41_He4_to_Ti45_reaclib, k_Ti45_to_He4_Ca41_derived,
        N, Ca42, -1, Ca43, -1, -1, Ca42, Ca43, k_Ca42_n_to_Ca43_reaclib, k_Ca43_to_n_Ca42_derived,
        H1, Ca42, -1, Sc43, -1, -1, Ca42, Sc43, k_Ca42_p_to_Sc43_reaclib, k_Sc43_to_p_Ca42_derived,
        He4, Ca42, -1, Ti46, -1, -1, Ca42, Ti46, k_Ca42_He4_to_Ti46_reaclib, k_Ti46_to_He4_Ca42_derived,
        N, Ca43, -1, Ca44, -1, -1, Ca43, Ca44, k_Ca43_n_to_Ca44_reaclib, k_Ca44_to_n_Ca43_derived,
        H1, Ca43, -1, Sc44, -1, -1, Ca43, Sc44, k_Ca43_p_to_Sc44_reaclib, k_Sc44_to_p_Ca43_derived,
        He4, Ca43, -1, Ti47, -1, -1, Ca43, Ti47, k_Ca43_He4_to_Ti47_reaclib, k_Ti47_to_He4_Ca43_derived,
        N, Ca44, -1, Ca45, -1, -1, Ca44, Ca45, k_Ca44_n_to_Ca45_reaclib, k_Ca45_to_n_Ca44_derived,
        H1, Ca44, -1, Sc45, -1, -1, Ca44, Sc45, k_Ca44_p_to_Sc45_reaclib, k_Sc45_to_p_Ca44_derived,
        He4, Ca44, -1, Ti48, -1, -1, Ca44, Ti48, k_Ca44_He4_to_Ti48_reaclib, k_Ti48_to_He4_Ca44_derived,
        N, Ca45, -1, Ca46, -1, -1, Ca45, Ca46, k_Ca45_n_to_Ca46_reaclib, k_Ca46_to_n_Ca45_derived,
        H1, Ca45, -1, Sc46, -1, -1, Ca45, Sc46, k_Ca45_p_to_Sc46_reaclib, k_Sc46_to_p_Ca45_derived,
        He4, Ca45, -1, Ti49, -1, -1, Ca45, Ti49, k_Ca45_He4_to_Ti49_reaclib, k_Ti49_to_He4_Ca45_derived,
        N, Ca46, -1, Ca47, -1, -1, Ca46, Ca47, k_Ca46_n_to_Ca47_reaclib, k_Ca47_to_n_Ca46_derived,
        H1, Ca46, -1, Sc47, -1, -1, Ca46, Sc47, k_Ca46_p_to_Sc47_reaclib, k_Sc47_to_p_Ca46_derived,
        He4, Ca46, -1, Ti50, -1, -1, Ca46, Ti50, k_Ca46_He4_to_Ti50_reaclib, k_Ti50_to_He4_Ca46_derived,
        N, Ca47, -1, Ca48, -1, -1, Ca47, Ca48, k_Ca47_n_to_Ca48_reaclib, k_Ca48_to_n_Ca47_derived,
        H1, Ca47, -1, Sc48, -1, -1, Ca47, Sc48, k_Ca47_p_to_Sc48_reaclib, k_Sc48_to_p_Ca47_derived,
        He4, Ca47, -1, Ti51, -1, -1, Ca47, Ti51, k_Ca47_He4_to_Ti51_reaclib, k_Ti51_to_He4_Ca47_derived,
        H1, Ca48, -1, Sc49, -1, -1, Ca48, Sc49, k_Ca48_p_to_Sc49_reaclib, k_Sc49_to_p_Ca48_derived,
        N, Sc43, -1, Sc44, -1, -1, Sc43, Sc44, k_Sc43_n_to_Sc44_reaclib, k_Sc44_to_n_Sc43_derived,
        H1, Sc43, -1, Ti44, -1, -1, Sc43, Ti44, k_Sc43_p_to_Ti44_reaclib, k_Ti44_to_p_Sc43_derived,
        He4, Sc43, -1, V47, -1, -1, Sc43, V47, k_Sc43_He4_to_V47_reaclib, k_V47_to_He4_Sc43_derived,
        N, Sc44, -1, Sc45, -1, -1, Sc44, Sc45, k_Sc44_n_to_Sc45_reaclib, k_Sc45_to_n_Sc44_derived,
        H1, Sc44, -1, Ti45, -1, -1, Sc44, Ti45, k_Sc44_p_to_Ti45_reaclib, k_Ti45_to_p_Sc44_derived,
        He4, Sc44, -1, V48, -1, -1, Sc44, V48, k_Sc44_He4_to_V48_reaclib, k_V48_to_He4_Sc44_derived,
        N, Sc45, -1, Sc46, -1, -1, Sc45, Sc46, k_Sc45_n_to_Sc46_reaclib, k_Sc46_to_n_Sc45_derived,
        H1, Sc45, -1, Ti46, -1, -1, Sc45, Ti46, k_Sc45_p_to_Ti46_reaclib, k_Ti46_to_p_Sc45_derived,
        He4, Sc45, -1, V49, -1, -1, Sc45, V49, k_Sc45_He4_to_V49_reaclib, k_V49_to_He4_Sc45_derived,
        N, Sc46, -1, Sc47, -1, -1, Sc46, Sc47, k_Sc46_n_to_Sc47_reaclib, k_Sc47_to_n_Sc46_derived,
        H1, Sc46, -1, Ti47, -1, -1, Sc46, Ti47, k_Sc46_p_to_Ti47_reaclib, k_Ti47_to_p_Sc46_derived,
        He4, Sc46, -1, V50, -1, -1, Sc46, V50, k_Sc46_He4_to_V50_reaclib, k_V50_to_He4_Sc46_derived,
        N, Sc47, -1, Sc48, -1, -1, Sc47, Sc48, k_Sc47_n_to_Sc48_reaclib, k_Sc48_to_n_Sc47_derived,
        H1, Sc47, -1, Ti48, -1, -1, Sc47, Ti48, k_Sc47_p_to_Ti48_reaclib, k_Ti48_to_p_Sc47_derived,
        He4, Sc47, -1, V51, -1, -1, Sc47, V51, k_Sc47_He4_to_V51_reaclib, k_V51_to_He4_Sc47_derived,
        N, Sc48, -1, Sc49, -1, -1, Sc48, Sc49, k_Sc48_n_to_Sc49_reaclib, k_Sc49_to_n_Sc48_derived,
        H1, Sc48, -1, Ti49, -1, -1, Sc48, Ti49, k_Sc48_p_to_Ti49_reaclib, k_Ti49_to_p_Sc48_derived,
        He4, Sc48, -1, V52, -1, -1, Sc48, V52, k_Sc48_He4_to_V52_reaclib, k_V52_to_He4_Sc48_derived,
        H1, Sc49, -1, Ti50, -1, -1, Sc49, Ti50, k_Sc49_p_to_Ti50_reaclib, k_Ti50_to_p_Sc49_derived,
        N, Ti44, -1, Ti45, -1, -1, Ti44, Ti45, k_Ti44_n_to_Ti45_reaclib, k_Ti45_to_n_Ti44_derived,
        He4, Ti44, -1, Cr48, -1, -1, Ti44, Cr48, k_Ti44_He4_to_Cr48_reaclib, k_Cr48_to_He4_Ti44_derived,
        N, Ti45, -1, Ti46, -1, -1, Ti45, Ti46, k_Ti45_n_to_Ti46_reaclib, k_Ti46_to_n_Ti45_derived,
        H1, Ti45, -1, V46, -1, -1, Ti45, V46, k_Ti45_p_to_V46_reaclib, k_V46_to_p_Ti45_derived,
        He4, Ti45, -1, Cr49, -1, -1, Ti45, Cr49, k_Ti45_He4_to_Cr49_reaclib, k_Cr49_to_He4_Ti45_derived,
        N, Ti46, -1, Ti47, -1, -1, Ti46, Ti47, k_Ti46_n_to_Ti47_reaclib, k_Ti47_to_n_Ti46_derived,
        H1, Ti46, -1, V47, -1, -1, Ti46, V47, k_Ti46_p_to_V47_reaclib, k_V47_to_p_Ti46_derived,
        He4, Ti46, -1, Cr50, -1, -1, Ti46, Cr50, k_Ti46_He4_to_Cr50_reaclib, k_Cr50_to_He4_Ti46_derived,
        N, Ti47, -1, Ti48, -1, -1, Ti47, Ti48, k_Ti47_n_to_Ti48_reaclib, k_Ti48_to_n_Ti47_derived,
        H1, Ti47, -1, V48, -1, -1, Ti47, V48, k_Ti47_p_to_V48_reaclib, k_V48_to_p_Ti47_derived,
        He4, Ti47, -1, Cr51, -1, -1, Ti47, Cr51, k_Ti47_He4_to_Cr51_reaclib, k_Cr51_to_He4_Ti47_derived,
        N, Ti48, -1, Ti49, -1, -1, Ti48, Ti49, k_Ti48_n_to_Ti49_reaclib, k_Ti49_to_n_Ti48_derived,
        H1, Ti48, -1, V49, -1, -1, Ti48, V49, k_Ti48_p_to_V49_reaclib, k_V49_to_p_Ti48_derived,
        He4, Ti48, -1, Cr52, -1, -1, Ti48, Cr52, k_Ti48_He4_to_Cr52_reaclib, k_Cr52_to_He4_Ti48_derived,
        N, Ti49, -1, Ti50, -1, -1, Ti49, Ti50, k_Ti49_n_to_Ti50_reaclib, k_Ti50_to_n_Ti49_derived,
        H1, Ti49, -1, V50, -1, -1, Ti49, V50, k_Ti49_p_to_V50_reaclib, k_V50_to_p_Ti49_derived,
        He4, Ti49, -1, Cr53, -1, -1, Ti49, Cr53, k_Ti49_He4_to_Cr53_reaclib, k_Cr53_to_He4_Ti49_derived,
        N, Ti50, -1, Ti51, -1, -1, Ti50, Ti51, k_Ti50_n_to_Ti51_reaclib, k_Ti51_to_n_Ti50_derived,
        H1, Ti50, -1, V51, -1, -1, Ti50, V51, k_Ti50_p_to_V51_reaclib, k_V51_to_p_Ti50_derived,
        He4, Ti50, -1, Cr54, -1, -1, Ti50, Cr54, k_Ti50_He4_to_Cr54_reaclib, k_Cr54_to_He4_Ti50_derived,
        H1, Ti51, -1, V52, -1, -1, Ti51, V52, k_Ti51_p_to_V52_reaclib, k_V52_to_p_Ti51_derived,
        N, V46, -1, V47, -1, -1, V46, V47, k_V46_n_to_V47_reaclib, k_V47_to_n_V46_derived,
        He4, V46, -1, Mn50, -1, -1, V46, Mn50, k_V46_He4_to_Mn50_reaclib, k_Mn50_to_He4_V46_derived,
        N, V47, -1, V48, -1, -1, V47, V48, k_V47_n_to_V48_reaclib, k_V48_to_n_V47_derived,
        H1, V47, -1, Cr48, -1, -1, V47, Cr48, k_V47_p_to_Cr48_reaclib, k_Cr48_to_p_V47_derived,
        He4, V47, -1, Mn51, -1, -1, V47, Mn51, k_V47_He4_to_Mn51_reaclib, k_Mn51_to_He4_V47_derived,
        N, V48, -1, V49, -1, -1, V48, V49, k_V48_n_to_V49_reaclib, k_V49_to_n_V48_derived,
        H1, V48, -1, Cr49, -1, -1, V48, Cr49, k_V48_p_to_Cr49_reaclib, k_Cr49_to_p_V48_derived,
        He4, V48, -1, Mn52, -1, -1, V48, Mn52, k_V48_He4_to_Mn52_reaclib, k_Mn52_to_He4_V48_derived,
        N, V49, -1, V50, -1, -1, V49, V50, k_V49_n_to_V50_reaclib, k_V50_to_n_V49_derived,
        H1, V49, -1, Cr50, -1, -1, V49, Cr50, k_V49_p_to_Cr50_reaclib, k_Cr50_to_p_V49_derived,
        He4, V49, -1, Mn53, -1, -1, V49, Mn53, k_V49_He4_to_Mn53_reaclib, k_Mn53_to_He4_V49_derived,
        N, V50, -1, V51, -1, -1, V50, V51, k_V50_n_to_V51_reaclib, k_V51_to_n_V50_derived,
        H1, V50, -1, Cr51, -1, -1, V50, Cr51, k_V50_p_to_Cr51_reaclib, k_Cr51_to_p_V50_derived,
        He4, V50, -1, Mn54, -1, -1, V50, Mn54, k_V50_He4_to_Mn54_reaclib, k_Mn54_to_He4_V50_derived,
        N, V51, -1, V52, -1, -1, V51, V52, k_V51_n_to_V52_reaclib, k_V52_to_n_V51_derived,
        H1, V51, -1, Cr52, -1, -1, V51, Cr52, k_V51_p_to_Cr52_reaclib, k_Cr52_to_p_V51_derived,
        He4, V51, -1, Mn55, -1, -1, V51, Mn55, k_V51_He4_to_Mn55_reaclib, k_Mn55_to_He4_V51_derived,
        H1, V52, -1, Cr53, -1, -1, V52, Cr53, k_V52_p_to_Cr53_reaclib, k_Cr53_to_p_V52_derived,
        N, Cr48, -1, Cr49, -1, -1, Cr48, Cr49, k_Cr48_n_to_Cr49_reaclib, k_Cr49_to_n_Cr48_derived,
        He4, Cr48, -1, Fe52, -1, -1, Cr48, Fe52, k_Cr48_He4_to_Fe52_reaclib, k_Fe52_to_He4_Cr48_derived,
        N, Cr49, -1, Cr50, -1, -1, Cr49, Cr50, k_Cr49_n_to_Cr50_reaclib, k_Cr50_to_n_Cr49_derived,
        H1, Cr49, -1, Mn50, -1, -1, Cr49, Mn50, k_Cr49_p_to_Mn50_reaclib, k_Mn50_to_p_Cr49_derived,
        He4, Cr49, -1, Fe53, -1, -1, Cr49, Fe53, k_Cr49_He4_to_Fe53_reaclib, k_Fe53_to_He4_Cr49_derived,
        N, Cr50, -1, Cr51, -1, -1, Cr50, Cr51, k_Cr50_n_to_Cr51_reaclib, k_Cr51_to_n_Cr50_derived,
        H1, Cr50, -1, Mn51, -1, -1, Cr50, Mn51, k_Cr50_p_to_Mn51_reaclib, k_Mn51_to_p_Cr50_derived,
        He4, Cr50, -1, Fe54, -1, -1, Cr50, Fe54, k_Cr50_He4_to_Fe54_reaclib, k_Fe54_to_He4_Cr50_derived,
        N, Cr51, -1, Cr52, -1, -1, Cr51, Cr52, k_Cr51_n_to_Cr52_reaclib, k_Cr52_to_n_Cr51_derived,
        H1, Cr51, -1, Mn52, -1, -1, Cr51, Mn52, k_Cr51_p_to_Mn52_reaclib, k_Mn52_to_p_Cr51_derived,
        He4, Cr51, -1, Fe55, -1, -1, Cr51, Fe55, k_Cr51_He4_to_Fe55_reaclib, k_Fe55_to_He4_Cr51_derived,
        N, Cr52, -1, Cr53, -1, -1, Cr52, Cr53, k_Cr52_n_to_Cr53_reaclib, k_Cr53_to_n_Cr52_derived,
        H1, Cr52, -1, Mn53, -1, -1, Cr52, Mn53, k_Cr52_p_to_Mn53_reaclib, k_Mn53_to_p_Cr52_derived,
        He4, Cr52, -1, Fe56, -1, -1, Cr52, Fe56, k_Cr52_He4_to_Fe56_reaclib, k_Fe56_to_He4_Cr52_derived,
        N, Cr53, -1, Cr54, -1, -1, Cr53, Cr54, k_Cr53_n_to_Cr54_reaclib, k_Cr54_to_n_Cr53_derived,
        H1, Cr53, -1, Mn54, -1, -1, Cr53, Mn54, k_Cr53_p_to_Mn54_reaclib, k_Mn54_to_p_Cr53_derived,
        He4, Cr53, -1, Fe57, -1, -1, Cr53, Fe57, k_Cr53_He4_to_Fe57_reaclib, k_Fe57_to_He4_Cr53_derived,
        H1, Cr54, -1, Mn55, -1, -1, Cr54, Mn55, k_Cr54_p_to_Mn55_reaclib, k_Mn55_to_p_Cr54_derived,
        He4, Cr54, -1, Fe58, -1, -1, Cr54, Fe58, k_Cr54_He4_to_Fe58_reaclib, k_Fe58_to_He4_Cr54_derived,
        N, Mn50, -1, Mn51, -1, -1, Mn50, Mn51, k_Mn50_n_to_Mn51_reaclib, k_Mn51_to_n_Mn50_derived,
        He4, Mn50, -1, Co54, -1, -1, Mn50, Co54, k_Mn50_He4_to_Co54_reaclib, k_Co54_to_He4_Mn50_derived,
        N, Mn51, -1, Mn52, -1, -1, Mn51, Mn52, k_Mn51_n_to_Mn52_reaclib, k_Mn52_to_n_Mn51_derived,
        H1, Mn51, -1, Fe52, -1, -1, Mn51, Fe52, k_Mn51_p_to_Fe52_reaclib, k_Fe52_to_p_Mn51_derived,
        He4, Mn51, -1, Co55, -1, -1, Mn51, Co55, k_Mn51_He4_to_Co55_reaclib, k_Co55_to_He4_Mn51_derived,
        N, Mn52, -1, Mn53, -1, -1, Mn52, Mn53, k_Mn52_n_to_Mn53_reaclib, k_Mn53_to_n_Mn52_derived,
        H1, Mn52, -1, Fe53, -1, -1, Mn52, Fe53, k_Mn52_p_to_Fe53_reaclib, k_Fe53_to_p_Mn52_derived,
        He4, Mn52, -1, Co56, -1, -1, Mn52, Co56, k_Mn52_He4_to_Co56_reaclib, k_Co56_to_He4_Mn52_derived,
        N, Mn53, -1, Mn54, -1, -1, Mn53, Mn54, k_Mn53_n_to_Mn54_reaclib, k_Mn54_to_n_Mn53_derived,
        H1, Mn53, -1, Fe54, -1, -1, Mn53, Fe54, k_Mn53_p_to_Fe54_reaclib, k_Fe54_to_p_Mn53_derived,
        He4, Mn53, -1, Co57, -1, -1, Mn53, Co57, k_Mn53_He4_to_Co57_reaclib, k_Co57_to_He4_Mn53_derived,
        N, Mn54, -1, Mn55, -1, -1, Mn54, Mn55, k_Mn54_n_to_Mn55_reaclib, k_Mn55_to_n_Mn54_derived,
        H1, Mn54, -1, Fe55, -1, -1, Mn54, Fe55, k_Mn54_p_to_Fe55_reaclib, k_Fe55_to_p_Mn54_derived,
        He4, Mn54, -1, Co58, -1, -1, Mn54, Co58, k_Mn54_He4_to_Co58_reaclib, k_Co58_to_He4_Mn54_derived,
        H1, Mn55, -1, Fe56, -1, -1, Mn55, Fe56, k_Mn55_p_to_Fe56_reaclib, k_Fe56_to_p_Mn55_derived,
        He4, Mn55, -1, Co59, -1, -1, Mn55, Co59, k_Mn55_He4_to_Co59_reaclib, k_Co59_to_He4_Mn55_derived,
        N, Fe52, -1, Fe53, -1, -1, Fe52, Fe53, k_Fe52_n_to_Fe53_reaclib, k_Fe53_to_n_Fe52_derived,
        H1, Fe52, -1, Co53, -1, -1, Fe52, Co53, k_Fe52_p_to_Co53_reaclib, k_Co53_to_p_Fe52_derived,
        He4, Fe52, -1, Ni56, -1, -1, Fe52, Ni56, k_Fe52_He4_to_Ni56_reaclib, k_Ni56_to_He4_Fe52_derived,
        N, Fe53, -1, Fe54, -1, -1, Fe53, Fe54, k_Fe53_n_to_Fe54_reaclib, k_Fe54_to_n_Fe53_derived,
        H1, Fe53, -1, Co54, -1, -1, Fe53, Co54, k_Fe53_p_to_Co54_reaclib, k_Co54_to_p_Fe53_derived,
        He4, Fe53, -1, Ni57, -1, -1, Fe53, Ni57, k_Fe53_He4_to_Ni57_reaclib, k_Ni57_to_He4_Fe53_derived,
        N, Fe54, -1, Fe55, -1, -1, Fe54, Fe55, k_Fe54_n_to_Fe55_reaclib, k_Fe55_to_n_Fe54_derived,
        H1, Fe54, -1, Co55, -1, -1, Fe54, Co55, k_Fe54_p_to_Co55_reaclib, k_Co55_to_p_Fe54_derived,
        He4, Fe54, -1, Ni58, -1, -1, Fe54, Ni58, k_Fe54_He4_to_Ni58_reaclib, k_Ni58_to_He4_Fe54_derived,
        N, Fe55, -1, Fe56, -1, -1, Fe55, Fe56, k_Fe55_n_to_Fe56_reaclib, k_Fe56_to_n_Fe55_derived,
        H1, Fe55, -1, Co56, -1, -1, Fe55, Co56, k_Fe55_p_to_Co56_reaclib, k_Co56_to_p_Fe55_derived,
        He4, Fe55, -1, Ni59, -1, -1, Fe55, Ni59, k_Fe55_He4_to_Ni59_reaclib, k_Ni59_to_He4_Fe55_derived,
        N, Fe56, -1, Fe57, -1, -1, Fe56, Fe57, k_Fe56_n_to_Fe57_reaclib, k_Fe57_to_n_Fe56_derived,
        H1, Fe56, -1, Co57, -1, -1, Fe56, Co57, k_Fe56_p_to_Co57_reaclib, k_Co57_to_p_Fe56_derived,
        He4, Fe56, -1, Ni60, -1, -1, Fe56, Ni60, k_Fe56_He4_to_Ni60_reaclib, k_Ni60_to_He4_Fe56_derived,
        N, Fe57, -1, Fe58, -1, -1, Fe57, Fe58, k_Fe57_n_to_Fe58_reaclib, k_Fe58_to_n_Fe57_derived,
        H1, Fe57, -1, Co58, -1, -1, Fe57, Co58, k_Fe57_p_to_Co58_reaclib, k_Co58_to_p_Fe57_derived,
        He4, Fe57, -1, Ni61, -1, -1, Fe57, Ni61, k_Fe57_He4_to_Ni61_reaclib, k_Ni61_to_He4_Fe57_derived,
        H1, Fe58, -1, Co59, -1, -1, Fe58, Co59, k_Fe58_p_to_Co59_reaclib, k_Co59_to_p_Fe58_derived,
        He4, Fe58, -1, Ni62, -1, -1, Fe58, Ni62, k_Fe58_He4_to_Ni62_reaclib, k_Ni62_to_He4_Fe58_derived,
        N, Co53, -1, Co54, -1, -1, Co53, Co54, k_Co53_n_to_Co54_reaclib, k_Co54_to_n_Co53_derived,
        He4, Co53, -1, Cu57, -1, -1, Co53, Cu57, k_Co53_He4_to_Cu57_reaclib, k_Cu57_to_He4_Co53_derived,
        N, Co54, -1, Co55, -1, -1, Co54, Co55, k_Co54_n_to_Co55_reaclib, k_Co55_to_n_Co54_derived,
        He4, Co54, -1, Cu58, -1, -1, Co54, Cu58, k_Co54_He4_to_Cu58_reaclib, k_Cu58_to_He4_Co54_derived,
        N, Co55, -1, Co56, -1, -1, Co55, Co56, k_Co55_n_to_Co56_reaclib, k_Co56_to_n_Co55_derived,
        H1, Co55, -1, Ni56, -1, -1, Co55, Ni56, k_Co55_p_to_Ni56_reaclib, k_Ni56_to_p_Co55_derived,
        He4, Co55, -1, Cu59, -1, -1, Co55, Cu59, k_Co55_He4_to_Cu59_reaclib, k_Cu59_to_He4_Co55_derived,
        N, Co56, -1, Co57, -1, -1, Co56, Co57, k_Co56_n_to_Co57_reaclib, k_Co57_to_n_Co56_derived,
        H1, Co56, -1, Ni57, -1, -1, Co56, Ni57, k_Co56_p_to_Ni57_reaclib, k_Ni57_to_p_Co56_derived,
        He4, Co56, -1, Cu60, -1, -1, Co56, Cu60, k_Co56_He4_to_Cu60_reaclib, k_Cu60_to_He4_Co56_derived,
        N, Co57, -1, Co58, -1, -1, Co57, Co58, k_Co57_n_to_Co58_reaclib, k_Co58_to_n_Co57_derived,
        H1, Co57, -1, Ni58, -1, -1, Co57, Ni58, k_Co57_p_to_Ni58_reaclib, k_Ni58_to_p_Co57_derived,
        He4, Co57, -1, Cu61, -1, -1, Co57, Cu61, k_Co57_He4_to_Cu61_reaclib, k_Cu61_to_He4_Co57_derived,
        N, Co58, -1, Co59, -1, -1, Co58, Co59, k_Co58_n_to_Co59_reaclib, k_Co59_to_n_Co58_derived,
        H1, Co58, -1, Ni59, -1, -1, Co58, Ni59, k_Co58_p_to_Ni59_reaclib, k_Ni59_to_p_Co58_derived,
        He4, Co58, -1, Cu62, -1, -1, Co58, Cu62, k_Co58_He4_to_Cu62_reaclib, k_Cu62_to_He4_Co58_derived,
        H1, Co59, -1, Ni60, -1, -1, Co59, Ni60, k_Co59_p_to_Ni60_reaclib, k_Ni60_to_p_Co59_derived,
        He4, Co59, -1, Cu63, -1, -1, Co59, Cu63, k_Co59_He4_to_Cu63_reaclib, k_Cu63_to_He4_Co59_derived,
        N, Ni56, -1, Ni57, -1, -1, Ni56, Ni57, k_Ni56_n_to_Ni57_reaclib, k_Ni57_to_n_Ni56_derived,
        H1, Ni56, -1, Cu57, -1, -1, Ni56, Cu57, k_Ni56_p_to_Cu57_reaclib, k_Cu57_to_p_Ni56_derived,
        He4, Ni56, -1, Zn60, -1, -1, Ni56, Zn60, k_Ni56_He4_to_Zn60_reaclib, k_Zn60_to_He4_Ni56_derived,
        N, Ni57, -1, Ni58, -1, -1, Ni57, Ni58, k_Ni57_n_to_Ni58_reaclib, k_Ni58_to_n_Ni57_derived,
        H1, Ni57, -1, Cu58, -1, -1, Ni57, Cu58, k_Ni57_p_to_Cu58_reaclib, k_Cu58_to_p_Ni57_derived,
        He4, Ni57, -1, Zn61, -1, -1, Ni57, Zn61, k_Ni57_He4_to_Zn61_reaclib, k_Zn61_to_He4_Ni57_derived,
        N, Ni58, -1, Ni59, -1, -1, Ni58, Ni59, k_Ni58_n_to_Ni59_reaclib, k_Ni59_to_n_Ni58_derived,
        H1, Ni58, -1, Cu59, -1, -1, Ni58, Cu59, k_Ni58_p_to_Cu59_reaclib, k_Cu59_to_p_Ni58_derived,
        He4, Ni58, -1, Zn62, -1, -1, Ni58, Zn62, k_Ni58_He4_to_Zn62_reaclib, k_Zn62_to_He4_Ni58_derived,
        N, Ni59, -1, Ni60, -1, -1, Ni59, Ni60, k_Ni59_n_to_Ni60_reaclib, k_Ni60_to_n_Ni59_derived,
        H1, Ni59, -1, Cu60, -1, -1, Ni59, Cu60, k_Ni59_p_to_Cu60_reaclib, k_Cu60_to_p_Ni59_derived,
        He4, Ni59, -1, Zn63, -1, -1, Ni59, Zn63, k_Ni59_He4_to_Zn63_reaclib, k_Zn63_to_He4_Ni59_derived,
        N, Ni60, -1, Ni61, -1, -1, Ni60, Ni61, k_Ni60_n_to_Ni61_reaclib, k_Ni61_to_n_Ni60_derived,
        H1, Ni60, -1, Cu61, -1, -1, Ni60, Cu61, k_Ni60_p_to_Cu61_reaclib, k_Cu61_to_p_Ni60_derived,
        He4, Ni60, -1, Zn64, -1, -1, Ni60, Zn64, k_Ni60_He4_to_Zn64_reaclib, k_Zn64_to_He4_Ni60_derived,
        N, Ni61, -1, Ni62, -1, -1, Ni61, Ni62, k_Ni61_n_to_Ni62_reaclib, k_Ni62_to_n_Ni61_derived,
        H1, Ni61, -1, Cu62, -1, -1, Ni61, Cu62, k_Ni61_p_to_Cu62_reaclib, k_Cu62_to_p_Ni61_derived,
        He4, Ni61, -1, Zn65, -1, -1, Ni61, Zn65, k_Ni61_He4_to_Zn65_reaclib, k_Zn65_to_He4_Ni61_derived,
        N, Ni62, -1, Ni63, -1, -1, Ni62, Ni63, k_Ni62_n_to_Ni63_reaclib, k_Ni63_to_n_Ni62_derived,
        H1, Ni62, -1, Cu63, -1, -1, Ni62, Cu63, k_Ni62_p_to_Cu63_reaclib, k_Cu63_to_p_Ni62_derived,
        He4, Ni62, -1, Zn66, -1, -1, Ni62, Zn66, k_Ni62_He4_to_Zn66_reaclib, k_Zn66_to_He4_Ni62_derived,
        N, Ni63, -1, Ni64, -1, -1, Ni63, Ni64, k_Ni63_n_to_Ni64_reaclib, k_Ni64_to_n_Ni63_derived,
        H1, Ni63, -1, Cu64, -1, -1, Ni63, Cu64, k_Ni63_p_to_Cu64_reaclib, k_Cu64_to_p_Ni63_derived,
        H1, Ni64, -1, Cu65, -1, -1, Ni64, Cu65, k_Ni64_p_to_Cu65_reaclib, k_Cu65_to_p_Ni64_derived,
        N, Cu57, -1, Cu58, -1, -1, Cu57, Cu58, k_Cu57_n_to_Cu58_reaclib, k_Cu58_to_n_Cu57_derived,
        N, Cu58, -1, Cu59, -1, -1, Cu58, Cu59, k_Cu58_n_to_Cu59_reaclib, k_Cu59_to_n_Cu58_derived,
        H1, Cu58, -1, Zn59, -1, -1, Cu58, Zn59, k_Cu58_p_to_Zn59_reaclib, k_Zn59_to_p_Cu58_derived,
        He4, Cu58, -1, Ga62, -1, -1, Cu58, Ga62, k_Cu58_He4_to_Ga62_reaclib, k_Ga62_to_He4_Cu58_derived,
        N, Cu59, -1, Cu60, -1, -1, Cu59, Cu60, k_Cu59_n_to_Cu60_reaclib, k_Cu60_to_n_Cu59_derived,
        H1, Cu59, -1, Zn60, -1, -1, Cu59, Zn60, k_Cu59_p_to_Zn60_reaclib, k_Zn60_to_p_Cu59_derived,
        He4, Cu59, -1, Ga63, -1, -1, Cu59, Ga63, k_Cu59_He4_to_Ga63_reaclib, k_Ga63_to_He4_Cu59_derived,
        N, Cu60, -1, Cu61, -1, -1, Cu60, Cu61, k_Cu60_n_to_Cu61_reaclib, k_Cu61_to_n_Cu60_derived,
        H1, Cu60, -1, Zn61, -1, -1, Cu60, Zn61, k_Cu60_p_to_Zn61_reaclib, k_Zn61_to_p_Cu60_derived,
        He4, Cu60, -1, Ga64, -1, -1, Cu60, Ga64, k_Cu60_He4_to_Ga64_reaclib, k_Ga64_to_He4_Cu60_derived,
        N, Cu61, -1, Cu62, -1, -1, Cu61, Cu62, k_Cu61_n_to_Cu62_reaclib, k_Cu62_to_n_Cu61_derived,
        H1, Cu61, -1, Zn62, -1, -1, Cu61, Zn62, k_Cu61_p_to_Zn62_reaclib, k_Zn62_to_p_Cu61_derived,
        N, Cu62, -1, Cu63, -1, -1, Cu62, Cu63, k_Cu62_n_to_Cu63_reaclib, k_Cu63_to_n_Cu62_derived,
        H1, Cu62, -1, Zn63, -1, -1, Cu62, Zn63, k_Cu62_p_to_Zn63_reaclib, k_Zn63_to_p_Cu62_derived,
        N, Cu63, -1, Cu64, -1, -1, Cu63, Cu64, k_Cu63_n_to_Cu64_reaclib, k_Cu64_to_n_Cu63_derived,
        H1, Cu63, -1, Zn64, -1, -1, Cu63, Zn64, k_Cu63_p_to_Zn64_reaclib, k_Zn64_to_p_Cu63_derived,
        N, Cu64, -1, Cu65, -1, -1, Cu64, Cu65, k_Cu64_n_to_Cu65_reaclib, k_Cu65_to_n_Cu64_derived,
        H1, Cu64, -1, Zn65, -1, -1, Cu64, Zn65, k_Cu64_p_to_Zn65_reaclib, k_Zn65_to_p_Cu64_derived,
        H1, Cu65, -1, Zn66, -1, -1, Cu65, Zn66, k_Cu65_p_to_Zn66_reaclib, k_Zn66_to_p_Cu65_derived,
        N, Zn59, -1, Zn60, -1, -1, Zn59, Zn60, k_Zn59_n_to_Zn60_reaclib, k_Zn60_to_n_Zn59_derived,
        He4, Zn59, -1, Ge63, -1, -1, Zn59, Ge63, k_Zn59_He4_to_Ge63_reaclib, k_Ge63_to_He4_Zn59_derived,
        N, Zn60, -1, Zn61, -1, -1, Zn60, Zn61, k_Zn60_n_to_Zn61_reaclib, k_Zn61_to_n_Zn60_derived,
        He4, Zn60, -1, Ge64, -1, -1, Zn60, Ge64, k_Zn60_He4_to_Ge64_reaclib, k_Ge64_to_He4_Zn60_derived,
        N, Zn61, -1, Zn62, -1, -1, Zn61, Zn62, k_Zn61_n_to_Zn62_reaclib, k_Zn62_to_n_Zn61_derived,
        H1, Zn61, -1, Ga62, -1, -1, Zn61, Ga62, k_Zn61_p_to_Ga62_reaclib, k_Ga62_to_p_Zn61_derived,
        N, Zn62, -1, Zn63, -1, -1, Zn62, Zn63, k_Zn62_n_to_Zn63_reaclib, k_Zn63_to_n_Zn62_derived,
        H1, Zn62, -1, Ga63, -1, -1, Zn62, Ga63, k_Zn62_p_to_Ga63_reaclib, k_Ga63_to_p_Zn62_derived,
        N, Zn63, -1, Zn64, -1, -1, Zn63, Zn64, k_Zn63_n_to_Zn64_reaclib, k_Zn64_to_n_Zn63_derived,
        H1, Zn63, -1, Ga64, -1, -1, Zn63, Ga64, k_Zn63_p_to_Ga64_reaclib, k_Ga64_to_p_Zn63_derived,
        N, Zn64, -1, Zn65, -1, -1, Zn64, Zn65, k_Zn64_n_to_Zn65_reaclib, k_Zn65_to_n_Zn64_derived,
        N, Zn65, -1, Zn66, -1, -1, Zn65, Zn66, k_Zn65_n_to_Zn66_reaclib, k_Zn66_to_n_Zn65_derived,
        N, Ga62, -1, Ga63, -1, -1, Ga62, Ga63, k_Ga62_n_to_Ga63_reaclib, k_Ga63_to_n_Ga62_derived,
        H1, Ga62, -1, Ge63, -1, -1, Ga62, Ge63, k_Ga62_p_to_Ge63_reaclib, k_Ge63_to_p_Ga62_derived,
        N, Ga63, -1, Ga64, -1, -1, Ga63, Ga64, k_Ga63_n_to_Ga64_reaclib, k_Ga64_to_n_Ga63_derived,
        H1, Ga63, -1, Ge64, -1, -1, Ga63, Ge64, k_Ga63_p_to_Ge64_reaclib, k_Ge64_to_p_Ga63_derived,
        N, Ge63, -1, Ge64, -1, -1, Ge63, Ge64, k_Ge63_n_to_Ge64_reaclib, k_Ge64_to_n_Ge63_derived,
        H2, He3, -1, H1, He4, -1, H2, He3, k_He3_d_to_p_He4_reaclib, k_He4_p_to_d_He3_derived,
        H1, Li6, -1, He4, He3, -1, Li6, He3, k_Li6_p_to_He4_He3_reaclib, k_He3_He4_to_p_Li6_derived,
        H1, Li7, -1, He4, He4, -1, Li7, -1, k_Li7_p_to_He4_He4_reaclib, k_He4_He4_to_p_Li7_derived,
        N, Be7, -1, H1, Li7, -1, Be7, Li7, k_Be7_n_to_p_Li7_reaclib, k_Li7_p_to_n_Be7_derived,
        N, Be7, -1, He4, He4, -1, Be7, -1, k_Be7_n_to_He4_He4_reaclib, k_He4_He4_to_n_Be7_derived,
        H1, Be9, -1, He4, Li6, -1, Be9, Li6, k_Be9_p_to_He4_Li6_reaclib, k_Li6_He4_to_p_Be9_derived,
        He4, Be9, -1, N, C12, -1, Be9, C12, k_Be9_He4_to_n_C12_reaclib, k_C12_n_to_He4_Be9_derived,
        N, B10, -1, He4, Li7, -1, B10, Li7, k_B10_n_to_He4_Li7_reaclib, k_Li7_He4_to_n_B10_derived,
        H1, B10, -1, He4, Be7, -1, B10, Be7, k_B10_p_to_He4_Be7_reaclib, k_Be7_He4_to_p_B10_derived,
        He4, B10, -1, N, N13, -1, B10, N13, k_B10_He4_to_n_N13_reaclib, k_N13_n_to_He4_B10_derived,
        He4, B10, -1, H1, C13, -1, B10, C13, k_B10_He4_to_p_C13_reaclib, k_C13_p_to_He4_B10_derived,
        He4, B11, -1, N, N14, -1, B11, N14, k_B11_He4_to_n_N14_reaclib, k_N14_n_to_He4_B11_derived,
        He4, B11, -1, H1, C14, -1, B11, C14, k_B11_He4_to_p_C14_reaclib, k_C14_p_to_He4_B11_derived,
        He4, C13, -1, N, O16, -1, C13, O16, k_C13_He4_to_n_O16_reaclib, k_O16_n_to_He4_C13_derived,
        He4, N13, -1, H1, O16, -1, N13, O16, k_N13_He4_to_p_O16_reaclib, k_O16_p_to_He4_N13_derived,
        N, N14, -1, H1, C14, -1, N14, C14, k_N14_n_to_p_C14_reaclib, k_C14_p_to_n_N14_derived,
        H1, N15, -1, He4, C12, -1, N15, C12, k_N15_p_to_He4_C12_reaclib, k_C12_He4_to_p_N15_derived,
        He4, O14, -1, H1, F17, -1, O14, F17, k_O14_He4_to_p_F17_reaclib, k_F17_p_to_He4_O14_derived,
        N, O15, -1, He4, C12, -1, O15, C12, k_O15_n_to_He4_C12_reaclib, k_C12_He4_to_n_O15_derived,
        N, O17, -1, He4, C14, -1, O17, C14, k_O17_n_to_He4_C14_reaclib, k_C14_He4_to_n_O17_derived,
        H1, O17, -1, He4, N14, -1, O17, N14, k_O17_p_to_He4_N14_reaclib, k_N14_He4_to_p_O17_derived,
        He4, O17, -1, N, Ne20, -1, O17, Ne20, k_O17_He4_to_n_Ne20_reaclib, k_Ne20_n_to_He4_O17_derived,
        H1, O18, -1, He4, N15, -1, O18, N15, k_O18_p_to_He4_N15_reaclib, k_N15_He4_to_p_O18_derived,
        N, F17, -1, H1, O17, -1, F17, O17, k_F17_n_to_p_O17_reaclib, k_O17_p_to_n_F17_derived,
        N, F18, -1, H1, O18, -1, F18, O18, k_F18_n_to_p_O18_reaclib, k_O18_p_to_n_F18_derived,
        N, F18, -1, He4, N15, -1, F18, N15, k_F18_n_to_He4_N15_reaclib, k_N15_He4_to_n_F18_derived,
        H1, F18, -1, He4, O15, -1, F18, O15, k_F18_p_to_He4_O15_reaclib, k_O15_He4_to_p_F18_derived,
        He4, F18, -1, H1, Ne21, -1, F18, Ne21, k_F18_He4_to_p_Ne21_reaclib, k_Ne21_p_to_He4_F18_derived,
        H1, F19, -1, He4, O16, -1, F19, O16, k_F19_p_to_He4_O16_reaclib, k_O16_He4_to_p_F19_derived,
        He4, F19, -1, H1, Ne22, -1, F19, Ne22, k_F19_He4_to_p_Ne22_reaclib, k_Ne22_p_to_He4_F19_derived,
        N, Ne18, -1, H1, F18, -1, Ne18, F18, k_Ne18_n_to_p_F18_reaclib, k_F18_p_to_n_Ne18_derived,
        N, Ne18, -1, He4, O15, -1, Ne18, O15, k_Ne18_n_to_He4_O15_reaclib, k_O15_He4_to_n_Ne18_derived,
        He4, Ne18, -1, H1, Na21, -1, Ne18, Na21, k_Ne18_He4_to_p_Na21_reaclib, k_Na21_p_to_He4_Ne18_derived,
        N, Ne19, -1, He4, O16, -1, Ne19, O16, k_Ne19_n_to_He4_O16_reaclib, k_O16_He4_to_n_Ne19_derived,
        He4, Ne19, -1, H1, Na22, -1, Ne19, Na22, k_Ne19_He4_to_p_Na22_reaclib, k_Na22_p_to_He4_Ne19_derived,
        He4, Ne21, -1, N, Mg24, -1, Ne21, Mg24, k_Ne21_He4_to_n_Mg24_reaclib, k_Mg24_n_to_He4_Ne21_derived,
        N, Na21, -1, H1, Ne21, -1, Na21, Ne21, k_Na21_n_to_p_Ne21_reaclib, k_Ne21_p_to_n_Na21_derived,
        N, Na21, -1, He4, F18, -1, Na21, F18, k_Na21_n_to_He4_F18_reaclib, k_F18_He4_to_n_Na21_derived,
        N, Na22, -1, H1, Ne22, -1, Na22, Ne22, k_Na22_n_to_p_Ne22_reaclib, k_Ne22_p_to_n_Na22_derived,
        N, Na22, -1, He4, F19, -1, Na22, F19, k_Na22_n_to_He4_F19_reaclib, k_F19_He4_to_n_Na22_derived,
        He4, Na22, -1, H1, Mg25, -1, Na22, Mg25, k_Na22_He4_to_p_Mg25_reaclib, k_Mg25_p_to_He4_Na22_derived,
        H1, Na23, -1, He4, Ne20, -1, Na23, Ne20, k_Na23_p_to_He4_Ne20_reaclib, k_Ne20_He4_to_p_Na23_derived,
        He4, Na23, -1, H1, Mg26, -1, Na23, Mg26, k_Na23_He4_to_p_Mg26_reaclib, k_Mg26_p_to_He4_Na23_derived,
        N, Mg23, -1, He4, Ne20, -1, Mg23, Ne20, k_Mg23_n_to_He4_Ne20_reaclib, k_Ne20_He4_to_n_Mg23_derived,
        He4, Mg23, -1, H1, Al26, -1, Mg23, Al26, k_Mg23_He4_to_p_Al26_reaclib, k_Al26_p_to_He4_Mg23_derived,
        He4, Mg25, -1, N, Si28, -1, Mg25, Si28, k_Mg25_He4_to_n_Si28_reaclib, k_Si28_n_to_He4_Mg25_derived,
        He4, Mg26, -1, N, Si29, -1, Mg26, Si29, k_Mg26_He4_to_n_Si29_reaclib, k_Si29_n_to_He4_Mg26_derived,
        N, Al25, -1, H1, Mg25, -1, Al25, Mg25, k_Al25_n_to_p_Mg25_reaclib, k_Mg25_p_to_n_Al25_derived,
        N, Al25, -1, He4, Na22, -1, Al25, Na22, k_Al25_n_to_He4_Na22_reaclib, k_Na22_He4_to_n_Al25_derived,
        He4, Al25, -1, H1, Si28, -1, Al25, Si28, k_Al25_He4_to_p_Si28_reaclib, k_Si28_p_to_He4_Al25_derived,
        N, Al26, -1, H1, Mg26, -1, Al26, Mg26, k_Al26_n_to_p_Mg26_reaclib, k_Mg26_p_to_n_Al26_derived,
        N, Al26, -1, He4, Na23, -1, Al26, Na23, k_Al26_n_to_He4_Na23_reaclib, k_Na23_He4_to_n_Al26_derived,
        He4, Al26, -1, H1, Si29, -1, Al26, Si29, k_Al26_He4_to_p_Si29_reaclib, k_Si29_p_to_He4_Al26_derived,
        H1, Al27, -1, He4, Mg24, -1, Al27, Mg24, k_Al27_p_to_He4_Mg24_reaclib, k_Mg24_He4_to_p_Al27_derived,
        He4, Al27, -1, H1, Si30, -1, Al27, Si30, k_Al27_He4_to_p_Si30_reaclib, k_Si30_p_to_He4_Al27_derived,
        H1, Si31, -1, N, P31, -1, Si31, P31, k_Si31_p_to_n_P31_reaclib, k_P31_n_to_p_Si31_derived,
        He4, Si31, -1, N, S34, -1, Si31, S34, k_Si31_He4_to_n_S34_reaclib, k_S34_n_to_He4_Si31_derived,
        N, P29, -1, H1, Si29, -1, P29, Si29, k_P29_n_to_p_Si29_reaclib, k_Si29_p_to_n_P29_derived,
        N, P29, -1, He4, Al26, -1, P29, Al26, k_P29_n_to_He4_Al26_reaclib, k_Al26_He4_to_n_P29_derived,
        He4, P29, -1, H1, S32, -1, P29, S32, k_P29_He4_to_p_S32_reaclib, k_S32_p_to_He4_P29_derived,
        N, P30, -1, H1, Si30, -1, P30, Si30, k_P30_n_to_p_Si30_reaclib, k_Si30_p_to_n_P30_derived,
        He4, P30, -1, H1, S33, -1, P30, S33, k_P30_He4_to_p_S33_reaclib, k_S33_p_to_He4_P30_derived,
        H1, P31, -1, He4, Si28, -1, P31, Si28, k_P31_p_to_He4_Si28_reaclib, k_Si28_He4_to_p_P31_derived,
        He4, P31, -1, H1, S34, -1, P31, S34, k_P31_He4_to_p_S34_reaclib, k_S34_p_to_He4_P31_derived,
        N, P32, -1, H1, Si32, -1, P32, Si32, k_P32_n_to_p_Si32_reaclib, k_Si32_p_to_n_P32_derived,
        H1, P32, -1, N, S32, -1, P32, S32, k_P32_p_to_n_S32_reaclib, k_S32_n_to_p_P32_derived,
        H1, P32, -1, He4, Si29, -1, P32, Si29, k_P32_p_to_He4_Si29_reaclib, k_Si29_He4_to_p_P32_derived,
        H1, P33, -1, He4, Si30, -1, P33, Si30, k_P33_p_to_He4_Si30_reaclib, k_Si30_He4_to_p_P33_derived,
        N, S32, -1, He4, Si29, -1, S32, Si29, k_S32_n_to_He4_Si29_reaclib, k_Si29_He4_to_n_S32_derived,
        N, S33, -1, H1, P33, -1, S33, P33, k_S33_n_to_p_P33_reaclib, k_P33_p_to_n_S33_derived,
        N, S33, -1, He4, Si30, -1, S33, Si30, k_S33_n_to_He4_Si30_reaclib, k_Si30_He4_to_n_S33_derived,
        N, S35, -1, He4, Si32, -1, S35, Si32, k_S35_n_to_He4_Si32_reaclib, k_Si32_He4_to_n_S35_derived,
        H1, S35, -1, He4, P32, -1, S35, P32, k_S35_p_to_He4_P32_reaclib, k_P32_He4_to_p_S35_derived,
        He4, S35, -1, N, Ar38, -1, S35, Ar38, k_S35_He4_to_n_Ar38_reaclib, k_Ar38_n_to_He4_S35_derived,
        H1, S36, -1, He4, P33, -1, S36, P33, k_S36_p_to_He4_P33_reaclib, k_P33_He4_to_p_S36_derived,
        N, Cl33, -1, H1, S33, -1, Cl33, S33, k_Cl33_n_to_p_S33_reaclib, k_S33_p_to_n_Cl33_derived,
        N, Cl33, -1, He4, P30, -1, Cl33, P30, k_Cl33_n_to_He4_P30_reaclib, k_P30_He4_to_n_Cl33_derived,
        He4, Cl33, -1, H1, Ar36, -1, Cl33, Ar36, k_Cl33_He4_to_p_Ar36_reaclib, k_Ar36_p_to_He4_Cl33_derived,
        N, Cl34, -1, H1, S34, -1, Cl34, S34, k_Cl34_n_to_p_S34_reaclib, k_S34_p_to_n_Cl34_derived,
        N, Cl34, -1, He4, P31, -1, Cl34, P31, k_Cl34_n_to_He4_P31_reaclib, k_P31_He4_to_n_Cl34_derived,
        He4, Cl34, -1, H1, Ar37, -1, Cl34, Ar37, k_Cl34_He4_to_p_Ar37_reaclib, k_Ar37_p_to_He4_Cl34_derived,
        N, Cl35, -1, H1, S35, -1, Cl35, S35, k_Cl35_n_to_p_S35_reaclib, k_S35_p_to_n_Cl35_derived,
        N, Cl35, -1, He4, P32, -1, Cl35, P32, k_Cl35_n_to_He4_P32_reaclib, k_P32_He4_to_n_Cl35_derived,
        H1, Cl35, -1, He4, S32, -1, Cl35, S32, k_Cl35_p_to_He4_S32_reaclib, k_S32_He4_to_p_Cl35_derived,
        He4, Cl35, -1, H1, Ar38, -1, Cl35, Ar38, k_Cl35_He4_to_p_Ar38_reaclib, k_Ar38_p_to_He4_Cl35_derived,
        N, Cl36, -1, H1, S36, -1, Cl36, S36, k_Cl36_n_to_p_S36_reaclib, k_S36_p_to_n_Cl36_derived,
        N, Cl36, -1, He4, P33, -1, Cl36, P33, k_Cl36_n_to_He4_P33_reaclib, k_P33_He4_to_n_Cl36_derived,
        H1, Cl36, -1, He4, S33, -1, Cl36, S33, k_Cl36_p_to_He4_S33_reaclib, k_S33_He4_to_p_Cl36_derived,
        N, Ar36, -1, H1, Cl36, -1, Ar36, Cl36, k_Ar36_n_to_p_Cl36_reaclib, k_Cl36_p_to_n_Ar36_derived,
        N, Ar36, -1, He4, S33, -1, Ar36, S33, k_Ar36_n_to_He4_S33_reaclib, k_S33_He4_to_n_Ar36_derived,
        N, Ar37, -1, H1, Cl37, -1, Ar37, Cl37, k_Ar37_n_to_p_Cl37_reaclib, k_Cl37_p_to_n_Ar37_derived,
        N, Ar39, -1, He4, S36, -1, Ar39, S36, k_Ar39_n_to_He4_S36_reaclib, k_S36_He4_to_n_Ar39_derived,
        H1, Ar39, -1, He4, Cl36, -1, Ar39, Cl36, k_Ar39_p_to_He4_Cl36_reaclib, k_Cl36_He4_to_p_Ar39_derived,
        H1, Ar40, -1, He4, Cl37, -1, Ar40, Cl37, k_Ar40_p_to_He4_Cl37_reaclib, k_Cl37_He4_to_p_Ar40_derived,
        N, K37, -1, H1, Ar37, -1, K37, Ar37, k_K37_n_to_p_Ar37_reaclib, k_Ar37_p_to_n_K37_derived,
        N, K37, -1, He4, Cl34, -1, K37, Cl34, k_K37_n_to_He4_Cl34_reaclib, k_Cl34_He4_to_n_K37_derived,
        He4, K37, -1, H1, Ca40, -1, K37, Ca40, k_K37_He4_to_p_Ca40_reaclib, k_Ca40_p_to_He4_K37_derived,
        N, K38, -1, H1, Ar38, -1, K38, Ar38, k_K38_n_to_p_Ar38_reaclib, k_Ar38_p_to_n_K38_derived,
        N, K38, -1, He4, Cl35, -1, K38, Cl35, k_K38_n_to_He4_Cl35_reaclib, k_Cl35_He4_to_n_K38_derived,
        He4, K38, -1, H1, Ca41, -1, K38, Ca41, k_K38_He4_to_p_Ca41_reaclib, k_Ca41_p_to_He4_K38_derived,
        N, K39, -1, H1, Ar39, -1, K39, Ar39, k_K39_n_to_p_Ar39_reaclib, k_Ar39_p_to_n_K39_derived,
        N, K39, -1, He4, Cl36, -1, K39, Cl36, k_K39_n_to_He4_Cl36_reaclib, k_Cl36_He4_to_n_K39_derived,
        H1, K39, -1, He4, Ar36, -1, K39, Ar36, k_K39_p_to_He4_Ar36_reaclib, k_Ar36_He4_to_p_K39_derived,
        N, K40, -1, H1, Ar40, -1, K40, Ar40, k_K40_n_to_p_Ar40_reaclib, k_Ar40_p_to_n_K40_derived,
        N, K40, -1, He4, Cl37, -1, K40, Cl37, k_K40_n_to_He4_Cl37_reaclib, k_Cl37_He4_to_n_K40_derived,
        H1, K40, -1, N, Ca40, -1, K40, Ca40, k_K40_p_to_n_Ca40_reaclib, k_Ca40_n_to_p_K40_derived,
        H1, K40, -1, He4, Ar37, -1, K40, Ar37, k_K40_p_to_He4_Ar37_reaclib, k_Ar37_He4_to_p_K40_derived,
        He4, K40, -1, H1, Ca43, -1, K40, Ca43, k_K40_He4_to_p_Ca43_reaclib, k_Ca43_p_to_He4_K40_derived,
        H1, K41, -1, He4, Ar38, -1, K41, Ar38, k_K41_p_to_He4_Ar38_reaclib, k_Ar38_He4_to_p_K41_derived,
        He4, K41, -1, H1, Ca44, -1, K41, Ca44, k_K41_He4_to_p_Ca44_reaclib, k_Ca44_p_to_He4_K41_derived,
        N, Ca40, -1, He4, Ar37, -1, Ca40, Ar37, k_Ca40_n_to_He4_Ar37_reaclib, k_Ar37_He4_to_n_Ca40_derived,
        N, Ca41, -1, He4, Ar38, -1, Ca41, Ar38, k_Ca41_n_to_He4_Ar38_reaclib, k_Ar38_He4_to_n_Ca41_derived,
        N, Ca42, -1, He4, Ar39, -1, Ca42, Ar39, k_Ca42_n_to_He4_Ar39_reaclib, k_Ar39_He4_to_n_Ca42_derived,
        H1, Ca42, -1, He4, K39, -1, Ca42, K39, k_Ca42_p_to_He4_K39_reaclib, k_K39_He4_to_p_Ca42_derived,
        N, Ca43, -1, He4, Ar40, -1, Ca43, Ar40, k_Ca43_n_to_He4_Ar40_reaclib, k_Ar40_He4_to_n_Ca43_derived,
        He4, Ca43, -1, N, Ti46, -1, Ca43, Ti46, k_Ca43_He4_to_n_Ti46_reaclib, k_Ti46_n_to_He4_Ca43_derived,
        He4, Ca45, -1, N, Ti48, -1, Ca45, Ti48, k_Ca45_He4_to_n_Ti48_reaclib, k_Ti48_n_to_He4_Ca45_derived,
        H1, Ca47, -1, N, Sc47, -1, Ca47, Sc47, k_Ca47_p_to_n_Sc47_reaclib, k_Sc47_n_to_p_Ca47_derived,
        He4, Ca47, -1, N, Ti50, -1, Ca47, Ti50, k_Ca47_He4_to_n_Ti50_reaclib, k_Ti50_n_to_He4_Ca47_derived,
        N, Sc43, -1, H1, Ca43, -1, Sc43, Ca43, k_Sc43_n_to_p_Ca43_reaclib, k_Ca43_p_to_n_Sc43_derived,
        N, Sc43, -1, He4, K40, -1, Sc43, K40, k_Sc43_n_to_He4_K40_reaclib, k_K40_He4_to_n_Sc43_derived,
        H1, Sc43, -1, He4, Ca40, -1, Sc43, Ca40, k_Sc43_p_to_He4_Ca40_reaclib, k_Ca40_He4_to_p_Sc43_derived,
        He4, Sc43, -1, H1, Ti46, -1, Sc43, Ti46, k_Sc43_He4_to_p_Ti46_reaclib, k_Ti46_p_to_He4_Sc43_derived,
        N, Sc44, -1, H1, Ca44, -1, Sc44, Ca44, k_Sc44_n_to_p_Ca44_reaclib, k_Ca44_p_to_n_Sc44_derived,
        H1, Sc44, -1, He4, Ca41, -1, Sc44, Ca41, k_Sc44_p_to_He4_Ca41_reaclib, k_Ca41_He4_to_p_Sc44_derived,
        He4, Sc44, -1, H1, Ti47, -1, Sc44, Ti47, k_Sc44_He4_to_p_Ti47_reaclib, k_Ti47_p_to_He4_Sc44_derived,
        N, Sc45, -1, H1, Ca45, -1, Sc45, Ca45, k_Sc45_n_to_p_Ca45_reaclib, k_Ca45_p_to_n_Sc45_derived,
        H1, Sc45, -1, He4, Ca42, -1, Sc45, Ca42, k_Sc45_p_to_He4_Ca42_reaclib, k_Ca42_He4_to_p_Sc45_derived,
        He4, Sc45, -1, H1, Ti48, -1, Sc45, Ti48, k_Sc45_He4_to_p_Ti48_reaclib, k_Ti48_p_to_He4_Sc45_derived,
        N, Sc46, -1, H1, Ca46, -1, Sc46, Ca46, k_Sc46_n_to_p_Ca46_reaclib, k_Ca46_p_to_n_Sc46_derived,
        H1, Sc46, -1, N, Ti46, -1, Sc46, Ti46, k_Sc46_p_to_n_Ti46_reaclib, k_Ti46_n_to_p_Sc46_derived,
        H1, Sc46, -1, He4, Ca43, -1, Sc46, Ca43, k_Sc46_p_to_He4_Ca43_reaclib, k_Ca43_He4_to_p_Sc46_derived,
        He4, Sc46, -1, N, V49, -1, Sc46, V49, k_Sc46_He4_to_n_V49_reaclib, k_V49_n_to_He4_Sc46_derived,
        He4, Sc46, -1, H1, Ti49, -1, Sc46, Ti49, k_Sc46_He4_to_p_Ti49_reaclib, k_Ti49_p_to_He4_Sc46_derived,
        H1, Sc47, -1, He4, Ca44, -1, Sc47, Ca44, k_Sc47_p_to_He4_Ca44_reaclib, k_Ca44_He4_to_p_Sc47_derived,
        He4, Sc47, -1, H1, Ti50, -1, Sc47, Ti50, k_Sc47_He4_to_p_Ti50_reaclib, k_Ti50_p_to_He4_Sc47_derived,
        N, Sc48, -1, H1, Ca48, -1, Sc48, Ca48, k_Sc48_n_to_p_Ca48_reaclib, k_Ca48_p_to_n_Sc48_derived,
        H1, Sc48, -1, N, Ti48, -1, Sc48, Ti48, k_Sc48_p_to_n_Ti48_reaclib, k_Ti48_n_to_p_Sc48_derived,
        H1, Sc48, -1, He4, Ca45, -1, Sc48, Ca45, k_Sc48_p_to_He4_Ca45_reaclib, k_Ca45_He4_to_p_Sc48_derived,
        He4, Sc48, -1, N, V51, -1, Sc48, V51, k_Sc48_He4_to_n_V51_reaclib, k_V51_n_to_He4_Sc48_derived,
        He4, Sc48, -1, H1, Ti51, -1, Sc48, Ti51, k_Sc48_He4_to_p_Ti51_reaclib, k_Ti51_p_to_He4_Sc48_derived,
        H1, Sc49, -1, N, Ti49, -1, Sc49, Ti49, k_Sc49_p_to_n_Ti49_reaclib, k_Ti49_n_to_p_Sc49_derived,
        H1, Sc49, -1, He4, Ca46, -1, Sc49, Ca46, k_Sc49_p_to_He4_Ca46_reaclib, k_Ca46_He4_to_p_Sc49_derived,
        N, Ti44, -1, H1, Sc44, -1, Ti44, Sc44, k_Ti44_n_to_p_Sc44_reaclib, k_Sc44_p_to_n_Ti44_derived,
        N, Ti44, -1, He4, Ca41, -1, Ti44, Ca41, k_Ti44_n_to_He4_Ca41_reaclib, k_Ca41_He4_to_n_Ti44_derived,
        N, Ti45, -1, H1, Sc45, -1, Ti45, Sc45, k_Ti45_n_to_p_Sc45_reaclib, k_Sc45_p_to_n_Ti45_derived,
        N, Ti45, -1, He4, Ca42, -1, Ti45, Ca42, k_Ti45_n_to_He4_Ca42_reaclib, k_Ca42_He4_to_n_Ti45_derived,
        He4, Ti45, -1, H1, V48, -1, Ti45, V48, k_Ti45_He4_to_p_V48_reaclib, k_V48_p_to_He4_Ti45_derived,
        N, Ti47, -1, H1, Sc47, -1, Ti47, Sc47, k_Ti47_n_to_p_Sc47_reaclib, k_Sc47_p_to_n_Ti47_derived,
        N, Ti47, -1, He4, Ca44, -1, Ti47, Ca44, k_Ti47_n_to_He4_Ca44_reaclib, k_Ca44_He4_to_n_Ti47_derived,
        N, Ti49, -1, He4, Ca46, -1, Ti49, Ca46, k_Ti49_n_to_He4_Ca46_reaclib, k_Ca46_He4_to_n_Ti49_derived,
        He4, Ti49, -1, N, Cr52, -1, Ti49, Cr52, k_Ti49_He4_to_n_Cr52_reaclib, k_Cr52_n_to_He4_Ti49_derived,
        N, Ti51, -1, He4, Ca48, -1, Ti51, Ca48, k_Ti51_n_to_He4_Ca48_reaclib, k_Ca48_He4_to_n_Ti51_derived,
        H1, Ti51, -1, N, V51, -1, Ti51, V51, k_Ti51_p_to_n_V51_reaclib, k_V51_n_to_p_Ti51_derived,
        He4, Ti51, -1, N, Cr54, -1, Ti51, Cr54, k_Ti51_He4_to_n_Cr54_reaclib, k_Cr54_n_to_He4_Ti51_derived,
        N, V46, -1, H1, Ti46, -1, V46, Ti46, k_V46_n_to_p_Ti46_reaclib, k_Ti46_p_to_n_V46_derived,
        N, V46, -1, He4, Sc43, -1, V46, Sc43, k_V46_n_to_He4_Sc43_reaclib, k_Sc43_He4_to_n_V46_derived,
        He4, V46, -1, H1, Cr49, -1, V46, Cr49, k_V46_He4_to_p_Cr49_reaclib, k_Cr49_p_to_He4_V46_derived,
        N, V47, -1, H1, Ti47, -1, V47, Ti47, k_V47_n_to_p_Ti47_reaclib, k_Ti47_p_to_n_V47_derived,
        N, V47, -1, He4, Sc44, -1, V47, Sc44, k_V47_n_to_He4_Sc44_reaclib, k_Sc44_He4_to_n_V47_derived,
        He4, V47, -1, H1, Cr50, -1, V47, Cr50, k_V47_He4_to_p_Cr50_reaclib, k_Cr50_p_to_He4_V47_derived,
        N, V48, -1, H1, Ti48, -1, V48, Ti48, k_V48_n_to_p_Ti48_reaclib, k_Ti48_p_to_n_V48_derived,
        N, V48, -1, He4, Sc45, -1, V48, Sc45, k_V48_n_to_He4_Sc45_reaclib, k_Sc45_He4_to_n_V48_derived,
        He4, V48, -1, H1, Cr51, -1, V48, Cr51, k_V48_He4_to_p_Cr51_reaclib, k_Cr51_p_to_He4_V48_derived,
        N, V49, -1, H1, Ti49, -1, V49, Ti49, k_V49_n_to_p_Ti49_reaclib, k_Ti49_p_to_n_V49_derived,
        H1, V49, -1, He4, Ti46, -1, V49, Ti46, k_V49_p_to_He4_Ti46_reaclib, k_Ti46_He4_to_p_V49_derived,
        He4, V49, -1, H1, Cr52, -1, V49, Cr52, k_V49_He4_to_p_Cr52_reaclib, k_Cr52_p_to_He4_V49_derived,
        N, V50, -1, H1, Ti50, -1, V50, Ti50, k_V50_n_to_p_Ti50_reaclib, k_Ti50_p_to_n_V50_derived,
        N, V50, -1, He4, Sc47, -1, V50, Sc47, k_V50_n_to_He4_Sc47_reaclib, k_Sc47_He4_to_n_V50_derived,
        H1, V50, -1, N, Cr50, -1, V50, Cr50, k_V50_p_to_n_Cr50_reaclib, k_Cr50_n_to_p_V50_derived,
        H1, V50, -1, He4, Ti47, -1, V50, Ti47, k_V50_p_to_He4_Ti47_reaclib, k_Ti47_He4_to_p_V50_derived,
        He4, V50, -1, H1, Cr53, -1, V50, Cr53, k_V50_He4_to_p_Cr53_reaclib, k_Cr53_p_to_He4_V50_derived,
        H1, V51, -1, He4, Ti48, -1, V51, Ti48, k_V51_p_to_He4_Ti48_reaclib, k_Ti48_He4_to_p_V51_derived,
        N, V52, -1, He4, Sc49, -1, V52, Sc49, k_V52_n_to_He4_Sc49_reaclib, k_Sc49_He4_to_n_V52_derived,
        H1, V52, -1, N, Cr52, -1, V52, Cr52, k_V52_p_to_n_Cr52_reaclib, k_Cr52_n_to_p_V52_derived,
        H1, V52, -1, He4, Ti49, -1, V52, Ti49, k_V52_p_to_He4_Ti49_reaclib, k_Ti49_He4_to_p_V52_derived,
        He4, V52, -1, N, Mn55, -1, V52, Mn55, k_V52_He4_to_n_Mn55_reaclib, k_Mn55_n_to_He4_V52_derived,
        N, Cr48, -1, H1, V48, -1, Cr48, V48, k_Cr48_n_to_p_V48_reaclib, k_V48_p_to_n_Cr48_derived,
        N, Cr48, -1, He4, Ti45, -1, Cr48, Ti45, k_Cr48_n_to_He4_Ti45_reaclib, k_Ti45_He4_to_n_Cr48_derived,
        He4, Cr48, -1, H1, Mn51, -1, Cr48, Mn51, k_Cr48_He4_to_p_Mn51_reaclib, k_Mn51_p_to_He4_Cr48_derived,
        N, Cr49, -1, H1, V49, -1, Cr49, V49, k_Cr49_n_to_p_V49_reaclib, k_V49_p_to_n_Cr49_derived,
        N, Cr49, -1, He4, Ti46, -1, Cr49, Ti46, k_Cr49_n_to_He4_Ti46_reaclib, k_Ti46_He4_to_n_Cr49_derived,
        He4, Cr49, -1, H1, Mn52, -1, Cr49, Mn52, k_Cr49_He4_to_p_Mn52_reaclib, k_Mn52_p_to_He4_Cr49_derived,
        N, Cr50, -1, He4, Ti47, -1, Cr50, Ti47, k_Cr50_n_to_He4_Ti47_reaclib, k_Ti47_He4_to_n_Cr50_derived,
        N, Cr51, -1, H1, V51, -1, Cr51, V51, k_Cr51_n_to_p_V51_reaclib, k_V51_p_to_n_Cr51_derived,
        N, Cr51, -1, He4, Ti48, -1, Cr51, Ti48, k_Cr51_n_to_He4_Ti48_reaclib, k_Ti48_He4_to_n_Cr51_derived,
        N, Cr53, -1, He4, Ti50, -1, Cr53, Ti50, k_Cr53_n_to_He4_Ti50_reaclib, k_Ti50_He4_to_n_Cr53_derived,
        H1, Cr54, -1, He4, V51, -1, Cr54, V51, k_Cr54_p_to_He4_V51_reaclib, k_V51_He4_to_p_Cr54_derived,
        N, Mn50, -1, H1, Cr50, -1, Mn50, Cr50, k_Mn50_n_to_p_Cr50_reaclib, k_Cr50_p_to_n_Mn50_derived,
        N, Mn50, -1, He4, V47, -1, Mn50, V47, k_Mn50_n_to_He4_V47_reaclib, k_V47_He4_to_n_Mn50_derived,
        He4, Mn50, -1, H1, Fe53, -1, Mn50, Fe53, k_Mn50_He4_to_p_Fe53_reaclib, k_Fe53_p_to_He4_Mn50_derived,
        N, Mn51, -1, H1, Cr51, -1, Mn51, Cr51, k_Mn51_n_to_p_Cr51_reaclib, k_Cr51_p_to_n_Mn51_derived,
        N, Mn51, -1, He4, V48, -1, Mn51, V48, k_Mn51_n_to_He4_V48_reaclib, k_V48_He4_to_n_Mn51_derived,
        He4, Mn51, -1, H1, Fe54, -1, Mn51, Fe54, k_Mn51_He4_to_p_Fe54_reaclib, k_Fe54_p_to_He4_Mn51_derived,
        N, Mn52, -1, H1, Cr52, -1, Mn52, Cr52, k_Mn52_n_to_p_Cr52_reaclib, k_Cr52_p_to_n_Mn52_derived,
        N, Mn52, -1, He4, V49, -1, Mn52, V49, k_Mn52_n_to_He4_V49_reaclib, k_V49_He4_to_n_Mn52_derived,
        He4, Mn52, -1, H1, Fe55, -1, Mn52, Fe55, k_Mn52_He4_to_p_Fe55_reaclib, k_Fe55_p_to_He4_Mn52_derived,
        N, Mn53, -1, H1, Cr53, -1, Mn53, Cr53, k_Mn53_n_to_p_Cr53_reaclib, k_Cr53_p_to_n_Mn53_derived,
        N, Mn53, -1, He4, V50, -1, Mn53, V50, k_Mn53_n_to_He4_V50_reaclib, k_V50_He4_to_n_Mn53_derived,
        H1, Mn53, -1, He4, Cr50, -1, Mn53, Cr50, k_Mn53_p_to_He4_Cr50_reaclib, k_Cr50_He4_to_p_Mn53_derived,
        He4, Mn53, -1, H1, Fe56, -1, Mn53, Fe56, k_Mn53_He4_to_p_Fe56_reaclib, k_Fe56_p_to_He4_Mn53_derived,
        N, Mn54, -1, H1, Cr54, -1, Mn54, Cr54, k_Mn54_n_to_p_Cr54_reaclib, k_Cr54_p_to_n_Mn54_derived,
        N, Mn54, -1, He4, V51, -1, Mn54, V51, k_Mn54_n_to_He4_V51_reaclib, k_V51_He4_to_n_Mn54_derived,
        H1, Mn54, -1, He4, Cr51, -1, Mn54, Cr51, k_Mn54_p_to_He4_Cr51_reaclib, k_Cr51_He4_to_p_Mn54_derived,
        H1, Mn55, -1, He4, Cr52, -1, Mn55, Cr52, k_Mn55_p_to_He4_Cr52_reaclib, k_Cr52_He4_to_p_Mn55_derived,
        N, Fe52, -1, H1, Mn52, -1, Fe52, Mn52, k_Fe52_n_to_p_Mn52_reaclib, k_Mn52_p_to_n_Fe52_derived,
        N, Fe52, -1, He4, Cr49, -1, Fe52, Cr49, k_Fe52_n_to_He4_Cr49_reaclib, k_Cr49_He4_to_n_Fe52_derived,
        He4, Fe52, -1, H1, Co55, -1, Fe52, Co55, k_Fe52_He4_to_p_Co55_reaclib, k_Co55_p_to_He4_Fe52_derived,
        N, Fe53, -1, H1, Mn53, -1, Fe53, Mn53, k_Fe53_n_to_p_Mn53_reaclib, k_Mn53_p_to_n_Fe53_derived,
        N, Fe53, -1, He4, Cr50, -1, Fe53, Cr50, k_Fe53_n_to_He4_Cr50_reaclib, k_Cr50_He4_to_n_Fe53_derived,
        He4, Fe53, -1, H1, Co56, -1, Fe53, Co56, k_Fe53_He4_to_p_Co56_reaclib, k_Co56_p_to_He4_Fe53_derived,
        N, Fe54, -1, H1, Mn54, -1, Fe54, Mn54, k_Fe54_n_to_p_Mn54_reaclib, k_Mn54_p_to_n_Fe54_derived,
        N, Fe54, -1, He4, Cr51, -1, Fe54, Cr51, k_Fe54_n_to_He4_Cr51_reaclib, k_Cr51_He4_to_n_Fe54_derived,
        N, Fe55, -1, H1, Mn55, -1, Fe55, Mn55, k_Fe55_n_to_p_Mn55_reaclib, k_Mn55_p_to_n_Fe55_derived,
        N, Fe55, -1, He4, Cr52, -1, Fe55, Cr52, k_Fe55_n_to_He4_Cr52_reaclib, k_Cr52_He4_to_n_Fe55_derived,
        N, Fe56, -1, He4, Cr53, -1, Fe56, Cr53, k_Fe56_n_to_He4_Cr53_reaclib, k_Cr53_He4_to_n_Fe56_derived,
        N, Fe57, -1, He4, Cr54, -1, Fe57, Cr54, k_Fe57_n_to_He4_Cr54_reaclib, k_Cr54_He4_to_n_Fe57_derived,
        H1, Fe57, -1, He4, Mn54, -1, Fe57, Mn54, k_Fe57_p_to_He4_Mn54_reaclib, k_Mn54_He4_to_p_Fe57_derived,
        H1, Fe58, -1, He4, Mn55, -1, Fe58, Mn55, k_Fe58_p_to_He4_Mn55_reaclib, k_Mn55_He4_to_p_Fe58_derived,
        N, Co53, -1, H1, Fe53, -1, Co53, Fe53, k_Co53_n_to_p_Fe53_reaclib, k_Fe53_p_to_n_Co53_derived,
        N, Co53, -1, He4, Mn50, -1, Co53, Mn50, k_Co53_n_to_He4_Mn50_reaclib, k_Mn50_He4_to_n_Co53_derived,
        He4, Co53, -1, H1, Ni56, -1, Co53, Ni56, k_Co53_He4_to_p_Ni56_reaclib, k_Ni56_p_to_He4_Co53_derived,
        N, Co54, -1, H1, Fe54, -1, Co54, Fe54, k_Co54_n_to_p_Fe54_reaclib, k_Fe54_p_to_n_Co54_derived,
        N, Co54, -1, He4, Mn51, -1, Co54, Mn51, k_Co54_n_to_He4_Mn51_reaclib, k_Mn51_He4_to_n_Co54_derived,
        He4, Co54, -1, H1, Ni57, -1, Co54, Ni57, k_Co54_He4_to_p_Ni57_reaclib, k_Ni57_p_to_He4_Co54_derived,
        N, Co55, -1, H1, Fe55, -1, Co55, Fe55, k_Co55_n_to_p_Fe55_reaclib, k_Fe55_p_to_n_Co55_derived,
        N, Co55, -1, He4, Mn52, -1, Co55, Mn52, k_Co55_n_to_He4_Mn52_reaclib, k_Mn52_He4_to_n_Co55_derived,
        He4, Co55, -1, H1, Ni58, -1, Co55, Ni58, k_Co55_He4_to_p_Ni58_reaclib, k_Ni58_p_to_He4_Co55_derived,
        N, Co56, -1, H1, Fe56, -1, Co56, Fe56, k_Co56_n_to_p_Fe56_reaclib, k_Fe56_p_to_n_Co56_derived,
        N, Co56, -1, He4, Mn53, -1, Co56, Mn53, k_Co56_n_to_He4_Mn53_reaclib, k_Mn53_He4_to_n_Co56_derived,
        He4, Co56, -1, H1, Ni59, -1, Co56, Ni59, k_Co56_He4_to_p_Ni59_reaclib, k_Ni59_p_to_He4_Co56_derived,
        N, Co57, -1, H1, Fe57, -1, Co57, Fe57, k_Co57_n_to_p_Fe57_reaclib, k_Fe57_p_to_n_Co57_derived,
        N, Co57, -1, He4, Mn54, -1, Co57, Mn54, k_Co57_n_to_He4_Mn54_reaclib, k_Mn54_He4_to_n_Co57_derived,
        H1, Co57, -1, He4, Fe54, -1, Co57, Fe54, k_Co57_p_to_He4_Fe54_reaclib, k_Fe54_He4_to_p_Co57_derived,
        He4, Co57, -1, H1, Ni60, -1, Co57, Ni60, k_Co57_He4_to_p_Ni60_reaclib, k_Ni60_p_to_He4_Co57_derived,
        N, Co58, -1, H1, Fe58, -1, Co58, Fe58, k_Co58_n_to_p_Fe58_reaclib, k_Fe58_p_to_n_Co58_derived,
        N, Co58, -1, He4, Mn55, -1, Co58, Mn55, k_Co58_n_to_He4_Mn55_reaclib, k_Mn55_He4_to_n_Co58_derived,
        H1, Co58, -1, He4, Fe55, -1, Co58, Fe55, k_Co58_p_to_He4_Fe55_reaclib, k_Fe55_He4_to_p_Co58_derived,
        H1, Co59, -1, He4, Fe56, -1, Co59, Fe56, k_Co59_p_to_He4_Fe56_reaclib, k_Fe56_He4_to_p_Co59_derived,
        N, Ni56, -1, H1, Co56, -1, Ni56, Co56, k_Ni56_n_to_p_Co56_reaclib, k_Co56_p_to_n_Ni56_derived,
        N, Ni56, -1, He4, Fe53, -1, Ni56, Fe53, k_Ni56_n_to_He4_Fe53_reaclib, k_Fe53_He4_to_n_Ni56_derived,
        N, Ni57, -1, H1, Co57, -1, Ni57, Co57, k_Ni57_n_to_p_Co57_reaclib, k_Co57_p_to_n_Ni57_derived,
        N, Ni57, -1, He4, Fe54, -1, Ni57, Fe54, k_Ni57_n_to_He4_Fe54_reaclib, k_Fe54_He4_to_n_Ni57_derived,
        N, Ni58, -1, H1, Co58, -1, Ni58, Co58, k_Ni58_n_to_p_Co58_reaclib, k_Co58_p_to_n_Ni58_derived,
        N, Ni58, -1, He4, Fe55, -1, Ni58, Fe55, k_Ni58_n_to_He4_Fe55_reaclib, k_Fe55_He4_to_n_Ni58_derived,
        N, Ni59, -1, H1, Co59, -1, Ni59, Co59, k_Ni59_n_to_p_Co59_reaclib, k_Co59_p_to_n_Ni59_derived,
        N, Ni59, -1, He4, Fe56, -1, Ni59, Fe56, k_Ni59_n_to_He4_Fe56_reaclib, k_Fe56_He4_to_n_Ni59_derived,
        N, Ni60, -1, He4, Fe57, -1, Ni60, Fe57, k_Ni60_n_to_He4_Fe57_reaclib, k_Fe57_He4_to_n_Ni60_derived,
        N, Ni61, -1, He4, Fe58, -1, Ni61, Fe58, k_Ni61_n_to_He4_Fe58_reaclib, k_Fe58_He4_to_n_Ni61_derived,
        H1, Ni61, -1, He4, Co58, -1, Ni61, Co58, k_Ni61_p_to_He4_Co58_reaclib, k_Co58_He4_to_p_Ni61_derived,
        H1, Ni62, -1, He4, Co59, -1, Ni62, Co59, k_Ni62_p_to_He4_Co59_reaclib, k_Co59_He4_to_p_Ni62_derived,
        N, Cu57, -1, H1, Ni57, -1, Cu57, Ni57, k_Cu57_n_to_p_Ni57_reaclib, k_Ni57_p_to_n_Cu57_derived,
        N, Cu57, -1, He4, Co54, -1, Cu57, Co54, k_Cu57_n_to_He4_Co54_reaclib, k_Co54_He4_to_n_Cu57_derived,
        He4, Cu57, -1, H1, Zn60, -1, Cu57, Zn60, k_Cu57_He4_to_p_Zn60_reaclib, k_Zn60_p_to_He4_Cu57_derived,
        N, Cu58, -1, H1, Ni58, -1, Cu58, Ni58, k_Cu58_n_to_p_Ni58_reaclib, k_Ni58_p_to_n_Cu58_derived,
        N, Cu58, -1, He4, Co55, -1, Cu58, Co55, k_Cu58_n_to_He4_Co55_reaclib, k_Co55_He4_to_n_Cu58_derived,
        N, Cu59, -1, H1, Ni59, -1, Cu59, Ni59, k_Cu59_n_to_p_Ni59_reaclib, k_Ni59_p_to_n_Cu59_derived,
        N, Cu59, -1, He4, Co56, -1, Cu59, Co56, k_Cu59_n_to_He4_Co56_reaclib, k_Co56_He4_to_n_Cu59_derived,
        H1, Cu59, -1, He4, Ni56, -1, Cu59, Ni56, k_Cu59_p_to_He4_Ni56_reaclib, k_Ni56_He4_to_p_Cu59_derived,
        N, Cu60, -1, H1, Ni60, -1, Cu60, Ni60, k_Cu60_n_to_p_Ni60_reaclib, k_Ni60_p_to_n_Cu60_derived,
        N, Cu60, -1, He4, Co57, -1, Cu60, Co57, k_Cu60_n_to_He4_Co57_reaclib, k_Co57_He4_to_n_Cu60_derived,
        H1, Cu60, -1, He4, Ni57, -1, Cu60, Ni57, k_Cu60_p_to_He4_Ni57_reaclib, k_Ni57_He4_to_p_Cu60_derived,
        N, Cu61, -1, H1, Ni61, -1, Cu61, Ni61, k_Cu61_n_to_p_Ni61_reaclib, k_Ni61_p_to_n_Cu61_derived,
        N, Cu61, -1, He4, Co58, -1, Cu61, Co58, k_Cu61_n_to_He4_Co58_reaclib, k_Co58_He4_to_n_Cu61_derived,
        H1, Cu61, -1, He4, Ni58, -1, Cu61, Ni58, k_Cu61_p_to_He4_Ni58_reaclib, k_Ni58_He4_to_p_Cu61_derived,
        N, Cu62, -1, H1, Ni62, -1, Cu62, Ni62, k_Cu62_n_to_p_Ni62_reaclib, k_Ni62_p_to_n_Cu62_derived,
        N, Cu62, -1, He4, Co59, -1, Cu62, Co59, k_Cu62_n_to_He4_Co59_reaclib, k_Co59_He4_to_n_Cu62_derived,
        H1, Cu62, -1, He4, Ni59, -1, Cu62, Ni59, k_Cu62_p_to_He4_Ni59_reaclib, k_Ni59_He4_to_p_Cu62_derived,
        N, Cu63, -1, H1, Ni63, -1, Cu63, Ni63, k_Cu63_n_to_p_Ni63_reaclib, k_Ni63_p_to_n_Cu63_derived,
        H1, Cu63, -1, He4, Ni60, -1, Cu63, Ni60, k_Cu63_p_to_He4_Ni60_reaclib, k_Ni60_He4_to_p_Cu63_derived,
        N, Cu64, -1, H1, Ni64, -1, Cu64, Ni64, k_Cu64_n_to_p_Ni64_reaclib, k_Ni64_p_to_n_Cu64_derived,
        H1, Cu64, -1, He4, Ni61, -1, Cu64, Ni61, k_Cu64_p_to_He4_Ni61_reaclib, k_Ni61_He4_to_p_Cu64_derived,
        H1, Cu65, -1, He4, Ni62, -1, Cu65, Ni62, k_Cu65_p_to_He4_Ni62_reaclib, k_Ni62_He4_to_p_Cu65_derived,
        N, Zn59, -1, H1, Cu59, -1, Zn59, Cu59, k_Zn59_n_to_p_Cu59_reaclib, k_Cu59_p_to_n_Zn59_derived,
        N, Zn59, -1, He4, Ni56, -1, Zn59, Ni56, k_Zn59_n_to_He4_Ni56_reaclib, k_Ni56_He4_to_n_Zn59_derived,
        N, Zn60, -1, H1, Cu60, -1, Zn60, Cu60, k_Zn60_n_to_p_Cu60_reaclib, k_Cu60_p_to_n_Zn60_derived,
        N, Zn60, -1, He4, Ni57, -1, Zn60, Ni57, k_Zn60_n_to_He4_Ni57_reaclib, k_Ni57_He4_to_n_Zn60_derived,
        N, Zn61, -1, H1, Cu61, -1, Zn61, Cu61, k_Zn61_n_to_p_Cu61_reaclib, k_Cu61_p_to_n_Zn61_derived,
        N, Zn61, -1, He4, Ni58, -1, Zn61, Ni58, k_Zn61_n_to_He4_Ni58_reaclib, k_Ni58_He4_to_n_Zn61_derived,
        H1, Zn61, -1, He4, Cu58, -1, Zn61, Cu58, k_Zn61_p_to_He4_Cu58_reaclib, k_Cu58_He4_to_p_Zn61_derived,
        N, Zn62, -1, H1, Cu62, -1, Zn62, Cu62, k_Zn62_n_to_p_Cu62_reaclib, k_Cu62_p_to_n_Zn62_derived,
        N, Zn62, -1, He4, Ni59, -1, Zn62, Ni59, k_Zn62_n_to_He4_Ni59_reaclib, k_Ni59_He4_to_n_Zn62_derived,
        H1, Zn62, -1, He4, Cu59, -1, Zn62, Cu59, k_Zn62_p_to_He4_Cu59_reaclib, k_Cu59_He4_to_p_Zn62_derived,
        N, Zn63, -1, H1, Cu63, -1, Zn63, Cu63, k_Zn63_n_to_p_Cu63_reaclib, k_Cu63_p_to_n_Zn63_derived,
        N, Zn63, -1, He4, Ni60, -1, Zn63, Ni60, k_Zn63_n_to_He4_Ni60_reaclib, k_Ni60_He4_to_n_Zn63_derived,
        H1, Zn63, -1, He4, Cu60, -1, Zn63, Cu60, k_Zn63_p_to_He4_Cu60_reaclib, k_Cu60_He4_to_p_Zn63_derived,
        N, Zn64, -1, H1, Cu64, -1, Zn64, Cu64, k_Zn64_n_to_p_Cu64_reaclib, k_Cu64_p_to_n_Zn64_derived,
        N, Zn64, -1, He4, Ni61, -1, Zn64, Ni61, k_Zn64_n_to_He4_Ni61_reaclib, k_Ni61_He4_to_n_Zn64_derived,
        H1, Zn64, -1, He4, Cu61, -1, Zn64, Cu61, k_Zn64_p_to_He4_Cu61_reaclib, k_Cu61_He4_to_p_Zn64_derived,
        N, Zn65, -1, H1, Cu65, -1, Zn65, Cu65, k_Zn65_n_to_p_Cu65_reaclib, k_Cu65_p_to_n_Zn65_derived,
        N, Zn65, -1, He4, Ni62, -1, Zn65, Ni62, k_Zn65_n_to_He4_Ni62_reaclib, k_Ni62_He4_to_n_Zn65_derived,
        H1, Zn65, -1, He4, Cu62, -1, Zn65, Cu62, k_Zn65_p_to_He4_Cu62_reaclib, k_Cu62_He4_to_p_Zn65_derived,
        N, Zn66, -1, He4, Ni63, -1, Zn66, Ni63, k_Zn66_n_to_He4_Ni63_reaclib, k_Ni63_He4_to_n_Zn66_derived,
        H1, Zn66, -1, He4, Cu63, -1, Zn66, Cu63, k_Zn66_p_to_He4_Cu63_reaclib, k_Cu63_He4_to_p_Zn66_derived,
        N, Ga62, -1, H1, Zn62, -1, Ga62, Zn62, k_Ga62_n_to_p_Zn62_reaclib, k_Zn62_p_to_n_Ga62_derived,
        N, Ga62, -1, He4, Cu59, -1, Ga62, Cu59, k_Ga62_n_to_He4_Cu59_reaclib, k_Cu59_He4_to_n_Ga62_derived,
        H1, Ga62, -1, He4, Zn59, -1, Ga62, Zn59, k_Ga62_p_to_He4_Zn59_reaclib, k_Zn59_He4_to_p_Ga62_derived,
        N, Ga63, -1, H1, Zn63, -1, Ga63, Zn63, k_Ga63_n_to_p_Zn63_reaclib, k_Zn63_p_to_n_Ga63_derived,
        N, Ga63, -1, He4, Cu60, -1, Ga63, Cu60, k_Ga63_n_to_He4_Cu60_reaclib, k_Cu60_He4_to_n_Ga63_derived,
        H1, Ga63, -1, He4, Zn60, -1, Ga63, Zn60, k_Ga63_p_to_He4_Zn60_reaclib, k_Zn60_He4_to_p_Ga63_derived,
        N, Ga64, -1, H1, Zn64, -1, Ga64, Zn64, k_Ga64_n_to_p_Zn64_reaclib, k_Zn64_p_to_n_Ga64_derived,
        N, Ga64, -1, He4, Cu61, -1, Ga64, Cu61, k_Ga64_n_to_He4_Cu61_reaclib, k_Cu61_He4_to_n_Ga64_derived,
        H1, Ga64, -1, He4, Zn61, -1, Ga64, Zn61, k_Ga64_p_to_He4_Zn61_reaclib, k_Zn61_He4_to_p_Ga64_derived,
        N, Ge63, -1, H1, Ga63, -1, Ge63, Ga63, k_Ge63_n_to_p_Ga63_reaclib, k_Ga63_p_to_n_Ge63_derived,
        N, Ge63, -1, He4, Zn60, -1, Ge63, Zn60, k_Ge63_n_to_He4_Zn60_reaclib, k_Zn60_He4_to_n_Ge63_derived,
        N, Ge64, -1, H1, Ga64, -1, Ge64, Ga64, k_Ge64_n_to_p_Ga64_reaclib, k_Ga64_p_to_n_Ge64_derived,
        N, Ge64, -1, He4, Zn61, -1, Ge64, Zn61, k_Ge64_n_to_He4_Zn61_reaclib, k_Zn61_He4_to_n_Ge64_derived,
        He4, He4, He4, C12, -1, -1, C12, -1, k_He4_He4_He4_to_C12_reaclib, k_C12_to_3He4_derived,
        N, N13, -1, H1, C13, -1, N13, C13, k_N13_n_to_p_C13_derived, k_C13_p_to_n_N13_reaclib,
        N, O14, -1, H1, N14, -1, O14, N14, k_O14_n_to_p_N14_derived, k_N14_p_to_n_O14_reaclib,
        N, O15, -1, H1, N15, -1, O15, N15, k_O15_n_to_p_N15_derived, k_N15_p_to_n_O15_reaclib,
        N, F17, -1, He4, N14, -1, F17, N14, k_F17_n_to_He4_N14_derived, k_N14_He4_to_n_F17_reaclib,
        He4, F17, -1, H1, Ne20, -1, F17, Ne20, k_F17_He4_to_p_Ne20_derived, k_Ne20_p_to_He4_F17_reaclib,
        N, Ne19, -1, H1, F19, -1, Ne19, F19, k_Ne19_n_to_p_F19_derived, k_F19_p_to_n_Ne19_reaclib,
        N, Ne21, -1, He4, O18, -1, Ne21, O18, k_Ne21_n_to_He4_O18_derived, k_O18_He4_to_n_Ne21_reaclib,
        He4, Na21, -1, H1, Mg24, -1, Na21, Mg24, k_Na21_He4_to_p_Mg24_derived, k_Mg24_p_to_He4_Na21_reaclib,
        N, Mg23, -1, H1, Na23, -1, Mg23, Na23, k_Mg23_n_to_p_Na23_derived, k_Na23_p_to_n_Mg23_reaclib,
        N, Mg25, -1, He4, Ne22, -1, Mg25, Ne22, k_Mg25_n_to_He4_Ne22_derived, k_Ne22_He4_to_n_Mg25_reaclib,
        N, P30, -1, He4, Al27, -1, P30, Al27, k_P30_n_to_He4_Al27_derived, k_Al27_He4_to_n_P30_reaclib,
        H1, Cl37, -1, He4, S34, -1, Cl37, S34, k_Cl37_p_to_He4_S34_derived, k_S34_He4_to_p_Cl37_reaclib,
        N, Ar37, -1, He4, S34, -1, Ar37, S34, k_Ar37_n_to_He4_S34_derived, k_S34_He4_to_n_Ar37_reaclib,
        N, Ca41, -1, H1, K41, -1, Ca41, K41, k_Ca41_n_to_p_K41_derived, k_K41_p_to_n_Ca41_reaclib,
        N, Sc44, -1, He4, K41, -1, Sc44, K41, k_Sc44_n_to_He4_K41_derived, k_K41_He4_to_n_Sc44_reaclib,
        H1, V47, -1, He4, Ti44, -1, V47, Ti44, k_V47_p_to_He4_Ti44_derived, k_Ti44_He4_to_p_V47_reaclib
    };
}
#endif

#endif
