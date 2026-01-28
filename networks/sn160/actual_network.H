#ifndef actual_network_H
#define actual_network_H

#include <string_view>

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Loop.H>

#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex::literals;

void actual_network_init();

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
        k_Be7_to_Li7_reaclib = 1,
        k_C14_to_N14_reaclib = 2,
        k_N13_to_C13_reaclib = 3,
        k_O14_to_N14_reaclib = 4,
        k_O15_to_N15_reaclib = 5,
        k_F17_to_O17_reaclib = 6,
        k_F18_to_O18_reaclib = 7,
        k_Ne18_to_F18_reaclib = 8,
        k_Ne19_to_F19_reaclib = 9,
        k_Zn59_to_Cu59_reaclib = 10,
        k_Ge63_to_Ga63_reaclib = 11,
        k_d_to_n_p_reaclib = 12,
        k_He3_to_p_d_reaclib = 13,
        k_He4_to_n_He3_reaclib = 14,
        k_He4_to_d_d_reaclib = 15,
        k_Li6_to_He4_d_reaclib = 16,
        k_Li7_to_n_Li6_reaclib = 17,
        k_Be7_to_p_Li6_reaclib = 18,
        k_Be7_to_He4_He3_reaclib = 19,
        k_B8_to_p_Be7_reaclib = 20,
        k_B8_to_He4_He4_reaclib = 21,
        k_B10_to_p_Be9_reaclib = 22,
        k_B10_to_He4_Li6_reaclib = 23,
        k_B11_to_n_B10_reaclib = 24,
        k_B11_to_He4_Li7_reaclib = 25,
        k_C12_to_p_B11_reaclib = 26,
        k_C13_to_n_C12_reaclib = 27,
        k_C14_to_n_C13_reaclib = 28,
        k_N13_to_p_C12_reaclib = 29,
        k_N14_to_n_N13_reaclib = 30,
        k_N14_to_p_C13_reaclib = 31,
        k_N15_to_n_N14_reaclib = 32,
        k_N15_to_p_C14_reaclib = 33,
        k_O14_to_p_N13_reaclib = 34,
        k_O15_to_n_O14_reaclib = 35,
        k_O15_to_p_N14_reaclib = 36,
        k_O16_to_n_O15_reaclib = 37,
        k_O16_to_p_N15_reaclib = 38,
        k_O16_to_He4_C12_reaclib = 39,
        k_O17_to_n_O16_reaclib = 40,
        k_O18_to_n_O17_reaclib = 41,
        k_O18_to_He4_C14_reaclib = 42,
        k_F17_to_p_O16_reaclib = 43,
        k_F18_to_n_F17_reaclib = 44,
        k_F18_to_p_O17_reaclib = 45,
        k_F18_to_He4_N14_reaclib = 46,
        k_F19_to_n_F18_reaclib = 47,
        k_F19_to_p_O18_reaclib = 48,
        k_F19_to_He4_N15_reaclib = 49,
        k_Ne18_to_p_F17_reaclib = 50,
        k_Ne18_to_He4_O14_reaclib = 51,
        k_Ne19_to_n_Ne18_reaclib = 52,
        k_Ne19_to_p_F18_reaclib = 53,
        k_Ne19_to_He4_O15_reaclib = 54,
        k_Ne20_to_n_Ne19_reaclib = 55,
        k_Ne20_to_p_F19_reaclib = 56,
        k_Ne20_to_He4_O16_reaclib = 57,
        k_Ne21_to_n_Ne20_reaclib = 58,
        k_Ne21_to_He4_O17_reaclib = 59,
        k_Ne22_to_n_Ne21_reaclib = 60,
        k_Ne22_to_He4_O18_reaclib = 61,
        k_Na21_to_p_Ne20_reaclib = 62,
        k_Na21_to_He4_F17_reaclib = 63,
        k_Na22_to_n_Na21_reaclib = 64,
        k_Na22_to_p_Ne21_reaclib = 65,
        k_Na22_to_He4_F18_reaclib = 66,
        k_Na23_to_n_Na22_reaclib = 67,
        k_Na23_to_p_Ne22_reaclib = 68,
        k_Na23_to_He4_F19_reaclib = 69,
        k_Mg23_to_p_Na22_reaclib = 70,
        k_Mg23_to_He4_Ne19_reaclib = 71,
        k_Mg24_to_n_Mg23_reaclib = 72,
        k_Mg24_to_p_Na23_reaclib = 73,
        k_Mg24_to_He4_Ne20_reaclib = 74,
        k_Mg25_to_n_Mg24_reaclib = 75,
        k_Mg25_to_He4_Ne21_reaclib = 76,
        k_Mg26_to_n_Mg25_reaclib = 77,
        k_Mg26_to_He4_Ne22_reaclib = 78,
        k_Al25_to_p_Mg24_reaclib = 79,
        k_Al25_to_He4_Na21_reaclib = 80,
        k_Al26_to_n_Al25_reaclib = 81,
        k_Al26_to_p_Mg25_reaclib = 82,
        k_Al26_to_He4_Na22_reaclib = 83,
        k_Al27_to_n_Al26_reaclib = 84,
        k_Al27_to_p_Mg26_reaclib = 85,
        k_Al27_to_He4_Na23_reaclib = 86,
        k_Si28_to_p_Al27_reaclib = 87,
        k_Si28_to_He4_Mg24_reaclib = 88,
        k_Si29_to_n_Si28_reaclib = 89,
        k_Si29_to_He4_Mg25_reaclib = 90,
        k_Si30_to_n_Si29_reaclib = 91,
        k_Si30_to_He4_Mg26_reaclib = 92,
        k_Si31_to_n_Si30_reaclib = 93,
        k_Si32_to_n_Si31_reaclib = 94,
        k_P29_to_p_Si28_reaclib = 95,
        k_P29_to_He4_Al25_reaclib = 96,
        k_P30_to_n_P29_reaclib = 97,
        k_P30_to_p_Si29_reaclib = 98,
        k_P30_to_He4_Al26_reaclib = 99,
        k_P31_to_n_P30_reaclib = 100,
        k_P31_to_p_Si30_reaclib = 101,
        k_P31_to_He4_Al27_reaclib = 102,
        k_P32_to_n_P31_reaclib = 103,
        k_P32_to_p_Si31_reaclib = 104,
        k_P33_to_n_P32_reaclib = 105,
        k_P33_to_p_Si32_reaclib = 106,
        k_S32_to_p_P31_reaclib = 107,
        k_S32_to_He4_Si28_reaclib = 108,
        k_S33_to_n_S32_reaclib = 109,
        k_S33_to_p_P32_reaclib = 110,
        k_S33_to_He4_Si29_reaclib = 111,
        k_S34_to_n_S33_reaclib = 112,
        k_S34_to_p_P33_reaclib = 113,
        k_S34_to_He4_Si30_reaclib = 114,
        k_S35_to_n_S34_reaclib = 115,
        k_S35_to_He4_Si31_reaclib = 116,
        k_S36_to_n_S35_reaclib = 117,
        k_S36_to_He4_Si32_reaclib = 118,
        k_Cl33_to_p_S32_reaclib = 119,
        k_Cl33_to_He4_P29_reaclib = 120,
        k_Cl34_to_n_Cl33_reaclib = 121,
        k_Cl34_to_p_S33_reaclib = 122,
        k_Cl34_to_He4_P30_reaclib = 123,
        k_Cl35_to_n_Cl34_reaclib = 124,
        k_Cl35_to_p_S34_reaclib = 125,
        k_Cl35_to_He4_P31_reaclib = 126,
        k_Cl36_to_n_Cl35_reaclib = 127,
        k_Cl36_to_p_S35_reaclib = 128,
        k_Cl36_to_He4_P32_reaclib = 129,
        k_Cl37_to_n_Cl36_reaclib = 130,
        k_Cl37_to_p_S36_reaclib = 131,
        k_Cl37_to_He4_P33_reaclib = 132,
        k_Ar36_to_p_Cl35_reaclib = 133,
        k_Ar36_to_He4_S32_reaclib = 134,
        k_Ar37_to_n_Ar36_reaclib = 135,
        k_Ar37_to_p_Cl36_reaclib = 136,
        k_Ar37_to_He4_S33_reaclib = 137,
        k_Ar38_to_n_Ar37_reaclib = 138,
        k_Ar38_to_p_Cl37_reaclib = 139,
        k_Ar38_to_He4_S34_reaclib = 140,
        k_Ar39_to_n_Ar38_reaclib = 141,
        k_Ar39_to_He4_S35_reaclib = 142,
        k_Ar40_to_n_Ar39_reaclib = 143,
        k_Ar40_to_He4_S36_reaclib = 144,
        k_K37_to_p_Ar36_reaclib = 145,
        k_K37_to_He4_Cl33_reaclib = 146,
        k_K38_to_n_K37_reaclib = 147,
        k_K38_to_p_Ar37_reaclib = 148,
        k_K38_to_He4_Cl34_reaclib = 149,
        k_K39_to_n_K38_reaclib = 150,
        k_K39_to_p_Ar38_reaclib = 151,
        k_K39_to_He4_Cl35_reaclib = 152,
        k_K40_to_n_K39_reaclib = 153,
        k_K40_to_p_Ar39_reaclib = 154,
        k_K40_to_He4_Cl36_reaclib = 155,
        k_K41_to_n_K40_reaclib = 156,
        k_K41_to_p_Ar40_reaclib = 157,
        k_K41_to_He4_Cl37_reaclib = 158,
        k_Ca40_to_p_K39_reaclib = 159,
        k_Ca40_to_He4_Ar36_reaclib = 160,
        k_Ca41_to_n_Ca40_reaclib = 161,
        k_Ca41_to_p_K40_reaclib = 162,
        k_Ca41_to_He4_Ar37_reaclib = 163,
        k_Ca42_to_n_Ca41_reaclib = 164,
        k_Ca42_to_p_K41_reaclib = 165,
        k_Ca42_to_He4_Ar38_reaclib = 166,
        k_Ca43_to_n_Ca42_reaclib = 167,
        k_Ca43_to_He4_Ar39_reaclib = 168,
        k_Ca44_to_n_Ca43_reaclib = 169,
        k_Ca44_to_He4_Ar40_reaclib = 170,
        k_Ca45_to_n_Ca44_reaclib = 171,
        k_Ca46_to_n_Ca45_reaclib = 172,
        k_Ca47_to_n_Ca46_reaclib = 173,
        k_Ca48_to_n_Ca47_reaclib = 174,
        k_Sc43_to_p_Ca42_reaclib = 175,
        k_Sc43_to_He4_K39_reaclib = 176,
        k_Sc44_to_n_Sc43_reaclib = 177,
        k_Sc44_to_p_Ca43_reaclib = 178,
        k_Sc44_to_He4_K40_reaclib = 179,
        k_Sc45_to_n_Sc44_reaclib = 180,
        k_Sc45_to_p_Ca44_reaclib = 181,
        k_Sc45_to_He4_K41_reaclib = 182,
        k_Sc46_to_n_Sc45_reaclib = 183,
        k_Sc46_to_p_Ca45_reaclib = 184,
        k_Sc47_to_n_Sc46_reaclib = 185,
        k_Sc47_to_p_Ca46_reaclib = 186,
        k_Sc48_to_n_Sc47_reaclib = 187,
        k_Sc48_to_p_Ca47_reaclib = 188,
        k_Sc49_to_n_Sc48_reaclib = 189,
        k_Sc49_to_p_Ca48_reaclib = 190,
        k_Ti44_to_p_Sc43_reaclib = 191,
        k_Ti44_to_He4_Ca40_reaclib = 192,
        k_Ti45_to_n_Ti44_reaclib = 193,
        k_Ti45_to_p_Sc44_reaclib = 194,
        k_Ti45_to_He4_Ca41_reaclib = 195,
        k_Ti46_to_n_Ti45_reaclib = 196,
        k_Ti46_to_p_Sc45_reaclib = 197,
        k_Ti46_to_He4_Ca42_reaclib = 198,
        k_Ti47_to_n_Ti46_reaclib = 199,
        k_Ti47_to_p_Sc46_reaclib = 200,
        k_Ti47_to_He4_Ca43_reaclib = 201,
        k_Ti48_to_n_Ti47_reaclib = 202,
        k_Ti48_to_p_Sc47_reaclib = 203,
        k_Ti48_to_He4_Ca44_reaclib = 204,
        k_Ti49_to_n_Ti48_reaclib = 205,
        k_Ti49_to_p_Sc48_reaclib = 206,
        k_Ti49_to_He4_Ca45_reaclib = 207,
        k_Ti50_to_n_Ti49_reaclib = 208,
        k_Ti50_to_p_Sc49_reaclib = 209,
        k_Ti50_to_He4_Ca46_reaclib = 210,
        k_Ti51_to_n_Ti50_reaclib = 211,
        k_Ti51_to_He4_Ca47_reaclib = 212,
        k_V46_to_p_Ti45_reaclib = 213,
        k_V47_to_n_V46_reaclib = 214,
        k_V47_to_p_Ti46_reaclib = 215,
        k_V47_to_He4_Sc43_reaclib = 216,
        k_V48_to_n_V47_reaclib = 217,
        k_V48_to_p_Ti47_reaclib = 218,
        k_V48_to_He4_Sc44_reaclib = 219,
        k_V49_to_n_V48_reaclib = 220,
        k_V49_to_p_Ti48_reaclib = 221,
        k_V49_to_He4_Sc45_reaclib = 222,
        k_V50_to_n_V49_reaclib = 223,
        k_V50_to_p_Ti49_reaclib = 224,
        k_V50_to_He4_Sc46_reaclib = 225,
        k_V51_to_n_V50_reaclib = 226,
        k_V51_to_p_Ti50_reaclib = 227,
        k_V51_to_He4_Sc47_reaclib = 228,
        k_V52_to_n_V51_reaclib = 229,
        k_V52_to_p_Ti51_reaclib = 230,
        k_V52_to_He4_Sc48_reaclib = 231,
        k_Cr48_to_p_V47_reaclib = 232,
        k_Cr48_to_He4_Ti44_reaclib = 233,
        k_Cr49_to_n_Cr48_reaclib = 234,
        k_Cr49_to_p_V48_reaclib = 235,
        k_Cr49_to_He4_Ti45_reaclib = 236,
        k_Cr50_to_n_Cr49_reaclib = 237,
        k_Cr50_to_p_V49_reaclib = 238,
        k_Cr50_to_He4_Ti46_reaclib = 239,
        k_Cr51_to_n_Cr50_reaclib = 240,
        k_Cr51_to_p_V50_reaclib = 241,
        k_Cr51_to_He4_Ti47_reaclib = 242,
        k_Cr52_to_n_Cr51_reaclib = 243,
        k_Cr52_to_p_V51_reaclib = 244,
        k_Cr52_to_He4_Ti48_reaclib = 245,
        k_Cr53_to_n_Cr52_reaclib = 246,
        k_Cr53_to_p_V52_reaclib = 247,
        k_Cr53_to_He4_Ti49_reaclib = 248,
        k_Cr54_to_n_Cr53_reaclib = 249,
        k_Cr54_to_He4_Ti50_reaclib = 250,
        k_Mn50_to_p_Cr49_reaclib = 251,
        k_Mn50_to_He4_V46_reaclib = 252,
        k_Mn51_to_n_Mn50_reaclib = 253,
        k_Mn51_to_p_Cr50_reaclib = 254,
        k_Mn51_to_He4_V47_reaclib = 255,
        k_Mn52_to_n_Mn51_reaclib = 256,
        k_Mn52_to_p_Cr51_reaclib = 257,
        k_Mn52_to_He4_V48_reaclib = 258,
        k_Mn53_to_n_Mn52_reaclib = 259,
        k_Mn53_to_p_Cr52_reaclib = 260,
        k_Mn53_to_He4_V49_reaclib = 261,
        k_Mn54_to_n_Mn53_reaclib = 262,
        k_Mn54_to_p_Cr53_reaclib = 263,
        k_Mn54_to_He4_V50_reaclib = 264,
        k_Mn55_to_n_Mn54_reaclib = 265,
        k_Mn55_to_p_Cr54_reaclib = 266,
        k_Mn55_to_He4_V51_reaclib = 267,
        k_Fe52_to_p_Mn51_reaclib = 268,
        k_Fe52_to_He4_Cr48_reaclib = 269,
        k_Fe53_to_n_Fe52_reaclib = 270,
        k_Fe53_to_p_Mn52_reaclib = 271,
        k_Fe53_to_He4_Cr49_reaclib = 272,
        k_Fe54_to_n_Fe53_reaclib = 273,
        k_Fe54_to_p_Mn53_reaclib = 274,
        k_Fe54_to_He4_Cr50_reaclib = 275,
        k_Fe55_to_n_Fe54_reaclib = 276,
        k_Fe55_to_p_Mn54_reaclib = 277,
        k_Fe55_to_He4_Cr51_reaclib = 278,
        k_Fe56_to_n_Fe55_reaclib = 279,
        k_Fe56_to_p_Mn55_reaclib = 280,
        k_Fe56_to_He4_Cr52_reaclib = 281,
        k_Fe57_to_n_Fe56_reaclib = 282,
        k_Fe57_to_He4_Cr53_reaclib = 283,
        k_Fe58_to_n_Fe57_reaclib = 284,
        k_Fe58_to_He4_Cr54_reaclib = 285,
        k_Co53_to_p_Fe52_reaclib = 286,
        k_Co54_to_n_Co53_reaclib = 287,
        k_Co54_to_p_Fe53_reaclib = 288,
        k_Co54_to_He4_Mn50_reaclib = 289,
        k_Co55_to_n_Co54_reaclib = 290,
        k_Co55_to_p_Fe54_reaclib = 291,
        k_Co55_to_He4_Mn51_reaclib = 292,
        k_Co56_to_n_Co55_reaclib = 293,
        k_Co56_to_p_Fe55_reaclib = 294,
        k_Co56_to_He4_Mn52_reaclib = 295,
        k_Co57_to_n_Co56_reaclib = 296,
        k_Co57_to_p_Fe56_reaclib = 297,
        k_Co57_to_He4_Mn53_reaclib = 298,
        k_Co58_to_n_Co57_reaclib = 299,
        k_Co58_to_p_Fe57_reaclib = 300,
        k_Co58_to_He4_Mn54_reaclib = 301,
        k_Co59_to_n_Co58_reaclib = 302,
        k_Co59_to_p_Fe58_reaclib = 303,
        k_Co59_to_He4_Mn55_reaclib = 304,
        k_Ni56_to_p_Co55_reaclib = 305,
        k_Ni56_to_He4_Fe52_reaclib = 306,
        k_Ni57_to_n_Ni56_reaclib = 307,
        k_Ni57_to_p_Co56_reaclib = 308,
        k_Ni57_to_He4_Fe53_reaclib = 309,
        k_Ni58_to_n_Ni57_reaclib = 310,
        k_Ni58_to_p_Co57_reaclib = 311,
        k_Ni58_to_He4_Fe54_reaclib = 312,
        k_Ni59_to_n_Ni58_reaclib = 313,
        k_Ni59_to_p_Co58_reaclib = 314,
        k_Ni59_to_He4_Fe55_reaclib = 315,
        k_Ni60_to_n_Ni59_reaclib = 316,
        k_Ni60_to_p_Co59_reaclib = 317,
        k_Ni60_to_He4_Fe56_reaclib = 318,
        k_Ni61_to_n_Ni60_reaclib = 319,
        k_Ni61_to_He4_Fe57_reaclib = 320,
        k_Ni62_to_n_Ni61_reaclib = 321,
        k_Ni62_to_He4_Fe58_reaclib = 322,
        k_Ni63_to_n_Ni62_reaclib = 323,
        k_Ni64_to_n_Ni63_reaclib = 324,
        k_Cu57_to_p_Ni56_reaclib = 325,
        k_Cu57_to_He4_Co53_reaclib = 326,
        k_Cu58_to_n_Cu57_reaclib = 327,
        k_Cu58_to_p_Ni57_reaclib = 328,
        k_Cu58_to_He4_Co54_reaclib = 329,
        k_Cu59_to_n_Cu58_reaclib = 330,
        k_Cu59_to_p_Ni58_reaclib = 331,
        k_Cu59_to_He4_Co55_reaclib = 332,
        k_Cu60_to_n_Cu59_reaclib = 333,
        k_Cu60_to_p_Ni59_reaclib = 334,
        k_Cu60_to_He4_Co56_reaclib = 335,
        k_Cu61_to_n_Cu60_reaclib = 336,
        k_Cu61_to_p_Ni60_reaclib = 337,
        k_Cu61_to_He4_Co57_reaclib = 338,
        k_Cu62_to_n_Cu61_reaclib = 339,
        k_Cu62_to_p_Ni61_reaclib = 340,
        k_Cu62_to_He4_Co58_reaclib = 341,
        k_Cu63_to_n_Cu62_reaclib = 342,
        k_Cu63_to_p_Ni62_reaclib = 343,
        k_Cu63_to_He4_Co59_reaclib = 344,
        k_Cu64_to_n_Cu63_reaclib = 345,
        k_Cu64_to_p_Ni63_reaclib = 346,
        k_Cu65_to_n_Cu64_reaclib = 347,
        k_Cu65_to_p_Ni64_reaclib = 348,
        k_Zn59_to_p_Cu58_reaclib = 349,
        k_Zn59_to_p_Ni58_reaclib = 350,
        k_Zn60_to_n_Zn59_reaclib = 351,
        k_Zn60_to_p_Cu59_reaclib = 352,
        k_Zn60_to_He4_Ni56_reaclib = 353,
        k_Zn61_to_n_Zn60_reaclib = 354,
        k_Zn61_to_p_Cu60_reaclib = 355,
        k_Zn61_to_He4_Ni57_reaclib = 356,
        k_Zn62_to_n_Zn61_reaclib = 357,
        k_Zn62_to_p_Cu61_reaclib = 358,
        k_Zn62_to_He4_Ni58_reaclib = 359,
        k_Zn63_to_n_Zn62_reaclib = 360,
        k_Zn63_to_p_Cu62_reaclib = 361,
        k_Zn63_to_He4_Ni59_reaclib = 362,
        k_Zn64_to_n_Zn63_reaclib = 363,
        k_Zn64_to_p_Cu63_reaclib = 364,
        k_Zn64_to_He4_Ni60_reaclib = 365,
        k_Zn65_to_n_Zn64_reaclib = 366,
        k_Zn65_to_p_Cu64_reaclib = 367,
        k_Zn65_to_He4_Ni61_reaclib = 368,
        k_Zn66_to_n_Zn65_reaclib = 369,
        k_Zn66_to_p_Cu65_reaclib = 370,
        k_Zn66_to_He4_Ni62_reaclib = 371,
        k_Ga62_to_p_Zn61_reaclib = 372,
        k_Ga62_to_He4_Cu58_reaclib = 373,
        k_Ga63_to_n_Ga62_reaclib = 374,
        k_Ga63_to_p_Zn62_reaclib = 375,
        k_Ga63_to_He4_Cu59_reaclib = 376,
        k_Ga64_to_n_Ga63_reaclib = 377,
        k_Ga64_to_p_Zn63_reaclib = 378,
        k_Ga64_to_He4_Cu60_reaclib = 379,
        k_Ge63_to_p_Ga62_reaclib = 380,
        k_Ge63_to_He4_Zn59_reaclib = 381,
        k_Ge64_to_n_Ge63_reaclib = 382,
        k_Ge64_to_p_Ga63_reaclib = 383,
        k_Ge64_to_He4_Zn60_reaclib = 384,
        k_Li6_to_n_p_He4_reaclib = 385,
        k_Be9_to_n_He4_He4_reaclib = 386,
        k_C12_to_He4_He4_He4_reaclib = 387,
        k_n_p_to_d_reaclib = 388,
        k_p_p_to_d_reaclib_bet_pos = 389,
        k_p_p_to_d_reaclib_electron_capture = 390,
        k_p_d_to_He3_reaclib = 391,
        k_d_d_to_He4_reaclib = 392,
        k_He4_d_to_Li6_reaclib = 393,
        k_n_He3_to_He4_reaclib = 394,
        k_p_He3_to_He4_reaclib = 395,
        k_He4_He3_to_Be7_reaclib = 396,
        k_n_Li6_to_Li7_reaclib = 397,
        k_p_Li6_to_Be7_reaclib = 398,
        k_He4_Li6_to_B10_reaclib = 399,
        k_He4_Li7_to_B11_reaclib = 400,
        k_p_Be7_to_B8_reaclib = 401,
        k_p_Be9_to_B10_reaclib = 402,
        k_n_B10_to_B11_reaclib = 403,
        k_p_B11_to_C12_reaclib = 404,
        k_n_C12_to_C13_reaclib = 405,
        k_p_C12_to_N13_reaclib = 406,
        k_He4_C12_to_O16_reaclib = 407,
        k_n_C13_to_C14_reaclib = 408,
        k_p_C13_to_N14_reaclib = 409,
        k_p_C14_to_N15_reaclib = 410,
        k_He4_C14_to_O18_reaclib = 411,
        k_n_N13_to_N14_reaclib = 412,
        k_p_N13_to_O14_reaclib = 413,
        k_n_N14_to_N15_reaclib = 414,
        k_p_N14_to_O15_reaclib = 415,
        k_He4_N14_to_F18_reaclib = 416,
        k_p_N15_to_O16_reaclib = 417,
        k_He4_N15_to_F19_reaclib = 418,
        k_n_O14_to_O15_reaclib = 419,
        k_He4_O14_to_Ne18_reaclib = 420,
        k_n_O15_to_O16_reaclib = 421,
        k_He4_O15_to_Ne19_reaclib = 422,
        k_n_O16_to_O17_reaclib = 423,
        k_p_O16_to_F17_reaclib = 424,
        k_He4_O16_to_Ne20_reaclib = 425,
        k_n_O17_to_O18_reaclib = 426,
        k_p_O17_to_F18_reaclib = 427,
        k_He4_O17_to_Ne21_reaclib = 428,
        k_p_O18_to_F19_reaclib = 429,
        k_He4_O18_to_Ne22_reaclib = 430,
        k_n_F17_to_F18_reaclib = 431,
        k_p_F17_to_Ne18_reaclib = 432,
        k_He4_F17_to_Na21_reaclib = 433,
        k_n_F18_to_F19_reaclib = 434,
        k_p_F18_to_Ne19_reaclib = 435,
        k_He4_F18_to_Na22_reaclib = 436,
        k_p_F19_to_Ne20_reaclib = 437,
        k_He4_F19_to_Na23_reaclib = 438,
        k_n_Ne18_to_Ne19_reaclib = 439,
        k_n_Ne19_to_Ne20_reaclib = 440,
        k_He4_Ne19_to_Mg23_reaclib = 441,
        k_n_Ne20_to_Ne21_reaclib = 442,
        k_p_Ne20_to_Na21_reaclib = 443,
        k_He4_Ne20_to_Mg24_reaclib = 444,
        k_n_Ne21_to_Ne22_reaclib = 445,
        k_p_Ne21_to_Na22_reaclib = 446,
        k_He4_Ne21_to_Mg25_reaclib = 447,
        k_p_Ne22_to_Na23_reaclib = 448,
        k_He4_Ne22_to_Mg26_reaclib = 449,
        k_n_Na21_to_Na22_reaclib = 450,
        k_He4_Na21_to_Al25_reaclib = 451,
        k_n_Na22_to_Na23_reaclib = 452,
        k_p_Na22_to_Mg23_reaclib = 453,
        k_He4_Na22_to_Al26_reaclib = 454,
        k_p_Na23_to_Mg24_reaclib = 455,
        k_He4_Na23_to_Al27_reaclib = 456,
        k_n_Mg23_to_Mg24_reaclib = 457,
        k_n_Mg24_to_Mg25_reaclib = 458,
        k_p_Mg24_to_Al25_reaclib = 459,
        k_He4_Mg24_to_Si28_reaclib = 460,
        k_n_Mg25_to_Mg26_reaclib = 461,
        k_p_Mg25_to_Al26_reaclib = 462,
        k_He4_Mg25_to_Si29_reaclib = 463,
        k_p_Mg26_to_Al27_reaclib = 464,
        k_He4_Mg26_to_Si30_reaclib = 465,
        k_n_Al25_to_Al26_reaclib = 466,
        k_He4_Al25_to_P29_reaclib = 467,
        k_n_Al26_to_Al27_reaclib = 468,
        k_He4_Al26_to_P30_reaclib = 469,
        k_p_Al27_to_Si28_reaclib = 470,
        k_He4_Al27_to_P31_reaclib = 471,
        k_n_Si28_to_Si29_reaclib = 472,
        k_p_Si28_to_P29_reaclib = 473,
        k_He4_Si28_to_S32_reaclib = 474,
        k_n_Si29_to_Si30_reaclib = 475,
        k_p_Si29_to_P30_reaclib = 476,
        k_He4_Si29_to_S33_reaclib = 477,
        k_n_Si30_to_Si31_reaclib = 478,
        k_p_Si30_to_P31_reaclib = 479,
        k_He4_Si30_to_S34_reaclib = 480,
        k_n_Si31_to_Si32_reaclib = 481,
        k_p_Si31_to_P32_reaclib = 482,
        k_He4_Si31_to_S35_reaclib = 483,
        k_p_Si32_to_P33_reaclib = 484,
        k_He4_Si32_to_S36_reaclib = 485,
        k_n_P29_to_P30_reaclib = 486,
        k_He4_P29_to_Cl33_reaclib = 487,
        k_n_P30_to_P31_reaclib = 488,
        k_He4_P30_to_Cl34_reaclib = 489,
        k_n_P31_to_P32_reaclib = 490,
        k_p_P31_to_S32_reaclib = 491,
        k_He4_P31_to_Cl35_reaclib = 492,
        k_n_P32_to_P33_reaclib = 493,
        k_p_P32_to_S33_reaclib = 494,
        k_He4_P32_to_Cl36_reaclib = 495,
        k_p_P33_to_S34_reaclib = 496,
        k_He4_P33_to_Cl37_reaclib = 497,
        k_n_S32_to_S33_reaclib = 498,
        k_p_S32_to_Cl33_reaclib = 499,
        k_He4_S32_to_Ar36_reaclib = 500,
        k_n_S33_to_S34_reaclib = 501,
        k_p_S33_to_Cl34_reaclib = 502,
        k_He4_S33_to_Ar37_reaclib = 503,
        k_n_S34_to_S35_reaclib = 504,
        k_p_S34_to_Cl35_reaclib = 505,
        k_He4_S34_to_Ar38_reaclib = 506,
        k_n_S35_to_S36_reaclib = 507,
        k_p_S35_to_Cl36_reaclib = 508,
        k_He4_S35_to_Ar39_reaclib = 509,
        k_p_S36_to_Cl37_reaclib = 510,
        k_He4_S36_to_Ar40_reaclib = 511,
        k_n_Cl33_to_Cl34_reaclib = 512,
        k_He4_Cl33_to_K37_reaclib = 513,
        k_n_Cl34_to_Cl35_reaclib = 514,
        k_He4_Cl34_to_K38_reaclib = 515,
        k_n_Cl35_to_Cl36_reaclib = 516,
        k_p_Cl35_to_Ar36_reaclib = 517,
        k_He4_Cl35_to_K39_reaclib = 518,
        k_n_Cl36_to_Cl37_reaclib = 519,
        k_p_Cl36_to_Ar37_reaclib = 520,
        k_He4_Cl36_to_K40_reaclib = 521,
        k_p_Cl37_to_Ar38_reaclib = 522,
        k_He4_Cl37_to_K41_reaclib = 523,
        k_n_Ar36_to_Ar37_reaclib = 524,
        k_p_Ar36_to_K37_reaclib = 525,
        k_He4_Ar36_to_Ca40_reaclib = 526,
        k_n_Ar37_to_Ar38_reaclib = 527,
        k_p_Ar37_to_K38_reaclib = 528,
        k_He4_Ar37_to_Ca41_reaclib = 529,
        k_n_Ar38_to_Ar39_reaclib = 530,
        k_p_Ar38_to_K39_reaclib = 531,
        k_He4_Ar38_to_Ca42_reaclib = 532,
        k_n_Ar39_to_Ar40_reaclib = 533,
        k_p_Ar39_to_K40_reaclib = 534,
        k_He4_Ar39_to_Ca43_reaclib = 535,
        k_p_Ar40_to_K41_reaclib = 536,
        k_He4_Ar40_to_Ca44_reaclib = 537,
        k_n_K37_to_K38_reaclib = 538,
        k_n_K38_to_K39_reaclib = 539,
        k_n_K39_to_K40_reaclib = 540,
        k_p_K39_to_Ca40_reaclib = 541,
        k_He4_K39_to_Sc43_reaclib = 542,
        k_n_K40_to_K41_reaclib = 543,
        k_p_K40_to_Ca41_reaclib = 544,
        k_He4_K40_to_Sc44_reaclib = 545,
        k_p_K41_to_Ca42_reaclib = 546,
        k_He4_K41_to_Sc45_reaclib = 547,
        k_n_Ca40_to_Ca41_reaclib = 548,
        k_He4_Ca40_to_Ti44_reaclib = 549,
        k_n_Ca41_to_Ca42_reaclib = 550,
        k_He4_Ca41_to_Ti45_reaclib = 551,
        k_n_Ca42_to_Ca43_reaclib = 552,
        k_p_Ca42_to_Sc43_reaclib = 553,
        k_He4_Ca42_to_Ti46_reaclib = 554,
        k_n_Ca43_to_Ca44_reaclib = 555,
        k_p_Ca43_to_Sc44_reaclib = 556,
        k_He4_Ca43_to_Ti47_reaclib = 557,
        k_n_Ca44_to_Ca45_reaclib = 558,
        k_p_Ca44_to_Sc45_reaclib = 559,
        k_He4_Ca44_to_Ti48_reaclib = 560,
        k_n_Ca45_to_Ca46_reaclib = 561,
        k_p_Ca45_to_Sc46_reaclib = 562,
        k_He4_Ca45_to_Ti49_reaclib = 563,
        k_n_Ca46_to_Ca47_reaclib = 564,
        k_p_Ca46_to_Sc47_reaclib = 565,
        k_He4_Ca46_to_Ti50_reaclib = 566,
        k_n_Ca47_to_Ca48_reaclib = 567,
        k_p_Ca47_to_Sc48_reaclib = 568,
        k_He4_Ca47_to_Ti51_reaclib = 569,
        k_p_Ca48_to_Sc49_reaclib = 570,
        k_n_Sc43_to_Sc44_reaclib = 571,
        k_p_Sc43_to_Ti44_reaclib = 572,
        k_He4_Sc43_to_V47_reaclib = 573,
        k_n_Sc44_to_Sc45_reaclib = 574,
        k_p_Sc44_to_Ti45_reaclib = 575,
        k_He4_Sc44_to_V48_reaclib = 576,
        k_n_Sc45_to_Sc46_reaclib = 577,
        k_p_Sc45_to_Ti46_reaclib = 578,
        k_He4_Sc45_to_V49_reaclib = 579,
        k_n_Sc46_to_Sc47_reaclib = 580,
        k_p_Sc46_to_Ti47_reaclib = 581,
        k_He4_Sc46_to_V50_reaclib = 582,
        k_n_Sc47_to_Sc48_reaclib = 583,
        k_p_Sc47_to_Ti48_reaclib = 584,
        k_He4_Sc47_to_V51_reaclib = 585,
        k_n_Sc48_to_Sc49_reaclib = 586,
        k_p_Sc48_to_Ti49_reaclib = 587,
        k_He4_Sc48_to_V52_reaclib = 588,
        k_p_Sc49_to_Ti50_reaclib = 589,
        k_n_Ti44_to_Ti45_reaclib = 590,
        k_He4_Ti44_to_Cr48_reaclib = 591,
        k_n_Ti45_to_Ti46_reaclib = 592,
        k_p_Ti45_to_V46_reaclib = 593,
        k_He4_Ti45_to_Cr49_reaclib = 594,
        k_n_Ti46_to_Ti47_reaclib = 595,
        k_p_Ti46_to_V47_reaclib = 596,
        k_He4_Ti46_to_Cr50_reaclib = 597,
        k_n_Ti47_to_Ti48_reaclib = 598,
        k_p_Ti47_to_V48_reaclib = 599,
        k_He4_Ti47_to_Cr51_reaclib = 600,
        k_n_Ti48_to_Ti49_reaclib = 601,
        k_p_Ti48_to_V49_reaclib = 602,
        k_He4_Ti48_to_Cr52_reaclib = 603,
        k_n_Ti49_to_Ti50_reaclib = 604,
        k_p_Ti49_to_V50_reaclib = 605,
        k_He4_Ti49_to_Cr53_reaclib = 606,
        k_n_Ti50_to_Ti51_reaclib = 607,
        k_p_Ti50_to_V51_reaclib = 608,
        k_He4_Ti50_to_Cr54_reaclib = 609,
        k_p_Ti51_to_V52_reaclib = 610,
        k_n_V46_to_V47_reaclib = 611,
        k_He4_V46_to_Mn50_reaclib = 612,
        k_n_V47_to_V48_reaclib = 613,
        k_p_V47_to_Cr48_reaclib = 614,
        k_He4_V47_to_Mn51_reaclib = 615,
        k_n_V48_to_V49_reaclib = 616,
        k_p_V48_to_Cr49_reaclib = 617,
        k_He4_V48_to_Mn52_reaclib = 618,
        k_n_V49_to_V50_reaclib = 619,
        k_p_V49_to_Cr50_reaclib = 620,
        k_He4_V49_to_Mn53_reaclib = 621,
        k_n_V50_to_V51_reaclib = 622,
        k_p_V50_to_Cr51_reaclib = 623,
        k_He4_V50_to_Mn54_reaclib = 624,
        k_n_V51_to_V52_reaclib = 625,
        k_p_V51_to_Cr52_reaclib = 626,
        k_He4_V51_to_Mn55_reaclib = 627,
        k_p_V52_to_Cr53_reaclib = 628,
        k_n_Cr48_to_Cr49_reaclib = 629,
        k_He4_Cr48_to_Fe52_reaclib = 630,
        k_n_Cr49_to_Cr50_reaclib = 631,
        k_p_Cr49_to_Mn50_reaclib = 632,
        k_He4_Cr49_to_Fe53_reaclib = 633,
        k_n_Cr50_to_Cr51_reaclib = 634,
        k_p_Cr50_to_Mn51_reaclib = 635,
        k_He4_Cr50_to_Fe54_reaclib = 636,
        k_n_Cr51_to_Cr52_reaclib = 637,
        k_p_Cr51_to_Mn52_reaclib = 638,
        k_He4_Cr51_to_Fe55_reaclib = 639,
        k_n_Cr52_to_Cr53_reaclib = 640,
        k_p_Cr52_to_Mn53_reaclib = 641,
        k_He4_Cr52_to_Fe56_reaclib = 642,
        k_n_Cr53_to_Cr54_reaclib = 643,
        k_p_Cr53_to_Mn54_reaclib = 644,
        k_He4_Cr53_to_Fe57_reaclib = 645,
        k_p_Cr54_to_Mn55_reaclib = 646,
        k_He4_Cr54_to_Fe58_reaclib = 647,
        k_n_Mn50_to_Mn51_reaclib = 648,
        k_He4_Mn50_to_Co54_reaclib = 649,
        k_n_Mn51_to_Mn52_reaclib = 650,
        k_p_Mn51_to_Fe52_reaclib = 651,
        k_He4_Mn51_to_Co55_reaclib = 652,
        k_n_Mn52_to_Mn53_reaclib = 653,
        k_p_Mn52_to_Fe53_reaclib = 654,
        k_He4_Mn52_to_Co56_reaclib = 655,
        k_n_Mn53_to_Mn54_reaclib = 656,
        k_p_Mn53_to_Fe54_reaclib = 657,
        k_He4_Mn53_to_Co57_reaclib = 658,
        k_n_Mn54_to_Mn55_reaclib = 659,
        k_p_Mn54_to_Fe55_reaclib = 660,
        k_He4_Mn54_to_Co58_reaclib = 661,
        k_p_Mn55_to_Fe56_reaclib = 662,
        k_He4_Mn55_to_Co59_reaclib = 663,
        k_n_Fe52_to_Fe53_reaclib = 664,
        k_p_Fe52_to_Co53_reaclib = 665,
        k_He4_Fe52_to_Ni56_reaclib = 666,
        k_n_Fe53_to_Fe54_reaclib = 667,
        k_p_Fe53_to_Co54_reaclib = 668,
        k_He4_Fe53_to_Ni57_reaclib = 669,
        k_n_Fe54_to_Fe55_reaclib = 670,
        k_p_Fe54_to_Co55_reaclib = 671,
        k_He4_Fe54_to_Ni58_reaclib = 672,
        k_n_Fe55_to_Fe56_reaclib = 673,
        k_p_Fe55_to_Co56_reaclib = 674,
        k_He4_Fe55_to_Ni59_reaclib = 675,
        k_n_Fe56_to_Fe57_reaclib = 676,
        k_p_Fe56_to_Co57_reaclib = 677,
        k_He4_Fe56_to_Ni60_reaclib = 678,
        k_n_Fe57_to_Fe58_reaclib = 679,
        k_p_Fe57_to_Co58_reaclib = 680,
        k_He4_Fe57_to_Ni61_reaclib = 681,
        k_p_Fe58_to_Co59_reaclib = 682,
        k_He4_Fe58_to_Ni62_reaclib = 683,
        k_n_Co53_to_Co54_reaclib = 684,
        k_He4_Co53_to_Cu57_reaclib = 685,
        k_n_Co54_to_Co55_reaclib = 686,
        k_He4_Co54_to_Cu58_reaclib = 687,
        k_n_Co55_to_Co56_reaclib = 688,
        k_p_Co55_to_Ni56_reaclib = 689,
        k_He4_Co55_to_Cu59_reaclib = 690,
        k_n_Co56_to_Co57_reaclib = 691,
        k_p_Co56_to_Ni57_reaclib = 692,
        k_He4_Co56_to_Cu60_reaclib = 693,
        k_n_Co57_to_Co58_reaclib = 694,
        k_p_Co57_to_Ni58_reaclib = 695,
        k_He4_Co57_to_Cu61_reaclib = 696,
        k_n_Co58_to_Co59_reaclib = 697,
        k_p_Co58_to_Ni59_reaclib = 698,
        k_He4_Co58_to_Cu62_reaclib = 699,
        k_p_Co59_to_Ni60_reaclib = 700,
        k_He4_Co59_to_Cu63_reaclib = 701,
        k_n_Ni56_to_Ni57_reaclib = 702,
        k_p_Ni56_to_Cu57_reaclib = 703,
        k_He4_Ni56_to_Zn60_reaclib = 704,
        k_n_Ni57_to_Ni58_reaclib = 705,
        k_p_Ni57_to_Cu58_reaclib = 706,
        k_He4_Ni57_to_Zn61_reaclib = 707,
        k_n_Ni58_to_Ni59_reaclib = 708,
        k_p_Ni58_to_Cu59_reaclib = 709,
        k_He4_Ni58_to_Zn62_reaclib = 710,
        k_n_Ni59_to_Ni60_reaclib = 711,
        k_p_Ni59_to_Cu60_reaclib = 712,
        k_He4_Ni59_to_Zn63_reaclib = 713,
        k_n_Ni60_to_Ni61_reaclib = 714,
        k_p_Ni60_to_Cu61_reaclib = 715,
        k_He4_Ni60_to_Zn64_reaclib = 716,
        k_n_Ni61_to_Ni62_reaclib = 717,
        k_p_Ni61_to_Cu62_reaclib = 718,
        k_He4_Ni61_to_Zn65_reaclib = 719,
        k_n_Ni62_to_Ni63_reaclib = 720,
        k_p_Ni62_to_Cu63_reaclib = 721,
        k_He4_Ni62_to_Zn66_reaclib = 722,
        k_n_Ni63_to_Ni64_reaclib = 723,
        k_p_Ni63_to_Cu64_reaclib = 724,
        k_p_Ni64_to_Cu65_reaclib = 725,
        k_n_Cu57_to_Cu58_reaclib = 726,
        k_n_Cu58_to_Cu59_reaclib = 727,
        k_p_Cu58_to_Zn59_reaclib = 728,
        k_He4_Cu58_to_Ga62_reaclib = 729,
        k_n_Cu59_to_Cu60_reaclib = 730,
        k_p_Cu59_to_Zn60_reaclib = 731,
        k_He4_Cu59_to_Ga63_reaclib = 732,
        k_n_Cu60_to_Cu61_reaclib = 733,
        k_p_Cu60_to_Zn61_reaclib = 734,
        k_He4_Cu60_to_Ga64_reaclib = 735,
        k_n_Cu61_to_Cu62_reaclib = 736,
        k_p_Cu61_to_Zn62_reaclib = 737,
        k_n_Cu62_to_Cu63_reaclib = 738,
        k_p_Cu62_to_Zn63_reaclib = 739,
        k_n_Cu63_to_Cu64_reaclib = 740,
        k_p_Cu63_to_Zn64_reaclib = 741,
        k_n_Cu64_to_Cu65_reaclib = 742,
        k_p_Cu64_to_Zn65_reaclib = 743,
        k_p_Cu65_to_Zn66_reaclib = 744,
        k_n_Zn59_to_Zn60_reaclib = 745,
        k_He4_Zn59_to_Ge63_reaclib = 746,
        k_n_Zn60_to_Zn61_reaclib = 747,
        k_He4_Zn60_to_Ge64_reaclib = 748,
        k_n_Zn61_to_Zn62_reaclib = 749,
        k_p_Zn61_to_Ga62_reaclib = 750,
        k_n_Zn62_to_Zn63_reaclib = 751,
        k_p_Zn62_to_Ga63_reaclib = 752,
        k_n_Zn63_to_Zn64_reaclib = 753,
        k_p_Zn63_to_Ga64_reaclib = 754,
        k_n_Zn64_to_Zn65_reaclib = 755,
        k_n_Zn65_to_Zn66_reaclib = 756,
        k_n_Ga62_to_Ga63_reaclib = 757,
        k_p_Ga62_to_Ge63_reaclib = 758,
        k_n_Ga63_to_Ga64_reaclib = 759,
        k_p_Ga63_to_Ge64_reaclib = 760,
        k_n_Ge63_to_Ge64_reaclib = 761,
        k_d_d_to_n_He3_reaclib = 762,
        k_n_He3_to_d_d_reaclib = 763,
        k_d_He3_to_p_He4_reaclib = 764,
        k_He4_He3_to_p_Li6_reaclib = 765,
        k_p_He4_to_d_He3_reaclib = 766,
        k_He4_He4_to_n_Be7_reaclib = 767,
        k_He4_He4_to_p_Li7_reaclib = 768,
        k_p_Li6_to_He4_He3_reaclib = 769,
        k_d_Li6_to_n_Be7_reaclib = 770,
        k_d_Li6_to_p_Li7_reaclib = 771,
        k_He4_Li6_to_p_Be9_reaclib = 772,
        k_p_Li7_to_n_Be7_reaclib = 773,
        k_p_Li7_to_d_Li6_reaclib = 774,
        k_p_Li7_to_He4_He4_reaclib = 775,
        k_He4_Li7_to_n_B10_reaclib = 776,
        k_n_Be7_to_p_Li7_reaclib = 777,
        k_n_Be7_to_d_Li6_reaclib = 778,
        k_n_Be7_to_He4_He4_reaclib = 779,
        k_He4_Be7_to_p_B10_reaclib = 780,
        k_p_Be9_to_He4_Li6_reaclib = 781,
        k_He4_Be9_to_n_C12_reaclib = 782,
        k_n_B10_to_He4_Li7_reaclib = 783,
        k_p_B10_to_He4_Be7_reaclib = 784,
        k_He4_B10_to_n_N13_reaclib = 785,
        k_He4_B10_to_p_C13_reaclib = 786,
        k_He4_B11_to_n_N14_reaclib = 787,
        k_He4_B11_to_p_C14_reaclib = 788,
        k_n_C12_to_He4_Be9_reaclib = 789,
        k_He4_C12_to_n_O15_reaclib = 790,
        k_He4_C12_to_p_N15_reaclib = 791,
        k_C12_C12_to_n_Mg23_reaclib = 792,
        k_C12_C12_to_p_Na23_reaclib = 793,
        k_C12_C12_to_He4_Ne20_reaclib = 794,
        k_p_C13_to_n_N13_reaclib = 795,
        k_p_C13_to_He4_B10_reaclib = 796,
        k_d_C13_to_n_N14_reaclib = 797,
        k_He4_C13_to_n_O16_reaclib = 798,
        k_p_C14_to_n_N14_reaclib = 799,
        k_p_C14_to_He4_B11_reaclib = 800,
        k_d_C14_to_n_N15_reaclib = 801,
        k_He4_C14_to_n_O17_reaclib = 802,
        k_n_N13_to_p_C13_reaclib = 803,
        k_n_N13_to_He4_B10_reaclib = 804,
        k_He4_N13_to_p_O16_reaclib = 805,
        k_n_N14_to_p_C14_reaclib = 806,
        k_n_N14_to_d_C13_reaclib = 807,
        k_n_N14_to_He4_B11_reaclib = 808,
        k_p_N14_to_n_O14_reaclib = 809,
        k_He4_N14_to_n_F17_reaclib = 810,
        k_He4_N14_to_p_O17_reaclib = 811,
        k_n_N15_to_d_C14_reaclib = 812,
        k_p_N15_to_n_O15_reaclib = 813,
        k_p_N15_to_He4_C12_reaclib = 814,
        k_He4_N15_to_n_F18_reaclib = 815,
        k_He4_N15_to_p_O18_reaclib = 816,
        k_n_O14_to_p_N14_reaclib = 817,
        k_He4_O14_to_p_F17_reaclib = 818,
        k_n_O15_to_p_N15_reaclib = 819,
        k_n_O15_to_He4_C12_reaclib = 820,
        k_He4_O15_to_n_Ne18_reaclib = 821,
        k_He4_O15_to_p_F18_reaclib = 822,
        k_n_O16_to_He4_C13_reaclib = 823,
        k_p_O16_to_He4_N13_reaclib = 824,
        k_He4_O16_to_n_Ne19_reaclib = 825,
        k_He4_O16_to_p_F19_reaclib = 826,
        k_C12_O16_to_p_Al27_reaclib = 827,
        k_C12_O16_to_He4_Mg24_reaclib = 828,
        k_O16_O16_to_p_P31_reaclib = 829,
        k_O16_O16_to_He4_Si28_reaclib = 830,
        k_n_O17_to_He4_C14_reaclib = 831,
        k_p_O17_to_n_F17_reaclib = 832,
        k_p_O17_to_He4_N14_reaclib = 833,
        k_He4_O17_to_n_Ne20_reaclib = 834,
        k_p_O18_to_n_F18_reaclib = 835,
        k_p_O18_to_He4_N15_reaclib = 836,
        k_He4_O18_to_n_Ne21_reaclib = 837,
        k_n_F17_to_p_O17_reaclib = 838,
        k_n_F17_to_He4_N14_reaclib = 839,
        k_p_F17_to_He4_O14_reaclib = 840,
        k_He4_F17_to_p_Ne20_reaclib = 841,
        k_n_F18_to_p_O18_reaclib = 842,
        k_n_F18_to_He4_N15_reaclib = 843,
        k_p_F18_to_n_Ne18_reaclib = 844,
        k_p_F18_to_He4_O15_reaclib = 845,
        k_He4_F18_to_n_Na21_reaclib = 846,
        k_He4_F18_to_p_Ne21_reaclib = 847,
        k_p_F19_to_n_Ne19_reaclib = 848,
        k_p_F19_to_He4_O16_reaclib = 849,
        k_He4_F19_to_n_Na22_reaclib = 850,
        k_He4_F19_to_p_Ne22_reaclib = 851,
        k_n_Ne18_to_p_F18_reaclib = 852,
        k_n_Ne18_to_He4_O15_reaclib = 853,
        k_He4_Ne18_to_p_Na21_reaclib = 854,
        k_n_Ne19_to_p_F19_reaclib = 855,
        k_n_Ne19_to_He4_O16_reaclib = 856,
        k_He4_Ne19_to_p_Na22_reaclib = 857,
        k_n_Ne20_to_He4_O17_reaclib = 858,
        k_p_Ne20_to_He4_F17_reaclib = 859,
        k_He4_Ne20_to_n_Mg23_reaclib = 860,
        k_He4_Ne20_to_p_Na23_reaclib = 861,
        k_He4_Ne20_to_C12_C12_reaclib = 862,
        k_C12_Ne20_to_p_P31_reaclib = 863,
        k_C12_Ne20_to_He4_Si28_reaclib = 864,
        k_n_Ne21_to_He4_O18_reaclib = 865,
        k_p_Ne21_to_n_Na21_reaclib = 866,
        k_p_Ne21_to_He4_F18_reaclib = 867,
        k_He4_Ne21_to_n_Mg24_reaclib = 868,
        k_p_Ne22_to_n_Na22_reaclib = 869,
        k_p_Ne22_to_He4_F19_reaclib = 870,
        k_He4_Ne22_to_n_Mg25_reaclib = 871,
        k_n_Na21_to_p_Ne21_reaclib = 872,
        k_n_Na21_to_He4_F18_reaclib = 873,
        k_p_Na21_to_He4_Ne18_reaclib = 874,
        k_He4_Na21_to_p_Mg24_reaclib = 875,
        k_n_Na22_to_p_Ne22_reaclib = 876,
        k_n_Na22_to_He4_F19_reaclib = 877,
        k_p_Na22_to_He4_Ne19_reaclib = 878,
        k_He4_Na22_to_n_Al25_reaclib = 879,
        k_He4_Na22_to_p_Mg25_reaclib = 880,
        k_p_Na23_to_n_Mg23_reaclib = 881,
        k_p_Na23_to_He4_Ne20_reaclib = 882,
        k_p_Na23_to_C12_C12_reaclib = 883,
        k_He4_Na23_to_n_Al26_reaclib = 884,
        k_He4_Na23_to_p_Mg26_reaclib = 885,
        k_n_Mg23_to_p_Na23_reaclib = 886,
        k_n_Mg23_to_He4_Ne20_reaclib = 887,
        k_n_Mg23_to_C12_C12_reaclib = 888,
        k_He4_Mg23_to_p_Al26_reaclib = 889,
        k_n_Mg24_to_He4_Ne21_reaclib = 890,
        k_p_Mg24_to_He4_Na21_reaclib = 891,
        k_He4_Mg24_to_p_Al27_reaclib = 892,
        k_He4_Mg24_to_C12_O16_reaclib = 893,
        k_n_Mg25_to_He4_Ne22_reaclib = 894,
        k_p_Mg25_to_n_Al25_reaclib = 895,
        k_p_Mg25_to_He4_Na22_reaclib = 896,
        k_He4_Mg25_to_n_Si28_reaclib = 897,
        k_p_Mg26_to_n_Al26_reaclib = 898,
        k_p_Mg26_to_He4_Na23_reaclib = 899,
        k_He4_Mg26_to_n_Si29_reaclib = 900,
        k_n_Al25_to_p_Mg25_reaclib = 901,
        k_n_Al25_to_He4_Na22_reaclib = 902,
        k_He4_Al25_to_p_Si28_reaclib = 903,
        k_n_Al26_to_p_Mg26_reaclib = 904,
        k_n_Al26_to_He4_Na23_reaclib = 905,
        k_p_Al26_to_He4_Mg23_reaclib = 906,
        k_He4_Al26_to_n_P29_reaclib = 907,
        k_He4_Al26_to_p_Si29_reaclib = 908,
        k_p_Al27_to_He4_Mg24_reaclib = 909,
        k_p_Al27_to_C12_O16_reaclib = 910,
        k_He4_Al27_to_n_P30_reaclib = 911,
        k_He4_Al27_to_p_Si30_reaclib = 912,
        k_n_Si28_to_He4_Mg25_reaclib = 913,
        k_p_Si28_to_He4_Al25_reaclib = 914,
        k_He4_Si28_to_p_P31_reaclib = 915,
        k_He4_Si28_to_C12_Ne20_reaclib = 916,
        k_He4_Si28_to_O16_O16_reaclib = 917,
        k_n_Si29_to_He4_Mg26_reaclib = 918,
        k_p_Si29_to_n_P29_reaclib = 919,
        k_p_Si29_to_He4_Al26_reaclib = 920,
        k_He4_Si29_to_n_S32_reaclib = 921,
        k_He4_Si29_to_p_P32_reaclib = 922,
        k_p_Si30_to_n_P30_reaclib = 923,
        k_p_Si30_to_He4_Al27_reaclib = 924,
        k_He4_Si30_to_n_S33_reaclib = 925,
        k_He4_Si30_to_p_P33_reaclib = 926,
        k_p_Si31_to_n_P31_reaclib = 927,
        k_He4_Si31_to_n_S34_reaclib = 928,
        k_p_Si32_to_n_P32_reaclib = 929,
        k_He4_Si32_to_n_S35_reaclib = 930,
        k_n_P29_to_p_Si29_reaclib = 931,
        k_n_P29_to_He4_Al26_reaclib = 932,
        k_He4_P29_to_p_S32_reaclib = 933,
        k_n_P30_to_p_Si30_reaclib = 934,
        k_n_P30_to_He4_Al27_reaclib = 935,
        k_He4_P30_to_n_Cl33_reaclib = 936,
        k_He4_P30_to_p_S33_reaclib = 937,
        k_n_P31_to_p_Si31_reaclib = 938,
        k_p_P31_to_He4_Si28_reaclib = 939,
        k_p_P31_to_C12_Ne20_reaclib = 940,
        k_p_P31_to_O16_O16_reaclib = 941,
        k_He4_P31_to_n_Cl34_reaclib = 942,
        k_He4_P31_to_p_S34_reaclib = 943,
        k_n_P32_to_p_Si32_reaclib = 944,
        k_p_P32_to_n_S32_reaclib = 945,
        k_p_P32_to_He4_Si29_reaclib = 946,
        k_He4_P32_to_n_Cl35_reaclib = 947,
        k_He4_P32_to_p_S35_reaclib = 948,
        k_p_P33_to_n_S33_reaclib = 949,
        k_p_P33_to_He4_Si30_reaclib = 950,
        k_He4_P33_to_n_Cl36_reaclib = 951,
        k_He4_P33_to_p_S36_reaclib = 952,
        k_n_S32_to_p_P32_reaclib = 953,
        k_n_S32_to_He4_Si29_reaclib = 954,
        k_p_S32_to_He4_P29_reaclib = 955,
        k_He4_S32_to_p_Cl35_reaclib = 956,
        k_n_S33_to_p_P33_reaclib = 957,
        k_n_S33_to_He4_Si30_reaclib = 958,
        k_p_S33_to_n_Cl33_reaclib = 959,
        k_p_S33_to_He4_P30_reaclib = 960,
        k_He4_S33_to_n_Ar36_reaclib = 961,
        k_He4_S33_to_p_Cl36_reaclib = 962,
        k_n_S34_to_He4_Si31_reaclib = 963,
        k_p_S34_to_n_Cl34_reaclib = 964,
        k_p_S34_to_He4_P31_reaclib = 965,
        k_He4_S34_to_n_Ar37_reaclib = 966,
        k_He4_S34_to_p_Cl37_reaclib = 967,
        k_n_S35_to_He4_Si32_reaclib = 968,
        k_p_S35_to_n_Cl35_reaclib = 969,
        k_p_S35_to_He4_P32_reaclib = 970,
        k_He4_S35_to_n_Ar38_reaclib = 971,
        k_p_S36_to_n_Cl36_reaclib = 972,
        k_p_S36_to_He4_P33_reaclib = 973,
        k_He4_S36_to_n_Ar39_reaclib = 974,
        k_n_Cl33_to_p_S33_reaclib = 975,
        k_n_Cl33_to_He4_P30_reaclib = 976,
        k_He4_Cl33_to_p_Ar36_reaclib = 977,
        k_n_Cl34_to_p_S34_reaclib = 978,
        k_n_Cl34_to_He4_P31_reaclib = 979,
        k_He4_Cl34_to_n_K37_reaclib = 980,
        k_He4_Cl34_to_p_Ar37_reaclib = 981,
        k_n_Cl35_to_p_S35_reaclib = 982,
        k_n_Cl35_to_He4_P32_reaclib = 983,
        k_p_Cl35_to_He4_S32_reaclib = 984,
        k_He4_Cl35_to_n_K38_reaclib = 985,
        k_He4_Cl35_to_p_Ar38_reaclib = 986,
        k_n_Cl36_to_p_S36_reaclib = 987,
        k_n_Cl36_to_He4_P33_reaclib = 988,
        k_p_Cl36_to_n_Ar36_reaclib = 989,
        k_p_Cl36_to_He4_S33_reaclib = 990,
        k_He4_Cl36_to_n_K39_reaclib = 991,
        k_He4_Cl36_to_p_Ar39_reaclib = 992,
        k_p_Cl37_to_n_Ar37_reaclib = 993,
        k_p_Cl37_to_He4_S34_reaclib = 994,
        k_He4_Cl37_to_n_K40_reaclib = 995,
        k_He4_Cl37_to_p_Ar40_reaclib = 996,
        k_n_Ar36_to_p_Cl36_reaclib = 997,
        k_n_Ar36_to_He4_S33_reaclib = 998,
        k_p_Ar36_to_He4_Cl33_reaclib = 999,
        k_He4_Ar36_to_p_K39_reaclib = 1000,
        k_n_Ar37_to_p_Cl37_reaclib = 1001,
        k_n_Ar37_to_He4_S34_reaclib = 1002,
        k_p_Ar37_to_n_K37_reaclib = 1003,
        k_p_Ar37_to_He4_Cl34_reaclib = 1004,
        k_He4_Ar37_to_n_Ca40_reaclib = 1005,
        k_He4_Ar37_to_p_K40_reaclib = 1006,
        k_n_Ar38_to_He4_S35_reaclib = 1007,
        k_p_Ar38_to_n_K38_reaclib = 1008,
        k_p_Ar38_to_He4_Cl35_reaclib = 1009,
        k_He4_Ar38_to_n_Ca41_reaclib = 1010,
        k_He4_Ar38_to_p_K41_reaclib = 1011,
        k_n_Ar39_to_He4_S36_reaclib = 1012,
        k_p_Ar39_to_n_K39_reaclib = 1013,
        k_p_Ar39_to_He4_Cl36_reaclib = 1014,
        k_He4_Ar39_to_n_Ca42_reaclib = 1015,
        k_p_Ar40_to_n_K40_reaclib = 1016,
        k_p_Ar40_to_He4_Cl37_reaclib = 1017,
        k_He4_Ar40_to_n_Ca43_reaclib = 1018,
        k_n_K37_to_p_Ar37_reaclib = 1019,
        k_n_K37_to_He4_Cl34_reaclib = 1020,
        k_He4_K37_to_p_Ca40_reaclib = 1021,
        k_n_K38_to_p_Ar38_reaclib = 1022,
        k_n_K38_to_He4_Cl35_reaclib = 1023,
        k_He4_K38_to_p_Ca41_reaclib = 1024,
        k_n_K39_to_p_Ar39_reaclib = 1025,
        k_n_K39_to_He4_Cl36_reaclib = 1026,
        k_p_K39_to_He4_Ar36_reaclib = 1027,
        k_He4_K39_to_p_Ca42_reaclib = 1028,
        k_n_K40_to_p_Ar40_reaclib = 1029,
        k_n_K40_to_He4_Cl37_reaclib = 1030,
        k_p_K40_to_n_Ca40_reaclib = 1031,
        k_p_K40_to_He4_Ar37_reaclib = 1032,
        k_He4_K40_to_n_Sc43_reaclib = 1033,
        k_He4_K40_to_p_Ca43_reaclib = 1034,
        k_p_K41_to_n_Ca41_reaclib = 1035,
        k_p_K41_to_He4_Ar38_reaclib = 1036,
        k_He4_K41_to_n_Sc44_reaclib = 1037,
        k_He4_K41_to_p_Ca44_reaclib = 1038,
        k_n_Ca40_to_p_K40_reaclib = 1039,
        k_n_Ca40_to_He4_Ar37_reaclib = 1040,
        k_p_Ca40_to_He4_K37_reaclib = 1041,
        k_He4_Ca40_to_p_Sc43_reaclib = 1042,
        k_n_Ca41_to_p_K41_reaclib = 1043,
        k_n_Ca41_to_He4_Ar38_reaclib = 1044,
        k_p_Ca41_to_He4_K38_reaclib = 1045,
        k_He4_Ca41_to_n_Ti44_reaclib = 1046,
        k_He4_Ca41_to_p_Sc44_reaclib = 1047,
        k_n_Ca42_to_He4_Ar39_reaclib = 1048,
        k_p_Ca42_to_He4_K39_reaclib = 1049,
        k_He4_Ca42_to_n_Ti45_reaclib = 1050,
        k_He4_Ca42_to_p_Sc45_reaclib = 1051,
        k_n_Ca43_to_He4_Ar40_reaclib = 1052,
        k_p_Ca43_to_n_Sc43_reaclib = 1053,
        k_p_Ca43_to_He4_K40_reaclib = 1054,
        k_He4_Ca43_to_n_Ti46_reaclib = 1055,
        k_He4_Ca43_to_p_Sc46_reaclib = 1056,
        k_p_Ca44_to_n_Sc44_reaclib = 1057,
        k_p_Ca44_to_He4_K41_reaclib = 1058,
        k_He4_Ca44_to_n_Ti47_reaclib = 1059,
        k_He4_Ca44_to_p_Sc47_reaclib = 1060,
        k_p_Ca45_to_n_Sc45_reaclib = 1061,
        k_He4_Ca45_to_n_Ti48_reaclib = 1062,
        k_He4_Ca45_to_p_Sc48_reaclib = 1063,
        k_p_Ca46_to_n_Sc46_reaclib = 1064,
        k_He4_Ca46_to_n_Ti49_reaclib = 1065,
        k_He4_Ca46_to_p_Sc49_reaclib = 1066,
        k_p_Ca47_to_n_Sc47_reaclib = 1067,
        k_He4_Ca47_to_n_Ti50_reaclib = 1068,
        k_p_Ca48_to_n_Sc48_reaclib = 1069,
        k_He4_Ca48_to_n_Ti51_reaclib = 1070,
        k_n_Sc43_to_p_Ca43_reaclib = 1071,
        k_n_Sc43_to_He4_K40_reaclib = 1072,
        k_p_Sc43_to_He4_Ca40_reaclib = 1073,
        k_He4_Sc43_to_n_V46_reaclib = 1074,
        k_He4_Sc43_to_p_Ti46_reaclib = 1075,
        k_n_Sc44_to_p_Ca44_reaclib = 1076,
        k_n_Sc44_to_He4_K41_reaclib = 1077,
        k_p_Sc44_to_n_Ti44_reaclib = 1078,
        k_p_Sc44_to_He4_Ca41_reaclib = 1079,
        k_He4_Sc44_to_n_V47_reaclib = 1080,
        k_He4_Sc44_to_p_Ti47_reaclib = 1081,
        k_n_Sc45_to_p_Ca45_reaclib = 1082,
        k_p_Sc45_to_n_Ti45_reaclib = 1083,
        k_p_Sc45_to_He4_Ca42_reaclib = 1084,
        k_He4_Sc45_to_n_V48_reaclib = 1085,
        k_He4_Sc45_to_p_Ti48_reaclib = 1086,
        k_n_Sc46_to_p_Ca46_reaclib = 1087,
        k_p_Sc46_to_n_Ti46_reaclib = 1088,
        k_p_Sc46_to_He4_Ca43_reaclib = 1089,
        k_He4_Sc46_to_n_V49_reaclib = 1090,
        k_He4_Sc46_to_p_Ti49_reaclib = 1091,
        k_n_Sc47_to_p_Ca47_reaclib = 1092,
        k_p_Sc47_to_n_Ti47_reaclib = 1093,
        k_p_Sc47_to_He4_Ca44_reaclib = 1094,
        k_He4_Sc47_to_n_V50_reaclib = 1095,
        k_He4_Sc47_to_p_Ti50_reaclib = 1096,
        k_n_Sc48_to_p_Ca48_reaclib = 1097,
        k_p_Sc48_to_n_Ti48_reaclib = 1098,
        k_p_Sc48_to_He4_Ca45_reaclib = 1099,
        k_He4_Sc48_to_n_V51_reaclib = 1100,
        k_He4_Sc48_to_p_Ti51_reaclib = 1101,
        k_p_Sc49_to_n_Ti49_reaclib = 1102,
        k_p_Sc49_to_He4_Ca46_reaclib = 1103,
        k_He4_Sc49_to_n_V52_reaclib = 1104,
        k_n_Ti44_to_p_Sc44_reaclib = 1105,
        k_n_Ti44_to_He4_Ca41_reaclib = 1106,
        k_He4_Ti44_to_p_V47_reaclib = 1107,
        k_n_Ti45_to_p_Sc45_reaclib = 1108,
        k_n_Ti45_to_He4_Ca42_reaclib = 1109,
        k_He4_Ti45_to_n_Cr48_reaclib = 1110,
        k_He4_Ti45_to_p_V48_reaclib = 1111,
        k_n_Ti46_to_p_Sc46_reaclib = 1112,
        k_n_Ti46_to_He4_Ca43_reaclib = 1113,
        k_p_Ti46_to_n_V46_reaclib = 1114,
        k_p_Ti46_to_He4_Sc43_reaclib = 1115,
        k_He4_Ti46_to_n_Cr49_reaclib = 1116,
        k_He4_Ti46_to_p_V49_reaclib = 1117,
        k_n_Ti47_to_p_Sc47_reaclib = 1118,
        k_n_Ti47_to_He4_Ca44_reaclib = 1119,
        k_p_Ti47_to_n_V47_reaclib = 1120,
        k_p_Ti47_to_He4_Sc44_reaclib = 1121,
        k_He4_Ti47_to_n_Cr50_reaclib = 1122,
        k_He4_Ti47_to_p_V50_reaclib = 1123,
        k_n_Ti48_to_p_Sc48_reaclib = 1124,
        k_n_Ti48_to_He4_Ca45_reaclib = 1125,
        k_p_Ti48_to_n_V48_reaclib = 1126,
        k_p_Ti48_to_He4_Sc45_reaclib = 1127,
        k_He4_Ti48_to_n_Cr51_reaclib = 1128,
        k_He4_Ti48_to_p_V51_reaclib = 1129,
        k_n_Ti49_to_p_Sc49_reaclib = 1130,
        k_n_Ti49_to_He4_Ca46_reaclib = 1131,
        k_p_Ti49_to_n_V49_reaclib = 1132,
        k_p_Ti49_to_He4_Sc46_reaclib = 1133,
        k_He4_Ti49_to_n_Cr52_reaclib = 1134,
        k_He4_Ti49_to_p_V52_reaclib = 1135,
        k_n_Ti50_to_He4_Ca47_reaclib = 1136,
        k_p_Ti50_to_n_V50_reaclib = 1137,
        k_p_Ti50_to_He4_Sc47_reaclib = 1138,
        k_He4_Ti50_to_n_Cr53_reaclib = 1139,
        k_n_Ti51_to_He4_Ca48_reaclib = 1140,
        k_p_Ti51_to_n_V51_reaclib = 1141,
        k_p_Ti51_to_He4_Sc48_reaclib = 1142,
        k_He4_Ti51_to_n_Cr54_reaclib = 1143,
        k_n_V46_to_p_Ti46_reaclib = 1144,
        k_n_V46_to_He4_Sc43_reaclib = 1145,
        k_He4_V46_to_p_Cr49_reaclib = 1146,
        k_n_V47_to_p_Ti47_reaclib = 1147,
        k_n_V47_to_He4_Sc44_reaclib = 1148,
        k_p_V47_to_He4_Ti44_reaclib = 1149,
        k_He4_V47_to_n_Mn50_reaclib = 1150,
        k_He4_V47_to_p_Cr50_reaclib = 1151,
        k_n_V48_to_p_Ti48_reaclib = 1152,
        k_n_V48_to_He4_Sc45_reaclib = 1153,
        k_p_V48_to_n_Cr48_reaclib = 1154,
        k_p_V48_to_He4_Ti45_reaclib = 1155,
        k_He4_V48_to_n_Mn51_reaclib = 1156,
        k_He4_V48_to_p_Cr51_reaclib = 1157,
        k_n_V49_to_p_Ti49_reaclib = 1158,
        k_n_V49_to_He4_Sc46_reaclib = 1159,
        k_p_V49_to_n_Cr49_reaclib = 1160,
        k_p_V49_to_He4_Ti46_reaclib = 1161,
        k_He4_V49_to_n_Mn52_reaclib = 1162,
        k_He4_V49_to_p_Cr52_reaclib = 1163,
        k_n_V50_to_p_Ti50_reaclib = 1164,
        k_n_V50_to_He4_Sc47_reaclib = 1165,
        k_p_V50_to_n_Cr50_reaclib = 1166,
        k_p_V50_to_He4_Ti47_reaclib = 1167,
        k_He4_V50_to_n_Mn53_reaclib = 1168,
        k_He4_V50_to_p_Cr53_reaclib = 1169,
        k_n_V51_to_p_Ti51_reaclib = 1170,
        k_n_V51_to_He4_Sc48_reaclib = 1171,
        k_p_V51_to_n_Cr51_reaclib = 1172,
        k_p_V51_to_He4_Ti48_reaclib = 1173,
        k_He4_V51_to_n_Mn54_reaclib = 1174,
        k_He4_V51_to_p_Cr54_reaclib = 1175,
        k_n_V52_to_He4_Sc49_reaclib = 1176,
        k_p_V52_to_n_Cr52_reaclib = 1177,
        k_p_V52_to_He4_Ti49_reaclib = 1178,
        k_He4_V52_to_n_Mn55_reaclib = 1179,
        k_n_Cr48_to_p_V48_reaclib = 1180,
        k_n_Cr48_to_He4_Ti45_reaclib = 1181,
        k_He4_Cr48_to_p_Mn51_reaclib = 1182,
        k_n_Cr49_to_p_V49_reaclib = 1183,
        k_n_Cr49_to_He4_Ti46_reaclib = 1184,
        k_p_Cr49_to_He4_V46_reaclib = 1185,
        k_He4_Cr49_to_n_Fe52_reaclib = 1186,
        k_He4_Cr49_to_p_Mn52_reaclib = 1187,
        k_n_Cr50_to_p_V50_reaclib = 1188,
        k_n_Cr50_to_He4_Ti47_reaclib = 1189,
        k_p_Cr50_to_n_Mn50_reaclib = 1190,
        k_p_Cr50_to_He4_V47_reaclib = 1191,
        k_He4_Cr50_to_n_Fe53_reaclib = 1192,
        k_He4_Cr50_to_p_Mn53_reaclib = 1193,
        k_n_Cr51_to_p_V51_reaclib = 1194,
        k_n_Cr51_to_He4_Ti48_reaclib = 1195,
        k_p_Cr51_to_n_Mn51_reaclib = 1196,
        k_p_Cr51_to_He4_V48_reaclib = 1197,
        k_He4_Cr51_to_n_Fe54_reaclib = 1198,
        k_He4_Cr51_to_p_Mn54_reaclib = 1199,
        k_n_Cr52_to_p_V52_reaclib = 1200,
        k_n_Cr52_to_He4_Ti49_reaclib = 1201,
        k_p_Cr52_to_n_Mn52_reaclib = 1202,
        k_p_Cr52_to_He4_V49_reaclib = 1203,
        k_He4_Cr52_to_n_Fe55_reaclib = 1204,
        k_He4_Cr52_to_p_Mn55_reaclib = 1205,
        k_n_Cr53_to_He4_Ti50_reaclib = 1206,
        k_p_Cr53_to_n_Mn53_reaclib = 1207,
        k_p_Cr53_to_He4_V50_reaclib = 1208,
        k_He4_Cr53_to_n_Fe56_reaclib = 1209,
        k_n_Cr54_to_He4_Ti51_reaclib = 1210,
        k_p_Cr54_to_n_Mn54_reaclib = 1211,
        k_p_Cr54_to_He4_V51_reaclib = 1212,
        k_He4_Cr54_to_n_Fe57_reaclib = 1213,
        k_n_Mn50_to_p_Cr50_reaclib = 1214,
        k_n_Mn50_to_He4_V47_reaclib = 1215,
        k_He4_Mn50_to_n_Co53_reaclib = 1216,
        k_He4_Mn50_to_p_Fe53_reaclib = 1217,
        k_n_Mn51_to_p_Cr51_reaclib = 1218,
        k_n_Mn51_to_He4_V48_reaclib = 1219,
        k_p_Mn51_to_He4_Cr48_reaclib = 1220,
        k_He4_Mn51_to_n_Co54_reaclib = 1221,
        k_He4_Mn51_to_p_Fe54_reaclib = 1222,
        k_n_Mn52_to_p_Cr52_reaclib = 1223,
        k_n_Mn52_to_He4_V49_reaclib = 1224,
        k_p_Mn52_to_n_Fe52_reaclib = 1225,
        k_p_Mn52_to_He4_Cr49_reaclib = 1226,
        k_He4_Mn52_to_n_Co55_reaclib = 1227,
        k_He4_Mn52_to_p_Fe55_reaclib = 1228,
        k_n_Mn53_to_p_Cr53_reaclib = 1229,
        k_n_Mn53_to_He4_V50_reaclib = 1230,
        k_p_Mn53_to_n_Fe53_reaclib = 1231,
        k_p_Mn53_to_He4_Cr50_reaclib = 1232,
        k_He4_Mn53_to_n_Co56_reaclib = 1233,
        k_He4_Mn53_to_p_Fe56_reaclib = 1234,
        k_n_Mn54_to_p_Cr54_reaclib = 1235,
        k_n_Mn54_to_He4_V51_reaclib = 1236,
        k_p_Mn54_to_n_Fe54_reaclib = 1237,
        k_p_Mn54_to_He4_Cr51_reaclib = 1238,
        k_He4_Mn54_to_n_Co57_reaclib = 1239,
        k_He4_Mn54_to_p_Fe57_reaclib = 1240,
        k_n_Mn55_to_He4_V52_reaclib = 1241,
        k_p_Mn55_to_n_Fe55_reaclib = 1242,
        k_p_Mn55_to_He4_Cr52_reaclib = 1243,
        k_He4_Mn55_to_n_Co58_reaclib = 1244,
        k_He4_Mn55_to_p_Fe58_reaclib = 1245,
        k_n_Fe52_to_p_Mn52_reaclib = 1246,
        k_n_Fe52_to_He4_Cr49_reaclib = 1247,
        k_He4_Fe52_to_p_Co55_reaclib = 1248,
        k_n_Fe53_to_p_Mn53_reaclib = 1249,
        k_n_Fe53_to_He4_Cr50_reaclib = 1250,
        k_p_Fe53_to_n_Co53_reaclib = 1251,
        k_p_Fe53_to_He4_Mn50_reaclib = 1252,
        k_He4_Fe53_to_n_Ni56_reaclib = 1253,
        k_He4_Fe53_to_p_Co56_reaclib = 1254,
        k_n_Fe54_to_p_Mn54_reaclib = 1255,
        k_n_Fe54_to_He4_Cr51_reaclib = 1256,
        k_p_Fe54_to_n_Co54_reaclib = 1257,
        k_p_Fe54_to_He4_Mn51_reaclib = 1258,
        k_He4_Fe54_to_n_Ni57_reaclib = 1259,
        k_He4_Fe54_to_p_Co57_reaclib = 1260,
        k_n_Fe55_to_p_Mn55_reaclib = 1261,
        k_n_Fe55_to_He4_Cr52_reaclib = 1262,
        k_p_Fe55_to_n_Co55_reaclib = 1263,
        k_p_Fe55_to_He4_Mn52_reaclib = 1264,
        k_He4_Fe55_to_n_Ni58_reaclib = 1265,
        k_He4_Fe55_to_p_Co58_reaclib = 1266,
        k_n_Fe56_to_He4_Cr53_reaclib = 1267,
        k_p_Fe56_to_n_Co56_reaclib = 1268,
        k_p_Fe56_to_He4_Mn53_reaclib = 1269,
        k_He4_Fe56_to_n_Ni59_reaclib = 1270,
        k_He4_Fe56_to_p_Co59_reaclib = 1271,
        k_n_Fe57_to_He4_Cr54_reaclib = 1272,
        k_p_Fe57_to_n_Co57_reaclib = 1273,
        k_p_Fe57_to_He4_Mn54_reaclib = 1274,
        k_He4_Fe57_to_n_Ni60_reaclib = 1275,
        k_p_Fe58_to_n_Co58_reaclib = 1276,
        k_p_Fe58_to_He4_Mn55_reaclib = 1277,
        k_He4_Fe58_to_n_Ni61_reaclib = 1278,
        k_n_Co53_to_p_Fe53_reaclib = 1279,
        k_n_Co53_to_He4_Mn50_reaclib = 1280,
        k_He4_Co53_to_p_Ni56_reaclib = 1281,
        k_n_Co54_to_p_Fe54_reaclib = 1282,
        k_n_Co54_to_He4_Mn51_reaclib = 1283,
        k_He4_Co54_to_n_Cu57_reaclib = 1284,
        k_He4_Co54_to_p_Ni57_reaclib = 1285,
        k_n_Co55_to_p_Fe55_reaclib = 1286,
        k_n_Co55_to_He4_Mn52_reaclib = 1287,
        k_p_Co55_to_He4_Fe52_reaclib = 1288,
        k_He4_Co55_to_n_Cu58_reaclib = 1289,
        k_He4_Co55_to_p_Ni58_reaclib = 1290,
        k_n_Co56_to_p_Fe56_reaclib = 1291,
        k_n_Co56_to_He4_Mn53_reaclib = 1292,
        k_p_Co56_to_n_Ni56_reaclib = 1293,
        k_p_Co56_to_He4_Fe53_reaclib = 1294,
        k_He4_Co56_to_n_Cu59_reaclib = 1295,
        k_He4_Co56_to_p_Ni59_reaclib = 1296,
        k_n_Co57_to_p_Fe57_reaclib = 1297,
        k_n_Co57_to_He4_Mn54_reaclib = 1298,
        k_p_Co57_to_n_Ni57_reaclib = 1299,
        k_p_Co57_to_He4_Fe54_reaclib = 1300,
        k_He4_Co57_to_n_Cu60_reaclib = 1301,
        k_He4_Co57_to_p_Ni60_reaclib = 1302,
        k_n_Co58_to_p_Fe58_reaclib = 1303,
        k_n_Co58_to_He4_Mn55_reaclib = 1304,
        k_p_Co58_to_n_Ni58_reaclib = 1305,
        k_p_Co58_to_He4_Fe55_reaclib = 1306,
        k_He4_Co58_to_n_Cu61_reaclib = 1307,
        k_He4_Co58_to_p_Ni61_reaclib = 1308,
        k_p_Co59_to_n_Ni59_reaclib = 1309,
        k_p_Co59_to_He4_Fe56_reaclib = 1310,
        k_He4_Co59_to_n_Cu62_reaclib = 1311,
        k_He4_Co59_to_p_Ni62_reaclib = 1312,
        k_n_Ni56_to_p_Co56_reaclib = 1313,
        k_n_Ni56_to_He4_Fe53_reaclib = 1314,
        k_p_Ni56_to_He4_Co53_reaclib = 1315,
        k_He4_Ni56_to_n_Zn59_reaclib = 1316,
        k_He4_Ni56_to_p_Cu59_reaclib = 1317,
        k_n_Ni57_to_p_Co57_reaclib = 1318,
        k_n_Ni57_to_He4_Fe54_reaclib = 1319,
        k_p_Ni57_to_n_Cu57_reaclib = 1320,
        k_p_Ni57_to_He4_Co54_reaclib = 1321,
        k_He4_Ni57_to_n_Zn60_reaclib = 1322,
        k_He4_Ni57_to_p_Cu60_reaclib = 1323,
        k_n_Ni58_to_p_Co58_reaclib = 1324,
        k_n_Ni58_to_He4_Fe55_reaclib = 1325,
        k_p_Ni58_to_n_Cu58_reaclib = 1326,
        k_p_Ni58_to_He4_Co55_reaclib = 1327,
        k_He4_Ni58_to_n_Zn61_reaclib = 1328,
        k_He4_Ni58_to_p_Cu61_reaclib = 1329,
        k_n_Ni59_to_p_Co59_reaclib = 1330,
        k_n_Ni59_to_He4_Fe56_reaclib = 1331,
        k_p_Ni59_to_n_Cu59_reaclib = 1332,
        k_p_Ni59_to_He4_Co56_reaclib = 1333,
        k_He4_Ni59_to_n_Zn62_reaclib = 1334,
        k_He4_Ni59_to_p_Cu62_reaclib = 1335,
        k_n_Ni60_to_He4_Fe57_reaclib = 1336,
        k_p_Ni60_to_n_Cu60_reaclib = 1337,
        k_p_Ni60_to_He4_Co57_reaclib = 1338,
        k_He4_Ni60_to_n_Zn63_reaclib = 1339,
        k_He4_Ni60_to_p_Cu63_reaclib = 1340,
        k_n_Ni61_to_He4_Fe58_reaclib = 1341,
        k_p_Ni61_to_n_Cu61_reaclib = 1342,
        k_p_Ni61_to_He4_Co58_reaclib = 1343,
        k_He4_Ni61_to_n_Zn64_reaclib = 1344,
        k_He4_Ni61_to_p_Cu64_reaclib = 1345,
        k_p_Ni62_to_n_Cu62_reaclib = 1346,
        k_p_Ni62_to_He4_Co59_reaclib = 1347,
        k_He4_Ni62_to_n_Zn65_reaclib = 1348,
        k_He4_Ni62_to_p_Cu65_reaclib = 1349,
        k_p_Ni63_to_n_Cu63_reaclib = 1350,
        k_He4_Ni63_to_n_Zn66_reaclib = 1351,
        k_p_Ni64_to_n_Cu64_reaclib = 1352,
        k_n_Cu57_to_p_Ni57_reaclib = 1353,
        k_n_Cu57_to_He4_Co54_reaclib = 1354,
        k_He4_Cu57_to_p_Zn60_reaclib = 1355,
        k_n_Cu58_to_p_Ni58_reaclib = 1356,
        k_n_Cu58_to_He4_Co55_reaclib = 1357,
        k_He4_Cu58_to_p_Zn61_reaclib = 1358,
        k_n_Cu59_to_p_Ni59_reaclib = 1359,
        k_n_Cu59_to_He4_Co56_reaclib = 1360,
        k_p_Cu59_to_n_Zn59_reaclib = 1361,
        k_p_Cu59_to_He4_Ni56_reaclib = 1362,
        k_He4_Cu59_to_n_Ga62_reaclib = 1363,
        k_He4_Cu59_to_p_Zn62_reaclib = 1364,
        k_n_Cu60_to_p_Ni60_reaclib = 1365,
        k_n_Cu60_to_He4_Co57_reaclib = 1366,
        k_p_Cu60_to_n_Zn60_reaclib = 1367,
        k_p_Cu60_to_He4_Ni57_reaclib = 1368,
        k_He4_Cu60_to_n_Ga63_reaclib = 1369,
        k_He4_Cu60_to_p_Zn63_reaclib = 1370,
        k_n_Cu61_to_p_Ni61_reaclib = 1371,
        k_n_Cu61_to_He4_Co58_reaclib = 1372,
        k_p_Cu61_to_n_Zn61_reaclib = 1373,
        k_p_Cu61_to_He4_Ni58_reaclib = 1374,
        k_He4_Cu61_to_n_Ga64_reaclib = 1375,
        k_He4_Cu61_to_p_Zn64_reaclib = 1376,
        k_n_Cu62_to_p_Ni62_reaclib = 1377,
        k_n_Cu62_to_He4_Co59_reaclib = 1378,
        k_p_Cu62_to_n_Zn62_reaclib = 1379,
        k_p_Cu62_to_He4_Ni59_reaclib = 1380,
        k_He4_Cu62_to_p_Zn65_reaclib = 1381,
        k_n_Cu63_to_p_Ni63_reaclib = 1382,
        k_p_Cu63_to_n_Zn63_reaclib = 1383,
        k_p_Cu63_to_He4_Ni60_reaclib = 1384,
        k_He4_Cu63_to_p_Zn66_reaclib = 1385,
        k_n_Cu64_to_p_Ni64_reaclib = 1386,
        k_p_Cu64_to_n_Zn64_reaclib = 1387,
        k_p_Cu64_to_He4_Ni61_reaclib = 1388,
        k_p_Cu65_to_n_Zn65_reaclib = 1389,
        k_p_Cu65_to_He4_Ni62_reaclib = 1390,
        k_n_Zn59_to_p_Cu59_reaclib = 1391,
        k_n_Zn59_to_He4_Ni56_reaclib = 1392,
        k_He4_Zn59_to_p_Ga62_reaclib = 1393,
        k_n_Zn60_to_p_Cu60_reaclib = 1394,
        k_n_Zn60_to_He4_Ni57_reaclib = 1395,
        k_p_Zn60_to_He4_Cu57_reaclib = 1396,
        k_He4_Zn60_to_n_Ge63_reaclib = 1397,
        k_He4_Zn60_to_p_Ga63_reaclib = 1398,
        k_n_Zn61_to_p_Cu61_reaclib = 1399,
        k_n_Zn61_to_He4_Ni58_reaclib = 1400,
        k_p_Zn61_to_He4_Cu58_reaclib = 1401,
        k_He4_Zn61_to_n_Ge64_reaclib = 1402,
        k_He4_Zn61_to_p_Ga64_reaclib = 1403,
        k_n_Zn62_to_p_Cu62_reaclib = 1404,
        k_n_Zn62_to_He4_Ni59_reaclib = 1405,
        k_p_Zn62_to_n_Ga62_reaclib = 1406,
        k_p_Zn62_to_He4_Cu59_reaclib = 1407,
        k_n_Zn63_to_p_Cu63_reaclib = 1408,
        k_n_Zn63_to_He4_Ni60_reaclib = 1409,
        k_p_Zn63_to_n_Ga63_reaclib = 1410,
        k_p_Zn63_to_He4_Cu60_reaclib = 1411,
        k_n_Zn64_to_p_Cu64_reaclib = 1412,
        k_n_Zn64_to_He4_Ni61_reaclib = 1413,
        k_p_Zn64_to_n_Ga64_reaclib = 1414,
        k_p_Zn64_to_He4_Cu61_reaclib = 1415,
        k_n_Zn65_to_p_Cu65_reaclib = 1416,
        k_n_Zn65_to_He4_Ni62_reaclib = 1417,
        k_p_Zn65_to_He4_Cu62_reaclib = 1418,
        k_n_Zn66_to_He4_Ni63_reaclib = 1419,
        k_p_Zn66_to_He4_Cu63_reaclib = 1420,
        k_n_Ga62_to_p_Zn62_reaclib = 1421,
        k_n_Ga62_to_He4_Cu59_reaclib = 1422,
        k_p_Ga62_to_He4_Zn59_reaclib = 1423,
        k_n_Ga63_to_p_Zn63_reaclib = 1424,
        k_n_Ga63_to_He4_Cu60_reaclib = 1425,
        k_p_Ga63_to_n_Ge63_reaclib = 1426,
        k_p_Ga63_to_He4_Zn60_reaclib = 1427,
        k_n_Ga64_to_p_Zn64_reaclib = 1428,
        k_n_Ga64_to_He4_Cu61_reaclib = 1429,
        k_p_Ga64_to_n_Ge64_reaclib = 1430,
        k_p_Ga64_to_He4_Zn61_reaclib = 1431,
        k_n_Ge63_to_p_Ga63_reaclib = 1432,
        k_n_Ge63_to_He4_Zn60_reaclib = 1433,
        k_n_Ge64_to_p_Ga64_reaclib = 1434,
        k_n_Ge64_to_He4_Zn61_reaclib = 1435,
        k_p_d_to_n_p_p_reaclib = 1436,
        k_He3_He3_to_p_p_He4_reaclib = 1437,
        k_d_Li7_to_n_He4_He4_reaclib = 1438,
        k_d_Be7_to_p_He4_He4_reaclib = 1439,
        k_p_Be9_to_d_He4_He4_reaclib = 1440,
        k_n_B8_to_p_He4_He4_reaclib = 1441,
        k_p_B11_to_He4_He4_He4_reaclib = 1442,
        k_He3_Li7_to_n_p_He4_He4_reaclib = 1443,
        k_He3_Be7_to_p_p_He4_He4_reaclib = 1444,
        k_p_Be9_to_n_p_He4_He4_reaclib = 1445,
        k_n_p_He4_to_Li6_reaclib = 1446,
        k_n_He4_He4_to_Be9_reaclib = 1447,
        k_He4_He4_He4_to_C12_reaclib = 1448,
        k_n_p_p_to_p_d_reaclib = 1449,
        k_p_p_He4_to_He3_He3_reaclib = 1450,
        k_n_He4_He4_to_d_Li7_reaclib = 1451,
        k_p_He4_He4_to_n_B8_reaclib = 1452,
        k_p_He4_He4_to_d_Be7_reaclib = 1453,
        k_d_He4_He4_to_p_Be9_reaclib = 1454,
        k_He4_He4_He4_to_p_B11_reaclib = 1455,
        k_n_p_He4_He4_to_He3_Li7_reaclib = 1456,
        k_n_p_He4_He4_to_p_Be9_reaclib = 1457,
        k_p_p_He4_He4_to_He3_Be7_reaclib = 1458,
        k_Na21_to_Ne21_weaktab = 1459,
        k_Ne21_to_Na21_weaktab = 1460,
        k_Na22_to_Ne22_weaktab = 1461,
        k_Ne22_to_Na22_weaktab = 1462,
        k_Mg23_to_Na23_weaktab = 1463,
        k_Na23_to_Mg23_weaktab = 1464,
        k_Al25_to_Mg25_weaktab = 1465,
        k_Mg25_to_Al25_weaktab = 1466,
        k_Al26_to_Mg26_weaktab = 1467,
        k_Mg26_to_Al26_weaktab = 1468,
        k_P29_to_Si29_weaktab = 1469,
        k_Si29_to_P29_weaktab = 1470,
        k_P30_to_Si30_weaktab = 1471,
        k_Si30_to_P30_weaktab = 1472,
        k_P31_to_Si31_weaktab = 1473,
        k_Si31_to_P31_weaktab = 1474,
        k_P32_to_S32_weaktab = 1475,
        k_P32_to_Si32_weaktab = 1476,
        k_S32_to_P32_weaktab = 1477,
        k_Si32_to_P32_weaktab = 1478,
        k_Cl33_to_S33_weaktab = 1479,
        k_P33_to_S33_weaktab = 1480,
        k_S33_to_Cl33_weaktab = 1481,
        k_S33_to_P33_weaktab = 1482,
        k_Cl34_to_S34_weaktab = 1483,
        k_S34_to_Cl34_weaktab = 1484,
        k_Cl35_to_S35_weaktab = 1485,
        k_S35_to_Cl35_weaktab = 1486,
        k_Ar36_to_Cl36_weaktab = 1487,
        k_Cl36_to_Ar36_weaktab = 1488,
        k_Cl36_to_S36_weaktab = 1489,
        k_S36_to_Cl36_weaktab = 1490,
        k_Ar37_to_Cl37_weaktab = 1491,
        k_Ar37_to_K37_weaktab = 1492,
        k_Cl37_to_Ar37_weaktab = 1493,
        k_K37_to_Ar37_weaktab = 1494,
        k_Ar38_to_K38_weaktab = 1495,
        k_K38_to_Ar38_weaktab = 1496,
        k_Ar39_to_K39_weaktab = 1497,
        k_K39_to_Ar39_weaktab = 1498,
        k_Ar40_to_K40_weaktab = 1499,
        k_Ca40_to_K40_weaktab = 1500,
        k_K40_to_Ar40_weaktab = 1501,
        k_K40_to_Ca40_weaktab = 1502,
        k_Ca41_to_K41_weaktab = 1503,
        k_K41_to_Ca41_weaktab = 1504,
        k_Ca43_to_Sc43_weaktab = 1505,
        k_Sc43_to_Ca43_weaktab = 1506,
        k_Ca44_to_Sc44_weaktab = 1507,
        k_Sc44_to_Ca44_weaktab = 1508,
        k_Sc44_to_Ti44_weaktab = 1509,
        k_Ti44_to_Sc44_weaktab = 1510,
        k_Co53_to_Fe53_weaktab = 1511,
        k_Fe53_to_Co53_weaktab = 1512,
        k_Cu57_to_Ni57_weaktab = 1513,
        k_Ni57_to_Cu57_weaktab = 1514,
        k_Ca45_to_Sc45_weaktab = 1515,
        k_Sc45_to_Ca45_weaktab = 1516,
        k_Sc45_to_Ti45_weaktab = 1517,
        k_Ti45_to_Sc45_weaktab = 1518,
        k_Ca46_to_Sc46_weaktab = 1519,
        k_Sc46_to_Ca46_weaktab = 1520,
        k_Sc46_to_Ti46_weaktab = 1521,
        k_Ti46_to_Sc46_weaktab = 1522,
        k_Ti46_to_V46_weaktab = 1523,
        k_V46_to_Ti46_weaktab = 1524,
        k_Ca47_to_Sc47_weaktab = 1525,
        k_Sc47_to_Ca47_weaktab = 1526,
        k_Sc47_to_Ti47_weaktab = 1527,
        k_Ti47_to_Sc47_weaktab = 1528,
        k_Ti47_to_V47_weaktab = 1529,
        k_V47_to_Ti47_weaktab = 1530,
        k_Ca48_to_Sc48_weaktab = 1531,
        k_Cr48_to_V48_weaktab = 1532,
        k_Sc48_to_Ca48_weaktab = 1533,
        k_Sc48_to_Ti48_weaktab = 1534,
        k_Ti48_to_Sc48_weaktab = 1535,
        k_Ti48_to_V48_weaktab = 1536,
        k_V48_to_Cr48_weaktab = 1537,
        k_V48_to_Ti48_weaktab = 1538,
        k_Cr49_to_V49_weaktab = 1539,
        k_Sc49_to_Ti49_weaktab = 1540,
        k_Ti49_to_Sc49_weaktab = 1541,
        k_Ti49_to_V49_weaktab = 1542,
        k_V49_to_Cr49_weaktab = 1543,
        k_V49_to_Ti49_weaktab = 1544,
        k_Cr50_to_Mn50_weaktab = 1545,
        k_Cr50_to_V50_weaktab = 1546,
        k_Mn50_to_Cr50_weaktab = 1547,
        k_Ti50_to_V50_weaktab = 1548,
        k_V50_to_Cr50_weaktab = 1549,
        k_V50_to_Ti50_weaktab = 1550,
        k_Cr51_to_Mn51_weaktab = 1551,
        k_Cr51_to_V51_weaktab = 1552,
        k_Mn51_to_Cr51_weaktab = 1553,
        k_Ti51_to_V51_weaktab = 1554,
        k_V51_to_Cr51_weaktab = 1555,
        k_V51_to_Ti51_weaktab = 1556,
        k_Cr52_to_Mn52_weaktab = 1557,
        k_Cr52_to_V52_weaktab = 1558,
        k_Fe52_to_Mn52_weaktab = 1559,
        k_Mn52_to_Cr52_weaktab = 1560,
        k_Mn52_to_Fe52_weaktab = 1561,
        k_V52_to_Cr52_weaktab = 1562,
        k_Cr53_to_Mn53_weaktab = 1563,
        k_Fe53_to_Mn53_weaktab = 1564,
        k_Mn53_to_Cr53_weaktab = 1565,
        k_Mn53_to_Fe53_weaktab = 1566,
        k_Co54_to_Fe54_weaktab = 1567,
        k_Cr54_to_Mn54_weaktab = 1568,
        k_Fe54_to_Co54_weaktab = 1569,
        k_Fe54_to_Mn54_weaktab = 1570,
        k_Mn54_to_Cr54_weaktab = 1571,
        k_Mn54_to_Fe54_weaktab = 1572,
        k_Co55_to_Fe55_weaktab = 1573,
        k_Fe55_to_Co55_weaktab = 1574,
        k_Fe55_to_Mn55_weaktab = 1575,
        k_Mn55_to_Fe55_weaktab = 1576,
        k_Co56_to_Fe56_weaktab = 1577,
        k_Co56_to_Ni56_weaktab = 1578,
        k_Fe56_to_Co56_weaktab = 1579,
        k_Ni56_to_Co56_weaktab = 1580,
        k_Co57_to_Fe57_weaktab = 1581,
        k_Co57_to_Ni57_weaktab = 1582,
        k_Fe57_to_Co57_weaktab = 1583,
        k_Ni57_to_Co57_weaktab = 1584,
        k_Co58_to_Fe58_weaktab = 1585,
        k_Co58_to_Ni58_weaktab = 1586,
        k_Cu58_to_Ni58_weaktab = 1587,
        k_Fe58_to_Co58_weaktab = 1588,
        k_Ni58_to_Co58_weaktab = 1589,
        k_Ni58_to_Cu58_weaktab = 1590,
        k_Co59_to_Ni59_weaktab = 1591,
        k_Cu59_to_Ni59_weaktab = 1592,
        k_Ni59_to_Co59_weaktab = 1593,
        k_Ni59_to_Cu59_weaktab = 1594,
        k_Cu60_to_Ni60_weaktab = 1595,
        k_Cu60_to_Zn60_weaktab = 1596,
        k_Ni60_to_Cu60_weaktab = 1597,
        k_Zn60_to_Cu60_weaktab = 1598,
        k_Cu61_to_Ni61_weaktab = 1599,
        k_Cu61_to_Zn61_weaktab = 1600,
        k_Ni61_to_Cu61_weaktab = 1601,
        k_Zn61_to_Cu61_weaktab = 1602,
        k_Cu62_to_Ni62_weaktab = 1603,
        k_Cu62_to_Zn62_weaktab = 1604,
        k_Ga62_to_Zn62_weaktab = 1605,
        k_Ni62_to_Cu62_weaktab = 1606,
        k_Zn62_to_Cu62_weaktab = 1607,
        k_Zn62_to_Ga62_weaktab = 1608,
        k_Cu63_to_Ni63_weaktab = 1609,
        k_Cu63_to_Zn63_weaktab = 1610,
        k_Ga63_to_Zn63_weaktab = 1611,
        k_Ni63_to_Cu63_weaktab = 1612,
        k_Zn63_to_Cu63_weaktab = 1613,
        k_Zn63_to_Ga63_weaktab = 1614,
        k_Cu64_to_Ni64_weaktab = 1615,
        k_Cu64_to_Zn64_weaktab = 1616,
        k_Ga64_to_Ge64_weaktab = 1617,
        k_Ga64_to_Zn64_weaktab = 1618,
        k_Ge64_to_Ga64_weaktab = 1619,
        k_Ni64_to_Cu64_weaktab = 1620,
        k_Zn64_to_Cu64_weaktab = 1621,
        k_Zn64_to_Ga64_weaktab = 1622,
        k_Cu65_to_Zn65_weaktab = 1623,
        k_Zn65_to_Cu65_weaktab = 1624,
        k_n_to_p_weaktab = 1625,
        k_p_to_n_weaktab = 1626,
        NumRates = k_p_to_n_weaktab
    };

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "Be7_to_Li7_reaclib",  // 1,
        "C14_to_N14_reaclib",  // 2,
        "N13_to_C13_reaclib",  // 3,
        "O14_to_N14_reaclib",  // 4,
        "O15_to_N15_reaclib",  // 5,
        "F17_to_O17_reaclib",  // 6,
        "F18_to_O18_reaclib",  // 7,
        "Ne18_to_F18_reaclib",  // 8,
        "Ne19_to_F19_reaclib",  // 9,
        "Zn59_to_Cu59_reaclib",  // 10,
        "Ge63_to_Ga63_reaclib",  // 11,
        "d_to_n_p_reaclib",  // 12,
        "He3_to_p_d_reaclib",  // 13,
        "He4_to_n_He3_reaclib",  // 14,
        "He4_to_d_d_reaclib",  // 15,
        "Li6_to_He4_d_reaclib",  // 16,
        "Li7_to_n_Li6_reaclib",  // 17,
        "Be7_to_p_Li6_reaclib",  // 18,
        "Be7_to_He4_He3_reaclib",  // 19,
        "B8_to_p_Be7_reaclib",  // 20,
        "B8_to_He4_He4_reaclib",  // 21,
        "B10_to_p_Be9_reaclib",  // 22,
        "B10_to_He4_Li6_reaclib",  // 23,
        "B11_to_n_B10_reaclib",  // 24,
        "B11_to_He4_Li7_reaclib",  // 25,
        "C12_to_p_B11_reaclib",  // 26,
        "C13_to_n_C12_reaclib",  // 27,
        "C14_to_n_C13_reaclib",  // 28,
        "N13_to_p_C12_reaclib",  // 29,
        "N14_to_n_N13_reaclib",  // 30,
        "N14_to_p_C13_reaclib",  // 31,
        "N15_to_n_N14_reaclib",  // 32,
        "N15_to_p_C14_reaclib",  // 33,
        "O14_to_p_N13_reaclib",  // 34,
        "O15_to_n_O14_reaclib",  // 35,
        "O15_to_p_N14_reaclib",  // 36,
        "O16_to_n_O15_reaclib",  // 37,
        "O16_to_p_N15_reaclib",  // 38,
        "O16_to_He4_C12_reaclib",  // 39,
        "O17_to_n_O16_reaclib",  // 40,
        "O18_to_n_O17_reaclib",  // 41,
        "O18_to_He4_C14_reaclib",  // 42,
        "F17_to_p_O16_reaclib",  // 43,
        "F18_to_n_F17_reaclib",  // 44,
        "F18_to_p_O17_reaclib",  // 45,
        "F18_to_He4_N14_reaclib",  // 46,
        "F19_to_n_F18_reaclib",  // 47,
        "F19_to_p_O18_reaclib",  // 48,
        "F19_to_He4_N15_reaclib",  // 49,
        "Ne18_to_p_F17_reaclib",  // 50,
        "Ne18_to_He4_O14_reaclib",  // 51,
        "Ne19_to_n_Ne18_reaclib",  // 52,
        "Ne19_to_p_F18_reaclib",  // 53,
        "Ne19_to_He4_O15_reaclib",  // 54,
        "Ne20_to_n_Ne19_reaclib",  // 55,
        "Ne20_to_p_F19_reaclib",  // 56,
        "Ne20_to_He4_O16_reaclib",  // 57,
        "Ne21_to_n_Ne20_reaclib",  // 58,
        "Ne21_to_He4_O17_reaclib",  // 59,
        "Ne22_to_n_Ne21_reaclib",  // 60,
        "Ne22_to_He4_O18_reaclib",  // 61,
        "Na21_to_p_Ne20_reaclib",  // 62,
        "Na21_to_He4_F17_reaclib",  // 63,
        "Na22_to_n_Na21_reaclib",  // 64,
        "Na22_to_p_Ne21_reaclib",  // 65,
        "Na22_to_He4_F18_reaclib",  // 66,
        "Na23_to_n_Na22_reaclib",  // 67,
        "Na23_to_p_Ne22_reaclib",  // 68,
        "Na23_to_He4_F19_reaclib",  // 69,
        "Mg23_to_p_Na22_reaclib",  // 70,
        "Mg23_to_He4_Ne19_reaclib",  // 71,
        "Mg24_to_n_Mg23_reaclib",  // 72,
        "Mg24_to_p_Na23_reaclib",  // 73,
        "Mg24_to_He4_Ne20_reaclib",  // 74,
        "Mg25_to_n_Mg24_reaclib",  // 75,
        "Mg25_to_He4_Ne21_reaclib",  // 76,
        "Mg26_to_n_Mg25_reaclib",  // 77,
        "Mg26_to_He4_Ne22_reaclib",  // 78,
        "Al25_to_p_Mg24_reaclib",  // 79,
        "Al25_to_He4_Na21_reaclib",  // 80,
        "Al26_to_n_Al25_reaclib",  // 81,
        "Al26_to_p_Mg25_reaclib",  // 82,
        "Al26_to_He4_Na22_reaclib",  // 83,
        "Al27_to_n_Al26_reaclib",  // 84,
        "Al27_to_p_Mg26_reaclib",  // 85,
        "Al27_to_He4_Na23_reaclib",  // 86,
        "Si28_to_p_Al27_reaclib",  // 87,
        "Si28_to_He4_Mg24_reaclib",  // 88,
        "Si29_to_n_Si28_reaclib",  // 89,
        "Si29_to_He4_Mg25_reaclib",  // 90,
        "Si30_to_n_Si29_reaclib",  // 91,
        "Si30_to_He4_Mg26_reaclib",  // 92,
        "Si31_to_n_Si30_reaclib",  // 93,
        "Si32_to_n_Si31_reaclib",  // 94,
        "P29_to_p_Si28_reaclib",  // 95,
        "P29_to_He4_Al25_reaclib",  // 96,
        "P30_to_n_P29_reaclib",  // 97,
        "P30_to_p_Si29_reaclib",  // 98,
        "P30_to_He4_Al26_reaclib",  // 99,
        "P31_to_n_P30_reaclib",  // 100,
        "P31_to_p_Si30_reaclib",  // 101,
        "P31_to_He4_Al27_reaclib",  // 102,
        "P32_to_n_P31_reaclib",  // 103,
        "P32_to_p_Si31_reaclib",  // 104,
        "P33_to_n_P32_reaclib",  // 105,
        "P33_to_p_Si32_reaclib",  // 106,
        "S32_to_p_P31_reaclib",  // 107,
        "S32_to_He4_Si28_reaclib",  // 108,
        "S33_to_n_S32_reaclib",  // 109,
        "S33_to_p_P32_reaclib",  // 110,
        "S33_to_He4_Si29_reaclib",  // 111,
        "S34_to_n_S33_reaclib",  // 112,
        "S34_to_p_P33_reaclib",  // 113,
        "S34_to_He4_Si30_reaclib",  // 114,
        "S35_to_n_S34_reaclib",  // 115,
        "S35_to_He4_Si31_reaclib",  // 116,
        "S36_to_n_S35_reaclib",  // 117,
        "S36_to_He4_Si32_reaclib",  // 118,
        "Cl33_to_p_S32_reaclib",  // 119,
        "Cl33_to_He4_P29_reaclib",  // 120,
        "Cl34_to_n_Cl33_reaclib",  // 121,
        "Cl34_to_p_S33_reaclib",  // 122,
        "Cl34_to_He4_P30_reaclib",  // 123,
        "Cl35_to_n_Cl34_reaclib",  // 124,
        "Cl35_to_p_S34_reaclib",  // 125,
        "Cl35_to_He4_P31_reaclib",  // 126,
        "Cl36_to_n_Cl35_reaclib",  // 127,
        "Cl36_to_p_S35_reaclib",  // 128,
        "Cl36_to_He4_P32_reaclib",  // 129,
        "Cl37_to_n_Cl36_reaclib",  // 130,
        "Cl37_to_p_S36_reaclib",  // 131,
        "Cl37_to_He4_P33_reaclib",  // 132,
        "Ar36_to_p_Cl35_reaclib",  // 133,
        "Ar36_to_He4_S32_reaclib",  // 134,
        "Ar37_to_n_Ar36_reaclib",  // 135,
        "Ar37_to_p_Cl36_reaclib",  // 136,
        "Ar37_to_He4_S33_reaclib",  // 137,
        "Ar38_to_n_Ar37_reaclib",  // 138,
        "Ar38_to_p_Cl37_reaclib",  // 139,
        "Ar38_to_He4_S34_reaclib",  // 140,
        "Ar39_to_n_Ar38_reaclib",  // 141,
        "Ar39_to_He4_S35_reaclib",  // 142,
        "Ar40_to_n_Ar39_reaclib",  // 143,
        "Ar40_to_He4_S36_reaclib",  // 144,
        "K37_to_p_Ar36_reaclib",  // 145,
        "K37_to_He4_Cl33_reaclib",  // 146,
        "K38_to_n_K37_reaclib",  // 147,
        "K38_to_p_Ar37_reaclib",  // 148,
        "K38_to_He4_Cl34_reaclib",  // 149,
        "K39_to_n_K38_reaclib",  // 150,
        "K39_to_p_Ar38_reaclib",  // 151,
        "K39_to_He4_Cl35_reaclib",  // 152,
        "K40_to_n_K39_reaclib",  // 153,
        "K40_to_p_Ar39_reaclib",  // 154,
        "K40_to_He4_Cl36_reaclib",  // 155,
        "K41_to_n_K40_reaclib",  // 156,
        "K41_to_p_Ar40_reaclib",  // 157,
        "K41_to_He4_Cl37_reaclib",  // 158,
        "Ca40_to_p_K39_reaclib",  // 159,
        "Ca40_to_He4_Ar36_reaclib",  // 160,
        "Ca41_to_n_Ca40_reaclib",  // 161,
        "Ca41_to_p_K40_reaclib",  // 162,
        "Ca41_to_He4_Ar37_reaclib",  // 163,
        "Ca42_to_n_Ca41_reaclib",  // 164,
        "Ca42_to_p_K41_reaclib",  // 165,
        "Ca42_to_He4_Ar38_reaclib",  // 166,
        "Ca43_to_n_Ca42_reaclib",  // 167,
        "Ca43_to_He4_Ar39_reaclib",  // 168,
        "Ca44_to_n_Ca43_reaclib",  // 169,
        "Ca44_to_He4_Ar40_reaclib",  // 170,
        "Ca45_to_n_Ca44_reaclib",  // 171,
        "Ca46_to_n_Ca45_reaclib",  // 172,
        "Ca47_to_n_Ca46_reaclib",  // 173,
        "Ca48_to_n_Ca47_reaclib",  // 174,
        "Sc43_to_p_Ca42_reaclib",  // 175,
        "Sc43_to_He4_K39_reaclib",  // 176,
        "Sc44_to_n_Sc43_reaclib",  // 177,
        "Sc44_to_p_Ca43_reaclib",  // 178,
        "Sc44_to_He4_K40_reaclib",  // 179,
        "Sc45_to_n_Sc44_reaclib",  // 180,
        "Sc45_to_p_Ca44_reaclib",  // 181,
        "Sc45_to_He4_K41_reaclib",  // 182,
        "Sc46_to_n_Sc45_reaclib",  // 183,
        "Sc46_to_p_Ca45_reaclib",  // 184,
        "Sc47_to_n_Sc46_reaclib",  // 185,
        "Sc47_to_p_Ca46_reaclib",  // 186,
        "Sc48_to_n_Sc47_reaclib",  // 187,
        "Sc48_to_p_Ca47_reaclib",  // 188,
        "Sc49_to_n_Sc48_reaclib",  // 189,
        "Sc49_to_p_Ca48_reaclib",  // 190,
        "Ti44_to_p_Sc43_reaclib",  // 191,
        "Ti44_to_He4_Ca40_reaclib",  // 192,
        "Ti45_to_n_Ti44_reaclib",  // 193,
        "Ti45_to_p_Sc44_reaclib",  // 194,
        "Ti45_to_He4_Ca41_reaclib",  // 195,
        "Ti46_to_n_Ti45_reaclib",  // 196,
        "Ti46_to_p_Sc45_reaclib",  // 197,
        "Ti46_to_He4_Ca42_reaclib",  // 198,
        "Ti47_to_n_Ti46_reaclib",  // 199,
        "Ti47_to_p_Sc46_reaclib",  // 200,
        "Ti47_to_He4_Ca43_reaclib",  // 201,
        "Ti48_to_n_Ti47_reaclib",  // 202,
        "Ti48_to_p_Sc47_reaclib",  // 203,
        "Ti48_to_He4_Ca44_reaclib",  // 204,
        "Ti49_to_n_Ti48_reaclib",  // 205,
        "Ti49_to_p_Sc48_reaclib",  // 206,
        "Ti49_to_He4_Ca45_reaclib",  // 207,
        "Ti50_to_n_Ti49_reaclib",  // 208,
        "Ti50_to_p_Sc49_reaclib",  // 209,
        "Ti50_to_He4_Ca46_reaclib",  // 210,
        "Ti51_to_n_Ti50_reaclib",  // 211,
        "Ti51_to_He4_Ca47_reaclib",  // 212,
        "V46_to_p_Ti45_reaclib",  // 213,
        "V47_to_n_V46_reaclib",  // 214,
        "V47_to_p_Ti46_reaclib",  // 215,
        "V47_to_He4_Sc43_reaclib",  // 216,
        "V48_to_n_V47_reaclib",  // 217,
        "V48_to_p_Ti47_reaclib",  // 218,
        "V48_to_He4_Sc44_reaclib",  // 219,
        "V49_to_n_V48_reaclib",  // 220,
        "V49_to_p_Ti48_reaclib",  // 221,
        "V49_to_He4_Sc45_reaclib",  // 222,
        "V50_to_n_V49_reaclib",  // 223,
        "V50_to_p_Ti49_reaclib",  // 224,
        "V50_to_He4_Sc46_reaclib",  // 225,
        "V51_to_n_V50_reaclib",  // 226,
        "V51_to_p_Ti50_reaclib",  // 227,
        "V51_to_He4_Sc47_reaclib",  // 228,
        "V52_to_n_V51_reaclib",  // 229,
        "V52_to_p_Ti51_reaclib",  // 230,
        "V52_to_He4_Sc48_reaclib",  // 231,
        "Cr48_to_p_V47_reaclib",  // 232,
        "Cr48_to_He4_Ti44_reaclib",  // 233,
        "Cr49_to_n_Cr48_reaclib",  // 234,
        "Cr49_to_p_V48_reaclib",  // 235,
        "Cr49_to_He4_Ti45_reaclib",  // 236,
        "Cr50_to_n_Cr49_reaclib",  // 237,
        "Cr50_to_p_V49_reaclib",  // 238,
        "Cr50_to_He4_Ti46_reaclib",  // 239,
        "Cr51_to_n_Cr50_reaclib",  // 240,
        "Cr51_to_p_V50_reaclib",  // 241,
        "Cr51_to_He4_Ti47_reaclib",  // 242,
        "Cr52_to_n_Cr51_reaclib",  // 243,
        "Cr52_to_p_V51_reaclib",  // 244,
        "Cr52_to_He4_Ti48_reaclib",  // 245,
        "Cr53_to_n_Cr52_reaclib",  // 246,
        "Cr53_to_p_V52_reaclib",  // 247,
        "Cr53_to_He4_Ti49_reaclib",  // 248,
        "Cr54_to_n_Cr53_reaclib",  // 249,
        "Cr54_to_He4_Ti50_reaclib",  // 250,
        "Mn50_to_p_Cr49_reaclib",  // 251,
        "Mn50_to_He4_V46_reaclib",  // 252,
        "Mn51_to_n_Mn50_reaclib",  // 253,
        "Mn51_to_p_Cr50_reaclib",  // 254,
        "Mn51_to_He4_V47_reaclib",  // 255,
        "Mn52_to_n_Mn51_reaclib",  // 256,
        "Mn52_to_p_Cr51_reaclib",  // 257,
        "Mn52_to_He4_V48_reaclib",  // 258,
        "Mn53_to_n_Mn52_reaclib",  // 259,
        "Mn53_to_p_Cr52_reaclib",  // 260,
        "Mn53_to_He4_V49_reaclib",  // 261,
        "Mn54_to_n_Mn53_reaclib",  // 262,
        "Mn54_to_p_Cr53_reaclib",  // 263,
        "Mn54_to_He4_V50_reaclib",  // 264,
        "Mn55_to_n_Mn54_reaclib",  // 265,
        "Mn55_to_p_Cr54_reaclib",  // 266,
        "Mn55_to_He4_V51_reaclib",  // 267,
        "Fe52_to_p_Mn51_reaclib",  // 268,
        "Fe52_to_He4_Cr48_reaclib",  // 269,
        "Fe53_to_n_Fe52_reaclib",  // 270,
        "Fe53_to_p_Mn52_reaclib",  // 271,
        "Fe53_to_He4_Cr49_reaclib",  // 272,
        "Fe54_to_n_Fe53_reaclib",  // 273,
        "Fe54_to_p_Mn53_reaclib",  // 274,
        "Fe54_to_He4_Cr50_reaclib",  // 275,
        "Fe55_to_n_Fe54_reaclib",  // 276,
        "Fe55_to_p_Mn54_reaclib",  // 277,
        "Fe55_to_He4_Cr51_reaclib",  // 278,
        "Fe56_to_n_Fe55_reaclib",  // 279,
        "Fe56_to_p_Mn55_reaclib",  // 280,
        "Fe56_to_He4_Cr52_reaclib",  // 281,
        "Fe57_to_n_Fe56_reaclib",  // 282,
        "Fe57_to_He4_Cr53_reaclib",  // 283,
        "Fe58_to_n_Fe57_reaclib",  // 284,
        "Fe58_to_He4_Cr54_reaclib",  // 285,
        "Co53_to_p_Fe52_reaclib",  // 286,
        "Co54_to_n_Co53_reaclib",  // 287,
        "Co54_to_p_Fe53_reaclib",  // 288,
        "Co54_to_He4_Mn50_reaclib",  // 289,
        "Co55_to_n_Co54_reaclib",  // 290,
        "Co55_to_p_Fe54_reaclib",  // 291,
        "Co55_to_He4_Mn51_reaclib",  // 292,
        "Co56_to_n_Co55_reaclib",  // 293,
        "Co56_to_p_Fe55_reaclib",  // 294,
        "Co56_to_He4_Mn52_reaclib",  // 295,
        "Co57_to_n_Co56_reaclib",  // 296,
        "Co57_to_p_Fe56_reaclib",  // 297,
        "Co57_to_He4_Mn53_reaclib",  // 298,
        "Co58_to_n_Co57_reaclib",  // 299,
        "Co58_to_p_Fe57_reaclib",  // 300,
        "Co58_to_He4_Mn54_reaclib",  // 301,
        "Co59_to_n_Co58_reaclib",  // 302,
        "Co59_to_p_Fe58_reaclib",  // 303,
        "Co59_to_He4_Mn55_reaclib",  // 304,
        "Ni56_to_p_Co55_reaclib",  // 305,
        "Ni56_to_He4_Fe52_reaclib",  // 306,
        "Ni57_to_n_Ni56_reaclib",  // 307,
        "Ni57_to_p_Co56_reaclib",  // 308,
        "Ni57_to_He4_Fe53_reaclib",  // 309,
        "Ni58_to_n_Ni57_reaclib",  // 310,
        "Ni58_to_p_Co57_reaclib",  // 311,
        "Ni58_to_He4_Fe54_reaclib",  // 312,
        "Ni59_to_n_Ni58_reaclib",  // 313,
        "Ni59_to_p_Co58_reaclib",  // 314,
        "Ni59_to_He4_Fe55_reaclib",  // 315,
        "Ni60_to_n_Ni59_reaclib",  // 316,
        "Ni60_to_p_Co59_reaclib",  // 317,
        "Ni60_to_He4_Fe56_reaclib",  // 318,
        "Ni61_to_n_Ni60_reaclib",  // 319,
        "Ni61_to_He4_Fe57_reaclib",  // 320,
        "Ni62_to_n_Ni61_reaclib",  // 321,
        "Ni62_to_He4_Fe58_reaclib",  // 322,
        "Ni63_to_n_Ni62_reaclib",  // 323,
        "Ni64_to_n_Ni63_reaclib",  // 324,
        "Cu57_to_p_Ni56_reaclib",  // 325,
        "Cu57_to_He4_Co53_reaclib",  // 326,
        "Cu58_to_n_Cu57_reaclib",  // 327,
        "Cu58_to_p_Ni57_reaclib",  // 328,
        "Cu58_to_He4_Co54_reaclib",  // 329,
        "Cu59_to_n_Cu58_reaclib",  // 330,
        "Cu59_to_p_Ni58_reaclib",  // 331,
        "Cu59_to_He4_Co55_reaclib",  // 332,
        "Cu60_to_n_Cu59_reaclib",  // 333,
        "Cu60_to_p_Ni59_reaclib",  // 334,
        "Cu60_to_He4_Co56_reaclib",  // 335,
        "Cu61_to_n_Cu60_reaclib",  // 336,
        "Cu61_to_p_Ni60_reaclib",  // 337,
        "Cu61_to_He4_Co57_reaclib",  // 338,
        "Cu62_to_n_Cu61_reaclib",  // 339,
        "Cu62_to_p_Ni61_reaclib",  // 340,
        "Cu62_to_He4_Co58_reaclib",  // 341,
        "Cu63_to_n_Cu62_reaclib",  // 342,
        "Cu63_to_p_Ni62_reaclib",  // 343,
        "Cu63_to_He4_Co59_reaclib",  // 344,
        "Cu64_to_n_Cu63_reaclib",  // 345,
        "Cu64_to_p_Ni63_reaclib",  // 346,
        "Cu65_to_n_Cu64_reaclib",  // 347,
        "Cu65_to_p_Ni64_reaclib",  // 348,
        "Zn59_to_p_Cu58_reaclib",  // 349,
        "Zn59_to_p_Ni58_reaclib",  // 350,
        "Zn60_to_n_Zn59_reaclib",  // 351,
        "Zn60_to_p_Cu59_reaclib",  // 352,
        "Zn60_to_He4_Ni56_reaclib",  // 353,
        "Zn61_to_n_Zn60_reaclib",  // 354,
        "Zn61_to_p_Cu60_reaclib",  // 355,
        "Zn61_to_He4_Ni57_reaclib",  // 356,
        "Zn62_to_n_Zn61_reaclib",  // 357,
        "Zn62_to_p_Cu61_reaclib",  // 358,
        "Zn62_to_He4_Ni58_reaclib",  // 359,
        "Zn63_to_n_Zn62_reaclib",  // 360,
        "Zn63_to_p_Cu62_reaclib",  // 361,
        "Zn63_to_He4_Ni59_reaclib",  // 362,
        "Zn64_to_n_Zn63_reaclib",  // 363,
        "Zn64_to_p_Cu63_reaclib",  // 364,
        "Zn64_to_He4_Ni60_reaclib",  // 365,
        "Zn65_to_n_Zn64_reaclib",  // 366,
        "Zn65_to_p_Cu64_reaclib",  // 367,
        "Zn65_to_He4_Ni61_reaclib",  // 368,
        "Zn66_to_n_Zn65_reaclib",  // 369,
        "Zn66_to_p_Cu65_reaclib",  // 370,
        "Zn66_to_He4_Ni62_reaclib",  // 371,
        "Ga62_to_p_Zn61_reaclib",  // 372,
        "Ga62_to_He4_Cu58_reaclib",  // 373,
        "Ga63_to_n_Ga62_reaclib",  // 374,
        "Ga63_to_p_Zn62_reaclib",  // 375,
        "Ga63_to_He4_Cu59_reaclib",  // 376,
        "Ga64_to_n_Ga63_reaclib",  // 377,
        "Ga64_to_p_Zn63_reaclib",  // 378,
        "Ga64_to_He4_Cu60_reaclib",  // 379,
        "Ge63_to_p_Ga62_reaclib",  // 380,
        "Ge63_to_He4_Zn59_reaclib",  // 381,
        "Ge64_to_n_Ge63_reaclib",  // 382,
        "Ge64_to_p_Ga63_reaclib",  // 383,
        "Ge64_to_He4_Zn60_reaclib",  // 384,
        "Li6_to_n_p_He4_reaclib",  // 385,
        "Be9_to_n_He4_He4_reaclib",  // 386,
        "C12_to_He4_He4_He4_reaclib",  // 387,
        "n_p_to_d_reaclib",  // 388,
        "p_p_to_d_reaclib_bet_pos",  // 389,
        "p_p_to_d_reaclib_electron_capture",  // 390,
        "p_d_to_He3_reaclib",  // 391,
        "d_d_to_He4_reaclib",  // 392,
        "He4_d_to_Li6_reaclib",  // 393,
        "n_He3_to_He4_reaclib",  // 394,
        "p_He3_to_He4_reaclib",  // 395,
        "He4_He3_to_Be7_reaclib",  // 396,
        "n_Li6_to_Li7_reaclib",  // 397,
        "p_Li6_to_Be7_reaclib",  // 398,
        "He4_Li6_to_B10_reaclib",  // 399,
        "He4_Li7_to_B11_reaclib",  // 400,
        "p_Be7_to_B8_reaclib",  // 401,
        "p_Be9_to_B10_reaclib",  // 402,
        "n_B10_to_B11_reaclib",  // 403,
        "p_B11_to_C12_reaclib",  // 404,
        "n_C12_to_C13_reaclib",  // 405,
        "p_C12_to_N13_reaclib",  // 406,
        "He4_C12_to_O16_reaclib",  // 407,
        "n_C13_to_C14_reaclib",  // 408,
        "p_C13_to_N14_reaclib",  // 409,
        "p_C14_to_N15_reaclib",  // 410,
        "He4_C14_to_O18_reaclib",  // 411,
        "n_N13_to_N14_reaclib",  // 412,
        "p_N13_to_O14_reaclib",  // 413,
        "n_N14_to_N15_reaclib",  // 414,
        "p_N14_to_O15_reaclib",  // 415,
        "He4_N14_to_F18_reaclib",  // 416,
        "p_N15_to_O16_reaclib",  // 417,
        "He4_N15_to_F19_reaclib",  // 418,
        "n_O14_to_O15_reaclib",  // 419,
        "He4_O14_to_Ne18_reaclib",  // 420,
        "n_O15_to_O16_reaclib",  // 421,
        "He4_O15_to_Ne19_reaclib",  // 422,
        "n_O16_to_O17_reaclib",  // 423,
        "p_O16_to_F17_reaclib",  // 424,
        "He4_O16_to_Ne20_reaclib",  // 425,
        "n_O17_to_O18_reaclib",  // 426,
        "p_O17_to_F18_reaclib",  // 427,
        "He4_O17_to_Ne21_reaclib",  // 428,
        "p_O18_to_F19_reaclib",  // 429,
        "He4_O18_to_Ne22_reaclib",  // 430,
        "n_F17_to_F18_reaclib",  // 431,
        "p_F17_to_Ne18_reaclib",  // 432,
        "He4_F17_to_Na21_reaclib",  // 433,
        "n_F18_to_F19_reaclib",  // 434,
        "p_F18_to_Ne19_reaclib",  // 435,
        "He4_F18_to_Na22_reaclib",  // 436,
        "p_F19_to_Ne20_reaclib",  // 437,
        "He4_F19_to_Na23_reaclib",  // 438,
        "n_Ne18_to_Ne19_reaclib",  // 439,
        "n_Ne19_to_Ne20_reaclib",  // 440,
        "He4_Ne19_to_Mg23_reaclib",  // 441,
        "n_Ne20_to_Ne21_reaclib",  // 442,
        "p_Ne20_to_Na21_reaclib",  // 443,
        "He4_Ne20_to_Mg24_reaclib",  // 444,
        "n_Ne21_to_Ne22_reaclib",  // 445,
        "p_Ne21_to_Na22_reaclib",  // 446,
        "He4_Ne21_to_Mg25_reaclib",  // 447,
        "p_Ne22_to_Na23_reaclib",  // 448,
        "He4_Ne22_to_Mg26_reaclib",  // 449,
        "n_Na21_to_Na22_reaclib",  // 450,
        "He4_Na21_to_Al25_reaclib",  // 451,
        "n_Na22_to_Na23_reaclib",  // 452,
        "p_Na22_to_Mg23_reaclib",  // 453,
        "He4_Na22_to_Al26_reaclib",  // 454,
        "p_Na23_to_Mg24_reaclib",  // 455,
        "He4_Na23_to_Al27_reaclib",  // 456,
        "n_Mg23_to_Mg24_reaclib",  // 457,
        "n_Mg24_to_Mg25_reaclib",  // 458,
        "p_Mg24_to_Al25_reaclib",  // 459,
        "He4_Mg24_to_Si28_reaclib",  // 460,
        "n_Mg25_to_Mg26_reaclib",  // 461,
        "p_Mg25_to_Al26_reaclib",  // 462,
        "He4_Mg25_to_Si29_reaclib",  // 463,
        "p_Mg26_to_Al27_reaclib",  // 464,
        "He4_Mg26_to_Si30_reaclib",  // 465,
        "n_Al25_to_Al26_reaclib",  // 466,
        "He4_Al25_to_P29_reaclib",  // 467,
        "n_Al26_to_Al27_reaclib",  // 468,
        "He4_Al26_to_P30_reaclib",  // 469,
        "p_Al27_to_Si28_reaclib",  // 470,
        "He4_Al27_to_P31_reaclib",  // 471,
        "n_Si28_to_Si29_reaclib",  // 472,
        "p_Si28_to_P29_reaclib",  // 473,
        "He4_Si28_to_S32_reaclib",  // 474,
        "n_Si29_to_Si30_reaclib",  // 475,
        "p_Si29_to_P30_reaclib",  // 476,
        "He4_Si29_to_S33_reaclib",  // 477,
        "n_Si30_to_Si31_reaclib",  // 478,
        "p_Si30_to_P31_reaclib",  // 479,
        "He4_Si30_to_S34_reaclib",  // 480,
        "n_Si31_to_Si32_reaclib",  // 481,
        "p_Si31_to_P32_reaclib",  // 482,
        "He4_Si31_to_S35_reaclib",  // 483,
        "p_Si32_to_P33_reaclib",  // 484,
        "He4_Si32_to_S36_reaclib",  // 485,
        "n_P29_to_P30_reaclib",  // 486,
        "He4_P29_to_Cl33_reaclib",  // 487,
        "n_P30_to_P31_reaclib",  // 488,
        "He4_P30_to_Cl34_reaclib",  // 489,
        "n_P31_to_P32_reaclib",  // 490,
        "p_P31_to_S32_reaclib",  // 491,
        "He4_P31_to_Cl35_reaclib",  // 492,
        "n_P32_to_P33_reaclib",  // 493,
        "p_P32_to_S33_reaclib",  // 494,
        "He4_P32_to_Cl36_reaclib",  // 495,
        "p_P33_to_S34_reaclib",  // 496,
        "He4_P33_to_Cl37_reaclib",  // 497,
        "n_S32_to_S33_reaclib",  // 498,
        "p_S32_to_Cl33_reaclib",  // 499,
        "He4_S32_to_Ar36_reaclib",  // 500,
        "n_S33_to_S34_reaclib",  // 501,
        "p_S33_to_Cl34_reaclib",  // 502,
        "He4_S33_to_Ar37_reaclib",  // 503,
        "n_S34_to_S35_reaclib",  // 504,
        "p_S34_to_Cl35_reaclib",  // 505,
        "He4_S34_to_Ar38_reaclib",  // 506,
        "n_S35_to_S36_reaclib",  // 507,
        "p_S35_to_Cl36_reaclib",  // 508,
        "He4_S35_to_Ar39_reaclib",  // 509,
        "p_S36_to_Cl37_reaclib",  // 510,
        "He4_S36_to_Ar40_reaclib",  // 511,
        "n_Cl33_to_Cl34_reaclib",  // 512,
        "He4_Cl33_to_K37_reaclib",  // 513,
        "n_Cl34_to_Cl35_reaclib",  // 514,
        "He4_Cl34_to_K38_reaclib",  // 515,
        "n_Cl35_to_Cl36_reaclib",  // 516,
        "p_Cl35_to_Ar36_reaclib",  // 517,
        "He4_Cl35_to_K39_reaclib",  // 518,
        "n_Cl36_to_Cl37_reaclib",  // 519,
        "p_Cl36_to_Ar37_reaclib",  // 520,
        "He4_Cl36_to_K40_reaclib",  // 521,
        "p_Cl37_to_Ar38_reaclib",  // 522,
        "He4_Cl37_to_K41_reaclib",  // 523,
        "n_Ar36_to_Ar37_reaclib",  // 524,
        "p_Ar36_to_K37_reaclib",  // 525,
        "He4_Ar36_to_Ca40_reaclib",  // 526,
        "n_Ar37_to_Ar38_reaclib",  // 527,
        "p_Ar37_to_K38_reaclib",  // 528,
        "He4_Ar37_to_Ca41_reaclib",  // 529,
        "n_Ar38_to_Ar39_reaclib",  // 530,
        "p_Ar38_to_K39_reaclib",  // 531,
        "He4_Ar38_to_Ca42_reaclib",  // 532,
        "n_Ar39_to_Ar40_reaclib",  // 533,
        "p_Ar39_to_K40_reaclib",  // 534,
        "He4_Ar39_to_Ca43_reaclib",  // 535,
        "p_Ar40_to_K41_reaclib",  // 536,
        "He4_Ar40_to_Ca44_reaclib",  // 537,
        "n_K37_to_K38_reaclib",  // 538,
        "n_K38_to_K39_reaclib",  // 539,
        "n_K39_to_K40_reaclib",  // 540,
        "p_K39_to_Ca40_reaclib",  // 541,
        "He4_K39_to_Sc43_reaclib",  // 542,
        "n_K40_to_K41_reaclib",  // 543,
        "p_K40_to_Ca41_reaclib",  // 544,
        "He4_K40_to_Sc44_reaclib",  // 545,
        "p_K41_to_Ca42_reaclib",  // 546,
        "He4_K41_to_Sc45_reaclib",  // 547,
        "n_Ca40_to_Ca41_reaclib",  // 548,
        "He4_Ca40_to_Ti44_reaclib",  // 549,
        "n_Ca41_to_Ca42_reaclib",  // 550,
        "He4_Ca41_to_Ti45_reaclib",  // 551,
        "n_Ca42_to_Ca43_reaclib",  // 552,
        "p_Ca42_to_Sc43_reaclib",  // 553,
        "He4_Ca42_to_Ti46_reaclib",  // 554,
        "n_Ca43_to_Ca44_reaclib",  // 555,
        "p_Ca43_to_Sc44_reaclib",  // 556,
        "He4_Ca43_to_Ti47_reaclib",  // 557,
        "n_Ca44_to_Ca45_reaclib",  // 558,
        "p_Ca44_to_Sc45_reaclib",  // 559,
        "He4_Ca44_to_Ti48_reaclib",  // 560,
        "n_Ca45_to_Ca46_reaclib",  // 561,
        "p_Ca45_to_Sc46_reaclib",  // 562,
        "He4_Ca45_to_Ti49_reaclib",  // 563,
        "n_Ca46_to_Ca47_reaclib",  // 564,
        "p_Ca46_to_Sc47_reaclib",  // 565,
        "He4_Ca46_to_Ti50_reaclib",  // 566,
        "n_Ca47_to_Ca48_reaclib",  // 567,
        "p_Ca47_to_Sc48_reaclib",  // 568,
        "He4_Ca47_to_Ti51_reaclib",  // 569,
        "p_Ca48_to_Sc49_reaclib",  // 570,
        "n_Sc43_to_Sc44_reaclib",  // 571,
        "p_Sc43_to_Ti44_reaclib",  // 572,
        "He4_Sc43_to_V47_reaclib",  // 573,
        "n_Sc44_to_Sc45_reaclib",  // 574,
        "p_Sc44_to_Ti45_reaclib",  // 575,
        "He4_Sc44_to_V48_reaclib",  // 576,
        "n_Sc45_to_Sc46_reaclib",  // 577,
        "p_Sc45_to_Ti46_reaclib",  // 578,
        "He4_Sc45_to_V49_reaclib",  // 579,
        "n_Sc46_to_Sc47_reaclib",  // 580,
        "p_Sc46_to_Ti47_reaclib",  // 581,
        "He4_Sc46_to_V50_reaclib",  // 582,
        "n_Sc47_to_Sc48_reaclib",  // 583,
        "p_Sc47_to_Ti48_reaclib",  // 584,
        "He4_Sc47_to_V51_reaclib",  // 585,
        "n_Sc48_to_Sc49_reaclib",  // 586,
        "p_Sc48_to_Ti49_reaclib",  // 587,
        "He4_Sc48_to_V52_reaclib",  // 588,
        "p_Sc49_to_Ti50_reaclib",  // 589,
        "n_Ti44_to_Ti45_reaclib",  // 590,
        "He4_Ti44_to_Cr48_reaclib",  // 591,
        "n_Ti45_to_Ti46_reaclib",  // 592,
        "p_Ti45_to_V46_reaclib",  // 593,
        "He4_Ti45_to_Cr49_reaclib",  // 594,
        "n_Ti46_to_Ti47_reaclib",  // 595,
        "p_Ti46_to_V47_reaclib",  // 596,
        "He4_Ti46_to_Cr50_reaclib",  // 597,
        "n_Ti47_to_Ti48_reaclib",  // 598,
        "p_Ti47_to_V48_reaclib",  // 599,
        "He4_Ti47_to_Cr51_reaclib",  // 600,
        "n_Ti48_to_Ti49_reaclib",  // 601,
        "p_Ti48_to_V49_reaclib",  // 602,
        "He4_Ti48_to_Cr52_reaclib",  // 603,
        "n_Ti49_to_Ti50_reaclib",  // 604,
        "p_Ti49_to_V50_reaclib",  // 605,
        "He4_Ti49_to_Cr53_reaclib",  // 606,
        "n_Ti50_to_Ti51_reaclib",  // 607,
        "p_Ti50_to_V51_reaclib",  // 608,
        "He4_Ti50_to_Cr54_reaclib",  // 609,
        "p_Ti51_to_V52_reaclib",  // 610,
        "n_V46_to_V47_reaclib",  // 611,
        "He4_V46_to_Mn50_reaclib",  // 612,
        "n_V47_to_V48_reaclib",  // 613,
        "p_V47_to_Cr48_reaclib",  // 614,
        "He4_V47_to_Mn51_reaclib",  // 615,
        "n_V48_to_V49_reaclib",  // 616,
        "p_V48_to_Cr49_reaclib",  // 617,
        "He4_V48_to_Mn52_reaclib",  // 618,
        "n_V49_to_V50_reaclib",  // 619,
        "p_V49_to_Cr50_reaclib",  // 620,
        "He4_V49_to_Mn53_reaclib",  // 621,
        "n_V50_to_V51_reaclib",  // 622,
        "p_V50_to_Cr51_reaclib",  // 623,
        "He4_V50_to_Mn54_reaclib",  // 624,
        "n_V51_to_V52_reaclib",  // 625,
        "p_V51_to_Cr52_reaclib",  // 626,
        "He4_V51_to_Mn55_reaclib",  // 627,
        "p_V52_to_Cr53_reaclib",  // 628,
        "n_Cr48_to_Cr49_reaclib",  // 629,
        "He4_Cr48_to_Fe52_reaclib",  // 630,
        "n_Cr49_to_Cr50_reaclib",  // 631,
        "p_Cr49_to_Mn50_reaclib",  // 632,
        "He4_Cr49_to_Fe53_reaclib",  // 633,
        "n_Cr50_to_Cr51_reaclib",  // 634,
        "p_Cr50_to_Mn51_reaclib",  // 635,
        "He4_Cr50_to_Fe54_reaclib",  // 636,
        "n_Cr51_to_Cr52_reaclib",  // 637,
        "p_Cr51_to_Mn52_reaclib",  // 638,
        "He4_Cr51_to_Fe55_reaclib",  // 639,
        "n_Cr52_to_Cr53_reaclib",  // 640,
        "p_Cr52_to_Mn53_reaclib",  // 641,
        "He4_Cr52_to_Fe56_reaclib",  // 642,
        "n_Cr53_to_Cr54_reaclib",  // 643,
        "p_Cr53_to_Mn54_reaclib",  // 644,
        "He4_Cr53_to_Fe57_reaclib",  // 645,
        "p_Cr54_to_Mn55_reaclib",  // 646,
        "He4_Cr54_to_Fe58_reaclib",  // 647,
        "n_Mn50_to_Mn51_reaclib",  // 648,
        "He4_Mn50_to_Co54_reaclib",  // 649,
        "n_Mn51_to_Mn52_reaclib",  // 650,
        "p_Mn51_to_Fe52_reaclib",  // 651,
        "He4_Mn51_to_Co55_reaclib",  // 652,
        "n_Mn52_to_Mn53_reaclib",  // 653,
        "p_Mn52_to_Fe53_reaclib",  // 654,
        "He4_Mn52_to_Co56_reaclib",  // 655,
        "n_Mn53_to_Mn54_reaclib",  // 656,
        "p_Mn53_to_Fe54_reaclib",  // 657,
        "He4_Mn53_to_Co57_reaclib",  // 658,
        "n_Mn54_to_Mn55_reaclib",  // 659,
        "p_Mn54_to_Fe55_reaclib",  // 660,
        "He4_Mn54_to_Co58_reaclib",  // 661,
        "p_Mn55_to_Fe56_reaclib",  // 662,
        "He4_Mn55_to_Co59_reaclib",  // 663,
        "n_Fe52_to_Fe53_reaclib",  // 664,
        "p_Fe52_to_Co53_reaclib",  // 665,
        "He4_Fe52_to_Ni56_reaclib",  // 666,
        "n_Fe53_to_Fe54_reaclib",  // 667,
        "p_Fe53_to_Co54_reaclib",  // 668,
        "He4_Fe53_to_Ni57_reaclib",  // 669,
        "n_Fe54_to_Fe55_reaclib",  // 670,
        "p_Fe54_to_Co55_reaclib",  // 671,
        "He4_Fe54_to_Ni58_reaclib",  // 672,
        "n_Fe55_to_Fe56_reaclib",  // 673,
        "p_Fe55_to_Co56_reaclib",  // 674,
        "He4_Fe55_to_Ni59_reaclib",  // 675,
        "n_Fe56_to_Fe57_reaclib",  // 676,
        "p_Fe56_to_Co57_reaclib",  // 677,
        "He4_Fe56_to_Ni60_reaclib",  // 678,
        "n_Fe57_to_Fe58_reaclib",  // 679,
        "p_Fe57_to_Co58_reaclib",  // 680,
        "He4_Fe57_to_Ni61_reaclib",  // 681,
        "p_Fe58_to_Co59_reaclib",  // 682,
        "He4_Fe58_to_Ni62_reaclib",  // 683,
        "n_Co53_to_Co54_reaclib",  // 684,
        "He4_Co53_to_Cu57_reaclib",  // 685,
        "n_Co54_to_Co55_reaclib",  // 686,
        "He4_Co54_to_Cu58_reaclib",  // 687,
        "n_Co55_to_Co56_reaclib",  // 688,
        "p_Co55_to_Ni56_reaclib",  // 689,
        "He4_Co55_to_Cu59_reaclib",  // 690,
        "n_Co56_to_Co57_reaclib",  // 691,
        "p_Co56_to_Ni57_reaclib",  // 692,
        "He4_Co56_to_Cu60_reaclib",  // 693,
        "n_Co57_to_Co58_reaclib",  // 694,
        "p_Co57_to_Ni58_reaclib",  // 695,
        "He4_Co57_to_Cu61_reaclib",  // 696,
        "n_Co58_to_Co59_reaclib",  // 697,
        "p_Co58_to_Ni59_reaclib",  // 698,
        "He4_Co58_to_Cu62_reaclib",  // 699,
        "p_Co59_to_Ni60_reaclib",  // 700,
        "He4_Co59_to_Cu63_reaclib",  // 701,
        "n_Ni56_to_Ni57_reaclib",  // 702,
        "p_Ni56_to_Cu57_reaclib",  // 703,
        "He4_Ni56_to_Zn60_reaclib",  // 704,
        "n_Ni57_to_Ni58_reaclib",  // 705,
        "p_Ni57_to_Cu58_reaclib",  // 706,
        "He4_Ni57_to_Zn61_reaclib",  // 707,
        "n_Ni58_to_Ni59_reaclib",  // 708,
        "p_Ni58_to_Cu59_reaclib",  // 709,
        "He4_Ni58_to_Zn62_reaclib",  // 710,
        "n_Ni59_to_Ni60_reaclib",  // 711,
        "p_Ni59_to_Cu60_reaclib",  // 712,
        "He4_Ni59_to_Zn63_reaclib",  // 713,
        "n_Ni60_to_Ni61_reaclib",  // 714,
        "p_Ni60_to_Cu61_reaclib",  // 715,
        "He4_Ni60_to_Zn64_reaclib",  // 716,
        "n_Ni61_to_Ni62_reaclib",  // 717,
        "p_Ni61_to_Cu62_reaclib",  // 718,
        "He4_Ni61_to_Zn65_reaclib",  // 719,
        "n_Ni62_to_Ni63_reaclib",  // 720,
        "p_Ni62_to_Cu63_reaclib",  // 721,
        "He4_Ni62_to_Zn66_reaclib",  // 722,
        "n_Ni63_to_Ni64_reaclib",  // 723,
        "p_Ni63_to_Cu64_reaclib",  // 724,
        "p_Ni64_to_Cu65_reaclib",  // 725,
        "n_Cu57_to_Cu58_reaclib",  // 726,
        "n_Cu58_to_Cu59_reaclib",  // 727,
        "p_Cu58_to_Zn59_reaclib",  // 728,
        "He4_Cu58_to_Ga62_reaclib",  // 729,
        "n_Cu59_to_Cu60_reaclib",  // 730,
        "p_Cu59_to_Zn60_reaclib",  // 731,
        "He4_Cu59_to_Ga63_reaclib",  // 732,
        "n_Cu60_to_Cu61_reaclib",  // 733,
        "p_Cu60_to_Zn61_reaclib",  // 734,
        "He4_Cu60_to_Ga64_reaclib",  // 735,
        "n_Cu61_to_Cu62_reaclib",  // 736,
        "p_Cu61_to_Zn62_reaclib",  // 737,
        "n_Cu62_to_Cu63_reaclib",  // 738,
        "p_Cu62_to_Zn63_reaclib",  // 739,
        "n_Cu63_to_Cu64_reaclib",  // 740,
        "p_Cu63_to_Zn64_reaclib",  // 741,
        "n_Cu64_to_Cu65_reaclib",  // 742,
        "p_Cu64_to_Zn65_reaclib",  // 743,
        "p_Cu65_to_Zn66_reaclib",  // 744,
        "n_Zn59_to_Zn60_reaclib",  // 745,
        "He4_Zn59_to_Ge63_reaclib",  // 746,
        "n_Zn60_to_Zn61_reaclib",  // 747,
        "He4_Zn60_to_Ge64_reaclib",  // 748,
        "n_Zn61_to_Zn62_reaclib",  // 749,
        "p_Zn61_to_Ga62_reaclib",  // 750,
        "n_Zn62_to_Zn63_reaclib",  // 751,
        "p_Zn62_to_Ga63_reaclib",  // 752,
        "n_Zn63_to_Zn64_reaclib",  // 753,
        "p_Zn63_to_Ga64_reaclib",  // 754,
        "n_Zn64_to_Zn65_reaclib",  // 755,
        "n_Zn65_to_Zn66_reaclib",  // 756,
        "n_Ga62_to_Ga63_reaclib",  // 757,
        "p_Ga62_to_Ge63_reaclib",  // 758,
        "n_Ga63_to_Ga64_reaclib",  // 759,
        "p_Ga63_to_Ge64_reaclib",  // 760,
        "n_Ge63_to_Ge64_reaclib",  // 761,
        "d_d_to_n_He3_reaclib",  // 762,
        "n_He3_to_d_d_reaclib",  // 763,
        "d_He3_to_p_He4_reaclib",  // 764,
        "He4_He3_to_p_Li6_reaclib",  // 765,
        "p_He4_to_d_He3_reaclib",  // 766,
        "He4_He4_to_n_Be7_reaclib",  // 767,
        "He4_He4_to_p_Li7_reaclib",  // 768,
        "p_Li6_to_He4_He3_reaclib",  // 769,
        "d_Li6_to_n_Be7_reaclib",  // 770,
        "d_Li6_to_p_Li7_reaclib",  // 771,
        "He4_Li6_to_p_Be9_reaclib",  // 772,
        "p_Li7_to_n_Be7_reaclib",  // 773,
        "p_Li7_to_d_Li6_reaclib",  // 774,
        "p_Li7_to_He4_He4_reaclib",  // 775,
        "He4_Li7_to_n_B10_reaclib",  // 776,
        "n_Be7_to_p_Li7_reaclib",  // 777,
        "n_Be7_to_d_Li6_reaclib",  // 778,
        "n_Be7_to_He4_He4_reaclib",  // 779,
        "He4_Be7_to_p_B10_reaclib",  // 780,
        "p_Be9_to_He4_Li6_reaclib",  // 781,
        "He4_Be9_to_n_C12_reaclib",  // 782,
        "n_B10_to_He4_Li7_reaclib",  // 783,
        "p_B10_to_He4_Be7_reaclib",  // 784,
        "He4_B10_to_n_N13_reaclib",  // 785,
        "He4_B10_to_p_C13_reaclib",  // 786,
        "He4_B11_to_n_N14_reaclib",  // 787,
        "He4_B11_to_p_C14_reaclib",  // 788,
        "n_C12_to_He4_Be9_reaclib",  // 789,
        "He4_C12_to_n_O15_reaclib",  // 790,
        "He4_C12_to_p_N15_reaclib",  // 791,
        "C12_C12_to_n_Mg23_reaclib",  // 792,
        "C12_C12_to_p_Na23_reaclib",  // 793,
        "C12_C12_to_He4_Ne20_reaclib",  // 794,
        "p_C13_to_n_N13_reaclib",  // 795,
        "p_C13_to_He4_B10_reaclib",  // 796,
        "d_C13_to_n_N14_reaclib",  // 797,
        "He4_C13_to_n_O16_reaclib",  // 798,
        "p_C14_to_n_N14_reaclib",  // 799,
        "p_C14_to_He4_B11_reaclib",  // 800,
        "d_C14_to_n_N15_reaclib",  // 801,
        "He4_C14_to_n_O17_reaclib",  // 802,
        "n_N13_to_p_C13_reaclib",  // 803,
        "n_N13_to_He4_B10_reaclib",  // 804,
        "He4_N13_to_p_O16_reaclib",  // 805,
        "n_N14_to_p_C14_reaclib",  // 806,
        "n_N14_to_d_C13_reaclib",  // 807,
        "n_N14_to_He4_B11_reaclib",  // 808,
        "p_N14_to_n_O14_reaclib",  // 809,
        "He4_N14_to_n_F17_reaclib",  // 810,
        "He4_N14_to_p_O17_reaclib",  // 811,
        "n_N15_to_d_C14_reaclib",  // 812,
        "p_N15_to_n_O15_reaclib",  // 813,
        "p_N15_to_He4_C12_reaclib",  // 814,
        "He4_N15_to_n_F18_reaclib",  // 815,
        "He4_N15_to_p_O18_reaclib",  // 816,
        "n_O14_to_p_N14_reaclib",  // 817,
        "He4_O14_to_p_F17_reaclib",  // 818,
        "n_O15_to_p_N15_reaclib",  // 819,
        "n_O15_to_He4_C12_reaclib",  // 820,
        "He4_O15_to_n_Ne18_reaclib",  // 821,
        "He4_O15_to_p_F18_reaclib",  // 822,
        "n_O16_to_He4_C13_reaclib",  // 823,
        "p_O16_to_He4_N13_reaclib",  // 824,
        "He4_O16_to_n_Ne19_reaclib",  // 825,
        "He4_O16_to_p_F19_reaclib",  // 826,
        "C12_O16_to_p_Al27_reaclib",  // 827,
        "C12_O16_to_He4_Mg24_reaclib",  // 828,
        "O16_O16_to_p_P31_reaclib",  // 829,
        "O16_O16_to_He4_Si28_reaclib",  // 830,
        "n_O17_to_He4_C14_reaclib",  // 831,
        "p_O17_to_n_F17_reaclib",  // 832,
        "p_O17_to_He4_N14_reaclib",  // 833,
        "He4_O17_to_n_Ne20_reaclib",  // 834,
        "p_O18_to_n_F18_reaclib",  // 835,
        "p_O18_to_He4_N15_reaclib",  // 836,
        "He4_O18_to_n_Ne21_reaclib",  // 837,
        "n_F17_to_p_O17_reaclib",  // 838,
        "n_F17_to_He4_N14_reaclib",  // 839,
        "p_F17_to_He4_O14_reaclib",  // 840,
        "He4_F17_to_p_Ne20_reaclib",  // 841,
        "n_F18_to_p_O18_reaclib",  // 842,
        "n_F18_to_He4_N15_reaclib",  // 843,
        "p_F18_to_n_Ne18_reaclib",  // 844,
        "p_F18_to_He4_O15_reaclib",  // 845,
        "He4_F18_to_n_Na21_reaclib",  // 846,
        "He4_F18_to_p_Ne21_reaclib",  // 847,
        "p_F19_to_n_Ne19_reaclib",  // 848,
        "p_F19_to_He4_O16_reaclib",  // 849,
        "He4_F19_to_n_Na22_reaclib",  // 850,
        "He4_F19_to_p_Ne22_reaclib",  // 851,
        "n_Ne18_to_p_F18_reaclib",  // 852,
        "n_Ne18_to_He4_O15_reaclib",  // 853,
        "He4_Ne18_to_p_Na21_reaclib",  // 854,
        "n_Ne19_to_p_F19_reaclib",  // 855,
        "n_Ne19_to_He4_O16_reaclib",  // 856,
        "He4_Ne19_to_p_Na22_reaclib",  // 857,
        "n_Ne20_to_He4_O17_reaclib",  // 858,
        "p_Ne20_to_He4_F17_reaclib",  // 859,
        "He4_Ne20_to_n_Mg23_reaclib",  // 860,
        "He4_Ne20_to_p_Na23_reaclib",  // 861,
        "He4_Ne20_to_C12_C12_reaclib",  // 862,
        "C12_Ne20_to_p_P31_reaclib",  // 863,
        "C12_Ne20_to_He4_Si28_reaclib",  // 864,
        "n_Ne21_to_He4_O18_reaclib",  // 865,
        "p_Ne21_to_n_Na21_reaclib",  // 866,
        "p_Ne21_to_He4_F18_reaclib",  // 867,
        "He4_Ne21_to_n_Mg24_reaclib",  // 868,
        "p_Ne22_to_n_Na22_reaclib",  // 869,
        "p_Ne22_to_He4_F19_reaclib",  // 870,
        "He4_Ne22_to_n_Mg25_reaclib",  // 871,
        "n_Na21_to_p_Ne21_reaclib",  // 872,
        "n_Na21_to_He4_F18_reaclib",  // 873,
        "p_Na21_to_He4_Ne18_reaclib",  // 874,
        "He4_Na21_to_p_Mg24_reaclib",  // 875,
        "n_Na22_to_p_Ne22_reaclib",  // 876,
        "n_Na22_to_He4_F19_reaclib",  // 877,
        "p_Na22_to_He4_Ne19_reaclib",  // 878,
        "He4_Na22_to_n_Al25_reaclib",  // 879,
        "He4_Na22_to_p_Mg25_reaclib",  // 880,
        "p_Na23_to_n_Mg23_reaclib",  // 881,
        "p_Na23_to_He4_Ne20_reaclib",  // 882,
        "p_Na23_to_C12_C12_reaclib",  // 883,
        "He4_Na23_to_n_Al26_reaclib",  // 884,
        "He4_Na23_to_p_Mg26_reaclib",  // 885,
        "n_Mg23_to_p_Na23_reaclib",  // 886,
        "n_Mg23_to_He4_Ne20_reaclib",  // 887,
        "n_Mg23_to_C12_C12_reaclib",  // 888,
        "He4_Mg23_to_p_Al26_reaclib",  // 889,
        "n_Mg24_to_He4_Ne21_reaclib",  // 890,
        "p_Mg24_to_He4_Na21_reaclib",  // 891,
        "He4_Mg24_to_p_Al27_reaclib",  // 892,
        "He4_Mg24_to_C12_O16_reaclib",  // 893,
        "n_Mg25_to_He4_Ne22_reaclib",  // 894,
        "p_Mg25_to_n_Al25_reaclib",  // 895,
        "p_Mg25_to_He4_Na22_reaclib",  // 896,
        "He4_Mg25_to_n_Si28_reaclib",  // 897,
        "p_Mg26_to_n_Al26_reaclib",  // 898,
        "p_Mg26_to_He4_Na23_reaclib",  // 899,
        "He4_Mg26_to_n_Si29_reaclib",  // 900,
        "n_Al25_to_p_Mg25_reaclib",  // 901,
        "n_Al25_to_He4_Na22_reaclib",  // 902,
        "He4_Al25_to_p_Si28_reaclib",  // 903,
        "n_Al26_to_p_Mg26_reaclib",  // 904,
        "n_Al26_to_He4_Na23_reaclib",  // 905,
        "p_Al26_to_He4_Mg23_reaclib",  // 906,
        "He4_Al26_to_n_P29_reaclib",  // 907,
        "He4_Al26_to_p_Si29_reaclib",  // 908,
        "p_Al27_to_He4_Mg24_reaclib",  // 909,
        "p_Al27_to_C12_O16_reaclib",  // 910,
        "He4_Al27_to_n_P30_reaclib",  // 911,
        "He4_Al27_to_p_Si30_reaclib",  // 912,
        "n_Si28_to_He4_Mg25_reaclib",  // 913,
        "p_Si28_to_He4_Al25_reaclib",  // 914,
        "He4_Si28_to_p_P31_reaclib",  // 915,
        "He4_Si28_to_C12_Ne20_reaclib",  // 916,
        "He4_Si28_to_O16_O16_reaclib",  // 917,
        "n_Si29_to_He4_Mg26_reaclib",  // 918,
        "p_Si29_to_n_P29_reaclib",  // 919,
        "p_Si29_to_He4_Al26_reaclib",  // 920,
        "He4_Si29_to_n_S32_reaclib",  // 921,
        "He4_Si29_to_p_P32_reaclib",  // 922,
        "p_Si30_to_n_P30_reaclib",  // 923,
        "p_Si30_to_He4_Al27_reaclib",  // 924,
        "He4_Si30_to_n_S33_reaclib",  // 925,
        "He4_Si30_to_p_P33_reaclib",  // 926,
        "p_Si31_to_n_P31_reaclib",  // 927,
        "He4_Si31_to_n_S34_reaclib",  // 928,
        "p_Si32_to_n_P32_reaclib",  // 929,
        "He4_Si32_to_n_S35_reaclib",  // 930,
        "n_P29_to_p_Si29_reaclib",  // 931,
        "n_P29_to_He4_Al26_reaclib",  // 932,
        "He4_P29_to_p_S32_reaclib",  // 933,
        "n_P30_to_p_Si30_reaclib",  // 934,
        "n_P30_to_He4_Al27_reaclib",  // 935,
        "He4_P30_to_n_Cl33_reaclib",  // 936,
        "He4_P30_to_p_S33_reaclib",  // 937,
        "n_P31_to_p_Si31_reaclib",  // 938,
        "p_P31_to_He4_Si28_reaclib",  // 939,
        "p_P31_to_C12_Ne20_reaclib",  // 940,
        "p_P31_to_O16_O16_reaclib",  // 941,
        "He4_P31_to_n_Cl34_reaclib",  // 942,
        "He4_P31_to_p_S34_reaclib",  // 943,
        "n_P32_to_p_Si32_reaclib",  // 944,
        "p_P32_to_n_S32_reaclib",  // 945,
        "p_P32_to_He4_Si29_reaclib",  // 946,
        "He4_P32_to_n_Cl35_reaclib",  // 947,
        "He4_P32_to_p_S35_reaclib",  // 948,
        "p_P33_to_n_S33_reaclib",  // 949,
        "p_P33_to_He4_Si30_reaclib",  // 950,
        "He4_P33_to_n_Cl36_reaclib",  // 951,
        "He4_P33_to_p_S36_reaclib",  // 952,
        "n_S32_to_p_P32_reaclib",  // 953,
        "n_S32_to_He4_Si29_reaclib",  // 954,
        "p_S32_to_He4_P29_reaclib",  // 955,
        "He4_S32_to_p_Cl35_reaclib",  // 956,
        "n_S33_to_p_P33_reaclib",  // 957,
        "n_S33_to_He4_Si30_reaclib",  // 958,
        "p_S33_to_n_Cl33_reaclib",  // 959,
        "p_S33_to_He4_P30_reaclib",  // 960,
        "He4_S33_to_n_Ar36_reaclib",  // 961,
        "He4_S33_to_p_Cl36_reaclib",  // 962,
        "n_S34_to_He4_Si31_reaclib",  // 963,
        "p_S34_to_n_Cl34_reaclib",  // 964,
        "p_S34_to_He4_P31_reaclib",  // 965,
        "He4_S34_to_n_Ar37_reaclib",  // 966,
        "He4_S34_to_p_Cl37_reaclib",  // 967,
        "n_S35_to_He4_Si32_reaclib",  // 968,
        "p_S35_to_n_Cl35_reaclib",  // 969,
        "p_S35_to_He4_P32_reaclib",  // 970,
        "He4_S35_to_n_Ar38_reaclib",  // 971,
        "p_S36_to_n_Cl36_reaclib",  // 972,
        "p_S36_to_He4_P33_reaclib",  // 973,
        "He4_S36_to_n_Ar39_reaclib",  // 974,
        "n_Cl33_to_p_S33_reaclib",  // 975,
        "n_Cl33_to_He4_P30_reaclib",  // 976,
        "He4_Cl33_to_p_Ar36_reaclib",  // 977,
        "n_Cl34_to_p_S34_reaclib",  // 978,
        "n_Cl34_to_He4_P31_reaclib",  // 979,
        "He4_Cl34_to_n_K37_reaclib",  // 980,
        "He4_Cl34_to_p_Ar37_reaclib",  // 981,
        "n_Cl35_to_p_S35_reaclib",  // 982,
        "n_Cl35_to_He4_P32_reaclib",  // 983,
        "p_Cl35_to_He4_S32_reaclib",  // 984,
        "He4_Cl35_to_n_K38_reaclib",  // 985,
        "He4_Cl35_to_p_Ar38_reaclib",  // 986,
        "n_Cl36_to_p_S36_reaclib",  // 987,
        "n_Cl36_to_He4_P33_reaclib",  // 988,
        "p_Cl36_to_n_Ar36_reaclib",  // 989,
        "p_Cl36_to_He4_S33_reaclib",  // 990,
        "He4_Cl36_to_n_K39_reaclib",  // 991,
        "He4_Cl36_to_p_Ar39_reaclib",  // 992,
        "p_Cl37_to_n_Ar37_reaclib",  // 993,
        "p_Cl37_to_He4_S34_reaclib",  // 994,
        "He4_Cl37_to_n_K40_reaclib",  // 995,
        "He4_Cl37_to_p_Ar40_reaclib",  // 996,
        "n_Ar36_to_p_Cl36_reaclib",  // 997,
        "n_Ar36_to_He4_S33_reaclib",  // 998,
        "p_Ar36_to_He4_Cl33_reaclib",  // 999,
        "He4_Ar36_to_p_K39_reaclib",  // 1000,
        "n_Ar37_to_p_Cl37_reaclib",  // 1001,
        "n_Ar37_to_He4_S34_reaclib",  // 1002,
        "p_Ar37_to_n_K37_reaclib",  // 1003,
        "p_Ar37_to_He4_Cl34_reaclib",  // 1004,
        "He4_Ar37_to_n_Ca40_reaclib",  // 1005,
        "He4_Ar37_to_p_K40_reaclib",  // 1006,
        "n_Ar38_to_He4_S35_reaclib",  // 1007,
        "p_Ar38_to_n_K38_reaclib",  // 1008,
        "p_Ar38_to_He4_Cl35_reaclib",  // 1009,
        "He4_Ar38_to_n_Ca41_reaclib",  // 1010,
        "He4_Ar38_to_p_K41_reaclib",  // 1011,
        "n_Ar39_to_He4_S36_reaclib",  // 1012,
        "p_Ar39_to_n_K39_reaclib",  // 1013,
        "p_Ar39_to_He4_Cl36_reaclib",  // 1014,
        "He4_Ar39_to_n_Ca42_reaclib",  // 1015,
        "p_Ar40_to_n_K40_reaclib",  // 1016,
        "p_Ar40_to_He4_Cl37_reaclib",  // 1017,
        "He4_Ar40_to_n_Ca43_reaclib",  // 1018,
        "n_K37_to_p_Ar37_reaclib",  // 1019,
        "n_K37_to_He4_Cl34_reaclib",  // 1020,
        "He4_K37_to_p_Ca40_reaclib",  // 1021,
        "n_K38_to_p_Ar38_reaclib",  // 1022,
        "n_K38_to_He4_Cl35_reaclib",  // 1023,
        "He4_K38_to_p_Ca41_reaclib",  // 1024,
        "n_K39_to_p_Ar39_reaclib",  // 1025,
        "n_K39_to_He4_Cl36_reaclib",  // 1026,
        "p_K39_to_He4_Ar36_reaclib",  // 1027,
        "He4_K39_to_p_Ca42_reaclib",  // 1028,
        "n_K40_to_p_Ar40_reaclib",  // 1029,
        "n_K40_to_He4_Cl37_reaclib",  // 1030,
        "p_K40_to_n_Ca40_reaclib",  // 1031,
        "p_K40_to_He4_Ar37_reaclib",  // 1032,
        "He4_K40_to_n_Sc43_reaclib",  // 1033,
        "He4_K40_to_p_Ca43_reaclib",  // 1034,
        "p_K41_to_n_Ca41_reaclib",  // 1035,
        "p_K41_to_He4_Ar38_reaclib",  // 1036,
        "He4_K41_to_n_Sc44_reaclib",  // 1037,
        "He4_K41_to_p_Ca44_reaclib",  // 1038,
        "n_Ca40_to_p_K40_reaclib",  // 1039,
        "n_Ca40_to_He4_Ar37_reaclib",  // 1040,
        "p_Ca40_to_He4_K37_reaclib",  // 1041,
        "He4_Ca40_to_p_Sc43_reaclib",  // 1042,
        "n_Ca41_to_p_K41_reaclib",  // 1043,
        "n_Ca41_to_He4_Ar38_reaclib",  // 1044,
        "p_Ca41_to_He4_K38_reaclib",  // 1045,
        "He4_Ca41_to_n_Ti44_reaclib",  // 1046,
        "He4_Ca41_to_p_Sc44_reaclib",  // 1047,
        "n_Ca42_to_He4_Ar39_reaclib",  // 1048,
        "p_Ca42_to_He4_K39_reaclib",  // 1049,
        "He4_Ca42_to_n_Ti45_reaclib",  // 1050,
        "He4_Ca42_to_p_Sc45_reaclib",  // 1051,
        "n_Ca43_to_He4_Ar40_reaclib",  // 1052,
        "p_Ca43_to_n_Sc43_reaclib",  // 1053,
        "p_Ca43_to_He4_K40_reaclib",  // 1054,
        "He4_Ca43_to_n_Ti46_reaclib",  // 1055,
        "He4_Ca43_to_p_Sc46_reaclib",  // 1056,
        "p_Ca44_to_n_Sc44_reaclib",  // 1057,
        "p_Ca44_to_He4_K41_reaclib",  // 1058,
        "He4_Ca44_to_n_Ti47_reaclib",  // 1059,
        "He4_Ca44_to_p_Sc47_reaclib",  // 1060,
        "p_Ca45_to_n_Sc45_reaclib",  // 1061,
        "He4_Ca45_to_n_Ti48_reaclib",  // 1062,
        "He4_Ca45_to_p_Sc48_reaclib",  // 1063,
        "p_Ca46_to_n_Sc46_reaclib",  // 1064,
        "He4_Ca46_to_n_Ti49_reaclib",  // 1065,
        "He4_Ca46_to_p_Sc49_reaclib",  // 1066,
        "p_Ca47_to_n_Sc47_reaclib",  // 1067,
        "He4_Ca47_to_n_Ti50_reaclib",  // 1068,
        "p_Ca48_to_n_Sc48_reaclib",  // 1069,
        "He4_Ca48_to_n_Ti51_reaclib",  // 1070,
        "n_Sc43_to_p_Ca43_reaclib",  // 1071,
        "n_Sc43_to_He4_K40_reaclib",  // 1072,
        "p_Sc43_to_He4_Ca40_reaclib",  // 1073,
        "He4_Sc43_to_n_V46_reaclib",  // 1074,
        "He4_Sc43_to_p_Ti46_reaclib",  // 1075,
        "n_Sc44_to_p_Ca44_reaclib",  // 1076,
        "n_Sc44_to_He4_K41_reaclib",  // 1077,
        "p_Sc44_to_n_Ti44_reaclib",  // 1078,
        "p_Sc44_to_He4_Ca41_reaclib",  // 1079,
        "He4_Sc44_to_n_V47_reaclib",  // 1080,
        "He4_Sc44_to_p_Ti47_reaclib",  // 1081,
        "n_Sc45_to_p_Ca45_reaclib",  // 1082,
        "p_Sc45_to_n_Ti45_reaclib",  // 1083,
        "p_Sc45_to_He4_Ca42_reaclib",  // 1084,
        "He4_Sc45_to_n_V48_reaclib",  // 1085,
        "He4_Sc45_to_p_Ti48_reaclib",  // 1086,
        "n_Sc46_to_p_Ca46_reaclib",  // 1087,
        "p_Sc46_to_n_Ti46_reaclib",  // 1088,
        "p_Sc46_to_He4_Ca43_reaclib",  // 1089,
        "He4_Sc46_to_n_V49_reaclib",  // 1090,
        "He4_Sc46_to_p_Ti49_reaclib",  // 1091,
        "n_Sc47_to_p_Ca47_reaclib",  // 1092,
        "p_Sc47_to_n_Ti47_reaclib",  // 1093,
        "p_Sc47_to_He4_Ca44_reaclib",  // 1094,
        "He4_Sc47_to_n_V50_reaclib",  // 1095,
        "He4_Sc47_to_p_Ti50_reaclib",  // 1096,
        "n_Sc48_to_p_Ca48_reaclib",  // 1097,
        "p_Sc48_to_n_Ti48_reaclib",  // 1098,
        "p_Sc48_to_He4_Ca45_reaclib",  // 1099,
        "He4_Sc48_to_n_V51_reaclib",  // 1100,
        "He4_Sc48_to_p_Ti51_reaclib",  // 1101,
        "p_Sc49_to_n_Ti49_reaclib",  // 1102,
        "p_Sc49_to_He4_Ca46_reaclib",  // 1103,
        "He4_Sc49_to_n_V52_reaclib",  // 1104,
        "n_Ti44_to_p_Sc44_reaclib",  // 1105,
        "n_Ti44_to_He4_Ca41_reaclib",  // 1106,
        "He4_Ti44_to_p_V47_reaclib",  // 1107,
        "n_Ti45_to_p_Sc45_reaclib",  // 1108,
        "n_Ti45_to_He4_Ca42_reaclib",  // 1109,
        "He4_Ti45_to_n_Cr48_reaclib",  // 1110,
        "He4_Ti45_to_p_V48_reaclib",  // 1111,
        "n_Ti46_to_p_Sc46_reaclib",  // 1112,
        "n_Ti46_to_He4_Ca43_reaclib",  // 1113,
        "p_Ti46_to_n_V46_reaclib",  // 1114,
        "p_Ti46_to_He4_Sc43_reaclib",  // 1115,
        "He4_Ti46_to_n_Cr49_reaclib",  // 1116,
        "He4_Ti46_to_p_V49_reaclib",  // 1117,
        "n_Ti47_to_p_Sc47_reaclib",  // 1118,
        "n_Ti47_to_He4_Ca44_reaclib",  // 1119,
        "p_Ti47_to_n_V47_reaclib",  // 1120,
        "p_Ti47_to_He4_Sc44_reaclib",  // 1121,
        "He4_Ti47_to_n_Cr50_reaclib",  // 1122,
        "He4_Ti47_to_p_V50_reaclib",  // 1123,
        "n_Ti48_to_p_Sc48_reaclib",  // 1124,
        "n_Ti48_to_He4_Ca45_reaclib",  // 1125,
        "p_Ti48_to_n_V48_reaclib",  // 1126,
        "p_Ti48_to_He4_Sc45_reaclib",  // 1127,
        "He4_Ti48_to_n_Cr51_reaclib",  // 1128,
        "He4_Ti48_to_p_V51_reaclib",  // 1129,
        "n_Ti49_to_p_Sc49_reaclib",  // 1130,
        "n_Ti49_to_He4_Ca46_reaclib",  // 1131,
        "p_Ti49_to_n_V49_reaclib",  // 1132,
        "p_Ti49_to_He4_Sc46_reaclib",  // 1133,
        "He4_Ti49_to_n_Cr52_reaclib",  // 1134,
        "He4_Ti49_to_p_V52_reaclib",  // 1135,
        "n_Ti50_to_He4_Ca47_reaclib",  // 1136,
        "p_Ti50_to_n_V50_reaclib",  // 1137,
        "p_Ti50_to_He4_Sc47_reaclib",  // 1138,
        "He4_Ti50_to_n_Cr53_reaclib",  // 1139,
        "n_Ti51_to_He4_Ca48_reaclib",  // 1140,
        "p_Ti51_to_n_V51_reaclib",  // 1141,
        "p_Ti51_to_He4_Sc48_reaclib",  // 1142,
        "He4_Ti51_to_n_Cr54_reaclib",  // 1143,
        "n_V46_to_p_Ti46_reaclib",  // 1144,
        "n_V46_to_He4_Sc43_reaclib",  // 1145,
        "He4_V46_to_p_Cr49_reaclib",  // 1146,
        "n_V47_to_p_Ti47_reaclib",  // 1147,
        "n_V47_to_He4_Sc44_reaclib",  // 1148,
        "p_V47_to_He4_Ti44_reaclib",  // 1149,
        "He4_V47_to_n_Mn50_reaclib",  // 1150,
        "He4_V47_to_p_Cr50_reaclib",  // 1151,
        "n_V48_to_p_Ti48_reaclib",  // 1152,
        "n_V48_to_He4_Sc45_reaclib",  // 1153,
        "p_V48_to_n_Cr48_reaclib",  // 1154,
        "p_V48_to_He4_Ti45_reaclib",  // 1155,
        "He4_V48_to_n_Mn51_reaclib",  // 1156,
        "He4_V48_to_p_Cr51_reaclib",  // 1157,
        "n_V49_to_p_Ti49_reaclib",  // 1158,
        "n_V49_to_He4_Sc46_reaclib",  // 1159,
        "p_V49_to_n_Cr49_reaclib",  // 1160,
        "p_V49_to_He4_Ti46_reaclib",  // 1161,
        "He4_V49_to_n_Mn52_reaclib",  // 1162,
        "He4_V49_to_p_Cr52_reaclib",  // 1163,
        "n_V50_to_p_Ti50_reaclib",  // 1164,
        "n_V50_to_He4_Sc47_reaclib",  // 1165,
        "p_V50_to_n_Cr50_reaclib",  // 1166,
        "p_V50_to_He4_Ti47_reaclib",  // 1167,
        "He4_V50_to_n_Mn53_reaclib",  // 1168,
        "He4_V50_to_p_Cr53_reaclib",  // 1169,
        "n_V51_to_p_Ti51_reaclib",  // 1170,
        "n_V51_to_He4_Sc48_reaclib",  // 1171,
        "p_V51_to_n_Cr51_reaclib",  // 1172,
        "p_V51_to_He4_Ti48_reaclib",  // 1173,
        "He4_V51_to_n_Mn54_reaclib",  // 1174,
        "He4_V51_to_p_Cr54_reaclib",  // 1175,
        "n_V52_to_He4_Sc49_reaclib",  // 1176,
        "p_V52_to_n_Cr52_reaclib",  // 1177,
        "p_V52_to_He4_Ti49_reaclib",  // 1178,
        "He4_V52_to_n_Mn55_reaclib",  // 1179,
        "n_Cr48_to_p_V48_reaclib",  // 1180,
        "n_Cr48_to_He4_Ti45_reaclib",  // 1181,
        "He4_Cr48_to_p_Mn51_reaclib",  // 1182,
        "n_Cr49_to_p_V49_reaclib",  // 1183,
        "n_Cr49_to_He4_Ti46_reaclib",  // 1184,
        "p_Cr49_to_He4_V46_reaclib",  // 1185,
        "He4_Cr49_to_n_Fe52_reaclib",  // 1186,
        "He4_Cr49_to_p_Mn52_reaclib",  // 1187,
        "n_Cr50_to_p_V50_reaclib",  // 1188,
        "n_Cr50_to_He4_Ti47_reaclib",  // 1189,
        "p_Cr50_to_n_Mn50_reaclib",  // 1190,
        "p_Cr50_to_He4_V47_reaclib",  // 1191,
        "He4_Cr50_to_n_Fe53_reaclib",  // 1192,
        "He4_Cr50_to_p_Mn53_reaclib",  // 1193,
        "n_Cr51_to_p_V51_reaclib",  // 1194,
        "n_Cr51_to_He4_Ti48_reaclib",  // 1195,
        "p_Cr51_to_n_Mn51_reaclib",  // 1196,
        "p_Cr51_to_He4_V48_reaclib",  // 1197,
        "He4_Cr51_to_n_Fe54_reaclib",  // 1198,
        "He4_Cr51_to_p_Mn54_reaclib",  // 1199,
        "n_Cr52_to_p_V52_reaclib",  // 1200,
        "n_Cr52_to_He4_Ti49_reaclib",  // 1201,
        "p_Cr52_to_n_Mn52_reaclib",  // 1202,
        "p_Cr52_to_He4_V49_reaclib",  // 1203,
        "He4_Cr52_to_n_Fe55_reaclib",  // 1204,
        "He4_Cr52_to_p_Mn55_reaclib",  // 1205,
        "n_Cr53_to_He4_Ti50_reaclib",  // 1206,
        "p_Cr53_to_n_Mn53_reaclib",  // 1207,
        "p_Cr53_to_He4_V50_reaclib",  // 1208,
        "He4_Cr53_to_n_Fe56_reaclib",  // 1209,
        "n_Cr54_to_He4_Ti51_reaclib",  // 1210,
        "p_Cr54_to_n_Mn54_reaclib",  // 1211,
        "p_Cr54_to_He4_V51_reaclib",  // 1212,
        "He4_Cr54_to_n_Fe57_reaclib",  // 1213,
        "n_Mn50_to_p_Cr50_reaclib",  // 1214,
        "n_Mn50_to_He4_V47_reaclib",  // 1215,
        "He4_Mn50_to_n_Co53_reaclib",  // 1216,
        "He4_Mn50_to_p_Fe53_reaclib",  // 1217,
        "n_Mn51_to_p_Cr51_reaclib",  // 1218,
        "n_Mn51_to_He4_V48_reaclib",  // 1219,
        "p_Mn51_to_He4_Cr48_reaclib",  // 1220,
        "He4_Mn51_to_n_Co54_reaclib",  // 1221,
        "He4_Mn51_to_p_Fe54_reaclib",  // 1222,
        "n_Mn52_to_p_Cr52_reaclib",  // 1223,
        "n_Mn52_to_He4_V49_reaclib",  // 1224,
        "p_Mn52_to_n_Fe52_reaclib",  // 1225,
        "p_Mn52_to_He4_Cr49_reaclib",  // 1226,
        "He4_Mn52_to_n_Co55_reaclib",  // 1227,
        "He4_Mn52_to_p_Fe55_reaclib",  // 1228,
        "n_Mn53_to_p_Cr53_reaclib",  // 1229,
        "n_Mn53_to_He4_V50_reaclib",  // 1230,
        "p_Mn53_to_n_Fe53_reaclib",  // 1231,
        "p_Mn53_to_He4_Cr50_reaclib",  // 1232,
        "He4_Mn53_to_n_Co56_reaclib",  // 1233,
        "He4_Mn53_to_p_Fe56_reaclib",  // 1234,
        "n_Mn54_to_p_Cr54_reaclib",  // 1235,
        "n_Mn54_to_He4_V51_reaclib",  // 1236,
        "p_Mn54_to_n_Fe54_reaclib",  // 1237,
        "p_Mn54_to_He4_Cr51_reaclib",  // 1238,
        "He4_Mn54_to_n_Co57_reaclib",  // 1239,
        "He4_Mn54_to_p_Fe57_reaclib",  // 1240,
        "n_Mn55_to_He4_V52_reaclib",  // 1241,
        "p_Mn55_to_n_Fe55_reaclib",  // 1242,
        "p_Mn55_to_He4_Cr52_reaclib",  // 1243,
        "He4_Mn55_to_n_Co58_reaclib",  // 1244,
        "He4_Mn55_to_p_Fe58_reaclib",  // 1245,
        "n_Fe52_to_p_Mn52_reaclib",  // 1246,
        "n_Fe52_to_He4_Cr49_reaclib",  // 1247,
        "He4_Fe52_to_p_Co55_reaclib",  // 1248,
        "n_Fe53_to_p_Mn53_reaclib",  // 1249,
        "n_Fe53_to_He4_Cr50_reaclib",  // 1250,
        "p_Fe53_to_n_Co53_reaclib",  // 1251,
        "p_Fe53_to_He4_Mn50_reaclib",  // 1252,
        "He4_Fe53_to_n_Ni56_reaclib",  // 1253,
        "He4_Fe53_to_p_Co56_reaclib",  // 1254,
        "n_Fe54_to_p_Mn54_reaclib",  // 1255,
        "n_Fe54_to_He4_Cr51_reaclib",  // 1256,
        "p_Fe54_to_n_Co54_reaclib",  // 1257,
        "p_Fe54_to_He4_Mn51_reaclib",  // 1258,
        "He4_Fe54_to_n_Ni57_reaclib",  // 1259,
        "He4_Fe54_to_p_Co57_reaclib",  // 1260,
        "n_Fe55_to_p_Mn55_reaclib",  // 1261,
        "n_Fe55_to_He4_Cr52_reaclib",  // 1262,
        "p_Fe55_to_n_Co55_reaclib",  // 1263,
        "p_Fe55_to_He4_Mn52_reaclib",  // 1264,
        "He4_Fe55_to_n_Ni58_reaclib",  // 1265,
        "He4_Fe55_to_p_Co58_reaclib",  // 1266,
        "n_Fe56_to_He4_Cr53_reaclib",  // 1267,
        "p_Fe56_to_n_Co56_reaclib",  // 1268,
        "p_Fe56_to_He4_Mn53_reaclib",  // 1269,
        "He4_Fe56_to_n_Ni59_reaclib",  // 1270,
        "He4_Fe56_to_p_Co59_reaclib",  // 1271,
        "n_Fe57_to_He4_Cr54_reaclib",  // 1272,
        "p_Fe57_to_n_Co57_reaclib",  // 1273,
        "p_Fe57_to_He4_Mn54_reaclib",  // 1274,
        "He4_Fe57_to_n_Ni60_reaclib",  // 1275,
        "p_Fe58_to_n_Co58_reaclib",  // 1276,
        "p_Fe58_to_He4_Mn55_reaclib",  // 1277,
        "He4_Fe58_to_n_Ni61_reaclib",  // 1278,
        "n_Co53_to_p_Fe53_reaclib",  // 1279,
        "n_Co53_to_He4_Mn50_reaclib",  // 1280,
        "He4_Co53_to_p_Ni56_reaclib",  // 1281,
        "n_Co54_to_p_Fe54_reaclib",  // 1282,
        "n_Co54_to_He4_Mn51_reaclib",  // 1283,
        "He4_Co54_to_n_Cu57_reaclib",  // 1284,
        "He4_Co54_to_p_Ni57_reaclib",  // 1285,
        "n_Co55_to_p_Fe55_reaclib",  // 1286,
        "n_Co55_to_He4_Mn52_reaclib",  // 1287,
        "p_Co55_to_He4_Fe52_reaclib",  // 1288,
        "He4_Co55_to_n_Cu58_reaclib",  // 1289,
        "He4_Co55_to_p_Ni58_reaclib",  // 1290,
        "n_Co56_to_p_Fe56_reaclib",  // 1291,
        "n_Co56_to_He4_Mn53_reaclib",  // 1292,
        "p_Co56_to_n_Ni56_reaclib",  // 1293,
        "p_Co56_to_He4_Fe53_reaclib",  // 1294,
        "He4_Co56_to_n_Cu59_reaclib",  // 1295,
        "He4_Co56_to_p_Ni59_reaclib",  // 1296,
        "n_Co57_to_p_Fe57_reaclib",  // 1297,
        "n_Co57_to_He4_Mn54_reaclib",  // 1298,
        "p_Co57_to_n_Ni57_reaclib",  // 1299,
        "p_Co57_to_He4_Fe54_reaclib",  // 1300,
        "He4_Co57_to_n_Cu60_reaclib",  // 1301,
        "He4_Co57_to_p_Ni60_reaclib",  // 1302,
        "n_Co58_to_p_Fe58_reaclib",  // 1303,
        "n_Co58_to_He4_Mn55_reaclib",  // 1304,
        "p_Co58_to_n_Ni58_reaclib",  // 1305,
        "p_Co58_to_He4_Fe55_reaclib",  // 1306,
        "He4_Co58_to_n_Cu61_reaclib",  // 1307,
        "He4_Co58_to_p_Ni61_reaclib",  // 1308,
        "p_Co59_to_n_Ni59_reaclib",  // 1309,
        "p_Co59_to_He4_Fe56_reaclib",  // 1310,
        "He4_Co59_to_n_Cu62_reaclib",  // 1311,
        "He4_Co59_to_p_Ni62_reaclib",  // 1312,
        "n_Ni56_to_p_Co56_reaclib",  // 1313,
        "n_Ni56_to_He4_Fe53_reaclib",  // 1314,
        "p_Ni56_to_He4_Co53_reaclib",  // 1315,
        "He4_Ni56_to_n_Zn59_reaclib",  // 1316,
        "He4_Ni56_to_p_Cu59_reaclib",  // 1317,
        "n_Ni57_to_p_Co57_reaclib",  // 1318,
        "n_Ni57_to_He4_Fe54_reaclib",  // 1319,
        "p_Ni57_to_n_Cu57_reaclib",  // 1320,
        "p_Ni57_to_He4_Co54_reaclib",  // 1321,
        "He4_Ni57_to_n_Zn60_reaclib",  // 1322,
        "He4_Ni57_to_p_Cu60_reaclib",  // 1323,
        "n_Ni58_to_p_Co58_reaclib",  // 1324,
        "n_Ni58_to_He4_Fe55_reaclib",  // 1325,
        "p_Ni58_to_n_Cu58_reaclib",  // 1326,
        "p_Ni58_to_He4_Co55_reaclib",  // 1327,
        "He4_Ni58_to_n_Zn61_reaclib",  // 1328,
        "He4_Ni58_to_p_Cu61_reaclib",  // 1329,
        "n_Ni59_to_p_Co59_reaclib",  // 1330,
        "n_Ni59_to_He4_Fe56_reaclib",  // 1331,
        "p_Ni59_to_n_Cu59_reaclib",  // 1332,
        "p_Ni59_to_He4_Co56_reaclib",  // 1333,
        "He4_Ni59_to_n_Zn62_reaclib",  // 1334,
        "He4_Ni59_to_p_Cu62_reaclib",  // 1335,
        "n_Ni60_to_He4_Fe57_reaclib",  // 1336,
        "p_Ni60_to_n_Cu60_reaclib",  // 1337,
        "p_Ni60_to_He4_Co57_reaclib",  // 1338,
        "He4_Ni60_to_n_Zn63_reaclib",  // 1339,
        "He4_Ni60_to_p_Cu63_reaclib",  // 1340,
        "n_Ni61_to_He4_Fe58_reaclib",  // 1341,
        "p_Ni61_to_n_Cu61_reaclib",  // 1342,
        "p_Ni61_to_He4_Co58_reaclib",  // 1343,
        "He4_Ni61_to_n_Zn64_reaclib",  // 1344,
        "He4_Ni61_to_p_Cu64_reaclib",  // 1345,
        "p_Ni62_to_n_Cu62_reaclib",  // 1346,
        "p_Ni62_to_He4_Co59_reaclib",  // 1347,
        "He4_Ni62_to_n_Zn65_reaclib",  // 1348,
        "He4_Ni62_to_p_Cu65_reaclib",  // 1349,
        "p_Ni63_to_n_Cu63_reaclib",  // 1350,
        "He4_Ni63_to_n_Zn66_reaclib",  // 1351,
        "p_Ni64_to_n_Cu64_reaclib",  // 1352,
        "n_Cu57_to_p_Ni57_reaclib",  // 1353,
        "n_Cu57_to_He4_Co54_reaclib",  // 1354,
        "He4_Cu57_to_p_Zn60_reaclib",  // 1355,
        "n_Cu58_to_p_Ni58_reaclib",  // 1356,
        "n_Cu58_to_He4_Co55_reaclib",  // 1357,
        "He4_Cu58_to_p_Zn61_reaclib",  // 1358,
        "n_Cu59_to_p_Ni59_reaclib",  // 1359,
        "n_Cu59_to_He4_Co56_reaclib",  // 1360,
        "p_Cu59_to_n_Zn59_reaclib",  // 1361,
        "p_Cu59_to_He4_Ni56_reaclib",  // 1362,
        "He4_Cu59_to_n_Ga62_reaclib",  // 1363,
        "He4_Cu59_to_p_Zn62_reaclib",  // 1364,
        "n_Cu60_to_p_Ni60_reaclib",  // 1365,
        "n_Cu60_to_He4_Co57_reaclib",  // 1366,
        "p_Cu60_to_n_Zn60_reaclib",  // 1367,
        "p_Cu60_to_He4_Ni57_reaclib",  // 1368,
        "He4_Cu60_to_n_Ga63_reaclib",  // 1369,
        "He4_Cu60_to_p_Zn63_reaclib",  // 1370,
        "n_Cu61_to_p_Ni61_reaclib",  // 1371,
        "n_Cu61_to_He4_Co58_reaclib",  // 1372,
        "p_Cu61_to_n_Zn61_reaclib",  // 1373,
        "p_Cu61_to_He4_Ni58_reaclib",  // 1374,
        "He4_Cu61_to_n_Ga64_reaclib",  // 1375,
        "He4_Cu61_to_p_Zn64_reaclib",  // 1376,
        "n_Cu62_to_p_Ni62_reaclib",  // 1377,
        "n_Cu62_to_He4_Co59_reaclib",  // 1378,
        "p_Cu62_to_n_Zn62_reaclib",  // 1379,
        "p_Cu62_to_He4_Ni59_reaclib",  // 1380,
        "He4_Cu62_to_p_Zn65_reaclib",  // 1381,
        "n_Cu63_to_p_Ni63_reaclib",  // 1382,
        "p_Cu63_to_n_Zn63_reaclib",  // 1383,
        "p_Cu63_to_He4_Ni60_reaclib",  // 1384,
        "He4_Cu63_to_p_Zn66_reaclib",  // 1385,
        "n_Cu64_to_p_Ni64_reaclib",  // 1386,
        "p_Cu64_to_n_Zn64_reaclib",  // 1387,
        "p_Cu64_to_He4_Ni61_reaclib",  // 1388,
        "p_Cu65_to_n_Zn65_reaclib",  // 1389,
        "p_Cu65_to_He4_Ni62_reaclib",  // 1390,
        "n_Zn59_to_p_Cu59_reaclib",  // 1391,
        "n_Zn59_to_He4_Ni56_reaclib",  // 1392,
        "He4_Zn59_to_p_Ga62_reaclib",  // 1393,
        "n_Zn60_to_p_Cu60_reaclib",  // 1394,
        "n_Zn60_to_He4_Ni57_reaclib",  // 1395,
        "p_Zn60_to_He4_Cu57_reaclib",  // 1396,
        "He4_Zn60_to_n_Ge63_reaclib",  // 1397,
        "He4_Zn60_to_p_Ga63_reaclib",  // 1398,
        "n_Zn61_to_p_Cu61_reaclib",  // 1399,
        "n_Zn61_to_He4_Ni58_reaclib",  // 1400,
        "p_Zn61_to_He4_Cu58_reaclib",  // 1401,
        "He4_Zn61_to_n_Ge64_reaclib",  // 1402,
        "He4_Zn61_to_p_Ga64_reaclib",  // 1403,
        "n_Zn62_to_p_Cu62_reaclib",  // 1404,
        "n_Zn62_to_He4_Ni59_reaclib",  // 1405,
        "p_Zn62_to_n_Ga62_reaclib",  // 1406,
        "p_Zn62_to_He4_Cu59_reaclib",  // 1407,
        "n_Zn63_to_p_Cu63_reaclib",  // 1408,
        "n_Zn63_to_He4_Ni60_reaclib",  // 1409,
        "p_Zn63_to_n_Ga63_reaclib",  // 1410,
        "p_Zn63_to_He4_Cu60_reaclib",  // 1411,
        "n_Zn64_to_p_Cu64_reaclib",  // 1412,
        "n_Zn64_to_He4_Ni61_reaclib",  // 1413,
        "p_Zn64_to_n_Ga64_reaclib",  // 1414,
        "p_Zn64_to_He4_Cu61_reaclib",  // 1415,
        "n_Zn65_to_p_Cu65_reaclib",  // 1416,
        "n_Zn65_to_He4_Ni62_reaclib",  // 1417,
        "p_Zn65_to_He4_Cu62_reaclib",  // 1418,
        "n_Zn66_to_He4_Ni63_reaclib",  // 1419,
        "p_Zn66_to_He4_Cu63_reaclib",  // 1420,
        "n_Ga62_to_p_Zn62_reaclib",  // 1421,
        "n_Ga62_to_He4_Cu59_reaclib",  // 1422,
        "p_Ga62_to_He4_Zn59_reaclib",  // 1423,
        "n_Ga63_to_p_Zn63_reaclib",  // 1424,
        "n_Ga63_to_He4_Cu60_reaclib",  // 1425,
        "p_Ga63_to_n_Ge63_reaclib",  // 1426,
        "p_Ga63_to_He4_Zn60_reaclib",  // 1427,
        "n_Ga64_to_p_Zn64_reaclib",  // 1428,
        "n_Ga64_to_He4_Cu61_reaclib",  // 1429,
        "p_Ga64_to_n_Ge64_reaclib",  // 1430,
        "p_Ga64_to_He4_Zn61_reaclib",  // 1431,
        "n_Ge63_to_p_Ga63_reaclib",  // 1432,
        "n_Ge63_to_He4_Zn60_reaclib",  // 1433,
        "n_Ge64_to_p_Ga64_reaclib",  // 1434,
        "n_Ge64_to_He4_Zn61_reaclib",  // 1435,
        "p_d_to_n_p_p_reaclib",  // 1436,
        "He3_He3_to_p_p_He4_reaclib",  // 1437,
        "d_Li7_to_n_He4_He4_reaclib",  // 1438,
        "d_Be7_to_p_He4_He4_reaclib",  // 1439,
        "p_Be9_to_d_He4_He4_reaclib",  // 1440,
        "n_B8_to_p_He4_He4_reaclib",  // 1441,
        "p_B11_to_He4_He4_He4_reaclib",  // 1442,
        "He3_Li7_to_n_p_He4_He4_reaclib",  // 1443,
        "He3_Be7_to_p_p_He4_He4_reaclib",  // 1444,
        "p_Be9_to_n_p_He4_He4_reaclib",  // 1445,
        "n_p_He4_to_Li6_reaclib",  // 1446,
        "n_He4_He4_to_Be9_reaclib",  // 1447,
        "He4_He4_He4_to_C12_reaclib",  // 1448,
        "n_p_p_to_p_d_reaclib",  // 1449,
        "p_p_He4_to_He3_He3_reaclib",  // 1450,
        "n_He4_He4_to_d_Li7_reaclib",  // 1451,
        "p_He4_He4_to_n_B8_reaclib",  // 1452,
        "p_He4_He4_to_d_Be7_reaclib",  // 1453,
        "d_He4_He4_to_p_Be9_reaclib",  // 1454,
        "He4_He4_He4_to_p_B11_reaclib",  // 1455,
        "n_p_He4_He4_to_He3_Li7_reaclib",  // 1456,
        "n_p_He4_He4_to_p_Be9_reaclib",  // 1457,
        "p_p_He4_He4_to_He3_Be7_reaclib",  // 1458,
        "Na21_to_Ne21_weaktab",  // 1459,
        "Ne21_to_Na21_weaktab",  // 1460,
        "Na22_to_Ne22_weaktab",  // 1461,
        "Ne22_to_Na22_weaktab",  // 1462,
        "Mg23_to_Na23_weaktab",  // 1463,
        "Na23_to_Mg23_weaktab",  // 1464,
        "Al25_to_Mg25_weaktab",  // 1465,
        "Mg25_to_Al25_weaktab",  // 1466,
        "Al26_to_Mg26_weaktab",  // 1467,
        "Mg26_to_Al26_weaktab",  // 1468,
        "P29_to_Si29_weaktab",  // 1469,
        "Si29_to_P29_weaktab",  // 1470,
        "P30_to_Si30_weaktab",  // 1471,
        "Si30_to_P30_weaktab",  // 1472,
        "P31_to_Si31_weaktab",  // 1473,
        "Si31_to_P31_weaktab",  // 1474,
        "P32_to_S32_weaktab",  // 1475,
        "P32_to_Si32_weaktab",  // 1476,
        "S32_to_P32_weaktab",  // 1477,
        "Si32_to_P32_weaktab",  // 1478,
        "Cl33_to_S33_weaktab",  // 1479,
        "P33_to_S33_weaktab",  // 1480,
        "S33_to_Cl33_weaktab",  // 1481,
        "S33_to_P33_weaktab",  // 1482,
        "Cl34_to_S34_weaktab",  // 1483,
        "S34_to_Cl34_weaktab",  // 1484,
        "Cl35_to_S35_weaktab",  // 1485,
        "S35_to_Cl35_weaktab",  // 1486,
        "Ar36_to_Cl36_weaktab",  // 1487,
        "Cl36_to_Ar36_weaktab",  // 1488,
        "Cl36_to_S36_weaktab",  // 1489,
        "S36_to_Cl36_weaktab",  // 1490,
        "Ar37_to_Cl37_weaktab",  // 1491,
        "Ar37_to_K37_weaktab",  // 1492,
        "Cl37_to_Ar37_weaktab",  // 1493,
        "K37_to_Ar37_weaktab",  // 1494,
        "Ar38_to_K38_weaktab",  // 1495,
        "K38_to_Ar38_weaktab",  // 1496,
        "Ar39_to_K39_weaktab",  // 1497,
        "K39_to_Ar39_weaktab",  // 1498,
        "Ar40_to_K40_weaktab",  // 1499,
        "Ca40_to_K40_weaktab",  // 1500,
        "K40_to_Ar40_weaktab",  // 1501,
        "K40_to_Ca40_weaktab",  // 1502,
        "Ca41_to_K41_weaktab",  // 1503,
        "K41_to_Ca41_weaktab",  // 1504,
        "Ca43_to_Sc43_weaktab",  // 1505,
        "Sc43_to_Ca43_weaktab",  // 1506,
        "Ca44_to_Sc44_weaktab",  // 1507,
        "Sc44_to_Ca44_weaktab",  // 1508,
        "Sc44_to_Ti44_weaktab",  // 1509,
        "Ti44_to_Sc44_weaktab",  // 1510,
        "Co53_to_Fe53_weaktab",  // 1511,
        "Fe53_to_Co53_weaktab",  // 1512,
        "Cu57_to_Ni57_weaktab",  // 1513,
        "Ni57_to_Cu57_weaktab",  // 1514,
        "Ca45_to_Sc45_weaktab",  // 1515,
        "Sc45_to_Ca45_weaktab",  // 1516,
        "Sc45_to_Ti45_weaktab",  // 1517,
        "Ti45_to_Sc45_weaktab",  // 1518,
        "Ca46_to_Sc46_weaktab",  // 1519,
        "Sc46_to_Ca46_weaktab",  // 1520,
        "Sc46_to_Ti46_weaktab",  // 1521,
        "Ti46_to_Sc46_weaktab",  // 1522,
        "Ti46_to_V46_weaktab",  // 1523,
        "V46_to_Ti46_weaktab",  // 1524,
        "Ca47_to_Sc47_weaktab",  // 1525,
        "Sc47_to_Ca47_weaktab",  // 1526,
        "Sc47_to_Ti47_weaktab",  // 1527,
        "Ti47_to_Sc47_weaktab",  // 1528,
        "Ti47_to_V47_weaktab",  // 1529,
        "V47_to_Ti47_weaktab",  // 1530,
        "Ca48_to_Sc48_weaktab",  // 1531,
        "Cr48_to_V48_weaktab",  // 1532,
        "Sc48_to_Ca48_weaktab",  // 1533,
        "Sc48_to_Ti48_weaktab",  // 1534,
        "Ti48_to_Sc48_weaktab",  // 1535,
        "Ti48_to_V48_weaktab",  // 1536,
        "V48_to_Cr48_weaktab",  // 1537,
        "V48_to_Ti48_weaktab",  // 1538,
        "Cr49_to_V49_weaktab",  // 1539,
        "Sc49_to_Ti49_weaktab",  // 1540,
        "Ti49_to_Sc49_weaktab",  // 1541,
        "Ti49_to_V49_weaktab",  // 1542,
        "V49_to_Cr49_weaktab",  // 1543,
        "V49_to_Ti49_weaktab",  // 1544,
        "Cr50_to_Mn50_weaktab",  // 1545,
        "Cr50_to_V50_weaktab",  // 1546,
        "Mn50_to_Cr50_weaktab",  // 1547,
        "Ti50_to_V50_weaktab",  // 1548,
        "V50_to_Cr50_weaktab",  // 1549,
        "V50_to_Ti50_weaktab",  // 1550,
        "Cr51_to_Mn51_weaktab",  // 1551,
        "Cr51_to_V51_weaktab",  // 1552,
        "Mn51_to_Cr51_weaktab",  // 1553,
        "Ti51_to_V51_weaktab",  // 1554,
        "V51_to_Cr51_weaktab",  // 1555,
        "V51_to_Ti51_weaktab",  // 1556,
        "Cr52_to_Mn52_weaktab",  // 1557,
        "Cr52_to_V52_weaktab",  // 1558,
        "Fe52_to_Mn52_weaktab",  // 1559,
        "Mn52_to_Cr52_weaktab",  // 1560,
        "Mn52_to_Fe52_weaktab",  // 1561,
        "V52_to_Cr52_weaktab",  // 1562,
        "Cr53_to_Mn53_weaktab",  // 1563,
        "Fe53_to_Mn53_weaktab",  // 1564,
        "Mn53_to_Cr53_weaktab",  // 1565,
        "Mn53_to_Fe53_weaktab",  // 1566,
        "Co54_to_Fe54_weaktab",  // 1567,
        "Cr54_to_Mn54_weaktab",  // 1568,
        "Fe54_to_Co54_weaktab",  // 1569,
        "Fe54_to_Mn54_weaktab",  // 1570,
        "Mn54_to_Cr54_weaktab",  // 1571,
        "Mn54_to_Fe54_weaktab",  // 1572,
        "Co55_to_Fe55_weaktab",  // 1573,
        "Fe55_to_Co55_weaktab",  // 1574,
        "Fe55_to_Mn55_weaktab",  // 1575,
        "Mn55_to_Fe55_weaktab",  // 1576,
        "Co56_to_Fe56_weaktab",  // 1577,
        "Co56_to_Ni56_weaktab",  // 1578,
        "Fe56_to_Co56_weaktab",  // 1579,
        "Ni56_to_Co56_weaktab",  // 1580,
        "Co57_to_Fe57_weaktab",  // 1581,
        "Co57_to_Ni57_weaktab",  // 1582,
        "Fe57_to_Co57_weaktab",  // 1583,
        "Ni57_to_Co57_weaktab",  // 1584,
        "Co58_to_Fe58_weaktab",  // 1585,
        "Co58_to_Ni58_weaktab",  // 1586,
        "Cu58_to_Ni58_weaktab",  // 1587,
        "Fe58_to_Co58_weaktab",  // 1588,
        "Ni58_to_Co58_weaktab",  // 1589,
        "Ni58_to_Cu58_weaktab",  // 1590,
        "Co59_to_Ni59_weaktab",  // 1591,
        "Cu59_to_Ni59_weaktab",  // 1592,
        "Ni59_to_Co59_weaktab",  // 1593,
        "Ni59_to_Cu59_weaktab",  // 1594,
        "Cu60_to_Ni60_weaktab",  // 1595,
        "Cu60_to_Zn60_weaktab",  // 1596,
        "Ni60_to_Cu60_weaktab",  // 1597,
        "Zn60_to_Cu60_weaktab",  // 1598,
        "Cu61_to_Ni61_weaktab",  // 1599,
        "Cu61_to_Zn61_weaktab",  // 1600,
        "Ni61_to_Cu61_weaktab",  // 1601,
        "Zn61_to_Cu61_weaktab",  // 1602,
        "Cu62_to_Ni62_weaktab",  // 1603,
        "Cu62_to_Zn62_weaktab",  // 1604,
        "Ga62_to_Zn62_weaktab",  // 1605,
        "Ni62_to_Cu62_weaktab",  // 1606,
        "Zn62_to_Cu62_weaktab",  // 1607,
        "Zn62_to_Ga62_weaktab",  // 1608,
        "Cu63_to_Ni63_weaktab",  // 1609,
        "Cu63_to_Zn63_weaktab",  // 1610,
        "Ga63_to_Zn63_weaktab",  // 1611,
        "Ni63_to_Cu63_weaktab",  // 1612,
        "Zn63_to_Cu63_weaktab",  // 1613,
        "Zn63_to_Ga63_weaktab",  // 1614,
        "Cu64_to_Ni64_weaktab",  // 1615,
        "Cu64_to_Zn64_weaktab",  // 1616,
        "Ga64_to_Ge64_weaktab",  // 1617,
        "Ga64_to_Zn64_weaktab",  // 1618,
        "Ge64_to_Ga64_weaktab",  // 1619,
        "Ni64_to_Cu64_weaktab",  // 1620,
        "Zn64_to_Cu64_weaktab",  // 1621,
        "Zn64_to_Ga64_weaktab",  // 1622,
        "Cu65_to_Zn65_weaktab",  // 1623,
        "Zn65_to_Cu65_weaktab",  // 1624,
        "n_to_p_weaktab",  // 1625,
        "p_to_n_weaktab"  // 1626,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    // A RatePair data table, 2D array of shape (NumNSERatePairs, 8).
    // Each row corresponds to a RatePair returned by _get_nse_rate_pairs().
    // Nuclei indices follow NetworkSpecies and rate indices follow NetworkRates.

    //  - Columns 1–3: 1-based index of reactant nuclei (forward rate)
    //  - Columns 4–6: 1-based index of product nuclei (forward rate)
    //  - Columns 7-8: 1-based index of non (n, H1, He4) nuclei.
    //  - Column  9:   1-based index of the forward rate
    //  - Column  10:  1-based index of the reverse rate

    constexpr int H1_index = 2;
    constexpr int N_index = 1;
    constexpr int He4_index = 5;
    constexpr int NumNSERatePairs = 695;

    inline AMREX_GPU_MANAGED amrex::Array2D<short, 1, NumNSERatePairs, 1, 10, amrex::Order::C> rate_pair_data {
        1, 2, -1, 3, -1, -1, 3, -1, 388, 12,  // fr: n_p_to_d_reaclib, rr: d_to_n_p_reaclib
        2, 3, -1, 4, -1, -1, 3, 4, 391, 13,  // fr: p_d_to_He3_reaclib, rr: He3_to_p_d_reaclib
        3, 3, -1, 5, -1, -1, 3, 3, 392, 15,  // fr: d_d_to_He4_reaclib, rr: He4_to_d_d_reaclib
        5, 3, -1, 6, -1, -1, 3, 6, 393, 16,  // fr: He4_d_to_Li6_reaclib, rr: Li6_to_He4_d_reaclib
        1, 4, -1, 5, -1, -1, 4, -1, 394, 14,  // fr: n_He3_to_He4_reaclib, rr: He4_to_n_He3_reaclib
        5, 4, -1, 8, -1, -1, 4, 8, 396, 19,  // fr: He4_He3_to_Be7_reaclib, rr: Be7_to_He4_He3_reaclib
        1, 6, -1, 7, -1, -1, 6, 7, 397, 17,  // fr: n_Li6_to_Li7_reaclib, rr: Li7_to_n_Li6_reaclib
        2, 6, -1, 8, -1, -1, 6, 8, 398, 18,  // fr: p_Li6_to_Be7_reaclib, rr: Be7_to_p_Li6_reaclib
        5, 6, -1, 11, -1, -1, 6, 11, 399, 23,  // fr: He4_Li6_to_B10_reaclib, rr: B10_to_He4_Li6_reaclib
        5, 7, -1, 12, -1, -1, 7, 12, 400, 25,  // fr: He4_Li7_to_B11_reaclib, rr: B11_to_He4_Li7_reaclib
        2, 8, -1, 10, -1, -1, 8, 10, 401, 20,  // fr: p_Be7_to_B8_reaclib, rr: B8_to_p_Be7_reaclib
        2, 9, -1, 11, -1, -1, 9, 11, 402, 22,  // fr: p_Be9_to_B10_reaclib, rr: B10_to_p_Be9_reaclib
        1, 11, -1, 12, -1, -1, 11, 12, 403, 24,  // fr: n_B10_to_B11_reaclib, rr: B11_to_n_B10_reaclib
        2, 12, -1, 13, -1, -1, 12, 13, 404, 26,  // fr: p_B11_to_C12_reaclib, rr: C12_to_p_B11_reaclib
        1, 13, -1, 14, -1, -1, 13, 14, 405, 27,  // fr: n_C12_to_C13_reaclib, rr: C13_to_n_C12_reaclib
        2, 13, -1, 16, -1, -1, 13, 16, 406, 29,  // fr: p_C12_to_N13_reaclib, rr: N13_to_p_C12_reaclib
        5, 13, -1, 21, -1, -1, 13, 21, 407, 39,  // fr: He4_C12_to_O16_reaclib, rr: O16_to_He4_C12_reaclib
        1, 14, -1, 15, -1, -1, 14, 15, 408, 28,  // fr: n_C13_to_C14_reaclib, rr: C14_to_n_C13_reaclib
        2, 14, -1, 17, -1, -1, 14, 17, 409, 31,  // fr: p_C13_to_N14_reaclib, rr: N14_to_p_C13_reaclib
        2, 15, -1, 18, -1, -1, 15, 18, 410, 33,  // fr: p_C14_to_N15_reaclib, rr: N15_to_p_C14_reaclib
        5, 15, -1, 23, -1, -1, 15, 23, 411, 42,  // fr: He4_C14_to_O18_reaclib, rr: O18_to_He4_C14_reaclib
        1, 16, -1, 17, -1, -1, 16, 17, 412, 30,  // fr: n_N13_to_N14_reaclib, rr: N14_to_n_N13_reaclib
        2, 16, -1, 19, -1, -1, 16, 19, 413, 34,  // fr: p_N13_to_O14_reaclib, rr: O14_to_p_N13_reaclib
        1, 17, -1, 18, -1, -1, 17, 18, 414, 32,  // fr: n_N14_to_N15_reaclib, rr: N15_to_n_N14_reaclib
        2, 17, -1, 20, -1, -1, 17, 20, 415, 36,  // fr: p_N14_to_O15_reaclib, rr: O15_to_p_N14_reaclib
        5, 17, -1, 25, -1, -1, 17, 25, 416, 46,  // fr: He4_N14_to_F18_reaclib, rr: F18_to_He4_N14_reaclib
        2, 18, -1, 21, -1, -1, 18, 21, 417, 38,  // fr: p_N15_to_O16_reaclib, rr: O16_to_p_N15_reaclib
        5, 18, -1, 26, -1, -1, 18, 26, 418, 49,  // fr: He4_N15_to_F19_reaclib, rr: F19_to_He4_N15_reaclib
        1, 19, -1, 20, -1, -1, 19, 20, 419, 35,  // fr: n_O14_to_O15_reaclib, rr: O15_to_n_O14_reaclib
        5, 19, -1, 27, -1, -1, 19, 27, 420, 51,  // fr: He4_O14_to_Ne18_reaclib, rr: Ne18_to_He4_O14_reaclib
        1, 20, -1, 21, -1, -1, 20, 21, 421, 37,  // fr: n_O15_to_O16_reaclib, rr: O16_to_n_O15_reaclib
        5, 20, -1, 28, -1, -1, 20, 28, 422, 54,  // fr: He4_O15_to_Ne19_reaclib, rr: Ne19_to_He4_O15_reaclib
        1, 21, -1, 22, -1, -1, 21, 22, 423, 40,  // fr: n_O16_to_O17_reaclib, rr: O17_to_n_O16_reaclib
        2, 21, -1, 24, -1, -1, 21, 24, 424, 43,  // fr: p_O16_to_F17_reaclib, rr: F17_to_p_O16_reaclib
        5, 21, -1, 29, -1, -1, 21, 29, 425, 57,  // fr: He4_O16_to_Ne20_reaclib, rr: Ne20_to_He4_O16_reaclib
        1, 22, -1, 23, -1, -1, 22, 23, 426, 41,  // fr: n_O17_to_O18_reaclib, rr: O18_to_n_O17_reaclib
        2, 22, -1, 25, -1, -1, 22, 25, 427, 45,  // fr: p_O17_to_F18_reaclib, rr: F18_to_p_O17_reaclib
        5, 22, -1, 30, -1, -1, 22, 30, 428, 59,  // fr: He4_O17_to_Ne21_reaclib, rr: Ne21_to_He4_O17_reaclib
        2, 23, -1, 26, -1, -1, 23, 26, 429, 48,  // fr: p_O18_to_F19_reaclib, rr: F19_to_p_O18_reaclib
        5, 23, -1, 31, -1, -1, 23, 31, 430, 61,  // fr: He4_O18_to_Ne22_reaclib, rr: Ne22_to_He4_O18_reaclib
        1, 24, -1, 25, -1, -1, 24, 25, 431, 44,  // fr: n_F17_to_F18_reaclib, rr: F18_to_n_F17_reaclib
        2, 24, -1, 27, -1, -1, 24, 27, 432, 50,  // fr: p_F17_to_Ne18_reaclib, rr: Ne18_to_p_F17_reaclib
        5, 24, -1, 32, -1, -1, 24, 32, 433, 63,  // fr: He4_F17_to_Na21_reaclib, rr: Na21_to_He4_F17_reaclib
        1, 25, -1, 26, -1, -1, 25, 26, 434, 47,  // fr: n_F18_to_F19_reaclib, rr: F19_to_n_F18_reaclib
        2, 25, -1, 28, -1, -1, 25, 28, 435, 53,  // fr: p_F18_to_Ne19_reaclib, rr: Ne19_to_p_F18_reaclib
        5, 25, -1, 33, -1, -1, 25, 33, 436, 66,  // fr: He4_F18_to_Na22_reaclib, rr: Na22_to_He4_F18_reaclib
        2, 26, -1, 29, -1, -1, 26, 29, 437, 56,  // fr: p_F19_to_Ne20_reaclib, rr: Ne20_to_p_F19_reaclib
        5, 26, -1, 34, -1, -1, 26, 34, 438, 69,  // fr: He4_F19_to_Na23_reaclib, rr: Na23_to_He4_F19_reaclib
        1, 27, -1, 28, -1, -1, 27, 28, 439, 52,  // fr: n_Ne18_to_Ne19_reaclib, rr: Ne19_to_n_Ne18_reaclib
        1, 28, -1, 29, -1, -1, 28, 29, 440, 55,  // fr: n_Ne19_to_Ne20_reaclib, rr: Ne20_to_n_Ne19_reaclib
        5, 28, -1, 35, -1, -1, 28, 35, 441, 71,  // fr: He4_Ne19_to_Mg23_reaclib, rr: Mg23_to_He4_Ne19_reaclib
        1, 29, -1, 30, -1, -1, 29, 30, 442, 58,  // fr: n_Ne20_to_Ne21_reaclib, rr: Ne21_to_n_Ne20_reaclib
        2, 29, -1, 32, -1, -1, 29, 32, 443, 62,  // fr: p_Ne20_to_Na21_reaclib, rr: Na21_to_p_Ne20_reaclib
        5, 29, -1, 36, -1, -1, 29, 36, 444, 74,  // fr: He4_Ne20_to_Mg24_reaclib, rr: Mg24_to_He4_Ne20_reaclib
        1, 30, -1, 31, -1, -1, 30, 31, 445, 60,  // fr: n_Ne21_to_Ne22_reaclib, rr: Ne22_to_n_Ne21_reaclib
        2, 30, -1, 33, -1, -1, 30, 33, 446, 65,  // fr: p_Ne21_to_Na22_reaclib, rr: Na22_to_p_Ne21_reaclib
        5, 30, -1, 37, -1, -1, 30, 37, 447, 76,  // fr: He4_Ne21_to_Mg25_reaclib, rr: Mg25_to_He4_Ne21_reaclib
        2, 31, -1, 34, -1, -1, 31, 34, 448, 68,  // fr: p_Ne22_to_Na23_reaclib, rr: Na23_to_p_Ne22_reaclib
        5, 31, -1, 38, -1, -1, 31, 38, 449, 78,  // fr: He4_Ne22_to_Mg26_reaclib, rr: Mg26_to_He4_Ne22_reaclib
        1, 32, -1, 33, -1, -1, 32, 33, 450, 64,  // fr: n_Na21_to_Na22_reaclib, rr: Na22_to_n_Na21_reaclib
        5, 32, -1, 39, -1, -1, 32, 39, 451, 80,  // fr: He4_Na21_to_Al25_reaclib, rr: Al25_to_He4_Na21_reaclib
        1, 33, -1, 34, -1, -1, 33, 34, 452, 67,  // fr: n_Na22_to_Na23_reaclib, rr: Na23_to_n_Na22_reaclib
        2, 33, -1, 35, -1, -1, 33, 35, 453, 70,  // fr: p_Na22_to_Mg23_reaclib, rr: Mg23_to_p_Na22_reaclib
        5, 33, -1, 40, -1, -1, 33, 40, 454, 83,  // fr: He4_Na22_to_Al26_reaclib, rr: Al26_to_He4_Na22_reaclib
        2, 34, -1, 36, -1, -1, 34, 36, 455, 73,  // fr: p_Na23_to_Mg24_reaclib, rr: Mg24_to_p_Na23_reaclib
        5, 34, -1, 41, -1, -1, 34, 41, 456, 86,  // fr: He4_Na23_to_Al27_reaclib, rr: Al27_to_He4_Na23_reaclib
        1, 35, -1, 36, -1, -1, 35, 36, 457, 72,  // fr: n_Mg23_to_Mg24_reaclib, rr: Mg24_to_n_Mg23_reaclib
        1, 36, -1, 37, -1, -1, 36, 37, 458, 75,  // fr: n_Mg24_to_Mg25_reaclib, rr: Mg25_to_n_Mg24_reaclib
        2, 36, -1, 39, -1, -1, 36, 39, 459, 79,  // fr: p_Mg24_to_Al25_reaclib, rr: Al25_to_p_Mg24_reaclib
        5, 36, -1, 42, -1, -1, 36, 42, 460, 88,  // fr: He4_Mg24_to_Si28_reaclib, rr: Si28_to_He4_Mg24_reaclib
        1, 37, -1, 38, -1, -1, 37, 38, 461, 77,  // fr: n_Mg25_to_Mg26_reaclib, rr: Mg26_to_n_Mg25_reaclib
        2, 37, -1, 40, -1, -1, 37, 40, 462, 82,  // fr: p_Mg25_to_Al26_reaclib, rr: Al26_to_p_Mg25_reaclib
        5, 37, -1, 43, -1, -1, 37, 43, 463, 90,  // fr: He4_Mg25_to_Si29_reaclib, rr: Si29_to_He4_Mg25_reaclib
        2, 38, -1, 41, -1, -1, 38, 41, 464, 85,  // fr: p_Mg26_to_Al27_reaclib, rr: Al27_to_p_Mg26_reaclib
        5, 38, -1, 44, -1, -1, 38, 44, 465, 92,  // fr: He4_Mg26_to_Si30_reaclib, rr: Si30_to_He4_Mg26_reaclib
        1, 39, -1, 40, -1, -1, 39, 40, 466, 81,  // fr: n_Al25_to_Al26_reaclib, rr: Al26_to_n_Al25_reaclib
        5, 39, -1, 47, -1, -1, 39, 47, 467, 96,  // fr: He4_Al25_to_P29_reaclib, rr: P29_to_He4_Al25_reaclib
        1, 40, -1, 41, -1, -1, 40, 41, 468, 84,  // fr: n_Al26_to_Al27_reaclib, rr: Al27_to_n_Al26_reaclib
        5, 40, -1, 48, -1, -1, 40, 48, 469, 99,  // fr: He4_Al26_to_P30_reaclib, rr: P30_to_He4_Al26_reaclib
        2, 41, -1, 42, -1, -1, 41, 42, 470, 87,  // fr: p_Al27_to_Si28_reaclib, rr: Si28_to_p_Al27_reaclib
        5, 41, -1, 49, -1, -1, 41, 49, 471, 102,  // fr: He4_Al27_to_P31_reaclib, rr: P31_to_He4_Al27_reaclib
        1, 42, -1, 43, -1, -1, 42, 43, 472, 89,  // fr: n_Si28_to_Si29_reaclib, rr: Si29_to_n_Si28_reaclib
        2, 42, -1, 47, -1, -1, 42, 47, 473, 95,  // fr: p_Si28_to_P29_reaclib, rr: P29_to_p_Si28_reaclib
        5, 42, -1, 52, -1, -1, 42, 52, 474, 108,  // fr: He4_Si28_to_S32_reaclib, rr: S32_to_He4_Si28_reaclib
        1, 43, -1, 44, -1, -1, 43, 44, 475, 91,  // fr: n_Si29_to_Si30_reaclib, rr: Si30_to_n_Si29_reaclib
        2, 43, -1, 48, -1, -1, 43, 48, 476, 98,  // fr: p_Si29_to_P30_reaclib, rr: P30_to_p_Si29_reaclib
        5, 43, -1, 53, -1, -1, 43, 53, 477, 111,  // fr: He4_Si29_to_S33_reaclib, rr: S33_to_He4_Si29_reaclib
        1, 44, -1, 45, -1, -1, 44, 45, 478, 93,  // fr: n_Si30_to_Si31_reaclib, rr: Si31_to_n_Si30_reaclib
        2, 44, -1, 49, -1, -1, 44, 49, 479, 101,  // fr: p_Si30_to_P31_reaclib, rr: P31_to_p_Si30_reaclib
        5, 44, -1, 54, -1, -1, 44, 54, 480, 114,  // fr: He4_Si30_to_S34_reaclib, rr: S34_to_He4_Si30_reaclib
        1, 45, -1, 46, -1, -1, 45, 46, 481, 94,  // fr: n_Si31_to_Si32_reaclib, rr: Si32_to_n_Si31_reaclib
        2, 45, -1, 50, -1, -1, 45, 50, 482, 104,  // fr: p_Si31_to_P32_reaclib, rr: P32_to_p_Si31_reaclib
        5, 45, -1, 55, -1, -1, 45, 55, 483, 116,  // fr: He4_Si31_to_S35_reaclib, rr: S35_to_He4_Si31_reaclib
        2, 46, -1, 51, -1, -1, 46, 51, 484, 106,  // fr: p_Si32_to_P33_reaclib, rr: P33_to_p_Si32_reaclib
        5, 46, -1, 56, -1, -1, 46, 56, 485, 118,  // fr: He4_Si32_to_S36_reaclib, rr: S36_to_He4_Si32_reaclib
        1, 47, -1, 48, -1, -1, 47, 48, 486, 97,  // fr: n_P29_to_P30_reaclib, rr: P30_to_n_P29_reaclib
        5, 47, -1, 57, -1, -1, 47, 57, 487, 120,  // fr: He4_P29_to_Cl33_reaclib, rr: Cl33_to_He4_P29_reaclib
        1, 48, -1, 49, -1, -1, 48, 49, 488, 100,  // fr: n_P30_to_P31_reaclib, rr: P31_to_n_P30_reaclib
        5, 48, -1, 58, -1, -1, 48, 58, 489, 123,  // fr: He4_P30_to_Cl34_reaclib, rr: Cl34_to_He4_P30_reaclib
        1, 49, -1, 50, -1, -1, 49, 50, 490, 103,  // fr: n_P31_to_P32_reaclib, rr: P32_to_n_P31_reaclib
        2, 49, -1, 52, -1, -1, 49, 52, 491, 107,  // fr: p_P31_to_S32_reaclib, rr: S32_to_p_P31_reaclib
        5, 49, -1, 59, -1, -1, 49, 59, 492, 126,  // fr: He4_P31_to_Cl35_reaclib, rr: Cl35_to_He4_P31_reaclib
        1, 50, -1, 51, -1, -1, 50, 51, 493, 105,  // fr: n_P32_to_P33_reaclib, rr: P33_to_n_P32_reaclib
        2, 50, -1, 53, -1, -1, 50, 53, 494, 110,  // fr: p_P32_to_S33_reaclib, rr: S33_to_p_P32_reaclib
        5, 50, -1, 60, -1, -1, 50, 60, 495, 129,  // fr: He4_P32_to_Cl36_reaclib, rr: Cl36_to_He4_P32_reaclib
        2, 51, -1, 54, -1, -1, 51, 54, 496, 113,  // fr: p_P33_to_S34_reaclib, rr: S34_to_p_P33_reaclib
        5, 51, -1, 61, -1, -1, 51, 61, 497, 132,  // fr: He4_P33_to_Cl37_reaclib, rr: Cl37_to_He4_P33_reaclib
        1, 52, -1, 53, -1, -1, 52, 53, 498, 109,  // fr: n_S32_to_S33_reaclib, rr: S33_to_n_S32_reaclib
        2, 52, -1, 57, -1, -1, 52, 57, 499, 119,  // fr: p_S32_to_Cl33_reaclib, rr: Cl33_to_p_S32_reaclib
        5, 52, -1, 62, -1, -1, 52, 62, 500, 134,  // fr: He4_S32_to_Ar36_reaclib, rr: Ar36_to_He4_S32_reaclib
        1, 53, -1, 54, -1, -1, 53, 54, 501, 112,  // fr: n_S33_to_S34_reaclib, rr: S34_to_n_S33_reaclib
        2, 53, -1, 58, -1, -1, 53, 58, 502, 122,  // fr: p_S33_to_Cl34_reaclib, rr: Cl34_to_p_S33_reaclib
        5, 53, -1, 63, -1, -1, 53, 63, 503, 137,  // fr: He4_S33_to_Ar37_reaclib, rr: Ar37_to_He4_S33_reaclib
        1, 54, -1, 55, -1, -1, 54, 55, 504, 115,  // fr: n_S34_to_S35_reaclib, rr: S35_to_n_S34_reaclib
        2, 54, -1, 59, -1, -1, 54, 59, 505, 125,  // fr: p_S34_to_Cl35_reaclib, rr: Cl35_to_p_S34_reaclib
        5, 54, -1, 64, -1, -1, 54, 64, 506, 140,  // fr: He4_S34_to_Ar38_reaclib, rr: Ar38_to_He4_S34_reaclib
        1, 55, -1, 56, -1, -1, 55, 56, 507, 117,  // fr: n_S35_to_S36_reaclib, rr: S36_to_n_S35_reaclib
        2, 55, -1, 60, -1, -1, 55, 60, 508, 128,  // fr: p_S35_to_Cl36_reaclib, rr: Cl36_to_p_S35_reaclib
        5, 55, -1, 65, -1, -1, 55, 65, 509, 142,  // fr: He4_S35_to_Ar39_reaclib, rr: Ar39_to_He4_S35_reaclib
        2, 56, -1, 61, -1, -1, 56, 61, 510, 131,  // fr: p_S36_to_Cl37_reaclib, rr: Cl37_to_p_S36_reaclib
        5, 56, -1, 66, -1, -1, 56, 66, 511, 144,  // fr: He4_S36_to_Ar40_reaclib, rr: Ar40_to_He4_S36_reaclib
        1, 57, -1, 58, -1, -1, 57, 58, 512, 121,  // fr: n_Cl33_to_Cl34_reaclib, rr: Cl34_to_n_Cl33_reaclib
        5, 57, -1, 67, -1, -1, 57, 67, 513, 146,  // fr: He4_Cl33_to_K37_reaclib, rr: K37_to_He4_Cl33_reaclib
        1, 58, -1, 59, -1, -1, 58, 59, 514, 124,  // fr: n_Cl34_to_Cl35_reaclib, rr: Cl35_to_n_Cl34_reaclib
        5, 58, -1, 68, -1, -1, 58, 68, 515, 149,  // fr: He4_Cl34_to_K38_reaclib, rr: K38_to_He4_Cl34_reaclib
        1, 59, -1, 60, -1, -1, 59, 60, 516, 127,  // fr: n_Cl35_to_Cl36_reaclib, rr: Cl36_to_n_Cl35_reaclib
        2, 59, -1, 62, -1, -1, 59, 62, 517, 133,  // fr: p_Cl35_to_Ar36_reaclib, rr: Ar36_to_p_Cl35_reaclib
        5, 59, -1, 69, -1, -1, 59, 69, 518, 152,  // fr: He4_Cl35_to_K39_reaclib, rr: K39_to_He4_Cl35_reaclib
        1, 60, -1, 61, -1, -1, 60, 61, 519, 130,  // fr: n_Cl36_to_Cl37_reaclib, rr: Cl37_to_n_Cl36_reaclib
        2, 60, -1, 63, -1, -1, 60, 63, 520, 136,  // fr: p_Cl36_to_Ar37_reaclib, rr: Ar37_to_p_Cl36_reaclib
        5, 60, -1, 70, -1, -1, 60, 70, 521, 155,  // fr: He4_Cl36_to_K40_reaclib, rr: K40_to_He4_Cl36_reaclib
        2, 61, -1, 64, -1, -1, 61, 64, 522, 139,  // fr: p_Cl37_to_Ar38_reaclib, rr: Ar38_to_p_Cl37_reaclib
        5, 61, -1, 71, -1, -1, 61, 71, 523, 158,  // fr: He4_Cl37_to_K41_reaclib, rr: K41_to_He4_Cl37_reaclib
        1, 62, -1, 63, -1, -1, 62, 63, 524, 135,  // fr: n_Ar36_to_Ar37_reaclib, rr: Ar37_to_n_Ar36_reaclib
        2, 62, -1, 67, -1, -1, 62, 67, 525, 145,  // fr: p_Ar36_to_K37_reaclib, rr: K37_to_p_Ar36_reaclib
        5, 62, -1, 72, -1, -1, 62, 72, 526, 160,  // fr: He4_Ar36_to_Ca40_reaclib, rr: Ca40_to_He4_Ar36_reaclib
        1, 63, -1, 64, -1, -1, 63, 64, 527, 138,  // fr: n_Ar37_to_Ar38_reaclib, rr: Ar38_to_n_Ar37_reaclib
        2, 63, -1, 68, -1, -1, 63, 68, 528, 148,  // fr: p_Ar37_to_K38_reaclib, rr: K38_to_p_Ar37_reaclib
        5, 63, -1, 73, -1, -1, 63, 73, 529, 163,  // fr: He4_Ar37_to_Ca41_reaclib, rr: Ca41_to_He4_Ar37_reaclib
        1, 64, -1, 65, -1, -1, 64, 65, 530, 141,  // fr: n_Ar38_to_Ar39_reaclib, rr: Ar39_to_n_Ar38_reaclib
        2, 64, -1, 69, -1, -1, 64, 69, 531, 151,  // fr: p_Ar38_to_K39_reaclib, rr: K39_to_p_Ar38_reaclib
        5, 64, -1, 74, -1, -1, 64, 74, 532, 166,  // fr: He4_Ar38_to_Ca42_reaclib, rr: Ca42_to_He4_Ar38_reaclib
        1, 65, -1, 66, -1, -1, 65, 66, 533, 143,  // fr: n_Ar39_to_Ar40_reaclib, rr: Ar40_to_n_Ar39_reaclib
        2, 65, -1, 70, -1, -1, 65, 70, 534, 154,  // fr: p_Ar39_to_K40_reaclib, rr: K40_to_p_Ar39_reaclib
        5, 65, -1, 75, -1, -1, 65, 75, 535, 168,  // fr: He4_Ar39_to_Ca43_reaclib, rr: Ca43_to_He4_Ar39_reaclib
        2, 66, -1, 71, -1, -1, 66, 71, 536, 157,  // fr: p_Ar40_to_K41_reaclib, rr: K41_to_p_Ar40_reaclib
        5, 66, -1, 76, -1, -1, 66, 76, 537, 170,  // fr: He4_Ar40_to_Ca44_reaclib, rr: Ca44_to_He4_Ar40_reaclib
        1, 67, -1, 68, -1, -1, 67, 68, 538, 147,  // fr: n_K37_to_K38_reaclib, rr: K38_to_n_K37_reaclib
        1, 68, -1, 69, -1, -1, 68, 69, 539, 150,  // fr: n_K38_to_K39_reaclib, rr: K39_to_n_K38_reaclib
        1, 69, -1, 70, -1, -1, 69, 70, 540, 153,  // fr: n_K39_to_K40_reaclib, rr: K40_to_n_K39_reaclib
        2, 69, -1, 72, -1, -1, 69, 72, 541, 159,  // fr: p_K39_to_Ca40_reaclib, rr: Ca40_to_p_K39_reaclib
        5, 69, -1, 81, -1, -1, 69, 81, 542, 176,  // fr: He4_K39_to_Sc43_reaclib, rr: Sc43_to_He4_K39_reaclib
        1, 70, -1, 71, -1, -1, 70, 71, 543, 156,  // fr: n_K40_to_K41_reaclib, rr: K41_to_n_K40_reaclib
        2, 70, -1, 73, -1, -1, 70, 73, 544, 162,  // fr: p_K40_to_Ca41_reaclib, rr: Ca41_to_p_K40_reaclib
        5, 70, -1, 82, -1, -1, 70, 82, 545, 179,  // fr: He4_K40_to_Sc44_reaclib, rr: Sc44_to_He4_K40_reaclib
        2, 71, -1, 74, -1, -1, 71, 74, 546, 165,  // fr: p_K41_to_Ca42_reaclib, rr: Ca42_to_p_K41_reaclib
        5, 71, -1, 83, -1, -1, 71, 83, 547, 182,  // fr: He4_K41_to_Sc45_reaclib, rr: Sc45_to_He4_K41_reaclib
        1, 72, -1, 73, -1, -1, 72, 73, 548, 161,  // fr: n_Ca40_to_Ca41_reaclib, rr: Ca41_to_n_Ca40_reaclib
        5, 72, -1, 88, -1, -1, 72, 88, 549, 192,  // fr: He4_Ca40_to_Ti44_reaclib, rr: Ti44_to_He4_Ca40_reaclib
        1, 73, -1, 74, -1, -1, 73, 74, 550, 164,  // fr: n_Ca41_to_Ca42_reaclib, rr: Ca42_to_n_Ca41_reaclib
        5, 73, -1, 89, -1, -1, 73, 89, 551, 195,  // fr: He4_Ca41_to_Ti45_reaclib, rr: Ti45_to_He4_Ca41_reaclib
        1, 74, -1, 75, -1, -1, 74, 75, 552, 167,  // fr: n_Ca42_to_Ca43_reaclib, rr: Ca43_to_n_Ca42_reaclib
        2, 74, -1, 81, -1, -1, 74, 81, 553, 175,  // fr: p_Ca42_to_Sc43_reaclib, rr: Sc43_to_p_Ca42_reaclib
        5, 74, -1, 90, -1, -1, 74, 90, 554, 198,  // fr: He4_Ca42_to_Ti46_reaclib, rr: Ti46_to_He4_Ca42_reaclib
        1, 75, -1, 76, -1, -1, 75, 76, 555, 169,  // fr: n_Ca43_to_Ca44_reaclib, rr: Ca44_to_n_Ca43_reaclib
        2, 75, -1, 82, -1, -1, 75, 82, 556, 178,  // fr: p_Ca43_to_Sc44_reaclib, rr: Sc44_to_p_Ca43_reaclib
        5, 75, -1, 91, -1, -1, 75, 91, 557, 201,  // fr: He4_Ca43_to_Ti47_reaclib, rr: Ti47_to_He4_Ca43_reaclib
        1, 76, -1, 77, -1, -1, 76, 77, 558, 171,  // fr: n_Ca44_to_Ca45_reaclib, rr: Ca45_to_n_Ca44_reaclib
        2, 76, -1, 83, -1, -1, 76, 83, 559, 181,  // fr: p_Ca44_to_Sc45_reaclib, rr: Sc45_to_p_Ca44_reaclib
        5, 76, -1, 92, -1, -1, 76, 92, 560, 204,  // fr: He4_Ca44_to_Ti48_reaclib, rr: Ti48_to_He4_Ca44_reaclib
        1, 77, -1, 78, -1, -1, 77, 78, 561, 172,  // fr: n_Ca45_to_Ca46_reaclib, rr: Ca46_to_n_Ca45_reaclib
        2, 77, -1, 84, -1, -1, 77, 84, 562, 184,  // fr: p_Ca45_to_Sc46_reaclib, rr: Sc46_to_p_Ca45_reaclib
        5, 77, -1, 93, -1, -1, 77, 93, 563, 207,  // fr: He4_Ca45_to_Ti49_reaclib, rr: Ti49_to_He4_Ca45_reaclib
        1, 78, -1, 79, -1, -1, 78, 79, 564, 173,  // fr: n_Ca46_to_Ca47_reaclib, rr: Ca47_to_n_Ca46_reaclib
        2, 78, -1, 85, -1, -1, 78, 85, 565, 186,  // fr: p_Ca46_to_Sc47_reaclib, rr: Sc47_to_p_Ca46_reaclib
        5, 78, -1, 94, -1, -1, 78, 94, 566, 210,  // fr: He4_Ca46_to_Ti50_reaclib, rr: Ti50_to_He4_Ca46_reaclib
        1, 79, -1, 80, -1, -1, 79, 80, 567, 174,  // fr: n_Ca47_to_Ca48_reaclib, rr: Ca48_to_n_Ca47_reaclib
        2, 79, -1, 86, -1, -1, 79, 86, 568, 188,  // fr: p_Ca47_to_Sc48_reaclib, rr: Sc48_to_p_Ca47_reaclib
        5, 79, -1, 95, -1, -1, 79, 95, 569, 212,  // fr: He4_Ca47_to_Ti51_reaclib, rr: Ti51_to_He4_Ca47_reaclib
        2, 80, -1, 87, -1, -1, 80, 87, 570, 190,  // fr: p_Ca48_to_Sc49_reaclib, rr: Sc49_to_p_Ca48_reaclib
        1, 81, -1, 82, -1, -1, 81, 82, 571, 177,  // fr: n_Sc43_to_Sc44_reaclib, rr: Sc44_to_n_Sc43_reaclib
        2, 81, -1, 88, -1, -1, 81, 88, 572, 191,  // fr: p_Sc43_to_Ti44_reaclib, rr: Ti44_to_p_Sc43_reaclib
        5, 81, -1, 97, -1, -1, 81, 97, 573, 216,  // fr: He4_Sc43_to_V47_reaclib, rr: V47_to_He4_Sc43_reaclib
        1, 82, -1, 83, -1, -1, 82, 83, 574, 180,  // fr: n_Sc44_to_Sc45_reaclib, rr: Sc45_to_n_Sc44_reaclib
        2, 82, -1, 89, -1, -1, 82, 89, 575, 194,  // fr: p_Sc44_to_Ti45_reaclib, rr: Ti45_to_p_Sc44_reaclib
        5, 82, -1, 98, -1, -1, 82, 98, 576, 219,  // fr: He4_Sc44_to_V48_reaclib, rr: V48_to_He4_Sc44_reaclib
        1, 83, -1, 84, -1, -1, 83, 84, 577, 183,  // fr: n_Sc45_to_Sc46_reaclib, rr: Sc46_to_n_Sc45_reaclib
        2, 83, -1, 90, -1, -1, 83, 90, 578, 197,  // fr: p_Sc45_to_Ti46_reaclib, rr: Ti46_to_p_Sc45_reaclib
        5, 83, -1, 99, -1, -1, 83, 99, 579, 222,  // fr: He4_Sc45_to_V49_reaclib, rr: V49_to_He4_Sc45_reaclib
        1, 84, -1, 85, -1, -1, 84, 85, 580, 185,  // fr: n_Sc46_to_Sc47_reaclib, rr: Sc47_to_n_Sc46_reaclib
        2, 84, -1, 91, -1, -1, 84, 91, 581, 200,  // fr: p_Sc46_to_Ti47_reaclib, rr: Ti47_to_p_Sc46_reaclib
        5, 84, -1, 100, -1, -1, 84, 100, 582, 225,  // fr: He4_Sc46_to_V50_reaclib, rr: V50_to_He4_Sc46_reaclib
        1, 85, -1, 86, -1, -1, 85, 86, 583, 187,  // fr: n_Sc47_to_Sc48_reaclib, rr: Sc48_to_n_Sc47_reaclib
        2, 85, -1, 92, -1, -1, 85, 92, 584, 203,  // fr: p_Sc47_to_Ti48_reaclib, rr: Ti48_to_p_Sc47_reaclib
        5, 85, -1, 101, -1, -1, 85, 101, 585, 228,  // fr: He4_Sc47_to_V51_reaclib, rr: V51_to_He4_Sc47_reaclib
        1, 86, -1, 87, -1, -1, 86, 87, 586, 189,  // fr: n_Sc48_to_Sc49_reaclib, rr: Sc49_to_n_Sc48_reaclib
        2, 86, -1, 93, -1, -1, 86, 93, 587, 206,  // fr: p_Sc48_to_Ti49_reaclib, rr: Ti49_to_p_Sc48_reaclib
        5, 86, -1, 102, -1, -1, 86, 102, 588, 231,  // fr: He4_Sc48_to_V52_reaclib, rr: V52_to_He4_Sc48_reaclib
        2, 87, -1, 94, -1, -1, 87, 94, 589, 209,  // fr: p_Sc49_to_Ti50_reaclib, rr: Ti50_to_p_Sc49_reaclib
        1, 88, -1, 89, -1, -1, 88, 89, 590, 193,  // fr: n_Ti44_to_Ti45_reaclib, rr: Ti45_to_n_Ti44_reaclib
        5, 88, -1, 103, -1, -1, 88, 103, 591, 233,  // fr: He4_Ti44_to_Cr48_reaclib, rr: Cr48_to_He4_Ti44_reaclib
        1, 89, -1, 90, -1, -1, 89, 90, 592, 196,  // fr: n_Ti45_to_Ti46_reaclib, rr: Ti46_to_n_Ti45_reaclib
        2, 89, -1, 96, -1, -1, 89, 96, 593, 213,  // fr: p_Ti45_to_V46_reaclib, rr: V46_to_p_Ti45_reaclib
        5, 89, -1, 104, -1, -1, 89, 104, 594, 236,  // fr: He4_Ti45_to_Cr49_reaclib, rr: Cr49_to_He4_Ti45_reaclib
        1, 90, -1, 91, -1, -1, 90, 91, 595, 199,  // fr: n_Ti46_to_Ti47_reaclib, rr: Ti47_to_n_Ti46_reaclib
        2, 90, -1, 97, -1, -1, 90, 97, 596, 215,  // fr: p_Ti46_to_V47_reaclib, rr: V47_to_p_Ti46_reaclib
        5, 90, -1, 105, -1, -1, 90, 105, 597, 239,  // fr: He4_Ti46_to_Cr50_reaclib, rr: Cr50_to_He4_Ti46_reaclib
        1, 91, -1, 92, -1, -1, 91, 92, 598, 202,  // fr: n_Ti47_to_Ti48_reaclib, rr: Ti48_to_n_Ti47_reaclib
        2, 91, -1, 98, -1, -1, 91, 98, 599, 218,  // fr: p_Ti47_to_V48_reaclib, rr: V48_to_p_Ti47_reaclib
        5, 91, -1, 106, -1, -1, 91, 106, 600, 242,  // fr: He4_Ti47_to_Cr51_reaclib, rr: Cr51_to_He4_Ti47_reaclib
        1, 92, -1, 93, -1, -1, 92, 93, 601, 205,  // fr: n_Ti48_to_Ti49_reaclib, rr: Ti49_to_n_Ti48_reaclib
        2, 92, -1, 99, -1, -1, 92, 99, 602, 221,  // fr: p_Ti48_to_V49_reaclib, rr: V49_to_p_Ti48_reaclib
        5, 92, -1, 107, -1, -1, 92, 107, 603, 245,  // fr: He4_Ti48_to_Cr52_reaclib, rr: Cr52_to_He4_Ti48_reaclib
        1, 93, -1, 94, -1, -1, 93, 94, 604, 208,  // fr: n_Ti49_to_Ti50_reaclib, rr: Ti50_to_n_Ti49_reaclib
        2, 93, -1, 100, -1, -1, 93, 100, 605, 224,  // fr: p_Ti49_to_V50_reaclib, rr: V50_to_p_Ti49_reaclib
        5, 93, -1, 108, -1, -1, 93, 108, 606, 248,  // fr: He4_Ti49_to_Cr53_reaclib, rr: Cr53_to_He4_Ti49_reaclib
        1, 94, -1, 95, -1, -1, 94, 95, 607, 211,  // fr: n_Ti50_to_Ti51_reaclib, rr: Ti51_to_n_Ti50_reaclib
        2, 94, -1, 101, -1, -1, 94, 101, 608, 227,  // fr: p_Ti50_to_V51_reaclib, rr: V51_to_p_Ti50_reaclib
        5, 94, -1, 109, -1, -1, 94, 109, 609, 250,  // fr: He4_Ti50_to_Cr54_reaclib, rr: Cr54_to_He4_Ti50_reaclib
        2, 95, -1, 102, -1, -1, 95, 102, 610, 230,  // fr: p_Ti51_to_V52_reaclib, rr: V52_to_p_Ti51_reaclib
        1, 96, -1, 97, -1, -1, 96, 97, 611, 214,  // fr: n_V46_to_V47_reaclib, rr: V47_to_n_V46_reaclib
        5, 96, -1, 110, -1, -1, 96, 110, 612, 252,  // fr: He4_V46_to_Mn50_reaclib, rr: Mn50_to_He4_V46_reaclib
        1, 97, -1, 98, -1, -1, 97, 98, 613, 217,  // fr: n_V47_to_V48_reaclib, rr: V48_to_n_V47_reaclib
        2, 97, -1, 103, -1, -1, 97, 103, 614, 232,  // fr: p_V47_to_Cr48_reaclib, rr: Cr48_to_p_V47_reaclib
        5, 97, -1, 111, -1, -1, 97, 111, 615, 255,  // fr: He4_V47_to_Mn51_reaclib, rr: Mn51_to_He4_V47_reaclib
        1, 98, -1, 99, -1, -1, 98, 99, 616, 220,  // fr: n_V48_to_V49_reaclib, rr: V49_to_n_V48_reaclib
        2, 98, -1, 104, -1, -1, 98, 104, 617, 235,  // fr: p_V48_to_Cr49_reaclib, rr: Cr49_to_p_V48_reaclib
        5, 98, -1, 112, -1, -1, 98, 112, 618, 258,  // fr: He4_V48_to_Mn52_reaclib, rr: Mn52_to_He4_V48_reaclib
        1, 99, -1, 100, -1, -1, 99, 100, 619, 223,  // fr: n_V49_to_V50_reaclib, rr: V50_to_n_V49_reaclib
        2, 99, -1, 105, -1, -1, 99, 105, 620, 238,  // fr: p_V49_to_Cr50_reaclib, rr: Cr50_to_p_V49_reaclib
        5, 99, -1, 113, -1, -1, 99, 113, 621, 261,  // fr: He4_V49_to_Mn53_reaclib, rr: Mn53_to_He4_V49_reaclib
        1, 100, -1, 101, -1, -1, 100, 101, 622, 226,  // fr: n_V50_to_V51_reaclib, rr: V51_to_n_V50_reaclib
        2, 100, -1, 106, -1, -1, 100, 106, 623, 241,  // fr: p_V50_to_Cr51_reaclib, rr: Cr51_to_p_V50_reaclib
        5, 100, -1, 114, -1, -1, 100, 114, 624, 264,  // fr: He4_V50_to_Mn54_reaclib, rr: Mn54_to_He4_V50_reaclib
        1, 101, -1, 102, -1, -1, 101, 102, 625, 229,  // fr: n_V51_to_V52_reaclib, rr: V52_to_n_V51_reaclib
        2, 101, -1, 107, -1, -1, 101, 107, 626, 244,  // fr: p_V51_to_Cr52_reaclib, rr: Cr52_to_p_V51_reaclib
        5, 101, -1, 115, -1, -1, 101, 115, 627, 267,  // fr: He4_V51_to_Mn55_reaclib, rr: Mn55_to_He4_V51_reaclib
        2, 102, -1, 108, -1, -1, 102, 108, 628, 247,  // fr: p_V52_to_Cr53_reaclib, rr: Cr53_to_p_V52_reaclib
        1, 103, -1, 104, -1, -1, 103, 104, 629, 234,  // fr: n_Cr48_to_Cr49_reaclib, rr: Cr49_to_n_Cr48_reaclib
        5, 103, -1, 116, -1, -1, 103, 116, 630, 269,  // fr: He4_Cr48_to_Fe52_reaclib, rr: Fe52_to_He4_Cr48_reaclib
        1, 104, -1, 105, -1, -1, 104, 105, 631, 237,  // fr: n_Cr49_to_Cr50_reaclib, rr: Cr50_to_n_Cr49_reaclib
        2, 104, -1, 110, -1, -1, 104, 110, 632, 251,  // fr: p_Cr49_to_Mn50_reaclib, rr: Mn50_to_p_Cr49_reaclib
        5, 104, -1, 117, -1, -1, 104, 117, 633, 272,  // fr: He4_Cr49_to_Fe53_reaclib, rr: Fe53_to_He4_Cr49_reaclib
        1, 105, -1, 106, -1, -1, 105, 106, 634, 240,  // fr: n_Cr50_to_Cr51_reaclib, rr: Cr51_to_n_Cr50_reaclib
        2, 105, -1, 111, -1, -1, 105, 111, 635, 254,  // fr: p_Cr50_to_Mn51_reaclib, rr: Mn51_to_p_Cr50_reaclib
        5, 105, -1, 118, -1, -1, 105, 118, 636, 275,  // fr: He4_Cr50_to_Fe54_reaclib, rr: Fe54_to_He4_Cr50_reaclib
        1, 106, -1, 107, -1, -1, 106, 107, 637, 243,  // fr: n_Cr51_to_Cr52_reaclib, rr: Cr52_to_n_Cr51_reaclib
        2, 106, -1, 112, -1, -1, 106, 112, 638, 257,  // fr: p_Cr51_to_Mn52_reaclib, rr: Mn52_to_p_Cr51_reaclib
        5, 106, -1, 119, -1, -1, 106, 119, 639, 278,  // fr: He4_Cr51_to_Fe55_reaclib, rr: Fe55_to_He4_Cr51_reaclib
        1, 107, -1, 108, -1, -1, 107, 108, 640, 246,  // fr: n_Cr52_to_Cr53_reaclib, rr: Cr53_to_n_Cr52_reaclib
        2, 107, -1, 113, -1, -1, 107, 113, 641, 260,  // fr: p_Cr52_to_Mn53_reaclib, rr: Mn53_to_p_Cr52_reaclib
        5, 107, -1, 120, -1, -1, 107, 120, 642, 281,  // fr: He4_Cr52_to_Fe56_reaclib, rr: Fe56_to_He4_Cr52_reaclib
        1, 108, -1, 109, -1, -1, 108, 109, 643, 249,  // fr: n_Cr53_to_Cr54_reaclib, rr: Cr54_to_n_Cr53_reaclib
        2, 108, -1, 114, -1, -1, 108, 114, 644, 263,  // fr: p_Cr53_to_Mn54_reaclib, rr: Mn54_to_p_Cr53_reaclib
        5, 108, -1, 121, -1, -1, 108, 121, 645, 283,  // fr: He4_Cr53_to_Fe57_reaclib, rr: Fe57_to_He4_Cr53_reaclib
        2, 109, -1, 115, -1, -1, 109, 115, 646, 266,  // fr: p_Cr54_to_Mn55_reaclib, rr: Mn55_to_p_Cr54_reaclib
        5, 109, -1, 122, -1, -1, 109, 122, 647, 285,  // fr: He4_Cr54_to_Fe58_reaclib, rr: Fe58_to_He4_Cr54_reaclib
        1, 110, -1, 111, -1, -1, 110, 111, 648, 253,  // fr: n_Mn50_to_Mn51_reaclib, rr: Mn51_to_n_Mn50_reaclib
        5, 110, -1, 124, -1, -1, 110, 124, 649, 289,  // fr: He4_Mn50_to_Co54_reaclib, rr: Co54_to_He4_Mn50_reaclib
        1, 111, -1, 112, -1, -1, 111, 112, 650, 256,  // fr: n_Mn51_to_Mn52_reaclib, rr: Mn52_to_n_Mn51_reaclib
        2, 111, -1, 116, -1, -1, 111, 116, 651, 268,  // fr: p_Mn51_to_Fe52_reaclib, rr: Fe52_to_p_Mn51_reaclib
        5, 111, -1, 125, -1, -1, 111, 125, 652, 292,  // fr: He4_Mn51_to_Co55_reaclib, rr: Co55_to_He4_Mn51_reaclib
        1, 112, -1, 113, -1, -1, 112, 113, 653, 259,  // fr: n_Mn52_to_Mn53_reaclib, rr: Mn53_to_n_Mn52_reaclib
        2, 112, -1, 117, -1, -1, 112, 117, 654, 271,  // fr: p_Mn52_to_Fe53_reaclib, rr: Fe53_to_p_Mn52_reaclib
        5, 112, -1, 126, -1, -1, 112, 126, 655, 295,  // fr: He4_Mn52_to_Co56_reaclib, rr: Co56_to_He4_Mn52_reaclib
        1, 113, -1, 114, -1, -1, 113, 114, 656, 262,  // fr: n_Mn53_to_Mn54_reaclib, rr: Mn54_to_n_Mn53_reaclib
        2, 113, -1, 118, -1, -1, 113, 118, 657, 274,  // fr: p_Mn53_to_Fe54_reaclib, rr: Fe54_to_p_Mn53_reaclib
        5, 113, -1, 127, -1, -1, 113, 127, 658, 298,  // fr: He4_Mn53_to_Co57_reaclib, rr: Co57_to_He4_Mn53_reaclib
        1, 114, -1, 115, -1, -1, 114, 115, 659, 265,  // fr: n_Mn54_to_Mn55_reaclib, rr: Mn55_to_n_Mn54_reaclib
        2, 114, -1, 119, -1, -1, 114, 119, 660, 277,  // fr: p_Mn54_to_Fe55_reaclib, rr: Fe55_to_p_Mn54_reaclib
        5, 114, -1, 128, -1, -1, 114, 128, 661, 301,  // fr: He4_Mn54_to_Co58_reaclib, rr: Co58_to_He4_Mn54_reaclib
        2, 115, -1, 120, -1, -1, 115, 120, 662, 280,  // fr: p_Mn55_to_Fe56_reaclib, rr: Fe56_to_p_Mn55_reaclib
        5, 115, -1, 129, -1, -1, 115, 129, 663, 304,  // fr: He4_Mn55_to_Co59_reaclib, rr: Co59_to_He4_Mn55_reaclib
        1, 116, -1, 117, -1, -1, 116, 117, 664, 270,  // fr: n_Fe52_to_Fe53_reaclib, rr: Fe53_to_n_Fe52_reaclib
        2, 116, -1, 123, -1, -1, 116, 123, 665, 286,  // fr: p_Fe52_to_Co53_reaclib, rr: Co53_to_p_Fe52_reaclib
        5, 116, -1, 130, -1, -1, 116, 130, 666, 306,  // fr: He4_Fe52_to_Ni56_reaclib, rr: Ni56_to_He4_Fe52_reaclib
        1, 117, -1, 118, -1, -1, 117, 118, 667, 273,  // fr: n_Fe53_to_Fe54_reaclib, rr: Fe54_to_n_Fe53_reaclib
        2, 117, -1, 124, -1, -1, 117, 124, 668, 288,  // fr: p_Fe53_to_Co54_reaclib, rr: Co54_to_p_Fe53_reaclib
        5, 117, -1, 131, -1, -1, 117, 131, 669, 309,  // fr: He4_Fe53_to_Ni57_reaclib, rr: Ni57_to_He4_Fe53_reaclib
        1, 118, -1, 119, -1, -1, 118, 119, 670, 276,  // fr: n_Fe54_to_Fe55_reaclib, rr: Fe55_to_n_Fe54_reaclib
        2, 118, -1, 125, -1, -1, 118, 125, 671, 291,  // fr: p_Fe54_to_Co55_reaclib, rr: Co55_to_p_Fe54_reaclib
        5, 118, -1, 132, -1, -1, 118, 132, 672, 312,  // fr: He4_Fe54_to_Ni58_reaclib, rr: Ni58_to_He4_Fe54_reaclib
        1, 119, -1, 120, -1, -1, 119, 120, 673, 279,  // fr: n_Fe55_to_Fe56_reaclib, rr: Fe56_to_n_Fe55_reaclib
        2, 119, -1, 126, -1, -1, 119, 126, 674, 294,  // fr: p_Fe55_to_Co56_reaclib, rr: Co56_to_p_Fe55_reaclib
        5, 119, -1, 133, -1, -1, 119, 133, 675, 315,  // fr: He4_Fe55_to_Ni59_reaclib, rr: Ni59_to_He4_Fe55_reaclib
        1, 120, -1, 121, -1, -1, 120, 121, 676, 282,  // fr: n_Fe56_to_Fe57_reaclib, rr: Fe57_to_n_Fe56_reaclib
        2, 120, -1, 127, -1, -1, 120, 127, 677, 297,  // fr: p_Fe56_to_Co57_reaclib, rr: Co57_to_p_Fe56_reaclib
        5, 120, -1, 134, -1, -1, 120, 134, 678, 318,  // fr: He4_Fe56_to_Ni60_reaclib, rr: Ni60_to_He4_Fe56_reaclib
        1, 121, -1, 122, -1, -1, 121, 122, 679, 284,  // fr: n_Fe57_to_Fe58_reaclib, rr: Fe58_to_n_Fe57_reaclib
        2, 121, -1, 128, -1, -1, 121, 128, 680, 300,  // fr: p_Fe57_to_Co58_reaclib, rr: Co58_to_p_Fe57_reaclib
        5, 121, -1, 135, -1, -1, 121, 135, 681, 320,  // fr: He4_Fe57_to_Ni61_reaclib, rr: Ni61_to_He4_Fe57_reaclib
        2, 122, -1, 129, -1, -1, 122, 129, 682, 303,  // fr: p_Fe58_to_Co59_reaclib, rr: Co59_to_p_Fe58_reaclib
        5, 122, -1, 136, -1, -1, 122, 136, 683, 322,  // fr: He4_Fe58_to_Ni62_reaclib, rr: Ni62_to_He4_Fe58_reaclib
        1, 123, -1, 124, -1, -1, 123, 124, 684, 287,  // fr: n_Co53_to_Co54_reaclib, rr: Co54_to_n_Co53_reaclib
        5, 123, -1, 139, -1, -1, 123, 139, 685, 326,  // fr: He4_Co53_to_Cu57_reaclib, rr: Cu57_to_He4_Co53_reaclib
        1, 124, -1, 125, -1, -1, 124, 125, 686, 290,  // fr: n_Co54_to_Co55_reaclib, rr: Co55_to_n_Co54_reaclib
        5, 124, -1, 140, -1, -1, 124, 140, 687, 329,  // fr: He4_Co54_to_Cu58_reaclib, rr: Cu58_to_He4_Co54_reaclib
        1, 125, -1, 126, -1, -1, 125, 126, 688, 293,  // fr: n_Co55_to_Co56_reaclib, rr: Co56_to_n_Co55_reaclib
        2, 125, -1, 130, -1, -1, 125, 130, 689, 305,  // fr: p_Co55_to_Ni56_reaclib, rr: Ni56_to_p_Co55_reaclib
        5, 125, -1, 141, -1, -1, 125, 141, 690, 332,  // fr: He4_Co55_to_Cu59_reaclib, rr: Cu59_to_He4_Co55_reaclib
        1, 126, -1, 127, -1, -1, 126, 127, 691, 296,  // fr: n_Co56_to_Co57_reaclib, rr: Co57_to_n_Co56_reaclib
        2, 126, -1, 131, -1, -1, 126, 131, 692, 308,  // fr: p_Co56_to_Ni57_reaclib, rr: Ni57_to_p_Co56_reaclib
        5, 126, -1, 142, -1, -1, 126, 142, 693, 335,  // fr: He4_Co56_to_Cu60_reaclib, rr: Cu60_to_He4_Co56_reaclib
        1, 127, -1, 128, -1, -1, 127, 128, 694, 299,  // fr: n_Co57_to_Co58_reaclib, rr: Co58_to_n_Co57_reaclib
        2, 127, -1, 132, -1, -1, 127, 132, 695, 311,  // fr: p_Co57_to_Ni58_reaclib, rr: Ni58_to_p_Co57_reaclib
        5, 127, -1, 143, -1, -1, 127, 143, 696, 338,  // fr: He4_Co57_to_Cu61_reaclib, rr: Cu61_to_He4_Co57_reaclib
        1, 128, -1, 129, -1, -1, 128, 129, 697, 302,  // fr: n_Co58_to_Co59_reaclib, rr: Co59_to_n_Co58_reaclib
        2, 128, -1, 133, -1, -1, 128, 133, 698, 314,  // fr: p_Co58_to_Ni59_reaclib, rr: Ni59_to_p_Co58_reaclib
        5, 128, -1, 144, -1, -1, 128, 144, 699, 341,  // fr: He4_Co58_to_Cu62_reaclib, rr: Cu62_to_He4_Co58_reaclib
        2, 129, -1, 134, -1, -1, 129, 134, 700, 317,  // fr: p_Co59_to_Ni60_reaclib, rr: Ni60_to_p_Co59_reaclib
        5, 129, -1, 145, -1, -1, 129, 145, 701, 344,  // fr: He4_Co59_to_Cu63_reaclib, rr: Cu63_to_He4_Co59_reaclib
        1, 130, -1, 131, -1, -1, 130, 131, 702, 307,  // fr: n_Ni56_to_Ni57_reaclib, rr: Ni57_to_n_Ni56_reaclib
        2, 130, -1, 139, -1, -1, 130, 139, 703, 325,  // fr: p_Ni56_to_Cu57_reaclib, rr: Cu57_to_p_Ni56_reaclib
        5, 130, -1, 149, -1, -1, 130, 149, 704, 353,  // fr: He4_Ni56_to_Zn60_reaclib, rr: Zn60_to_He4_Ni56_reaclib
        1, 131, -1, 132, -1, -1, 131, 132, 705, 310,  // fr: n_Ni57_to_Ni58_reaclib, rr: Ni58_to_n_Ni57_reaclib
        2, 131, -1, 140, -1, -1, 131, 140, 706, 328,  // fr: p_Ni57_to_Cu58_reaclib, rr: Cu58_to_p_Ni57_reaclib
        5, 131, -1, 150, -1, -1, 131, 150, 707, 356,  // fr: He4_Ni57_to_Zn61_reaclib, rr: Zn61_to_He4_Ni57_reaclib
        1, 132, -1, 133, -1, -1, 132, 133, 708, 313,  // fr: n_Ni58_to_Ni59_reaclib, rr: Ni59_to_n_Ni58_reaclib
        2, 132, -1, 141, -1, -1, 132, 141, 709, 331,  // fr: p_Ni58_to_Cu59_reaclib, rr: Cu59_to_p_Ni58_reaclib
        5, 132, -1, 151, -1, -1, 132, 151, 710, 359,  // fr: He4_Ni58_to_Zn62_reaclib, rr: Zn62_to_He4_Ni58_reaclib
        1, 133, -1, 134, -1, -1, 133, 134, 711, 316,  // fr: n_Ni59_to_Ni60_reaclib, rr: Ni60_to_n_Ni59_reaclib
        2, 133, -1, 142, -1, -1, 133, 142, 712, 334,  // fr: p_Ni59_to_Cu60_reaclib, rr: Cu60_to_p_Ni59_reaclib
        5, 133, -1, 152, -1, -1, 133, 152, 713, 362,  // fr: He4_Ni59_to_Zn63_reaclib, rr: Zn63_to_He4_Ni59_reaclib
        1, 134, -1, 135, -1, -1, 134, 135, 714, 319,  // fr: n_Ni60_to_Ni61_reaclib, rr: Ni61_to_n_Ni60_reaclib
        2, 134, -1, 143, -1, -1, 134, 143, 715, 337,  // fr: p_Ni60_to_Cu61_reaclib, rr: Cu61_to_p_Ni60_reaclib
        5, 134, -1, 153, -1, -1, 134, 153, 716, 365,  // fr: He4_Ni60_to_Zn64_reaclib, rr: Zn64_to_He4_Ni60_reaclib
        1, 135, -1, 136, -1, -1, 135, 136, 717, 321,  // fr: n_Ni61_to_Ni62_reaclib, rr: Ni62_to_n_Ni61_reaclib
        2, 135, -1, 144, -1, -1, 135, 144, 718, 340,  // fr: p_Ni61_to_Cu62_reaclib, rr: Cu62_to_p_Ni61_reaclib
        5, 135, -1, 154, -1, -1, 135, 154, 719, 368,  // fr: He4_Ni61_to_Zn65_reaclib, rr: Zn65_to_He4_Ni61_reaclib
        1, 136, -1, 137, -1, -1, 136, 137, 720, 323,  // fr: n_Ni62_to_Ni63_reaclib, rr: Ni63_to_n_Ni62_reaclib
        2, 136, -1, 145, -1, -1, 136, 145, 721, 343,  // fr: p_Ni62_to_Cu63_reaclib, rr: Cu63_to_p_Ni62_reaclib
        5, 136, -1, 155, -1, -1, 136, 155, 722, 371,  // fr: He4_Ni62_to_Zn66_reaclib, rr: Zn66_to_He4_Ni62_reaclib
        1, 137, -1, 138, -1, -1, 137, 138, 723, 324,  // fr: n_Ni63_to_Ni64_reaclib, rr: Ni64_to_n_Ni63_reaclib
        2, 137, -1, 146, -1, -1, 137, 146, 724, 346,  // fr: p_Ni63_to_Cu64_reaclib, rr: Cu64_to_p_Ni63_reaclib
        2, 138, -1, 147, -1, -1, 138, 147, 725, 348,  // fr: p_Ni64_to_Cu65_reaclib, rr: Cu65_to_p_Ni64_reaclib
        1, 139, -1, 140, -1, -1, 139, 140, 726, 327,  // fr: n_Cu57_to_Cu58_reaclib, rr: Cu58_to_n_Cu57_reaclib
        1, 140, -1, 141, -1, -1, 140, 141, 727, 330,  // fr: n_Cu58_to_Cu59_reaclib, rr: Cu59_to_n_Cu58_reaclib
        2, 140, -1, 148, -1, -1, 140, 148, 728, 349,  // fr: p_Cu58_to_Zn59_reaclib, rr: Zn59_to_p_Cu58_reaclib
        5, 140, -1, 156, -1, -1, 140, 156, 729, 373,  // fr: He4_Cu58_to_Ga62_reaclib, rr: Ga62_to_He4_Cu58_reaclib
        1, 141, -1, 142, -1, -1, 141, 142, 730, 333,  // fr: n_Cu59_to_Cu60_reaclib, rr: Cu60_to_n_Cu59_reaclib
        2, 141, -1, 149, -1, -1, 141, 149, 731, 352,  // fr: p_Cu59_to_Zn60_reaclib, rr: Zn60_to_p_Cu59_reaclib
        5, 141, -1, 157, -1, -1, 141, 157, 732, 376,  // fr: He4_Cu59_to_Ga63_reaclib, rr: Ga63_to_He4_Cu59_reaclib
        1, 142, -1, 143, -1, -1, 142, 143, 733, 336,  // fr: n_Cu60_to_Cu61_reaclib, rr: Cu61_to_n_Cu60_reaclib
        2, 142, -1, 150, -1, -1, 142, 150, 734, 355,  // fr: p_Cu60_to_Zn61_reaclib, rr: Zn61_to_p_Cu60_reaclib
        5, 142, -1, 158, -1, -1, 142, 158, 735, 379,  // fr: He4_Cu60_to_Ga64_reaclib, rr: Ga64_to_He4_Cu60_reaclib
        1, 143, -1, 144, -1, -1, 143, 144, 736, 339,  // fr: n_Cu61_to_Cu62_reaclib, rr: Cu62_to_n_Cu61_reaclib
        2, 143, -1, 151, -1, -1, 143, 151, 737, 358,  // fr: p_Cu61_to_Zn62_reaclib, rr: Zn62_to_p_Cu61_reaclib
        1, 144, -1, 145, -1, -1, 144, 145, 738, 342,  // fr: n_Cu62_to_Cu63_reaclib, rr: Cu63_to_n_Cu62_reaclib
        2, 144, -1, 152, -1, -1, 144, 152, 739, 361,  // fr: p_Cu62_to_Zn63_reaclib, rr: Zn63_to_p_Cu62_reaclib
        1, 145, -1, 146, -1, -1, 145, 146, 740, 345,  // fr: n_Cu63_to_Cu64_reaclib, rr: Cu64_to_n_Cu63_reaclib
        2, 145, -1, 153, -1, -1, 145, 153, 741, 364,  // fr: p_Cu63_to_Zn64_reaclib, rr: Zn64_to_p_Cu63_reaclib
        1, 146, -1, 147, -1, -1, 146, 147, 742, 347,  // fr: n_Cu64_to_Cu65_reaclib, rr: Cu65_to_n_Cu64_reaclib
        2, 146, -1, 154, -1, -1, 146, 154, 743, 367,  // fr: p_Cu64_to_Zn65_reaclib, rr: Zn65_to_p_Cu64_reaclib
        2, 147, -1, 155, -1, -1, 147, 155, 744, 370,  // fr: p_Cu65_to_Zn66_reaclib, rr: Zn66_to_p_Cu65_reaclib
        1, 148, -1, 149, -1, -1, 148, 149, 745, 351,  // fr: n_Zn59_to_Zn60_reaclib, rr: Zn60_to_n_Zn59_reaclib
        5, 148, -1, 159, -1, -1, 148, 159, 746, 381,  // fr: He4_Zn59_to_Ge63_reaclib, rr: Ge63_to_He4_Zn59_reaclib
        1, 149, -1, 150, -1, -1, 149, 150, 747, 354,  // fr: n_Zn60_to_Zn61_reaclib, rr: Zn61_to_n_Zn60_reaclib
        5, 149, -1, 160, -1, -1, 149, 160, 748, 384,  // fr: He4_Zn60_to_Ge64_reaclib, rr: Ge64_to_He4_Zn60_reaclib
        1, 150, -1, 151, -1, -1, 150, 151, 749, 357,  // fr: n_Zn61_to_Zn62_reaclib, rr: Zn62_to_n_Zn61_reaclib
        2, 150, -1, 156, -1, -1, 150, 156, 750, 372,  // fr: p_Zn61_to_Ga62_reaclib, rr: Ga62_to_p_Zn61_reaclib
        1, 151, -1, 152, -1, -1, 151, 152, 751, 360,  // fr: n_Zn62_to_Zn63_reaclib, rr: Zn63_to_n_Zn62_reaclib
        2, 151, -1, 157, -1, -1, 151, 157, 752, 375,  // fr: p_Zn62_to_Ga63_reaclib, rr: Ga63_to_p_Zn62_reaclib
        1, 152, -1, 153, -1, -1, 152, 153, 753, 363,  // fr: n_Zn63_to_Zn64_reaclib, rr: Zn64_to_n_Zn63_reaclib
        2, 152, -1, 158, -1, -1, 152, 158, 754, 378,  // fr: p_Zn63_to_Ga64_reaclib, rr: Ga64_to_p_Zn63_reaclib
        1, 153, -1, 154, -1, -1, 153, 154, 755, 366,  // fr: n_Zn64_to_Zn65_reaclib, rr: Zn65_to_n_Zn64_reaclib
        1, 154, -1, 155, -1, -1, 154, 155, 756, 369,  // fr: n_Zn65_to_Zn66_reaclib, rr: Zn66_to_n_Zn65_reaclib
        1, 156, -1, 157, -1, -1, 156, 157, 757, 374,  // fr: n_Ga62_to_Ga63_reaclib, rr: Ga63_to_n_Ga62_reaclib
        2, 156, -1, 159, -1, -1, 156, 159, 758, 380,  // fr: p_Ga62_to_Ge63_reaclib, rr: Ge63_to_p_Ga62_reaclib
        1, 157, -1, 158, -1, -1, 157, 158, 759, 377,  // fr: n_Ga63_to_Ga64_reaclib, rr: Ga64_to_n_Ga63_reaclib
        2, 157, -1, 160, -1, -1, 157, 160, 760, 383,  // fr: p_Ga63_to_Ge64_reaclib, rr: Ge64_to_p_Ga63_reaclib
        1, 159, -1, 160, -1, -1, 159, 160, 761, 382,  // fr: n_Ge63_to_Ge64_reaclib, rr: Ge64_to_n_Ge63_reaclib
        3, 4, -1, 2, 5, -1, 3, 4, 764, 766,  // fr: d_He3_to_p_He4_reaclib, rr: p_He4_to_d_He3_reaclib
        2, 6, -1, 5, 4, -1, 6, 4, 769, 765,  // fr: p_Li6_to_He4_He3_reaclib, rr: He4_He3_to_p_Li6_reaclib
        2, 7, -1, 5, 5, -1, 7, -1, 775, 768,  // fr: p_Li7_to_He4_He4_reaclib, rr: He4_He4_to_p_Li7_reaclib
        1, 8, -1, 2, 7, -1, 8, 7, 777, 773,  // fr: n_Be7_to_p_Li7_reaclib, rr: p_Li7_to_n_Be7_reaclib
        1, 8, -1, 5, 5, -1, 8, -1, 779, 767,  // fr: n_Be7_to_He4_He4_reaclib, rr: He4_He4_to_n_Be7_reaclib
        2, 9, -1, 5, 6, -1, 9, 6, 781, 772,  // fr: p_Be9_to_He4_Li6_reaclib, rr: He4_Li6_to_p_Be9_reaclib
        5, 9, -1, 1, 13, -1, 9, 13, 782, 789,  // fr: He4_Be9_to_n_C12_reaclib, rr: n_C12_to_He4_Be9_reaclib
        1, 11, -1, 5, 7, -1, 11, 7, 783, 776,  // fr: n_B10_to_He4_Li7_reaclib, rr: He4_Li7_to_n_B10_reaclib
        2, 11, -1, 5, 8, -1, 11, 8, 784, 780,  // fr: p_B10_to_He4_Be7_reaclib, rr: He4_Be7_to_p_B10_reaclib
        5, 11, -1, 1, 16, -1, 11, 16, 785, 804,  // fr: He4_B10_to_n_N13_reaclib, rr: n_N13_to_He4_B10_reaclib
        5, 11, -1, 2, 14, -1, 11, 14, 786, 796,  // fr: He4_B10_to_p_C13_reaclib, rr: p_C13_to_He4_B10_reaclib
        5, 12, -1, 1, 17, -1, 12, 17, 787, 808,  // fr: He4_B11_to_n_N14_reaclib, rr: n_N14_to_He4_B11_reaclib
        5, 12, -1, 2, 15, -1, 12, 15, 788, 800,  // fr: He4_B11_to_p_C14_reaclib, rr: p_C14_to_He4_B11_reaclib
        5, 14, -1, 1, 21, -1, 14, 21, 798, 823,  // fr: He4_C13_to_n_O16_reaclib, rr: n_O16_to_He4_C13_reaclib
        1, 16, -1, 2, 14, -1, 16, 14, 803, 795,  // fr: n_N13_to_p_C13_reaclib, rr: p_C13_to_n_N13_reaclib
        5, 16, -1, 2, 21, -1, 16, 21, 805, 824,  // fr: He4_N13_to_p_O16_reaclib, rr: p_O16_to_He4_N13_reaclib
        1, 17, -1, 2, 15, -1, 17, 15, 806, 799,  // fr: n_N14_to_p_C14_reaclib, rr: p_C14_to_n_N14_reaclib
        2, 18, -1, 5, 13, -1, 18, 13, 814, 791,  // fr: p_N15_to_He4_C12_reaclib, rr: He4_C12_to_p_N15_reaclib
        1, 19, -1, 2, 17, -1, 19, 17, 817, 809,  // fr: n_O14_to_p_N14_reaclib, rr: p_N14_to_n_O14_reaclib
        5, 19, -1, 2, 24, -1, 19, 24, 818, 840,  // fr: He4_O14_to_p_F17_reaclib, rr: p_F17_to_He4_O14_reaclib
        1, 20, -1, 2, 18, -1, 20, 18, 819, 813,  // fr: n_O15_to_p_N15_reaclib, rr: p_N15_to_n_O15_reaclib
        1, 20, -1, 5, 13, -1, 20, 13, 820, 790,  // fr: n_O15_to_He4_C12_reaclib, rr: He4_C12_to_n_O15_reaclib
        1, 22, -1, 5, 15, -1, 22, 15, 831, 802,  // fr: n_O17_to_He4_C14_reaclib, rr: He4_C14_to_n_O17_reaclib
        2, 22, -1, 5, 17, -1, 22, 17, 833, 811,  // fr: p_O17_to_He4_N14_reaclib, rr: He4_N14_to_p_O17_reaclib
        5, 22, -1, 1, 29, -1, 22, 29, 834, 858,  // fr: He4_O17_to_n_Ne20_reaclib, rr: n_Ne20_to_He4_O17_reaclib
        2, 23, -1, 5, 18, -1, 23, 18, 836, 816,  // fr: p_O18_to_He4_N15_reaclib, rr: He4_N15_to_p_O18_reaclib
        1, 24, -1, 2, 22, -1, 24, 22, 838, 832,  // fr: n_F17_to_p_O17_reaclib, rr: p_O17_to_n_F17_reaclib
        1, 24, -1, 5, 17, -1, 24, 17, 839, 810,  // fr: n_F17_to_He4_N14_reaclib, rr: He4_N14_to_n_F17_reaclib
        5, 24, -1, 2, 29, -1, 24, 29, 841, 859,  // fr: He4_F17_to_p_Ne20_reaclib, rr: p_Ne20_to_He4_F17_reaclib
        1, 25, -1, 2, 23, -1, 25, 23, 842, 835,  // fr: n_F18_to_p_O18_reaclib, rr: p_O18_to_n_F18_reaclib
        1, 25, -1, 5, 18, -1, 25, 18, 843, 815,  // fr: n_F18_to_He4_N15_reaclib, rr: He4_N15_to_n_F18_reaclib
        2, 25, -1, 5, 20, -1, 25, 20, 845, 822,  // fr: p_F18_to_He4_O15_reaclib, rr: He4_O15_to_p_F18_reaclib
        5, 25, -1, 2, 30, -1, 25, 30, 847, 867,  // fr: He4_F18_to_p_Ne21_reaclib, rr: p_Ne21_to_He4_F18_reaclib
        2, 26, -1, 5, 21, -1, 26, 21, 849, 826,  // fr: p_F19_to_He4_O16_reaclib, rr: He4_O16_to_p_F19_reaclib
        5, 26, -1, 2, 31, -1, 26, 31, 851, 870,  // fr: He4_F19_to_p_Ne22_reaclib, rr: p_Ne22_to_He4_F19_reaclib
        1, 27, -1, 2, 25, -1, 27, 25, 852, 844,  // fr: n_Ne18_to_p_F18_reaclib, rr: p_F18_to_n_Ne18_reaclib
        1, 27, -1, 5, 20, -1, 27, 20, 853, 821,  // fr: n_Ne18_to_He4_O15_reaclib, rr: He4_O15_to_n_Ne18_reaclib
        5, 27, -1, 2, 32, -1, 27, 32, 854, 874,  // fr: He4_Ne18_to_p_Na21_reaclib, rr: p_Na21_to_He4_Ne18_reaclib
        1, 28, -1, 2, 26, -1, 28, 26, 855, 848,  // fr: n_Ne19_to_p_F19_reaclib, rr: p_F19_to_n_Ne19_reaclib
        1, 28, -1, 5, 21, -1, 28, 21, 856, 825,  // fr: n_Ne19_to_He4_O16_reaclib, rr: He4_O16_to_n_Ne19_reaclib
        5, 28, -1, 2, 33, -1, 28, 33, 857, 878,  // fr: He4_Ne19_to_p_Na22_reaclib, rr: p_Na22_to_He4_Ne19_reaclib
        1, 30, -1, 5, 23, -1, 30, 23, 865, 837,  // fr: n_Ne21_to_He4_O18_reaclib, rr: He4_O18_to_n_Ne21_reaclib
        5, 30, -1, 1, 36, -1, 30, 36, 868, 890,  // fr: He4_Ne21_to_n_Mg24_reaclib, rr: n_Mg24_to_He4_Ne21_reaclib
        1, 32, -1, 2, 30, -1, 32, 30, 872, 866,  // fr: n_Na21_to_p_Ne21_reaclib, rr: p_Ne21_to_n_Na21_reaclib
        1, 32, -1, 5, 25, -1, 32, 25, 873, 846,  // fr: n_Na21_to_He4_F18_reaclib, rr: He4_F18_to_n_Na21_reaclib
        5, 32, -1, 2, 36, -1, 32, 36, 875, 891,  // fr: He4_Na21_to_p_Mg24_reaclib, rr: p_Mg24_to_He4_Na21_reaclib
        1, 33, -1, 2, 31, -1, 33, 31, 876, 869,  // fr: n_Na22_to_p_Ne22_reaclib, rr: p_Ne22_to_n_Na22_reaclib
        1, 33, -1, 5, 26, -1, 33, 26, 877, 850,  // fr: n_Na22_to_He4_F19_reaclib, rr: He4_F19_to_n_Na22_reaclib
        5, 33, -1, 2, 37, -1, 33, 37, 880, 896,  // fr: He4_Na22_to_p_Mg25_reaclib, rr: p_Mg25_to_He4_Na22_reaclib
        2, 34, -1, 5, 29, -1, 34, 29, 882, 861,  // fr: p_Na23_to_He4_Ne20_reaclib, rr: He4_Ne20_to_p_Na23_reaclib
        5, 34, -1, 2, 38, -1, 34, 38, 885, 899,  // fr: He4_Na23_to_p_Mg26_reaclib, rr: p_Mg26_to_He4_Na23_reaclib
        1, 35, -1, 2, 34, -1, 35, 34, 886, 881,  // fr: n_Mg23_to_p_Na23_reaclib, rr: p_Na23_to_n_Mg23_reaclib
        1, 35, -1, 5, 29, -1, 35, 29, 887, 860,  // fr: n_Mg23_to_He4_Ne20_reaclib, rr: He4_Ne20_to_n_Mg23_reaclib
        5, 35, -1, 2, 40, -1, 35, 40, 889, 906,  // fr: He4_Mg23_to_p_Al26_reaclib, rr: p_Al26_to_He4_Mg23_reaclib
        1, 37, -1, 5, 31, -1, 37, 31, 894, 871,  // fr: n_Mg25_to_He4_Ne22_reaclib, rr: He4_Ne22_to_n_Mg25_reaclib
        5, 37, -1, 1, 42, -1, 37, 42, 897, 913,  // fr: He4_Mg25_to_n_Si28_reaclib, rr: n_Si28_to_He4_Mg25_reaclib
        5, 38, -1, 1, 43, -1, 38, 43, 900, 918,  // fr: He4_Mg26_to_n_Si29_reaclib, rr: n_Si29_to_He4_Mg26_reaclib
        1, 39, -1, 2, 37, -1, 39, 37, 901, 895,  // fr: n_Al25_to_p_Mg25_reaclib, rr: p_Mg25_to_n_Al25_reaclib
        1, 39, -1, 5, 33, -1, 39, 33, 902, 879,  // fr: n_Al25_to_He4_Na22_reaclib, rr: He4_Na22_to_n_Al25_reaclib
        5, 39, -1, 2, 42, -1, 39, 42, 903, 914,  // fr: He4_Al25_to_p_Si28_reaclib, rr: p_Si28_to_He4_Al25_reaclib
        1, 40, -1, 2, 38, -1, 40, 38, 904, 898,  // fr: n_Al26_to_p_Mg26_reaclib, rr: p_Mg26_to_n_Al26_reaclib
        1, 40, -1, 5, 34, -1, 40, 34, 905, 884,  // fr: n_Al26_to_He4_Na23_reaclib, rr: He4_Na23_to_n_Al26_reaclib
        5, 40, -1, 2, 43, -1, 40, 43, 908, 920,  // fr: He4_Al26_to_p_Si29_reaclib, rr: p_Si29_to_He4_Al26_reaclib
        2, 41, -1, 5, 36, -1, 41, 36, 909, 892,  // fr: p_Al27_to_He4_Mg24_reaclib, rr: He4_Mg24_to_p_Al27_reaclib
        5, 41, -1, 2, 44, -1, 41, 44, 912, 924,  // fr: He4_Al27_to_p_Si30_reaclib, rr: p_Si30_to_He4_Al27_reaclib
        2, 45, -1, 1, 49, -1, 45, 49, 927, 938,  // fr: p_Si31_to_n_P31_reaclib, rr: n_P31_to_p_Si31_reaclib
        5, 45, -1, 1, 54, -1, 45, 54, 928, 963,  // fr: He4_Si31_to_n_S34_reaclib, rr: n_S34_to_He4_Si31_reaclib
        1, 47, -1, 2, 43, -1, 47, 43, 931, 919,  // fr: n_P29_to_p_Si29_reaclib, rr: p_Si29_to_n_P29_reaclib
        1, 47, -1, 5, 40, -1, 47, 40, 932, 907,  // fr: n_P29_to_He4_Al26_reaclib, rr: He4_Al26_to_n_P29_reaclib
        5, 47, -1, 2, 52, -1, 47, 52, 933, 955,  // fr: He4_P29_to_p_S32_reaclib, rr: p_S32_to_He4_P29_reaclib
        1, 48, -1, 2, 44, -1, 48, 44, 934, 923,  // fr: n_P30_to_p_Si30_reaclib, rr: p_Si30_to_n_P30_reaclib
        1, 48, -1, 5, 41, -1, 48, 41, 935, 911,  // fr: n_P30_to_He4_Al27_reaclib, rr: He4_Al27_to_n_P30_reaclib
        5, 48, -1, 2, 53, -1, 48, 53, 937, 960,  // fr: He4_P30_to_p_S33_reaclib, rr: p_S33_to_He4_P30_reaclib
        2, 49, -1, 5, 42, -1, 49, 42, 939, 915,  // fr: p_P31_to_He4_Si28_reaclib, rr: He4_Si28_to_p_P31_reaclib
        5, 49, -1, 2, 54, -1, 49, 54, 943, 965,  // fr: He4_P31_to_p_S34_reaclib, rr: p_S34_to_He4_P31_reaclib
        1, 50, -1, 2, 46, -1, 50, 46, 944, 929,  // fr: n_P32_to_p_Si32_reaclib, rr: p_Si32_to_n_P32_reaclib
        2, 50, -1, 1, 52, -1, 50, 52, 945, 953,  // fr: p_P32_to_n_S32_reaclib, rr: n_S32_to_p_P32_reaclib
        2, 50, -1, 5, 43, -1, 50, 43, 946, 922,  // fr: p_P32_to_He4_Si29_reaclib, rr: He4_Si29_to_p_P32_reaclib
        2, 51, -1, 5, 44, -1, 51, 44, 950, 926,  // fr: p_P33_to_He4_Si30_reaclib, rr: He4_Si30_to_p_P33_reaclib
        1, 52, -1, 5, 43, -1, 52, 43, 954, 921,  // fr: n_S32_to_He4_Si29_reaclib, rr: He4_Si29_to_n_S32_reaclib
        1, 53, -1, 2, 51, -1, 53, 51, 957, 949,  // fr: n_S33_to_p_P33_reaclib, rr: p_P33_to_n_S33_reaclib
        1, 53, -1, 5, 44, -1, 53, 44, 958, 925,  // fr: n_S33_to_He4_Si30_reaclib, rr: He4_Si30_to_n_S33_reaclib
        1, 55, -1, 5, 46, -1, 55, 46, 968, 930,  // fr: n_S35_to_He4_Si32_reaclib, rr: He4_Si32_to_n_S35_reaclib
        2, 55, -1, 5, 50, -1, 55, 50, 970, 948,  // fr: p_S35_to_He4_P32_reaclib, rr: He4_P32_to_p_S35_reaclib
        5, 55, -1, 1, 64, -1, 55, 64, 971, 1007,  // fr: He4_S35_to_n_Ar38_reaclib, rr: n_Ar38_to_He4_S35_reaclib
        2, 56, -1, 5, 51, -1, 56, 51, 973, 952,  // fr: p_S36_to_He4_P33_reaclib, rr: He4_P33_to_p_S36_reaclib
        1, 57, -1, 2, 53, -1, 57, 53, 975, 959,  // fr: n_Cl33_to_p_S33_reaclib, rr: p_S33_to_n_Cl33_reaclib
        1, 57, -1, 5, 48, -1, 57, 48, 976, 936,  // fr: n_Cl33_to_He4_P30_reaclib, rr: He4_P30_to_n_Cl33_reaclib
        5, 57, -1, 2, 62, -1, 57, 62, 977, 999,  // fr: He4_Cl33_to_p_Ar36_reaclib, rr: p_Ar36_to_He4_Cl33_reaclib
        1, 58, -1, 2, 54, -1, 58, 54, 978, 964,  // fr: n_Cl34_to_p_S34_reaclib, rr: p_S34_to_n_Cl34_reaclib
        1, 58, -1, 5, 49, -1, 58, 49, 979, 942,  // fr: n_Cl34_to_He4_P31_reaclib, rr: He4_P31_to_n_Cl34_reaclib
        5, 58, -1, 2, 63, -1, 58, 63, 981, 1004,  // fr: He4_Cl34_to_p_Ar37_reaclib, rr: p_Ar37_to_He4_Cl34_reaclib
        1, 59, -1, 2, 55, -1, 59, 55, 982, 969,  // fr: n_Cl35_to_p_S35_reaclib, rr: p_S35_to_n_Cl35_reaclib
        1, 59, -1, 5, 50, -1, 59, 50, 983, 947,  // fr: n_Cl35_to_He4_P32_reaclib, rr: He4_P32_to_n_Cl35_reaclib
        2, 59, -1, 5, 52, -1, 59, 52, 984, 956,  // fr: p_Cl35_to_He4_S32_reaclib, rr: He4_S32_to_p_Cl35_reaclib
        5, 59, -1, 2, 64, -1, 59, 64, 986, 1009,  // fr: He4_Cl35_to_p_Ar38_reaclib, rr: p_Ar38_to_He4_Cl35_reaclib
        1, 60, -1, 2, 56, -1, 60, 56, 987, 972,  // fr: n_Cl36_to_p_S36_reaclib, rr: p_S36_to_n_Cl36_reaclib
        1, 60, -1, 5, 51, -1, 60, 51, 988, 951,  // fr: n_Cl36_to_He4_P33_reaclib, rr: He4_P33_to_n_Cl36_reaclib
        2, 60, -1, 5, 53, -1, 60, 53, 990, 962,  // fr: p_Cl36_to_He4_S33_reaclib, rr: He4_S33_to_p_Cl36_reaclib
        2, 61, -1, 5, 54, -1, 61, 54, 994, 967,  // fr: p_Cl37_to_He4_S34_reaclib, rr: He4_S34_to_p_Cl37_reaclib
        1, 62, -1, 2, 60, -1, 62, 60, 997, 989,  // fr: n_Ar36_to_p_Cl36_reaclib, rr: p_Cl36_to_n_Ar36_reaclib
        1, 62, -1, 5, 53, -1, 62, 53, 998, 961,  // fr: n_Ar36_to_He4_S33_reaclib, rr: He4_S33_to_n_Ar36_reaclib
        1, 63, -1, 2, 61, -1, 63, 61, 1001, 993,  // fr: n_Ar37_to_p_Cl37_reaclib, rr: p_Cl37_to_n_Ar37_reaclib
        1, 63, -1, 5, 54, -1, 63, 54, 1002, 966,  // fr: n_Ar37_to_He4_S34_reaclib, rr: He4_S34_to_n_Ar37_reaclib
        1, 65, -1, 5, 56, -1, 65, 56, 1012, 974,  // fr: n_Ar39_to_He4_S36_reaclib, rr: He4_S36_to_n_Ar39_reaclib
        2, 65, -1, 5, 60, -1, 65, 60, 1014, 992,  // fr: p_Ar39_to_He4_Cl36_reaclib, rr: He4_Cl36_to_p_Ar39_reaclib
        2, 66, -1, 5, 61, -1, 66, 61, 1017, 996,  // fr: p_Ar40_to_He4_Cl37_reaclib, rr: He4_Cl37_to_p_Ar40_reaclib
        1, 67, -1, 2, 63, -1, 67, 63, 1019, 1003,  // fr: n_K37_to_p_Ar37_reaclib, rr: p_Ar37_to_n_K37_reaclib
        1, 67, -1, 5, 58, -1, 67, 58, 1020, 980,  // fr: n_K37_to_He4_Cl34_reaclib, rr: He4_Cl34_to_n_K37_reaclib
        5, 67, -1, 2, 72, -1, 67, 72, 1021, 1041,  // fr: He4_K37_to_p_Ca40_reaclib, rr: p_Ca40_to_He4_K37_reaclib
        1, 68, -1, 2, 64, -1, 68, 64, 1022, 1008,  // fr: n_K38_to_p_Ar38_reaclib, rr: p_Ar38_to_n_K38_reaclib
        1, 68, -1, 5, 59, -1, 68, 59, 1023, 985,  // fr: n_K38_to_He4_Cl35_reaclib, rr: He4_Cl35_to_n_K38_reaclib
        5, 68, -1, 2, 73, -1, 68, 73, 1024, 1045,  // fr: He4_K38_to_p_Ca41_reaclib, rr: p_Ca41_to_He4_K38_reaclib
        1, 69, -1, 2, 65, -1, 69, 65, 1025, 1013,  // fr: n_K39_to_p_Ar39_reaclib, rr: p_Ar39_to_n_K39_reaclib
        1, 69, -1, 5, 60, -1, 69, 60, 1026, 991,  // fr: n_K39_to_He4_Cl36_reaclib, rr: He4_Cl36_to_n_K39_reaclib
        2, 69, -1, 5, 62, -1, 69, 62, 1027, 1000,  // fr: p_K39_to_He4_Ar36_reaclib, rr: He4_Ar36_to_p_K39_reaclib
        1, 70, -1, 2, 66, -1, 70, 66, 1029, 1016,  // fr: n_K40_to_p_Ar40_reaclib, rr: p_Ar40_to_n_K40_reaclib
        1, 70, -1, 5, 61, -1, 70, 61, 1030, 995,  // fr: n_K40_to_He4_Cl37_reaclib, rr: He4_Cl37_to_n_K40_reaclib
        2, 70, -1, 1, 72, -1, 70, 72, 1031, 1039,  // fr: p_K40_to_n_Ca40_reaclib, rr: n_Ca40_to_p_K40_reaclib
        2, 70, -1, 5, 63, -1, 70, 63, 1032, 1006,  // fr: p_K40_to_He4_Ar37_reaclib, rr: He4_Ar37_to_p_K40_reaclib
        5, 70, -1, 2, 75, -1, 70, 75, 1034, 1054,  // fr: He4_K40_to_p_Ca43_reaclib, rr: p_Ca43_to_He4_K40_reaclib
        2, 71, -1, 5, 64, -1, 71, 64, 1036, 1011,  // fr: p_K41_to_He4_Ar38_reaclib, rr: He4_Ar38_to_p_K41_reaclib
        5, 71, -1, 2, 76, -1, 71, 76, 1038, 1058,  // fr: He4_K41_to_p_Ca44_reaclib, rr: p_Ca44_to_He4_K41_reaclib
        1, 72, -1, 5, 63, -1, 72, 63, 1040, 1005,  // fr: n_Ca40_to_He4_Ar37_reaclib, rr: He4_Ar37_to_n_Ca40_reaclib
        1, 73, -1, 2, 71, -1, 73, 71, 1043, 1035,  // fr: n_Ca41_to_p_K41_reaclib, rr: p_K41_to_n_Ca41_reaclib
        1, 73, -1, 5, 64, -1, 73, 64, 1044, 1010,  // fr: n_Ca41_to_He4_Ar38_reaclib, rr: He4_Ar38_to_n_Ca41_reaclib
        1, 74, -1, 5, 65, -1, 74, 65, 1048, 1015,  // fr: n_Ca42_to_He4_Ar39_reaclib, rr: He4_Ar39_to_n_Ca42_reaclib
        2, 74, -1, 5, 69, -1, 74, 69, 1049, 1028,  // fr: p_Ca42_to_He4_K39_reaclib, rr: He4_K39_to_p_Ca42_reaclib
        1, 75, -1, 5, 66, -1, 75, 66, 1052, 1018,  // fr: n_Ca43_to_He4_Ar40_reaclib, rr: He4_Ar40_to_n_Ca43_reaclib
        5, 75, -1, 1, 90, -1, 75, 90, 1055, 1113,  // fr: He4_Ca43_to_n_Ti46_reaclib, rr: n_Ti46_to_He4_Ca43_reaclib
        5, 77, -1, 1, 92, -1, 77, 92, 1062, 1125,  // fr: He4_Ca45_to_n_Ti48_reaclib, rr: n_Ti48_to_He4_Ca45_reaclib
        2, 79, -1, 1, 85, -1, 79, 85, 1067, 1092,  // fr: p_Ca47_to_n_Sc47_reaclib, rr: n_Sc47_to_p_Ca47_reaclib
        5, 79, -1, 1, 94, -1, 79, 94, 1068, 1136,  // fr: He4_Ca47_to_n_Ti50_reaclib, rr: n_Ti50_to_He4_Ca47_reaclib
        1, 81, -1, 2, 75, -1, 81, 75, 1071, 1053,  // fr: n_Sc43_to_p_Ca43_reaclib, rr: p_Ca43_to_n_Sc43_reaclib
        1, 81, -1, 5, 70, -1, 81, 70, 1072, 1033,  // fr: n_Sc43_to_He4_K40_reaclib, rr: He4_K40_to_n_Sc43_reaclib
        2, 81, -1, 5, 72, -1, 81, 72, 1073, 1042,  // fr: p_Sc43_to_He4_Ca40_reaclib, rr: He4_Ca40_to_p_Sc43_reaclib
        5, 81, -1, 2, 90, -1, 81, 90, 1075, 1115,  // fr: He4_Sc43_to_p_Ti46_reaclib, rr: p_Ti46_to_He4_Sc43_reaclib
        1, 82, -1, 2, 76, -1, 82, 76, 1076, 1057,  // fr: n_Sc44_to_p_Ca44_reaclib, rr: p_Ca44_to_n_Sc44_reaclib
        1, 82, -1, 5, 71, -1, 82, 71, 1077, 1037,  // fr: n_Sc44_to_He4_K41_reaclib, rr: He4_K41_to_n_Sc44_reaclib
        2, 82, -1, 5, 73, -1, 82, 73, 1079, 1047,  // fr: p_Sc44_to_He4_Ca41_reaclib, rr: He4_Ca41_to_p_Sc44_reaclib
        5, 82, -1, 2, 91, -1, 82, 91, 1081, 1121,  // fr: He4_Sc44_to_p_Ti47_reaclib, rr: p_Ti47_to_He4_Sc44_reaclib
        1, 83, -1, 2, 77, -1, 83, 77, 1082, 1061,  // fr: n_Sc45_to_p_Ca45_reaclib, rr: p_Ca45_to_n_Sc45_reaclib
        2, 83, -1, 5, 74, -1, 83, 74, 1084, 1051,  // fr: p_Sc45_to_He4_Ca42_reaclib, rr: He4_Ca42_to_p_Sc45_reaclib
        5, 83, -1, 2, 92, -1, 83, 92, 1086, 1127,  // fr: He4_Sc45_to_p_Ti48_reaclib, rr: p_Ti48_to_He4_Sc45_reaclib
        1, 84, -1, 2, 78, -1, 84, 78, 1087, 1064,  // fr: n_Sc46_to_p_Ca46_reaclib, rr: p_Ca46_to_n_Sc46_reaclib
        2, 84, -1, 1, 90, -1, 84, 90, 1088, 1112,  // fr: p_Sc46_to_n_Ti46_reaclib, rr: n_Ti46_to_p_Sc46_reaclib
        2, 84, -1, 5, 75, -1, 84, 75, 1089, 1056,  // fr: p_Sc46_to_He4_Ca43_reaclib, rr: He4_Ca43_to_p_Sc46_reaclib
        5, 84, -1, 1, 99, -1, 84, 99, 1090, 1159,  // fr: He4_Sc46_to_n_V49_reaclib, rr: n_V49_to_He4_Sc46_reaclib
        5, 84, -1, 2, 93, -1, 84, 93, 1091, 1133,  // fr: He4_Sc46_to_p_Ti49_reaclib, rr: p_Ti49_to_He4_Sc46_reaclib
        2, 85, -1, 5, 76, -1, 85, 76, 1094, 1060,  // fr: p_Sc47_to_He4_Ca44_reaclib, rr: He4_Ca44_to_p_Sc47_reaclib
        5, 85, -1, 2, 94, -1, 85, 94, 1096, 1138,  // fr: He4_Sc47_to_p_Ti50_reaclib, rr: p_Ti50_to_He4_Sc47_reaclib
        1, 86, -1, 2, 80, -1, 86, 80, 1097, 1069,  // fr: n_Sc48_to_p_Ca48_reaclib, rr: p_Ca48_to_n_Sc48_reaclib
        2, 86, -1, 1, 92, -1, 86, 92, 1098, 1124,  // fr: p_Sc48_to_n_Ti48_reaclib, rr: n_Ti48_to_p_Sc48_reaclib
        2, 86, -1, 5, 77, -1, 86, 77, 1099, 1063,  // fr: p_Sc48_to_He4_Ca45_reaclib, rr: He4_Ca45_to_p_Sc48_reaclib
        5, 86, -1, 1, 101, -1, 86, 101, 1100, 1171,  // fr: He4_Sc48_to_n_V51_reaclib, rr: n_V51_to_He4_Sc48_reaclib
        5, 86, -1, 2, 95, -1, 86, 95, 1101, 1142,  // fr: He4_Sc48_to_p_Ti51_reaclib, rr: p_Ti51_to_He4_Sc48_reaclib
        2, 87, -1, 1, 93, -1, 87, 93, 1102, 1130,  // fr: p_Sc49_to_n_Ti49_reaclib, rr: n_Ti49_to_p_Sc49_reaclib
        2, 87, -1, 5, 78, -1, 87, 78, 1103, 1066,  // fr: p_Sc49_to_He4_Ca46_reaclib, rr: He4_Ca46_to_p_Sc49_reaclib
        1, 88, -1, 2, 82, -1, 88, 82, 1105, 1078,  // fr: n_Ti44_to_p_Sc44_reaclib, rr: p_Sc44_to_n_Ti44_reaclib
        1, 88, -1, 5, 73, -1, 88, 73, 1106, 1046,  // fr: n_Ti44_to_He4_Ca41_reaclib, rr: He4_Ca41_to_n_Ti44_reaclib
        1, 89, -1, 2, 83, -1, 89, 83, 1108, 1083,  // fr: n_Ti45_to_p_Sc45_reaclib, rr: p_Sc45_to_n_Ti45_reaclib
        1, 89, -1, 5, 74, -1, 89, 74, 1109, 1050,  // fr: n_Ti45_to_He4_Ca42_reaclib, rr: He4_Ca42_to_n_Ti45_reaclib
        5, 89, -1, 2, 98, -1, 89, 98, 1111, 1155,  // fr: He4_Ti45_to_p_V48_reaclib, rr: p_V48_to_He4_Ti45_reaclib
        1, 91, -1, 2, 85, -1, 91, 85, 1118, 1093,  // fr: n_Ti47_to_p_Sc47_reaclib, rr: p_Sc47_to_n_Ti47_reaclib
        1, 91, -1, 5, 76, -1, 91, 76, 1119, 1059,  // fr: n_Ti47_to_He4_Ca44_reaclib, rr: He4_Ca44_to_n_Ti47_reaclib
        1, 93, -1, 5, 78, -1, 93, 78, 1131, 1065,  // fr: n_Ti49_to_He4_Ca46_reaclib, rr: He4_Ca46_to_n_Ti49_reaclib
        5, 93, -1, 1, 107, -1, 93, 107, 1134, 1201,  // fr: He4_Ti49_to_n_Cr52_reaclib, rr: n_Cr52_to_He4_Ti49_reaclib
        1, 95, -1, 5, 80, -1, 95, 80, 1140, 1070,  // fr: n_Ti51_to_He4_Ca48_reaclib, rr: He4_Ca48_to_n_Ti51_reaclib
        2, 95, -1, 1, 101, -1, 95, 101, 1141, 1170,  // fr: p_Ti51_to_n_V51_reaclib, rr: n_V51_to_p_Ti51_reaclib
        5, 95, -1, 1, 109, -1, 95, 109, 1143, 1210,  // fr: He4_Ti51_to_n_Cr54_reaclib, rr: n_Cr54_to_He4_Ti51_reaclib
        1, 96, -1, 2, 90, -1, 96, 90, 1144, 1114,  // fr: n_V46_to_p_Ti46_reaclib, rr: p_Ti46_to_n_V46_reaclib
        1, 96, -1, 5, 81, -1, 96, 81, 1145, 1074,  // fr: n_V46_to_He4_Sc43_reaclib, rr: He4_Sc43_to_n_V46_reaclib
        5, 96, -1, 2, 104, -1, 96, 104, 1146, 1185,  // fr: He4_V46_to_p_Cr49_reaclib, rr: p_Cr49_to_He4_V46_reaclib
        1, 97, -1, 2, 91, -1, 97, 91, 1147, 1120,  // fr: n_V47_to_p_Ti47_reaclib, rr: p_Ti47_to_n_V47_reaclib
        1, 97, -1, 5, 82, -1, 97, 82, 1148, 1080,  // fr: n_V47_to_He4_Sc44_reaclib, rr: He4_Sc44_to_n_V47_reaclib
        2, 97, -1, 5, 88, -1, 97, 88, 1149, 1107,  // fr: p_V47_to_He4_Ti44_reaclib, rr: He4_Ti44_to_p_V47_reaclib
        5, 97, -1, 2, 105, -1, 97, 105, 1151, 1191,  // fr: He4_V47_to_p_Cr50_reaclib, rr: p_Cr50_to_He4_V47_reaclib
        1, 98, -1, 2, 92, -1, 98, 92, 1152, 1126,  // fr: n_V48_to_p_Ti48_reaclib, rr: p_Ti48_to_n_V48_reaclib
        1, 98, -1, 5, 83, -1, 98, 83, 1153, 1085,  // fr: n_V48_to_He4_Sc45_reaclib, rr: He4_Sc45_to_n_V48_reaclib
        5, 98, -1, 2, 106, -1, 98, 106, 1157, 1197,  // fr: He4_V48_to_p_Cr51_reaclib, rr: p_Cr51_to_He4_V48_reaclib
        1, 99, -1, 2, 93, -1, 99, 93, 1158, 1132,  // fr: n_V49_to_p_Ti49_reaclib, rr: p_Ti49_to_n_V49_reaclib
        2, 99, -1, 5, 90, -1, 99, 90, 1161, 1117,  // fr: p_V49_to_He4_Ti46_reaclib, rr: He4_Ti46_to_p_V49_reaclib
        5, 99, -1, 2, 107, -1, 99, 107, 1163, 1203,  // fr: He4_V49_to_p_Cr52_reaclib, rr: p_Cr52_to_He4_V49_reaclib
        1, 100, -1, 2, 94, -1, 100, 94, 1164, 1137,  // fr: n_V50_to_p_Ti50_reaclib, rr: p_Ti50_to_n_V50_reaclib
        1, 100, -1, 5, 85, -1, 100, 85, 1165, 1095,  // fr: n_V50_to_He4_Sc47_reaclib, rr: He4_Sc47_to_n_V50_reaclib
        2, 100, -1, 1, 105, -1, 100, 105, 1166, 1188,  // fr: p_V50_to_n_Cr50_reaclib, rr: n_Cr50_to_p_V50_reaclib
        2, 100, -1, 5, 91, -1, 100, 91, 1167, 1123,  // fr: p_V50_to_He4_Ti47_reaclib, rr: He4_Ti47_to_p_V50_reaclib
        5, 100, -1, 2, 108, -1, 100, 108, 1169, 1208,  // fr: He4_V50_to_p_Cr53_reaclib, rr: p_Cr53_to_He4_V50_reaclib
        2, 101, -1, 5, 92, -1, 101, 92, 1173, 1129,  // fr: p_V51_to_He4_Ti48_reaclib, rr: He4_Ti48_to_p_V51_reaclib
        1, 102, -1, 5, 87, -1, 102, 87, 1176, 1104,  // fr: n_V52_to_He4_Sc49_reaclib, rr: He4_Sc49_to_n_V52_reaclib
        2, 102, -1, 1, 107, -1, 102, 107, 1177, 1200,  // fr: p_V52_to_n_Cr52_reaclib, rr: n_Cr52_to_p_V52_reaclib
        2, 102, -1, 5, 93, -1, 102, 93, 1178, 1135,  // fr: p_V52_to_He4_Ti49_reaclib, rr: He4_Ti49_to_p_V52_reaclib
        5, 102, -1, 1, 115, -1, 102, 115, 1179, 1241,  // fr: He4_V52_to_n_Mn55_reaclib, rr: n_Mn55_to_He4_V52_reaclib
        1, 103, -1, 2, 98, -1, 103, 98, 1180, 1154,  // fr: n_Cr48_to_p_V48_reaclib, rr: p_V48_to_n_Cr48_reaclib
        1, 103, -1, 5, 89, -1, 103, 89, 1181, 1110,  // fr: n_Cr48_to_He4_Ti45_reaclib, rr: He4_Ti45_to_n_Cr48_reaclib
        5, 103, -1, 2, 111, -1, 103, 111, 1182, 1220,  // fr: He4_Cr48_to_p_Mn51_reaclib, rr: p_Mn51_to_He4_Cr48_reaclib
        1, 104, -1, 2, 99, -1, 104, 99, 1183, 1160,  // fr: n_Cr49_to_p_V49_reaclib, rr: p_V49_to_n_Cr49_reaclib
        1, 104, -1, 5, 90, -1, 104, 90, 1184, 1116,  // fr: n_Cr49_to_He4_Ti46_reaclib, rr: He4_Ti46_to_n_Cr49_reaclib
        5, 104, -1, 2, 112, -1, 104, 112, 1187, 1226,  // fr: He4_Cr49_to_p_Mn52_reaclib, rr: p_Mn52_to_He4_Cr49_reaclib
        1, 105, -1, 5, 91, -1, 105, 91, 1189, 1122,  // fr: n_Cr50_to_He4_Ti47_reaclib, rr: He4_Ti47_to_n_Cr50_reaclib
        1, 106, -1, 2, 101, -1, 106, 101, 1194, 1172,  // fr: n_Cr51_to_p_V51_reaclib, rr: p_V51_to_n_Cr51_reaclib
        1, 106, -1, 5, 92, -1, 106, 92, 1195, 1128,  // fr: n_Cr51_to_He4_Ti48_reaclib, rr: He4_Ti48_to_n_Cr51_reaclib
        1, 108, -1, 5, 94, -1, 108, 94, 1206, 1139,  // fr: n_Cr53_to_He4_Ti50_reaclib, rr: He4_Ti50_to_n_Cr53_reaclib
        2, 109, -1, 5, 101, -1, 109, 101, 1212, 1175,  // fr: p_Cr54_to_He4_V51_reaclib, rr: He4_V51_to_p_Cr54_reaclib
        1, 110, -1, 2, 105, -1, 110, 105, 1214, 1190,  // fr: n_Mn50_to_p_Cr50_reaclib, rr: p_Cr50_to_n_Mn50_reaclib
        1, 110, -1, 5, 97, -1, 110, 97, 1215, 1150,  // fr: n_Mn50_to_He4_V47_reaclib, rr: He4_V47_to_n_Mn50_reaclib
        5, 110, -1, 2, 117, -1, 110, 117, 1217, 1252,  // fr: He4_Mn50_to_p_Fe53_reaclib, rr: p_Fe53_to_He4_Mn50_reaclib
        1, 111, -1, 2, 106, -1, 111, 106, 1218, 1196,  // fr: n_Mn51_to_p_Cr51_reaclib, rr: p_Cr51_to_n_Mn51_reaclib
        1, 111, -1, 5, 98, -1, 111, 98, 1219, 1156,  // fr: n_Mn51_to_He4_V48_reaclib, rr: He4_V48_to_n_Mn51_reaclib
        5, 111, -1, 2, 118, -1, 111, 118, 1222, 1258,  // fr: He4_Mn51_to_p_Fe54_reaclib, rr: p_Fe54_to_He4_Mn51_reaclib
        1, 112, -1, 2, 107, -1, 112, 107, 1223, 1202,  // fr: n_Mn52_to_p_Cr52_reaclib, rr: p_Cr52_to_n_Mn52_reaclib
        1, 112, -1, 5, 99, -1, 112, 99, 1224, 1162,  // fr: n_Mn52_to_He4_V49_reaclib, rr: He4_V49_to_n_Mn52_reaclib
        5, 112, -1, 2, 119, -1, 112, 119, 1228, 1264,  // fr: He4_Mn52_to_p_Fe55_reaclib, rr: p_Fe55_to_He4_Mn52_reaclib
        1, 113, -1, 2, 108, -1, 113, 108, 1229, 1207,  // fr: n_Mn53_to_p_Cr53_reaclib, rr: p_Cr53_to_n_Mn53_reaclib
        1, 113, -1, 5, 100, -1, 113, 100, 1230, 1168,  // fr: n_Mn53_to_He4_V50_reaclib, rr: He4_V50_to_n_Mn53_reaclib
        2, 113, -1, 5, 105, -1, 113, 105, 1232, 1193,  // fr: p_Mn53_to_He4_Cr50_reaclib, rr: He4_Cr50_to_p_Mn53_reaclib
        5, 113, -1, 2, 120, -1, 113, 120, 1234, 1269,  // fr: He4_Mn53_to_p_Fe56_reaclib, rr: p_Fe56_to_He4_Mn53_reaclib
        1, 114, -1, 2, 109, -1, 114, 109, 1235, 1211,  // fr: n_Mn54_to_p_Cr54_reaclib, rr: p_Cr54_to_n_Mn54_reaclib
        1, 114, -1, 5, 101, -1, 114, 101, 1236, 1174,  // fr: n_Mn54_to_He4_V51_reaclib, rr: He4_V51_to_n_Mn54_reaclib
        2, 114, -1, 5, 106, -1, 114, 106, 1238, 1199,  // fr: p_Mn54_to_He4_Cr51_reaclib, rr: He4_Cr51_to_p_Mn54_reaclib
        2, 115, -1, 5, 107, -1, 115, 107, 1243, 1205,  // fr: p_Mn55_to_He4_Cr52_reaclib, rr: He4_Cr52_to_p_Mn55_reaclib
        1, 116, -1, 2, 112, -1, 116, 112, 1246, 1225,  // fr: n_Fe52_to_p_Mn52_reaclib, rr: p_Mn52_to_n_Fe52_reaclib
        1, 116, -1, 5, 104, -1, 116, 104, 1247, 1186,  // fr: n_Fe52_to_He4_Cr49_reaclib, rr: He4_Cr49_to_n_Fe52_reaclib
        5, 116, -1, 2, 125, -1, 116, 125, 1248, 1288,  // fr: He4_Fe52_to_p_Co55_reaclib, rr: p_Co55_to_He4_Fe52_reaclib
        1, 117, -1, 2, 113, -1, 117, 113, 1249, 1231,  // fr: n_Fe53_to_p_Mn53_reaclib, rr: p_Mn53_to_n_Fe53_reaclib
        1, 117, -1, 5, 105, -1, 117, 105, 1250, 1192,  // fr: n_Fe53_to_He4_Cr50_reaclib, rr: He4_Cr50_to_n_Fe53_reaclib
        5, 117, -1, 2, 126, -1, 117, 126, 1254, 1294,  // fr: He4_Fe53_to_p_Co56_reaclib, rr: p_Co56_to_He4_Fe53_reaclib
        1, 118, -1, 2, 114, -1, 118, 114, 1255, 1237,  // fr: n_Fe54_to_p_Mn54_reaclib, rr: p_Mn54_to_n_Fe54_reaclib
        1, 118, -1, 5, 106, -1, 118, 106, 1256, 1198,  // fr: n_Fe54_to_He4_Cr51_reaclib, rr: He4_Cr51_to_n_Fe54_reaclib
        1, 119, -1, 2, 115, -1, 119, 115, 1261, 1242,  // fr: n_Fe55_to_p_Mn55_reaclib, rr: p_Mn55_to_n_Fe55_reaclib
        1, 119, -1, 5, 107, -1, 119, 107, 1262, 1204,  // fr: n_Fe55_to_He4_Cr52_reaclib, rr: He4_Cr52_to_n_Fe55_reaclib
        1, 120, -1, 5, 108, -1, 120, 108, 1267, 1209,  // fr: n_Fe56_to_He4_Cr53_reaclib, rr: He4_Cr53_to_n_Fe56_reaclib
        1, 121, -1, 5, 109, -1, 121, 109, 1272, 1213,  // fr: n_Fe57_to_He4_Cr54_reaclib, rr: He4_Cr54_to_n_Fe57_reaclib
        2, 121, -1, 5, 114, -1, 121, 114, 1274, 1240,  // fr: p_Fe57_to_He4_Mn54_reaclib, rr: He4_Mn54_to_p_Fe57_reaclib
        2, 122, -1, 5, 115, -1, 122, 115, 1277, 1245,  // fr: p_Fe58_to_He4_Mn55_reaclib, rr: He4_Mn55_to_p_Fe58_reaclib
        1, 123, -1, 2, 117, -1, 123, 117, 1279, 1251,  // fr: n_Co53_to_p_Fe53_reaclib, rr: p_Fe53_to_n_Co53_reaclib
        1, 123, -1, 5, 110, -1, 123, 110, 1280, 1216,  // fr: n_Co53_to_He4_Mn50_reaclib, rr: He4_Mn50_to_n_Co53_reaclib
        5, 123, -1, 2, 130, -1, 123, 130, 1281, 1315,  // fr: He4_Co53_to_p_Ni56_reaclib, rr: p_Ni56_to_He4_Co53_reaclib
        1, 124, -1, 2, 118, -1, 124, 118, 1282, 1257,  // fr: n_Co54_to_p_Fe54_reaclib, rr: p_Fe54_to_n_Co54_reaclib
        1, 124, -1, 5, 111, -1, 124, 111, 1283, 1221,  // fr: n_Co54_to_He4_Mn51_reaclib, rr: He4_Mn51_to_n_Co54_reaclib
        5, 124, -1, 2, 131, -1, 124, 131, 1285, 1321,  // fr: He4_Co54_to_p_Ni57_reaclib, rr: p_Ni57_to_He4_Co54_reaclib
        1, 125, -1, 2, 119, -1, 125, 119, 1286, 1263,  // fr: n_Co55_to_p_Fe55_reaclib, rr: p_Fe55_to_n_Co55_reaclib
        1, 125, -1, 5, 112, -1, 125, 112, 1287, 1227,  // fr: n_Co55_to_He4_Mn52_reaclib, rr: He4_Mn52_to_n_Co55_reaclib
        5, 125, -1, 2, 132, -1, 125, 132, 1290, 1327,  // fr: He4_Co55_to_p_Ni58_reaclib, rr: p_Ni58_to_He4_Co55_reaclib
        1, 126, -1, 2, 120, -1, 126, 120, 1291, 1268,  // fr: n_Co56_to_p_Fe56_reaclib, rr: p_Fe56_to_n_Co56_reaclib
        1, 126, -1, 5, 113, -1, 126, 113, 1292, 1233,  // fr: n_Co56_to_He4_Mn53_reaclib, rr: He4_Mn53_to_n_Co56_reaclib
        5, 126, -1, 2, 133, -1, 126, 133, 1296, 1333,  // fr: He4_Co56_to_p_Ni59_reaclib, rr: p_Ni59_to_He4_Co56_reaclib
        1, 127, -1, 2, 121, -1, 127, 121, 1297, 1273,  // fr: n_Co57_to_p_Fe57_reaclib, rr: p_Fe57_to_n_Co57_reaclib
        1, 127, -1, 5, 114, -1, 127, 114, 1298, 1239,  // fr: n_Co57_to_He4_Mn54_reaclib, rr: He4_Mn54_to_n_Co57_reaclib
        2, 127, -1, 5, 118, -1, 127, 118, 1300, 1260,  // fr: p_Co57_to_He4_Fe54_reaclib, rr: He4_Fe54_to_p_Co57_reaclib
        5, 127, -1, 2, 134, -1, 127, 134, 1302, 1338,  // fr: He4_Co57_to_p_Ni60_reaclib, rr: p_Ni60_to_He4_Co57_reaclib
        1, 128, -1, 2, 122, -1, 128, 122, 1303, 1276,  // fr: n_Co58_to_p_Fe58_reaclib, rr: p_Fe58_to_n_Co58_reaclib
        1, 128, -1, 5, 115, -1, 128, 115, 1304, 1244,  // fr: n_Co58_to_He4_Mn55_reaclib, rr: He4_Mn55_to_n_Co58_reaclib
        2, 128, -1, 5, 119, -1, 128, 119, 1306, 1266,  // fr: p_Co58_to_He4_Fe55_reaclib, rr: He4_Fe55_to_p_Co58_reaclib
        2, 129, -1, 5, 120, -1, 129, 120, 1310, 1271,  // fr: p_Co59_to_He4_Fe56_reaclib, rr: He4_Fe56_to_p_Co59_reaclib
        1, 130, -1, 2, 126, -1, 130, 126, 1313, 1293,  // fr: n_Ni56_to_p_Co56_reaclib, rr: p_Co56_to_n_Ni56_reaclib
        1, 130, -1, 5, 117, -1, 130, 117, 1314, 1253,  // fr: n_Ni56_to_He4_Fe53_reaclib, rr: He4_Fe53_to_n_Ni56_reaclib
        1, 131, -1, 2, 127, -1, 131, 127, 1318, 1299,  // fr: n_Ni57_to_p_Co57_reaclib, rr: p_Co57_to_n_Ni57_reaclib
        1, 131, -1, 5, 118, -1, 131, 118, 1319, 1259,  // fr: n_Ni57_to_He4_Fe54_reaclib, rr: He4_Fe54_to_n_Ni57_reaclib
        1, 132, -1, 2, 128, -1, 132, 128, 1324, 1305,  // fr: n_Ni58_to_p_Co58_reaclib, rr: p_Co58_to_n_Ni58_reaclib
        1, 132, -1, 5, 119, -1, 132, 119, 1325, 1265,  // fr: n_Ni58_to_He4_Fe55_reaclib, rr: He4_Fe55_to_n_Ni58_reaclib
        1, 133, -1, 2, 129, -1, 133, 129, 1330, 1309,  // fr: n_Ni59_to_p_Co59_reaclib, rr: p_Co59_to_n_Ni59_reaclib
        1, 133, -1, 5, 120, -1, 133, 120, 1331, 1270,  // fr: n_Ni59_to_He4_Fe56_reaclib, rr: He4_Fe56_to_n_Ni59_reaclib
        1, 134, -1, 5, 121, -1, 134, 121, 1336, 1275,  // fr: n_Ni60_to_He4_Fe57_reaclib, rr: He4_Fe57_to_n_Ni60_reaclib
        1, 135, -1, 5, 122, -1, 135, 122, 1341, 1278,  // fr: n_Ni61_to_He4_Fe58_reaclib, rr: He4_Fe58_to_n_Ni61_reaclib
        2, 135, -1, 5, 128, -1, 135, 128, 1343, 1308,  // fr: p_Ni61_to_He4_Co58_reaclib, rr: He4_Co58_to_p_Ni61_reaclib
        2, 136, -1, 5, 129, -1, 136, 129, 1347, 1312,  // fr: p_Ni62_to_He4_Co59_reaclib, rr: He4_Co59_to_p_Ni62_reaclib
        1, 139, -1, 2, 131, -1, 139, 131, 1353, 1320,  // fr: n_Cu57_to_p_Ni57_reaclib, rr: p_Ni57_to_n_Cu57_reaclib
        1, 139, -1, 5, 124, -1, 139, 124, 1354, 1284,  // fr: n_Cu57_to_He4_Co54_reaclib, rr: He4_Co54_to_n_Cu57_reaclib
        5, 139, -1, 2, 149, -1, 139, 149, 1355, 1396,  // fr: He4_Cu57_to_p_Zn60_reaclib, rr: p_Zn60_to_He4_Cu57_reaclib
        1, 140, -1, 2, 132, -1, 140, 132, 1356, 1326,  // fr: n_Cu58_to_p_Ni58_reaclib, rr: p_Ni58_to_n_Cu58_reaclib
        1, 140, -1, 5, 125, -1, 140, 125, 1357, 1289,  // fr: n_Cu58_to_He4_Co55_reaclib, rr: He4_Co55_to_n_Cu58_reaclib
        1, 141, -1, 2, 133, -1, 141, 133, 1359, 1332,  // fr: n_Cu59_to_p_Ni59_reaclib, rr: p_Ni59_to_n_Cu59_reaclib
        1, 141, -1, 5, 126, -1, 141, 126, 1360, 1295,  // fr: n_Cu59_to_He4_Co56_reaclib, rr: He4_Co56_to_n_Cu59_reaclib
        2, 141, -1, 5, 130, -1, 141, 130, 1362, 1317,  // fr: p_Cu59_to_He4_Ni56_reaclib, rr: He4_Ni56_to_p_Cu59_reaclib
        1, 142, -1, 2, 134, -1, 142, 134, 1365, 1337,  // fr: n_Cu60_to_p_Ni60_reaclib, rr: p_Ni60_to_n_Cu60_reaclib
        1, 142, -1, 5, 127, -1, 142, 127, 1366, 1301,  // fr: n_Cu60_to_He4_Co57_reaclib, rr: He4_Co57_to_n_Cu60_reaclib
        2, 142, -1, 5, 131, -1, 142, 131, 1368, 1323,  // fr: p_Cu60_to_He4_Ni57_reaclib, rr: He4_Ni57_to_p_Cu60_reaclib
        1, 143, -1, 2, 135, -1, 143, 135, 1371, 1342,  // fr: n_Cu61_to_p_Ni61_reaclib, rr: p_Ni61_to_n_Cu61_reaclib
        1, 143, -1, 5, 128, -1, 143, 128, 1372, 1307,  // fr: n_Cu61_to_He4_Co58_reaclib, rr: He4_Co58_to_n_Cu61_reaclib
        2, 143, -1, 5, 132, -1, 143, 132, 1374, 1329,  // fr: p_Cu61_to_He4_Ni58_reaclib, rr: He4_Ni58_to_p_Cu61_reaclib
        1, 144, -1, 2, 136, -1, 144, 136, 1377, 1346,  // fr: n_Cu62_to_p_Ni62_reaclib, rr: p_Ni62_to_n_Cu62_reaclib
        1, 144, -1, 5, 129, -1, 144, 129, 1378, 1311,  // fr: n_Cu62_to_He4_Co59_reaclib, rr: He4_Co59_to_n_Cu62_reaclib
        2, 144, -1, 5, 133, -1, 144, 133, 1380, 1335,  // fr: p_Cu62_to_He4_Ni59_reaclib, rr: He4_Ni59_to_p_Cu62_reaclib
        1, 145, -1, 2, 137, -1, 145, 137, 1382, 1350,  // fr: n_Cu63_to_p_Ni63_reaclib, rr: p_Ni63_to_n_Cu63_reaclib
        2, 145, -1, 5, 134, -1, 145, 134, 1384, 1340,  // fr: p_Cu63_to_He4_Ni60_reaclib, rr: He4_Ni60_to_p_Cu63_reaclib
        1, 146, -1, 2, 138, -1, 146, 138, 1386, 1352,  // fr: n_Cu64_to_p_Ni64_reaclib, rr: p_Ni64_to_n_Cu64_reaclib
        2, 146, -1, 5, 135, -1, 146, 135, 1388, 1345,  // fr: p_Cu64_to_He4_Ni61_reaclib, rr: He4_Ni61_to_p_Cu64_reaclib
        2, 147, -1, 5, 136, -1, 147, 136, 1390, 1349,  // fr: p_Cu65_to_He4_Ni62_reaclib, rr: He4_Ni62_to_p_Cu65_reaclib
        1, 148, -1, 2, 141, -1, 148, 141, 1391, 1361,  // fr: n_Zn59_to_p_Cu59_reaclib, rr: p_Cu59_to_n_Zn59_reaclib
        1, 148, -1, 5, 130, -1, 148, 130, 1392, 1316,  // fr: n_Zn59_to_He4_Ni56_reaclib, rr: He4_Ni56_to_n_Zn59_reaclib
        1, 149, -1, 2, 142, -1, 149, 142, 1394, 1367,  // fr: n_Zn60_to_p_Cu60_reaclib, rr: p_Cu60_to_n_Zn60_reaclib
        1, 149, -1, 5, 131, -1, 149, 131, 1395, 1322,  // fr: n_Zn60_to_He4_Ni57_reaclib, rr: He4_Ni57_to_n_Zn60_reaclib
        1, 150, -1, 2, 143, -1, 150, 143, 1399, 1373,  // fr: n_Zn61_to_p_Cu61_reaclib, rr: p_Cu61_to_n_Zn61_reaclib
        1, 150, -1, 5, 132, -1, 150, 132, 1400, 1328,  // fr: n_Zn61_to_He4_Ni58_reaclib, rr: He4_Ni58_to_n_Zn61_reaclib
        2, 150, -1, 5, 140, -1, 150, 140, 1401, 1358,  // fr: p_Zn61_to_He4_Cu58_reaclib, rr: He4_Cu58_to_p_Zn61_reaclib
        1, 151, -1, 2, 144, -1, 151, 144, 1404, 1379,  // fr: n_Zn62_to_p_Cu62_reaclib, rr: p_Cu62_to_n_Zn62_reaclib
        1, 151, -1, 5, 133, -1, 151, 133, 1405, 1334,  // fr: n_Zn62_to_He4_Ni59_reaclib, rr: He4_Ni59_to_n_Zn62_reaclib
        2, 151, -1, 5, 141, -1, 151, 141, 1407, 1364,  // fr: p_Zn62_to_He4_Cu59_reaclib, rr: He4_Cu59_to_p_Zn62_reaclib
        1, 152, -1, 2, 145, -1, 152, 145, 1408, 1383,  // fr: n_Zn63_to_p_Cu63_reaclib, rr: p_Cu63_to_n_Zn63_reaclib
        1, 152, -1, 5, 134, -1, 152, 134, 1409, 1339,  // fr: n_Zn63_to_He4_Ni60_reaclib, rr: He4_Ni60_to_n_Zn63_reaclib
        2, 152, -1, 5, 142, -1, 152, 142, 1411, 1370,  // fr: p_Zn63_to_He4_Cu60_reaclib, rr: He4_Cu60_to_p_Zn63_reaclib
        1, 153, -1, 2, 146, -1, 153, 146, 1412, 1387,  // fr: n_Zn64_to_p_Cu64_reaclib, rr: p_Cu64_to_n_Zn64_reaclib
        1, 153, -1, 5, 135, -1, 153, 135, 1413, 1344,  // fr: n_Zn64_to_He4_Ni61_reaclib, rr: He4_Ni61_to_n_Zn64_reaclib
        2, 153, -1, 5, 143, -1, 153, 143, 1415, 1376,  // fr: p_Zn64_to_He4_Cu61_reaclib, rr: He4_Cu61_to_p_Zn64_reaclib
        1, 154, -1, 2, 147, -1, 154, 147, 1416, 1389,  // fr: n_Zn65_to_p_Cu65_reaclib, rr: p_Cu65_to_n_Zn65_reaclib
        1, 154, -1, 5, 136, -1, 154, 136, 1417, 1348,  // fr: n_Zn65_to_He4_Ni62_reaclib, rr: He4_Ni62_to_n_Zn65_reaclib
        2, 154, -1, 5, 144, -1, 154, 144, 1418, 1381,  // fr: p_Zn65_to_He4_Cu62_reaclib, rr: He4_Cu62_to_p_Zn65_reaclib
        1, 155, -1, 5, 137, -1, 155, 137, 1419, 1351,  // fr: n_Zn66_to_He4_Ni63_reaclib, rr: He4_Ni63_to_n_Zn66_reaclib
        2, 155, -1, 5, 145, -1, 155, 145, 1420, 1385,  // fr: p_Zn66_to_He4_Cu63_reaclib, rr: He4_Cu63_to_p_Zn66_reaclib
        1, 156, -1, 2, 151, -1, 156, 151, 1421, 1406,  // fr: n_Ga62_to_p_Zn62_reaclib, rr: p_Zn62_to_n_Ga62_reaclib
        1, 156, -1, 5, 141, -1, 156, 141, 1422, 1363,  // fr: n_Ga62_to_He4_Cu59_reaclib, rr: He4_Cu59_to_n_Ga62_reaclib
        2, 156, -1, 5, 148, -1, 156, 148, 1423, 1393,  // fr: p_Ga62_to_He4_Zn59_reaclib, rr: He4_Zn59_to_p_Ga62_reaclib
        1, 157, -1, 2, 152, -1, 157, 152, 1424, 1410,  // fr: n_Ga63_to_p_Zn63_reaclib, rr: p_Zn63_to_n_Ga63_reaclib
        1, 157, -1, 5, 142, -1, 157, 142, 1425, 1369,  // fr: n_Ga63_to_He4_Cu60_reaclib, rr: He4_Cu60_to_n_Ga63_reaclib
        2, 157, -1, 5, 149, -1, 157, 149, 1427, 1398,  // fr: p_Ga63_to_He4_Zn60_reaclib, rr: He4_Zn60_to_p_Ga63_reaclib
        1, 158, -1, 2, 153, -1, 158, 153, 1428, 1414,  // fr: n_Ga64_to_p_Zn64_reaclib, rr: p_Zn64_to_n_Ga64_reaclib
        1, 158, -1, 5, 143, -1, 158, 143, 1429, 1375,  // fr: n_Ga64_to_He4_Cu61_reaclib, rr: He4_Cu61_to_n_Ga64_reaclib
        2, 158, -1, 5, 150, -1, 158, 150, 1431, 1403,  // fr: p_Ga64_to_He4_Zn61_reaclib, rr: He4_Zn61_to_p_Ga64_reaclib
        1, 159, -1, 2, 157, -1, 159, 157, 1432, 1426,  // fr: n_Ge63_to_p_Ga63_reaclib, rr: p_Ga63_to_n_Ge63_reaclib
        1, 159, -1, 5, 149, -1, 159, 149, 1433, 1397,  // fr: n_Ge63_to_He4_Zn60_reaclib, rr: He4_Zn60_to_n_Ge63_reaclib
        1, 160, -1, 2, 158, -1, 160, 158, 1434, 1430,  // fr: n_Ge64_to_p_Ga64_reaclib, rr: p_Ga64_to_n_Ge64_reaclib
        1, 160, -1, 5, 150, -1, 160, 150, 1435, 1402,  // fr: n_Ge64_to_He4_Zn61_reaclib, rr: He4_Zn61_to_n_Ge64_reaclib
        5, 5, 5, 13, -1, -1, 13, -1, 1448, 387  // fr: He4_He4_He4_to_C12_reaclib, rr: C12_to_He4_He4_He4_reaclib
    };
}
#endif

#endif
