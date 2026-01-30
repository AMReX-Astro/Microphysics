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
        N, H1, -1, H2, -1, -1, H2, -1, k_n_p_to_d_reaclib, k_d_to_n_p_reaclib,
        H1, H2, -1, He3, -1, -1, H2, He3, k_p_d_to_He3_reaclib, k_He3_to_p_d_reaclib,
        H2, H2, -1, He4, -1, -1, H2, H2, k_d_d_to_He4_reaclib, k_He4_to_d_d_reaclib,
        He4, H2, -1, Li6, -1, -1, H2, Li6, k_He4_d_to_Li6_reaclib, k_Li6_to_He4_d_reaclib,
        N, He3, -1, He4, -1, -1, He3, -1, k_n_He3_to_He4_reaclib, k_He4_to_n_He3_reaclib,
        He4, He3, -1, Be7, -1, -1, He3, Be7, k_He4_He3_to_Be7_reaclib, k_Be7_to_He4_He3_reaclib,
        N, Li6, -1, Li7, -1, -1, Li6, Li7, k_n_Li6_to_Li7_reaclib, k_Li7_to_n_Li6_reaclib,
        H1, Li6, -1, Be7, -1, -1, Li6, Be7, k_p_Li6_to_Be7_reaclib, k_Be7_to_p_Li6_reaclib,
        He4, Li6, -1, B10, -1, -1, Li6, B10, k_He4_Li6_to_B10_reaclib, k_B10_to_He4_Li6_reaclib,
        He4, Li7, -1, B11, -1, -1, Li7, B11, k_He4_Li7_to_B11_reaclib, k_B11_to_He4_Li7_reaclib,
        H1, Be7, -1, B8, -1, -1, Be7, B8, k_p_Be7_to_B8_reaclib, k_B8_to_p_Be7_reaclib,
        H1, Be9, -1, B10, -1, -1, Be9, B10, k_p_Be9_to_B10_reaclib, k_B10_to_p_Be9_reaclib,
        N, B10, -1, B11, -1, -1, B10, B11, k_n_B10_to_B11_reaclib, k_B11_to_n_B10_reaclib,
        H1, B11, -1, C12, -1, -1, B11, C12, k_p_B11_to_C12_reaclib, k_C12_to_p_B11_reaclib,
        N, C12, -1, C13, -1, -1, C12, C13, k_n_C12_to_C13_reaclib, k_C13_to_n_C12_reaclib,
        H1, C12, -1, N13, -1, -1, C12, N13, k_p_C12_to_N13_reaclib, k_N13_to_p_C12_reaclib,
        He4, C12, -1, O16, -1, -1, C12, O16, k_He4_C12_to_O16_reaclib, k_O16_to_He4_C12_reaclib,
        N, C13, -1, C14, -1, -1, C13, C14, k_n_C13_to_C14_reaclib, k_C14_to_n_C13_reaclib,
        H1, C13, -1, N14, -1, -1, C13, N14, k_p_C13_to_N14_reaclib, k_N14_to_p_C13_reaclib,
        H1, C14, -1, N15, -1, -1, C14, N15, k_p_C14_to_N15_reaclib, k_N15_to_p_C14_reaclib,
        He4, C14, -1, O18, -1, -1, C14, O18, k_He4_C14_to_O18_reaclib, k_O18_to_He4_C14_reaclib,
        N, N13, -1, N14, -1, -1, N13, N14, k_n_N13_to_N14_reaclib, k_N14_to_n_N13_reaclib,
        H1, N13, -1, O14, -1, -1, N13, O14, k_p_N13_to_O14_reaclib, k_O14_to_p_N13_reaclib,
        N, N14, -1, N15, -1, -1, N14, N15, k_n_N14_to_N15_reaclib, k_N15_to_n_N14_reaclib,
        H1, N14, -1, O15, -1, -1, N14, O15, k_p_N14_to_O15_reaclib, k_O15_to_p_N14_reaclib,
        He4, N14, -1, F18, -1, -1, N14, F18, k_He4_N14_to_F18_reaclib, k_F18_to_He4_N14_reaclib,
        H1, N15, -1, O16, -1, -1, N15, O16, k_p_N15_to_O16_reaclib, k_O16_to_p_N15_reaclib,
        He4, N15, -1, F19, -1, -1, N15, F19, k_He4_N15_to_F19_reaclib, k_F19_to_He4_N15_reaclib,
        N, O14, -1, O15, -1, -1, O14, O15, k_n_O14_to_O15_reaclib, k_O15_to_n_O14_reaclib,
        He4, O14, -1, Ne18, -1, -1, O14, Ne18, k_He4_O14_to_Ne18_reaclib, k_Ne18_to_He4_O14_reaclib,
        N, O15, -1, O16, -1, -1, O15, O16, k_n_O15_to_O16_reaclib, k_O16_to_n_O15_reaclib,
        He4, O15, -1, Ne19, -1, -1, O15, Ne19, k_He4_O15_to_Ne19_reaclib, k_Ne19_to_He4_O15_reaclib,
        N, O16, -1, O17, -1, -1, O16, O17, k_n_O16_to_O17_reaclib, k_O17_to_n_O16_reaclib,
        H1, O16, -1, F17, -1, -1, O16, F17, k_p_O16_to_F17_reaclib, k_F17_to_p_O16_reaclib,
        He4, O16, -1, Ne20, -1, -1, O16, Ne20, k_He4_O16_to_Ne20_reaclib, k_Ne20_to_He4_O16_reaclib,
        N, O17, -1, O18, -1, -1, O17, O18, k_n_O17_to_O18_reaclib, k_O18_to_n_O17_reaclib,
        H1, O17, -1, F18, -1, -1, O17, F18, k_p_O17_to_F18_reaclib, k_F18_to_p_O17_reaclib,
        He4, O17, -1, Ne21, -1, -1, O17, Ne21, k_He4_O17_to_Ne21_reaclib, k_Ne21_to_He4_O17_reaclib,
        H1, O18, -1, F19, -1, -1, O18, F19, k_p_O18_to_F19_reaclib, k_F19_to_p_O18_reaclib,
        He4, O18, -1, Ne22, -1, -1, O18, Ne22, k_He4_O18_to_Ne22_reaclib, k_Ne22_to_He4_O18_reaclib,
        N, F17, -1, F18, -1, -1, F17, F18, k_n_F17_to_F18_reaclib, k_F18_to_n_F17_reaclib,
        H1, F17, -1, Ne18, -1, -1, F17, Ne18, k_p_F17_to_Ne18_reaclib, k_Ne18_to_p_F17_reaclib,
        He4, F17, -1, Na21, -1, -1, F17, Na21, k_He4_F17_to_Na21_reaclib, k_Na21_to_He4_F17_reaclib,
        N, F18, -1, F19, -1, -1, F18, F19, k_n_F18_to_F19_reaclib, k_F19_to_n_F18_reaclib,
        H1, F18, -1, Ne19, -1, -1, F18, Ne19, k_p_F18_to_Ne19_reaclib, k_Ne19_to_p_F18_reaclib,
        He4, F18, -1, Na22, -1, -1, F18, Na22, k_He4_F18_to_Na22_reaclib, k_Na22_to_He4_F18_reaclib,
        H1, F19, -1, Ne20, -1, -1, F19, Ne20, k_p_F19_to_Ne20_reaclib, k_Ne20_to_p_F19_reaclib,
        He4, F19, -1, Na23, -1, -1, F19, Na23, k_He4_F19_to_Na23_reaclib, k_Na23_to_He4_F19_reaclib,
        N, Ne18, -1, Ne19, -1, -1, Ne18, Ne19, k_n_Ne18_to_Ne19_reaclib, k_Ne19_to_n_Ne18_reaclib,
        N, Ne19, -1, Ne20, -1, -1, Ne19, Ne20, k_n_Ne19_to_Ne20_reaclib, k_Ne20_to_n_Ne19_reaclib,
        He4, Ne19, -1, Mg23, -1, -1, Ne19, Mg23, k_He4_Ne19_to_Mg23_reaclib, k_Mg23_to_He4_Ne19_reaclib,
        N, Ne20, -1, Ne21, -1, -1, Ne20, Ne21, k_n_Ne20_to_Ne21_reaclib, k_Ne21_to_n_Ne20_reaclib,
        H1, Ne20, -1, Na21, -1, -1, Ne20, Na21, k_p_Ne20_to_Na21_reaclib, k_Na21_to_p_Ne20_reaclib,
        He4, Ne20, -1, Mg24, -1, -1, Ne20, Mg24, k_He4_Ne20_to_Mg24_reaclib, k_Mg24_to_He4_Ne20_reaclib,
        N, Ne21, -1, Ne22, -1, -1, Ne21, Ne22, k_n_Ne21_to_Ne22_reaclib, k_Ne22_to_n_Ne21_reaclib,
        H1, Ne21, -1, Na22, -1, -1, Ne21, Na22, k_p_Ne21_to_Na22_reaclib, k_Na22_to_p_Ne21_reaclib,
        He4, Ne21, -1, Mg25, -1, -1, Ne21, Mg25, k_He4_Ne21_to_Mg25_reaclib, k_Mg25_to_He4_Ne21_reaclib,
        H1, Ne22, -1, Na23, -1, -1, Ne22, Na23, k_p_Ne22_to_Na23_reaclib, k_Na23_to_p_Ne22_reaclib,
        He4, Ne22, -1, Mg26, -1, -1, Ne22, Mg26, k_He4_Ne22_to_Mg26_reaclib, k_Mg26_to_He4_Ne22_reaclib,
        N, Na21, -1, Na22, -1, -1, Na21, Na22, k_n_Na21_to_Na22_reaclib, k_Na22_to_n_Na21_reaclib,
        He4, Na21, -1, Al25, -1, -1, Na21, Al25, k_He4_Na21_to_Al25_reaclib, k_Al25_to_He4_Na21_reaclib,
        N, Na22, -1, Na23, -1, -1, Na22, Na23, k_n_Na22_to_Na23_reaclib, k_Na23_to_n_Na22_reaclib,
        H1, Na22, -1, Mg23, -1, -1, Na22, Mg23, k_p_Na22_to_Mg23_reaclib, k_Mg23_to_p_Na22_reaclib,
        He4, Na22, -1, Al26, -1, -1, Na22, Al26, k_He4_Na22_to_Al26_reaclib, k_Al26_to_He4_Na22_reaclib,
        H1, Na23, -1, Mg24, -1, -1, Na23, Mg24, k_p_Na23_to_Mg24_reaclib, k_Mg24_to_p_Na23_reaclib,
        He4, Na23, -1, Al27, -1, -1, Na23, Al27, k_He4_Na23_to_Al27_reaclib, k_Al27_to_He4_Na23_reaclib,
        N, Mg23, -1, Mg24, -1, -1, Mg23, Mg24, k_n_Mg23_to_Mg24_reaclib, k_Mg24_to_n_Mg23_reaclib,
        N, Mg24, -1, Mg25, -1, -1, Mg24, Mg25, k_n_Mg24_to_Mg25_reaclib, k_Mg25_to_n_Mg24_reaclib,
        H1, Mg24, -1, Al25, -1, -1, Mg24, Al25, k_p_Mg24_to_Al25_reaclib, k_Al25_to_p_Mg24_reaclib,
        He4, Mg24, -1, Si28, -1, -1, Mg24, Si28, k_He4_Mg24_to_Si28_reaclib, k_Si28_to_He4_Mg24_reaclib,
        N, Mg25, -1, Mg26, -1, -1, Mg25, Mg26, k_n_Mg25_to_Mg26_reaclib, k_Mg26_to_n_Mg25_reaclib,
        H1, Mg25, -1, Al26, -1, -1, Mg25, Al26, k_p_Mg25_to_Al26_reaclib, k_Al26_to_p_Mg25_reaclib,
        He4, Mg25, -1, Si29, -1, -1, Mg25, Si29, k_He4_Mg25_to_Si29_reaclib, k_Si29_to_He4_Mg25_reaclib,
        H1, Mg26, -1, Al27, -1, -1, Mg26, Al27, k_p_Mg26_to_Al27_reaclib, k_Al27_to_p_Mg26_reaclib,
        He4, Mg26, -1, Si30, -1, -1, Mg26, Si30, k_He4_Mg26_to_Si30_reaclib, k_Si30_to_He4_Mg26_reaclib,
        N, Al25, -1, Al26, -1, -1, Al25, Al26, k_n_Al25_to_Al26_reaclib, k_Al26_to_n_Al25_reaclib,
        He4, Al25, -1, P29, -1, -1, Al25, P29, k_He4_Al25_to_P29_reaclib, k_P29_to_He4_Al25_reaclib,
        N, Al26, -1, Al27, -1, -1, Al26, Al27, k_n_Al26_to_Al27_reaclib, k_Al27_to_n_Al26_reaclib,
        He4, Al26, -1, P30, -1, -1, Al26, P30, k_He4_Al26_to_P30_reaclib, k_P30_to_He4_Al26_reaclib,
        H1, Al27, -1, Si28, -1, -1, Al27, Si28, k_p_Al27_to_Si28_reaclib, k_Si28_to_p_Al27_reaclib,
        He4, Al27, -1, P31, -1, -1, Al27, P31, k_He4_Al27_to_P31_reaclib, k_P31_to_He4_Al27_reaclib,
        N, Si28, -1, Si29, -1, -1, Si28, Si29, k_n_Si28_to_Si29_reaclib, k_Si29_to_n_Si28_reaclib,
        H1, Si28, -1, P29, -1, -1, Si28, P29, k_p_Si28_to_P29_reaclib, k_P29_to_p_Si28_reaclib,
        He4, Si28, -1, S32, -1, -1, Si28, S32, k_He4_Si28_to_S32_reaclib, k_S32_to_He4_Si28_reaclib,
        N, Si29, -1, Si30, -1, -1, Si29, Si30, k_n_Si29_to_Si30_reaclib, k_Si30_to_n_Si29_reaclib,
        H1, Si29, -1, P30, -1, -1, Si29, P30, k_p_Si29_to_P30_reaclib, k_P30_to_p_Si29_reaclib,
        He4, Si29, -1, S33, -1, -1, Si29, S33, k_He4_Si29_to_S33_reaclib, k_S33_to_He4_Si29_reaclib,
        N, Si30, -1, Si31, -1, -1, Si30, Si31, k_n_Si30_to_Si31_reaclib, k_Si31_to_n_Si30_reaclib,
        H1, Si30, -1, P31, -1, -1, Si30, P31, k_p_Si30_to_P31_reaclib, k_P31_to_p_Si30_reaclib,
        He4, Si30, -1, S34, -1, -1, Si30, S34, k_He4_Si30_to_S34_reaclib, k_S34_to_He4_Si30_reaclib,
        N, Si31, -1, Si32, -1, -1, Si31, Si32, k_n_Si31_to_Si32_reaclib, k_Si32_to_n_Si31_reaclib,
        H1, Si31, -1, P32, -1, -1, Si31, P32, k_p_Si31_to_P32_reaclib, k_P32_to_p_Si31_reaclib,
        He4, Si31, -1, S35, -1, -1, Si31, S35, k_He4_Si31_to_S35_reaclib, k_S35_to_He4_Si31_reaclib,
        H1, Si32, -1, P33, -1, -1, Si32, P33, k_p_Si32_to_P33_reaclib, k_P33_to_p_Si32_reaclib,
        He4, Si32, -1, S36, -1, -1, Si32, S36, k_He4_Si32_to_S36_reaclib, k_S36_to_He4_Si32_reaclib,
        N, P29, -1, P30, -1, -1, P29, P30, k_n_P29_to_P30_reaclib, k_P30_to_n_P29_reaclib,
        He4, P29, -1, Cl33, -1, -1, P29, Cl33, k_He4_P29_to_Cl33_reaclib, k_Cl33_to_He4_P29_reaclib,
        N, P30, -1, P31, -1, -1, P30, P31, k_n_P30_to_P31_reaclib, k_P31_to_n_P30_reaclib,
        He4, P30, -1, Cl34, -1, -1, P30, Cl34, k_He4_P30_to_Cl34_reaclib, k_Cl34_to_He4_P30_reaclib,
        N, P31, -1, P32, -1, -1, P31, P32, k_n_P31_to_P32_reaclib, k_P32_to_n_P31_reaclib,
        H1, P31, -1, S32, -1, -1, P31, S32, k_p_P31_to_S32_reaclib, k_S32_to_p_P31_reaclib,
        He4, P31, -1, Cl35, -1, -1, P31, Cl35, k_He4_P31_to_Cl35_reaclib, k_Cl35_to_He4_P31_reaclib,
        N, P32, -1, P33, -1, -1, P32, P33, k_n_P32_to_P33_reaclib, k_P33_to_n_P32_reaclib,
        H1, P32, -1, S33, -1, -1, P32, S33, k_p_P32_to_S33_reaclib, k_S33_to_p_P32_reaclib,
        He4, P32, -1, Cl36, -1, -1, P32, Cl36, k_He4_P32_to_Cl36_reaclib, k_Cl36_to_He4_P32_reaclib,
        H1, P33, -1, S34, -1, -1, P33, S34, k_p_P33_to_S34_reaclib, k_S34_to_p_P33_reaclib,
        He4, P33, -1, Cl37, -1, -1, P33, Cl37, k_He4_P33_to_Cl37_reaclib, k_Cl37_to_He4_P33_reaclib,
        N, S32, -1, S33, -1, -1, S32, S33, k_n_S32_to_S33_reaclib, k_S33_to_n_S32_reaclib,
        H1, S32, -1, Cl33, -1, -1, S32, Cl33, k_p_S32_to_Cl33_reaclib, k_Cl33_to_p_S32_reaclib,
        He4, S32, -1, Ar36, -1, -1, S32, Ar36, k_He4_S32_to_Ar36_reaclib, k_Ar36_to_He4_S32_reaclib,
        N, S33, -1, S34, -1, -1, S33, S34, k_n_S33_to_S34_reaclib, k_S34_to_n_S33_reaclib,
        H1, S33, -1, Cl34, -1, -1, S33, Cl34, k_p_S33_to_Cl34_reaclib, k_Cl34_to_p_S33_reaclib,
        He4, S33, -1, Ar37, -1, -1, S33, Ar37, k_He4_S33_to_Ar37_reaclib, k_Ar37_to_He4_S33_reaclib,
        N, S34, -1, S35, -1, -1, S34, S35, k_n_S34_to_S35_reaclib, k_S35_to_n_S34_reaclib,
        H1, S34, -1, Cl35, -1, -1, S34, Cl35, k_p_S34_to_Cl35_reaclib, k_Cl35_to_p_S34_reaclib,
        He4, S34, -1, Ar38, -1, -1, S34, Ar38, k_He4_S34_to_Ar38_reaclib, k_Ar38_to_He4_S34_reaclib,
        N, S35, -1, S36, -1, -1, S35, S36, k_n_S35_to_S36_reaclib, k_S36_to_n_S35_reaclib,
        H1, S35, -1, Cl36, -1, -1, S35, Cl36, k_p_S35_to_Cl36_reaclib, k_Cl36_to_p_S35_reaclib,
        He4, S35, -1, Ar39, -1, -1, S35, Ar39, k_He4_S35_to_Ar39_reaclib, k_Ar39_to_He4_S35_reaclib,
        H1, S36, -1, Cl37, -1, -1, S36, Cl37, k_p_S36_to_Cl37_reaclib, k_Cl37_to_p_S36_reaclib,
        He4, S36, -1, Ar40, -1, -1, S36, Ar40, k_He4_S36_to_Ar40_reaclib, k_Ar40_to_He4_S36_reaclib,
        N, Cl33, -1, Cl34, -1, -1, Cl33, Cl34, k_n_Cl33_to_Cl34_reaclib, k_Cl34_to_n_Cl33_reaclib,
        He4, Cl33, -1, K37, -1, -1, Cl33, K37, k_He4_Cl33_to_K37_reaclib, k_K37_to_He4_Cl33_reaclib,
        N, Cl34, -1, Cl35, -1, -1, Cl34, Cl35, k_n_Cl34_to_Cl35_reaclib, k_Cl35_to_n_Cl34_reaclib,
        He4, Cl34, -1, K38, -1, -1, Cl34, K38, k_He4_Cl34_to_K38_reaclib, k_K38_to_He4_Cl34_reaclib,
        N, Cl35, -1, Cl36, -1, -1, Cl35, Cl36, k_n_Cl35_to_Cl36_reaclib, k_Cl36_to_n_Cl35_reaclib,
        H1, Cl35, -1, Ar36, -1, -1, Cl35, Ar36, k_p_Cl35_to_Ar36_reaclib, k_Ar36_to_p_Cl35_reaclib,
        He4, Cl35, -1, K39, -1, -1, Cl35, K39, k_He4_Cl35_to_K39_reaclib, k_K39_to_He4_Cl35_reaclib,
        N, Cl36, -1, Cl37, -1, -1, Cl36, Cl37, k_n_Cl36_to_Cl37_reaclib, k_Cl37_to_n_Cl36_reaclib,
        H1, Cl36, -1, Ar37, -1, -1, Cl36, Ar37, k_p_Cl36_to_Ar37_reaclib, k_Ar37_to_p_Cl36_reaclib,
        He4, Cl36, -1, K40, -1, -1, Cl36, K40, k_He4_Cl36_to_K40_reaclib, k_K40_to_He4_Cl36_reaclib,
        H1, Cl37, -1, Ar38, -1, -1, Cl37, Ar38, k_p_Cl37_to_Ar38_reaclib, k_Ar38_to_p_Cl37_reaclib,
        He4, Cl37, -1, K41, -1, -1, Cl37, K41, k_He4_Cl37_to_K41_reaclib, k_K41_to_He4_Cl37_reaclib,
        N, Ar36, -1, Ar37, -1, -1, Ar36, Ar37, k_n_Ar36_to_Ar37_reaclib, k_Ar37_to_n_Ar36_reaclib,
        H1, Ar36, -1, K37, -1, -1, Ar36, K37, k_p_Ar36_to_K37_reaclib, k_K37_to_p_Ar36_reaclib,
        He4, Ar36, -1, Ca40, -1, -1, Ar36, Ca40, k_He4_Ar36_to_Ca40_reaclib, k_Ca40_to_He4_Ar36_reaclib,
        N, Ar37, -1, Ar38, -1, -1, Ar37, Ar38, k_n_Ar37_to_Ar38_reaclib, k_Ar38_to_n_Ar37_reaclib,
        H1, Ar37, -1, K38, -1, -1, Ar37, K38, k_p_Ar37_to_K38_reaclib, k_K38_to_p_Ar37_reaclib,
        He4, Ar37, -1, Ca41, -1, -1, Ar37, Ca41, k_He4_Ar37_to_Ca41_reaclib, k_Ca41_to_He4_Ar37_reaclib,
        N, Ar38, -1, Ar39, -1, -1, Ar38, Ar39, k_n_Ar38_to_Ar39_reaclib, k_Ar39_to_n_Ar38_reaclib,
        H1, Ar38, -1, K39, -1, -1, Ar38, K39, k_p_Ar38_to_K39_reaclib, k_K39_to_p_Ar38_reaclib,
        He4, Ar38, -1, Ca42, -1, -1, Ar38, Ca42, k_He4_Ar38_to_Ca42_reaclib, k_Ca42_to_He4_Ar38_reaclib,
        N, Ar39, -1, Ar40, -1, -1, Ar39, Ar40, k_n_Ar39_to_Ar40_reaclib, k_Ar40_to_n_Ar39_reaclib,
        H1, Ar39, -1, K40, -1, -1, Ar39, K40, k_p_Ar39_to_K40_reaclib, k_K40_to_p_Ar39_reaclib,
        He4, Ar39, -1, Ca43, -1, -1, Ar39, Ca43, k_He4_Ar39_to_Ca43_reaclib, k_Ca43_to_He4_Ar39_reaclib,
        H1, Ar40, -1, K41, -1, -1, Ar40, K41, k_p_Ar40_to_K41_reaclib, k_K41_to_p_Ar40_reaclib,
        He4, Ar40, -1, Ca44, -1, -1, Ar40, Ca44, k_He4_Ar40_to_Ca44_reaclib, k_Ca44_to_He4_Ar40_reaclib,
        N, K37, -1, K38, -1, -1, K37, K38, k_n_K37_to_K38_reaclib, k_K38_to_n_K37_reaclib,
        N, K38, -1, K39, -1, -1, K38, K39, k_n_K38_to_K39_reaclib, k_K39_to_n_K38_reaclib,
        N, K39, -1, K40, -1, -1, K39, K40, k_n_K39_to_K40_reaclib, k_K40_to_n_K39_reaclib,
        H1, K39, -1, Ca40, -1, -1, K39, Ca40, k_p_K39_to_Ca40_reaclib, k_Ca40_to_p_K39_reaclib,
        He4, K39, -1, Sc43, -1, -1, K39, Sc43, k_He4_K39_to_Sc43_reaclib, k_Sc43_to_He4_K39_reaclib,
        N, K40, -1, K41, -1, -1, K40, K41, k_n_K40_to_K41_reaclib, k_K41_to_n_K40_reaclib,
        H1, K40, -1, Ca41, -1, -1, K40, Ca41, k_p_K40_to_Ca41_reaclib, k_Ca41_to_p_K40_reaclib,
        He4, K40, -1, Sc44, -1, -1, K40, Sc44, k_He4_K40_to_Sc44_reaclib, k_Sc44_to_He4_K40_reaclib,
        H1, K41, -1, Ca42, -1, -1, K41, Ca42, k_p_K41_to_Ca42_reaclib, k_Ca42_to_p_K41_reaclib,
        He4, K41, -1, Sc45, -1, -1, K41, Sc45, k_He4_K41_to_Sc45_reaclib, k_Sc45_to_He4_K41_reaclib,
        N, Ca40, -1, Ca41, -1, -1, Ca40, Ca41, k_n_Ca40_to_Ca41_reaclib, k_Ca41_to_n_Ca40_reaclib,
        He4, Ca40, -1, Ti44, -1, -1, Ca40, Ti44, k_He4_Ca40_to_Ti44_reaclib, k_Ti44_to_He4_Ca40_reaclib,
        N, Ca41, -1, Ca42, -1, -1, Ca41, Ca42, k_n_Ca41_to_Ca42_reaclib, k_Ca42_to_n_Ca41_reaclib,
        He4, Ca41, -1, Ti45, -1, -1, Ca41, Ti45, k_He4_Ca41_to_Ti45_reaclib, k_Ti45_to_He4_Ca41_reaclib,
        N, Ca42, -1, Ca43, -1, -1, Ca42, Ca43, k_n_Ca42_to_Ca43_reaclib, k_Ca43_to_n_Ca42_reaclib,
        H1, Ca42, -1, Sc43, -1, -1, Ca42, Sc43, k_p_Ca42_to_Sc43_reaclib, k_Sc43_to_p_Ca42_reaclib,
        He4, Ca42, -1, Ti46, -1, -1, Ca42, Ti46, k_He4_Ca42_to_Ti46_reaclib, k_Ti46_to_He4_Ca42_reaclib,
        N, Ca43, -1, Ca44, -1, -1, Ca43, Ca44, k_n_Ca43_to_Ca44_reaclib, k_Ca44_to_n_Ca43_reaclib,
        H1, Ca43, -1, Sc44, -1, -1, Ca43, Sc44, k_p_Ca43_to_Sc44_reaclib, k_Sc44_to_p_Ca43_reaclib,
        He4, Ca43, -1, Ti47, -1, -1, Ca43, Ti47, k_He4_Ca43_to_Ti47_reaclib, k_Ti47_to_He4_Ca43_reaclib,
        N, Ca44, -1, Ca45, -1, -1, Ca44, Ca45, k_n_Ca44_to_Ca45_reaclib, k_Ca45_to_n_Ca44_reaclib,
        H1, Ca44, -1, Sc45, -1, -1, Ca44, Sc45, k_p_Ca44_to_Sc45_reaclib, k_Sc45_to_p_Ca44_reaclib,
        He4, Ca44, -1, Ti48, -1, -1, Ca44, Ti48, k_He4_Ca44_to_Ti48_reaclib, k_Ti48_to_He4_Ca44_reaclib,
        N, Ca45, -1, Ca46, -1, -1, Ca45, Ca46, k_n_Ca45_to_Ca46_reaclib, k_Ca46_to_n_Ca45_reaclib,
        H1, Ca45, -1, Sc46, -1, -1, Ca45, Sc46, k_p_Ca45_to_Sc46_reaclib, k_Sc46_to_p_Ca45_reaclib,
        He4, Ca45, -1, Ti49, -1, -1, Ca45, Ti49, k_He4_Ca45_to_Ti49_reaclib, k_Ti49_to_He4_Ca45_reaclib,
        N, Ca46, -1, Ca47, -1, -1, Ca46, Ca47, k_n_Ca46_to_Ca47_reaclib, k_Ca47_to_n_Ca46_reaclib,
        H1, Ca46, -1, Sc47, -1, -1, Ca46, Sc47, k_p_Ca46_to_Sc47_reaclib, k_Sc47_to_p_Ca46_reaclib,
        He4, Ca46, -1, Ti50, -1, -1, Ca46, Ti50, k_He4_Ca46_to_Ti50_reaclib, k_Ti50_to_He4_Ca46_reaclib,
        N, Ca47, -1, Ca48, -1, -1, Ca47, Ca48, k_n_Ca47_to_Ca48_reaclib, k_Ca48_to_n_Ca47_reaclib,
        H1, Ca47, -1, Sc48, -1, -1, Ca47, Sc48, k_p_Ca47_to_Sc48_reaclib, k_Sc48_to_p_Ca47_reaclib,
        He4, Ca47, -1, Ti51, -1, -1, Ca47, Ti51, k_He4_Ca47_to_Ti51_reaclib, k_Ti51_to_He4_Ca47_reaclib,
        H1, Ca48, -1, Sc49, -1, -1, Ca48, Sc49, k_p_Ca48_to_Sc49_reaclib, k_Sc49_to_p_Ca48_reaclib,
        N, Sc43, -1, Sc44, -1, -1, Sc43, Sc44, k_n_Sc43_to_Sc44_reaclib, k_Sc44_to_n_Sc43_reaclib,
        H1, Sc43, -1, Ti44, -1, -1, Sc43, Ti44, k_p_Sc43_to_Ti44_reaclib, k_Ti44_to_p_Sc43_reaclib,
        He4, Sc43, -1, V47, -1, -1, Sc43, V47, k_He4_Sc43_to_V47_reaclib, k_V47_to_He4_Sc43_reaclib,
        N, Sc44, -1, Sc45, -1, -1, Sc44, Sc45, k_n_Sc44_to_Sc45_reaclib, k_Sc45_to_n_Sc44_reaclib,
        H1, Sc44, -1, Ti45, -1, -1, Sc44, Ti45, k_p_Sc44_to_Ti45_reaclib, k_Ti45_to_p_Sc44_reaclib,
        He4, Sc44, -1, V48, -1, -1, Sc44, V48, k_He4_Sc44_to_V48_reaclib, k_V48_to_He4_Sc44_reaclib,
        N, Sc45, -1, Sc46, -1, -1, Sc45, Sc46, k_n_Sc45_to_Sc46_reaclib, k_Sc46_to_n_Sc45_reaclib,
        H1, Sc45, -1, Ti46, -1, -1, Sc45, Ti46, k_p_Sc45_to_Ti46_reaclib, k_Ti46_to_p_Sc45_reaclib,
        He4, Sc45, -1, V49, -1, -1, Sc45, V49, k_He4_Sc45_to_V49_reaclib, k_V49_to_He4_Sc45_reaclib,
        N, Sc46, -1, Sc47, -1, -1, Sc46, Sc47, k_n_Sc46_to_Sc47_reaclib, k_Sc47_to_n_Sc46_reaclib,
        H1, Sc46, -1, Ti47, -1, -1, Sc46, Ti47, k_p_Sc46_to_Ti47_reaclib, k_Ti47_to_p_Sc46_reaclib,
        He4, Sc46, -1, V50, -1, -1, Sc46, V50, k_He4_Sc46_to_V50_reaclib, k_V50_to_He4_Sc46_reaclib,
        N, Sc47, -1, Sc48, -1, -1, Sc47, Sc48, k_n_Sc47_to_Sc48_reaclib, k_Sc48_to_n_Sc47_reaclib,
        H1, Sc47, -1, Ti48, -1, -1, Sc47, Ti48, k_p_Sc47_to_Ti48_reaclib, k_Ti48_to_p_Sc47_reaclib,
        He4, Sc47, -1, V51, -1, -1, Sc47, V51, k_He4_Sc47_to_V51_reaclib, k_V51_to_He4_Sc47_reaclib,
        N, Sc48, -1, Sc49, -1, -1, Sc48, Sc49, k_n_Sc48_to_Sc49_reaclib, k_Sc49_to_n_Sc48_reaclib,
        H1, Sc48, -1, Ti49, -1, -1, Sc48, Ti49, k_p_Sc48_to_Ti49_reaclib, k_Ti49_to_p_Sc48_reaclib,
        He4, Sc48, -1, V52, -1, -1, Sc48, V52, k_He4_Sc48_to_V52_reaclib, k_V52_to_He4_Sc48_reaclib,
        H1, Sc49, -1, Ti50, -1, -1, Sc49, Ti50, k_p_Sc49_to_Ti50_reaclib, k_Ti50_to_p_Sc49_reaclib,
        N, Ti44, -1, Ti45, -1, -1, Ti44, Ti45, k_n_Ti44_to_Ti45_reaclib, k_Ti45_to_n_Ti44_reaclib,
        He4, Ti44, -1, Cr48, -1, -1, Ti44, Cr48, k_He4_Ti44_to_Cr48_reaclib, k_Cr48_to_He4_Ti44_reaclib,
        N, Ti45, -1, Ti46, -1, -1, Ti45, Ti46, k_n_Ti45_to_Ti46_reaclib, k_Ti46_to_n_Ti45_reaclib,
        H1, Ti45, -1, V46, -1, -1, Ti45, V46, k_p_Ti45_to_V46_reaclib, k_V46_to_p_Ti45_reaclib,
        He4, Ti45, -1, Cr49, -1, -1, Ti45, Cr49, k_He4_Ti45_to_Cr49_reaclib, k_Cr49_to_He4_Ti45_reaclib,
        N, Ti46, -1, Ti47, -1, -1, Ti46, Ti47, k_n_Ti46_to_Ti47_reaclib, k_Ti47_to_n_Ti46_reaclib,
        H1, Ti46, -1, V47, -1, -1, Ti46, V47, k_p_Ti46_to_V47_reaclib, k_V47_to_p_Ti46_reaclib,
        He4, Ti46, -1, Cr50, -1, -1, Ti46, Cr50, k_He4_Ti46_to_Cr50_reaclib, k_Cr50_to_He4_Ti46_reaclib,
        N, Ti47, -1, Ti48, -1, -1, Ti47, Ti48, k_n_Ti47_to_Ti48_reaclib, k_Ti48_to_n_Ti47_reaclib,
        H1, Ti47, -1, V48, -1, -1, Ti47, V48, k_p_Ti47_to_V48_reaclib, k_V48_to_p_Ti47_reaclib,
        He4, Ti47, -1, Cr51, -1, -1, Ti47, Cr51, k_He4_Ti47_to_Cr51_reaclib, k_Cr51_to_He4_Ti47_reaclib,
        N, Ti48, -1, Ti49, -1, -1, Ti48, Ti49, k_n_Ti48_to_Ti49_reaclib, k_Ti49_to_n_Ti48_reaclib,
        H1, Ti48, -1, V49, -1, -1, Ti48, V49, k_p_Ti48_to_V49_reaclib, k_V49_to_p_Ti48_reaclib,
        He4, Ti48, -1, Cr52, -1, -1, Ti48, Cr52, k_He4_Ti48_to_Cr52_reaclib, k_Cr52_to_He4_Ti48_reaclib,
        N, Ti49, -1, Ti50, -1, -1, Ti49, Ti50, k_n_Ti49_to_Ti50_reaclib, k_Ti50_to_n_Ti49_reaclib,
        H1, Ti49, -1, V50, -1, -1, Ti49, V50, k_p_Ti49_to_V50_reaclib, k_V50_to_p_Ti49_reaclib,
        He4, Ti49, -1, Cr53, -1, -1, Ti49, Cr53, k_He4_Ti49_to_Cr53_reaclib, k_Cr53_to_He4_Ti49_reaclib,
        N, Ti50, -1, Ti51, -1, -1, Ti50, Ti51, k_n_Ti50_to_Ti51_reaclib, k_Ti51_to_n_Ti50_reaclib,
        H1, Ti50, -1, V51, -1, -1, Ti50, V51, k_p_Ti50_to_V51_reaclib, k_V51_to_p_Ti50_reaclib,
        He4, Ti50, -1, Cr54, -1, -1, Ti50, Cr54, k_He4_Ti50_to_Cr54_reaclib, k_Cr54_to_He4_Ti50_reaclib,
        H1, Ti51, -1, V52, -1, -1, Ti51, V52, k_p_Ti51_to_V52_reaclib, k_V52_to_p_Ti51_reaclib,
        N, V46, -1, V47, -1, -1, V46, V47, k_n_V46_to_V47_reaclib, k_V47_to_n_V46_reaclib,
        He4, V46, -1, Mn50, -1, -1, V46, Mn50, k_He4_V46_to_Mn50_reaclib, k_Mn50_to_He4_V46_reaclib,
        N, V47, -1, V48, -1, -1, V47, V48, k_n_V47_to_V48_reaclib, k_V48_to_n_V47_reaclib,
        H1, V47, -1, Cr48, -1, -1, V47, Cr48, k_p_V47_to_Cr48_reaclib, k_Cr48_to_p_V47_reaclib,
        He4, V47, -1, Mn51, -1, -1, V47, Mn51, k_He4_V47_to_Mn51_reaclib, k_Mn51_to_He4_V47_reaclib,
        N, V48, -1, V49, -1, -1, V48, V49, k_n_V48_to_V49_reaclib, k_V49_to_n_V48_reaclib,
        H1, V48, -1, Cr49, -1, -1, V48, Cr49, k_p_V48_to_Cr49_reaclib, k_Cr49_to_p_V48_reaclib,
        He4, V48, -1, Mn52, -1, -1, V48, Mn52, k_He4_V48_to_Mn52_reaclib, k_Mn52_to_He4_V48_reaclib,
        N, V49, -1, V50, -1, -1, V49, V50, k_n_V49_to_V50_reaclib, k_V50_to_n_V49_reaclib,
        H1, V49, -1, Cr50, -1, -1, V49, Cr50, k_p_V49_to_Cr50_reaclib, k_Cr50_to_p_V49_reaclib,
        He4, V49, -1, Mn53, -1, -1, V49, Mn53, k_He4_V49_to_Mn53_reaclib, k_Mn53_to_He4_V49_reaclib,
        N, V50, -1, V51, -1, -1, V50, V51, k_n_V50_to_V51_reaclib, k_V51_to_n_V50_reaclib,
        H1, V50, -1, Cr51, -1, -1, V50, Cr51, k_p_V50_to_Cr51_reaclib, k_Cr51_to_p_V50_reaclib,
        He4, V50, -1, Mn54, -1, -1, V50, Mn54, k_He4_V50_to_Mn54_reaclib, k_Mn54_to_He4_V50_reaclib,
        N, V51, -1, V52, -1, -1, V51, V52, k_n_V51_to_V52_reaclib, k_V52_to_n_V51_reaclib,
        H1, V51, -1, Cr52, -1, -1, V51, Cr52, k_p_V51_to_Cr52_reaclib, k_Cr52_to_p_V51_reaclib,
        He4, V51, -1, Mn55, -1, -1, V51, Mn55, k_He4_V51_to_Mn55_reaclib, k_Mn55_to_He4_V51_reaclib,
        H1, V52, -1, Cr53, -1, -1, V52, Cr53, k_p_V52_to_Cr53_reaclib, k_Cr53_to_p_V52_reaclib,
        N, Cr48, -1, Cr49, -1, -1, Cr48, Cr49, k_n_Cr48_to_Cr49_reaclib, k_Cr49_to_n_Cr48_reaclib,
        He4, Cr48, -1, Fe52, -1, -1, Cr48, Fe52, k_He4_Cr48_to_Fe52_reaclib, k_Fe52_to_He4_Cr48_reaclib,
        N, Cr49, -1, Cr50, -1, -1, Cr49, Cr50, k_n_Cr49_to_Cr50_reaclib, k_Cr50_to_n_Cr49_reaclib,
        H1, Cr49, -1, Mn50, -1, -1, Cr49, Mn50, k_p_Cr49_to_Mn50_reaclib, k_Mn50_to_p_Cr49_reaclib,
        He4, Cr49, -1, Fe53, -1, -1, Cr49, Fe53, k_He4_Cr49_to_Fe53_reaclib, k_Fe53_to_He4_Cr49_reaclib,
        N, Cr50, -1, Cr51, -1, -1, Cr50, Cr51, k_n_Cr50_to_Cr51_reaclib, k_Cr51_to_n_Cr50_reaclib,
        H1, Cr50, -1, Mn51, -1, -1, Cr50, Mn51, k_p_Cr50_to_Mn51_reaclib, k_Mn51_to_p_Cr50_reaclib,
        He4, Cr50, -1, Fe54, -1, -1, Cr50, Fe54, k_He4_Cr50_to_Fe54_reaclib, k_Fe54_to_He4_Cr50_reaclib,
        N, Cr51, -1, Cr52, -1, -1, Cr51, Cr52, k_n_Cr51_to_Cr52_reaclib, k_Cr52_to_n_Cr51_reaclib,
        H1, Cr51, -1, Mn52, -1, -1, Cr51, Mn52, k_p_Cr51_to_Mn52_reaclib, k_Mn52_to_p_Cr51_reaclib,
        He4, Cr51, -1, Fe55, -1, -1, Cr51, Fe55, k_He4_Cr51_to_Fe55_reaclib, k_Fe55_to_He4_Cr51_reaclib,
        N, Cr52, -1, Cr53, -1, -1, Cr52, Cr53, k_n_Cr52_to_Cr53_reaclib, k_Cr53_to_n_Cr52_reaclib,
        H1, Cr52, -1, Mn53, -1, -1, Cr52, Mn53, k_p_Cr52_to_Mn53_reaclib, k_Mn53_to_p_Cr52_reaclib,
        He4, Cr52, -1, Fe56, -1, -1, Cr52, Fe56, k_He4_Cr52_to_Fe56_reaclib, k_Fe56_to_He4_Cr52_reaclib,
        N, Cr53, -1, Cr54, -1, -1, Cr53, Cr54, k_n_Cr53_to_Cr54_reaclib, k_Cr54_to_n_Cr53_reaclib,
        H1, Cr53, -1, Mn54, -1, -1, Cr53, Mn54, k_p_Cr53_to_Mn54_reaclib, k_Mn54_to_p_Cr53_reaclib,
        He4, Cr53, -1, Fe57, -1, -1, Cr53, Fe57, k_He4_Cr53_to_Fe57_reaclib, k_Fe57_to_He4_Cr53_reaclib,
        H1, Cr54, -1, Mn55, -1, -1, Cr54, Mn55, k_p_Cr54_to_Mn55_reaclib, k_Mn55_to_p_Cr54_reaclib,
        He4, Cr54, -1, Fe58, -1, -1, Cr54, Fe58, k_He4_Cr54_to_Fe58_reaclib, k_Fe58_to_He4_Cr54_reaclib,
        N, Mn50, -1, Mn51, -1, -1, Mn50, Mn51, k_n_Mn50_to_Mn51_reaclib, k_Mn51_to_n_Mn50_reaclib,
        He4, Mn50, -1, Co54, -1, -1, Mn50, Co54, k_He4_Mn50_to_Co54_reaclib, k_Co54_to_He4_Mn50_reaclib,
        N, Mn51, -1, Mn52, -1, -1, Mn51, Mn52, k_n_Mn51_to_Mn52_reaclib, k_Mn52_to_n_Mn51_reaclib,
        H1, Mn51, -1, Fe52, -1, -1, Mn51, Fe52, k_p_Mn51_to_Fe52_reaclib, k_Fe52_to_p_Mn51_reaclib,
        He4, Mn51, -1, Co55, -1, -1, Mn51, Co55, k_He4_Mn51_to_Co55_reaclib, k_Co55_to_He4_Mn51_reaclib,
        N, Mn52, -1, Mn53, -1, -1, Mn52, Mn53, k_n_Mn52_to_Mn53_reaclib, k_Mn53_to_n_Mn52_reaclib,
        H1, Mn52, -1, Fe53, -1, -1, Mn52, Fe53, k_p_Mn52_to_Fe53_reaclib, k_Fe53_to_p_Mn52_reaclib,
        He4, Mn52, -1, Co56, -1, -1, Mn52, Co56, k_He4_Mn52_to_Co56_reaclib, k_Co56_to_He4_Mn52_reaclib,
        N, Mn53, -1, Mn54, -1, -1, Mn53, Mn54, k_n_Mn53_to_Mn54_reaclib, k_Mn54_to_n_Mn53_reaclib,
        H1, Mn53, -1, Fe54, -1, -1, Mn53, Fe54, k_p_Mn53_to_Fe54_reaclib, k_Fe54_to_p_Mn53_reaclib,
        He4, Mn53, -1, Co57, -1, -1, Mn53, Co57, k_He4_Mn53_to_Co57_reaclib, k_Co57_to_He4_Mn53_reaclib,
        N, Mn54, -1, Mn55, -1, -1, Mn54, Mn55, k_n_Mn54_to_Mn55_reaclib, k_Mn55_to_n_Mn54_reaclib,
        H1, Mn54, -1, Fe55, -1, -1, Mn54, Fe55, k_p_Mn54_to_Fe55_reaclib, k_Fe55_to_p_Mn54_reaclib,
        He4, Mn54, -1, Co58, -1, -1, Mn54, Co58, k_He4_Mn54_to_Co58_reaclib, k_Co58_to_He4_Mn54_reaclib,
        H1, Mn55, -1, Fe56, -1, -1, Mn55, Fe56, k_p_Mn55_to_Fe56_reaclib, k_Fe56_to_p_Mn55_reaclib,
        He4, Mn55, -1, Co59, -1, -1, Mn55, Co59, k_He4_Mn55_to_Co59_reaclib, k_Co59_to_He4_Mn55_reaclib,
        N, Fe52, -1, Fe53, -1, -1, Fe52, Fe53, k_n_Fe52_to_Fe53_reaclib, k_Fe53_to_n_Fe52_reaclib,
        H1, Fe52, -1, Co53, -1, -1, Fe52, Co53, k_p_Fe52_to_Co53_reaclib, k_Co53_to_p_Fe52_reaclib,
        He4, Fe52, -1, Ni56, -1, -1, Fe52, Ni56, k_He4_Fe52_to_Ni56_reaclib, k_Ni56_to_He4_Fe52_reaclib,
        N, Fe53, -1, Fe54, -1, -1, Fe53, Fe54, k_n_Fe53_to_Fe54_reaclib, k_Fe54_to_n_Fe53_reaclib,
        H1, Fe53, -1, Co54, -1, -1, Fe53, Co54, k_p_Fe53_to_Co54_reaclib, k_Co54_to_p_Fe53_reaclib,
        He4, Fe53, -1, Ni57, -1, -1, Fe53, Ni57, k_He4_Fe53_to_Ni57_reaclib, k_Ni57_to_He4_Fe53_reaclib,
        N, Fe54, -1, Fe55, -1, -1, Fe54, Fe55, k_n_Fe54_to_Fe55_reaclib, k_Fe55_to_n_Fe54_reaclib,
        H1, Fe54, -1, Co55, -1, -1, Fe54, Co55, k_p_Fe54_to_Co55_reaclib, k_Co55_to_p_Fe54_reaclib,
        He4, Fe54, -1, Ni58, -1, -1, Fe54, Ni58, k_He4_Fe54_to_Ni58_reaclib, k_Ni58_to_He4_Fe54_reaclib,
        N, Fe55, -1, Fe56, -1, -1, Fe55, Fe56, k_n_Fe55_to_Fe56_reaclib, k_Fe56_to_n_Fe55_reaclib,
        H1, Fe55, -1, Co56, -1, -1, Fe55, Co56, k_p_Fe55_to_Co56_reaclib, k_Co56_to_p_Fe55_reaclib,
        He4, Fe55, -1, Ni59, -1, -1, Fe55, Ni59, k_He4_Fe55_to_Ni59_reaclib, k_Ni59_to_He4_Fe55_reaclib,
        N, Fe56, -1, Fe57, -1, -1, Fe56, Fe57, k_n_Fe56_to_Fe57_reaclib, k_Fe57_to_n_Fe56_reaclib,
        H1, Fe56, -1, Co57, -1, -1, Fe56, Co57, k_p_Fe56_to_Co57_reaclib, k_Co57_to_p_Fe56_reaclib,
        He4, Fe56, -1, Ni60, -1, -1, Fe56, Ni60, k_He4_Fe56_to_Ni60_reaclib, k_Ni60_to_He4_Fe56_reaclib,
        N, Fe57, -1, Fe58, -1, -1, Fe57, Fe58, k_n_Fe57_to_Fe58_reaclib, k_Fe58_to_n_Fe57_reaclib,
        H1, Fe57, -1, Co58, -1, -1, Fe57, Co58, k_p_Fe57_to_Co58_reaclib, k_Co58_to_p_Fe57_reaclib,
        He4, Fe57, -1, Ni61, -1, -1, Fe57, Ni61, k_He4_Fe57_to_Ni61_reaclib, k_Ni61_to_He4_Fe57_reaclib,
        H1, Fe58, -1, Co59, -1, -1, Fe58, Co59, k_p_Fe58_to_Co59_reaclib, k_Co59_to_p_Fe58_reaclib,
        He4, Fe58, -1, Ni62, -1, -1, Fe58, Ni62, k_He4_Fe58_to_Ni62_reaclib, k_Ni62_to_He4_Fe58_reaclib,
        N, Co53, -1, Co54, -1, -1, Co53, Co54, k_n_Co53_to_Co54_reaclib, k_Co54_to_n_Co53_reaclib,
        He4, Co53, -1, Cu57, -1, -1, Co53, Cu57, k_He4_Co53_to_Cu57_reaclib, k_Cu57_to_He4_Co53_reaclib,
        N, Co54, -1, Co55, -1, -1, Co54, Co55, k_n_Co54_to_Co55_reaclib, k_Co55_to_n_Co54_reaclib,
        He4, Co54, -1, Cu58, -1, -1, Co54, Cu58, k_He4_Co54_to_Cu58_reaclib, k_Cu58_to_He4_Co54_reaclib,
        N, Co55, -1, Co56, -1, -1, Co55, Co56, k_n_Co55_to_Co56_reaclib, k_Co56_to_n_Co55_reaclib,
        H1, Co55, -1, Ni56, -1, -1, Co55, Ni56, k_p_Co55_to_Ni56_reaclib, k_Ni56_to_p_Co55_reaclib,
        He4, Co55, -1, Cu59, -1, -1, Co55, Cu59, k_He4_Co55_to_Cu59_reaclib, k_Cu59_to_He4_Co55_reaclib,
        N, Co56, -1, Co57, -1, -1, Co56, Co57, k_n_Co56_to_Co57_reaclib, k_Co57_to_n_Co56_reaclib,
        H1, Co56, -1, Ni57, -1, -1, Co56, Ni57, k_p_Co56_to_Ni57_reaclib, k_Ni57_to_p_Co56_reaclib,
        He4, Co56, -1, Cu60, -1, -1, Co56, Cu60, k_He4_Co56_to_Cu60_reaclib, k_Cu60_to_He4_Co56_reaclib,
        N, Co57, -1, Co58, -1, -1, Co57, Co58, k_n_Co57_to_Co58_reaclib, k_Co58_to_n_Co57_reaclib,
        H1, Co57, -1, Ni58, -1, -1, Co57, Ni58, k_p_Co57_to_Ni58_reaclib, k_Ni58_to_p_Co57_reaclib,
        He4, Co57, -1, Cu61, -1, -1, Co57, Cu61, k_He4_Co57_to_Cu61_reaclib, k_Cu61_to_He4_Co57_reaclib,
        N, Co58, -1, Co59, -1, -1, Co58, Co59, k_n_Co58_to_Co59_reaclib, k_Co59_to_n_Co58_reaclib,
        H1, Co58, -1, Ni59, -1, -1, Co58, Ni59, k_p_Co58_to_Ni59_reaclib, k_Ni59_to_p_Co58_reaclib,
        He4, Co58, -1, Cu62, -1, -1, Co58, Cu62, k_He4_Co58_to_Cu62_reaclib, k_Cu62_to_He4_Co58_reaclib,
        H1, Co59, -1, Ni60, -1, -1, Co59, Ni60, k_p_Co59_to_Ni60_reaclib, k_Ni60_to_p_Co59_reaclib,
        He4, Co59, -1, Cu63, -1, -1, Co59, Cu63, k_He4_Co59_to_Cu63_reaclib, k_Cu63_to_He4_Co59_reaclib,
        N, Ni56, -1, Ni57, -1, -1, Ni56, Ni57, k_n_Ni56_to_Ni57_reaclib, k_Ni57_to_n_Ni56_reaclib,
        H1, Ni56, -1, Cu57, -1, -1, Ni56, Cu57, k_p_Ni56_to_Cu57_reaclib, k_Cu57_to_p_Ni56_reaclib,
        He4, Ni56, -1, Zn60, -1, -1, Ni56, Zn60, k_He4_Ni56_to_Zn60_reaclib, k_Zn60_to_He4_Ni56_reaclib,
        N, Ni57, -1, Ni58, -1, -1, Ni57, Ni58, k_n_Ni57_to_Ni58_reaclib, k_Ni58_to_n_Ni57_reaclib,
        H1, Ni57, -1, Cu58, -1, -1, Ni57, Cu58, k_p_Ni57_to_Cu58_reaclib, k_Cu58_to_p_Ni57_reaclib,
        He4, Ni57, -1, Zn61, -1, -1, Ni57, Zn61, k_He4_Ni57_to_Zn61_reaclib, k_Zn61_to_He4_Ni57_reaclib,
        N, Ni58, -1, Ni59, -1, -1, Ni58, Ni59, k_n_Ni58_to_Ni59_reaclib, k_Ni59_to_n_Ni58_reaclib,
        H1, Ni58, -1, Cu59, -1, -1, Ni58, Cu59, k_p_Ni58_to_Cu59_reaclib, k_Cu59_to_p_Ni58_reaclib,
        He4, Ni58, -1, Zn62, -1, -1, Ni58, Zn62, k_He4_Ni58_to_Zn62_reaclib, k_Zn62_to_He4_Ni58_reaclib,
        N, Ni59, -1, Ni60, -1, -1, Ni59, Ni60, k_n_Ni59_to_Ni60_reaclib, k_Ni60_to_n_Ni59_reaclib,
        H1, Ni59, -1, Cu60, -1, -1, Ni59, Cu60, k_p_Ni59_to_Cu60_reaclib, k_Cu60_to_p_Ni59_reaclib,
        He4, Ni59, -1, Zn63, -1, -1, Ni59, Zn63, k_He4_Ni59_to_Zn63_reaclib, k_Zn63_to_He4_Ni59_reaclib,
        N, Ni60, -1, Ni61, -1, -1, Ni60, Ni61, k_n_Ni60_to_Ni61_reaclib, k_Ni61_to_n_Ni60_reaclib,
        H1, Ni60, -1, Cu61, -1, -1, Ni60, Cu61, k_p_Ni60_to_Cu61_reaclib, k_Cu61_to_p_Ni60_reaclib,
        He4, Ni60, -1, Zn64, -1, -1, Ni60, Zn64, k_He4_Ni60_to_Zn64_reaclib, k_Zn64_to_He4_Ni60_reaclib,
        N, Ni61, -1, Ni62, -1, -1, Ni61, Ni62, k_n_Ni61_to_Ni62_reaclib, k_Ni62_to_n_Ni61_reaclib,
        H1, Ni61, -1, Cu62, -1, -1, Ni61, Cu62, k_p_Ni61_to_Cu62_reaclib, k_Cu62_to_p_Ni61_reaclib,
        He4, Ni61, -1, Zn65, -1, -1, Ni61, Zn65, k_He4_Ni61_to_Zn65_reaclib, k_Zn65_to_He4_Ni61_reaclib,
        N, Ni62, -1, Ni63, -1, -1, Ni62, Ni63, k_n_Ni62_to_Ni63_reaclib, k_Ni63_to_n_Ni62_reaclib,
        H1, Ni62, -1, Cu63, -1, -1, Ni62, Cu63, k_p_Ni62_to_Cu63_reaclib, k_Cu63_to_p_Ni62_reaclib,
        He4, Ni62, -1, Zn66, -1, -1, Ni62, Zn66, k_He4_Ni62_to_Zn66_reaclib, k_Zn66_to_He4_Ni62_reaclib,
        N, Ni63, -1, Ni64, -1, -1, Ni63, Ni64, k_n_Ni63_to_Ni64_reaclib, k_Ni64_to_n_Ni63_reaclib,
        H1, Ni63, -1, Cu64, -1, -1, Ni63, Cu64, k_p_Ni63_to_Cu64_reaclib, k_Cu64_to_p_Ni63_reaclib,
        H1, Ni64, -1, Cu65, -1, -1, Ni64, Cu65, k_p_Ni64_to_Cu65_reaclib, k_Cu65_to_p_Ni64_reaclib,
        N, Cu57, -1, Cu58, -1, -1, Cu57, Cu58, k_n_Cu57_to_Cu58_reaclib, k_Cu58_to_n_Cu57_reaclib,
        N, Cu58, -1, Cu59, -1, -1, Cu58, Cu59, k_n_Cu58_to_Cu59_reaclib, k_Cu59_to_n_Cu58_reaclib,
        H1, Cu58, -1, Zn59, -1, -1, Cu58, Zn59, k_p_Cu58_to_Zn59_reaclib, k_Zn59_to_p_Cu58_reaclib,
        He4, Cu58, -1, Ga62, -1, -1, Cu58, Ga62, k_He4_Cu58_to_Ga62_reaclib, k_Ga62_to_He4_Cu58_reaclib,
        N, Cu59, -1, Cu60, -1, -1, Cu59, Cu60, k_n_Cu59_to_Cu60_reaclib, k_Cu60_to_n_Cu59_reaclib,
        H1, Cu59, -1, Zn60, -1, -1, Cu59, Zn60, k_p_Cu59_to_Zn60_reaclib, k_Zn60_to_p_Cu59_reaclib,
        He4, Cu59, -1, Ga63, -1, -1, Cu59, Ga63, k_He4_Cu59_to_Ga63_reaclib, k_Ga63_to_He4_Cu59_reaclib,
        N, Cu60, -1, Cu61, -1, -1, Cu60, Cu61, k_n_Cu60_to_Cu61_reaclib, k_Cu61_to_n_Cu60_reaclib,
        H1, Cu60, -1, Zn61, -1, -1, Cu60, Zn61, k_p_Cu60_to_Zn61_reaclib, k_Zn61_to_p_Cu60_reaclib,
        He4, Cu60, -1, Ga64, -1, -1, Cu60, Ga64, k_He4_Cu60_to_Ga64_reaclib, k_Ga64_to_He4_Cu60_reaclib,
        N, Cu61, -1, Cu62, -1, -1, Cu61, Cu62, k_n_Cu61_to_Cu62_reaclib, k_Cu62_to_n_Cu61_reaclib,
        H1, Cu61, -1, Zn62, -1, -1, Cu61, Zn62, k_p_Cu61_to_Zn62_reaclib, k_Zn62_to_p_Cu61_reaclib,
        N, Cu62, -1, Cu63, -1, -1, Cu62, Cu63, k_n_Cu62_to_Cu63_reaclib, k_Cu63_to_n_Cu62_reaclib,
        H1, Cu62, -1, Zn63, -1, -1, Cu62, Zn63, k_p_Cu62_to_Zn63_reaclib, k_Zn63_to_p_Cu62_reaclib,
        N, Cu63, -1, Cu64, -1, -1, Cu63, Cu64, k_n_Cu63_to_Cu64_reaclib, k_Cu64_to_n_Cu63_reaclib,
        H1, Cu63, -1, Zn64, -1, -1, Cu63, Zn64, k_p_Cu63_to_Zn64_reaclib, k_Zn64_to_p_Cu63_reaclib,
        N, Cu64, -1, Cu65, -1, -1, Cu64, Cu65, k_n_Cu64_to_Cu65_reaclib, k_Cu65_to_n_Cu64_reaclib,
        H1, Cu64, -1, Zn65, -1, -1, Cu64, Zn65, k_p_Cu64_to_Zn65_reaclib, k_Zn65_to_p_Cu64_reaclib,
        H1, Cu65, -1, Zn66, -1, -1, Cu65, Zn66, k_p_Cu65_to_Zn66_reaclib, k_Zn66_to_p_Cu65_reaclib,
        N, Zn59, -1, Zn60, -1, -1, Zn59, Zn60, k_n_Zn59_to_Zn60_reaclib, k_Zn60_to_n_Zn59_reaclib,
        He4, Zn59, -1, Ge63, -1, -1, Zn59, Ge63, k_He4_Zn59_to_Ge63_reaclib, k_Ge63_to_He4_Zn59_reaclib,
        N, Zn60, -1, Zn61, -1, -1, Zn60, Zn61, k_n_Zn60_to_Zn61_reaclib, k_Zn61_to_n_Zn60_reaclib,
        He4, Zn60, -1, Ge64, -1, -1, Zn60, Ge64, k_He4_Zn60_to_Ge64_reaclib, k_Ge64_to_He4_Zn60_reaclib,
        N, Zn61, -1, Zn62, -1, -1, Zn61, Zn62, k_n_Zn61_to_Zn62_reaclib, k_Zn62_to_n_Zn61_reaclib,
        H1, Zn61, -1, Ga62, -1, -1, Zn61, Ga62, k_p_Zn61_to_Ga62_reaclib, k_Ga62_to_p_Zn61_reaclib,
        N, Zn62, -1, Zn63, -1, -1, Zn62, Zn63, k_n_Zn62_to_Zn63_reaclib, k_Zn63_to_n_Zn62_reaclib,
        H1, Zn62, -1, Ga63, -1, -1, Zn62, Ga63, k_p_Zn62_to_Ga63_reaclib, k_Ga63_to_p_Zn62_reaclib,
        N, Zn63, -1, Zn64, -1, -1, Zn63, Zn64, k_n_Zn63_to_Zn64_reaclib, k_Zn64_to_n_Zn63_reaclib,
        H1, Zn63, -1, Ga64, -1, -1, Zn63, Ga64, k_p_Zn63_to_Ga64_reaclib, k_Ga64_to_p_Zn63_reaclib,
        N, Zn64, -1, Zn65, -1, -1, Zn64, Zn65, k_n_Zn64_to_Zn65_reaclib, k_Zn65_to_n_Zn64_reaclib,
        N, Zn65, -1, Zn66, -1, -1, Zn65, Zn66, k_n_Zn65_to_Zn66_reaclib, k_Zn66_to_n_Zn65_reaclib,
        N, Ga62, -1, Ga63, -1, -1, Ga62, Ga63, k_n_Ga62_to_Ga63_reaclib, k_Ga63_to_n_Ga62_reaclib,
        H1, Ga62, -1, Ge63, -1, -1, Ga62, Ge63, k_p_Ga62_to_Ge63_reaclib, k_Ge63_to_p_Ga62_reaclib,
        N, Ga63, -1, Ga64, -1, -1, Ga63, Ga64, k_n_Ga63_to_Ga64_reaclib, k_Ga64_to_n_Ga63_reaclib,
        H1, Ga63, -1, Ge64, -1, -1, Ga63, Ge64, k_p_Ga63_to_Ge64_reaclib, k_Ge64_to_p_Ga63_reaclib,
        N, Ge63, -1, Ge64, -1, -1, Ge63, Ge64, k_n_Ge63_to_Ge64_reaclib, k_Ge64_to_n_Ge63_reaclib,
        H2, He3, -1, H1, He4, -1, H2, He3, k_d_He3_to_p_He4_reaclib, k_p_He4_to_d_He3_reaclib,
        H1, Li6, -1, He4, He3, -1, Li6, He3, k_p_Li6_to_He4_He3_reaclib, k_He4_He3_to_p_Li6_reaclib,
        H1, Li7, -1, He4, He4, -1, Li7, -1, k_p_Li7_to_He4_He4_reaclib, k_He4_He4_to_p_Li7_reaclib,
        N, Be7, -1, H1, Li7, -1, Be7, Li7, k_n_Be7_to_p_Li7_reaclib, k_p_Li7_to_n_Be7_reaclib,
        N, Be7, -1, He4, He4, -1, Be7, -1, k_n_Be7_to_He4_He4_reaclib, k_He4_He4_to_n_Be7_reaclib,
        H1, Be9, -1, He4, Li6, -1, Be9, Li6, k_p_Be9_to_He4_Li6_reaclib, k_He4_Li6_to_p_Be9_reaclib,
        He4, Be9, -1, N, C12, -1, Be9, C12, k_He4_Be9_to_n_C12_reaclib, k_n_C12_to_He4_Be9_reaclib,
        N, B10, -1, He4, Li7, -1, B10, Li7, k_n_B10_to_He4_Li7_reaclib, k_He4_Li7_to_n_B10_reaclib,
        H1, B10, -1, He4, Be7, -1, B10, Be7, k_p_B10_to_He4_Be7_reaclib, k_He4_Be7_to_p_B10_reaclib,
        He4, B10, -1, N, N13, -1, B10, N13, k_He4_B10_to_n_N13_reaclib, k_n_N13_to_He4_B10_reaclib,
        He4, B10, -1, H1, C13, -1, B10, C13, k_He4_B10_to_p_C13_reaclib, k_p_C13_to_He4_B10_reaclib,
        He4, B11, -1, N, N14, -1, B11, N14, k_He4_B11_to_n_N14_reaclib, k_n_N14_to_He4_B11_reaclib,
        He4, B11, -1, H1, C14, -1, B11, C14, k_He4_B11_to_p_C14_reaclib, k_p_C14_to_He4_B11_reaclib,
        He4, C13, -1, N, O16, -1, C13, O16, k_He4_C13_to_n_O16_reaclib, k_n_O16_to_He4_C13_reaclib,
        N, N13, -1, H1, C13, -1, N13, C13, k_n_N13_to_p_C13_reaclib, k_p_C13_to_n_N13_reaclib,
        He4, N13, -1, H1, O16, -1, N13, O16, k_He4_N13_to_p_O16_reaclib, k_p_O16_to_He4_N13_reaclib,
        N, N14, -1, H1, C14, -1, N14, C14, k_n_N14_to_p_C14_reaclib, k_p_C14_to_n_N14_reaclib,
        H1, N15, -1, He4, C12, -1, N15, C12, k_p_N15_to_He4_C12_reaclib, k_He4_C12_to_p_N15_reaclib,
        N, O14, -1, H1, N14, -1, O14, N14, k_n_O14_to_p_N14_reaclib, k_p_N14_to_n_O14_reaclib,
        He4, O14, -1, H1, F17, -1, O14, F17, k_He4_O14_to_p_F17_reaclib, k_p_F17_to_He4_O14_reaclib,
        N, O15, -1, H1, N15, -1, O15, N15, k_n_O15_to_p_N15_reaclib, k_p_N15_to_n_O15_reaclib,
        N, O15, -1, He4, C12, -1, O15, C12, k_n_O15_to_He4_C12_reaclib, k_He4_C12_to_n_O15_reaclib,
        N, O17, -1, He4, C14, -1, O17, C14, k_n_O17_to_He4_C14_reaclib, k_He4_C14_to_n_O17_reaclib,
        H1, O17, -1, He4, N14, -1, O17, N14, k_p_O17_to_He4_N14_reaclib, k_He4_N14_to_p_O17_reaclib,
        He4, O17, -1, N, Ne20, -1, O17, Ne20, k_He4_O17_to_n_Ne20_reaclib, k_n_Ne20_to_He4_O17_reaclib,
        H1, O18, -1, He4, N15, -1, O18, N15, k_p_O18_to_He4_N15_reaclib, k_He4_N15_to_p_O18_reaclib,
        N, F17, -1, H1, O17, -1, F17, O17, k_n_F17_to_p_O17_reaclib, k_p_O17_to_n_F17_reaclib,
        N, F17, -1, He4, N14, -1, F17, N14, k_n_F17_to_He4_N14_reaclib, k_He4_N14_to_n_F17_reaclib,
        He4, F17, -1, H1, Ne20, -1, F17, Ne20, k_He4_F17_to_p_Ne20_reaclib, k_p_Ne20_to_He4_F17_reaclib,
        N, F18, -1, H1, O18, -1, F18, O18, k_n_F18_to_p_O18_reaclib, k_p_O18_to_n_F18_reaclib,
        N, F18, -1, He4, N15, -1, F18, N15, k_n_F18_to_He4_N15_reaclib, k_He4_N15_to_n_F18_reaclib,
        H1, F18, -1, He4, O15, -1, F18, O15, k_p_F18_to_He4_O15_reaclib, k_He4_O15_to_p_F18_reaclib,
        He4, F18, -1, H1, Ne21, -1, F18, Ne21, k_He4_F18_to_p_Ne21_reaclib, k_p_Ne21_to_He4_F18_reaclib,
        H1, F19, -1, He4, O16, -1, F19, O16, k_p_F19_to_He4_O16_reaclib, k_He4_O16_to_p_F19_reaclib,
        He4, F19, -1, H1, Ne22, -1, F19, Ne22, k_He4_F19_to_p_Ne22_reaclib, k_p_Ne22_to_He4_F19_reaclib,
        N, Ne18, -1, H1, F18, -1, Ne18, F18, k_n_Ne18_to_p_F18_reaclib, k_p_F18_to_n_Ne18_reaclib,
        N, Ne18, -1, He4, O15, -1, Ne18, O15, k_n_Ne18_to_He4_O15_reaclib, k_He4_O15_to_n_Ne18_reaclib,
        He4, Ne18, -1, H1, Na21, -1, Ne18, Na21, k_He4_Ne18_to_p_Na21_reaclib, k_p_Na21_to_He4_Ne18_reaclib,
        N, Ne19, -1, H1, F19, -1, Ne19, F19, k_n_Ne19_to_p_F19_reaclib, k_p_F19_to_n_Ne19_reaclib,
        N, Ne19, -1, He4, O16, -1, Ne19, O16, k_n_Ne19_to_He4_O16_reaclib, k_He4_O16_to_n_Ne19_reaclib,
        He4, Ne19, -1, H1, Na22, -1, Ne19, Na22, k_He4_Ne19_to_p_Na22_reaclib, k_p_Na22_to_He4_Ne19_reaclib,
        N, Ne21, -1, He4, O18, -1, Ne21, O18, k_n_Ne21_to_He4_O18_reaclib, k_He4_O18_to_n_Ne21_reaclib,
        He4, Ne21, -1, N, Mg24, -1, Ne21, Mg24, k_He4_Ne21_to_n_Mg24_reaclib, k_n_Mg24_to_He4_Ne21_reaclib,
        N, Na21, -1, H1, Ne21, -1, Na21, Ne21, k_n_Na21_to_p_Ne21_reaclib, k_p_Ne21_to_n_Na21_reaclib,
        N, Na21, -1, He4, F18, -1, Na21, F18, k_n_Na21_to_He4_F18_reaclib, k_He4_F18_to_n_Na21_reaclib,
        He4, Na21, -1, H1, Mg24, -1, Na21, Mg24, k_He4_Na21_to_p_Mg24_reaclib, k_p_Mg24_to_He4_Na21_reaclib,
        N, Na22, -1, H1, Ne22, -1, Na22, Ne22, k_n_Na22_to_p_Ne22_reaclib, k_p_Ne22_to_n_Na22_reaclib,
        N, Na22, -1, He4, F19, -1, Na22, F19, k_n_Na22_to_He4_F19_reaclib, k_He4_F19_to_n_Na22_reaclib,
        He4, Na22, -1, H1, Mg25, -1, Na22, Mg25, k_He4_Na22_to_p_Mg25_reaclib, k_p_Mg25_to_He4_Na22_reaclib,
        H1, Na23, -1, He4, Ne20, -1, Na23, Ne20, k_p_Na23_to_He4_Ne20_reaclib, k_He4_Ne20_to_p_Na23_reaclib,
        He4, Na23, -1, H1, Mg26, -1, Na23, Mg26, k_He4_Na23_to_p_Mg26_reaclib, k_p_Mg26_to_He4_Na23_reaclib,
        N, Mg23, -1, H1, Na23, -1, Mg23, Na23, k_n_Mg23_to_p_Na23_reaclib, k_p_Na23_to_n_Mg23_reaclib,
        N, Mg23, -1, He4, Ne20, -1, Mg23, Ne20, k_n_Mg23_to_He4_Ne20_reaclib, k_He4_Ne20_to_n_Mg23_reaclib,
        He4, Mg23, -1, H1, Al26, -1, Mg23, Al26, k_He4_Mg23_to_p_Al26_reaclib, k_p_Al26_to_He4_Mg23_reaclib,
        N, Mg25, -1, He4, Ne22, -1, Mg25, Ne22, k_n_Mg25_to_He4_Ne22_reaclib, k_He4_Ne22_to_n_Mg25_reaclib,
        He4, Mg25, -1, N, Si28, -1, Mg25, Si28, k_He4_Mg25_to_n_Si28_reaclib, k_n_Si28_to_He4_Mg25_reaclib,
        He4, Mg26, -1, N, Si29, -1, Mg26, Si29, k_He4_Mg26_to_n_Si29_reaclib, k_n_Si29_to_He4_Mg26_reaclib,
        N, Al25, -1, H1, Mg25, -1, Al25, Mg25, k_n_Al25_to_p_Mg25_reaclib, k_p_Mg25_to_n_Al25_reaclib,
        N, Al25, -1, He4, Na22, -1, Al25, Na22, k_n_Al25_to_He4_Na22_reaclib, k_He4_Na22_to_n_Al25_reaclib,
        He4, Al25, -1, H1, Si28, -1, Al25, Si28, k_He4_Al25_to_p_Si28_reaclib, k_p_Si28_to_He4_Al25_reaclib,
        N, Al26, -1, H1, Mg26, -1, Al26, Mg26, k_n_Al26_to_p_Mg26_reaclib, k_p_Mg26_to_n_Al26_reaclib,
        N, Al26, -1, He4, Na23, -1, Al26, Na23, k_n_Al26_to_He4_Na23_reaclib, k_He4_Na23_to_n_Al26_reaclib,
        He4, Al26, -1, H1, Si29, -1, Al26, Si29, k_He4_Al26_to_p_Si29_reaclib, k_p_Si29_to_He4_Al26_reaclib,
        H1, Al27, -1, He4, Mg24, -1, Al27, Mg24, k_p_Al27_to_He4_Mg24_reaclib, k_He4_Mg24_to_p_Al27_reaclib,
        He4, Al27, -1, H1, Si30, -1, Al27, Si30, k_He4_Al27_to_p_Si30_reaclib, k_p_Si30_to_He4_Al27_reaclib,
        H1, Si31, -1, N, P31, -1, Si31, P31, k_p_Si31_to_n_P31_reaclib, k_n_P31_to_p_Si31_reaclib,
        He4, Si31, -1, N, S34, -1, Si31, S34, k_He4_Si31_to_n_S34_reaclib, k_n_S34_to_He4_Si31_reaclib,
        N, P29, -1, H1, Si29, -1, P29, Si29, k_n_P29_to_p_Si29_reaclib, k_p_Si29_to_n_P29_reaclib,
        N, P29, -1, He4, Al26, -1, P29, Al26, k_n_P29_to_He4_Al26_reaclib, k_He4_Al26_to_n_P29_reaclib,
        He4, P29, -1, H1, S32, -1, P29, S32, k_He4_P29_to_p_S32_reaclib, k_p_S32_to_He4_P29_reaclib,
        N, P30, -1, H1, Si30, -1, P30, Si30, k_n_P30_to_p_Si30_reaclib, k_p_Si30_to_n_P30_reaclib,
        N, P30, -1, He4, Al27, -1, P30, Al27, k_n_P30_to_He4_Al27_reaclib, k_He4_Al27_to_n_P30_reaclib,
        He4, P30, -1, H1, S33, -1, P30, S33, k_He4_P30_to_p_S33_reaclib, k_p_S33_to_He4_P30_reaclib,
        H1, P31, -1, He4, Si28, -1, P31, Si28, k_p_P31_to_He4_Si28_reaclib, k_He4_Si28_to_p_P31_reaclib,
        He4, P31, -1, H1, S34, -1, P31, S34, k_He4_P31_to_p_S34_reaclib, k_p_S34_to_He4_P31_reaclib,
        N, P32, -1, H1, Si32, -1, P32, Si32, k_n_P32_to_p_Si32_reaclib, k_p_Si32_to_n_P32_reaclib,
        H1, P32, -1, N, S32, -1, P32, S32, k_p_P32_to_n_S32_reaclib, k_n_S32_to_p_P32_reaclib,
        H1, P32, -1, He4, Si29, -1, P32, Si29, k_p_P32_to_He4_Si29_reaclib, k_He4_Si29_to_p_P32_reaclib,
        H1, P33, -1, He4, Si30, -1, P33, Si30, k_p_P33_to_He4_Si30_reaclib, k_He4_Si30_to_p_P33_reaclib,
        N, S32, -1, He4, Si29, -1, S32, Si29, k_n_S32_to_He4_Si29_reaclib, k_He4_Si29_to_n_S32_reaclib,
        N, S33, -1, H1, P33, -1, S33, P33, k_n_S33_to_p_P33_reaclib, k_p_P33_to_n_S33_reaclib,
        N, S33, -1, He4, Si30, -1, S33, Si30, k_n_S33_to_He4_Si30_reaclib, k_He4_Si30_to_n_S33_reaclib,
        N, S35, -1, He4, Si32, -1, S35, Si32, k_n_S35_to_He4_Si32_reaclib, k_He4_Si32_to_n_S35_reaclib,
        H1, S35, -1, He4, P32, -1, S35, P32, k_p_S35_to_He4_P32_reaclib, k_He4_P32_to_p_S35_reaclib,
        He4, S35, -1, N, Ar38, -1, S35, Ar38, k_He4_S35_to_n_Ar38_reaclib, k_n_Ar38_to_He4_S35_reaclib,
        H1, S36, -1, He4, P33, -1, S36, P33, k_p_S36_to_He4_P33_reaclib, k_He4_P33_to_p_S36_reaclib,
        N, Cl33, -1, H1, S33, -1, Cl33, S33, k_n_Cl33_to_p_S33_reaclib, k_p_S33_to_n_Cl33_reaclib,
        N, Cl33, -1, He4, P30, -1, Cl33, P30, k_n_Cl33_to_He4_P30_reaclib, k_He4_P30_to_n_Cl33_reaclib,
        He4, Cl33, -1, H1, Ar36, -1, Cl33, Ar36, k_He4_Cl33_to_p_Ar36_reaclib, k_p_Ar36_to_He4_Cl33_reaclib,
        N, Cl34, -1, H1, S34, -1, Cl34, S34, k_n_Cl34_to_p_S34_reaclib, k_p_S34_to_n_Cl34_reaclib,
        N, Cl34, -1, He4, P31, -1, Cl34, P31, k_n_Cl34_to_He4_P31_reaclib, k_He4_P31_to_n_Cl34_reaclib,
        He4, Cl34, -1, H1, Ar37, -1, Cl34, Ar37, k_He4_Cl34_to_p_Ar37_reaclib, k_p_Ar37_to_He4_Cl34_reaclib,
        N, Cl35, -1, H1, S35, -1, Cl35, S35, k_n_Cl35_to_p_S35_reaclib, k_p_S35_to_n_Cl35_reaclib,
        N, Cl35, -1, He4, P32, -1, Cl35, P32, k_n_Cl35_to_He4_P32_reaclib, k_He4_P32_to_n_Cl35_reaclib,
        H1, Cl35, -1, He4, S32, -1, Cl35, S32, k_p_Cl35_to_He4_S32_reaclib, k_He4_S32_to_p_Cl35_reaclib,
        He4, Cl35, -1, H1, Ar38, -1, Cl35, Ar38, k_He4_Cl35_to_p_Ar38_reaclib, k_p_Ar38_to_He4_Cl35_reaclib,
        N, Cl36, -1, H1, S36, -1, Cl36, S36, k_n_Cl36_to_p_S36_reaclib, k_p_S36_to_n_Cl36_reaclib,
        N, Cl36, -1, He4, P33, -1, Cl36, P33, k_n_Cl36_to_He4_P33_reaclib, k_He4_P33_to_n_Cl36_reaclib,
        H1, Cl36, -1, He4, S33, -1, Cl36, S33, k_p_Cl36_to_He4_S33_reaclib, k_He4_S33_to_p_Cl36_reaclib,
        H1, Cl37, -1, He4, S34, -1, Cl37, S34, k_p_Cl37_to_He4_S34_reaclib, k_He4_S34_to_p_Cl37_reaclib,
        N, Ar36, -1, H1, Cl36, -1, Ar36, Cl36, k_n_Ar36_to_p_Cl36_reaclib, k_p_Cl36_to_n_Ar36_reaclib,
        N, Ar36, -1, He4, S33, -1, Ar36, S33, k_n_Ar36_to_He4_S33_reaclib, k_He4_S33_to_n_Ar36_reaclib,
        N, Ar37, -1, H1, Cl37, -1, Ar37, Cl37, k_n_Ar37_to_p_Cl37_reaclib, k_p_Cl37_to_n_Ar37_reaclib,
        N, Ar37, -1, He4, S34, -1, Ar37, S34, k_n_Ar37_to_He4_S34_reaclib, k_He4_S34_to_n_Ar37_reaclib,
        N, Ar39, -1, He4, S36, -1, Ar39, S36, k_n_Ar39_to_He4_S36_reaclib, k_He4_S36_to_n_Ar39_reaclib,
        H1, Ar39, -1, He4, Cl36, -1, Ar39, Cl36, k_p_Ar39_to_He4_Cl36_reaclib, k_He4_Cl36_to_p_Ar39_reaclib,
        H1, Ar40, -1, He4, Cl37, -1, Ar40, Cl37, k_p_Ar40_to_He4_Cl37_reaclib, k_He4_Cl37_to_p_Ar40_reaclib,
        N, K37, -1, H1, Ar37, -1, K37, Ar37, k_n_K37_to_p_Ar37_reaclib, k_p_Ar37_to_n_K37_reaclib,
        N, K37, -1, He4, Cl34, -1, K37, Cl34, k_n_K37_to_He4_Cl34_reaclib, k_He4_Cl34_to_n_K37_reaclib,
        He4, K37, -1, H1, Ca40, -1, K37, Ca40, k_He4_K37_to_p_Ca40_reaclib, k_p_Ca40_to_He4_K37_reaclib,
        N, K38, -1, H1, Ar38, -1, K38, Ar38, k_n_K38_to_p_Ar38_reaclib, k_p_Ar38_to_n_K38_reaclib,
        N, K38, -1, He4, Cl35, -1, K38, Cl35, k_n_K38_to_He4_Cl35_reaclib, k_He4_Cl35_to_n_K38_reaclib,
        He4, K38, -1, H1, Ca41, -1, K38, Ca41, k_He4_K38_to_p_Ca41_reaclib, k_p_Ca41_to_He4_K38_reaclib,
        N, K39, -1, H1, Ar39, -1, K39, Ar39, k_n_K39_to_p_Ar39_reaclib, k_p_Ar39_to_n_K39_reaclib,
        N, K39, -1, He4, Cl36, -1, K39, Cl36, k_n_K39_to_He4_Cl36_reaclib, k_He4_Cl36_to_n_K39_reaclib,
        H1, K39, -1, He4, Ar36, -1, K39, Ar36, k_p_K39_to_He4_Ar36_reaclib, k_He4_Ar36_to_p_K39_reaclib,
        N, K40, -1, H1, Ar40, -1, K40, Ar40, k_n_K40_to_p_Ar40_reaclib, k_p_Ar40_to_n_K40_reaclib,
        N, K40, -1, He4, Cl37, -1, K40, Cl37, k_n_K40_to_He4_Cl37_reaclib, k_He4_Cl37_to_n_K40_reaclib,
        H1, K40, -1, N, Ca40, -1, K40, Ca40, k_p_K40_to_n_Ca40_reaclib, k_n_Ca40_to_p_K40_reaclib,
        H1, K40, -1, He4, Ar37, -1, K40, Ar37, k_p_K40_to_He4_Ar37_reaclib, k_He4_Ar37_to_p_K40_reaclib,
        He4, K40, -1, H1, Ca43, -1, K40, Ca43, k_He4_K40_to_p_Ca43_reaclib, k_p_Ca43_to_He4_K40_reaclib,
        H1, K41, -1, He4, Ar38, -1, K41, Ar38, k_p_K41_to_He4_Ar38_reaclib, k_He4_Ar38_to_p_K41_reaclib,
        He4, K41, -1, H1, Ca44, -1, K41, Ca44, k_He4_K41_to_p_Ca44_reaclib, k_p_Ca44_to_He4_K41_reaclib,
        N, Ca40, -1, He4, Ar37, -1, Ca40, Ar37, k_n_Ca40_to_He4_Ar37_reaclib, k_He4_Ar37_to_n_Ca40_reaclib,
        N, Ca41, -1, H1, K41, -1, Ca41, K41, k_n_Ca41_to_p_K41_reaclib, k_p_K41_to_n_Ca41_reaclib,
        N, Ca41, -1, He4, Ar38, -1, Ca41, Ar38, k_n_Ca41_to_He4_Ar38_reaclib, k_He4_Ar38_to_n_Ca41_reaclib,
        N, Ca42, -1, He4, Ar39, -1, Ca42, Ar39, k_n_Ca42_to_He4_Ar39_reaclib, k_He4_Ar39_to_n_Ca42_reaclib,
        H1, Ca42, -1, He4, K39, -1, Ca42, K39, k_p_Ca42_to_He4_K39_reaclib, k_He4_K39_to_p_Ca42_reaclib,
        N, Ca43, -1, He4, Ar40, -1, Ca43, Ar40, k_n_Ca43_to_He4_Ar40_reaclib, k_He4_Ar40_to_n_Ca43_reaclib,
        He4, Ca43, -1, N, Ti46, -1, Ca43, Ti46, k_He4_Ca43_to_n_Ti46_reaclib, k_n_Ti46_to_He4_Ca43_reaclib,
        He4, Ca45, -1, N, Ti48, -1, Ca45, Ti48, k_He4_Ca45_to_n_Ti48_reaclib, k_n_Ti48_to_He4_Ca45_reaclib,
        H1, Ca47, -1, N, Sc47, -1, Ca47, Sc47, k_p_Ca47_to_n_Sc47_reaclib, k_n_Sc47_to_p_Ca47_reaclib,
        He4, Ca47, -1, N, Ti50, -1, Ca47, Ti50, k_He4_Ca47_to_n_Ti50_reaclib, k_n_Ti50_to_He4_Ca47_reaclib,
        N, Sc43, -1, H1, Ca43, -1, Sc43, Ca43, k_n_Sc43_to_p_Ca43_reaclib, k_p_Ca43_to_n_Sc43_reaclib,
        N, Sc43, -1, He4, K40, -1, Sc43, K40, k_n_Sc43_to_He4_K40_reaclib, k_He4_K40_to_n_Sc43_reaclib,
        H1, Sc43, -1, He4, Ca40, -1, Sc43, Ca40, k_p_Sc43_to_He4_Ca40_reaclib, k_He4_Ca40_to_p_Sc43_reaclib,
        He4, Sc43, -1, H1, Ti46, -1, Sc43, Ti46, k_He4_Sc43_to_p_Ti46_reaclib, k_p_Ti46_to_He4_Sc43_reaclib,
        N, Sc44, -1, H1, Ca44, -1, Sc44, Ca44, k_n_Sc44_to_p_Ca44_reaclib, k_p_Ca44_to_n_Sc44_reaclib,
        N, Sc44, -1, He4, K41, -1, Sc44, K41, k_n_Sc44_to_He4_K41_reaclib, k_He4_K41_to_n_Sc44_reaclib,
        H1, Sc44, -1, He4, Ca41, -1, Sc44, Ca41, k_p_Sc44_to_He4_Ca41_reaclib, k_He4_Ca41_to_p_Sc44_reaclib,
        He4, Sc44, -1, H1, Ti47, -1, Sc44, Ti47, k_He4_Sc44_to_p_Ti47_reaclib, k_p_Ti47_to_He4_Sc44_reaclib,
        N, Sc45, -1, H1, Ca45, -1, Sc45, Ca45, k_n_Sc45_to_p_Ca45_reaclib, k_p_Ca45_to_n_Sc45_reaclib,
        H1, Sc45, -1, He4, Ca42, -1, Sc45, Ca42, k_p_Sc45_to_He4_Ca42_reaclib, k_He4_Ca42_to_p_Sc45_reaclib,
        He4, Sc45, -1, H1, Ti48, -1, Sc45, Ti48, k_He4_Sc45_to_p_Ti48_reaclib, k_p_Ti48_to_He4_Sc45_reaclib,
        N, Sc46, -1, H1, Ca46, -1, Sc46, Ca46, k_n_Sc46_to_p_Ca46_reaclib, k_p_Ca46_to_n_Sc46_reaclib,
        H1, Sc46, -1, N, Ti46, -1, Sc46, Ti46, k_p_Sc46_to_n_Ti46_reaclib, k_n_Ti46_to_p_Sc46_reaclib,
        H1, Sc46, -1, He4, Ca43, -1, Sc46, Ca43, k_p_Sc46_to_He4_Ca43_reaclib, k_He4_Ca43_to_p_Sc46_reaclib,
        He4, Sc46, -1, N, V49, -1, Sc46, V49, k_He4_Sc46_to_n_V49_reaclib, k_n_V49_to_He4_Sc46_reaclib,
        He4, Sc46, -1, H1, Ti49, -1, Sc46, Ti49, k_He4_Sc46_to_p_Ti49_reaclib, k_p_Ti49_to_He4_Sc46_reaclib,
        H1, Sc47, -1, He4, Ca44, -1, Sc47, Ca44, k_p_Sc47_to_He4_Ca44_reaclib, k_He4_Ca44_to_p_Sc47_reaclib,
        He4, Sc47, -1, H1, Ti50, -1, Sc47, Ti50, k_He4_Sc47_to_p_Ti50_reaclib, k_p_Ti50_to_He4_Sc47_reaclib,
        N, Sc48, -1, H1, Ca48, -1, Sc48, Ca48, k_n_Sc48_to_p_Ca48_reaclib, k_p_Ca48_to_n_Sc48_reaclib,
        H1, Sc48, -1, N, Ti48, -1, Sc48, Ti48, k_p_Sc48_to_n_Ti48_reaclib, k_n_Ti48_to_p_Sc48_reaclib,
        H1, Sc48, -1, He4, Ca45, -1, Sc48, Ca45, k_p_Sc48_to_He4_Ca45_reaclib, k_He4_Ca45_to_p_Sc48_reaclib,
        He4, Sc48, -1, N, V51, -1, Sc48, V51, k_He4_Sc48_to_n_V51_reaclib, k_n_V51_to_He4_Sc48_reaclib,
        He4, Sc48, -1, H1, Ti51, -1, Sc48, Ti51, k_He4_Sc48_to_p_Ti51_reaclib, k_p_Ti51_to_He4_Sc48_reaclib,
        H1, Sc49, -1, N, Ti49, -1, Sc49, Ti49, k_p_Sc49_to_n_Ti49_reaclib, k_n_Ti49_to_p_Sc49_reaclib,
        H1, Sc49, -1, He4, Ca46, -1, Sc49, Ca46, k_p_Sc49_to_He4_Ca46_reaclib, k_He4_Ca46_to_p_Sc49_reaclib,
        N, Ti44, -1, H1, Sc44, -1, Ti44, Sc44, k_n_Ti44_to_p_Sc44_reaclib, k_p_Sc44_to_n_Ti44_reaclib,
        N, Ti44, -1, He4, Ca41, -1, Ti44, Ca41, k_n_Ti44_to_He4_Ca41_reaclib, k_He4_Ca41_to_n_Ti44_reaclib,
        N, Ti45, -1, H1, Sc45, -1, Ti45, Sc45, k_n_Ti45_to_p_Sc45_reaclib, k_p_Sc45_to_n_Ti45_reaclib,
        N, Ti45, -1, He4, Ca42, -1, Ti45, Ca42, k_n_Ti45_to_He4_Ca42_reaclib, k_He4_Ca42_to_n_Ti45_reaclib,
        He4, Ti45, -1, H1, V48, -1, Ti45, V48, k_He4_Ti45_to_p_V48_reaclib, k_p_V48_to_He4_Ti45_reaclib,
        N, Ti47, -1, H1, Sc47, -1, Ti47, Sc47, k_n_Ti47_to_p_Sc47_reaclib, k_p_Sc47_to_n_Ti47_reaclib,
        N, Ti47, -1, He4, Ca44, -1, Ti47, Ca44, k_n_Ti47_to_He4_Ca44_reaclib, k_He4_Ca44_to_n_Ti47_reaclib,
        N, Ti49, -1, He4, Ca46, -1, Ti49, Ca46, k_n_Ti49_to_He4_Ca46_reaclib, k_He4_Ca46_to_n_Ti49_reaclib,
        He4, Ti49, -1, N, Cr52, -1, Ti49, Cr52, k_He4_Ti49_to_n_Cr52_reaclib, k_n_Cr52_to_He4_Ti49_reaclib,
        N, Ti51, -1, He4, Ca48, -1, Ti51, Ca48, k_n_Ti51_to_He4_Ca48_reaclib, k_He4_Ca48_to_n_Ti51_reaclib,
        H1, Ti51, -1, N, V51, -1, Ti51, V51, k_p_Ti51_to_n_V51_reaclib, k_n_V51_to_p_Ti51_reaclib,
        He4, Ti51, -1, N, Cr54, -1, Ti51, Cr54, k_He4_Ti51_to_n_Cr54_reaclib, k_n_Cr54_to_He4_Ti51_reaclib,
        N, V46, -1, H1, Ti46, -1, V46, Ti46, k_n_V46_to_p_Ti46_reaclib, k_p_Ti46_to_n_V46_reaclib,
        N, V46, -1, He4, Sc43, -1, V46, Sc43, k_n_V46_to_He4_Sc43_reaclib, k_He4_Sc43_to_n_V46_reaclib,
        He4, V46, -1, H1, Cr49, -1, V46, Cr49, k_He4_V46_to_p_Cr49_reaclib, k_p_Cr49_to_He4_V46_reaclib,
        N, V47, -1, H1, Ti47, -1, V47, Ti47, k_n_V47_to_p_Ti47_reaclib, k_p_Ti47_to_n_V47_reaclib,
        N, V47, -1, He4, Sc44, -1, V47, Sc44, k_n_V47_to_He4_Sc44_reaclib, k_He4_Sc44_to_n_V47_reaclib,
        H1, V47, -1, He4, Ti44, -1, V47, Ti44, k_p_V47_to_He4_Ti44_reaclib, k_He4_Ti44_to_p_V47_reaclib,
        He4, V47, -1, H1, Cr50, -1, V47, Cr50, k_He4_V47_to_p_Cr50_reaclib, k_p_Cr50_to_He4_V47_reaclib,
        N, V48, -1, H1, Ti48, -1, V48, Ti48, k_n_V48_to_p_Ti48_reaclib, k_p_Ti48_to_n_V48_reaclib,
        N, V48, -1, He4, Sc45, -1, V48, Sc45, k_n_V48_to_He4_Sc45_reaclib, k_He4_Sc45_to_n_V48_reaclib,
        He4, V48, -1, H1, Cr51, -1, V48, Cr51, k_He4_V48_to_p_Cr51_reaclib, k_p_Cr51_to_He4_V48_reaclib,
        N, V49, -1, H1, Ti49, -1, V49, Ti49, k_n_V49_to_p_Ti49_reaclib, k_p_Ti49_to_n_V49_reaclib,
        H1, V49, -1, He4, Ti46, -1, V49, Ti46, k_p_V49_to_He4_Ti46_reaclib, k_He4_Ti46_to_p_V49_reaclib,
        He4, V49, -1, H1, Cr52, -1, V49, Cr52, k_He4_V49_to_p_Cr52_reaclib, k_p_Cr52_to_He4_V49_reaclib,
        N, V50, -1, H1, Ti50, -1, V50, Ti50, k_n_V50_to_p_Ti50_reaclib, k_p_Ti50_to_n_V50_reaclib,
        N, V50, -1, He4, Sc47, -1, V50, Sc47, k_n_V50_to_He4_Sc47_reaclib, k_He4_Sc47_to_n_V50_reaclib,
        H1, V50, -1, N, Cr50, -1, V50, Cr50, k_p_V50_to_n_Cr50_reaclib, k_n_Cr50_to_p_V50_reaclib,
        H1, V50, -1, He4, Ti47, -1, V50, Ti47, k_p_V50_to_He4_Ti47_reaclib, k_He4_Ti47_to_p_V50_reaclib,
        He4, V50, -1, H1, Cr53, -1, V50, Cr53, k_He4_V50_to_p_Cr53_reaclib, k_p_Cr53_to_He4_V50_reaclib,
        H1, V51, -1, He4, Ti48, -1, V51, Ti48, k_p_V51_to_He4_Ti48_reaclib, k_He4_Ti48_to_p_V51_reaclib,
        N, V52, -1, He4, Sc49, -1, V52, Sc49, k_n_V52_to_He4_Sc49_reaclib, k_He4_Sc49_to_n_V52_reaclib,
        H1, V52, -1, N, Cr52, -1, V52, Cr52, k_p_V52_to_n_Cr52_reaclib, k_n_Cr52_to_p_V52_reaclib,
        H1, V52, -1, He4, Ti49, -1, V52, Ti49, k_p_V52_to_He4_Ti49_reaclib, k_He4_Ti49_to_p_V52_reaclib,
        He4, V52, -1, N, Mn55, -1, V52, Mn55, k_He4_V52_to_n_Mn55_reaclib, k_n_Mn55_to_He4_V52_reaclib,
        N, Cr48, -1, H1, V48, -1, Cr48, V48, k_n_Cr48_to_p_V48_reaclib, k_p_V48_to_n_Cr48_reaclib,
        N, Cr48, -1, He4, Ti45, -1, Cr48, Ti45, k_n_Cr48_to_He4_Ti45_reaclib, k_He4_Ti45_to_n_Cr48_reaclib,
        He4, Cr48, -1, H1, Mn51, -1, Cr48, Mn51, k_He4_Cr48_to_p_Mn51_reaclib, k_p_Mn51_to_He4_Cr48_reaclib,
        N, Cr49, -1, H1, V49, -1, Cr49, V49, k_n_Cr49_to_p_V49_reaclib, k_p_V49_to_n_Cr49_reaclib,
        N, Cr49, -1, He4, Ti46, -1, Cr49, Ti46, k_n_Cr49_to_He4_Ti46_reaclib, k_He4_Ti46_to_n_Cr49_reaclib,
        He4, Cr49, -1, H1, Mn52, -1, Cr49, Mn52, k_He4_Cr49_to_p_Mn52_reaclib, k_p_Mn52_to_He4_Cr49_reaclib,
        N, Cr50, -1, He4, Ti47, -1, Cr50, Ti47, k_n_Cr50_to_He4_Ti47_reaclib, k_He4_Ti47_to_n_Cr50_reaclib,
        N, Cr51, -1, H1, V51, -1, Cr51, V51, k_n_Cr51_to_p_V51_reaclib, k_p_V51_to_n_Cr51_reaclib,
        N, Cr51, -1, He4, Ti48, -1, Cr51, Ti48, k_n_Cr51_to_He4_Ti48_reaclib, k_He4_Ti48_to_n_Cr51_reaclib,
        N, Cr53, -1, He4, Ti50, -1, Cr53, Ti50, k_n_Cr53_to_He4_Ti50_reaclib, k_He4_Ti50_to_n_Cr53_reaclib,
        H1, Cr54, -1, He4, V51, -1, Cr54, V51, k_p_Cr54_to_He4_V51_reaclib, k_He4_V51_to_p_Cr54_reaclib,
        N, Mn50, -1, H1, Cr50, -1, Mn50, Cr50, k_n_Mn50_to_p_Cr50_reaclib, k_p_Cr50_to_n_Mn50_reaclib,
        N, Mn50, -1, He4, V47, -1, Mn50, V47, k_n_Mn50_to_He4_V47_reaclib, k_He4_V47_to_n_Mn50_reaclib,
        He4, Mn50, -1, H1, Fe53, -1, Mn50, Fe53, k_He4_Mn50_to_p_Fe53_reaclib, k_p_Fe53_to_He4_Mn50_reaclib,
        N, Mn51, -1, H1, Cr51, -1, Mn51, Cr51, k_n_Mn51_to_p_Cr51_reaclib, k_p_Cr51_to_n_Mn51_reaclib,
        N, Mn51, -1, He4, V48, -1, Mn51, V48, k_n_Mn51_to_He4_V48_reaclib, k_He4_V48_to_n_Mn51_reaclib,
        He4, Mn51, -1, H1, Fe54, -1, Mn51, Fe54, k_He4_Mn51_to_p_Fe54_reaclib, k_p_Fe54_to_He4_Mn51_reaclib,
        N, Mn52, -1, H1, Cr52, -1, Mn52, Cr52, k_n_Mn52_to_p_Cr52_reaclib, k_p_Cr52_to_n_Mn52_reaclib,
        N, Mn52, -1, He4, V49, -1, Mn52, V49, k_n_Mn52_to_He4_V49_reaclib, k_He4_V49_to_n_Mn52_reaclib,
        He4, Mn52, -1, H1, Fe55, -1, Mn52, Fe55, k_He4_Mn52_to_p_Fe55_reaclib, k_p_Fe55_to_He4_Mn52_reaclib,
        N, Mn53, -1, H1, Cr53, -1, Mn53, Cr53, k_n_Mn53_to_p_Cr53_reaclib, k_p_Cr53_to_n_Mn53_reaclib,
        N, Mn53, -1, He4, V50, -1, Mn53, V50, k_n_Mn53_to_He4_V50_reaclib, k_He4_V50_to_n_Mn53_reaclib,
        H1, Mn53, -1, He4, Cr50, -1, Mn53, Cr50, k_p_Mn53_to_He4_Cr50_reaclib, k_He4_Cr50_to_p_Mn53_reaclib,
        He4, Mn53, -1, H1, Fe56, -1, Mn53, Fe56, k_He4_Mn53_to_p_Fe56_reaclib, k_p_Fe56_to_He4_Mn53_reaclib,
        N, Mn54, -1, H1, Cr54, -1, Mn54, Cr54, k_n_Mn54_to_p_Cr54_reaclib, k_p_Cr54_to_n_Mn54_reaclib,
        N, Mn54, -1, He4, V51, -1, Mn54, V51, k_n_Mn54_to_He4_V51_reaclib, k_He4_V51_to_n_Mn54_reaclib,
        H1, Mn54, -1, He4, Cr51, -1, Mn54, Cr51, k_p_Mn54_to_He4_Cr51_reaclib, k_He4_Cr51_to_p_Mn54_reaclib,
        H1, Mn55, -1, He4, Cr52, -1, Mn55, Cr52, k_p_Mn55_to_He4_Cr52_reaclib, k_He4_Cr52_to_p_Mn55_reaclib,
        N, Fe52, -1, H1, Mn52, -1, Fe52, Mn52, k_n_Fe52_to_p_Mn52_reaclib, k_p_Mn52_to_n_Fe52_reaclib,
        N, Fe52, -1, He4, Cr49, -1, Fe52, Cr49, k_n_Fe52_to_He4_Cr49_reaclib, k_He4_Cr49_to_n_Fe52_reaclib,
        He4, Fe52, -1, H1, Co55, -1, Fe52, Co55, k_He4_Fe52_to_p_Co55_reaclib, k_p_Co55_to_He4_Fe52_reaclib,
        N, Fe53, -1, H1, Mn53, -1, Fe53, Mn53, k_n_Fe53_to_p_Mn53_reaclib, k_p_Mn53_to_n_Fe53_reaclib,
        N, Fe53, -1, He4, Cr50, -1, Fe53, Cr50, k_n_Fe53_to_He4_Cr50_reaclib, k_He4_Cr50_to_n_Fe53_reaclib,
        He4, Fe53, -1, H1, Co56, -1, Fe53, Co56, k_He4_Fe53_to_p_Co56_reaclib, k_p_Co56_to_He4_Fe53_reaclib,
        N, Fe54, -1, H1, Mn54, -1, Fe54, Mn54, k_n_Fe54_to_p_Mn54_reaclib, k_p_Mn54_to_n_Fe54_reaclib,
        N, Fe54, -1, He4, Cr51, -1, Fe54, Cr51, k_n_Fe54_to_He4_Cr51_reaclib, k_He4_Cr51_to_n_Fe54_reaclib,
        N, Fe55, -1, H1, Mn55, -1, Fe55, Mn55, k_n_Fe55_to_p_Mn55_reaclib, k_p_Mn55_to_n_Fe55_reaclib,
        N, Fe55, -1, He4, Cr52, -1, Fe55, Cr52, k_n_Fe55_to_He4_Cr52_reaclib, k_He4_Cr52_to_n_Fe55_reaclib,
        N, Fe56, -1, He4, Cr53, -1, Fe56, Cr53, k_n_Fe56_to_He4_Cr53_reaclib, k_He4_Cr53_to_n_Fe56_reaclib,
        N, Fe57, -1, He4, Cr54, -1, Fe57, Cr54, k_n_Fe57_to_He4_Cr54_reaclib, k_He4_Cr54_to_n_Fe57_reaclib,
        H1, Fe57, -1, He4, Mn54, -1, Fe57, Mn54, k_p_Fe57_to_He4_Mn54_reaclib, k_He4_Mn54_to_p_Fe57_reaclib,
        H1, Fe58, -1, He4, Mn55, -1, Fe58, Mn55, k_p_Fe58_to_He4_Mn55_reaclib, k_He4_Mn55_to_p_Fe58_reaclib,
        N, Co53, -1, H1, Fe53, -1, Co53, Fe53, k_n_Co53_to_p_Fe53_reaclib, k_p_Fe53_to_n_Co53_reaclib,
        N, Co53, -1, He4, Mn50, -1, Co53, Mn50, k_n_Co53_to_He4_Mn50_reaclib, k_He4_Mn50_to_n_Co53_reaclib,
        He4, Co53, -1, H1, Ni56, -1, Co53, Ni56, k_He4_Co53_to_p_Ni56_reaclib, k_p_Ni56_to_He4_Co53_reaclib,
        N, Co54, -1, H1, Fe54, -1, Co54, Fe54, k_n_Co54_to_p_Fe54_reaclib, k_p_Fe54_to_n_Co54_reaclib,
        N, Co54, -1, He4, Mn51, -1, Co54, Mn51, k_n_Co54_to_He4_Mn51_reaclib, k_He4_Mn51_to_n_Co54_reaclib,
        He4, Co54, -1, H1, Ni57, -1, Co54, Ni57, k_He4_Co54_to_p_Ni57_reaclib, k_p_Ni57_to_He4_Co54_reaclib,
        N, Co55, -1, H1, Fe55, -1, Co55, Fe55, k_n_Co55_to_p_Fe55_reaclib, k_p_Fe55_to_n_Co55_reaclib,
        N, Co55, -1, He4, Mn52, -1, Co55, Mn52, k_n_Co55_to_He4_Mn52_reaclib, k_He4_Mn52_to_n_Co55_reaclib,
        He4, Co55, -1, H1, Ni58, -1, Co55, Ni58, k_He4_Co55_to_p_Ni58_reaclib, k_p_Ni58_to_He4_Co55_reaclib,
        N, Co56, -1, H1, Fe56, -1, Co56, Fe56, k_n_Co56_to_p_Fe56_reaclib, k_p_Fe56_to_n_Co56_reaclib,
        N, Co56, -1, He4, Mn53, -1, Co56, Mn53, k_n_Co56_to_He4_Mn53_reaclib, k_He4_Mn53_to_n_Co56_reaclib,
        He4, Co56, -1, H1, Ni59, -1, Co56, Ni59, k_He4_Co56_to_p_Ni59_reaclib, k_p_Ni59_to_He4_Co56_reaclib,
        N, Co57, -1, H1, Fe57, -1, Co57, Fe57, k_n_Co57_to_p_Fe57_reaclib, k_p_Fe57_to_n_Co57_reaclib,
        N, Co57, -1, He4, Mn54, -1, Co57, Mn54, k_n_Co57_to_He4_Mn54_reaclib, k_He4_Mn54_to_n_Co57_reaclib,
        H1, Co57, -1, He4, Fe54, -1, Co57, Fe54, k_p_Co57_to_He4_Fe54_reaclib, k_He4_Fe54_to_p_Co57_reaclib,
        He4, Co57, -1, H1, Ni60, -1, Co57, Ni60, k_He4_Co57_to_p_Ni60_reaclib, k_p_Ni60_to_He4_Co57_reaclib,
        N, Co58, -1, H1, Fe58, -1, Co58, Fe58, k_n_Co58_to_p_Fe58_reaclib, k_p_Fe58_to_n_Co58_reaclib,
        N, Co58, -1, He4, Mn55, -1, Co58, Mn55, k_n_Co58_to_He4_Mn55_reaclib, k_He4_Mn55_to_n_Co58_reaclib,
        H1, Co58, -1, He4, Fe55, -1, Co58, Fe55, k_p_Co58_to_He4_Fe55_reaclib, k_He4_Fe55_to_p_Co58_reaclib,
        H1, Co59, -1, He4, Fe56, -1, Co59, Fe56, k_p_Co59_to_He4_Fe56_reaclib, k_He4_Fe56_to_p_Co59_reaclib,
        N, Ni56, -1, H1, Co56, -1, Ni56, Co56, k_n_Ni56_to_p_Co56_reaclib, k_p_Co56_to_n_Ni56_reaclib,
        N, Ni56, -1, He4, Fe53, -1, Ni56, Fe53, k_n_Ni56_to_He4_Fe53_reaclib, k_He4_Fe53_to_n_Ni56_reaclib,
        N, Ni57, -1, H1, Co57, -1, Ni57, Co57, k_n_Ni57_to_p_Co57_reaclib, k_p_Co57_to_n_Ni57_reaclib,
        N, Ni57, -1, He4, Fe54, -1, Ni57, Fe54, k_n_Ni57_to_He4_Fe54_reaclib, k_He4_Fe54_to_n_Ni57_reaclib,
        N, Ni58, -1, H1, Co58, -1, Ni58, Co58, k_n_Ni58_to_p_Co58_reaclib, k_p_Co58_to_n_Ni58_reaclib,
        N, Ni58, -1, He4, Fe55, -1, Ni58, Fe55, k_n_Ni58_to_He4_Fe55_reaclib, k_He4_Fe55_to_n_Ni58_reaclib,
        N, Ni59, -1, H1, Co59, -1, Ni59, Co59, k_n_Ni59_to_p_Co59_reaclib, k_p_Co59_to_n_Ni59_reaclib,
        N, Ni59, -1, He4, Fe56, -1, Ni59, Fe56, k_n_Ni59_to_He4_Fe56_reaclib, k_He4_Fe56_to_n_Ni59_reaclib,
        N, Ni60, -1, He4, Fe57, -1, Ni60, Fe57, k_n_Ni60_to_He4_Fe57_reaclib, k_He4_Fe57_to_n_Ni60_reaclib,
        N, Ni61, -1, He4, Fe58, -1, Ni61, Fe58, k_n_Ni61_to_He4_Fe58_reaclib, k_He4_Fe58_to_n_Ni61_reaclib,
        H1, Ni61, -1, He4, Co58, -1, Ni61, Co58, k_p_Ni61_to_He4_Co58_reaclib, k_He4_Co58_to_p_Ni61_reaclib,
        H1, Ni62, -1, He4, Co59, -1, Ni62, Co59, k_p_Ni62_to_He4_Co59_reaclib, k_He4_Co59_to_p_Ni62_reaclib,
        N, Cu57, -1, H1, Ni57, -1, Cu57, Ni57, k_n_Cu57_to_p_Ni57_reaclib, k_p_Ni57_to_n_Cu57_reaclib,
        N, Cu57, -1, He4, Co54, -1, Cu57, Co54, k_n_Cu57_to_He4_Co54_reaclib, k_He4_Co54_to_n_Cu57_reaclib,
        He4, Cu57, -1, H1, Zn60, -1, Cu57, Zn60, k_He4_Cu57_to_p_Zn60_reaclib, k_p_Zn60_to_He4_Cu57_reaclib,
        N, Cu58, -1, H1, Ni58, -1, Cu58, Ni58, k_n_Cu58_to_p_Ni58_reaclib, k_p_Ni58_to_n_Cu58_reaclib,
        N, Cu58, -1, He4, Co55, -1, Cu58, Co55, k_n_Cu58_to_He4_Co55_reaclib, k_He4_Co55_to_n_Cu58_reaclib,
        N, Cu59, -1, H1, Ni59, -1, Cu59, Ni59, k_n_Cu59_to_p_Ni59_reaclib, k_p_Ni59_to_n_Cu59_reaclib,
        N, Cu59, -1, He4, Co56, -1, Cu59, Co56, k_n_Cu59_to_He4_Co56_reaclib, k_He4_Co56_to_n_Cu59_reaclib,
        H1, Cu59, -1, He4, Ni56, -1, Cu59, Ni56, k_p_Cu59_to_He4_Ni56_reaclib, k_He4_Ni56_to_p_Cu59_reaclib,
        N, Cu60, -1, H1, Ni60, -1, Cu60, Ni60, k_n_Cu60_to_p_Ni60_reaclib, k_p_Ni60_to_n_Cu60_reaclib,
        N, Cu60, -1, He4, Co57, -1, Cu60, Co57, k_n_Cu60_to_He4_Co57_reaclib, k_He4_Co57_to_n_Cu60_reaclib,
        H1, Cu60, -1, He4, Ni57, -1, Cu60, Ni57, k_p_Cu60_to_He4_Ni57_reaclib, k_He4_Ni57_to_p_Cu60_reaclib,
        N, Cu61, -1, H1, Ni61, -1, Cu61, Ni61, k_n_Cu61_to_p_Ni61_reaclib, k_p_Ni61_to_n_Cu61_reaclib,
        N, Cu61, -1, He4, Co58, -1, Cu61, Co58, k_n_Cu61_to_He4_Co58_reaclib, k_He4_Co58_to_n_Cu61_reaclib,
        H1, Cu61, -1, He4, Ni58, -1, Cu61, Ni58, k_p_Cu61_to_He4_Ni58_reaclib, k_He4_Ni58_to_p_Cu61_reaclib,
        N, Cu62, -1, H1, Ni62, -1, Cu62, Ni62, k_n_Cu62_to_p_Ni62_reaclib, k_p_Ni62_to_n_Cu62_reaclib,
        N, Cu62, -1, He4, Co59, -1, Cu62, Co59, k_n_Cu62_to_He4_Co59_reaclib, k_He4_Co59_to_n_Cu62_reaclib,
        H1, Cu62, -1, He4, Ni59, -1, Cu62, Ni59, k_p_Cu62_to_He4_Ni59_reaclib, k_He4_Ni59_to_p_Cu62_reaclib,
        N, Cu63, -1, H1, Ni63, -1, Cu63, Ni63, k_n_Cu63_to_p_Ni63_reaclib, k_p_Ni63_to_n_Cu63_reaclib,
        H1, Cu63, -1, He4, Ni60, -1, Cu63, Ni60, k_p_Cu63_to_He4_Ni60_reaclib, k_He4_Ni60_to_p_Cu63_reaclib,
        N, Cu64, -1, H1, Ni64, -1, Cu64, Ni64, k_n_Cu64_to_p_Ni64_reaclib, k_p_Ni64_to_n_Cu64_reaclib,
        H1, Cu64, -1, He4, Ni61, -1, Cu64, Ni61, k_p_Cu64_to_He4_Ni61_reaclib, k_He4_Ni61_to_p_Cu64_reaclib,
        H1, Cu65, -1, He4, Ni62, -1, Cu65, Ni62, k_p_Cu65_to_He4_Ni62_reaclib, k_He4_Ni62_to_p_Cu65_reaclib,
        N, Zn59, -1, H1, Cu59, -1, Zn59, Cu59, k_n_Zn59_to_p_Cu59_reaclib, k_p_Cu59_to_n_Zn59_reaclib,
        N, Zn59, -1, He4, Ni56, -1, Zn59, Ni56, k_n_Zn59_to_He4_Ni56_reaclib, k_He4_Ni56_to_n_Zn59_reaclib,
        N, Zn60, -1, H1, Cu60, -1, Zn60, Cu60, k_n_Zn60_to_p_Cu60_reaclib, k_p_Cu60_to_n_Zn60_reaclib,
        N, Zn60, -1, He4, Ni57, -1, Zn60, Ni57, k_n_Zn60_to_He4_Ni57_reaclib, k_He4_Ni57_to_n_Zn60_reaclib,
        N, Zn61, -1, H1, Cu61, -1, Zn61, Cu61, k_n_Zn61_to_p_Cu61_reaclib, k_p_Cu61_to_n_Zn61_reaclib,
        N, Zn61, -1, He4, Ni58, -1, Zn61, Ni58, k_n_Zn61_to_He4_Ni58_reaclib, k_He4_Ni58_to_n_Zn61_reaclib,
        H1, Zn61, -1, He4, Cu58, -1, Zn61, Cu58, k_p_Zn61_to_He4_Cu58_reaclib, k_He4_Cu58_to_p_Zn61_reaclib,
        N, Zn62, -1, H1, Cu62, -1, Zn62, Cu62, k_n_Zn62_to_p_Cu62_reaclib, k_p_Cu62_to_n_Zn62_reaclib,
        N, Zn62, -1, He4, Ni59, -1, Zn62, Ni59, k_n_Zn62_to_He4_Ni59_reaclib, k_He4_Ni59_to_n_Zn62_reaclib,
        H1, Zn62, -1, He4, Cu59, -1, Zn62, Cu59, k_p_Zn62_to_He4_Cu59_reaclib, k_He4_Cu59_to_p_Zn62_reaclib,
        N, Zn63, -1, H1, Cu63, -1, Zn63, Cu63, k_n_Zn63_to_p_Cu63_reaclib, k_p_Cu63_to_n_Zn63_reaclib,
        N, Zn63, -1, He4, Ni60, -1, Zn63, Ni60, k_n_Zn63_to_He4_Ni60_reaclib, k_He4_Ni60_to_n_Zn63_reaclib,
        H1, Zn63, -1, He4, Cu60, -1, Zn63, Cu60, k_p_Zn63_to_He4_Cu60_reaclib, k_He4_Cu60_to_p_Zn63_reaclib,
        N, Zn64, -1, H1, Cu64, -1, Zn64, Cu64, k_n_Zn64_to_p_Cu64_reaclib, k_p_Cu64_to_n_Zn64_reaclib,
        N, Zn64, -1, He4, Ni61, -1, Zn64, Ni61, k_n_Zn64_to_He4_Ni61_reaclib, k_He4_Ni61_to_n_Zn64_reaclib,
        H1, Zn64, -1, He4, Cu61, -1, Zn64, Cu61, k_p_Zn64_to_He4_Cu61_reaclib, k_He4_Cu61_to_p_Zn64_reaclib,
        N, Zn65, -1, H1, Cu65, -1, Zn65, Cu65, k_n_Zn65_to_p_Cu65_reaclib, k_p_Cu65_to_n_Zn65_reaclib,
        N, Zn65, -1, He4, Ni62, -1, Zn65, Ni62, k_n_Zn65_to_He4_Ni62_reaclib, k_He4_Ni62_to_n_Zn65_reaclib,
        H1, Zn65, -1, He4, Cu62, -1, Zn65, Cu62, k_p_Zn65_to_He4_Cu62_reaclib, k_He4_Cu62_to_p_Zn65_reaclib,
        N, Zn66, -1, He4, Ni63, -1, Zn66, Ni63, k_n_Zn66_to_He4_Ni63_reaclib, k_He4_Ni63_to_n_Zn66_reaclib,
        H1, Zn66, -1, He4, Cu63, -1, Zn66, Cu63, k_p_Zn66_to_He4_Cu63_reaclib, k_He4_Cu63_to_p_Zn66_reaclib,
        N, Ga62, -1, H1, Zn62, -1, Ga62, Zn62, k_n_Ga62_to_p_Zn62_reaclib, k_p_Zn62_to_n_Ga62_reaclib,
        N, Ga62, -1, He4, Cu59, -1, Ga62, Cu59, k_n_Ga62_to_He4_Cu59_reaclib, k_He4_Cu59_to_n_Ga62_reaclib,
        H1, Ga62, -1, He4, Zn59, -1, Ga62, Zn59, k_p_Ga62_to_He4_Zn59_reaclib, k_He4_Zn59_to_p_Ga62_reaclib,
        N, Ga63, -1, H1, Zn63, -1, Ga63, Zn63, k_n_Ga63_to_p_Zn63_reaclib, k_p_Zn63_to_n_Ga63_reaclib,
        N, Ga63, -1, He4, Cu60, -1, Ga63, Cu60, k_n_Ga63_to_He4_Cu60_reaclib, k_He4_Cu60_to_n_Ga63_reaclib,
        H1, Ga63, -1, He4, Zn60, -1, Ga63, Zn60, k_p_Ga63_to_He4_Zn60_reaclib, k_He4_Zn60_to_p_Ga63_reaclib,
        N, Ga64, -1, H1, Zn64, -1, Ga64, Zn64, k_n_Ga64_to_p_Zn64_reaclib, k_p_Zn64_to_n_Ga64_reaclib,
        N, Ga64, -1, He4, Cu61, -1, Ga64, Cu61, k_n_Ga64_to_He4_Cu61_reaclib, k_He4_Cu61_to_n_Ga64_reaclib,
        H1, Ga64, -1, He4, Zn61, -1, Ga64, Zn61, k_p_Ga64_to_He4_Zn61_reaclib, k_He4_Zn61_to_p_Ga64_reaclib,
        N, Ge63, -1, H1, Ga63, -1, Ge63, Ga63, k_n_Ge63_to_p_Ga63_reaclib, k_p_Ga63_to_n_Ge63_reaclib,
        N, Ge63, -1, He4, Zn60, -1, Ge63, Zn60, k_n_Ge63_to_He4_Zn60_reaclib, k_He4_Zn60_to_n_Ge63_reaclib,
        N, Ge64, -1, H1, Ga64, -1, Ge64, Ga64, k_n_Ge64_to_p_Ga64_reaclib, k_p_Ga64_to_n_Ge64_reaclib,
        N, Ge64, -1, He4, Zn61, -1, Ge64, Zn61, k_n_Ge64_to_He4_Zn61_reaclib, k_He4_Zn61_to_n_Ge64_reaclib,
        He4, He4, He4, C12, -1, -1, C12, -1, k_He4_He4_He4_to_C12_reaclib, k_C12_to_He4_He4_He4_reaclib
    };
}
#endif

#endif
