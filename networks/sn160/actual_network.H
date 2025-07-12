#ifndef actual_network_H
#define actual_network_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Loop.H>
#include <string_view>

#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex;

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
        constexpr_for<1, NumSpec+1>([&] (auto n) {
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

        constexpr_for<1, NumSpec+1>([&] (auto n) {
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
        k_Be7_to_Li7_weak_electron_capture = 1,
        k_C14_to_N14_weak_wc12 = 2,
        k_N13_to_C13_weak_wc12 = 3,
        k_O14_to_N14_weak_wc12 = 4,
        k_O15_to_N15_weak_wc12 = 5,
        k_F17_to_O17_weak_wc12 = 6,
        k_F18_to_O18_weak_wc12 = 7,
        k_Ne18_to_F18_weak_wc12 = 8,
        k_Ne19_to_F19_weak_wc12 = 9,
        k_Zn59_to_Cu59_weak_wc12 = 10,
        k_Ge63_to_Ga63_weak_wc12 = 11,
        k_d_to_n_p = 12,
        k_He3_to_p_d = 13,
        k_He4_to_n_He3 = 14,
        k_He4_to_d_d = 15,
        k_Li6_to_He4_d = 16,
        k_Li7_to_n_Li6 = 17,
        k_Be7_to_p_Li6 = 18,
        k_Be7_to_He4_He3 = 19,
        k_B8_to_p_Be7 = 20,
        k_B8_to_He4_He4_weak_wc12 = 21,
        k_B10_to_p_Be9 = 22,
        k_B10_to_He4_Li6 = 23,
        k_B11_to_n_B10 = 24,
        k_B11_to_He4_Li7 = 25,
        k_C12_to_p_B11 = 26,
        k_C13_to_n_C12 = 27,
        k_C14_to_n_C13 = 28,
        k_N13_to_p_C12 = 29,
        k_N14_to_n_N13 = 30,
        k_N14_to_p_C13 = 31,
        k_N15_to_n_N14 = 32,
        k_N15_to_p_C14 = 33,
        k_O14_to_p_N13 = 34,
        k_O15_to_n_O14 = 35,
        k_O15_to_p_N14 = 36,
        k_O16_to_n_O15 = 37,
        k_O16_to_p_N15 = 38,
        k_O16_to_He4_C12 = 39,
        k_O17_to_n_O16 = 40,
        k_O18_to_n_O17 = 41,
        k_O18_to_He4_C14 = 42,
        k_F17_to_p_O16 = 43,
        k_F18_to_n_F17 = 44,
        k_F18_to_p_O17 = 45,
        k_F18_to_He4_N14 = 46,
        k_F19_to_n_F18 = 47,
        k_F19_to_p_O18 = 48,
        k_F19_to_He4_N15 = 49,
        k_Ne18_to_p_F17 = 50,
        k_Ne18_to_He4_O14 = 51,
        k_Ne19_to_n_Ne18 = 52,
        k_Ne19_to_p_F18 = 53,
        k_Ne19_to_He4_O15 = 54,
        k_Ne20_to_n_Ne19 = 55,
        k_Ne20_to_p_F19 = 56,
        k_Ne20_to_He4_O16 = 57,
        k_Ne21_to_n_Ne20 = 58,
        k_Ne21_to_He4_O17 = 59,
        k_Ne22_to_n_Ne21 = 60,
        k_Ne22_to_He4_O18 = 61,
        k_Na21_to_p_Ne20 = 62,
        k_Na21_to_He4_F17 = 63,
        k_Na22_to_n_Na21 = 64,
        k_Na22_to_p_Ne21 = 65,
        k_Na22_to_He4_F18 = 66,
        k_Na23_to_n_Na22 = 67,
        k_Na23_to_p_Ne22 = 68,
        k_Na23_to_He4_F19 = 69,
        k_Mg23_to_p_Na22 = 70,
        k_Mg23_to_He4_Ne19 = 71,
        k_Mg24_to_n_Mg23 = 72,
        k_Mg24_to_p_Na23 = 73,
        k_Mg24_to_He4_Ne20 = 74,
        k_Mg25_to_n_Mg24 = 75,
        k_Mg25_to_He4_Ne21 = 76,
        k_Mg26_to_n_Mg25 = 77,
        k_Mg26_to_He4_Ne22 = 78,
        k_Al25_to_p_Mg24 = 79,
        k_Al25_to_He4_Na21 = 80,
        k_Al26_to_n_Al25 = 81,
        k_Al26_to_p_Mg25 = 82,
        k_Al26_to_He4_Na22 = 83,
        k_Al27_to_n_Al26 = 84,
        k_Al27_to_p_Mg26 = 85,
        k_Al27_to_He4_Na23 = 86,
        k_Si28_to_p_Al27 = 87,
        k_Si28_to_He4_Mg24 = 88,
        k_Si29_to_n_Si28 = 89,
        k_Si29_to_He4_Mg25 = 90,
        k_Si30_to_n_Si29 = 91,
        k_Si30_to_He4_Mg26 = 92,
        k_Si31_to_n_Si30 = 93,
        k_Si32_to_n_Si31 = 94,
        k_P29_to_p_Si28 = 95,
        k_P29_to_He4_Al25 = 96,
        k_P30_to_n_P29 = 97,
        k_P30_to_p_Si29 = 98,
        k_P30_to_He4_Al26 = 99,
        k_P31_to_n_P30 = 100,
        k_P31_to_p_Si30 = 101,
        k_P31_to_He4_Al27 = 102,
        k_P32_to_n_P31 = 103,
        k_P32_to_p_Si31 = 104,
        k_P33_to_n_P32 = 105,
        k_P33_to_p_Si32 = 106,
        k_S32_to_p_P31 = 107,
        k_S32_to_He4_Si28 = 108,
        k_S33_to_n_S32 = 109,
        k_S33_to_p_P32 = 110,
        k_S33_to_He4_Si29 = 111,
        k_S34_to_n_S33 = 112,
        k_S34_to_p_P33 = 113,
        k_S34_to_He4_Si30 = 114,
        k_S35_to_n_S34 = 115,
        k_S35_to_He4_Si31 = 116,
        k_S36_to_n_S35 = 117,
        k_S36_to_He4_Si32 = 118,
        k_Cl33_to_p_S32 = 119,
        k_Cl33_to_He4_P29 = 120,
        k_Cl34_to_n_Cl33 = 121,
        k_Cl34_to_p_S33 = 122,
        k_Cl34_to_He4_P30 = 123,
        k_Cl35_to_n_Cl34 = 124,
        k_Cl35_to_p_S34 = 125,
        k_Cl35_to_He4_P31 = 126,
        k_Cl36_to_n_Cl35 = 127,
        k_Cl36_to_p_S35 = 128,
        k_Cl36_to_He4_P32 = 129,
        k_Cl37_to_n_Cl36 = 130,
        k_Cl37_to_p_S36 = 131,
        k_Cl37_to_He4_P33 = 132,
        k_Ar36_to_p_Cl35 = 133,
        k_Ar36_to_He4_S32 = 134,
        k_Ar37_to_n_Ar36 = 135,
        k_Ar37_to_p_Cl36 = 136,
        k_Ar37_to_He4_S33 = 137,
        k_Ar38_to_n_Ar37 = 138,
        k_Ar38_to_p_Cl37 = 139,
        k_Ar38_to_He4_S34 = 140,
        k_Ar39_to_n_Ar38 = 141,
        k_Ar39_to_He4_S35 = 142,
        k_Ar40_to_n_Ar39 = 143,
        k_Ar40_to_He4_S36 = 144,
        k_K37_to_p_Ar36 = 145,
        k_K37_to_He4_Cl33 = 146,
        k_K38_to_n_K37 = 147,
        k_K38_to_p_Ar37 = 148,
        k_K38_to_He4_Cl34 = 149,
        k_K39_to_n_K38 = 150,
        k_K39_to_p_Ar38 = 151,
        k_K39_to_He4_Cl35 = 152,
        k_K40_to_n_K39 = 153,
        k_K40_to_p_Ar39 = 154,
        k_K40_to_He4_Cl36 = 155,
        k_K41_to_n_K40 = 156,
        k_K41_to_p_Ar40 = 157,
        k_K41_to_He4_Cl37 = 158,
        k_Ca40_to_p_K39 = 159,
        k_Ca40_to_He4_Ar36 = 160,
        k_Ca41_to_n_Ca40 = 161,
        k_Ca41_to_p_K40 = 162,
        k_Ca41_to_He4_Ar37 = 163,
        k_Ca42_to_n_Ca41 = 164,
        k_Ca42_to_p_K41 = 165,
        k_Ca42_to_He4_Ar38 = 166,
        k_Ca43_to_n_Ca42 = 167,
        k_Ca43_to_He4_Ar39 = 168,
        k_Ca44_to_n_Ca43 = 169,
        k_Ca44_to_He4_Ar40 = 170,
        k_Ca45_to_n_Ca44 = 171,
        k_Ca46_to_n_Ca45 = 172,
        k_Ca47_to_n_Ca46 = 173,
        k_Ca48_to_n_Ca47 = 174,
        k_Sc43_to_p_Ca42 = 175,
        k_Sc43_to_He4_K39 = 176,
        k_Sc44_to_n_Sc43 = 177,
        k_Sc44_to_p_Ca43 = 178,
        k_Sc44_to_He4_K40 = 179,
        k_Sc45_to_n_Sc44 = 180,
        k_Sc45_to_p_Ca44 = 181,
        k_Sc45_to_He4_K41 = 182,
        k_Sc46_to_n_Sc45 = 183,
        k_Sc46_to_p_Ca45 = 184,
        k_Sc47_to_n_Sc46 = 185,
        k_Sc47_to_p_Ca46 = 186,
        k_Sc48_to_n_Sc47 = 187,
        k_Sc48_to_p_Ca47 = 188,
        k_Sc49_to_n_Sc48 = 189,
        k_Sc49_to_p_Ca48 = 190,
        k_Ti44_to_p_Sc43 = 191,
        k_Ti44_to_He4_Ca40 = 192,
        k_Ti45_to_n_Ti44 = 193,
        k_Ti45_to_p_Sc44 = 194,
        k_Ti45_to_He4_Ca41 = 195,
        k_Ti46_to_n_Ti45 = 196,
        k_Ti46_to_p_Sc45 = 197,
        k_Ti46_to_He4_Ca42 = 198,
        k_Ti47_to_n_Ti46 = 199,
        k_Ti47_to_p_Sc46 = 200,
        k_Ti47_to_He4_Ca43 = 201,
        k_Ti48_to_n_Ti47 = 202,
        k_Ti48_to_p_Sc47 = 203,
        k_Ti48_to_He4_Ca44 = 204,
        k_Ti49_to_n_Ti48 = 205,
        k_Ti49_to_p_Sc48 = 206,
        k_Ti49_to_He4_Ca45 = 207,
        k_Ti50_to_n_Ti49 = 208,
        k_Ti50_to_p_Sc49 = 209,
        k_Ti50_to_He4_Ca46 = 210,
        k_Ti51_to_n_Ti50 = 211,
        k_Ti51_to_He4_Ca47 = 212,
        k_V46_to_p_Ti45 = 213,
        k_V47_to_n_V46 = 214,
        k_V47_to_p_Ti46 = 215,
        k_V47_to_He4_Sc43 = 216,
        k_V48_to_n_V47 = 217,
        k_V48_to_p_Ti47 = 218,
        k_V48_to_He4_Sc44 = 219,
        k_V49_to_n_V48 = 220,
        k_V49_to_p_Ti48 = 221,
        k_V49_to_He4_Sc45 = 222,
        k_V50_to_n_V49 = 223,
        k_V50_to_p_Ti49 = 224,
        k_V50_to_He4_Sc46 = 225,
        k_V51_to_n_V50 = 226,
        k_V51_to_p_Ti50 = 227,
        k_V51_to_He4_Sc47 = 228,
        k_V52_to_n_V51 = 229,
        k_V52_to_p_Ti51 = 230,
        k_V52_to_He4_Sc48 = 231,
        k_Cr48_to_p_V47 = 232,
        k_Cr48_to_He4_Ti44 = 233,
        k_Cr49_to_n_Cr48 = 234,
        k_Cr49_to_p_V48 = 235,
        k_Cr49_to_He4_Ti45 = 236,
        k_Cr50_to_n_Cr49 = 237,
        k_Cr50_to_p_V49 = 238,
        k_Cr50_to_He4_Ti46 = 239,
        k_Cr51_to_n_Cr50 = 240,
        k_Cr51_to_p_V50 = 241,
        k_Cr51_to_He4_Ti47 = 242,
        k_Cr52_to_n_Cr51 = 243,
        k_Cr52_to_p_V51 = 244,
        k_Cr52_to_He4_Ti48 = 245,
        k_Cr53_to_n_Cr52 = 246,
        k_Cr53_to_p_V52 = 247,
        k_Cr53_to_He4_Ti49 = 248,
        k_Cr54_to_n_Cr53 = 249,
        k_Cr54_to_He4_Ti50 = 250,
        k_Mn50_to_p_Cr49 = 251,
        k_Mn50_to_He4_V46 = 252,
        k_Mn51_to_n_Mn50 = 253,
        k_Mn51_to_p_Cr50 = 254,
        k_Mn51_to_He4_V47 = 255,
        k_Mn52_to_n_Mn51 = 256,
        k_Mn52_to_p_Cr51 = 257,
        k_Mn52_to_He4_V48 = 258,
        k_Mn53_to_n_Mn52 = 259,
        k_Mn53_to_p_Cr52 = 260,
        k_Mn53_to_He4_V49 = 261,
        k_Mn54_to_n_Mn53 = 262,
        k_Mn54_to_p_Cr53 = 263,
        k_Mn54_to_He4_V50 = 264,
        k_Mn55_to_n_Mn54 = 265,
        k_Mn55_to_p_Cr54 = 266,
        k_Mn55_to_He4_V51 = 267,
        k_Fe52_to_p_Mn51 = 268,
        k_Fe52_to_He4_Cr48 = 269,
        k_Fe53_to_n_Fe52 = 270,
        k_Fe53_to_p_Mn52 = 271,
        k_Fe53_to_He4_Cr49 = 272,
        k_Fe54_to_n_Fe53 = 273,
        k_Fe54_to_p_Mn53 = 274,
        k_Fe54_to_He4_Cr50 = 275,
        k_Fe55_to_n_Fe54 = 276,
        k_Fe55_to_p_Mn54 = 277,
        k_Fe55_to_He4_Cr51 = 278,
        k_Fe56_to_n_Fe55 = 279,
        k_Fe56_to_p_Mn55 = 280,
        k_Fe56_to_He4_Cr52 = 281,
        k_Fe57_to_n_Fe56 = 282,
        k_Fe57_to_He4_Cr53 = 283,
        k_Fe58_to_n_Fe57 = 284,
        k_Fe58_to_He4_Cr54 = 285,
        k_Co53_to_p_Fe52 = 286,
        k_Co54_to_n_Co53 = 287,
        k_Co54_to_p_Fe53 = 288,
        k_Co54_to_He4_Mn50 = 289,
        k_Co55_to_n_Co54 = 290,
        k_Co55_to_p_Fe54 = 291,
        k_Co55_to_He4_Mn51 = 292,
        k_Co56_to_n_Co55 = 293,
        k_Co56_to_p_Fe55 = 294,
        k_Co56_to_He4_Mn52 = 295,
        k_Co57_to_n_Co56 = 296,
        k_Co57_to_p_Fe56 = 297,
        k_Co57_to_He4_Mn53 = 298,
        k_Co58_to_n_Co57 = 299,
        k_Co58_to_p_Fe57 = 300,
        k_Co58_to_He4_Mn54 = 301,
        k_Co59_to_n_Co58 = 302,
        k_Co59_to_p_Fe58 = 303,
        k_Co59_to_He4_Mn55 = 304,
        k_Ni56_to_p_Co55 = 305,
        k_Ni56_to_He4_Fe52 = 306,
        k_Ni57_to_n_Ni56 = 307,
        k_Ni57_to_p_Co56 = 308,
        k_Ni57_to_He4_Fe53 = 309,
        k_Ni58_to_n_Ni57 = 310,
        k_Ni58_to_p_Co57 = 311,
        k_Ni58_to_He4_Fe54 = 312,
        k_Ni59_to_n_Ni58 = 313,
        k_Ni59_to_p_Co58 = 314,
        k_Ni59_to_He4_Fe55 = 315,
        k_Ni60_to_n_Ni59 = 316,
        k_Ni60_to_p_Co59 = 317,
        k_Ni60_to_He4_Fe56 = 318,
        k_Ni61_to_n_Ni60 = 319,
        k_Ni61_to_He4_Fe57 = 320,
        k_Ni62_to_n_Ni61 = 321,
        k_Ni62_to_He4_Fe58 = 322,
        k_Ni63_to_n_Ni62 = 323,
        k_Ni64_to_n_Ni63 = 324,
        k_Cu57_to_p_Ni56 = 325,
        k_Cu57_to_He4_Co53 = 326,
        k_Cu58_to_n_Cu57 = 327,
        k_Cu58_to_p_Ni57 = 328,
        k_Cu58_to_He4_Co54 = 329,
        k_Cu59_to_n_Cu58 = 330,
        k_Cu59_to_p_Ni58 = 331,
        k_Cu59_to_He4_Co55 = 332,
        k_Cu60_to_n_Cu59 = 333,
        k_Cu60_to_p_Ni59 = 334,
        k_Cu60_to_He4_Co56 = 335,
        k_Cu61_to_n_Cu60 = 336,
        k_Cu61_to_p_Ni60 = 337,
        k_Cu61_to_He4_Co57 = 338,
        k_Cu62_to_n_Cu61 = 339,
        k_Cu62_to_p_Ni61 = 340,
        k_Cu62_to_He4_Co58 = 341,
        k_Cu63_to_n_Cu62 = 342,
        k_Cu63_to_p_Ni62 = 343,
        k_Cu63_to_He4_Co59 = 344,
        k_Cu64_to_n_Cu63 = 345,
        k_Cu64_to_p_Ni63 = 346,
        k_Cu65_to_n_Cu64 = 347,
        k_Cu65_to_p_Ni64 = 348,
        k_Zn59_to_p_Cu58 = 349,
        k_Zn59_to_p_Ni58_weak_wc12 = 350,
        k_Zn60_to_n_Zn59 = 351,
        k_Zn60_to_p_Cu59 = 352,
        k_Zn60_to_He4_Ni56 = 353,
        k_Zn61_to_n_Zn60 = 354,
        k_Zn61_to_p_Cu60 = 355,
        k_Zn61_to_He4_Ni57 = 356,
        k_Zn62_to_n_Zn61 = 357,
        k_Zn62_to_p_Cu61 = 358,
        k_Zn62_to_He4_Ni58 = 359,
        k_Zn63_to_n_Zn62 = 360,
        k_Zn63_to_p_Cu62 = 361,
        k_Zn63_to_He4_Ni59 = 362,
        k_Zn64_to_n_Zn63 = 363,
        k_Zn64_to_p_Cu63 = 364,
        k_Zn64_to_He4_Ni60 = 365,
        k_Zn65_to_n_Zn64 = 366,
        k_Zn65_to_p_Cu64 = 367,
        k_Zn65_to_He4_Ni61 = 368,
        k_Zn66_to_n_Zn65 = 369,
        k_Zn66_to_p_Cu65 = 370,
        k_Zn66_to_He4_Ni62 = 371,
        k_Ga62_to_p_Zn61 = 372,
        k_Ga62_to_He4_Cu58 = 373,
        k_Ga63_to_n_Ga62 = 374,
        k_Ga63_to_p_Zn62 = 375,
        k_Ga63_to_He4_Cu59 = 376,
        k_Ga64_to_n_Ga63 = 377,
        k_Ga64_to_p_Zn63 = 378,
        k_Ga64_to_He4_Cu60 = 379,
        k_Ge63_to_p_Ga62 = 380,
        k_Ge63_to_He4_Zn59 = 381,
        k_Ge64_to_n_Ge63 = 382,
        k_Ge64_to_p_Ga63 = 383,
        k_Ge64_to_He4_Zn60 = 384,
        k_Li6_to_n_p_He4 = 385,
        k_Be9_to_n_He4_He4 = 386,
        k_C12_to_He4_He4_He4 = 387,
        k_n_p_to_d = 388,
        k_p_p_to_d_weak_bet_pos_ = 389,
        k_p_p_to_d_weak_electron_capture = 390,
        k_p_d_to_He3 = 391,
        k_d_d_to_He4 = 392,
        k_He4_d_to_Li6 = 393,
        k_n_He3_to_He4 = 394,
        k_p_He3_to_He4_weak_bet_pos_ = 395,
        k_He4_He3_to_Be7 = 396,
        k_n_Li6_to_Li7 = 397,
        k_p_Li6_to_Be7 = 398,
        k_He4_Li6_to_B10 = 399,
        k_He4_Li7_to_B11 = 400,
        k_p_Be7_to_B8 = 401,
        k_p_Be9_to_B10 = 402,
        k_n_B10_to_B11 = 403,
        k_p_B11_to_C12 = 404,
        k_n_C12_to_C13 = 405,
        k_p_C12_to_N13 = 406,
        k_He4_C12_to_O16 = 407,
        k_n_C13_to_C14 = 408,
        k_p_C13_to_N14 = 409,
        k_p_C14_to_N15 = 410,
        k_He4_C14_to_O18 = 411,
        k_n_N13_to_N14 = 412,
        k_p_N13_to_O14 = 413,
        k_n_N14_to_N15 = 414,
        k_p_N14_to_O15 = 415,
        k_He4_N14_to_F18 = 416,
        k_p_N15_to_O16 = 417,
        k_He4_N15_to_F19 = 418,
        k_n_O14_to_O15 = 419,
        k_He4_O14_to_Ne18 = 420,
        k_n_O15_to_O16 = 421,
        k_He4_O15_to_Ne19 = 422,
        k_n_O16_to_O17 = 423,
        k_p_O16_to_F17 = 424,
        k_He4_O16_to_Ne20 = 425,
        k_n_O17_to_O18 = 426,
        k_p_O17_to_F18 = 427,
        k_He4_O17_to_Ne21 = 428,
        k_p_O18_to_F19 = 429,
        k_He4_O18_to_Ne22 = 430,
        k_n_F17_to_F18 = 431,
        k_p_F17_to_Ne18 = 432,
        k_He4_F17_to_Na21 = 433,
        k_n_F18_to_F19 = 434,
        k_p_F18_to_Ne19 = 435,
        k_He4_F18_to_Na22 = 436,
        k_p_F19_to_Ne20 = 437,
        k_He4_F19_to_Na23 = 438,
        k_n_Ne18_to_Ne19 = 439,
        k_n_Ne19_to_Ne20 = 440,
        k_He4_Ne19_to_Mg23 = 441,
        k_n_Ne20_to_Ne21 = 442,
        k_p_Ne20_to_Na21 = 443,
        k_He4_Ne20_to_Mg24 = 444,
        k_n_Ne21_to_Ne22 = 445,
        k_p_Ne21_to_Na22 = 446,
        k_He4_Ne21_to_Mg25 = 447,
        k_p_Ne22_to_Na23 = 448,
        k_He4_Ne22_to_Mg26 = 449,
        k_n_Na21_to_Na22 = 450,
        k_He4_Na21_to_Al25 = 451,
        k_n_Na22_to_Na23 = 452,
        k_p_Na22_to_Mg23 = 453,
        k_He4_Na22_to_Al26 = 454,
        k_p_Na23_to_Mg24 = 455,
        k_He4_Na23_to_Al27 = 456,
        k_n_Mg23_to_Mg24 = 457,
        k_n_Mg24_to_Mg25 = 458,
        k_p_Mg24_to_Al25 = 459,
        k_He4_Mg24_to_Si28 = 460,
        k_n_Mg25_to_Mg26 = 461,
        k_p_Mg25_to_Al26 = 462,
        k_He4_Mg25_to_Si29 = 463,
        k_p_Mg26_to_Al27 = 464,
        k_He4_Mg26_to_Si30 = 465,
        k_n_Al25_to_Al26 = 466,
        k_He4_Al25_to_P29 = 467,
        k_n_Al26_to_Al27 = 468,
        k_He4_Al26_to_P30 = 469,
        k_p_Al27_to_Si28 = 470,
        k_He4_Al27_to_P31 = 471,
        k_n_Si28_to_Si29 = 472,
        k_p_Si28_to_P29 = 473,
        k_He4_Si28_to_S32 = 474,
        k_n_Si29_to_Si30 = 475,
        k_p_Si29_to_P30 = 476,
        k_He4_Si29_to_S33 = 477,
        k_n_Si30_to_Si31 = 478,
        k_p_Si30_to_P31 = 479,
        k_He4_Si30_to_S34 = 480,
        k_n_Si31_to_Si32 = 481,
        k_p_Si31_to_P32 = 482,
        k_He4_Si31_to_S35 = 483,
        k_p_Si32_to_P33 = 484,
        k_He4_Si32_to_S36 = 485,
        k_n_P29_to_P30 = 486,
        k_He4_P29_to_Cl33 = 487,
        k_n_P30_to_P31 = 488,
        k_He4_P30_to_Cl34 = 489,
        k_n_P31_to_P32 = 490,
        k_p_P31_to_S32 = 491,
        k_He4_P31_to_Cl35 = 492,
        k_n_P32_to_P33 = 493,
        k_p_P32_to_S33 = 494,
        k_He4_P32_to_Cl36 = 495,
        k_p_P33_to_S34 = 496,
        k_He4_P33_to_Cl37 = 497,
        k_n_S32_to_S33 = 498,
        k_p_S32_to_Cl33 = 499,
        k_He4_S32_to_Ar36 = 500,
        k_n_S33_to_S34 = 501,
        k_p_S33_to_Cl34 = 502,
        k_He4_S33_to_Ar37 = 503,
        k_n_S34_to_S35 = 504,
        k_p_S34_to_Cl35 = 505,
        k_He4_S34_to_Ar38 = 506,
        k_n_S35_to_S36 = 507,
        k_p_S35_to_Cl36 = 508,
        k_He4_S35_to_Ar39 = 509,
        k_p_S36_to_Cl37 = 510,
        k_He4_S36_to_Ar40 = 511,
        k_n_Cl33_to_Cl34 = 512,
        k_He4_Cl33_to_K37 = 513,
        k_n_Cl34_to_Cl35 = 514,
        k_He4_Cl34_to_K38 = 515,
        k_n_Cl35_to_Cl36 = 516,
        k_p_Cl35_to_Ar36 = 517,
        k_He4_Cl35_to_K39 = 518,
        k_n_Cl36_to_Cl37 = 519,
        k_p_Cl36_to_Ar37 = 520,
        k_He4_Cl36_to_K40 = 521,
        k_p_Cl37_to_Ar38 = 522,
        k_He4_Cl37_to_K41 = 523,
        k_n_Ar36_to_Ar37 = 524,
        k_p_Ar36_to_K37 = 525,
        k_He4_Ar36_to_Ca40 = 526,
        k_n_Ar37_to_Ar38 = 527,
        k_p_Ar37_to_K38 = 528,
        k_He4_Ar37_to_Ca41 = 529,
        k_n_Ar38_to_Ar39 = 530,
        k_p_Ar38_to_K39 = 531,
        k_He4_Ar38_to_Ca42 = 532,
        k_n_Ar39_to_Ar40 = 533,
        k_p_Ar39_to_K40 = 534,
        k_He4_Ar39_to_Ca43 = 535,
        k_p_Ar40_to_K41 = 536,
        k_He4_Ar40_to_Ca44 = 537,
        k_n_K37_to_K38 = 538,
        k_n_K38_to_K39 = 539,
        k_n_K39_to_K40 = 540,
        k_p_K39_to_Ca40 = 541,
        k_He4_K39_to_Sc43 = 542,
        k_n_K40_to_K41 = 543,
        k_p_K40_to_Ca41 = 544,
        k_He4_K40_to_Sc44 = 545,
        k_p_K41_to_Ca42 = 546,
        k_He4_K41_to_Sc45 = 547,
        k_n_Ca40_to_Ca41 = 548,
        k_He4_Ca40_to_Ti44 = 549,
        k_n_Ca41_to_Ca42 = 550,
        k_He4_Ca41_to_Ti45 = 551,
        k_n_Ca42_to_Ca43 = 552,
        k_p_Ca42_to_Sc43 = 553,
        k_He4_Ca42_to_Ti46 = 554,
        k_n_Ca43_to_Ca44 = 555,
        k_p_Ca43_to_Sc44 = 556,
        k_He4_Ca43_to_Ti47 = 557,
        k_n_Ca44_to_Ca45 = 558,
        k_p_Ca44_to_Sc45 = 559,
        k_He4_Ca44_to_Ti48 = 560,
        k_n_Ca45_to_Ca46 = 561,
        k_p_Ca45_to_Sc46 = 562,
        k_He4_Ca45_to_Ti49 = 563,
        k_n_Ca46_to_Ca47 = 564,
        k_p_Ca46_to_Sc47 = 565,
        k_He4_Ca46_to_Ti50 = 566,
        k_n_Ca47_to_Ca48 = 567,
        k_p_Ca47_to_Sc48 = 568,
        k_He4_Ca47_to_Ti51 = 569,
        k_p_Ca48_to_Sc49 = 570,
        k_n_Sc43_to_Sc44 = 571,
        k_p_Sc43_to_Ti44 = 572,
        k_He4_Sc43_to_V47 = 573,
        k_n_Sc44_to_Sc45 = 574,
        k_p_Sc44_to_Ti45 = 575,
        k_He4_Sc44_to_V48 = 576,
        k_n_Sc45_to_Sc46 = 577,
        k_p_Sc45_to_Ti46 = 578,
        k_He4_Sc45_to_V49 = 579,
        k_n_Sc46_to_Sc47 = 580,
        k_p_Sc46_to_Ti47 = 581,
        k_He4_Sc46_to_V50 = 582,
        k_n_Sc47_to_Sc48 = 583,
        k_p_Sc47_to_Ti48 = 584,
        k_He4_Sc47_to_V51 = 585,
        k_n_Sc48_to_Sc49 = 586,
        k_p_Sc48_to_Ti49 = 587,
        k_He4_Sc48_to_V52 = 588,
        k_p_Sc49_to_Ti50 = 589,
        k_n_Ti44_to_Ti45 = 590,
        k_He4_Ti44_to_Cr48 = 591,
        k_n_Ti45_to_Ti46 = 592,
        k_p_Ti45_to_V46 = 593,
        k_He4_Ti45_to_Cr49 = 594,
        k_n_Ti46_to_Ti47 = 595,
        k_p_Ti46_to_V47 = 596,
        k_He4_Ti46_to_Cr50 = 597,
        k_n_Ti47_to_Ti48 = 598,
        k_p_Ti47_to_V48 = 599,
        k_He4_Ti47_to_Cr51 = 600,
        k_n_Ti48_to_Ti49 = 601,
        k_p_Ti48_to_V49 = 602,
        k_He4_Ti48_to_Cr52 = 603,
        k_n_Ti49_to_Ti50 = 604,
        k_p_Ti49_to_V50 = 605,
        k_He4_Ti49_to_Cr53 = 606,
        k_n_Ti50_to_Ti51 = 607,
        k_p_Ti50_to_V51 = 608,
        k_He4_Ti50_to_Cr54 = 609,
        k_p_Ti51_to_V52 = 610,
        k_n_V46_to_V47 = 611,
        k_He4_V46_to_Mn50 = 612,
        k_n_V47_to_V48 = 613,
        k_p_V47_to_Cr48 = 614,
        k_He4_V47_to_Mn51 = 615,
        k_n_V48_to_V49 = 616,
        k_p_V48_to_Cr49 = 617,
        k_He4_V48_to_Mn52 = 618,
        k_n_V49_to_V50 = 619,
        k_p_V49_to_Cr50 = 620,
        k_He4_V49_to_Mn53 = 621,
        k_n_V50_to_V51 = 622,
        k_p_V50_to_Cr51 = 623,
        k_He4_V50_to_Mn54 = 624,
        k_n_V51_to_V52 = 625,
        k_p_V51_to_Cr52 = 626,
        k_He4_V51_to_Mn55 = 627,
        k_p_V52_to_Cr53 = 628,
        k_n_Cr48_to_Cr49 = 629,
        k_He4_Cr48_to_Fe52 = 630,
        k_n_Cr49_to_Cr50 = 631,
        k_p_Cr49_to_Mn50 = 632,
        k_He4_Cr49_to_Fe53 = 633,
        k_n_Cr50_to_Cr51 = 634,
        k_p_Cr50_to_Mn51 = 635,
        k_He4_Cr50_to_Fe54 = 636,
        k_n_Cr51_to_Cr52 = 637,
        k_p_Cr51_to_Mn52 = 638,
        k_He4_Cr51_to_Fe55 = 639,
        k_n_Cr52_to_Cr53 = 640,
        k_p_Cr52_to_Mn53 = 641,
        k_He4_Cr52_to_Fe56 = 642,
        k_n_Cr53_to_Cr54 = 643,
        k_p_Cr53_to_Mn54 = 644,
        k_He4_Cr53_to_Fe57 = 645,
        k_p_Cr54_to_Mn55 = 646,
        k_He4_Cr54_to_Fe58 = 647,
        k_n_Mn50_to_Mn51 = 648,
        k_He4_Mn50_to_Co54 = 649,
        k_n_Mn51_to_Mn52 = 650,
        k_p_Mn51_to_Fe52 = 651,
        k_He4_Mn51_to_Co55 = 652,
        k_n_Mn52_to_Mn53 = 653,
        k_p_Mn52_to_Fe53 = 654,
        k_He4_Mn52_to_Co56 = 655,
        k_n_Mn53_to_Mn54 = 656,
        k_p_Mn53_to_Fe54 = 657,
        k_He4_Mn53_to_Co57 = 658,
        k_n_Mn54_to_Mn55 = 659,
        k_p_Mn54_to_Fe55 = 660,
        k_He4_Mn54_to_Co58 = 661,
        k_p_Mn55_to_Fe56 = 662,
        k_He4_Mn55_to_Co59 = 663,
        k_n_Fe52_to_Fe53 = 664,
        k_p_Fe52_to_Co53 = 665,
        k_He4_Fe52_to_Ni56 = 666,
        k_n_Fe53_to_Fe54 = 667,
        k_p_Fe53_to_Co54 = 668,
        k_He4_Fe53_to_Ni57 = 669,
        k_n_Fe54_to_Fe55 = 670,
        k_p_Fe54_to_Co55 = 671,
        k_He4_Fe54_to_Ni58 = 672,
        k_n_Fe55_to_Fe56 = 673,
        k_p_Fe55_to_Co56 = 674,
        k_He4_Fe55_to_Ni59 = 675,
        k_n_Fe56_to_Fe57 = 676,
        k_p_Fe56_to_Co57 = 677,
        k_He4_Fe56_to_Ni60 = 678,
        k_n_Fe57_to_Fe58 = 679,
        k_p_Fe57_to_Co58 = 680,
        k_He4_Fe57_to_Ni61 = 681,
        k_p_Fe58_to_Co59 = 682,
        k_He4_Fe58_to_Ni62 = 683,
        k_n_Co53_to_Co54 = 684,
        k_He4_Co53_to_Cu57 = 685,
        k_n_Co54_to_Co55 = 686,
        k_He4_Co54_to_Cu58 = 687,
        k_n_Co55_to_Co56 = 688,
        k_p_Co55_to_Ni56 = 689,
        k_He4_Co55_to_Cu59 = 690,
        k_n_Co56_to_Co57 = 691,
        k_p_Co56_to_Ni57 = 692,
        k_He4_Co56_to_Cu60 = 693,
        k_n_Co57_to_Co58 = 694,
        k_p_Co57_to_Ni58 = 695,
        k_He4_Co57_to_Cu61 = 696,
        k_n_Co58_to_Co59 = 697,
        k_p_Co58_to_Ni59 = 698,
        k_He4_Co58_to_Cu62 = 699,
        k_p_Co59_to_Ni60 = 700,
        k_He4_Co59_to_Cu63 = 701,
        k_n_Ni56_to_Ni57 = 702,
        k_p_Ni56_to_Cu57 = 703,
        k_He4_Ni56_to_Zn60 = 704,
        k_n_Ni57_to_Ni58 = 705,
        k_p_Ni57_to_Cu58 = 706,
        k_He4_Ni57_to_Zn61 = 707,
        k_n_Ni58_to_Ni59 = 708,
        k_p_Ni58_to_Cu59 = 709,
        k_He4_Ni58_to_Zn62 = 710,
        k_n_Ni59_to_Ni60 = 711,
        k_p_Ni59_to_Cu60 = 712,
        k_He4_Ni59_to_Zn63 = 713,
        k_n_Ni60_to_Ni61 = 714,
        k_p_Ni60_to_Cu61 = 715,
        k_He4_Ni60_to_Zn64 = 716,
        k_n_Ni61_to_Ni62 = 717,
        k_p_Ni61_to_Cu62 = 718,
        k_He4_Ni61_to_Zn65 = 719,
        k_n_Ni62_to_Ni63 = 720,
        k_p_Ni62_to_Cu63 = 721,
        k_He4_Ni62_to_Zn66 = 722,
        k_n_Ni63_to_Ni64 = 723,
        k_p_Ni63_to_Cu64 = 724,
        k_p_Ni64_to_Cu65 = 725,
        k_n_Cu57_to_Cu58 = 726,
        k_n_Cu58_to_Cu59 = 727,
        k_p_Cu58_to_Zn59 = 728,
        k_He4_Cu58_to_Ga62 = 729,
        k_n_Cu59_to_Cu60 = 730,
        k_p_Cu59_to_Zn60 = 731,
        k_He4_Cu59_to_Ga63 = 732,
        k_n_Cu60_to_Cu61 = 733,
        k_p_Cu60_to_Zn61 = 734,
        k_He4_Cu60_to_Ga64 = 735,
        k_n_Cu61_to_Cu62 = 736,
        k_p_Cu61_to_Zn62 = 737,
        k_n_Cu62_to_Cu63 = 738,
        k_p_Cu62_to_Zn63 = 739,
        k_n_Cu63_to_Cu64 = 740,
        k_p_Cu63_to_Zn64 = 741,
        k_n_Cu64_to_Cu65 = 742,
        k_p_Cu64_to_Zn65 = 743,
        k_p_Cu65_to_Zn66 = 744,
        k_n_Zn59_to_Zn60 = 745,
        k_He4_Zn59_to_Ge63 = 746,
        k_n_Zn60_to_Zn61 = 747,
        k_He4_Zn60_to_Ge64 = 748,
        k_n_Zn61_to_Zn62 = 749,
        k_p_Zn61_to_Ga62 = 750,
        k_n_Zn62_to_Zn63 = 751,
        k_p_Zn62_to_Ga63 = 752,
        k_n_Zn63_to_Zn64 = 753,
        k_p_Zn63_to_Ga64 = 754,
        k_n_Zn64_to_Zn65 = 755,
        k_n_Zn65_to_Zn66 = 756,
        k_n_Ga62_to_Ga63 = 757,
        k_p_Ga62_to_Ge63 = 758,
        k_n_Ga63_to_Ga64 = 759,
        k_p_Ga63_to_Ge64 = 760,
        k_n_Ge63_to_Ge64 = 761,
        k_d_d_to_n_He3 = 762,
        k_n_He3_to_d_d = 763,
        k_d_He3_to_p_He4 = 764,
        k_He4_He3_to_p_Li6 = 765,
        k_p_He4_to_d_He3 = 766,
        k_He4_He4_to_n_Be7 = 767,
        k_He4_He4_to_p_Li7 = 768,
        k_p_Li6_to_He4_He3 = 769,
        k_d_Li6_to_n_Be7 = 770,
        k_d_Li6_to_p_Li7 = 771,
        k_He4_Li6_to_p_Be9 = 772,
        k_p_Li7_to_n_Be7 = 773,
        k_p_Li7_to_d_Li6 = 774,
        k_p_Li7_to_He4_He4 = 775,
        k_He4_Li7_to_n_B10 = 776,
        k_n_Be7_to_p_Li7 = 777,
        k_n_Be7_to_d_Li6 = 778,
        k_n_Be7_to_He4_He4 = 779,
        k_He4_Be7_to_p_B10 = 780,
        k_p_Be9_to_He4_Li6 = 781,
        k_He4_Be9_to_n_C12 = 782,
        k_n_B10_to_He4_Li7 = 783,
        k_p_B10_to_He4_Be7 = 784,
        k_He4_B10_to_n_N13 = 785,
        k_He4_B10_to_p_C13 = 786,
        k_He4_B11_to_n_N14 = 787,
        k_He4_B11_to_p_C14 = 788,
        k_n_C12_to_He4_Be9 = 789,
        k_He4_C12_to_n_O15 = 790,
        k_He4_C12_to_p_N15 = 791,
        k_C12_C12_to_n_Mg23 = 792,
        k_C12_C12_to_p_Na23 = 793,
        k_C12_C12_to_He4_Ne20 = 794,
        k_p_C13_to_n_N13 = 795,
        k_p_C13_to_He4_B10 = 796,
        k_d_C13_to_n_N14 = 797,
        k_He4_C13_to_n_O16 = 798,
        k_p_C14_to_n_N14 = 799,
        k_p_C14_to_He4_B11 = 800,
        k_d_C14_to_n_N15 = 801,
        k_He4_C14_to_n_O17 = 802,
        k_n_N13_to_p_C13 = 803,
        k_n_N13_to_He4_B10 = 804,
        k_He4_N13_to_p_O16 = 805,
        k_n_N14_to_p_C14 = 806,
        k_n_N14_to_d_C13 = 807,
        k_n_N14_to_He4_B11 = 808,
        k_p_N14_to_n_O14 = 809,
        k_He4_N14_to_n_F17 = 810,
        k_He4_N14_to_p_O17 = 811,
        k_n_N15_to_d_C14 = 812,
        k_p_N15_to_n_O15 = 813,
        k_p_N15_to_He4_C12 = 814,
        k_He4_N15_to_n_F18 = 815,
        k_He4_N15_to_p_O18 = 816,
        k_n_O14_to_p_N14 = 817,
        k_He4_O14_to_p_F17 = 818,
        k_n_O15_to_p_N15 = 819,
        k_n_O15_to_He4_C12 = 820,
        k_He4_O15_to_n_Ne18 = 821,
        k_He4_O15_to_p_F18 = 822,
        k_n_O16_to_He4_C13 = 823,
        k_p_O16_to_He4_N13 = 824,
        k_He4_O16_to_n_Ne19 = 825,
        k_He4_O16_to_p_F19 = 826,
        k_C12_O16_to_p_Al27 = 827,
        k_C12_O16_to_He4_Mg24 = 828,
        k_O16_O16_to_p_P31 = 829,
        k_O16_O16_to_He4_Si28 = 830,
        k_n_O17_to_He4_C14 = 831,
        k_p_O17_to_n_F17 = 832,
        k_p_O17_to_He4_N14 = 833,
        k_He4_O17_to_n_Ne20 = 834,
        k_p_O18_to_n_F18 = 835,
        k_p_O18_to_He4_N15 = 836,
        k_He4_O18_to_n_Ne21 = 837,
        k_n_F17_to_p_O17 = 838,
        k_n_F17_to_He4_N14 = 839,
        k_p_F17_to_He4_O14 = 840,
        k_He4_F17_to_p_Ne20 = 841,
        k_n_F18_to_p_O18 = 842,
        k_n_F18_to_He4_N15 = 843,
        k_p_F18_to_n_Ne18 = 844,
        k_p_F18_to_He4_O15 = 845,
        k_He4_F18_to_n_Na21 = 846,
        k_He4_F18_to_p_Ne21 = 847,
        k_p_F19_to_n_Ne19 = 848,
        k_p_F19_to_He4_O16 = 849,
        k_He4_F19_to_n_Na22 = 850,
        k_He4_F19_to_p_Ne22 = 851,
        k_n_Ne18_to_p_F18 = 852,
        k_n_Ne18_to_He4_O15 = 853,
        k_He4_Ne18_to_p_Na21 = 854,
        k_n_Ne19_to_p_F19 = 855,
        k_n_Ne19_to_He4_O16 = 856,
        k_He4_Ne19_to_p_Na22 = 857,
        k_n_Ne20_to_He4_O17 = 858,
        k_p_Ne20_to_He4_F17 = 859,
        k_He4_Ne20_to_n_Mg23 = 860,
        k_He4_Ne20_to_p_Na23 = 861,
        k_He4_Ne20_to_C12_C12 = 862,
        k_C12_Ne20_to_p_P31 = 863,
        k_C12_Ne20_to_He4_Si28 = 864,
        k_n_Ne21_to_He4_O18 = 865,
        k_p_Ne21_to_n_Na21 = 866,
        k_p_Ne21_to_He4_F18 = 867,
        k_He4_Ne21_to_n_Mg24 = 868,
        k_p_Ne22_to_n_Na22 = 869,
        k_p_Ne22_to_He4_F19 = 870,
        k_He4_Ne22_to_n_Mg25 = 871,
        k_n_Na21_to_p_Ne21 = 872,
        k_n_Na21_to_He4_F18 = 873,
        k_p_Na21_to_He4_Ne18 = 874,
        k_He4_Na21_to_p_Mg24 = 875,
        k_n_Na22_to_p_Ne22 = 876,
        k_n_Na22_to_He4_F19 = 877,
        k_p_Na22_to_He4_Ne19 = 878,
        k_He4_Na22_to_n_Al25 = 879,
        k_He4_Na22_to_p_Mg25 = 880,
        k_p_Na23_to_n_Mg23 = 881,
        k_p_Na23_to_He4_Ne20 = 882,
        k_p_Na23_to_C12_C12 = 883,
        k_He4_Na23_to_n_Al26 = 884,
        k_He4_Na23_to_p_Mg26 = 885,
        k_n_Mg23_to_p_Na23 = 886,
        k_n_Mg23_to_He4_Ne20 = 887,
        k_n_Mg23_to_C12_C12 = 888,
        k_He4_Mg23_to_p_Al26 = 889,
        k_n_Mg24_to_He4_Ne21 = 890,
        k_p_Mg24_to_He4_Na21 = 891,
        k_He4_Mg24_to_p_Al27 = 892,
        k_He4_Mg24_to_C12_O16 = 893,
        k_n_Mg25_to_He4_Ne22 = 894,
        k_p_Mg25_to_n_Al25 = 895,
        k_p_Mg25_to_He4_Na22 = 896,
        k_He4_Mg25_to_n_Si28 = 897,
        k_p_Mg26_to_n_Al26 = 898,
        k_p_Mg26_to_He4_Na23 = 899,
        k_He4_Mg26_to_n_Si29 = 900,
        k_n_Al25_to_p_Mg25 = 901,
        k_n_Al25_to_He4_Na22 = 902,
        k_He4_Al25_to_p_Si28 = 903,
        k_n_Al26_to_p_Mg26 = 904,
        k_n_Al26_to_He4_Na23 = 905,
        k_p_Al26_to_He4_Mg23 = 906,
        k_He4_Al26_to_n_P29 = 907,
        k_He4_Al26_to_p_Si29 = 908,
        k_p_Al27_to_He4_Mg24 = 909,
        k_p_Al27_to_C12_O16 = 910,
        k_He4_Al27_to_n_P30 = 911,
        k_He4_Al27_to_p_Si30 = 912,
        k_n_Si28_to_He4_Mg25 = 913,
        k_p_Si28_to_He4_Al25 = 914,
        k_He4_Si28_to_p_P31 = 915,
        k_He4_Si28_to_C12_Ne20 = 916,
        k_He4_Si28_to_O16_O16 = 917,
        k_n_Si29_to_He4_Mg26 = 918,
        k_p_Si29_to_n_P29 = 919,
        k_p_Si29_to_He4_Al26 = 920,
        k_He4_Si29_to_n_S32 = 921,
        k_He4_Si29_to_p_P32 = 922,
        k_p_Si30_to_n_P30 = 923,
        k_p_Si30_to_He4_Al27 = 924,
        k_He4_Si30_to_n_S33 = 925,
        k_He4_Si30_to_p_P33 = 926,
        k_p_Si31_to_n_P31 = 927,
        k_He4_Si31_to_n_S34 = 928,
        k_p_Si32_to_n_P32 = 929,
        k_He4_Si32_to_n_S35 = 930,
        k_n_P29_to_p_Si29 = 931,
        k_n_P29_to_He4_Al26 = 932,
        k_He4_P29_to_p_S32 = 933,
        k_n_P30_to_p_Si30 = 934,
        k_n_P30_to_He4_Al27 = 935,
        k_He4_P30_to_n_Cl33 = 936,
        k_He4_P30_to_p_S33 = 937,
        k_n_P31_to_p_Si31 = 938,
        k_p_P31_to_He4_Si28 = 939,
        k_p_P31_to_C12_Ne20 = 940,
        k_p_P31_to_O16_O16 = 941,
        k_He4_P31_to_n_Cl34 = 942,
        k_He4_P31_to_p_S34 = 943,
        k_n_P32_to_p_Si32 = 944,
        k_p_P32_to_n_S32 = 945,
        k_p_P32_to_He4_Si29 = 946,
        k_He4_P32_to_n_Cl35 = 947,
        k_He4_P32_to_p_S35 = 948,
        k_p_P33_to_n_S33 = 949,
        k_p_P33_to_He4_Si30 = 950,
        k_He4_P33_to_n_Cl36 = 951,
        k_He4_P33_to_p_S36 = 952,
        k_n_S32_to_p_P32 = 953,
        k_n_S32_to_He4_Si29 = 954,
        k_p_S32_to_He4_P29 = 955,
        k_He4_S32_to_p_Cl35 = 956,
        k_n_S33_to_p_P33 = 957,
        k_n_S33_to_He4_Si30 = 958,
        k_p_S33_to_n_Cl33 = 959,
        k_p_S33_to_He4_P30 = 960,
        k_He4_S33_to_n_Ar36 = 961,
        k_He4_S33_to_p_Cl36 = 962,
        k_n_S34_to_He4_Si31 = 963,
        k_p_S34_to_n_Cl34 = 964,
        k_p_S34_to_He4_P31 = 965,
        k_He4_S34_to_n_Ar37 = 966,
        k_He4_S34_to_p_Cl37 = 967,
        k_n_S35_to_He4_Si32 = 968,
        k_p_S35_to_n_Cl35 = 969,
        k_p_S35_to_He4_P32 = 970,
        k_He4_S35_to_n_Ar38 = 971,
        k_p_S36_to_n_Cl36 = 972,
        k_p_S36_to_He4_P33 = 973,
        k_He4_S36_to_n_Ar39 = 974,
        k_n_Cl33_to_p_S33 = 975,
        k_n_Cl33_to_He4_P30 = 976,
        k_He4_Cl33_to_p_Ar36 = 977,
        k_n_Cl34_to_p_S34 = 978,
        k_n_Cl34_to_He4_P31 = 979,
        k_He4_Cl34_to_n_K37 = 980,
        k_He4_Cl34_to_p_Ar37 = 981,
        k_n_Cl35_to_p_S35 = 982,
        k_n_Cl35_to_He4_P32 = 983,
        k_p_Cl35_to_He4_S32 = 984,
        k_He4_Cl35_to_n_K38 = 985,
        k_He4_Cl35_to_p_Ar38 = 986,
        k_n_Cl36_to_p_S36 = 987,
        k_n_Cl36_to_He4_P33 = 988,
        k_p_Cl36_to_n_Ar36 = 989,
        k_p_Cl36_to_He4_S33 = 990,
        k_He4_Cl36_to_n_K39 = 991,
        k_He4_Cl36_to_p_Ar39 = 992,
        k_p_Cl37_to_n_Ar37 = 993,
        k_p_Cl37_to_He4_S34 = 994,
        k_He4_Cl37_to_n_K40 = 995,
        k_He4_Cl37_to_p_Ar40 = 996,
        k_n_Ar36_to_p_Cl36 = 997,
        k_n_Ar36_to_He4_S33 = 998,
        k_p_Ar36_to_He4_Cl33 = 999,
        k_He4_Ar36_to_p_K39 = 1000,
        k_n_Ar37_to_p_Cl37 = 1001,
        k_n_Ar37_to_He4_S34 = 1002,
        k_p_Ar37_to_n_K37 = 1003,
        k_p_Ar37_to_He4_Cl34 = 1004,
        k_He4_Ar37_to_n_Ca40 = 1005,
        k_He4_Ar37_to_p_K40 = 1006,
        k_n_Ar38_to_He4_S35 = 1007,
        k_p_Ar38_to_n_K38 = 1008,
        k_p_Ar38_to_He4_Cl35 = 1009,
        k_He4_Ar38_to_n_Ca41 = 1010,
        k_He4_Ar38_to_p_K41 = 1011,
        k_n_Ar39_to_He4_S36 = 1012,
        k_p_Ar39_to_n_K39 = 1013,
        k_p_Ar39_to_He4_Cl36 = 1014,
        k_He4_Ar39_to_n_Ca42 = 1015,
        k_p_Ar40_to_n_K40 = 1016,
        k_p_Ar40_to_He4_Cl37 = 1017,
        k_He4_Ar40_to_n_Ca43 = 1018,
        k_n_K37_to_p_Ar37 = 1019,
        k_n_K37_to_He4_Cl34 = 1020,
        k_He4_K37_to_p_Ca40 = 1021,
        k_n_K38_to_p_Ar38 = 1022,
        k_n_K38_to_He4_Cl35 = 1023,
        k_He4_K38_to_p_Ca41 = 1024,
        k_n_K39_to_p_Ar39 = 1025,
        k_n_K39_to_He4_Cl36 = 1026,
        k_p_K39_to_He4_Ar36 = 1027,
        k_He4_K39_to_p_Ca42 = 1028,
        k_n_K40_to_p_Ar40 = 1029,
        k_n_K40_to_He4_Cl37 = 1030,
        k_p_K40_to_n_Ca40 = 1031,
        k_p_K40_to_He4_Ar37 = 1032,
        k_He4_K40_to_n_Sc43 = 1033,
        k_He4_K40_to_p_Ca43 = 1034,
        k_p_K41_to_n_Ca41 = 1035,
        k_p_K41_to_He4_Ar38 = 1036,
        k_He4_K41_to_n_Sc44 = 1037,
        k_He4_K41_to_p_Ca44 = 1038,
        k_n_Ca40_to_p_K40 = 1039,
        k_n_Ca40_to_He4_Ar37 = 1040,
        k_p_Ca40_to_He4_K37 = 1041,
        k_He4_Ca40_to_p_Sc43 = 1042,
        k_n_Ca41_to_p_K41 = 1043,
        k_n_Ca41_to_He4_Ar38 = 1044,
        k_p_Ca41_to_He4_K38 = 1045,
        k_He4_Ca41_to_n_Ti44 = 1046,
        k_He4_Ca41_to_p_Sc44 = 1047,
        k_n_Ca42_to_He4_Ar39 = 1048,
        k_p_Ca42_to_He4_K39 = 1049,
        k_He4_Ca42_to_n_Ti45 = 1050,
        k_He4_Ca42_to_p_Sc45 = 1051,
        k_n_Ca43_to_He4_Ar40 = 1052,
        k_p_Ca43_to_n_Sc43 = 1053,
        k_p_Ca43_to_He4_K40 = 1054,
        k_He4_Ca43_to_n_Ti46 = 1055,
        k_He4_Ca43_to_p_Sc46 = 1056,
        k_p_Ca44_to_n_Sc44 = 1057,
        k_p_Ca44_to_He4_K41 = 1058,
        k_He4_Ca44_to_n_Ti47 = 1059,
        k_He4_Ca44_to_p_Sc47 = 1060,
        k_p_Ca45_to_n_Sc45 = 1061,
        k_He4_Ca45_to_n_Ti48 = 1062,
        k_He4_Ca45_to_p_Sc48 = 1063,
        k_p_Ca46_to_n_Sc46 = 1064,
        k_He4_Ca46_to_n_Ti49 = 1065,
        k_He4_Ca46_to_p_Sc49 = 1066,
        k_p_Ca47_to_n_Sc47 = 1067,
        k_He4_Ca47_to_n_Ti50 = 1068,
        k_p_Ca48_to_n_Sc48 = 1069,
        k_He4_Ca48_to_n_Ti51 = 1070,
        k_n_Sc43_to_p_Ca43 = 1071,
        k_n_Sc43_to_He4_K40 = 1072,
        k_p_Sc43_to_He4_Ca40 = 1073,
        k_He4_Sc43_to_n_V46 = 1074,
        k_He4_Sc43_to_p_Ti46 = 1075,
        k_n_Sc44_to_p_Ca44 = 1076,
        k_n_Sc44_to_He4_K41 = 1077,
        k_p_Sc44_to_n_Ti44 = 1078,
        k_p_Sc44_to_He4_Ca41 = 1079,
        k_He4_Sc44_to_n_V47 = 1080,
        k_He4_Sc44_to_p_Ti47 = 1081,
        k_n_Sc45_to_p_Ca45 = 1082,
        k_p_Sc45_to_n_Ti45 = 1083,
        k_p_Sc45_to_He4_Ca42 = 1084,
        k_He4_Sc45_to_n_V48 = 1085,
        k_He4_Sc45_to_p_Ti48 = 1086,
        k_n_Sc46_to_p_Ca46 = 1087,
        k_p_Sc46_to_n_Ti46 = 1088,
        k_p_Sc46_to_He4_Ca43 = 1089,
        k_He4_Sc46_to_n_V49 = 1090,
        k_He4_Sc46_to_p_Ti49 = 1091,
        k_n_Sc47_to_p_Ca47 = 1092,
        k_p_Sc47_to_n_Ti47 = 1093,
        k_p_Sc47_to_He4_Ca44 = 1094,
        k_He4_Sc47_to_n_V50 = 1095,
        k_He4_Sc47_to_p_Ti50 = 1096,
        k_n_Sc48_to_p_Ca48 = 1097,
        k_p_Sc48_to_n_Ti48 = 1098,
        k_p_Sc48_to_He4_Ca45 = 1099,
        k_He4_Sc48_to_n_V51 = 1100,
        k_He4_Sc48_to_p_Ti51 = 1101,
        k_p_Sc49_to_n_Ti49 = 1102,
        k_p_Sc49_to_He4_Ca46 = 1103,
        k_He4_Sc49_to_n_V52 = 1104,
        k_n_Ti44_to_p_Sc44 = 1105,
        k_n_Ti44_to_He4_Ca41 = 1106,
        k_He4_Ti44_to_p_V47 = 1107,
        k_n_Ti45_to_p_Sc45 = 1108,
        k_n_Ti45_to_He4_Ca42 = 1109,
        k_He4_Ti45_to_n_Cr48 = 1110,
        k_He4_Ti45_to_p_V48 = 1111,
        k_n_Ti46_to_p_Sc46 = 1112,
        k_n_Ti46_to_He4_Ca43 = 1113,
        k_p_Ti46_to_n_V46 = 1114,
        k_p_Ti46_to_He4_Sc43 = 1115,
        k_He4_Ti46_to_n_Cr49 = 1116,
        k_He4_Ti46_to_p_V49 = 1117,
        k_n_Ti47_to_p_Sc47 = 1118,
        k_n_Ti47_to_He4_Ca44 = 1119,
        k_p_Ti47_to_n_V47 = 1120,
        k_p_Ti47_to_He4_Sc44 = 1121,
        k_He4_Ti47_to_n_Cr50 = 1122,
        k_He4_Ti47_to_p_V50 = 1123,
        k_n_Ti48_to_p_Sc48 = 1124,
        k_n_Ti48_to_He4_Ca45 = 1125,
        k_p_Ti48_to_n_V48 = 1126,
        k_p_Ti48_to_He4_Sc45 = 1127,
        k_He4_Ti48_to_n_Cr51 = 1128,
        k_He4_Ti48_to_p_V51 = 1129,
        k_n_Ti49_to_p_Sc49 = 1130,
        k_n_Ti49_to_He4_Ca46 = 1131,
        k_p_Ti49_to_n_V49 = 1132,
        k_p_Ti49_to_He4_Sc46 = 1133,
        k_He4_Ti49_to_n_Cr52 = 1134,
        k_He4_Ti49_to_p_V52 = 1135,
        k_n_Ti50_to_He4_Ca47 = 1136,
        k_p_Ti50_to_n_V50 = 1137,
        k_p_Ti50_to_He4_Sc47 = 1138,
        k_He4_Ti50_to_n_Cr53 = 1139,
        k_n_Ti51_to_He4_Ca48 = 1140,
        k_p_Ti51_to_n_V51 = 1141,
        k_p_Ti51_to_He4_Sc48 = 1142,
        k_He4_Ti51_to_n_Cr54 = 1143,
        k_n_V46_to_p_Ti46 = 1144,
        k_n_V46_to_He4_Sc43 = 1145,
        k_He4_V46_to_p_Cr49 = 1146,
        k_n_V47_to_p_Ti47 = 1147,
        k_n_V47_to_He4_Sc44 = 1148,
        k_p_V47_to_He4_Ti44 = 1149,
        k_He4_V47_to_n_Mn50 = 1150,
        k_He4_V47_to_p_Cr50 = 1151,
        k_n_V48_to_p_Ti48 = 1152,
        k_n_V48_to_He4_Sc45 = 1153,
        k_p_V48_to_n_Cr48 = 1154,
        k_p_V48_to_He4_Ti45 = 1155,
        k_He4_V48_to_n_Mn51 = 1156,
        k_He4_V48_to_p_Cr51 = 1157,
        k_n_V49_to_p_Ti49 = 1158,
        k_n_V49_to_He4_Sc46 = 1159,
        k_p_V49_to_n_Cr49 = 1160,
        k_p_V49_to_He4_Ti46 = 1161,
        k_He4_V49_to_n_Mn52 = 1162,
        k_He4_V49_to_p_Cr52 = 1163,
        k_n_V50_to_p_Ti50 = 1164,
        k_n_V50_to_He4_Sc47 = 1165,
        k_p_V50_to_n_Cr50 = 1166,
        k_p_V50_to_He4_Ti47 = 1167,
        k_He4_V50_to_n_Mn53 = 1168,
        k_He4_V50_to_p_Cr53 = 1169,
        k_n_V51_to_p_Ti51 = 1170,
        k_n_V51_to_He4_Sc48 = 1171,
        k_p_V51_to_n_Cr51 = 1172,
        k_p_V51_to_He4_Ti48 = 1173,
        k_He4_V51_to_n_Mn54 = 1174,
        k_He4_V51_to_p_Cr54 = 1175,
        k_n_V52_to_He4_Sc49 = 1176,
        k_p_V52_to_n_Cr52 = 1177,
        k_p_V52_to_He4_Ti49 = 1178,
        k_He4_V52_to_n_Mn55 = 1179,
        k_n_Cr48_to_p_V48 = 1180,
        k_n_Cr48_to_He4_Ti45 = 1181,
        k_He4_Cr48_to_p_Mn51 = 1182,
        k_n_Cr49_to_p_V49 = 1183,
        k_n_Cr49_to_He4_Ti46 = 1184,
        k_p_Cr49_to_He4_V46 = 1185,
        k_He4_Cr49_to_n_Fe52 = 1186,
        k_He4_Cr49_to_p_Mn52 = 1187,
        k_n_Cr50_to_p_V50 = 1188,
        k_n_Cr50_to_He4_Ti47 = 1189,
        k_p_Cr50_to_n_Mn50 = 1190,
        k_p_Cr50_to_He4_V47 = 1191,
        k_He4_Cr50_to_n_Fe53 = 1192,
        k_He4_Cr50_to_p_Mn53 = 1193,
        k_n_Cr51_to_p_V51 = 1194,
        k_n_Cr51_to_He4_Ti48 = 1195,
        k_p_Cr51_to_n_Mn51 = 1196,
        k_p_Cr51_to_He4_V48 = 1197,
        k_He4_Cr51_to_n_Fe54 = 1198,
        k_He4_Cr51_to_p_Mn54 = 1199,
        k_n_Cr52_to_p_V52 = 1200,
        k_n_Cr52_to_He4_Ti49 = 1201,
        k_p_Cr52_to_n_Mn52 = 1202,
        k_p_Cr52_to_He4_V49 = 1203,
        k_He4_Cr52_to_n_Fe55 = 1204,
        k_He4_Cr52_to_p_Mn55 = 1205,
        k_n_Cr53_to_He4_Ti50 = 1206,
        k_p_Cr53_to_n_Mn53 = 1207,
        k_p_Cr53_to_He4_V50 = 1208,
        k_He4_Cr53_to_n_Fe56 = 1209,
        k_n_Cr54_to_He4_Ti51 = 1210,
        k_p_Cr54_to_n_Mn54 = 1211,
        k_p_Cr54_to_He4_V51 = 1212,
        k_He4_Cr54_to_n_Fe57 = 1213,
        k_n_Mn50_to_p_Cr50 = 1214,
        k_n_Mn50_to_He4_V47 = 1215,
        k_He4_Mn50_to_n_Co53 = 1216,
        k_He4_Mn50_to_p_Fe53 = 1217,
        k_n_Mn51_to_p_Cr51 = 1218,
        k_n_Mn51_to_He4_V48 = 1219,
        k_p_Mn51_to_He4_Cr48 = 1220,
        k_He4_Mn51_to_n_Co54 = 1221,
        k_He4_Mn51_to_p_Fe54 = 1222,
        k_n_Mn52_to_p_Cr52 = 1223,
        k_n_Mn52_to_He4_V49 = 1224,
        k_p_Mn52_to_n_Fe52 = 1225,
        k_p_Mn52_to_He4_Cr49 = 1226,
        k_He4_Mn52_to_n_Co55 = 1227,
        k_He4_Mn52_to_p_Fe55 = 1228,
        k_n_Mn53_to_p_Cr53 = 1229,
        k_n_Mn53_to_He4_V50 = 1230,
        k_p_Mn53_to_n_Fe53 = 1231,
        k_p_Mn53_to_He4_Cr50 = 1232,
        k_He4_Mn53_to_n_Co56 = 1233,
        k_He4_Mn53_to_p_Fe56 = 1234,
        k_n_Mn54_to_p_Cr54 = 1235,
        k_n_Mn54_to_He4_V51 = 1236,
        k_p_Mn54_to_n_Fe54 = 1237,
        k_p_Mn54_to_He4_Cr51 = 1238,
        k_He4_Mn54_to_n_Co57 = 1239,
        k_He4_Mn54_to_p_Fe57 = 1240,
        k_n_Mn55_to_He4_V52 = 1241,
        k_p_Mn55_to_n_Fe55 = 1242,
        k_p_Mn55_to_He4_Cr52 = 1243,
        k_He4_Mn55_to_n_Co58 = 1244,
        k_He4_Mn55_to_p_Fe58 = 1245,
        k_n_Fe52_to_p_Mn52 = 1246,
        k_n_Fe52_to_He4_Cr49 = 1247,
        k_He4_Fe52_to_p_Co55 = 1248,
        k_n_Fe53_to_p_Mn53 = 1249,
        k_n_Fe53_to_He4_Cr50 = 1250,
        k_p_Fe53_to_n_Co53 = 1251,
        k_p_Fe53_to_He4_Mn50 = 1252,
        k_He4_Fe53_to_n_Ni56 = 1253,
        k_He4_Fe53_to_p_Co56 = 1254,
        k_n_Fe54_to_p_Mn54 = 1255,
        k_n_Fe54_to_He4_Cr51 = 1256,
        k_p_Fe54_to_n_Co54 = 1257,
        k_p_Fe54_to_He4_Mn51 = 1258,
        k_He4_Fe54_to_n_Ni57 = 1259,
        k_He4_Fe54_to_p_Co57 = 1260,
        k_n_Fe55_to_p_Mn55 = 1261,
        k_n_Fe55_to_He4_Cr52 = 1262,
        k_p_Fe55_to_n_Co55 = 1263,
        k_p_Fe55_to_He4_Mn52 = 1264,
        k_He4_Fe55_to_n_Ni58 = 1265,
        k_He4_Fe55_to_p_Co58 = 1266,
        k_n_Fe56_to_He4_Cr53 = 1267,
        k_p_Fe56_to_n_Co56 = 1268,
        k_p_Fe56_to_He4_Mn53 = 1269,
        k_He4_Fe56_to_n_Ni59 = 1270,
        k_He4_Fe56_to_p_Co59 = 1271,
        k_n_Fe57_to_He4_Cr54 = 1272,
        k_p_Fe57_to_n_Co57 = 1273,
        k_p_Fe57_to_He4_Mn54 = 1274,
        k_He4_Fe57_to_n_Ni60 = 1275,
        k_p_Fe58_to_n_Co58 = 1276,
        k_p_Fe58_to_He4_Mn55 = 1277,
        k_He4_Fe58_to_n_Ni61 = 1278,
        k_n_Co53_to_p_Fe53 = 1279,
        k_n_Co53_to_He4_Mn50 = 1280,
        k_He4_Co53_to_p_Ni56 = 1281,
        k_n_Co54_to_p_Fe54 = 1282,
        k_n_Co54_to_He4_Mn51 = 1283,
        k_He4_Co54_to_n_Cu57 = 1284,
        k_He4_Co54_to_p_Ni57 = 1285,
        k_n_Co55_to_p_Fe55 = 1286,
        k_n_Co55_to_He4_Mn52 = 1287,
        k_p_Co55_to_He4_Fe52 = 1288,
        k_He4_Co55_to_n_Cu58 = 1289,
        k_He4_Co55_to_p_Ni58 = 1290,
        k_n_Co56_to_p_Fe56 = 1291,
        k_n_Co56_to_He4_Mn53 = 1292,
        k_p_Co56_to_n_Ni56 = 1293,
        k_p_Co56_to_He4_Fe53 = 1294,
        k_He4_Co56_to_n_Cu59 = 1295,
        k_He4_Co56_to_p_Ni59 = 1296,
        k_n_Co57_to_p_Fe57 = 1297,
        k_n_Co57_to_He4_Mn54 = 1298,
        k_p_Co57_to_n_Ni57 = 1299,
        k_p_Co57_to_He4_Fe54 = 1300,
        k_He4_Co57_to_n_Cu60 = 1301,
        k_He4_Co57_to_p_Ni60 = 1302,
        k_n_Co58_to_p_Fe58 = 1303,
        k_n_Co58_to_He4_Mn55 = 1304,
        k_p_Co58_to_n_Ni58 = 1305,
        k_p_Co58_to_He4_Fe55 = 1306,
        k_He4_Co58_to_n_Cu61 = 1307,
        k_He4_Co58_to_p_Ni61 = 1308,
        k_p_Co59_to_n_Ni59 = 1309,
        k_p_Co59_to_He4_Fe56 = 1310,
        k_He4_Co59_to_n_Cu62 = 1311,
        k_He4_Co59_to_p_Ni62 = 1312,
        k_n_Ni56_to_p_Co56 = 1313,
        k_n_Ni56_to_He4_Fe53 = 1314,
        k_p_Ni56_to_He4_Co53 = 1315,
        k_He4_Ni56_to_n_Zn59 = 1316,
        k_He4_Ni56_to_p_Cu59 = 1317,
        k_n_Ni57_to_p_Co57 = 1318,
        k_n_Ni57_to_He4_Fe54 = 1319,
        k_p_Ni57_to_n_Cu57 = 1320,
        k_p_Ni57_to_He4_Co54 = 1321,
        k_He4_Ni57_to_n_Zn60 = 1322,
        k_He4_Ni57_to_p_Cu60 = 1323,
        k_n_Ni58_to_p_Co58 = 1324,
        k_n_Ni58_to_He4_Fe55 = 1325,
        k_p_Ni58_to_n_Cu58 = 1326,
        k_p_Ni58_to_He4_Co55 = 1327,
        k_He4_Ni58_to_n_Zn61 = 1328,
        k_He4_Ni58_to_p_Cu61 = 1329,
        k_n_Ni59_to_p_Co59 = 1330,
        k_n_Ni59_to_He4_Fe56 = 1331,
        k_p_Ni59_to_n_Cu59 = 1332,
        k_p_Ni59_to_He4_Co56 = 1333,
        k_He4_Ni59_to_n_Zn62 = 1334,
        k_He4_Ni59_to_p_Cu62 = 1335,
        k_n_Ni60_to_He4_Fe57 = 1336,
        k_p_Ni60_to_n_Cu60 = 1337,
        k_p_Ni60_to_He4_Co57 = 1338,
        k_He4_Ni60_to_n_Zn63 = 1339,
        k_He4_Ni60_to_p_Cu63 = 1340,
        k_n_Ni61_to_He4_Fe58 = 1341,
        k_p_Ni61_to_n_Cu61 = 1342,
        k_p_Ni61_to_He4_Co58 = 1343,
        k_He4_Ni61_to_n_Zn64 = 1344,
        k_He4_Ni61_to_p_Cu64 = 1345,
        k_p_Ni62_to_n_Cu62 = 1346,
        k_p_Ni62_to_He4_Co59 = 1347,
        k_He4_Ni62_to_n_Zn65 = 1348,
        k_He4_Ni62_to_p_Cu65 = 1349,
        k_p_Ni63_to_n_Cu63 = 1350,
        k_He4_Ni63_to_n_Zn66 = 1351,
        k_p_Ni64_to_n_Cu64 = 1352,
        k_n_Cu57_to_p_Ni57 = 1353,
        k_n_Cu57_to_He4_Co54 = 1354,
        k_He4_Cu57_to_p_Zn60 = 1355,
        k_n_Cu58_to_p_Ni58 = 1356,
        k_n_Cu58_to_He4_Co55 = 1357,
        k_He4_Cu58_to_p_Zn61 = 1358,
        k_n_Cu59_to_p_Ni59 = 1359,
        k_n_Cu59_to_He4_Co56 = 1360,
        k_p_Cu59_to_n_Zn59 = 1361,
        k_p_Cu59_to_He4_Ni56 = 1362,
        k_He4_Cu59_to_n_Ga62 = 1363,
        k_He4_Cu59_to_p_Zn62 = 1364,
        k_n_Cu60_to_p_Ni60 = 1365,
        k_n_Cu60_to_He4_Co57 = 1366,
        k_p_Cu60_to_n_Zn60 = 1367,
        k_p_Cu60_to_He4_Ni57 = 1368,
        k_He4_Cu60_to_n_Ga63 = 1369,
        k_He4_Cu60_to_p_Zn63 = 1370,
        k_n_Cu61_to_p_Ni61 = 1371,
        k_n_Cu61_to_He4_Co58 = 1372,
        k_p_Cu61_to_n_Zn61 = 1373,
        k_p_Cu61_to_He4_Ni58 = 1374,
        k_He4_Cu61_to_n_Ga64 = 1375,
        k_He4_Cu61_to_p_Zn64 = 1376,
        k_n_Cu62_to_p_Ni62 = 1377,
        k_n_Cu62_to_He4_Co59 = 1378,
        k_p_Cu62_to_n_Zn62 = 1379,
        k_p_Cu62_to_He4_Ni59 = 1380,
        k_He4_Cu62_to_p_Zn65 = 1381,
        k_n_Cu63_to_p_Ni63 = 1382,
        k_p_Cu63_to_n_Zn63 = 1383,
        k_p_Cu63_to_He4_Ni60 = 1384,
        k_He4_Cu63_to_p_Zn66 = 1385,
        k_n_Cu64_to_p_Ni64 = 1386,
        k_p_Cu64_to_n_Zn64 = 1387,
        k_p_Cu64_to_He4_Ni61 = 1388,
        k_p_Cu65_to_n_Zn65 = 1389,
        k_p_Cu65_to_He4_Ni62 = 1390,
        k_n_Zn59_to_p_Cu59 = 1391,
        k_n_Zn59_to_He4_Ni56 = 1392,
        k_He4_Zn59_to_p_Ga62 = 1393,
        k_n_Zn60_to_p_Cu60 = 1394,
        k_n_Zn60_to_He4_Ni57 = 1395,
        k_p_Zn60_to_He4_Cu57 = 1396,
        k_He4_Zn60_to_n_Ge63 = 1397,
        k_He4_Zn60_to_p_Ga63 = 1398,
        k_n_Zn61_to_p_Cu61 = 1399,
        k_n_Zn61_to_He4_Ni58 = 1400,
        k_p_Zn61_to_He4_Cu58 = 1401,
        k_He4_Zn61_to_n_Ge64 = 1402,
        k_He4_Zn61_to_p_Ga64 = 1403,
        k_n_Zn62_to_p_Cu62 = 1404,
        k_n_Zn62_to_He4_Ni59 = 1405,
        k_p_Zn62_to_n_Ga62 = 1406,
        k_p_Zn62_to_He4_Cu59 = 1407,
        k_n_Zn63_to_p_Cu63 = 1408,
        k_n_Zn63_to_He4_Ni60 = 1409,
        k_p_Zn63_to_n_Ga63 = 1410,
        k_p_Zn63_to_He4_Cu60 = 1411,
        k_n_Zn64_to_p_Cu64 = 1412,
        k_n_Zn64_to_He4_Ni61 = 1413,
        k_p_Zn64_to_n_Ga64 = 1414,
        k_p_Zn64_to_He4_Cu61 = 1415,
        k_n_Zn65_to_p_Cu65 = 1416,
        k_n_Zn65_to_He4_Ni62 = 1417,
        k_p_Zn65_to_He4_Cu62 = 1418,
        k_n_Zn66_to_He4_Ni63 = 1419,
        k_p_Zn66_to_He4_Cu63 = 1420,
        k_n_Ga62_to_p_Zn62 = 1421,
        k_n_Ga62_to_He4_Cu59 = 1422,
        k_p_Ga62_to_He4_Zn59 = 1423,
        k_n_Ga63_to_p_Zn63 = 1424,
        k_n_Ga63_to_He4_Cu60 = 1425,
        k_p_Ga63_to_n_Ge63 = 1426,
        k_p_Ga63_to_He4_Zn60 = 1427,
        k_n_Ga64_to_p_Zn64 = 1428,
        k_n_Ga64_to_He4_Cu61 = 1429,
        k_p_Ga64_to_n_Ge64 = 1430,
        k_p_Ga64_to_He4_Zn61 = 1431,
        k_n_Ge63_to_p_Ga63 = 1432,
        k_n_Ge63_to_He4_Zn60 = 1433,
        k_n_Ge64_to_p_Ga64 = 1434,
        k_n_Ge64_to_He4_Zn61 = 1435,
        k_p_d_to_n_p_p = 1436,
        k_He3_He3_to_p_p_He4 = 1437,
        k_d_Li7_to_n_He4_He4 = 1438,
        k_d_Be7_to_p_He4_He4 = 1439,
        k_p_Be9_to_d_He4_He4 = 1440,
        k_n_B8_to_p_He4_He4 = 1441,
        k_p_B11_to_He4_He4_He4 = 1442,
        k_He3_Li7_to_n_p_He4_He4 = 1443,
        k_He3_Be7_to_p_p_He4_He4 = 1444,
        k_p_Be9_to_n_p_He4_He4 = 1445,
        k_n_p_He4_to_Li6 = 1446,
        k_n_He4_He4_to_Be9 = 1447,
        k_He4_He4_He4_to_C12 = 1448,
        k_n_p_p_to_p_d = 1449,
        k_p_p_He4_to_He3_He3 = 1450,
        k_n_He4_He4_to_d_Li7 = 1451,
        k_p_He4_He4_to_n_B8 = 1452,
        k_p_He4_He4_to_d_Be7 = 1453,
        k_d_He4_He4_to_p_Be9 = 1454,
        k_He4_He4_He4_to_p_B11 = 1455,
        k_n_p_He4_He4_to_He3_Li7 = 1456,
        k_n_p_He4_He4_to_p_Be9 = 1457,
        k_p_p_He4_He4_to_He3_Be7 = 1458,
        k_Al25_to_Mg25 = 1459,
        k_Al26_to_Mg26 = 1460,
        k_Ar36_to_Cl36 = 1461,
        k_Ar37_to_Cl37 = 1462,
        k_Ar37_to_K37 = 1463,
        k_Ar38_to_K38 = 1464,
        k_Ar39_to_K39 = 1465,
        k_Ar40_to_K40 = 1466,
        k_Ca40_to_K40 = 1467,
        k_Ca41_to_K41 = 1468,
        k_Ca43_to_Sc43 = 1469,
        k_Ca44_to_Sc44 = 1470,
        k_Cl33_to_S33 = 1471,
        k_Cl34_to_S34 = 1472,
        k_Cl35_to_S35 = 1473,
        k_Cl36_to_Ar36 = 1474,
        k_Cl36_to_S36 = 1475,
        k_Cl37_to_Ar37 = 1476,
        k_Co53_to_Fe53 = 1477,
        k_Cu57_to_Ni57 = 1478,
        k_Fe53_to_Co53 = 1479,
        k_K37_to_Ar37 = 1480,
        k_K38_to_Ar38 = 1481,
        k_K39_to_Ar39 = 1482,
        k_K40_to_Ar40 = 1483,
        k_K40_to_Ca40 = 1484,
        k_K41_to_Ca41 = 1485,
        k_Mg23_to_Na23 = 1486,
        k_Mg25_to_Al25 = 1487,
        k_Mg26_to_Al26 = 1488,
        k_Na21_to_Ne21 = 1489,
        k_Na22_to_Ne22 = 1490,
        k_Na23_to_Mg23 = 1491,
        k_Ne21_to_Na21 = 1492,
        k_Ne22_to_Na22 = 1493,
        k_Ni57_to_Cu57 = 1494,
        k_P29_to_Si29 = 1495,
        k_P30_to_Si30 = 1496,
        k_P31_to_Si31 = 1497,
        k_P32_to_S32 = 1498,
        k_P32_to_Si32 = 1499,
        k_P33_to_S33 = 1500,
        k_S32_to_P32 = 1501,
        k_S33_to_Cl33 = 1502,
        k_S33_to_P33 = 1503,
        k_S34_to_Cl34 = 1504,
        k_S35_to_Cl35 = 1505,
        k_S36_to_Cl36 = 1506,
        k_Sc43_to_Ca43 = 1507,
        k_Sc44_to_Ca44 = 1508,
        k_Sc44_to_Ti44 = 1509,
        k_Si29_to_P29 = 1510,
        k_Si30_to_P30 = 1511,
        k_Si31_to_P31 = 1512,
        k_Si32_to_P32 = 1513,
        k_Ti44_to_Sc44 = 1514,
        k_Ca45_to_Sc45 = 1515,
        k_Ca46_to_Sc46 = 1516,
        k_Ca47_to_Sc47 = 1517,
        k_Ca48_to_Sc48 = 1518,
        k_Co54_to_Fe54 = 1519,
        k_Co55_to_Fe55 = 1520,
        k_Co56_to_Fe56 = 1521,
        k_Co56_to_Ni56 = 1522,
        k_Co57_to_Fe57 = 1523,
        k_Co57_to_Ni57 = 1524,
        k_Co58_to_Fe58 = 1525,
        k_Co58_to_Ni58 = 1526,
        k_Co59_to_Ni59 = 1527,
        k_Cr48_to_V48 = 1528,
        k_Cr49_to_V49 = 1529,
        k_Cr50_to_Mn50 = 1530,
        k_Cr50_to_V50 = 1531,
        k_Cr51_to_Mn51 = 1532,
        k_Cr51_to_V51 = 1533,
        k_Cr52_to_Mn52 = 1534,
        k_Cr52_to_V52 = 1535,
        k_Cr53_to_Mn53 = 1536,
        k_Cr54_to_Mn54 = 1537,
        k_Cu58_to_Ni58 = 1538,
        k_Cu59_to_Ni59 = 1539,
        k_Cu60_to_Ni60 = 1540,
        k_Cu60_to_Zn60 = 1541,
        k_Cu61_to_Ni61 = 1542,
        k_Cu61_to_Zn61 = 1543,
        k_Cu62_to_Ni62 = 1544,
        k_Cu62_to_Zn62 = 1545,
        k_Cu63_to_Ni63 = 1546,
        k_Cu63_to_Zn63 = 1547,
        k_Cu64_to_Ni64 = 1548,
        k_Cu64_to_Zn64 = 1549,
        k_Cu65_to_Zn65 = 1550,
        k_Fe52_to_Mn52 = 1551,
        k_Fe53_to_Mn53 = 1552,
        k_Fe54_to_Co54 = 1553,
        k_Fe54_to_Mn54 = 1554,
        k_Fe55_to_Co55 = 1555,
        k_Fe55_to_Mn55 = 1556,
        k_Fe56_to_Co56 = 1557,
        k_Fe57_to_Co57 = 1558,
        k_Fe58_to_Co58 = 1559,
        k_Ga62_to_Zn62 = 1560,
        k_Ga63_to_Zn63 = 1561,
        k_Ga64_to_Ge64 = 1562,
        k_Ga64_to_Zn64 = 1563,
        k_Ge64_to_Ga64 = 1564,
        k_Mn50_to_Cr50 = 1565,
        k_Mn51_to_Cr51 = 1566,
        k_Mn52_to_Cr52 = 1567,
        k_Mn52_to_Fe52 = 1568,
        k_Mn53_to_Cr53 = 1569,
        k_Mn53_to_Fe53 = 1570,
        k_Mn54_to_Cr54 = 1571,
        k_Mn54_to_Fe54 = 1572,
        k_Mn55_to_Fe55 = 1573,
        k_n_to_p = 1574,
        k_Ni56_to_Co56 = 1575,
        k_Ni57_to_Co57 = 1576,
        k_Ni58_to_Co58 = 1577,
        k_Ni58_to_Cu58 = 1578,
        k_Ni59_to_Co59 = 1579,
        k_Ni59_to_Cu59 = 1580,
        k_Ni60_to_Cu60 = 1581,
        k_Ni61_to_Cu61 = 1582,
        k_Ni62_to_Cu62 = 1583,
        k_Ni63_to_Cu63 = 1584,
        k_Ni64_to_Cu64 = 1585,
        k_p_to_n = 1586,
        k_Sc45_to_Ca45 = 1587,
        k_Sc45_to_Ti45 = 1588,
        k_Sc46_to_Ca46 = 1589,
        k_Sc46_to_Ti46 = 1590,
        k_Sc47_to_Ca47 = 1591,
        k_Sc47_to_Ti47 = 1592,
        k_Sc48_to_Ca48 = 1593,
        k_Sc48_to_Ti48 = 1594,
        k_Sc49_to_Ti49 = 1595,
        k_Ti45_to_Sc45 = 1596,
        k_Ti46_to_Sc46 = 1597,
        k_Ti46_to_V46 = 1598,
        k_Ti47_to_Sc47 = 1599,
        k_Ti47_to_V47 = 1600,
        k_Ti48_to_Sc48 = 1601,
        k_Ti48_to_V48 = 1602,
        k_Ti49_to_Sc49 = 1603,
        k_Ti49_to_V49 = 1604,
        k_Ti50_to_V50 = 1605,
        k_Ti51_to_V51 = 1606,
        k_V46_to_Ti46 = 1607,
        k_V47_to_Ti47 = 1608,
        k_V48_to_Cr48 = 1609,
        k_V48_to_Ti48 = 1610,
        k_V49_to_Cr49 = 1611,
        k_V49_to_Ti49 = 1612,
        k_V50_to_Cr50 = 1613,
        k_V50_to_Ti50 = 1614,
        k_V51_to_Cr51 = 1615,
        k_V51_to_Ti51 = 1616,
        k_V52_to_Cr52 = 1617,
        k_Zn60_to_Cu60 = 1618,
        k_Zn61_to_Cu61 = 1619,
        k_Zn62_to_Cu62 = 1620,
        k_Zn62_to_Ga62 = 1621,
        k_Zn63_to_Cu63 = 1622,
        k_Zn63_to_Ga63 = 1623,
        k_Zn64_to_Cu64 = 1624,
        k_Zn64_to_Ga64 = 1625,
        k_Zn65_to_Cu65 = 1626,
        NumRates = k_Zn65_to_Cu65
    };

    // number of reaclib rates

    const int NrateReaclib = 1458;

    // number of tabular rates

    const int NrateTabular = 168;

    // rate names -- note: the rates are 1-based, not zero-based, so we pad
    // this vector with rate_names[0] = "" so the indices line up with the
    // NetworkRates enum

    static const std::vector<std::string> rate_names = {
        "",   // 0
        "Be7_to_Li7_weak_electron_capture",  // 1,
        "C14_to_N14_weak_wc12",  // 2,
        "N13_to_C13_weak_wc12",  // 3,
        "O14_to_N14_weak_wc12",  // 4,
        "O15_to_N15_weak_wc12",  // 5,
        "F17_to_O17_weak_wc12",  // 6,
        "F18_to_O18_weak_wc12",  // 7,
        "Ne18_to_F18_weak_wc12",  // 8,
        "Ne19_to_F19_weak_wc12",  // 9,
        "Zn59_to_Cu59_weak_wc12",  // 10,
        "Ge63_to_Ga63_weak_wc12",  // 11,
        "d_to_n_p",  // 12,
        "He3_to_p_d",  // 13,
        "He4_to_n_He3",  // 14,
        "He4_to_d_d",  // 15,
        "Li6_to_He4_d",  // 16,
        "Li7_to_n_Li6",  // 17,
        "Be7_to_p_Li6",  // 18,
        "Be7_to_He4_He3",  // 19,
        "B8_to_p_Be7",  // 20,
        "B8_to_He4_He4_weak_wc12",  // 21,
        "B10_to_p_Be9",  // 22,
        "B10_to_He4_Li6",  // 23,
        "B11_to_n_B10",  // 24,
        "B11_to_He4_Li7",  // 25,
        "C12_to_p_B11",  // 26,
        "C13_to_n_C12",  // 27,
        "C14_to_n_C13",  // 28,
        "N13_to_p_C12",  // 29,
        "N14_to_n_N13",  // 30,
        "N14_to_p_C13",  // 31,
        "N15_to_n_N14",  // 32,
        "N15_to_p_C14",  // 33,
        "O14_to_p_N13",  // 34,
        "O15_to_n_O14",  // 35,
        "O15_to_p_N14",  // 36,
        "O16_to_n_O15",  // 37,
        "O16_to_p_N15",  // 38,
        "O16_to_He4_C12",  // 39,
        "O17_to_n_O16",  // 40,
        "O18_to_n_O17",  // 41,
        "O18_to_He4_C14",  // 42,
        "F17_to_p_O16",  // 43,
        "F18_to_n_F17",  // 44,
        "F18_to_p_O17",  // 45,
        "F18_to_He4_N14",  // 46,
        "F19_to_n_F18",  // 47,
        "F19_to_p_O18",  // 48,
        "F19_to_He4_N15",  // 49,
        "Ne18_to_p_F17",  // 50,
        "Ne18_to_He4_O14",  // 51,
        "Ne19_to_n_Ne18",  // 52,
        "Ne19_to_p_F18",  // 53,
        "Ne19_to_He4_O15",  // 54,
        "Ne20_to_n_Ne19",  // 55,
        "Ne20_to_p_F19",  // 56,
        "Ne20_to_He4_O16",  // 57,
        "Ne21_to_n_Ne20",  // 58,
        "Ne21_to_He4_O17",  // 59,
        "Ne22_to_n_Ne21",  // 60,
        "Ne22_to_He4_O18",  // 61,
        "Na21_to_p_Ne20",  // 62,
        "Na21_to_He4_F17",  // 63,
        "Na22_to_n_Na21",  // 64,
        "Na22_to_p_Ne21",  // 65,
        "Na22_to_He4_F18",  // 66,
        "Na23_to_n_Na22",  // 67,
        "Na23_to_p_Ne22",  // 68,
        "Na23_to_He4_F19",  // 69,
        "Mg23_to_p_Na22",  // 70,
        "Mg23_to_He4_Ne19",  // 71,
        "Mg24_to_n_Mg23",  // 72,
        "Mg24_to_p_Na23",  // 73,
        "Mg24_to_He4_Ne20",  // 74,
        "Mg25_to_n_Mg24",  // 75,
        "Mg25_to_He4_Ne21",  // 76,
        "Mg26_to_n_Mg25",  // 77,
        "Mg26_to_He4_Ne22",  // 78,
        "Al25_to_p_Mg24",  // 79,
        "Al25_to_He4_Na21",  // 80,
        "Al26_to_n_Al25",  // 81,
        "Al26_to_p_Mg25",  // 82,
        "Al26_to_He4_Na22",  // 83,
        "Al27_to_n_Al26",  // 84,
        "Al27_to_p_Mg26",  // 85,
        "Al27_to_He4_Na23",  // 86,
        "Si28_to_p_Al27",  // 87,
        "Si28_to_He4_Mg24",  // 88,
        "Si29_to_n_Si28",  // 89,
        "Si29_to_He4_Mg25",  // 90,
        "Si30_to_n_Si29",  // 91,
        "Si30_to_He4_Mg26",  // 92,
        "Si31_to_n_Si30",  // 93,
        "Si32_to_n_Si31",  // 94,
        "P29_to_p_Si28",  // 95,
        "P29_to_He4_Al25",  // 96,
        "P30_to_n_P29",  // 97,
        "P30_to_p_Si29",  // 98,
        "P30_to_He4_Al26",  // 99,
        "P31_to_n_P30",  // 100,
        "P31_to_p_Si30",  // 101,
        "P31_to_He4_Al27",  // 102,
        "P32_to_n_P31",  // 103,
        "P32_to_p_Si31",  // 104,
        "P33_to_n_P32",  // 105,
        "P33_to_p_Si32",  // 106,
        "S32_to_p_P31",  // 107,
        "S32_to_He4_Si28",  // 108,
        "S33_to_n_S32",  // 109,
        "S33_to_p_P32",  // 110,
        "S33_to_He4_Si29",  // 111,
        "S34_to_n_S33",  // 112,
        "S34_to_p_P33",  // 113,
        "S34_to_He4_Si30",  // 114,
        "S35_to_n_S34",  // 115,
        "S35_to_He4_Si31",  // 116,
        "S36_to_n_S35",  // 117,
        "S36_to_He4_Si32",  // 118,
        "Cl33_to_p_S32",  // 119,
        "Cl33_to_He4_P29",  // 120,
        "Cl34_to_n_Cl33",  // 121,
        "Cl34_to_p_S33",  // 122,
        "Cl34_to_He4_P30",  // 123,
        "Cl35_to_n_Cl34",  // 124,
        "Cl35_to_p_S34",  // 125,
        "Cl35_to_He4_P31",  // 126,
        "Cl36_to_n_Cl35",  // 127,
        "Cl36_to_p_S35",  // 128,
        "Cl36_to_He4_P32",  // 129,
        "Cl37_to_n_Cl36",  // 130,
        "Cl37_to_p_S36",  // 131,
        "Cl37_to_He4_P33",  // 132,
        "Ar36_to_p_Cl35",  // 133,
        "Ar36_to_He4_S32",  // 134,
        "Ar37_to_n_Ar36",  // 135,
        "Ar37_to_p_Cl36",  // 136,
        "Ar37_to_He4_S33",  // 137,
        "Ar38_to_n_Ar37",  // 138,
        "Ar38_to_p_Cl37",  // 139,
        "Ar38_to_He4_S34",  // 140,
        "Ar39_to_n_Ar38",  // 141,
        "Ar39_to_He4_S35",  // 142,
        "Ar40_to_n_Ar39",  // 143,
        "Ar40_to_He4_S36",  // 144,
        "K37_to_p_Ar36",  // 145,
        "K37_to_He4_Cl33",  // 146,
        "K38_to_n_K37",  // 147,
        "K38_to_p_Ar37",  // 148,
        "K38_to_He4_Cl34",  // 149,
        "K39_to_n_K38",  // 150,
        "K39_to_p_Ar38",  // 151,
        "K39_to_He4_Cl35",  // 152,
        "K40_to_n_K39",  // 153,
        "K40_to_p_Ar39",  // 154,
        "K40_to_He4_Cl36",  // 155,
        "K41_to_n_K40",  // 156,
        "K41_to_p_Ar40",  // 157,
        "K41_to_He4_Cl37",  // 158,
        "Ca40_to_p_K39",  // 159,
        "Ca40_to_He4_Ar36",  // 160,
        "Ca41_to_n_Ca40",  // 161,
        "Ca41_to_p_K40",  // 162,
        "Ca41_to_He4_Ar37",  // 163,
        "Ca42_to_n_Ca41",  // 164,
        "Ca42_to_p_K41",  // 165,
        "Ca42_to_He4_Ar38",  // 166,
        "Ca43_to_n_Ca42",  // 167,
        "Ca43_to_He4_Ar39",  // 168,
        "Ca44_to_n_Ca43",  // 169,
        "Ca44_to_He4_Ar40",  // 170,
        "Ca45_to_n_Ca44",  // 171,
        "Ca46_to_n_Ca45",  // 172,
        "Ca47_to_n_Ca46",  // 173,
        "Ca48_to_n_Ca47",  // 174,
        "Sc43_to_p_Ca42",  // 175,
        "Sc43_to_He4_K39",  // 176,
        "Sc44_to_n_Sc43",  // 177,
        "Sc44_to_p_Ca43",  // 178,
        "Sc44_to_He4_K40",  // 179,
        "Sc45_to_n_Sc44",  // 180,
        "Sc45_to_p_Ca44",  // 181,
        "Sc45_to_He4_K41",  // 182,
        "Sc46_to_n_Sc45",  // 183,
        "Sc46_to_p_Ca45",  // 184,
        "Sc47_to_n_Sc46",  // 185,
        "Sc47_to_p_Ca46",  // 186,
        "Sc48_to_n_Sc47",  // 187,
        "Sc48_to_p_Ca47",  // 188,
        "Sc49_to_n_Sc48",  // 189,
        "Sc49_to_p_Ca48",  // 190,
        "Ti44_to_p_Sc43",  // 191,
        "Ti44_to_He4_Ca40",  // 192,
        "Ti45_to_n_Ti44",  // 193,
        "Ti45_to_p_Sc44",  // 194,
        "Ti45_to_He4_Ca41",  // 195,
        "Ti46_to_n_Ti45",  // 196,
        "Ti46_to_p_Sc45",  // 197,
        "Ti46_to_He4_Ca42",  // 198,
        "Ti47_to_n_Ti46",  // 199,
        "Ti47_to_p_Sc46",  // 200,
        "Ti47_to_He4_Ca43",  // 201,
        "Ti48_to_n_Ti47",  // 202,
        "Ti48_to_p_Sc47",  // 203,
        "Ti48_to_He4_Ca44",  // 204,
        "Ti49_to_n_Ti48",  // 205,
        "Ti49_to_p_Sc48",  // 206,
        "Ti49_to_He4_Ca45",  // 207,
        "Ti50_to_n_Ti49",  // 208,
        "Ti50_to_p_Sc49",  // 209,
        "Ti50_to_He4_Ca46",  // 210,
        "Ti51_to_n_Ti50",  // 211,
        "Ti51_to_He4_Ca47",  // 212,
        "V46_to_p_Ti45",  // 213,
        "V47_to_n_V46",  // 214,
        "V47_to_p_Ti46",  // 215,
        "V47_to_He4_Sc43",  // 216,
        "V48_to_n_V47",  // 217,
        "V48_to_p_Ti47",  // 218,
        "V48_to_He4_Sc44",  // 219,
        "V49_to_n_V48",  // 220,
        "V49_to_p_Ti48",  // 221,
        "V49_to_He4_Sc45",  // 222,
        "V50_to_n_V49",  // 223,
        "V50_to_p_Ti49",  // 224,
        "V50_to_He4_Sc46",  // 225,
        "V51_to_n_V50",  // 226,
        "V51_to_p_Ti50",  // 227,
        "V51_to_He4_Sc47",  // 228,
        "V52_to_n_V51",  // 229,
        "V52_to_p_Ti51",  // 230,
        "V52_to_He4_Sc48",  // 231,
        "Cr48_to_p_V47",  // 232,
        "Cr48_to_He4_Ti44",  // 233,
        "Cr49_to_n_Cr48",  // 234,
        "Cr49_to_p_V48",  // 235,
        "Cr49_to_He4_Ti45",  // 236,
        "Cr50_to_n_Cr49",  // 237,
        "Cr50_to_p_V49",  // 238,
        "Cr50_to_He4_Ti46",  // 239,
        "Cr51_to_n_Cr50",  // 240,
        "Cr51_to_p_V50",  // 241,
        "Cr51_to_He4_Ti47",  // 242,
        "Cr52_to_n_Cr51",  // 243,
        "Cr52_to_p_V51",  // 244,
        "Cr52_to_He4_Ti48",  // 245,
        "Cr53_to_n_Cr52",  // 246,
        "Cr53_to_p_V52",  // 247,
        "Cr53_to_He4_Ti49",  // 248,
        "Cr54_to_n_Cr53",  // 249,
        "Cr54_to_He4_Ti50",  // 250,
        "Mn50_to_p_Cr49",  // 251,
        "Mn50_to_He4_V46",  // 252,
        "Mn51_to_n_Mn50",  // 253,
        "Mn51_to_p_Cr50",  // 254,
        "Mn51_to_He4_V47",  // 255,
        "Mn52_to_n_Mn51",  // 256,
        "Mn52_to_p_Cr51",  // 257,
        "Mn52_to_He4_V48",  // 258,
        "Mn53_to_n_Mn52",  // 259,
        "Mn53_to_p_Cr52",  // 260,
        "Mn53_to_He4_V49",  // 261,
        "Mn54_to_n_Mn53",  // 262,
        "Mn54_to_p_Cr53",  // 263,
        "Mn54_to_He4_V50",  // 264,
        "Mn55_to_n_Mn54",  // 265,
        "Mn55_to_p_Cr54",  // 266,
        "Mn55_to_He4_V51",  // 267,
        "Fe52_to_p_Mn51",  // 268,
        "Fe52_to_He4_Cr48",  // 269,
        "Fe53_to_n_Fe52",  // 270,
        "Fe53_to_p_Mn52",  // 271,
        "Fe53_to_He4_Cr49",  // 272,
        "Fe54_to_n_Fe53",  // 273,
        "Fe54_to_p_Mn53",  // 274,
        "Fe54_to_He4_Cr50",  // 275,
        "Fe55_to_n_Fe54",  // 276,
        "Fe55_to_p_Mn54",  // 277,
        "Fe55_to_He4_Cr51",  // 278,
        "Fe56_to_n_Fe55",  // 279,
        "Fe56_to_p_Mn55",  // 280,
        "Fe56_to_He4_Cr52",  // 281,
        "Fe57_to_n_Fe56",  // 282,
        "Fe57_to_He4_Cr53",  // 283,
        "Fe58_to_n_Fe57",  // 284,
        "Fe58_to_He4_Cr54",  // 285,
        "Co53_to_p_Fe52",  // 286,
        "Co54_to_n_Co53",  // 287,
        "Co54_to_p_Fe53",  // 288,
        "Co54_to_He4_Mn50",  // 289,
        "Co55_to_n_Co54",  // 290,
        "Co55_to_p_Fe54",  // 291,
        "Co55_to_He4_Mn51",  // 292,
        "Co56_to_n_Co55",  // 293,
        "Co56_to_p_Fe55",  // 294,
        "Co56_to_He4_Mn52",  // 295,
        "Co57_to_n_Co56",  // 296,
        "Co57_to_p_Fe56",  // 297,
        "Co57_to_He4_Mn53",  // 298,
        "Co58_to_n_Co57",  // 299,
        "Co58_to_p_Fe57",  // 300,
        "Co58_to_He4_Mn54",  // 301,
        "Co59_to_n_Co58",  // 302,
        "Co59_to_p_Fe58",  // 303,
        "Co59_to_He4_Mn55",  // 304,
        "Ni56_to_p_Co55",  // 305,
        "Ni56_to_He4_Fe52",  // 306,
        "Ni57_to_n_Ni56",  // 307,
        "Ni57_to_p_Co56",  // 308,
        "Ni57_to_He4_Fe53",  // 309,
        "Ni58_to_n_Ni57",  // 310,
        "Ni58_to_p_Co57",  // 311,
        "Ni58_to_He4_Fe54",  // 312,
        "Ni59_to_n_Ni58",  // 313,
        "Ni59_to_p_Co58",  // 314,
        "Ni59_to_He4_Fe55",  // 315,
        "Ni60_to_n_Ni59",  // 316,
        "Ni60_to_p_Co59",  // 317,
        "Ni60_to_He4_Fe56",  // 318,
        "Ni61_to_n_Ni60",  // 319,
        "Ni61_to_He4_Fe57",  // 320,
        "Ni62_to_n_Ni61",  // 321,
        "Ni62_to_He4_Fe58",  // 322,
        "Ni63_to_n_Ni62",  // 323,
        "Ni64_to_n_Ni63",  // 324,
        "Cu57_to_p_Ni56",  // 325,
        "Cu57_to_He4_Co53",  // 326,
        "Cu58_to_n_Cu57",  // 327,
        "Cu58_to_p_Ni57",  // 328,
        "Cu58_to_He4_Co54",  // 329,
        "Cu59_to_n_Cu58",  // 330,
        "Cu59_to_p_Ni58",  // 331,
        "Cu59_to_He4_Co55",  // 332,
        "Cu60_to_n_Cu59",  // 333,
        "Cu60_to_p_Ni59",  // 334,
        "Cu60_to_He4_Co56",  // 335,
        "Cu61_to_n_Cu60",  // 336,
        "Cu61_to_p_Ni60",  // 337,
        "Cu61_to_He4_Co57",  // 338,
        "Cu62_to_n_Cu61",  // 339,
        "Cu62_to_p_Ni61",  // 340,
        "Cu62_to_He4_Co58",  // 341,
        "Cu63_to_n_Cu62",  // 342,
        "Cu63_to_p_Ni62",  // 343,
        "Cu63_to_He4_Co59",  // 344,
        "Cu64_to_n_Cu63",  // 345,
        "Cu64_to_p_Ni63",  // 346,
        "Cu65_to_n_Cu64",  // 347,
        "Cu65_to_p_Ni64",  // 348,
        "Zn59_to_p_Cu58",  // 349,
        "Zn59_to_p_Ni58_weak_wc12",  // 350,
        "Zn60_to_n_Zn59",  // 351,
        "Zn60_to_p_Cu59",  // 352,
        "Zn60_to_He4_Ni56",  // 353,
        "Zn61_to_n_Zn60",  // 354,
        "Zn61_to_p_Cu60",  // 355,
        "Zn61_to_He4_Ni57",  // 356,
        "Zn62_to_n_Zn61",  // 357,
        "Zn62_to_p_Cu61",  // 358,
        "Zn62_to_He4_Ni58",  // 359,
        "Zn63_to_n_Zn62",  // 360,
        "Zn63_to_p_Cu62",  // 361,
        "Zn63_to_He4_Ni59",  // 362,
        "Zn64_to_n_Zn63",  // 363,
        "Zn64_to_p_Cu63",  // 364,
        "Zn64_to_He4_Ni60",  // 365,
        "Zn65_to_n_Zn64",  // 366,
        "Zn65_to_p_Cu64",  // 367,
        "Zn65_to_He4_Ni61",  // 368,
        "Zn66_to_n_Zn65",  // 369,
        "Zn66_to_p_Cu65",  // 370,
        "Zn66_to_He4_Ni62",  // 371,
        "Ga62_to_p_Zn61",  // 372,
        "Ga62_to_He4_Cu58",  // 373,
        "Ga63_to_n_Ga62",  // 374,
        "Ga63_to_p_Zn62",  // 375,
        "Ga63_to_He4_Cu59",  // 376,
        "Ga64_to_n_Ga63",  // 377,
        "Ga64_to_p_Zn63",  // 378,
        "Ga64_to_He4_Cu60",  // 379,
        "Ge63_to_p_Ga62",  // 380,
        "Ge63_to_He4_Zn59",  // 381,
        "Ge64_to_n_Ge63",  // 382,
        "Ge64_to_p_Ga63",  // 383,
        "Ge64_to_He4_Zn60",  // 384,
        "Li6_to_n_p_He4",  // 385,
        "Be9_to_n_He4_He4",  // 386,
        "C12_to_He4_He4_He4",  // 387,
        "n_p_to_d",  // 388,
        "p_p_to_d_weak_bet_pos_",  // 389,
        "p_p_to_d_weak_electron_capture",  // 390,
        "p_d_to_He3",  // 391,
        "d_d_to_He4",  // 392,
        "He4_d_to_Li6",  // 393,
        "n_He3_to_He4",  // 394,
        "p_He3_to_He4_weak_bet_pos_",  // 395,
        "He4_He3_to_Be7",  // 396,
        "n_Li6_to_Li7",  // 397,
        "p_Li6_to_Be7",  // 398,
        "He4_Li6_to_B10",  // 399,
        "He4_Li7_to_B11",  // 400,
        "p_Be7_to_B8",  // 401,
        "p_Be9_to_B10",  // 402,
        "n_B10_to_B11",  // 403,
        "p_B11_to_C12",  // 404,
        "n_C12_to_C13",  // 405,
        "p_C12_to_N13",  // 406,
        "He4_C12_to_O16",  // 407,
        "n_C13_to_C14",  // 408,
        "p_C13_to_N14",  // 409,
        "p_C14_to_N15",  // 410,
        "He4_C14_to_O18",  // 411,
        "n_N13_to_N14",  // 412,
        "p_N13_to_O14",  // 413,
        "n_N14_to_N15",  // 414,
        "p_N14_to_O15",  // 415,
        "He4_N14_to_F18",  // 416,
        "p_N15_to_O16",  // 417,
        "He4_N15_to_F19",  // 418,
        "n_O14_to_O15",  // 419,
        "He4_O14_to_Ne18",  // 420,
        "n_O15_to_O16",  // 421,
        "He4_O15_to_Ne19",  // 422,
        "n_O16_to_O17",  // 423,
        "p_O16_to_F17",  // 424,
        "He4_O16_to_Ne20",  // 425,
        "n_O17_to_O18",  // 426,
        "p_O17_to_F18",  // 427,
        "He4_O17_to_Ne21",  // 428,
        "p_O18_to_F19",  // 429,
        "He4_O18_to_Ne22",  // 430,
        "n_F17_to_F18",  // 431,
        "p_F17_to_Ne18",  // 432,
        "He4_F17_to_Na21",  // 433,
        "n_F18_to_F19",  // 434,
        "p_F18_to_Ne19",  // 435,
        "He4_F18_to_Na22",  // 436,
        "p_F19_to_Ne20",  // 437,
        "He4_F19_to_Na23",  // 438,
        "n_Ne18_to_Ne19",  // 439,
        "n_Ne19_to_Ne20",  // 440,
        "He4_Ne19_to_Mg23",  // 441,
        "n_Ne20_to_Ne21",  // 442,
        "p_Ne20_to_Na21",  // 443,
        "He4_Ne20_to_Mg24",  // 444,
        "n_Ne21_to_Ne22",  // 445,
        "p_Ne21_to_Na22",  // 446,
        "He4_Ne21_to_Mg25",  // 447,
        "p_Ne22_to_Na23",  // 448,
        "He4_Ne22_to_Mg26",  // 449,
        "n_Na21_to_Na22",  // 450,
        "He4_Na21_to_Al25",  // 451,
        "n_Na22_to_Na23",  // 452,
        "p_Na22_to_Mg23",  // 453,
        "He4_Na22_to_Al26",  // 454,
        "p_Na23_to_Mg24",  // 455,
        "He4_Na23_to_Al27",  // 456,
        "n_Mg23_to_Mg24",  // 457,
        "n_Mg24_to_Mg25",  // 458,
        "p_Mg24_to_Al25",  // 459,
        "He4_Mg24_to_Si28",  // 460,
        "n_Mg25_to_Mg26",  // 461,
        "p_Mg25_to_Al26",  // 462,
        "He4_Mg25_to_Si29",  // 463,
        "p_Mg26_to_Al27",  // 464,
        "He4_Mg26_to_Si30",  // 465,
        "n_Al25_to_Al26",  // 466,
        "He4_Al25_to_P29",  // 467,
        "n_Al26_to_Al27",  // 468,
        "He4_Al26_to_P30",  // 469,
        "p_Al27_to_Si28",  // 470,
        "He4_Al27_to_P31",  // 471,
        "n_Si28_to_Si29",  // 472,
        "p_Si28_to_P29",  // 473,
        "He4_Si28_to_S32",  // 474,
        "n_Si29_to_Si30",  // 475,
        "p_Si29_to_P30",  // 476,
        "He4_Si29_to_S33",  // 477,
        "n_Si30_to_Si31",  // 478,
        "p_Si30_to_P31",  // 479,
        "He4_Si30_to_S34",  // 480,
        "n_Si31_to_Si32",  // 481,
        "p_Si31_to_P32",  // 482,
        "He4_Si31_to_S35",  // 483,
        "p_Si32_to_P33",  // 484,
        "He4_Si32_to_S36",  // 485,
        "n_P29_to_P30",  // 486,
        "He4_P29_to_Cl33",  // 487,
        "n_P30_to_P31",  // 488,
        "He4_P30_to_Cl34",  // 489,
        "n_P31_to_P32",  // 490,
        "p_P31_to_S32",  // 491,
        "He4_P31_to_Cl35",  // 492,
        "n_P32_to_P33",  // 493,
        "p_P32_to_S33",  // 494,
        "He4_P32_to_Cl36",  // 495,
        "p_P33_to_S34",  // 496,
        "He4_P33_to_Cl37",  // 497,
        "n_S32_to_S33",  // 498,
        "p_S32_to_Cl33",  // 499,
        "He4_S32_to_Ar36",  // 500,
        "n_S33_to_S34",  // 501,
        "p_S33_to_Cl34",  // 502,
        "He4_S33_to_Ar37",  // 503,
        "n_S34_to_S35",  // 504,
        "p_S34_to_Cl35",  // 505,
        "He4_S34_to_Ar38",  // 506,
        "n_S35_to_S36",  // 507,
        "p_S35_to_Cl36",  // 508,
        "He4_S35_to_Ar39",  // 509,
        "p_S36_to_Cl37",  // 510,
        "He4_S36_to_Ar40",  // 511,
        "n_Cl33_to_Cl34",  // 512,
        "He4_Cl33_to_K37",  // 513,
        "n_Cl34_to_Cl35",  // 514,
        "He4_Cl34_to_K38",  // 515,
        "n_Cl35_to_Cl36",  // 516,
        "p_Cl35_to_Ar36",  // 517,
        "He4_Cl35_to_K39",  // 518,
        "n_Cl36_to_Cl37",  // 519,
        "p_Cl36_to_Ar37",  // 520,
        "He4_Cl36_to_K40",  // 521,
        "p_Cl37_to_Ar38",  // 522,
        "He4_Cl37_to_K41",  // 523,
        "n_Ar36_to_Ar37",  // 524,
        "p_Ar36_to_K37",  // 525,
        "He4_Ar36_to_Ca40",  // 526,
        "n_Ar37_to_Ar38",  // 527,
        "p_Ar37_to_K38",  // 528,
        "He4_Ar37_to_Ca41",  // 529,
        "n_Ar38_to_Ar39",  // 530,
        "p_Ar38_to_K39",  // 531,
        "He4_Ar38_to_Ca42",  // 532,
        "n_Ar39_to_Ar40",  // 533,
        "p_Ar39_to_K40",  // 534,
        "He4_Ar39_to_Ca43",  // 535,
        "p_Ar40_to_K41",  // 536,
        "He4_Ar40_to_Ca44",  // 537,
        "n_K37_to_K38",  // 538,
        "n_K38_to_K39",  // 539,
        "n_K39_to_K40",  // 540,
        "p_K39_to_Ca40",  // 541,
        "He4_K39_to_Sc43",  // 542,
        "n_K40_to_K41",  // 543,
        "p_K40_to_Ca41",  // 544,
        "He4_K40_to_Sc44",  // 545,
        "p_K41_to_Ca42",  // 546,
        "He4_K41_to_Sc45",  // 547,
        "n_Ca40_to_Ca41",  // 548,
        "He4_Ca40_to_Ti44",  // 549,
        "n_Ca41_to_Ca42",  // 550,
        "He4_Ca41_to_Ti45",  // 551,
        "n_Ca42_to_Ca43",  // 552,
        "p_Ca42_to_Sc43",  // 553,
        "He4_Ca42_to_Ti46",  // 554,
        "n_Ca43_to_Ca44",  // 555,
        "p_Ca43_to_Sc44",  // 556,
        "He4_Ca43_to_Ti47",  // 557,
        "n_Ca44_to_Ca45",  // 558,
        "p_Ca44_to_Sc45",  // 559,
        "He4_Ca44_to_Ti48",  // 560,
        "n_Ca45_to_Ca46",  // 561,
        "p_Ca45_to_Sc46",  // 562,
        "He4_Ca45_to_Ti49",  // 563,
        "n_Ca46_to_Ca47",  // 564,
        "p_Ca46_to_Sc47",  // 565,
        "He4_Ca46_to_Ti50",  // 566,
        "n_Ca47_to_Ca48",  // 567,
        "p_Ca47_to_Sc48",  // 568,
        "He4_Ca47_to_Ti51",  // 569,
        "p_Ca48_to_Sc49",  // 570,
        "n_Sc43_to_Sc44",  // 571,
        "p_Sc43_to_Ti44",  // 572,
        "He4_Sc43_to_V47",  // 573,
        "n_Sc44_to_Sc45",  // 574,
        "p_Sc44_to_Ti45",  // 575,
        "He4_Sc44_to_V48",  // 576,
        "n_Sc45_to_Sc46",  // 577,
        "p_Sc45_to_Ti46",  // 578,
        "He4_Sc45_to_V49",  // 579,
        "n_Sc46_to_Sc47",  // 580,
        "p_Sc46_to_Ti47",  // 581,
        "He4_Sc46_to_V50",  // 582,
        "n_Sc47_to_Sc48",  // 583,
        "p_Sc47_to_Ti48",  // 584,
        "He4_Sc47_to_V51",  // 585,
        "n_Sc48_to_Sc49",  // 586,
        "p_Sc48_to_Ti49",  // 587,
        "He4_Sc48_to_V52",  // 588,
        "p_Sc49_to_Ti50",  // 589,
        "n_Ti44_to_Ti45",  // 590,
        "He4_Ti44_to_Cr48",  // 591,
        "n_Ti45_to_Ti46",  // 592,
        "p_Ti45_to_V46",  // 593,
        "He4_Ti45_to_Cr49",  // 594,
        "n_Ti46_to_Ti47",  // 595,
        "p_Ti46_to_V47",  // 596,
        "He4_Ti46_to_Cr50",  // 597,
        "n_Ti47_to_Ti48",  // 598,
        "p_Ti47_to_V48",  // 599,
        "He4_Ti47_to_Cr51",  // 600,
        "n_Ti48_to_Ti49",  // 601,
        "p_Ti48_to_V49",  // 602,
        "He4_Ti48_to_Cr52",  // 603,
        "n_Ti49_to_Ti50",  // 604,
        "p_Ti49_to_V50",  // 605,
        "He4_Ti49_to_Cr53",  // 606,
        "n_Ti50_to_Ti51",  // 607,
        "p_Ti50_to_V51",  // 608,
        "He4_Ti50_to_Cr54",  // 609,
        "p_Ti51_to_V52",  // 610,
        "n_V46_to_V47",  // 611,
        "He4_V46_to_Mn50",  // 612,
        "n_V47_to_V48",  // 613,
        "p_V47_to_Cr48",  // 614,
        "He4_V47_to_Mn51",  // 615,
        "n_V48_to_V49",  // 616,
        "p_V48_to_Cr49",  // 617,
        "He4_V48_to_Mn52",  // 618,
        "n_V49_to_V50",  // 619,
        "p_V49_to_Cr50",  // 620,
        "He4_V49_to_Mn53",  // 621,
        "n_V50_to_V51",  // 622,
        "p_V50_to_Cr51",  // 623,
        "He4_V50_to_Mn54",  // 624,
        "n_V51_to_V52",  // 625,
        "p_V51_to_Cr52",  // 626,
        "He4_V51_to_Mn55",  // 627,
        "p_V52_to_Cr53",  // 628,
        "n_Cr48_to_Cr49",  // 629,
        "He4_Cr48_to_Fe52",  // 630,
        "n_Cr49_to_Cr50",  // 631,
        "p_Cr49_to_Mn50",  // 632,
        "He4_Cr49_to_Fe53",  // 633,
        "n_Cr50_to_Cr51",  // 634,
        "p_Cr50_to_Mn51",  // 635,
        "He4_Cr50_to_Fe54",  // 636,
        "n_Cr51_to_Cr52",  // 637,
        "p_Cr51_to_Mn52",  // 638,
        "He4_Cr51_to_Fe55",  // 639,
        "n_Cr52_to_Cr53",  // 640,
        "p_Cr52_to_Mn53",  // 641,
        "He4_Cr52_to_Fe56",  // 642,
        "n_Cr53_to_Cr54",  // 643,
        "p_Cr53_to_Mn54",  // 644,
        "He4_Cr53_to_Fe57",  // 645,
        "p_Cr54_to_Mn55",  // 646,
        "He4_Cr54_to_Fe58",  // 647,
        "n_Mn50_to_Mn51",  // 648,
        "He4_Mn50_to_Co54",  // 649,
        "n_Mn51_to_Mn52",  // 650,
        "p_Mn51_to_Fe52",  // 651,
        "He4_Mn51_to_Co55",  // 652,
        "n_Mn52_to_Mn53",  // 653,
        "p_Mn52_to_Fe53",  // 654,
        "He4_Mn52_to_Co56",  // 655,
        "n_Mn53_to_Mn54",  // 656,
        "p_Mn53_to_Fe54",  // 657,
        "He4_Mn53_to_Co57",  // 658,
        "n_Mn54_to_Mn55",  // 659,
        "p_Mn54_to_Fe55",  // 660,
        "He4_Mn54_to_Co58",  // 661,
        "p_Mn55_to_Fe56",  // 662,
        "He4_Mn55_to_Co59",  // 663,
        "n_Fe52_to_Fe53",  // 664,
        "p_Fe52_to_Co53",  // 665,
        "He4_Fe52_to_Ni56",  // 666,
        "n_Fe53_to_Fe54",  // 667,
        "p_Fe53_to_Co54",  // 668,
        "He4_Fe53_to_Ni57",  // 669,
        "n_Fe54_to_Fe55",  // 670,
        "p_Fe54_to_Co55",  // 671,
        "He4_Fe54_to_Ni58",  // 672,
        "n_Fe55_to_Fe56",  // 673,
        "p_Fe55_to_Co56",  // 674,
        "He4_Fe55_to_Ni59",  // 675,
        "n_Fe56_to_Fe57",  // 676,
        "p_Fe56_to_Co57",  // 677,
        "He4_Fe56_to_Ni60",  // 678,
        "n_Fe57_to_Fe58",  // 679,
        "p_Fe57_to_Co58",  // 680,
        "He4_Fe57_to_Ni61",  // 681,
        "p_Fe58_to_Co59",  // 682,
        "He4_Fe58_to_Ni62",  // 683,
        "n_Co53_to_Co54",  // 684,
        "He4_Co53_to_Cu57",  // 685,
        "n_Co54_to_Co55",  // 686,
        "He4_Co54_to_Cu58",  // 687,
        "n_Co55_to_Co56",  // 688,
        "p_Co55_to_Ni56",  // 689,
        "He4_Co55_to_Cu59",  // 690,
        "n_Co56_to_Co57",  // 691,
        "p_Co56_to_Ni57",  // 692,
        "He4_Co56_to_Cu60",  // 693,
        "n_Co57_to_Co58",  // 694,
        "p_Co57_to_Ni58",  // 695,
        "He4_Co57_to_Cu61",  // 696,
        "n_Co58_to_Co59",  // 697,
        "p_Co58_to_Ni59",  // 698,
        "He4_Co58_to_Cu62",  // 699,
        "p_Co59_to_Ni60",  // 700,
        "He4_Co59_to_Cu63",  // 701,
        "n_Ni56_to_Ni57",  // 702,
        "p_Ni56_to_Cu57",  // 703,
        "He4_Ni56_to_Zn60",  // 704,
        "n_Ni57_to_Ni58",  // 705,
        "p_Ni57_to_Cu58",  // 706,
        "He4_Ni57_to_Zn61",  // 707,
        "n_Ni58_to_Ni59",  // 708,
        "p_Ni58_to_Cu59",  // 709,
        "He4_Ni58_to_Zn62",  // 710,
        "n_Ni59_to_Ni60",  // 711,
        "p_Ni59_to_Cu60",  // 712,
        "He4_Ni59_to_Zn63",  // 713,
        "n_Ni60_to_Ni61",  // 714,
        "p_Ni60_to_Cu61",  // 715,
        "He4_Ni60_to_Zn64",  // 716,
        "n_Ni61_to_Ni62",  // 717,
        "p_Ni61_to_Cu62",  // 718,
        "He4_Ni61_to_Zn65",  // 719,
        "n_Ni62_to_Ni63",  // 720,
        "p_Ni62_to_Cu63",  // 721,
        "He4_Ni62_to_Zn66",  // 722,
        "n_Ni63_to_Ni64",  // 723,
        "p_Ni63_to_Cu64",  // 724,
        "p_Ni64_to_Cu65",  // 725,
        "n_Cu57_to_Cu58",  // 726,
        "n_Cu58_to_Cu59",  // 727,
        "p_Cu58_to_Zn59",  // 728,
        "He4_Cu58_to_Ga62",  // 729,
        "n_Cu59_to_Cu60",  // 730,
        "p_Cu59_to_Zn60",  // 731,
        "He4_Cu59_to_Ga63",  // 732,
        "n_Cu60_to_Cu61",  // 733,
        "p_Cu60_to_Zn61",  // 734,
        "He4_Cu60_to_Ga64",  // 735,
        "n_Cu61_to_Cu62",  // 736,
        "p_Cu61_to_Zn62",  // 737,
        "n_Cu62_to_Cu63",  // 738,
        "p_Cu62_to_Zn63",  // 739,
        "n_Cu63_to_Cu64",  // 740,
        "p_Cu63_to_Zn64",  // 741,
        "n_Cu64_to_Cu65",  // 742,
        "p_Cu64_to_Zn65",  // 743,
        "p_Cu65_to_Zn66",  // 744,
        "n_Zn59_to_Zn60",  // 745,
        "He4_Zn59_to_Ge63",  // 746,
        "n_Zn60_to_Zn61",  // 747,
        "He4_Zn60_to_Ge64",  // 748,
        "n_Zn61_to_Zn62",  // 749,
        "p_Zn61_to_Ga62",  // 750,
        "n_Zn62_to_Zn63",  // 751,
        "p_Zn62_to_Ga63",  // 752,
        "n_Zn63_to_Zn64",  // 753,
        "p_Zn63_to_Ga64",  // 754,
        "n_Zn64_to_Zn65",  // 755,
        "n_Zn65_to_Zn66",  // 756,
        "n_Ga62_to_Ga63",  // 757,
        "p_Ga62_to_Ge63",  // 758,
        "n_Ga63_to_Ga64",  // 759,
        "p_Ga63_to_Ge64",  // 760,
        "n_Ge63_to_Ge64",  // 761,
        "d_d_to_n_He3",  // 762,
        "n_He3_to_d_d",  // 763,
        "d_He3_to_p_He4",  // 764,
        "He4_He3_to_p_Li6",  // 765,
        "p_He4_to_d_He3",  // 766,
        "He4_He4_to_n_Be7",  // 767,
        "He4_He4_to_p_Li7",  // 768,
        "p_Li6_to_He4_He3",  // 769,
        "d_Li6_to_n_Be7",  // 770,
        "d_Li6_to_p_Li7",  // 771,
        "He4_Li6_to_p_Be9",  // 772,
        "p_Li7_to_n_Be7",  // 773,
        "p_Li7_to_d_Li6",  // 774,
        "p_Li7_to_He4_He4",  // 775,
        "He4_Li7_to_n_B10",  // 776,
        "n_Be7_to_p_Li7",  // 777,
        "n_Be7_to_d_Li6",  // 778,
        "n_Be7_to_He4_He4",  // 779,
        "He4_Be7_to_p_B10",  // 780,
        "p_Be9_to_He4_Li6",  // 781,
        "He4_Be9_to_n_C12",  // 782,
        "n_B10_to_He4_Li7",  // 783,
        "p_B10_to_He4_Be7",  // 784,
        "He4_B10_to_n_N13",  // 785,
        "He4_B10_to_p_C13",  // 786,
        "He4_B11_to_n_N14",  // 787,
        "He4_B11_to_p_C14",  // 788,
        "n_C12_to_He4_Be9",  // 789,
        "He4_C12_to_n_O15",  // 790,
        "He4_C12_to_p_N15",  // 791,
        "C12_C12_to_n_Mg23",  // 792,
        "C12_C12_to_p_Na23",  // 793,
        "C12_C12_to_He4_Ne20",  // 794,
        "p_C13_to_n_N13",  // 795,
        "p_C13_to_He4_B10",  // 796,
        "d_C13_to_n_N14",  // 797,
        "He4_C13_to_n_O16",  // 798,
        "p_C14_to_n_N14",  // 799,
        "p_C14_to_He4_B11",  // 800,
        "d_C14_to_n_N15",  // 801,
        "He4_C14_to_n_O17",  // 802,
        "n_N13_to_p_C13",  // 803,
        "n_N13_to_He4_B10",  // 804,
        "He4_N13_to_p_O16",  // 805,
        "n_N14_to_p_C14",  // 806,
        "n_N14_to_d_C13",  // 807,
        "n_N14_to_He4_B11",  // 808,
        "p_N14_to_n_O14",  // 809,
        "He4_N14_to_n_F17",  // 810,
        "He4_N14_to_p_O17",  // 811,
        "n_N15_to_d_C14",  // 812,
        "p_N15_to_n_O15",  // 813,
        "p_N15_to_He4_C12",  // 814,
        "He4_N15_to_n_F18",  // 815,
        "He4_N15_to_p_O18",  // 816,
        "n_O14_to_p_N14",  // 817,
        "He4_O14_to_p_F17",  // 818,
        "n_O15_to_p_N15",  // 819,
        "n_O15_to_He4_C12",  // 820,
        "He4_O15_to_n_Ne18",  // 821,
        "He4_O15_to_p_F18",  // 822,
        "n_O16_to_He4_C13",  // 823,
        "p_O16_to_He4_N13",  // 824,
        "He4_O16_to_n_Ne19",  // 825,
        "He4_O16_to_p_F19",  // 826,
        "C12_O16_to_p_Al27",  // 827,
        "C12_O16_to_He4_Mg24",  // 828,
        "O16_O16_to_p_P31",  // 829,
        "O16_O16_to_He4_Si28",  // 830,
        "n_O17_to_He4_C14",  // 831,
        "p_O17_to_n_F17",  // 832,
        "p_O17_to_He4_N14",  // 833,
        "He4_O17_to_n_Ne20",  // 834,
        "p_O18_to_n_F18",  // 835,
        "p_O18_to_He4_N15",  // 836,
        "He4_O18_to_n_Ne21",  // 837,
        "n_F17_to_p_O17",  // 838,
        "n_F17_to_He4_N14",  // 839,
        "p_F17_to_He4_O14",  // 840,
        "He4_F17_to_p_Ne20",  // 841,
        "n_F18_to_p_O18",  // 842,
        "n_F18_to_He4_N15",  // 843,
        "p_F18_to_n_Ne18",  // 844,
        "p_F18_to_He4_O15",  // 845,
        "He4_F18_to_n_Na21",  // 846,
        "He4_F18_to_p_Ne21",  // 847,
        "p_F19_to_n_Ne19",  // 848,
        "p_F19_to_He4_O16",  // 849,
        "He4_F19_to_n_Na22",  // 850,
        "He4_F19_to_p_Ne22",  // 851,
        "n_Ne18_to_p_F18",  // 852,
        "n_Ne18_to_He4_O15",  // 853,
        "He4_Ne18_to_p_Na21",  // 854,
        "n_Ne19_to_p_F19",  // 855,
        "n_Ne19_to_He4_O16",  // 856,
        "He4_Ne19_to_p_Na22",  // 857,
        "n_Ne20_to_He4_O17",  // 858,
        "p_Ne20_to_He4_F17",  // 859,
        "He4_Ne20_to_n_Mg23",  // 860,
        "He4_Ne20_to_p_Na23",  // 861,
        "He4_Ne20_to_C12_C12",  // 862,
        "C12_Ne20_to_p_P31",  // 863,
        "C12_Ne20_to_He4_Si28",  // 864,
        "n_Ne21_to_He4_O18",  // 865,
        "p_Ne21_to_n_Na21",  // 866,
        "p_Ne21_to_He4_F18",  // 867,
        "He4_Ne21_to_n_Mg24",  // 868,
        "p_Ne22_to_n_Na22",  // 869,
        "p_Ne22_to_He4_F19",  // 870,
        "He4_Ne22_to_n_Mg25",  // 871,
        "n_Na21_to_p_Ne21",  // 872,
        "n_Na21_to_He4_F18",  // 873,
        "p_Na21_to_He4_Ne18",  // 874,
        "He4_Na21_to_p_Mg24",  // 875,
        "n_Na22_to_p_Ne22",  // 876,
        "n_Na22_to_He4_F19",  // 877,
        "p_Na22_to_He4_Ne19",  // 878,
        "He4_Na22_to_n_Al25",  // 879,
        "He4_Na22_to_p_Mg25",  // 880,
        "p_Na23_to_n_Mg23",  // 881,
        "p_Na23_to_He4_Ne20",  // 882,
        "p_Na23_to_C12_C12",  // 883,
        "He4_Na23_to_n_Al26",  // 884,
        "He4_Na23_to_p_Mg26",  // 885,
        "n_Mg23_to_p_Na23",  // 886,
        "n_Mg23_to_He4_Ne20",  // 887,
        "n_Mg23_to_C12_C12",  // 888,
        "He4_Mg23_to_p_Al26",  // 889,
        "n_Mg24_to_He4_Ne21",  // 890,
        "p_Mg24_to_He4_Na21",  // 891,
        "He4_Mg24_to_p_Al27",  // 892,
        "He4_Mg24_to_C12_O16",  // 893,
        "n_Mg25_to_He4_Ne22",  // 894,
        "p_Mg25_to_n_Al25",  // 895,
        "p_Mg25_to_He4_Na22",  // 896,
        "He4_Mg25_to_n_Si28",  // 897,
        "p_Mg26_to_n_Al26",  // 898,
        "p_Mg26_to_He4_Na23",  // 899,
        "He4_Mg26_to_n_Si29",  // 900,
        "n_Al25_to_p_Mg25",  // 901,
        "n_Al25_to_He4_Na22",  // 902,
        "He4_Al25_to_p_Si28",  // 903,
        "n_Al26_to_p_Mg26",  // 904,
        "n_Al26_to_He4_Na23",  // 905,
        "p_Al26_to_He4_Mg23",  // 906,
        "He4_Al26_to_n_P29",  // 907,
        "He4_Al26_to_p_Si29",  // 908,
        "p_Al27_to_He4_Mg24",  // 909,
        "p_Al27_to_C12_O16",  // 910,
        "He4_Al27_to_n_P30",  // 911,
        "He4_Al27_to_p_Si30",  // 912,
        "n_Si28_to_He4_Mg25",  // 913,
        "p_Si28_to_He4_Al25",  // 914,
        "He4_Si28_to_p_P31",  // 915,
        "He4_Si28_to_C12_Ne20",  // 916,
        "He4_Si28_to_O16_O16",  // 917,
        "n_Si29_to_He4_Mg26",  // 918,
        "p_Si29_to_n_P29",  // 919,
        "p_Si29_to_He4_Al26",  // 920,
        "He4_Si29_to_n_S32",  // 921,
        "He4_Si29_to_p_P32",  // 922,
        "p_Si30_to_n_P30",  // 923,
        "p_Si30_to_He4_Al27",  // 924,
        "He4_Si30_to_n_S33",  // 925,
        "He4_Si30_to_p_P33",  // 926,
        "p_Si31_to_n_P31",  // 927,
        "He4_Si31_to_n_S34",  // 928,
        "p_Si32_to_n_P32",  // 929,
        "He4_Si32_to_n_S35",  // 930,
        "n_P29_to_p_Si29",  // 931,
        "n_P29_to_He4_Al26",  // 932,
        "He4_P29_to_p_S32",  // 933,
        "n_P30_to_p_Si30",  // 934,
        "n_P30_to_He4_Al27",  // 935,
        "He4_P30_to_n_Cl33",  // 936,
        "He4_P30_to_p_S33",  // 937,
        "n_P31_to_p_Si31",  // 938,
        "p_P31_to_He4_Si28",  // 939,
        "p_P31_to_C12_Ne20",  // 940,
        "p_P31_to_O16_O16",  // 941,
        "He4_P31_to_n_Cl34",  // 942,
        "He4_P31_to_p_S34",  // 943,
        "n_P32_to_p_Si32",  // 944,
        "p_P32_to_n_S32",  // 945,
        "p_P32_to_He4_Si29",  // 946,
        "He4_P32_to_n_Cl35",  // 947,
        "He4_P32_to_p_S35",  // 948,
        "p_P33_to_n_S33",  // 949,
        "p_P33_to_He4_Si30",  // 950,
        "He4_P33_to_n_Cl36",  // 951,
        "He4_P33_to_p_S36",  // 952,
        "n_S32_to_p_P32",  // 953,
        "n_S32_to_He4_Si29",  // 954,
        "p_S32_to_He4_P29",  // 955,
        "He4_S32_to_p_Cl35",  // 956,
        "n_S33_to_p_P33",  // 957,
        "n_S33_to_He4_Si30",  // 958,
        "p_S33_to_n_Cl33",  // 959,
        "p_S33_to_He4_P30",  // 960,
        "He4_S33_to_n_Ar36",  // 961,
        "He4_S33_to_p_Cl36",  // 962,
        "n_S34_to_He4_Si31",  // 963,
        "p_S34_to_n_Cl34",  // 964,
        "p_S34_to_He4_P31",  // 965,
        "He4_S34_to_n_Ar37",  // 966,
        "He4_S34_to_p_Cl37",  // 967,
        "n_S35_to_He4_Si32",  // 968,
        "p_S35_to_n_Cl35",  // 969,
        "p_S35_to_He4_P32",  // 970,
        "He4_S35_to_n_Ar38",  // 971,
        "p_S36_to_n_Cl36",  // 972,
        "p_S36_to_He4_P33",  // 973,
        "He4_S36_to_n_Ar39",  // 974,
        "n_Cl33_to_p_S33",  // 975,
        "n_Cl33_to_He4_P30",  // 976,
        "He4_Cl33_to_p_Ar36",  // 977,
        "n_Cl34_to_p_S34",  // 978,
        "n_Cl34_to_He4_P31",  // 979,
        "He4_Cl34_to_n_K37",  // 980,
        "He4_Cl34_to_p_Ar37",  // 981,
        "n_Cl35_to_p_S35",  // 982,
        "n_Cl35_to_He4_P32",  // 983,
        "p_Cl35_to_He4_S32",  // 984,
        "He4_Cl35_to_n_K38",  // 985,
        "He4_Cl35_to_p_Ar38",  // 986,
        "n_Cl36_to_p_S36",  // 987,
        "n_Cl36_to_He4_P33",  // 988,
        "p_Cl36_to_n_Ar36",  // 989,
        "p_Cl36_to_He4_S33",  // 990,
        "He4_Cl36_to_n_K39",  // 991,
        "He4_Cl36_to_p_Ar39",  // 992,
        "p_Cl37_to_n_Ar37",  // 993,
        "p_Cl37_to_He4_S34",  // 994,
        "He4_Cl37_to_n_K40",  // 995,
        "He4_Cl37_to_p_Ar40",  // 996,
        "n_Ar36_to_p_Cl36",  // 997,
        "n_Ar36_to_He4_S33",  // 998,
        "p_Ar36_to_He4_Cl33",  // 999,
        "He4_Ar36_to_p_K39",  // 1000,
        "n_Ar37_to_p_Cl37",  // 1001,
        "n_Ar37_to_He4_S34",  // 1002,
        "p_Ar37_to_n_K37",  // 1003,
        "p_Ar37_to_He4_Cl34",  // 1004,
        "He4_Ar37_to_n_Ca40",  // 1005,
        "He4_Ar37_to_p_K40",  // 1006,
        "n_Ar38_to_He4_S35",  // 1007,
        "p_Ar38_to_n_K38",  // 1008,
        "p_Ar38_to_He4_Cl35",  // 1009,
        "He4_Ar38_to_n_Ca41",  // 1010,
        "He4_Ar38_to_p_K41",  // 1011,
        "n_Ar39_to_He4_S36",  // 1012,
        "p_Ar39_to_n_K39",  // 1013,
        "p_Ar39_to_He4_Cl36",  // 1014,
        "He4_Ar39_to_n_Ca42",  // 1015,
        "p_Ar40_to_n_K40",  // 1016,
        "p_Ar40_to_He4_Cl37",  // 1017,
        "He4_Ar40_to_n_Ca43",  // 1018,
        "n_K37_to_p_Ar37",  // 1019,
        "n_K37_to_He4_Cl34",  // 1020,
        "He4_K37_to_p_Ca40",  // 1021,
        "n_K38_to_p_Ar38",  // 1022,
        "n_K38_to_He4_Cl35",  // 1023,
        "He4_K38_to_p_Ca41",  // 1024,
        "n_K39_to_p_Ar39",  // 1025,
        "n_K39_to_He4_Cl36",  // 1026,
        "p_K39_to_He4_Ar36",  // 1027,
        "He4_K39_to_p_Ca42",  // 1028,
        "n_K40_to_p_Ar40",  // 1029,
        "n_K40_to_He4_Cl37",  // 1030,
        "p_K40_to_n_Ca40",  // 1031,
        "p_K40_to_He4_Ar37",  // 1032,
        "He4_K40_to_n_Sc43",  // 1033,
        "He4_K40_to_p_Ca43",  // 1034,
        "p_K41_to_n_Ca41",  // 1035,
        "p_K41_to_He4_Ar38",  // 1036,
        "He4_K41_to_n_Sc44",  // 1037,
        "He4_K41_to_p_Ca44",  // 1038,
        "n_Ca40_to_p_K40",  // 1039,
        "n_Ca40_to_He4_Ar37",  // 1040,
        "p_Ca40_to_He4_K37",  // 1041,
        "He4_Ca40_to_p_Sc43",  // 1042,
        "n_Ca41_to_p_K41",  // 1043,
        "n_Ca41_to_He4_Ar38",  // 1044,
        "p_Ca41_to_He4_K38",  // 1045,
        "He4_Ca41_to_n_Ti44",  // 1046,
        "He4_Ca41_to_p_Sc44",  // 1047,
        "n_Ca42_to_He4_Ar39",  // 1048,
        "p_Ca42_to_He4_K39",  // 1049,
        "He4_Ca42_to_n_Ti45",  // 1050,
        "He4_Ca42_to_p_Sc45",  // 1051,
        "n_Ca43_to_He4_Ar40",  // 1052,
        "p_Ca43_to_n_Sc43",  // 1053,
        "p_Ca43_to_He4_K40",  // 1054,
        "He4_Ca43_to_n_Ti46",  // 1055,
        "He4_Ca43_to_p_Sc46",  // 1056,
        "p_Ca44_to_n_Sc44",  // 1057,
        "p_Ca44_to_He4_K41",  // 1058,
        "He4_Ca44_to_n_Ti47",  // 1059,
        "He4_Ca44_to_p_Sc47",  // 1060,
        "p_Ca45_to_n_Sc45",  // 1061,
        "He4_Ca45_to_n_Ti48",  // 1062,
        "He4_Ca45_to_p_Sc48",  // 1063,
        "p_Ca46_to_n_Sc46",  // 1064,
        "He4_Ca46_to_n_Ti49",  // 1065,
        "He4_Ca46_to_p_Sc49",  // 1066,
        "p_Ca47_to_n_Sc47",  // 1067,
        "He4_Ca47_to_n_Ti50",  // 1068,
        "p_Ca48_to_n_Sc48",  // 1069,
        "He4_Ca48_to_n_Ti51",  // 1070,
        "n_Sc43_to_p_Ca43",  // 1071,
        "n_Sc43_to_He4_K40",  // 1072,
        "p_Sc43_to_He4_Ca40",  // 1073,
        "He4_Sc43_to_n_V46",  // 1074,
        "He4_Sc43_to_p_Ti46",  // 1075,
        "n_Sc44_to_p_Ca44",  // 1076,
        "n_Sc44_to_He4_K41",  // 1077,
        "p_Sc44_to_n_Ti44",  // 1078,
        "p_Sc44_to_He4_Ca41",  // 1079,
        "He4_Sc44_to_n_V47",  // 1080,
        "He4_Sc44_to_p_Ti47",  // 1081,
        "n_Sc45_to_p_Ca45",  // 1082,
        "p_Sc45_to_n_Ti45",  // 1083,
        "p_Sc45_to_He4_Ca42",  // 1084,
        "He4_Sc45_to_n_V48",  // 1085,
        "He4_Sc45_to_p_Ti48",  // 1086,
        "n_Sc46_to_p_Ca46",  // 1087,
        "p_Sc46_to_n_Ti46",  // 1088,
        "p_Sc46_to_He4_Ca43",  // 1089,
        "He4_Sc46_to_n_V49",  // 1090,
        "He4_Sc46_to_p_Ti49",  // 1091,
        "n_Sc47_to_p_Ca47",  // 1092,
        "p_Sc47_to_n_Ti47",  // 1093,
        "p_Sc47_to_He4_Ca44",  // 1094,
        "He4_Sc47_to_n_V50",  // 1095,
        "He4_Sc47_to_p_Ti50",  // 1096,
        "n_Sc48_to_p_Ca48",  // 1097,
        "p_Sc48_to_n_Ti48",  // 1098,
        "p_Sc48_to_He4_Ca45",  // 1099,
        "He4_Sc48_to_n_V51",  // 1100,
        "He4_Sc48_to_p_Ti51",  // 1101,
        "p_Sc49_to_n_Ti49",  // 1102,
        "p_Sc49_to_He4_Ca46",  // 1103,
        "He4_Sc49_to_n_V52",  // 1104,
        "n_Ti44_to_p_Sc44",  // 1105,
        "n_Ti44_to_He4_Ca41",  // 1106,
        "He4_Ti44_to_p_V47",  // 1107,
        "n_Ti45_to_p_Sc45",  // 1108,
        "n_Ti45_to_He4_Ca42",  // 1109,
        "He4_Ti45_to_n_Cr48",  // 1110,
        "He4_Ti45_to_p_V48",  // 1111,
        "n_Ti46_to_p_Sc46",  // 1112,
        "n_Ti46_to_He4_Ca43",  // 1113,
        "p_Ti46_to_n_V46",  // 1114,
        "p_Ti46_to_He4_Sc43",  // 1115,
        "He4_Ti46_to_n_Cr49",  // 1116,
        "He4_Ti46_to_p_V49",  // 1117,
        "n_Ti47_to_p_Sc47",  // 1118,
        "n_Ti47_to_He4_Ca44",  // 1119,
        "p_Ti47_to_n_V47",  // 1120,
        "p_Ti47_to_He4_Sc44",  // 1121,
        "He4_Ti47_to_n_Cr50",  // 1122,
        "He4_Ti47_to_p_V50",  // 1123,
        "n_Ti48_to_p_Sc48",  // 1124,
        "n_Ti48_to_He4_Ca45",  // 1125,
        "p_Ti48_to_n_V48",  // 1126,
        "p_Ti48_to_He4_Sc45",  // 1127,
        "He4_Ti48_to_n_Cr51",  // 1128,
        "He4_Ti48_to_p_V51",  // 1129,
        "n_Ti49_to_p_Sc49",  // 1130,
        "n_Ti49_to_He4_Ca46",  // 1131,
        "p_Ti49_to_n_V49",  // 1132,
        "p_Ti49_to_He4_Sc46",  // 1133,
        "He4_Ti49_to_n_Cr52",  // 1134,
        "He4_Ti49_to_p_V52",  // 1135,
        "n_Ti50_to_He4_Ca47",  // 1136,
        "p_Ti50_to_n_V50",  // 1137,
        "p_Ti50_to_He4_Sc47",  // 1138,
        "He4_Ti50_to_n_Cr53",  // 1139,
        "n_Ti51_to_He4_Ca48",  // 1140,
        "p_Ti51_to_n_V51",  // 1141,
        "p_Ti51_to_He4_Sc48",  // 1142,
        "He4_Ti51_to_n_Cr54",  // 1143,
        "n_V46_to_p_Ti46",  // 1144,
        "n_V46_to_He4_Sc43",  // 1145,
        "He4_V46_to_p_Cr49",  // 1146,
        "n_V47_to_p_Ti47",  // 1147,
        "n_V47_to_He4_Sc44",  // 1148,
        "p_V47_to_He4_Ti44",  // 1149,
        "He4_V47_to_n_Mn50",  // 1150,
        "He4_V47_to_p_Cr50",  // 1151,
        "n_V48_to_p_Ti48",  // 1152,
        "n_V48_to_He4_Sc45",  // 1153,
        "p_V48_to_n_Cr48",  // 1154,
        "p_V48_to_He4_Ti45",  // 1155,
        "He4_V48_to_n_Mn51",  // 1156,
        "He4_V48_to_p_Cr51",  // 1157,
        "n_V49_to_p_Ti49",  // 1158,
        "n_V49_to_He4_Sc46",  // 1159,
        "p_V49_to_n_Cr49",  // 1160,
        "p_V49_to_He4_Ti46",  // 1161,
        "He4_V49_to_n_Mn52",  // 1162,
        "He4_V49_to_p_Cr52",  // 1163,
        "n_V50_to_p_Ti50",  // 1164,
        "n_V50_to_He4_Sc47",  // 1165,
        "p_V50_to_n_Cr50",  // 1166,
        "p_V50_to_He4_Ti47",  // 1167,
        "He4_V50_to_n_Mn53",  // 1168,
        "He4_V50_to_p_Cr53",  // 1169,
        "n_V51_to_p_Ti51",  // 1170,
        "n_V51_to_He4_Sc48",  // 1171,
        "p_V51_to_n_Cr51",  // 1172,
        "p_V51_to_He4_Ti48",  // 1173,
        "He4_V51_to_n_Mn54",  // 1174,
        "He4_V51_to_p_Cr54",  // 1175,
        "n_V52_to_He4_Sc49",  // 1176,
        "p_V52_to_n_Cr52",  // 1177,
        "p_V52_to_He4_Ti49",  // 1178,
        "He4_V52_to_n_Mn55",  // 1179,
        "n_Cr48_to_p_V48",  // 1180,
        "n_Cr48_to_He4_Ti45",  // 1181,
        "He4_Cr48_to_p_Mn51",  // 1182,
        "n_Cr49_to_p_V49",  // 1183,
        "n_Cr49_to_He4_Ti46",  // 1184,
        "p_Cr49_to_He4_V46",  // 1185,
        "He4_Cr49_to_n_Fe52",  // 1186,
        "He4_Cr49_to_p_Mn52",  // 1187,
        "n_Cr50_to_p_V50",  // 1188,
        "n_Cr50_to_He4_Ti47",  // 1189,
        "p_Cr50_to_n_Mn50",  // 1190,
        "p_Cr50_to_He4_V47",  // 1191,
        "He4_Cr50_to_n_Fe53",  // 1192,
        "He4_Cr50_to_p_Mn53",  // 1193,
        "n_Cr51_to_p_V51",  // 1194,
        "n_Cr51_to_He4_Ti48",  // 1195,
        "p_Cr51_to_n_Mn51",  // 1196,
        "p_Cr51_to_He4_V48",  // 1197,
        "He4_Cr51_to_n_Fe54",  // 1198,
        "He4_Cr51_to_p_Mn54",  // 1199,
        "n_Cr52_to_p_V52",  // 1200,
        "n_Cr52_to_He4_Ti49",  // 1201,
        "p_Cr52_to_n_Mn52",  // 1202,
        "p_Cr52_to_He4_V49",  // 1203,
        "He4_Cr52_to_n_Fe55",  // 1204,
        "He4_Cr52_to_p_Mn55",  // 1205,
        "n_Cr53_to_He4_Ti50",  // 1206,
        "p_Cr53_to_n_Mn53",  // 1207,
        "p_Cr53_to_He4_V50",  // 1208,
        "He4_Cr53_to_n_Fe56",  // 1209,
        "n_Cr54_to_He4_Ti51",  // 1210,
        "p_Cr54_to_n_Mn54",  // 1211,
        "p_Cr54_to_He4_V51",  // 1212,
        "He4_Cr54_to_n_Fe57",  // 1213,
        "n_Mn50_to_p_Cr50",  // 1214,
        "n_Mn50_to_He4_V47",  // 1215,
        "He4_Mn50_to_n_Co53",  // 1216,
        "He4_Mn50_to_p_Fe53",  // 1217,
        "n_Mn51_to_p_Cr51",  // 1218,
        "n_Mn51_to_He4_V48",  // 1219,
        "p_Mn51_to_He4_Cr48",  // 1220,
        "He4_Mn51_to_n_Co54",  // 1221,
        "He4_Mn51_to_p_Fe54",  // 1222,
        "n_Mn52_to_p_Cr52",  // 1223,
        "n_Mn52_to_He4_V49",  // 1224,
        "p_Mn52_to_n_Fe52",  // 1225,
        "p_Mn52_to_He4_Cr49",  // 1226,
        "He4_Mn52_to_n_Co55",  // 1227,
        "He4_Mn52_to_p_Fe55",  // 1228,
        "n_Mn53_to_p_Cr53",  // 1229,
        "n_Mn53_to_He4_V50",  // 1230,
        "p_Mn53_to_n_Fe53",  // 1231,
        "p_Mn53_to_He4_Cr50",  // 1232,
        "He4_Mn53_to_n_Co56",  // 1233,
        "He4_Mn53_to_p_Fe56",  // 1234,
        "n_Mn54_to_p_Cr54",  // 1235,
        "n_Mn54_to_He4_V51",  // 1236,
        "p_Mn54_to_n_Fe54",  // 1237,
        "p_Mn54_to_He4_Cr51",  // 1238,
        "He4_Mn54_to_n_Co57",  // 1239,
        "He4_Mn54_to_p_Fe57",  // 1240,
        "n_Mn55_to_He4_V52",  // 1241,
        "p_Mn55_to_n_Fe55",  // 1242,
        "p_Mn55_to_He4_Cr52",  // 1243,
        "He4_Mn55_to_n_Co58",  // 1244,
        "He4_Mn55_to_p_Fe58",  // 1245,
        "n_Fe52_to_p_Mn52",  // 1246,
        "n_Fe52_to_He4_Cr49",  // 1247,
        "He4_Fe52_to_p_Co55",  // 1248,
        "n_Fe53_to_p_Mn53",  // 1249,
        "n_Fe53_to_He4_Cr50",  // 1250,
        "p_Fe53_to_n_Co53",  // 1251,
        "p_Fe53_to_He4_Mn50",  // 1252,
        "He4_Fe53_to_n_Ni56",  // 1253,
        "He4_Fe53_to_p_Co56",  // 1254,
        "n_Fe54_to_p_Mn54",  // 1255,
        "n_Fe54_to_He4_Cr51",  // 1256,
        "p_Fe54_to_n_Co54",  // 1257,
        "p_Fe54_to_He4_Mn51",  // 1258,
        "He4_Fe54_to_n_Ni57",  // 1259,
        "He4_Fe54_to_p_Co57",  // 1260,
        "n_Fe55_to_p_Mn55",  // 1261,
        "n_Fe55_to_He4_Cr52",  // 1262,
        "p_Fe55_to_n_Co55",  // 1263,
        "p_Fe55_to_He4_Mn52",  // 1264,
        "He4_Fe55_to_n_Ni58",  // 1265,
        "He4_Fe55_to_p_Co58",  // 1266,
        "n_Fe56_to_He4_Cr53",  // 1267,
        "p_Fe56_to_n_Co56",  // 1268,
        "p_Fe56_to_He4_Mn53",  // 1269,
        "He4_Fe56_to_n_Ni59",  // 1270,
        "He4_Fe56_to_p_Co59",  // 1271,
        "n_Fe57_to_He4_Cr54",  // 1272,
        "p_Fe57_to_n_Co57",  // 1273,
        "p_Fe57_to_He4_Mn54",  // 1274,
        "He4_Fe57_to_n_Ni60",  // 1275,
        "p_Fe58_to_n_Co58",  // 1276,
        "p_Fe58_to_He4_Mn55",  // 1277,
        "He4_Fe58_to_n_Ni61",  // 1278,
        "n_Co53_to_p_Fe53",  // 1279,
        "n_Co53_to_He4_Mn50",  // 1280,
        "He4_Co53_to_p_Ni56",  // 1281,
        "n_Co54_to_p_Fe54",  // 1282,
        "n_Co54_to_He4_Mn51",  // 1283,
        "He4_Co54_to_n_Cu57",  // 1284,
        "He4_Co54_to_p_Ni57",  // 1285,
        "n_Co55_to_p_Fe55",  // 1286,
        "n_Co55_to_He4_Mn52",  // 1287,
        "p_Co55_to_He4_Fe52",  // 1288,
        "He4_Co55_to_n_Cu58",  // 1289,
        "He4_Co55_to_p_Ni58",  // 1290,
        "n_Co56_to_p_Fe56",  // 1291,
        "n_Co56_to_He4_Mn53",  // 1292,
        "p_Co56_to_n_Ni56",  // 1293,
        "p_Co56_to_He4_Fe53",  // 1294,
        "He4_Co56_to_n_Cu59",  // 1295,
        "He4_Co56_to_p_Ni59",  // 1296,
        "n_Co57_to_p_Fe57",  // 1297,
        "n_Co57_to_He4_Mn54",  // 1298,
        "p_Co57_to_n_Ni57",  // 1299,
        "p_Co57_to_He4_Fe54",  // 1300,
        "He4_Co57_to_n_Cu60",  // 1301,
        "He4_Co57_to_p_Ni60",  // 1302,
        "n_Co58_to_p_Fe58",  // 1303,
        "n_Co58_to_He4_Mn55",  // 1304,
        "p_Co58_to_n_Ni58",  // 1305,
        "p_Co58_to_He4_Fe55",  // 1306,
        "He4_Co58_to_n_Cu61",  // 1307,
        "He4_Co58_to_p_Ni61",  // 1308,
        "p_Co59_to_n_Ni59",  // 1309,
        "p_Co59_to_He4_Fe56",  // 1310,
        "He4_Co59_to_n_Cu62",  // 1311,
        "He4_Co59_to_p_Ni62",  // 1312,
        "n_Ni56_to_p_Co56",  // 1313,
        "n_Ni56_to_He4_Fe53",  // 1314,
        "p_Ni56_to_He4_Co53",  // 1315,
        "He4_Ni56_to_n_Zn59",  // 1316,
        "He4_Ni56_to_p_Cu59",  // 1317,
        "n_Ni57_to_p_Co57",  // 1318,
        "n_Ni57_to_He4_Fe54",  // 1319,
        "p_Ni57_to_n_Cu57",  // 1320,
        "p_Ni57_to_He4_Co54",  // 1321,
        "He4_Ni57_to_n_Zn60",  // 1322,
        "He4_Ni57_to_p_Cu60",  // 1323,
        "n_Ni58_to_p_Co58",  // 1324,
        "n_Ni58_to_He4_Fe55",  // 1325,
        "p_Ni58_to_n_Cu58",  // 1326,
        "p_Ni58_to_He4_Co55",  // 1327,
        "He4_Ni58_to_n_Zn61",  // 1328,
        "He4_Ni58_to_p_Cu61",  // 1329,
        "n_Ni59_to_p_Co59",  // 1330,
        "n_Ni59_to_He4_Fe56",  // 1331,
        "p_Ni59_to_n_Cu59",  // 1332,
        "p_Ni59_to_He4_Co56",  // 1333,
        "He4_Ni59_to_n_Zn62",  // 1334,
        "He4_Ni59_to_p_Cu62",  // 1335,
        "n_Ni60_to_He4_Fe57",  // 1336,
        "p_Ni60_to_n_Cu60",  // 1337,
        "p_Ni60_to_He4_Co57",  // 1338,
        "He4_Ni60_to_n_Zn63",  // 1339,
        "He4_Ni60_to_p_Cu63",  // 1340,
        "n_Ni61_to_He4_Fe58",  // 1341,
        "p_Ni61_to_n_Cu61",  // 1342,
        "p_Ni61_to_He4_Co58",  // 1343,
        "He4_Ni61_to_n_Zn64",  // 1344,
        "He4_Ni61_to_p_Cu64",  // 1345,
        "p_Ni62_to_n_Cu62",  // 1346,
        "p_Ni62_to_He4_Co59",  // 1347,
        "He4_Ni62_to_n_Zn65",  // 1348,
        "He4_Ni62_to_p_Cu65",  // 1349,
        "p_Ni63_to_n_Cu63",  // 1350,
        "He4_Ni63_to_n_Zn66",  // 1351,
        "p_Ni64_to_n_Cu64",  // 1352,
        "n_Cu57_to_p_Ni57",  // 1353,
        "n_Cu57_to_He4_Co54",  // 1354,
        "He4_Cu57_to_p_Zn60",  // 1355,
        "n_Cu58_to_p_Ni58",  // 1356,
        "n_Cu58_to_He4_Co55",  // 1357,
        "He4_Cu58_to_p_Zn61",  // 1358,
        "n_Cu59_to_p_Ni59",  // 1359,
        "n_Cu59_to_He4_Co56",  // 1360,
        "p_Cu59_to_n_Zn59",  // 1361,
        "p_Cu59_to_He4_Ni56",  // 1362,
        "He4_Cu59_to_n_Ga62",  // 1363,
        "He4_Cu59_to_p_Zn62",  // 1364,
        "n_Cu60_to_p_Ni60",  // 1365,
        "n_Cu60_to_He4_Co57",  // 1366,
        "p_Cu60_to_n_Zn60",  // 1367,
        "p_Cu60_to_He4_Ni57",  // 1368,
        "He4_Cu60_to_n_Ga63",  // 1369,
        "He4_Cu60_to_p_Zn63",  // 1370,
        "n_Cu61_to_p_Ni61",  // 1371,
        "n_Cu61_to_He4_Co58",  // 1372,
        "p_Cu61_to_n_Zn61",  // 1373,
        "p_Cu61_to_He4_Ni58",  // 1374,
        "He4_Cu61_to_n_Ga64",  // 1375,
        "He4_Cu61_to_p_Zn64",  // 1376,
        "n_Cu62_to_p_Ni62",  // 1377,
        "n_Cu62_to_He4_Co59",  // 1378,
        "p_Cu62_to_n_Zn62",  // 1379,
        "p_Cu62_to_He4_Ni59",  // 1380,
        "He4_Cu62_to_p_Zn65",  // 1381,
        "n_Cu63_to_p_Ni63",  // 1382,
        "p_Cu63_to_n_Zn63",  // 1383,
        "p_Cu63_to_He4_Ni60",  // 1384,
        "He4_Cu63_to_p_Zn66",  // 1385,
        "n_Cu64_to_p_Ni64",  // 1386,
        "p_Cu64_to_n_Zn64",  // 1387,
        "p_Cu64_to_He4_Ni61",  // 1388,
        "p_Cu65_to_n_Zn65",  // 1389,
        "p_Cu65_to_He4_Ni62",  // 1390,
        "n_Zn59_to_p_Cu59",  // 1391,
        "n_Zn59_to_He4_Ni56",  // 1392,
        "He4_Zn59_to_p_Ga62",  // 1393,
        "n_Zn60_to_p_Cu60",  // 1394,
        "n_Zn60_to_He4_Ni57",  // 1395,
        "p_Zn60_to_He4_Cu57",  // 1396,
        "He4_Zn60_to_n_Ge63",  // 1397,
        "He4_Zn60_to_p_Ga63",  // 1398,
        "n_Zn61_to_p_Cu61",  // 1399,
        "n_Zn61_to_He4_Ni58",  // 1400,
        "p_Zn61_to_He4_Cu58",  // 1401,
        "He4_Zn61_to_n_Ge64",  // 1402,
        "He4_Zn61_to_p_Ga64",  // 1403,
        "n_Zn62_to_p_Cu62",  // 1404,
        "n_Zn62_to_He4_Ni59",  // 1405,
        "p_Zn62_to_n_Ga62",  // 1406,
        "p_Zn62_to_He4_Cu59",  // 1407,
        "n_Zn63_to_p_Cu63",  // 1408,
        "n_Zn63_to_He4_Ni60",  // 1409,
        "p_Zn63_to_n_Ga63",  // 1410,
        "p_Zn63_to_He4_Cu60",  // 1411,
        "n_Zn64_to_p_Cu64",  // 1412,
        "n_Zn64_to_He4_Ni61",  // 1413,
        "p_Zn64_to_n_Ga64",  // 1414,
        "p_Zn64_to_He4_Cu61",  // 1415,
        "n_Zn65_to_p_Cu65",  // 1416,
        "n_Zn65_to_He4_Ni62",  // 1417,
        "p_Zn65_to_He4_Cu62",  // 1418,
        "n_Zn66_to_He4_Ni63",  // 1419,
        "p_Zn66_to_He4_Cu63",  // 1420,
        "n_Ga62_to_p_Zn62",  // 1421,
        "n_Ga62_to_He4_Cu59",  // 1422,
        "p_Ga62_to_He4_Zn59",  // 1423,
        "n_Ga63_to_p_Zn63",  // 1424,
        "n_Ga63_to_He4_Cu60",  // 1425,
        "p_Ga63_to_n_Ge63",  // 1426,
        "p_Ga63_to_He4_Zn60",  // 1427,
        "n_Ga64_to_p_Zn64",  // 1428,
        "n_Ga64_to_He4_Cu61",  // 1429,
        "p_Ga64_to_n_Ge64",  // 1430,
        "p_Ga64_to_He4_Zn61",  // 1431,
        "n_Ge63_to_p_Ga63",  // 1432,
        "n_Ge63_to_He4_Zn60",  // 1433,
        "n_Ge64_to_p_Ga64",  // 1434,
        "n_Ge64_to_He4_Zn61",  // 1435,
        "p_d_to_n_p_p",  // 1436,
        "He3_He3_to_p_p_He4",  // 1437,
        "d_Li7_to_n_He4_He4",  // 1438,
        "d_Be7_to_p_He4_He4",  // 1439,
        "p_Be9_to_d_He4_He4",  // 1440,
        "n_B8_to_p_He4_He4",  // 1441,
        "p_B11_to_He4_He4_He4",  // 1442,
        "He3_Li7_to_n_p_He4_He4",  // 1443,
        "He3_Be7_to_p_p_He4_He4",  // 1444,
        "p_Be9_to_n_p_He4_He4",  // 1445,
        "n_p_He4_to_Li6",  // 1446,
        "n_He4_He4_to_Be9",  // 1447,
        "He4_He4_He4_to_C12",  // 1448,
        "n_p_p_to_p_d",  // 1449,
        "p_p_He4_to_He3_He3",  // 1450,
        "n_He4_He4_to_d_Li7",  // 1451,
        "p_He4_He4_to_n_B8",  // 1452,
        "p_He4_He4_to_d_Be7",  // 1453,
        "d_He4_He4_to_p_Be9",  // 1454,
        "He4_He4_He4_to_p_B11",  // 1455,
        "n_p_He4_He4_to_He3_Li7",  // 1456,
        "n_p_He4_He4_to_p_Be9",  // 1457,
        "p_p_He4_He4_to_He3_Be7",  // 1458,
        "Al25_to_Mg25",  // 1459,
        "Al26_to_Mg26",  // 1460,
        "Ar36_to_Cl36",  // 1461,
        "Ar37_to_Cl37",  // 1462,
        "Ar37_to_K37",  // 1463,
        "Ar38_to_K38",  // 1464,
        "Ar39_to_K39",  // 1465,
        "Ar40_to_K40",  // 1466,
        "Ca40_to_K40",  // 1467,
        "Ca41_to_K41",  // 1468,
        "Ca43_to_Sc43",  // 1469,
        "Ca44_to_Sc44",  // 1470,
        "Cl33_to_S33",  // 1471,
        "Cl34_to_S34",  // 1472,
        "Cl35_to_S35",  // 1473,
        "Cl36_to_Ar36",  // 1474,
        "Cl36_to_S36",  // 1475,
        "Cl37_to_Ar37",  // 1476,
        "Co53_to_Fe53",  // 1477,
        "Cu57_to_Ni57",  // 1478,
        "Fe53_to_Co53",  // 1479,
        "K37_to_Ar37",  // 1480,
        "K38_to_Ar38",  // 1481,
        "K39_to_Ar39",  // 1482,
        "K40_to_Ar40",  // 1483,
        "K40_to_Ca40",  // 1484,
        "K41_to_Ca41",  // 1485,
        "Mg23_to_Na23",  // 1486,
        "Mg25_to_Al25",  // 1487,
        "Mg26_to_Al26",  // 1488,
        "Na21_to_Ne21",  // 1489,
        "Na22_to_Ne22",  // 1490,
        "Na23_to_Mg23",  // 1491,
        "Ne21_to_Na21",  // 1492,
        "Ne22_to_Na22",  // 1493,
        "Ni57_to_Cu57",  // 1494,
        "P29_to_Si29",  // 1495,
        "P30_to_Si30",  // 1496,
        "P31_to_Si31",  // 1497,
        "P32_to_S32",  // 1498,
        "P32_to_Si32",  // 1499,
        "P33_to_S33",  // 1500,
        "S32_to_P32",  // 1501,
        "S33_to_Cl33",  // 1502,
        "S33_to_P33",  // 1503,
        "S34_to_Cl34",  // 1504,
        "S35_to_Cl35",  // 1505,
        "S36_to_Cl36",  // 1506,
        "Sc43_to_Ca43",  // 1507,
        "Sc44_to_Ca44",  // 1508,
        "Sc44_to_Ti44",  // 1509,
        "Si29_to_P29",  // 1510,
        "Si30_to_P30",  // 1511,
        "Si31_to_P31",  // 1512,
        "Si32_to_P32",  // 1513,
        "Ti44_to_Sc44",  // 1514,
        "Ca45_to_Sc45",  // 1515,
        "Ca46_to_Sc46",  // 1516,
        "Ca47_to_Sc47",  // 1517,
        "Ca48_to_Sc48",  // 1518,
        "Co54_to_Fe54",  // 1519,
        "Co55_to_Fe55",  // 1520,
        "Co56_to_Fe56",  // 1521,
        "Co56_to_Ni56",  // 1522,
        "Co57_to_Fe57",  // 1523,
        "Co57_to_Ni57",  // 1524,
        "Co58_to_Fe58",  // 1525,
        "Co58_to_Ni58",  // 1526,
        "Co59_to_Ni59",  // 1527,
        "Cr48_to_V48",  // 1528,
        "Cr49_to_V49",  // 1529,
        "Cr50_to_Mn50",  // 1530,
        "Cr50_to_V50",  // 1531,
        "Cr51_to_Mn51",  // 1532,
        "Cr51_to_V51",  // 1533,
        "Cr52_to_Mn52",  // 1534,
        "Cr52_to_V52",  // 1535,
        "Cr53_to_Mn53",  // 1536,
        "Cr54_to_Mn54",  // 1537,
        "Cu58_to_Ni58",  // 1538,
        "Cu59_to_Ni59",  // 1539,
        "Cu60_to_Ni60",  // 1540,
        "Cu60_to_Zn60",  // 1541,
        "Cu61_to_Ni61",  // 1542,
        "Cu61_to_Zn61",  // 1543,
        "Cu62_to_Ni62",  // 1544,
        "Cu62_to_Zn62",  // 1545,
        "Cu63_to_Ni63",  // 1546,
        "Cu63_to_Zn63",  // 1547,
        "Cu64_to_Ni64",  // 1548,
        "Cu64_to_Zn64",  // 1549,
        "Cu65_to_Zn65",  // 1550,
        "Fe52_to_Mn52",  // 1551,
        "Fe53_to_Mn53",  // 1552,
        "Fe54_to_Co54",  // 1553,
        "Fe54_to_Mn54",  // 1554,
        "Fe55_to_Co55",  // 1555,
        "Fe55_to_Mn55",  // 1556,
        "Fe56_to_Co56",  // 1557,
        "Fe57_to_Co57",  // 1558,
        "Fe58_to_Co58",  // 1559,
        "Ga62_to_Zn62",  // 1560,
        "Ga63_to_Zn63",  // 1561,
        "Ga64_to_Ge64",  // 1562,
        "Ga64_to_Zn64",  // 1563,
        "Ge64_to_Ga64",  // 1564,
        "Mn50_to_Cr50",  // 1565,
        "Mn51_to_Cr51",  // 1566,
        "Mn52_to_Cr52",  // 1567,
        "Mn52_to_Fe52",  // 1568,
        "Mn53_to_Cr53",  // 1569,
        "Mn53_to_Fe53",  // 1570,
        "Mn54_to_Cr54",  // 1571,
        "Mn54_to_Fe54",  // 1572,
        "Mn55_to_Fe55",  // 1573,
        "n_to_p",  // 1574,
        "Ni56_to_Co56",  // 1575,
        "Ni57_to_Co57",  // 1576,
        "Ni58_to_Co58",  // 1577,
        "Ni58_to_Cu58",  // 1578,
        "Ni59_to_Co59",  // 1579,
        "Ni59_to_Cu59",  // 1580,
        "Ni60_to_Cu60",  // 1581,
        "Ni61_to_Cu61",  // 1582,
        "Ni62_to_Cu62",  // 1583,
        "Ni63_to_Cu63",  // 1584,
        "Ni64_to_Cu64",  // 1585,
        "p_to_n",  // 1586,
        "Sc45_to_Ca45",  // 1587,
        "Sc45_to_Ti45",  // 1588,
        "Sc46_to_Ca46",  // 1589,
        "Sc46_to_Ti46",  // 1590,
        "Sc47_to_Ca47",  // 1591,
        "Sc47_to_Ti47",  // 1592,
        "Sc48_to_Ca48",  // 1593,
        "Sc48_to_Ti48",  // 1594,
        "Sc49_to_Ti49",  // 1595,
        "Ti45_to_Sc45",  // 1596,
        "Ti46_to_Sc46",  // 1597,
        "Ti46_to_V46",  // 1598,
        "Ti47_to_Sc47",  // 1599,
        "Ti47_to_V47",  // 1600,
        "Ti48_to_Sc48",  // 1601,
        "Ti48_to_V48",  // 1602,
        "Ti49_to_Sc49",  // 1603,
        "Ti49_to_V49",  // 1604,
        "Ti50_to_V50",  // 1605,
        "Ti51_to_V51",  // 1606,
        "V46_to_Ti46",  // 1607,
        "V47_to_Ti47",  // 1608,
        "V48_to_Cr48",  // 1609,
        "V48_to_Ti48",  // 1610,
        "V49_to_Cr49",  // 1611,
        "V49_to_Ti49",  // 1612,
        "V50_to_Cr50",  // 1613,
        "V50_to_Ti50",  // 1614,
        "V51_to_Cr51",  // 1615,
        "V51_to_Ti51",  // 1616,
        "V52_to_Cr52",  // 1617,
        "Zn60_to_Cu60",  // 1618,
        "Zn61_to_Cu61",  // 1619,
        "Zn62_to_Cu62",  // 1620,
        "Zn62_to_Ga62",  // 1621,
        "Zn63_to_Cu63",  // 1622,
        "Zn63_to_Ga63",  // 1623,
        "Zn64_to_Cu64",  // 1624,
        "Zn64_to_Ga64",  // 1625,
        "Zn65_to_Cu65"  // 1626,
    };

}

#ifdef NSE_NET
namespace NSE_INDEX
{
    constexpr int H1_index = 1;
    constexpr int N_index = 0;
    constexpr int He4_index = 4;

    // Each row corresponds to the rate in NetworkRates enum
    // First 3 row indices for reactants, followed by 3 product indices
    // last index is the corresponding reverse rate index.

    extern AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices;
}
#endif

#endif
