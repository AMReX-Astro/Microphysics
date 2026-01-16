#ifndef REACT_UTIL_H
#define REACT_UTIL_H

#include <network.H>
#include <extern_parameters.H>

using namespace amrex::literals;

struct init_t {
    int nprim;
    int is1;
    int is2;
    int is3;

    int n1;
    int n2;
    int n3;

    amrex::Real Xp_min;
    amrex::Real Xp_max;
};

AMREX_INLINE
init_t setup_composition(const int nz) {

    // get the primary species indices
    init_t comp_data;

    comp_data.nprim = 0;

    // we absolutely require primary_species_1 to be defined
    comp_data.is1 = network_spec_index(unit_test_rp::primary_species_1);
    if (comp_data.is1 >= 0) {
        comp_data.nprim++;
    } else {
        amrex::Error("Error: invalid primary_species_1");
    }

    // we'll check the next species, but if it is not valid,
    // we'll just move on
    comp_data.is2 = network_spec_index(unit_test_rp::primary_species_2);
    if (comp_data.is2 >= 0) {
        comp_data.nprim++;

        // only consider primary_species_3 if primary_species_2
        // was defined
        comp_data.is3 = network_spec_index(unit_test_rp::primary_species_3);
        if (comp_data.is3 >= 0) {
            comp_data.nprim++;
        }
    }

    if (comp_data.nprim == 0) {
        amrex::Error("ERROR: no primary species set");
    }

    // figure out how many zones to allocate to the each of the primary
    // species and the extrema for the primary species

    if (comp_data.nprim == 1) {
        comp_data.n1 = nz;
        comp_data.n2 = 0;
        comp_data.n3 = 0;

        comp_data.Xp_min = 0.2_rt;
        comp_data.Xp_max = 0.9_rt;

    } else if (comp_data.nprim == 2) {
        comp_data.n1 = nz/2;
        comp_data.n2 = nz - comp_data.n1;
        comp_data.n3 = 0;

        comp_data.Xp_min = 0.2_rt;
        comp_data.Xp_max = 0.8_rt;

    } else if (comp_data.nprim == 3) {
        comp_data.n1 = nz/3;
        comp_data.n2 = nz/3;
        comp_data.n3 = nz - comp_data.n1 - comp_data.n2;

        comp_data.Xp_min = 0.2_rt;
        comp_data.Xp_max = 0.7_rt;

    }

    return comp_data;
}


AMREX_INLINE AMREX_GPU_HOST_DEVICE
void get_xn(const int k, const init_t& cd, amrex::Real *xn_zone, bool uniform_composition=false) {

    for (int n = 0; n < NumSpec; n++) {
        xn_zone[n] = 0.0_rt;
    }

    if (uniform_composition) {

        for (int n = 0; n < NumSpec; ++n) {
            xn_zone[n] = 1.0_rt / NumSpec;
        }

    }
    else {

        if (k < cd.n1) {
            if (cd.nprim >= 2) xn_zone[cd.is2] = cd.Xp_min/2;
            if (cd.nprim >= 3) xn_zone[cd.is3] = cd.Xp_min/2;

            amrex::Real dX = (cd.Xp_max - cd.Xp_min) /
                amrex::max((cd.n1 - 1), 1);
            xn_zone[cd.is1] = cd.Xp_min + k * dX;

        } else if (cd.nprim >= 2 && k < cd.n1 + cd.n2) {
            xn_zone[cd.is1] = cd.Xp_min/2;
            if (cd.nprim >= 3) {
                xn_zone[cd.is3] = cd.Xp_min/2;
            }
            amrex::Real dX = (cd.Xp_max - cd.Xp_min) /
                amrex::max((cd.n2 - 1), 1);
            xn_zone[cd.is2] = cd.Xp_min + (k - cd.n1) * dX;

        } else {
            xn_zone[cd.is1] = cd.Xp_min/2;
            xn_zone[cd.is2] = cd.Xp_min/2;

            amrex::Real dX = (cd.Xp_max - cd.Xp_min) /
                amrex::max((cd.n3 - 1), 1);
            xn_zone[cd.is3] = cd.Xp_min + (k - (cd.n1 + cd.n2)) * dX;

        }

        amrex::Real excess = 0.0_rt;
        for (int n = 0; n < NumSpec; n++) {
            excess += xn_zone[n];
        }
        excess = 1.0_rt - excess;

        for (int n = 0; n < NumSpec; n++) {
            if (n == cd.is1 ||
                (cd.nprim >= 2 && n == cd.is2) ||
                (cd.nprim >= 3 && n == cd.is3)) {
                continue;
            }

            xn_zone[n] = excess / (NumSpec - cd.nprim);
        }

    }


    // normalize -- just in case

    amrex::Real sum_X = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        sum_X += xn_zone[n];
    }
    for (int n = 0; n < NumSpec; n++) {
        xn_zone[n] = xn_zone[n] / sum_X;
    }

}

///
/// return the value of the runtime parameter XN, where N is
/// the integer index passed in
///

AMREX_INLINE AMREX_GPU_HOST_DEVICE
amrex::Real get_xn(const int index, bool uniform_composition=false) {

    amrex::Real mass_fraction{};

    if (uniform_composition) {
        mass_fraction = 1.0_rt / NumSpec;
        return mass_fraction;
    }

    switch (index) {

    case 1:
        mass_fraction = unit_test_rp::X1;
        break;
    case 2:
        mass_fraction = unit_test_rp::X2;
        break;
    case 3:
        mass_fraction = unit_test_rp::X3;
        break;
    case 4:
        mass_fraction = unit_test_rp::X4;
        break;
    case 5:
        mass_fraction = unit_test_rp::X5;
        break;
    case 6:
        mass_fraction = unit_test_rp::X6;
        break;
    case 7:
        mass_fraction = unit_test_rp::X7;
        break;
    case 8:
        mass_fraction = unit_test_rp::X8;
        break;
    case 9:
        mass_fraction = unit_test_rp::X9;
        break;
    case 10:
        mass_fraction = unit_test_rp::X10;
        break;
    case 11:
        mass_fraction = unit_test_rp::X11;
        break;
    case 12:
        mass_fraction = unit_test_rp::X12;
        break;
    case 13:
        mass_fraction = unit_test_rp::X13;
        break;
    case 14:
        mass_fraction = unit_test_rp::X14;
        break;
    case 15:
        mass_fraction = unit_test_rp::X15;
        break;
    case 16:
        mass_fraction = unit_test_rp::X16;
        break;
    case 17:
        mass_fraction = unit_test_rp::X17;
        break;
    case 18:
        mass_fraction = unit_test_rp::X18;
        break;
    case 19:
        mass_fraction = unit_test_rp::X19;
        break;
    case 20:
        mass_fraction = unit_test_rp::X20;
        break;
    case 21:
        mass_fraction = unit_test_rp::X21;
        break;
    case 22:
        mass_fraction = unit_test_rp::X22;
        break;
    case 23:
        mass_fraction = unit_test_rp::X23;
        break;
    case 24:
        mass_fraction = unit_test_rp::X24;
        break;
    case 25:
        mass_fraction = unit_test_rp::X25;
        break;
    case 26:
        mass_fraction = unit_test_rp::X26;
        break;
    case 27:
        mass_fraction = unit_test_rp::X27;
        break;
    case 28:
        mass_fraction = unit_test_rp::X28;
        break;
    case 29:
        mass_fraction = unit_test_rp::X29;
        break;
    case 30:
        mass_fraction = unit_test_rp::X30;
        break;
    case 31:
        mass_fraction = unit_test_rp::X31;
        break;
    case 32:
        mass_fraction = unit_test_rp::X32;
        break;
    case 33:
        mass_fraction = unit_test_rp::X33;
        break;
    case 34:
        mass_fraction = unit_test_rp::X34;
        break;
    case 35:
        mass_fraction = unit_test_rp::X35;
        break;
    case 36:
        mass_fraction = unit_test_rp::X36;
        break;
    case 37:
        mass_fraction = unit_test_rp::X37;
        break;
    case 38:
        mass_fraction = unit_test_rp::X38;
        break;
    case 39:
        mass_fraction = unit_test_rp::X39;
        break;
    case 40:
        mass_fraction = unit_test_rp::X40;
        break;
    case 41:
        mass_fraction = unit_test_rp::X41;
        break;
    case 42:
        mass_fraction = unit_test_rp::X42;
        break;
    case 43:
        mass_fraction = unit_test_rp::X43;
        break;
    case 44:
        mass_fraction = unit_test_rp::X44;
        break;
    case 45:
        mass_fraction = unit_test_rp::X45;
        break;
    case 46:
        mass_fraction = unit_test_rp::X46;
        break;
    case 47:
        mass_fraction = unit_test_rp::X47;
        break;
    case 48:
        mass_fraction = unit_test_rp::X48;
        break;
    case 49:
        mass_fraction = unit_test_rp::X49;
        break;
    case 50:
        mass_fraction = unit_test_rp::X50;
        break;
    case 51:
        mass_fraction = unit_test_rp::X51;
        break;
    case 52:
        mass_fraction = unit_test_rp::X52;
        break;
    case 53:
        mass_fraction = unit_test_rp::X53;
        break;
    case 54:
        mass_fraction = unit_test_rp::X54;
        break;
    case 55:
        mass_fraction = unit_test_rp::X55;
        break;
    case 56:
        mass_fraction = unit_test_rp::X56;
        break;
    case 57:
        mass_fraction = unit_test_rp::X57;
        break;
    case 58:
        mass_fraction = unit_test_rp::X58;
        break;
    case 59:
        mass_fraction = unit_test_rp::X59;
        break;
    case 60:
        mass_fraction = unit_test_rp::X60;
        break;
    case 61:
        mass_fraction = unit_test_rp::X61;
        break;
    case 62:
        mass_fraction = unit_test_rp::X62;
        break;
    case 63:
        mass_fraction = unit_test_rp::X63;
        break;
    case 64:
        mass_fraction = unit_test_rp::X64;
        break;
    case 65:
        mass_fraction = unit_test_rp::X65;
        break;
    case 66:
        mass_fraction = unit_test_rp::X66;
        break;
    case 67:
        mass_fraction = unit_test_rp::X67;
        break;
    case 68:
        mass_fraction = unit_test_rp::X68;
        break;
    case 69:
        mass_fraction = unit_test_rp::X69;
        break;
    case 70:
        mass_fraction = unit_test_rp::X70;
        break;
    case 71:
        mass_fraction = unit_test_rp::X71;
        break;
    case 72:
        mass_fraction = unit_test_rp::X72;
        break;
    case 73:
        mass_fraction = unit_test_rp::X73;
        break;
    case 74:
        mass_fraction = unit_test_rp::X74;
        break;
    case 75:
        mass_fraction = unit_test_rp::X75;
        break;
    case 76:
        mass_fraction = unit_test_rp::X76;
        break;
    case 77:
        mass_fraction = unit_test_rp::X77;
        break;
    case 78:
        mass_fraction = unit_test_rp::X78;
        break;
    case 79:
        mass_fraction = unit_test_rp::X79;
        break;
    case 80:
        mass_fraction = unit_test_rp::X80;
        break;
    case 81:
        mass_fraction = unit_test_rp::X81;
        break;
    case 82:
        mass_fraction = unit_test_rp::X82;
        break;
    case 83:
        mass_fraction = unit_test_rp::X83;
        break;
    case 84:
        mass_fraction = unit_test_rp::X84;
        break;
    case 85:
        mass_fraction = unit_test_rp::X85;
        break;
    case 86:
        mass_fraction = unit_test_rp::X86;
        break;
    case 87:
        mass_fraction = unit_test_rp::X87;
        break;
    case 88:
        mass_fraction = unit_test_rp::X88;
        break;
    case 89:
        mass_fraction = unit_test_rp::X89;
        break;
    case 90:
        mass_fraction = unit_test_rp::X90;
        break;
    case 91:
        mass_fraction = unit_test_rp::X91;
        break;
    case 92:
        mass_fraction = unit_test_rp::X92;
        break;
    case 93:
        mass_fraction = unit_test_rp::X93;
        break;
    case 94:
        mass_fraction = unit_test_rp::X94;
        break;
    case 95:
        mass_fraction = unit_test_rp::X95;
        break;
    case 96:
        mass_fraction = unit_test_rp::X96;
        break;
    case 97:
        mass_fraction = unit_test_rp::X97;
        break;
    case 98:
        mass_fraction = unit_test_rp::X98;
        break;
    case 99:
        mass_fraction = unit_test_rp::X99;
        break;
    case 100:
        mass_fraction = unit_test_rp::X100;
        break;
    case 101:
        mass_fraction = unit_test_rp::X101;
        break;
    case 102:
        mass_fraction = unit_test_rp::X102;
        break;
    case 103:
        mass_fraction = unit_test_rp::X103;
        break;
    case 104:
        mass_fraction = unit_test_rp::X104;
        break;
    case 105:
        mass_fraction = unit_test_rp::X105;
        break;
    case 106:
        mass_fraction = unit_test_rp::X106;
        break;
    case 107:
        mass_fraction = unit_test_rp::X107;
        break;
    case 108:
        mass_fraction = unit_test_rp::X108;
        break;
    case 109:
        mass_fraction = unit_test_rp::X109;
        break;
    case 110:
        mass_fraction = unit_test_rp::X110;
        break;
    case 111:
        mass_fraction = unit_test_rp::X111;
        break;
    case 112:
        mass_fraction = unit_test_rp::X112;
        break;
    case 113:
        mass_fraction = unit_test_rp::X113;
        break;
    case 114:
        mass_fraction = unit_test_rp::X114;
        break;
    case 115:
        mass_fraction = unit_test_rp::X115;
        break;
    case 116:
        mass_fraction = unit_test_rp::X116;
        break;
    case 117:
        mass_fraction = unit_test_rp::X117;
        break;
    case 118:
        mass_fraction = unit_test_rp::X118;
        break;
    case 119:
        mass_fraction = unit_test_rp::X119;
        break;
    case 120:
        mass_fraction = unit_test_rp::X120;
        break;
    case 121:
        mass_fraction = unit_test_rp::X121;
        break;
    case 122:
        mass_fraction = unit_test_rp::X122;
        break;
    case 123:
        mass_fraction = unit_test_rp::X123;
        break;
    case 124:
        mass_fraction = unit_test_rp::X124;
        break;
    case 125:
        mass_fraction = unit_test_rp::X125;
        break;
    case 126:
        mass_fraction = unit_test_rp::X126;
        break;
    case 127:
        mass_fraction = unit_test_rp::X127;
        break;
    case 128:
        mass_fraction = unit_test_rp::X128;
        break;
    case 129:
        mass_fraction = unit_test_rp::X129;
        break;
    case 130:
        mass_fraction = unit_test_rp::X130;
        break;
    case 131:
        mass_fraction = unit_test_rp::X131;
        break;
    case 132:
        mass_fraction = unit_test_rp::X132;
        break;
    case 133:
        mass_fraction = unit_test_rp::X133;
        break;
    case 134:
        mass_fraction = unit_test_rp::X134;
        break;
    case 135:
        mass_fraction = unit_test_rp::X135;
        break;
    case 136:
        mass_fraction = unit_test_rp::X136;
        break;
    case 137:
        mass_fraction = unit_test_rp::X137;
        break;
    case 138:
        mass_fraction = unit_test_rp::X138;
        break;
    case 139:
        mass_fraction = unit_test_rp::X139;
        break;
    case 140:
        mass_fraction = unit_test_rp::X140;
        break;
    case 141:
        mass_fraction = unit_test_rp::X141;
        break;
    case 142:
        mass_fraction = unit_test_rp::X142;
        break;
    case 143:
        mass_fraction = unit_test_rp::X143;
        break;
    case 144:
        mass_fraction = unit_test_rp::X144;
        break;
    case 145:
        mass_fraction = unit_test_rp::X145;
        break;
    case 146:
        mass_fraction = unit_test_rp::X146;
        break;
    case 147:
        mass_fraction = unit_test_rp::X147;
        break;
    case 148:
        mass_fraction = unit_test_rp::X148;
        break;
    case 149:
        mass_fraction = unit_test_rp::X149;
        break;
    case 150:
        mass_fraction = unit_test_rp::X150;
        break;
    case 151:
        mass_fraction = unit_test_rp::X151;
        break;
    case 152:
        mass_fraction = unit_test_rp::X152;
        break;
    case 153:
        mass_fraction = unit_test_rp::X153;
        break;
    case 154:
        mass_fraction = unit_test_rp::X154;
        break;
    case 155:
        mass_fraction = unit_test_rp::X155;
        break;
    case 156:
        mass_fraction = unit_test_rp::X156;
        break;
    case 157:
        mass_fraction = unit_test_rp::X157;
        break;
    case 158:
        mass_fraction = unit_test_rp::X158;
        break;
    case 159:
        mass_fraction = unit_test_rp::X159;
        break;
    case 160:
        mass_fraction = unit_test_rp::X160;
        break;
    case 161:
        mass_fraction = unit_test_rp::X161;
        break;
    case 162:
        mass_fraction = unit_test_rp::X162;
        break;
    case 163:
        mass_fraction = unit_test_rp::X163;
        break;
    case 164:
        mass_fraction = unit_test_rp::X164;
        break;
    case 165:
        mass_fraction = unit_test_rp::X165;
        break;
    case 166:
        mass_fraction = unit_test_rp::X166;
        break;
    case 167:
        mass_fraction = unit_test_rp::X167;
        break;
    case 168:
        mass_fraction = unit_test_rp::X168;
        break;
    case 169:
        mass_fraction = unit_test_rp::X169;
        break;
    case 170:
        mass_fraction = unit_test_rp::X170;
        break;
    case 171:
        mass_fraction = unit_test_rp::X171;
        break;
    case 172:
        mass_fraction = unit_test_rp::X172;
        break;
    case 173:
        mass_fraction = unit_test_rp::X173;
        break;
    case 174:
        mass_fraction = unit_test_rp::X174;
        break;
    case 175:
        mass_fraction = unit_test_rp::X175;
        break;
    case 176:
        mass_fraction = unit_test_rp::X176;
        break;
    case 177:
        mass_fraction = unit_test_rp::X177;
        break;
    case 178:
        mass_fraction = unit_test_rp::X178;
        break;
    case 179:
        mass_fraction = unit_test_rp::X179;
        break;
    case 180:
        mass_fraction = unit_test_rp::X180;
        break;
    case 181:
        mass_fraction = unit_test_rp::X181;
        break;
    case 182:
        mass_fraction = unit_test_rp::X182;
        break;
    case 183:
        mass_fraction = unit_test_rp::X183;
        break;
    case 184:
        mass_fraction = unit_test_rp::X184;
        break;
    case 185:
        mass_fraction = unit_test_rp::X185;
        break;
    case 186:
        mass_fraction = unit_test_rp::X186;
        break;
    case 187:
        mass_fraction = unit_test_rp::X187;
        break;
    case 188:
        mass_fraction = unit_test_rp::X188;
        break;
    case 189:
        mass_fraction = unit_test_rp::X189;
        break;
    case 190:
        mass_fraction = unit_test_rp::X190;
        break;
    case 191:
        mass_fraction = unit_test_rp::X191;
        break;
    case 192:
        mass_fraction = unit_test_rp::X192;
        break;
    case 193:
        mass_fraction = unit_test_rp::X193;
        break;
    case 194:
        mass_fraction = unit_test_rp::X194;
        break;
    case 195:
        mass_fraction = unit_test_rp::X195;
        break;
    case 196:
        mass_fraction = unit_test_rp::X196;
        break;
    case 197:
        mass_fraction = unit_test_rp::X197;
        break;
    case 198:
        mass_fraction = unit_test_rp::X198;
        break;
    case 199:
        mass_fraction = unit_test_rp::X199;
        break;
    case 200:
        mass_fraction = unit_test_rp::X200;
        break;
    case 201:
        mass_fraction = unit_test_rp::X201;
        break;
    case 202:
        mass_fraction = unit_test_rp::X202;
        break;
    case 203:
        mass_fraction = unit_test_rp::X203;
        break;
    case 204:
        mass_fraction = unit_test_rp::X204;
        break;
    case 205:
        mass_fraction = unit_test_rp::X205;
        break;
    case 206:
        mass_fraction = unit_test_rp::X206;
        break;
    case 207:
        mass_fraction = unit_test_rp::X207;
        break;
    case 208:
        mass_fraction = unit_test_rp::X208;
        break;
    case 209:
        mass_fraction = unit_test_rp::X209;
        break;
    case 210:
        mass_fraction = unit_test_rp::X210;
        break;
    case 211:
        mass_fraction = unit_test_rp::X211;
        break;
    case 212:
        mass_fraction = unit_test_rp::X212;
        break;
    case 213:
        mass_fraction = unit_test_rp::X213;
        break;
    case 214:
        mass_fraction = unit_test_rp::X214;
        break;
    case 215:
        mass_fraction = unit_test_rp::X215;
        break;
    case 216:
        mass_fraction = unit_test_rp::X216;
        break;
    case 217:
        mass_fraction = unit_test_rp::X217;
        break;
    case 218:
        mass_fraction = unit_test_rp::X218;
        break;
    case 219:
        mass_fraction = unit_test_rp::X219;
        break;
    case 220:
        mass_fraction = unit_test_rp::X220;
        break;
    case 221:
        mass_fraction = unit_test_rp::X221;
        break;
    case 222:
        mass_fraction = unit_test_rp::X222;
        break;
    case 223:
        mass_fraction = unit_test_rp::X223;
        break;
    case 224:
        mass_fraction = unit_test_rp::X224;
        break;
    case 225:
        mass_fraction = unit_test_rp::X225;
        break;
    case 226:
        mass_fraction = unit_test_rp::X226;
        break;
    case 227:
        mass_fraction = unit_test_rp::X227;
        break;
    case 228:
        mass_fraction = unit_test_rp::X228;
        break;
    case 229:
        mass_fraction = unit_test_rp::X229;
        break;
    case 230:
        mass_fraction = unit_test_rp::X230;
        break;
    case 231:
        mass_fraction = unit_test_rp::X231;
        break;
    case 232:
        mass_fraction = unit_test_rp::X232;
        break;
    case 233:
        mass_fraction = unit_test_rp::X233;
        break;
    case 234:
        mass_fraction = unit_test_rp::X234;
        break;
    case 235:
        mass_fraction = unit_test_rp::X235;
        break;
    case 236:
        mass_fraction = unit_test_rp::X236;
        break;
    case 237:
        mass_fraction = unit_test_rp::X237;
        break;
    case 238:
        mass_fraction = unit_test_rp::X238;
        break;
    case 239:
        mass_fraction = unit_test_rp::X239;
        break;
    case 240:
        mass_fraction = unit_test_rp::X240;
        break;
    case 241:
        mass_fraction = unit_test_rp::X241;
        break;
    case 242:
        mass_fraction = unit_test_rp::X242;
        break;
    case 243:
        mass_fraction = unit_test_rp::X243;
        break;
    case 244:
        mass_fraction = unit_test_rp::X244;
        break;
    case 245:
        mass_fraction = unit_test_rp::X245;
        break;
    case 246:
        mass_fraction = unit_test_rp::X246;
        break;
    case 247:
        mass_fraction = unit_test_rp::X247;
        break;
    case 248:
        mass_fraction = unit_test_rp::X248;
        break;
    case 249:
        mass_fraction = unit_test_rp::X249;
        break;
    case 250:
        mass_fraction = unit_test_rp::X250;
        break;
    case 251:
        mass_fraction = unit_test_rp::X251;
        break;
    case 252:
        mass_fraction = unit_test_rp::X252;
        break;
    case 253:
        mass_fraction = unit_test_rp::X253;
        break;
    case 254:
        mass_fraction = unit_test_rp::X254;
        break;
    case 255:
        mass_fraction = unit_test_rp::X255;
        break;
    case 256:
        mass_fraction = unit_test_rp::X256;
        break;
    case 257:
        mass_fraction = unit_test_rp::X257;
        break;
    case 258:
        mass_fraction = unit_test_rp::X258;
        break;
    case 259:
        mass_fraction = unit_test_rp::X259;
        break;
    case 260:
        mass_fraction = unit_test_rp::X260;
        break;
    case 261:
        mass_fraction = unit_test_rp::X261;
        break;
    case 262:
        mass_fraction = unit_test_rp::X262;
        break;
    case 263:
        mass_fraction = unit_test_rp::X263;
        break;
    case 264:
        mass_fraction = unit_test_rp::X264;
        break;
    case 265:
        mass_fraction = unit_test_rp::X265;
        break;
    case 266:
        mass_fraction = unit_test_rp::X266;
        break;
    case 267:
        mass_fraction = unit_test_rp::X267;
        break;
    case 268:
        mass_fraction = unit_test_rp::X268;
        break;
    case 269:
        mass_fraction = unit_test_rp::X269;
        break;
    case 270:
        mass_fraction = unit_test_rp::X270;
        break;
    case 271:
        mass_fraction = unit_test_rp::X271;
        break;
    case 272:
        mass_fraction = unit_test_rp::X272;
        break;
    case 273:
        mass_fraction = unit_test_rp::X273;
        break;
    case 274:
        mass_fraction = unit_test_rp::X274;
        break;
    case 275:
        mass_fraction = unit_test_rp::X275;
        break;
    case 276:
        mass_fraction = unit_test_rp::X276;
        break;
    case 277:
        mass_fraction = unit_test_rp::X277;
        break;
    case 278:
        mass_fraction = unit_test_rp::X278;
        break;
    case 279:
        mass_fraction = unit_test_rp::X279;
        break;
    case 280:
        mass_fraction = unit_test_rp::X280;
        break;
    case 281:
        mass_fraction = unit_test_rp::X281;
        break;
    case 282:
        mass_fraction = unit_test_rp::X282;
        break;
    case 283:
        mass_fraction = unit_test_rp::X283;
        break;
    case 284:
        mass_fraction = unit_test_rp::X284;
        break;
    case 285:
        mass_fraction = unit_test_rp::X285;
        break;
    case 286:
        mass_fraction = unit_test_rp::X286;
        break;
    case 287:
        mass_fraction = unit_test_rp::X287;
        break;
    case 288:
        mass_fraction = unit_test_rp::X288;
        break;
    case 289:
        mass_fraction = unit_test_rp::X289;
        break;
    case 290:
        mass_fraction = unit_test_rp::X290;
        break;
    case 291:
        mass_fraction = unit_test_rp::X291;
        break;
    case 292:
        mass_fraction = unit_test_rp::X292;
        break;
    case 293:
        mass_fraction = unit_test_rp::X293;
        break;
    case 294:
        mass_fraction = unit_test_rp::X294;
        break;
    case 295:
        mass_fraction = unit_test_rp::X295;
        break;
    case 296:
        mass_fraction = unit_test_rp::X296;
        break;
    case 297:
        mass_fraction = unit_test_rp::X297;
        break;
    case 298:
        mass_fraction = unit_test_rp::X298;
        break;
    case 299:
        mass_fraction = unit_test_rp::X299;
        break;
    case 300:
        mass_fraction = unit_test_rp::X300;
        break;
    case 301:
        mass_fraction = unit_test_rp::X301;
        break;
    case 302:
        mass_fraction = unit_test_rp::X302;
        break;
    case 303:
        mass_fraction = unit_test_rp::X303;
        break;
    case 304:
        mass_fraction = unit_test_rp::X304;
        break;
    case 305:
        mass_fraction = unit_test_rp::X305;
        break;
    case 306:
        mass_fraction = unit_test_rp::X306;
        break;
    case 307:
        mass_fraction = unit_test_rp::X307;
        break;
    case 308:
        mass_fraction = unit_test_rp::X308;
        break;
    case 309:
        mass_fraction = unit_test_rp::X309;
        break;
    case 310:
        mass_fraction = unit_test_rp::X310;
        break;
    case 311:
        mass_fraction = unit_test_rp::X311;
        break;
    case 312:
        mass_fraction = unit_test_rp::X312;
        break;
    case 313:
        mass_fraction = unit_test_rp::X313;
        break;
    case 314:
        mass_fraction = unit_test_rp::X314;
        break;
    case 315:
        mass_fraction = unit_test_rp::X315;
        break;
    case 316:
        mass_fraction = unit_test_rp::X316;
        break;
    case 317:
        mass_fraction = unit_test_rp::X317;
        break;
    case 318:
        mass_fraction = unit_test_rp::X318;
        break;
    case 319:
        mass_fraction = unit_test_rp::X319;
        break;
    case 320:
        mass_fraction = unit_test_rp::X320;
        break;
    case 321:
        mass_fraction = unit_test_rp::X321;
        break;
    case 322:
        mass_fraction = unit_test_rp::X322;
        break;
    case 323:
        mass_fraction = unit_test_rp::X323;
        break;
    case 324:
        mass_fraction = unit_test_rp::X324;
        break;
    case 325:
        mass_fraction = unit_test_rp::X325;
        break;
    case 326:
        mass_fraction = unit_test_rp::X326;
        break;
    case 327:
        mass_fraction = unit_test_rp::X327;
        break;
    case 328:
        mass_fraction = unit_test_rp::X328;
        break;
    case 329:
        mass_fraction = unit_test_rp::X329;
        break;
    case 330:
        mass_fraction = unit_test_rp::X330;
        break;
    case 331:
        mass_fraction = unit_test_rp::X331;
        break;
    case 332:
        mass_fraction = unit_test_rp::X332;
        break;
    case 333:
        mass_fraction = unit_test_rp::X333;
        break;
    case 334:
        mass_fraction = unit_test_rp::X334;
        break;
    case 335:
        mass_fraction = unit_test_rp::X335;
        break;
    case 336:
        mass_fraction = unit_test_rp::X336;
        break;
    case 337:
        mass_fraction = unit_test_rp::X337;
        break;
    case 338:
        mass_fraction = unit_test_rp::X338;
        break;
    case 339:
        mass_fraction = unit_test_rp::X339;
        break;
    case 340:
        mass_fraction = unit_test_rp::X340;
        break;
    case 341:
        mass_fraction = unit_test_rp::X341;
        break;
    case 342:
        mass_fraction = unit_test_rp::X342;
        break;
    case 343:
        mass_fraction = unit_test_rp::X343;
        break;
    case 344:
        mass_fraction = unit_test_rp::X344;
        break;
    case 345:
        mass_fraction = unit_test_rp::X345;
        break;
    case 346:
        mass_fraction = unit_test_rp::X346;
        break;
    case 347:
        mass_fraction = unit_test_rp::X347;
        break;
    case 348:
        mass_fraction = unit_test_rp::X348;
        break;
    case 349:
        mass_fraction = unit_test_rp::X349;
        break;
    case 350:
        mass_fraction = unit_test_rp::X350;
        break;
    case 351:
        mass_fraction = unit_test_rp::X351;
        break;
    case 352:
        mass_fraction = unit_test_rp::X352;
        break;
    case 353:
        mass_fraction = unit_test_rp::X353;
        break;
    case 354:
        mass_fraction = unit_test_rp::X354;
        break;
    case 355:
        mass_fraction = unit_test_rp::X355;
        break;
    case 356:
        mass_fraction = unit_test_rp::X356;
        break;
    case 357:
        mass_fraction = unit_test_rp::X357;
        break;
    case 358:
        mass_fraction = unit_test_rp::X358;
        break;
    case 359:
        mass_fraction = unit_test_rp::X359;
        break;
    case 360:
        mass_fraction = unit_test_rp::X360;
        break;
    case 361:
        mass_fraction = unit_test_rp::X361;
        break;
    case 362:
        mass_fraction = unit_test_rp::X362;
        break;
    case 363:
        mass_fraction = unit_test_rp::X363;
        break;
    case 364:
        mass_fraction = unit_test_rp::X364;
        break;
    case 365:
        mass_fraction = unit_test_rp::X365;
        break;
    case 366:
        mass_fraction = unit_test_rp::X366;
        break;
    case 367:
        mass_fraction = unit_test_rp::X367;
        break;
    case 368:
        mass_fraction = unit_test_rp::X368;
        break;
    case 369:
        mass_fraction = unit_test_rp::X369;
        break;
    case 370:
        mass_fraction = unit_test_rp::X370;
        break;
    case 371:
        mass_fraction = unit_test_rp::X371;
        break;
    case 372:
        mass_fraction = unit_test_rp::X372;
        break;
    case 373:
        mass_fraction = unit_test_rp::X373;
        break;
    case 374:
        mass_fraction = unit_test_rp::X374;
        break;
    case 375:
        mass_fraction = unit_test_rp::X375;
        break;
    case 376:
        mass_fraction = unit_test_rp::X376;
        break;
    case 377:
        mass_fraction = unit_test_rp::X377;
        break;
    case 378:
        mass_fraction = unit_test_rp::X378;
        break;
    case 379:
        mass_fraction = unit_test_rp::X379;
        break;
    case 380:
        mass_fraction = unit_test_rp::X380;
        break;
    case 381:
        mass_fraction = unit_test_rp::X381;
        break;
    case 382:
        mass_fraction = unit_test_rp::X382;
        break;
    case 383:
        mass_fraction = unit_test_rp::X383;
        break;
    case 384:
        mass_fraction = unit_test_rp::X384;
        break;
    case 385:
        mass_fraction = unit_test_rp::X385;
        break;
    case 386:
        mass_fraction = unit_test_rp::X386;
        break;
    case 387:
        mass_fraction = unit_test_rp::X387;
        break;
    case 388:
        mass_fraction = unit_test_rp::X388;
        break;
    case 389:
        mass_fraction = unit_test_rp::X389;
        break;
    case 390:
        mass_fraction = unit_test_rp::X390;
        break;
    case 391:
        mass_fraction = unit_test_rp::X391;
        break;
    case 392:
        mass_fraction = unit_test_rp::X392;
        break;
    case 393:
        mass_fraction = unit_test_rp::X393;
        break;
    case 394:
        mass_fraction = unit_test_rp::X394;
        break;
    case 395:
        mass_fraction = unit_test_rp::X395;
        break;
    case 396:
        mass_fraction = unit_test_rp::X396;
        break;
    case 397:
        mass_fraction = unit_test_rp::X397;
        break;
    case 398:
        mass_fraction = unit_test_rp::X398;
        break;
    case 399:
        mass_fraction = unit_test_rp::X399;
        break;
    case 400:
        mass_fraction = unit_test_rp::X400;
        break;
    }

    return mass_fraction;

}

#endif
