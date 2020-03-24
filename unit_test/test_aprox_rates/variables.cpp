#include <variables.H>
#include <network.H>

plot_t init_variables() {

    plot_t p;

    constexpr int n_tests = 52;

    p.irho = p.next_index(1);
    p.itemp = p.next_index(1);
    p.ini56 = p.next_index(1);

    p.ic12ag = p.next_index(4);
    p.ic12ag_deboer17 = p.next_index(4);
    p.itriplealf = p.next_index(4);
    p.ic12c12 = p.next_index(4);
    p.ic12o16 = p.next_index(4);
    p.io16o16 = p.next_index(4);
    p.io16ag = p.next_index(4);
    p.ine20ag = p.next_index(4);
    p.img24ag = p.next_index(4);
    p.img24ap = p.next_index(4);
    p.ial27pg = p.next_index(4);
    p.ial27pg_old = p.next_index(4);
    p.isi28ag = p.next_index(4);
    p.isi28ap = p.next_index(4);
    p.ip31pg = p.next_index(4);
    p.is32ag = p.next_index(4);
    p.is32ap = p.next_index(4);
    p.icl35pg = p.next_index(4);
    p.iar36ag = p.next_index(4);
    p.iar36ap = p.next_index(4);
    p.ik39pg = p.next_index(4);
    p.ica40ag = p.next_index(4);
    p.ica40ap = p.next_index(4);
    p.isc43pg = p.next_index(4);
    p.iti44ag = p.next_index(4);
    p.iti44ap = p.next_index(4);
    p.iv47pg = p.next_index(4);
    p.icr48ag = p.next_index(4);
    p.icr48ap = p.next_index(4);
    p.imn51pg = p.next_index(4);
    p.ife52ag = p.next_index(4);
    p.ife52ap = p.next_index(4);
    p.ico55pg = p.next_index(4);
    p.ipp = p.next_index(4);
    p.ipng = p.next_index(4);
    p.idpg = p.next_index(4);
    p.ihe3ng = p.next_index(4);
    p.ihe3he3= p.next_index(4);
    p.ihe3he4 = p.next_index(4);
    p.ic12pg = p.next_index(4);
    p.in14pg = p.next_index(4);
    p.in15pg = p.next_index(4);
    p.in15pa = p.next_index(4);
    p.io16pg = p.next_index(4);
    p.in14ag = p.next_index(4);
    p.ife52ng = p.next_index(4);
    p.ife53ng = p.next_index(4);
    p.ife54ng = p.next_index(4);
    p.ife54pg = p.next_index(4);
    p.ife54ap = p.next_index(4);
    p.ife55ng = p.next_index(4);
    p.ife56pg = p.next_index(4);

    // langanke and ecapnuc are different so not included in n_tests 
    p.ilanganke = p.next_index(2);
    p.iecapnuc = p.next_index(4);


    p.names.resize(p.n_plot_comps);

    p.names[p.irho] = "density";
    p.names[p.itemp] = "temperature";
    p.names[p.ini56] = "X(ni56)";

    p.names[p.ic12ag] = "c12ag";
    p.names[p.ic12ag_deboer17] = "c12ag_deboer17";
    p.names[p.itriplealf] = "triplealf";
    p.names[p.ic12c12] = "c12c12";
    p.names[p.ic12o16] = "c12o16";
    p.names[p.io16o16] = "o16o16";
    p.names[p.io16ag] = "o16ag";

    p.names[p.ine20ag] = "ne20ag";
    p.names[p.img24ag] = "mg24ag";
    p.names[p.img24ap] = "mg24ap";
    p.names[p.ial27pg] = "al27pg";
    p.names[p.ial27pg_old] = "al27pg_old";
    p.names[p.isi28ag] = "si28ag";

    p.names[p.isi28ap] = "si28ap";
    p.names[p.ip31pg] = "p31pg";
    p.names[p.is32ag] = "s32ag";
    p.names[p.is32ap] = "s32ap";
    p.names[p.icl35pg] = "cl35pg";
    p.names[p.iar36ag] = "ar36ag";
    p.names[p.iar36ap] = "ar36ap";
    p.names[p.ik39pg] = "k39pg";
    p.names[p.ica40ag] = "ca40ag";
    p.names[p.ica40ap] = "ca40ap";
    p.names[p.isc43pg] = "sc43pg";
    p.names[p.iti44ag] = "ti44ag";
    p.names[p.iti44ap] = "ti44ap";
    p.names[p.iv47pg] = "v47pg";
    p.names[p.icr48ag] = "cr48ag";
    p.names[p.icr48ap] = "cr48ap";
    p.names[p.imn51pg] = "mn51pg";
    p.names[p.ife52ag] = "fe52ag";
    p.names[p.ife52ap] = "fe52ap";
    p.names[p.ico55pg] = "co55pg";
    p.names[p.ipp] = "pp";
    p.names[p.ipng] = "png";
    p.names[p.idpg] = "dpg";
    p.names[p.ihe3ng] = "he3ng";
    p.names[p.ihe3he3] = "he3he3";
    p.names[p.ihe3he4] = "he3he4";
    p.names[p.ic12pg] = "c12pg";
    p.names[p.in14pg] = "n14pg";
    p.names[p.in15pg] = "n15pg";
    p.names[p.in15pa] = "n15pa";
    p.names[p.io16pg] = "o16pg";
    p.names[p.in14ag] = "n14ag";
    p.names[p.ife52ng] = "fe52ng";
    p.names[p.ife53ng] = "fe53ng";
    p.names[p.ife54ng] = "fe54ng";
    p.names[p.ife54pg] = "fe54pg";
    p.names[p.ife54ap] = "fe54ap";
    p.names[p.ife55ng] = "fe55ng";
    p.names[p.ife56pg] = "fe56pg";

    for (auto i = 0; i < n_tests; ++i) {
        for (auto j = 1; j < 4; ++j) {
            p.names[3 + i*4 + j] = p.names[3 + i*4];
        }
        p.names[3 + i*4] += "_fr";
        p.names[3 + i*4 + 1] += "_dfrdt";
        p.names[3 + i*4 + 2] += "_rr";
        p.names[3 + i*4 + 3] += "_drrdt";
    }

    p.names[p.ilanganke] = "langanke_rn56ec";
    p.names[p.ilanganke+1] = "langanke_sn56ec";

    p.names[p.iecapnuc] = "ecapnuc_rpen";
    p.names[p.iecapnuc+1] = "ecapnuc_rnep";
    p.names[p.iecapnuc+2] = "ecapnuc_spenc";
    p.names[p.iecapnuc+3] = "ecapnuc_snepc";

  return p;
}

