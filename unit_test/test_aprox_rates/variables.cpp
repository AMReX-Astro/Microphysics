#include <variables.H>
#include <network.H>

plot_t init_variables() {

    plot_t p;

    constexpr int n_tests = 52;

    p.irho = p.next_index(1);
    p.itemp = p.next_index(1);

    p.ic12ag = p.next_index(4);
    p.ic12ag_deboer17 = p.next_index(4);
    p.itriplealf = p.next_index(4);
    p.ic12c12 = p.next_index(4);

    p.names.resize(p.n_plot_comps);

    p.names[p.irho] = "density";
    p.names[p.itemp] = "temperature";

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
            p.names[2 + i*4 + j] = p.names[2 + i*4];
        }
        p.names[2 + i*4] += "_fr";
        p.names[2 + i*4 + 1] += "_dfrdt";
        p.names[2 + i*4 + 2] += "_rr";
        p.names[2 + i*4 + 3] += "_drrdt";
    }

  return p;
}

