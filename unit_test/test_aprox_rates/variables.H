#include <vector>
#include <string>

#ifndef _variables_H_
#define _variables_H_

#include <AMReX_Vector.H>

class plot_t {

public:

    int irho = -1;
    int itemp = -1;
    int ini56 = -1;
    int ic12ag = -1;
    int ic12ag_deboer17 = -1;
    int itriplealf = -1;
    int ic12c12 = -1;
    int ic12o16 = -1;
    int io16o16 = -1;
    int io16ag = -1;
    int ine20ag = -1;
    int img24ag = -1;
    int img24ap = -1;
    int ial27pg = -1;
    int ial27pg_old = -1;
    int isi28ag = -1;
    int isi28ap = -1;
    int ip31pg = -1;
    int is32ag = -1;
    int is32ap = -1;
    int icl35pg = -1;
    int iar36ag = -1;
    int iar36ap = -1;
    int ik39pg = -1;
    int ica40ag = -1;
    int ica40ap = -1;
    int isc43pg = -1;
    int iti44ag = -1;
    int iti44ap = -1;
    int iv47pg = -1;
    int icr48ag = -1;
    int icr48ap = -1;
    int imn51pg = -1;
    int ife52ag = -1;
    int ife52ap = -1;
    int ico55pg = -1;
    int ipp = -1;
    int ipng = -1;
    int idpg = -1;
    int ihe3ng = -1;
    int ihe3he3 = -1;
    int ihe3he4 = -1;
    int ic12pg = -1;
    int in14pg = -1;
    int in15pg = -1;
    int in15pa = -1;
    int io16pg = -1;
    int in14ag = -1;
    int ife52ng = -1;
    int ife53ng = -1;
    int ife54ng = -1;
    int ife54pg = -1;
    int ife54ap = -1;
    int ife55ng = -1;
    int ife56pg = -1;
    int ilanganke = -1;
    int iecapnuc = -1;

    int n_plot_comps = 0;

    amrex::Vector<std::string> names;

    int next_index(const int num) {
        int next = n_plot_comps;
        n_plot_comps += num;
        return next;
    }

};

plot_t init_variables();


#endif
