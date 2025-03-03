#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>

#include <variables.H>
#include <network.H>
#include <eos.H>
#include <aprox_rates.H>
#include <microphysics_autodiff.H>

#include <cmath>

using namespace amrex;
using namespace unit_test_rp;

void aprox_rates_test(const Box& bx,
                      const amrex::Real dlogrho, const amrex::Real dlogT, const amrex::Real dNi,
                      const plot_t& vars,
                      Array4<amrex::Real> const sp) {

  const int ini56 = network_spec_index("nickel-56");

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
  {

    eos_extra_t eos_state;

    autodiff::dual temp_zone = std::pow(10.0_rt, std::log10(temp_min) + static_cast<amrex::Real>(j)*dlogT);
    autodiff::seed(temp_zone);
    eos_state.T = autodiff::val(temp_zone);

    amrex::Real dens_zone = std::pow(10.0, std::log10(dens_min) + static_cast<amrex::Real>(i)*dlogrho);
    eos_state.rho = dens_zone;

    amrex::Real ni56 = amrex::Real(k) * dNi;
    for(auto comp = 0; comp < NumSpec; ++comp) {
      eos_state.xn[comp] = 1.0_rt - ni56 / amrex::Real(NumSpec - 1);
    }

    eos_state.xn[ini56] = ni56;

    // call the EOS using rho, T
    eos(eos_input_rt, eos_state);

    auto tf = get_tfactors(temp_zone);

    // store state
    sp(i, j, k, vars.irho) = dens_zone;
    sp(i, j, k, vars.itemp) = autodiff::val(temp_zone);
    sp(i, j, k, vars.ini56) = ni56;

    autodiff::dual fr;
    autodiff::dual rr;

    rate_c12ag(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ic12ag) = autodiff::val(fr);
    sp(i, j, k, vars.ic12ag+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ic12ag+2) = autodiff::val(rr);
    sp(i, j, k, vars.ic12ag+3) = autodiff::derivative(rr);

    rate_triplealf(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.itriplealf) = autodiff::val(fr);
    sp(i, j, k, vars.itriplealf+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.itriplealf+2) = autodiff::val(rr);
    sp(i, j, k, vars.itriplealf+3) = autodiff::derivative(rr);

    rate_c12c12(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ic12c12) = autodiff::val(fr);
    sp(i, j, k, vars.ic12c12+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ic12c12+2) = autodiff::val(rr);
    sp(i, j, k, vars.ic12c12+3) = autodiff::derivative(rr);

    rate_c12o16(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ic12o16) = autodiff::val(fr);
    sp(i, j, k, vars.ic12o16+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ic12o16+2) = autodiff::val(rr);
    sp(i, j, k, vars.ic12o16+3) = autodiff::derivative(rr);

    rate_o16o16(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.io16o16) = autodiff::val(fr);
    sp(i, j, k, vars.io16o16+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.io16o16+2) = autodiff::val(rr);
    sp(i, j, k, vars.io16o16+3) = autodiff::derivative(rr);

    rate_o16ag(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.io16ag) = autodiff::val(fr);
    sp(i, j, k, vars.io16ag+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.io16ag+2) = autodiff::val(rr);
    sp(i, j, k, vars.io16ag+3) = autodiff::derivative(rr);

    rate_ne20ag(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ine20ag) = autodiff::val(fr);
    sp(i, j, k, vars.ine20ag+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ine20ag+2) = autodiff::val(rr);
    sp(i, j, k, vars.ine20ag+3) = autodiff::derivative(rr);

    rate_mg24ag(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.img24ag) = autodiff::val(fr);
    sp(i, j, k, vars.img24ag+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.img24ag+2) = autodiff::val(rr);
    sp(i, j, k, vars.img24ag+3) = autodiff::derivative(rr);

    rate_mg24ap(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.img24ap) = autodiff::val(fr);
    sp(i, j, k, vars.img24ap+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.img24ap+2) = autodiff::val(rr);
    sp(i, j, k, vars.img24ap+3) = autodiff::derivative(rr);

    rate_al27pg(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ial27pg) = autodiff::val(fr);
    sp(i, j, k, vars.ial27pg+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ial27pg+2) = autodiff::val(rr);
    sp(i, j, k, vars.ial27pg+3) = autodiff::derivative(rr);

    rate_al27pg_old(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ial27pg_old) = autodiff::val(fr);
    sp(i, j, k, vars.ial27pg_old+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ial27pg_old+2) = autodiff::val(rr);
    sp(i, j, k, vars.ial27pg_old+3) = autodiff::derivative(rr);

    rate_si28ag(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.isi28ag) = autodiff::val(fr);
    sp(i, j, k, vars.isi28ag+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.isi28ag+2) = autodiff::val(rr);
    sp(i, j, k, vars.isi28ag+3) = autodiff::derivative(rr);

    rate_si28ap(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.isi28ap) = autodiff::val(fr);
    sp(i, j, k, vars.isi28ap+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.isi28ap+2) = autodiff::val(rr);
    sp(i, j, k, vars.isi28ap+3) = autodiff::derivative(rr);

    rate_p31pg(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ip31pg) = autodiff::val(fr);
    sp(i, j, k, vars.ip31pg+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ip31pg+2) = autodiff::val(rr);
    sp(i, j, k, vars.ip31pg+3) = autodiff::derivative(rr);

    rate_s32ag(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.is32ag) = autodiff::val(fr);
    sp(i, j, k, vars.is32ag+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.is32ag+2) = autodiff::val(rr);
    sp(i, j, k, vars.is32ag+3) = autodiff::derivative(rr);

    rate_s32ap(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.is32ap) = autodiff::val(fr);
    sp(i, j, k, vars.is32ap+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.is32ap+2) = autodiff::val(rr);
    sp(i, j, k, vars.is32ap+3) = autodiff::derivative(rr);

    rate_cl35pg(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.icl35pg) = autodiff::val(fr);
    sp(i, j, k, vars.icl35pg+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.icl35pg+2) = autodiff::val(rr);
    sp(i, j, k, vars.icl35pg+3) = autodiff::derivative(rr);

    rate_ar36ag(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.iar36ag) = autodiff::val(fr);
    sp(i, j, k, vars.iar36ag+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.iar36ag+2) = autodiff::val(rr);
    sp(i, j, k, vars.iar36ag+3) = autodiff::derivative(rr);

    rate_ar36ap(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.iar36ap) = autodiff::val(fr);
    sp(i, j, k, vars.iar36ap+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.iar36ap+2) = autodiff::val(rr);
    sp(i, j, k, vars.iar36ap+3) = autodiff::derivative(rr);

    rate_k39pg(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ik39pg) = autodiff::val(fr);
    sp(i, j, k, vars.ik39pg+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ik39pg+2) = autodiff::val(rr);
    sp(i, j, k, vars.ik39pg+3) = autodiff::derivative(rr);

    rate_ca40ag(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ica40ag) = autodiff::val(fr);
    sp(i, j, k, vars.ica40ag+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ica40ag+2) = autodiff::val(rr);
    sp(i, j, k, vars.ica40ag+3) = autodiff::derivative(rr);

    rate_ca40ap(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ica40ap) = autodiff::val(fr);
    sp(i, j, k, vars.ica40ap+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ica40ap+2) = autodiff::val(rr);
    sp(i, j, k, vars.ica40ap+3) = autodiff::derivative(rr);

    rate_sc43pg(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.isc43pg) = autodiff::val(fr);
    sp(i, j, k, vars.isc43pg+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.isc43pg+2) = autodiff::val(rr);
    sp(i, j, k, vars.isc43pg+3) = autodiff::derivative(rr);

    rate_ti44ag(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.iti44ag) = autodiff::val(fr);
    sp(i, j, k, vars.iti44ag+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.iti44ag+2) = autodiff::val(rr);
    sp(i, j, k, vars.iti44ag+3) = autodiff::derivative(rr);

    rate_ti44ap(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.iti44ap) = autodiff::val(fr);
    sp(i, j, k, vars.iti44ap+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.iti44ap+2) = autodiff::val(rr);
    sp(i, j, k, vars.iti44ap+3) = autodiff::derivative(rr);

    rate_v47pg(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.iv47pg) = autodiff::val(fr);
    sp(i, j, k, vars.iv47pg+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.iv47pg+2) = autodiff::val(rr);
    sp(i, j, k, vars.iv47pg+3) = autodiff::derivative(rr);

    rate_cr48ag(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.icr48ag) = autodiff::val(fr);
    sp(i, j, k, vars.icr48ag+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.icr48ag+2) = autodiff::val(rr);
    sp(i, j, k, vars.icr48ag+3) = autodiff::derivative(rr);

    rate_cr48ap(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.icr48ap) = autodiff::val(fr);
    sp(i, j, k, vars.icr48ap+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.icr48ap+2) = autodiff::val(rr);
    sp(i, j, k, vars.icr48ap+3) = autodiff::derivative(rr);

    rate_mn51pg(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.imn51pg) = autodiff::val(fr);
    sp(i, j, k, vars.imn51pg+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.imn51pg+2) = autodiff::val(rr);
    sp(i, j, k, vars.imn51pg+3) = autodiff::derivative(rr);

    rate_fe52ag(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ife52ag) = autodiff::val(fr);
    sp(i, j, k, vars.ife52ag+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ife52ag+2) = autodiff::val(rr);
    sp(i, j, k, vars.ife52ag+3) = autodiff::derivative(rr);

    rate_fe52ap(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ife52ap) = autodiff::val(fr);
    sp(i, j, k, vars.ife52ap+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ife52ap+2) = autodiff::val(rr);
    sp(i, j, k, vars.ife52ap+3) = autodiff::derivative(rr);

    rate_co55pg(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ico55pg) = autodiff::val(fr);
    sp(i, j, k, vars.ico55pg+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ico55pg+2) = autodiff::val(rr);
    sp(i, j, k, vars.ico55pg+3) = autodiff::derivative(rr);

    rate_pp(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ipp) = autodiff::val(fr);
    sp(i, j, k, vars.ipp+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ipp+2) = autodiff::val(rr);
    sp(i, j, k, vars.ipp+3) = autodiff::derivative(rr);

    rate_png(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ipng) = autodiff::val(fr);
    sp(i, j, k, vars.ipng+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ipng+2) = autodiff::val(rr);
    sp(i, j, k, vars.ipng+3) = autodiff::derivative(rr);

    rate_dpg(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.idpg) = autodiff::val(fr);
    sp(i, j, k, vars.idpg+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.idpg+2) = autodiff::val(rr);
    sp(i, j, k, vars.idpg+3) = autodiff::derivative(rr);

    rate_he3ng(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ihe3ng) = autodiff::val(fr);
    sp(i, j, k, vars.ihe3ng+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ihe3ng+2) = autodiff::val(rr);
    sp(i, j, k, vars.ihe3ng+3) = autodiff::derivative(rr);

    rate_he3he3(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ihe3he3) = autodiff::val(fr);
    sp(i, j, k, vars.ihe3he3+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ihe3he3+2) = autodiff::val(rr);
    sp(i, j, k, vars.ihe3he3+3) = autodiff::derivative(rr);

    rate_he3he4(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ihe3he4) = autodiff::val(fr);
    sp(i, j, k, vars.ihe3he4+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ihe3he4+2) = autodiff::val(rr);
    sp(i, j, k, vars.ihe3he4+3) = autodiff::derivative(rr);

    rate_c12pg(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ic12pg) = autodiff::val(fr);
    sp(i, j, k, vars.ic12pg+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ic12pg+2) = autodiff::val(rr);
    sp(i, j, k, vars.ic12pg+3) = autodiff::derivative(rr);

    rate_n14pg(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.in14pg) = autodiff::val(fr);
    sp(i, j, k, vars.in14pg+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.in14pg+2) = autodiff::val(rr);
    sp(i, j, k, vars.in14pg+3) = autodiff::derivative(rr);

    rate_n15pg(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.in15pg) = autodiff::val(fr);
    sp(i, j, k, vars.in15pg+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.in15pg+2) = autodiff::val(rr);
    sp(i, j, k, vars.in15pg+3) = autodiff::derivative(rr);

    rate_n15pa(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.in15pa) = autodiff::val(fr);
    sp(i, j, k, vars.in15pa+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.in15pa+2) = autodiff::val(rr);
    sp(i, j, k, vars.in15pa+3) = autodiff::derivative(rr);

    rate_o16pg(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.io16pg) = autodiff::val(fr);
    sp(i, j, k, vars.io16pg+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.io16pg+2) = autodiff::val(rr);
    sp(i, j, k, vars.io16pg+3) = autodiff::derivative(rr);

    rate_n14ag(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.in14ag) = autodiff::val(fr);
    sp(i, j, k, vars.in14ag+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.in14ag+2) = autodiff::val(rr);
    sp(i, j, k, vars.in14ag+3) = autodiff::derivative(rr);

    rate_fe52ng(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ife52ng) = autodiff::val(fr);
    sp(i, j, k, vars.ife52ng+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ife52ng+2) = autodiff::val(rr);
    sp(i, j, k, vars.ife52ng+3) = autodiff::derivative(rr);

    rate_fe53ng(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ife53ng) = autodiff::val(fr);
    sp(i, j, k, vars.ife53ng+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ife53ng+2) = autodiff::val(rr);
    sp(i, j, k, vars.ife53ng+3) = autodiff::derivative(rr);

    rate_fe54ng(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ife54ng) = autodiff::val(fr);
    sp(i, j, k, vars.ife54ng+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ife54ng+2) = autodiff::val(rr);
    sp(i, j, k, vars.ife54ng+3) = autodiff::derivative(rr);

    rate_fe54pg(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ife54pg) = autodiff::val(fr);
    sp(i, j, k, vars.ife54pg+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ife54pg+2) = autodiff::val(rr);
    sp(i, j, k, vars.ife54pg+3) = autodiff::derivative(rr);

    rate_fe54ap(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ife54ap) = autodiff::val(fr);
    sp(i, j, k, vars.ife54ap+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ife54ap+2) = autodiff::val(rr);
    sp(i, j, k, vars.ife54ap+3) = autodiff::derivative(rr);

    rate_fe55ng(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ife55ng) = autodiff::val(fr);
    sp(i, j, k, vars.ife55ng+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ife55ng+2) = autodiff::val(rr);
    sp(i, j, k, vars.ife55ng+3) = autodiff::derivative(rr);

    rate_fe56pg(tf, dens_zone, fr, rr);

    sp(i, j, k, vars.ife56pg) = autodiff::val(fr);
    sp(i, j, k, vars.ife56pg+1) = autodiff::derivative(fr);
    sp(i, j, k, vars.ife56pg+2) = autodiff::val(rr);
    sp(i, j, k, vars.ife56pg+3) = autodiff::derivative(rr);


    amrex::Real rn56ec;
    amrex::Real sn56ec;
    langanke(autodiff::val(temp_zone), dens_zone, eos_state.xn[ini56],
             eos_state.y_e, rn56ec, sn56ec);

    sp(i, j, k, vars.ilanganke) = rn56ec;
    sp(i, j, k, vars.ilanganke+1) = sn56ec;

    amrex::Real rpen;
    amrex::Real rnep;
    amrex::Real spenc;
    amrex::Real snepc;
    ecapnuc(eos_state.eta, dens_zone, rpen, rnep, spenc, snepc);

    sp(i, j, k, vars.iecapnuc) = rpen;
    sp(i, j, k, vars.iecapnuc+1) = rnep;
    sp(i, j, k, vars.iecapnuc+2) = spenc;
    sp(i, j, k, vars.iecapnuc+3) = snepc;
  });

}


void aprox_rates_extra_c12ag(const Box& bx,
                             const amrex::Real dlogrho, const amrex::Real dlogT, const amrex::Real dNi,
                             const plot_t& vars,
                             Array4<amrex::Real> const sp) {

    const int ini56 = network_spec_index("nickel-56");

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

        eos_extra_t eos_state;

        autodiff::dual temp_zone = std::pow(10.0_rt, std::log10(temp_min) + static_cast<amrex::Real>(j)*dlogT);
        autodiff::seed(temp_zone);
        eos_state.T = autodiff::val(temp_zone);

        amrex::Real dens_zone = std::pow(10.0, std::log10(dens_min) + static_cast<amrex::Real>(i)*dlogrho);
        eos_state.rho = dens_zone;

        amrex::Real ni56 = amrex::Real(k) * dNi;
        for(auto comp = 0; comp < NumSpec; ++comp) {
            eos_state.xn[comp] = 1.0_rt - ni56 / amrex::Real(NumSpec - 1);
        }

        eos_state.xn[ini56] = ni56;

        // call the EOS using rho, T
        eos(eos_input_rt, eos_state);

        auto tf = get_tfactors(temp_zone);

        autodiff::dual fr;
        autodiff::dual rr;

        // override to get the other rate
        use_c12ag_deboer17 = 1;

        rate_c12ag(tf, dens_zone, fr, rr);

        sp(i, j, k, vars.ic12ag_deboer17) = autodiff::val(fr);
        sp(i, j, k, vars.ic12ag_deboer17+1) = autodiff::derivative(fr);
        sp(i, j, k, vars.ic12ag_deboer17+2) = autodiff::val(rr);
        sp(i, j, k, vars.ic12ag_deboer17+3) = autodiff::derivative(rr);

    });
}
