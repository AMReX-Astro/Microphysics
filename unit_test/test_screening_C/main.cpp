#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>


using namespace amrex;

#include "test_screen.H"
#include "test_screen_F.H"
#include "AMReX_buildInfo.H"

#include <network.H>
#include <eos.H>
#include <screen.H>
#include <screening.H>

#include <variables.H>

#include <cmath>

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    main_main();

    amrex::Finalize();
    return 0;
}

void main_main ()
{

    // AMREX_SPACEDIM: number of dimensions
    int n_cell, max_grid_size;
    Vector<int> bc_lo(AMREX_SPACEDIM,0);
    Vector<int> bc_hi(AMREX_SPACEDIM,0);

    // inputs parameters
    {
        // ParmParse is way of reading inputs from the inputs file
        ParmParse pp;

        // We need to get n_cell from the inputs file - this is the
        // number of cells on each side of a square (or cubic) domain.
        pp.get("n_cell", n_cell);

        // The domain is broken into boxes of size max_grid_size
        max_grid_size = 32;
        pp.query("max_grid_size", max_grid_size);

    }

    Vector<int> is_periodic(AMREX_SPACEDIM,0);
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
      is_periodic[idim] = 1;
    }

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
        IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);

        // Break up boxarray "ba" into chunks no larger than
        // "max_grid_size" along a direction
        ba.maxSize(max_grid_size);

        // This defines the physical box, [0, 1] in each direction.
        RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                         {AMREX_D_DECL(1.0, 1.0, 1.0)});

        // This defines a Geometry object
        geom.define(domain, &real_box,
                    CoordSys::cartesian, is_periodic.data());
    }

    // Nghost = number of ghost cells for each array
    int Nghost = 0;

    // do the runtime parameter initializations and microphysics inits
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "reading extern runtime parameters ..." << std::endl;
    }

    ParmParse ppa("amr");

    std::string probin_file = "probin";

    ppa.query("probin_file", probin_file);

    const int probin_file_length = probin_file.length();
    Vector<int> probin_file_name(probin_file_length);

    for (int i = 0; i < probin_file_length; i++)
      probin_file_name[i] = probin_file[i];

    init_unit_test(probin_file_name.dataPtr(), &probin_file_length);

    init_extern_parameters();

    eos_init();

    screening_init();

    auto vars = init_variables();

    // time = starting time in the simulation
    Real time = 0.0;

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // we allocate our main multifabs
    MultiFab state(ba, dm, vars.n_plot_comps, Nghost);

    // Initialize the state to zero; we will fill
    // it in below in do_eos.
    state.setVal(0.0);

    // What time is it now?  We'll use this to compute total run time.
    Real strt_time = ParallelDescriptor::second();

    const int ih1 = network_spec_index("hydrogen-1");
    if (ih1 < 0) amrex::Error("Error: ih1 not found");

    const int ihe3 = network_spec_index("helium-3");
    if (ihe3 < 0) amrex::Error("Error: ihe3 not found");

    const int ihe4 = network_spec_index("helium-4");
    if (ihe4 < 0) amrex::Error("Error: ihe4 not found");

    const int ic12 = network_spec_index("carbon-12");
    if (ic12 < 0) amrex::Error("Error: ic12 not found");

    const int in14 = network_spec_index("nitrogen-14");
    if (in14 < 0) amrex::Error("Error: in14 not found");

    const int io16 = network_spec_index("oxygen-16");
    if (io16 < 0) amrex::Error("Error: io16 not found");

    const int ine20 = network_spec_index("neon-20");
    if (ine20 < 0) amrex::Error("Error: ine20 not found");

    const int img24 = network_spec_index("magnesium-24");
    if (img24 < 0) amrex::Error("Error: img24 not found");

    const int isi28 = network_spec_index("silicon-28");
    if (isi28 < 0) amrex::Error("Error: isi28 not found");

    const int is32 = network_spec_index("sulfur-32");
    if (is32 < 0) amrex::Error("Error: is32 not found");

    const int iar36 = network_spec_index("argon-36");    
    if (iar36 < 0) amrex::Error("Error: iar36 not found");

    const int ica40 = network_spec_index("calcium-40");
    if (ica40 < 0) amrex::Error("Error: ica40 not found");

    const int iti44 = network_spec_index("titanium-44");
    if (iti44 < 0) amrex::Error("Error: iti44 not found");

    const int icr48 = network_spec_index("chromium-48");
    if (icr48 < 0) amrex::Error("Error: icr48 not found");

    const int ife52 = network_spec_index("iron-52");
    if (ife52 < 0) amrex::Error("Error: ife52 not found");

    const int ife54 = network_spec_index("iron-54");
    if (ife54 < 0) amrex::Error("Error: ife54 not found");

    const int ife56 = network_spec_index("iron-56");
    if (ife56 < 0) amrex::Error("Error: ife56 not found");


    // initialize the screening factors
    int jscr_init = 0;
    add_screening_factor(jscr_init, zion[ihe4],aion[ihe4],zion[ihe4],aion[ihe4]);

    jscr_init++;
    add_screening_factor(jscr_init, zion[ihe4],aion[ihe4],4.0e0_rt,8.0e0_rt);

    jscr_init++;
    add_screening_factor(jscr_init, zion[ic12],aion[ic12],zion[ihe4],aion[ihe4]);

    jscr_init++;
    add_screening_factor(jscr_init, zion[ic12],aion[ic12],zion[ic12],aion[ic12]);

    jscr_init++;
    add_screening_factor(jscr_init, zion[ic12],aion[ic12],zion[io16],aion[io16]);

    jscr_init++;
    add_screening_factor(jscr_init, zion[io16],aion[io16],zion[io16],aion[io16]);

    jscr_init++;
    add_screening_factor(jscr_init, zion[io16],aion[io16],zion[ihe4],aion[ihe4]);

    jscr_init++;
    add_screening_factor(jscr_init, zion[ine20],aion[ine20],zion[ihe4],aion[ihe4]);

    jscr_init++;
    add_screening_factor(jscr_init, zion[img24],aion[img24],zion[ihe4],aion[ihe4]);

    jscr_init++;
    add_screening_factor(jscr_init, 13.0e0_rt,27.0e0_rt,1.0e0_rt,1.0e0_rt);

    jscr_init++;
    add_screening_factor(jscr_init, zion[isi28],aion[isi28],zion[ihe4],aion[ihe4]);

    jscr_init++;
    add_screening_factor(jscr_init, 15.0e0_rt,31.0e0_rt,1.0e0_rt,1.0e0_rt);

    jscr_init++;
    add_screening_factor(jscr_init, zion[is32],aion[is32],zion[ihe4],aion[ihe4]);

    jscr_init++;
    add_screening_factor(jscr_init, 17.0e0_rt,35.0e0_rt,1.0e0_rt,1.0e0_rt);

    jscr_init++;
    add_screening_factor(jscr_init, zion[iar36],aion[iar36],zion[ihe4],aion[ihe4]);

    jscr_init++;
    add_screening_factor(jscr_init, 19.0e0_rt,39.0e0_rt,1.0e0_rt,1.0e0_rt);

    jscr_init++;
    add_screening_factor(jscr_init, zion[ica40],aion[ica40],zion[ihe4],aion[ihe4]);

    jscr_init++;
    add_screening_factor(jscr_init, 21.0e0_rt,43.0e0_rt,1.0e0_rt,1.0e0_rt);

    jscr_init++;
    add_screening_factor(jscr_init, zion[iti44],aion[iti44],zion[ihe4],aion[ihe4]);

    jscr_init++;
    add_screening_factor(jscr_init, 23.0e0_rt,47.0e0_rt,1.0e0_rt,1.0e0_rt);

    jscr_init++;
    add_screening_factor(jscr_init, zion[icr48],aion[icr48],zion[ihe4],aion[ihe4]);

    jscr_init++;
    add_screening_factor(jscr_init, 25.0e0_rt,51.0e0_rt,1.0e0_rt,1.0e0_rt);

    jscr_init++;
    add_screening_factor(jscr_init, zion[ife52],aion[ife52],zion[ihe4],aion[ihe4]);

    jscr_init++;
    add_screening_factor(jscr_init, 27.0e0_rt,55.0e0_rt,1.0e0_rt,1.0e0_rt);

    jscr_init++;
    add_screening_factor(jscr_init, zion[ife54],aion[ife54],1.0e0_rt,1.0e0_rt);

    jscr_init++;
    add_screening_factor(jscr_init, zion[ife54],aion[ife54],zion[ihe4],aion[ihe4]);

    jscr_init++;
    add_screening_factor(jscr_init, zion[ife56],aion[ife56],1.0e0_rt,1.0e0_rt);

    jscr_init++;
    add_screening_factor(jscr_init, 1.0e0_rt,2.0e0_rt,zion[ih1],aion[ih1]);

    jscr_init++;
    add_screening_factor(jscr_init, zion[ih1],aion[ih1],zion[ih1],aion[ih1]);

    jscr_init++;
    add_screening_factor(jscr_init, zion[ihe3],aion[ihe3],zion[ihe3],aion[ihe3]);

    jscr_init++;
    add_screening_factor(jscr_init, zion[ihe3],aion[ihe3],zion[ihe4],aion[ihe4]);

    jscr_init++;
    add_screening_factor(jscr_init, zion[ic12],aion[ic12],zion[ih1],aion[ih1]);

    jscr_init++;
    add_screening_factor(jscr_init, zion[in14],aion[in14],zion[ih1],aion[ih1]);

    jscr_init++;
    add_screening_factor(jscr_init, zion[io16],aion[io16],zion[ih1],aion[ih1]);

    jscr_init++;
    add_screening_factor(jscr_init, zion[in14],aion[in14],zion[ihe4],aion[ihe4]);


    Real dlogrho = 0.0e0_rt;
    Real dlogT   = 0.0e0_rt;
    Real dmetal  = 0.0e0_rt;

    if (n_cell > 1) {
        dlogrho = (std::log10(dens_max) - std::log10(dens_min))/(n_cell - 1);
        dlogT   = (std::log10(temp_max) - std::log10(temp_min))/(n_cell - 1);
        dmetal  = (metalicity_max  - 0.0)/(n_cell - 1);
    }

    // Initialize the state and compute the different thermodynamics
    // by inverting the EOS
    for ( MFIter mfi(state); mfi.isValid(); ++mfi )
    {
      const Box& bx = mfi.validbox();

      Array4<Real> const sp = state.array(mfi);

      AMREX_PARALLEL_FOR_3D(bx, i, j, k,
      {

        // set the composition -- approximately solar
        Real metalicity = 0.0 + static_cast<Real> (k) * dmetal;

        Real xn[NumSpec];
        Real ymass[NumSpec];

        for (int n = 0; n < NumSpec; n++) {
          xn[n] = metalicity/(NumSpec - 2);
        }
        xn[ih1] = 0.75 - 0.5*metalicity;
        xn[ihe4] = 0.25 - 0.5*metalicity;

        for (int n = 0; n < NumSpec; n++) {
          ymass[n] = xn[n] / aion[n];
        }

        Real temp_zone = std::pow(10.0, std::log10(temp_min) + static_cast<Real>(j)*dlogT);

        Real dens_zone = std::pow(10.0, std::log10(dens_min) + static_cast<Real>(i)*dlogrho);

        // store default state
        sp(i, j, k, vars.irho) = dens_zone;
        sp(i, j, k, vars.itemp) = temp_zone;
        for (int n = 0; n < NumSpec; n++) {
          sp(i, j, k, vars.ispec+n) = xn[n];
        }

        plasma_state_t pstate;
        fill_plasma_state(pstate, temp_zone, dens_zone, ymass);

        Real sc1a;
        Real sc1adt;
        Real sc1add;

        // 3-alpha
        int jscr = 0;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_he4_he4) = sc1a;
        sp(i, j, k, vars.iscn_he4_he4_dt) = sc1adt;

        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_he4_be8) = sc1a;
        sp(i, j, k, vars.iscn_he4_be8_dt) = sc1adt;

        // c12(a,g)o16
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_c12_he4) = sc1a;
        sp(i, j, k, vars.iscn_c12_he4_dt) = sc1adt;

        // c12 + c12
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_c12_c12) = sc1a;
        sp(i, j, k, vars.iscn_c12_c12_dt) = sc1adt;

        // c12 + o16
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_c12_o16) = sc1a;
        sp(i, j, k, vars.iscn_c12_o16_dt) = sc1adt;

        // o16 + o16
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_o16_o16) = sc1a;
        sp(i, j, k, vars.iscn_o16_o16_dt) = sc1adt;

        // o16 + ne20
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_o16_he4) = sc1a;
        sp(i, j, k, vars.iscn_o16_he4_dt) = sc1adt;

        // ne20(a,g)mg24
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_ne20_he4) = sc1a;
        sp(i, j, k, vars.iscn_ne20_he4_dt) = sc1adt;

        // mg24(a,g)si28
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_mg24_he4) = sc1a;
        sp(i, j, k, vars.iscn_mg24_he4_dt) = sc1adt;

        // al27(p,g)si28
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_al27_p) = sc1a;
        sp(i, j, k, vars.iscn_al27_p_dt) = sc1adt;

        // si28 tp s32
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_si28_he4) = sc1a;
        sp(i, j, k, vars.iscn_si28_he4_dt) = sc1adt;

        // p31(p,g)s32
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_p31_p) = sc1a;
        sp(i, j, k, vars.iscn_p31_p_dt) = sc1adt;

        // s32 to ar36
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_s32_he4) = sc1a;
        sp(i, j, k, vars.iscn_s32_he4_dt) = sc1adt;

        // cl35(p,g)ar36
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_cl35_p) = sc1a;
        sp(i, j, k, vars.iscn_cl35_p_dt) = sc1adt;

        // ar36 to ca40
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_ar36_he4) = sc1a;
        sp(i, j, k, vars.iscn_ar36_he4_dt) = sc1adt;

        // k39(p,g)ca40
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_k39_p) = sc1a;
        sp(i, j, k, vars.iscn_k39_p_dt) = sc1adt;

        // ca40 to ti44
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_ca40_he4) = sc1a;
        sp(i, j, k, vars.iscn_ca40_he4_dt) = sc1adt;

        // sc43(p,g)ti44
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_sc43_p) = sc1a;
        sp(i, j, k, vars.iscn_sc43_p_dt) = sc1adt;

        // ti44 to cr48
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_ti44_he4) = sc1a;
        sp(i, j, k, vars.iscn_ti44_he4_dt) = sc1adt;

        // v47(p,g)cr48
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_v47_p) = sc1a;
        sp(i, j, k, vars.iscn_v47_p_dt) = sc1adt;

        // cr48 to fe52
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_cr48_he4) = sc1a;
        sp(i, j, k, vars.iscn_cr48_he4_dt) = sc1adt;

        // mn51(p,g)fe52
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_mn51_p) = sc1a;
        sp(i, j, k, vars.iscn_mn51_p_dt) = sc1adt;

        // fe to ni
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_fe52_he4) = sc1a;
        sp(i, j, k, vars.iscn_fe52_he4_dt) = sc1adt;

        // co55(p,g)ni56
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_co55_p) = sc1a;
        sp(i, j, k, vars.iscn_co55_p_dt) = sc1adt;

        // fe54(p,g)co55
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_fe54_p) = sc1a;
        sp(i, j, k, vars.iscn_fe54_p_dt) = sc1adt;

        // fe54(a,p)co57
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_fe54_he4) = sc1a;
        sp(i, j, k, vars.iscn_fe54_he4_dt) = sc1adt;

        // fe56(p,g)co57
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_fe56_p) = sc1a;
        sp(i, j, k, vars.iscn_fe56_p_dt) = sc1adt;

        // d + p
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_d_p) = sc1a;
        sp(i, j, k, vars.iscn_d_p_dt) = sc1adt;

        // pp
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_p_p) = sc1a;
        sp(i, j, k, vars.iscn_p_p_dt) = sc1adt;

        // he3 + he3
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_he3_he3) = sc1a;
        sp(i, j, k, vars.iscn_he3_he3_dt) = sc1adt;

        // he3 + he4
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_he3_he4) = sc1a;
        sp(i, j, k, vars.iscn_he3_he4_dt) = sc1adt;

        // c12(p,g)n13
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_c12_p) = sc1a;
        sp(i, j, k, vars.iscn_c12_p_dt) = sc1adt;

        // n14(p,g)o15
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_n14_p) = sc1a;
        sp(i, j, k, vars.iscn_n14_p_dt) = sc1adt;

        // o16(p,g)f17
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_o16_p) = sc1a;
        sp(i, j, k, vars.iscn_o16_p_dt) = sc1adt;

        // n14(a,g)f18
        jscr++;
        screen5(pstate, jscr, sc1a, sc1adt, sc1add);
        sp(i, j, k, vars.iscn_n14_he4) = sc1a;
        sp(i, j, k, vars.iscn_n14_he4_dt) = sc1adt;


      });

    }

    // Call the timer again and compute the maximum difference between
    // the start time and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time, IOProc);


    std::string name = "test_screening_C";

    // Write a plotfile
    WriteSingleLevelPlotfile(name, state, vars.names, geom, time, 0);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;

}
