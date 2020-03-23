#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>


using namespace amrex;

#include "test_aprox_rates.H"
#include "test_aprox_rates_F.H"
#include "AMReX_buildInfo.H"

#include <network.H>
#include <eos.H>
#include <variables.H>
#include <aprox_rates.H>

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

    rates_init();

    auto vars = init_variables();

    // time = starting time in the simulation
    Real time = 0.0;

    // How Boxes are distributed among MPI processes
    DistributionMapping dm(ba);

    // we allocate our main multifabs
    MultiFab state(ba, dm, vars.n_plot_comps, Nghost);

    // Initialize the state to zero; we will fill
    // it in below in do_eos.
    state.setVal(0.0);

    // What time is it now?  We'll use this to compute total run time.
    Real strt_time = ParallelDescriptor::second();

    const int ih1 = network_spec_index("hydrogen-1");
    const int ihe4 = network_spec_index("helium-4");

    Real dlogrho = 0.0e0_rt;
    Real dlogT = 0.0e0_rt;

    if (n_cell > 1) {
        dlogrho = (log10(dens_max) - log10(dens_min))/(n_cell - 1);
        dlogT   = (log10(temp_max) - log10(temp_min))/(n_cell - 1);
    }

    // Initialize the state and compute the different thermodynamics
    // by inverting the EOS
    for ( MFIter mfi(state); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        Array4<Real> const sp = state.array(mfi);

        AMREX_PARALLEL_FOR_3D(bx, i, j, k,
        {
            Real temp_zone = std::pow(10.0_rt, log10(temp_min) + static_cast<Real>(j)*dlogT);
            // eos_state.T = temp_zone;

            Real dens_zone = std::pow(10.0, log10(dens_min) + static_cast<Real>(i)*dlogrho);
            // eos_state.rho = dens_zone;

            auto tf = get_tfactors(temp_zone);

            // store state
            sp(i, j, k, vars.irho) = dens_zone;
            sp(i, j, k, vars.itemp) = temp_zone;

            Real fr;
            Real dfrdt;
            Real rr;
            Real drrdt;

            rate_c12ag(tf, dens_zone, fr, dfrdt, rr, drrdt);

            sp(i, j, k, vars.ifr) = fr;
            sp(i, j, k, vars.idfrdt) = dfrdt;
            sp(i, j, k, vars.irr) = rr;
            sp(i, j, k, vars.idrrdt) = drrdt;

            rate_c12ag_deboer17(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_triplealf(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_c12c12(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_c12o16(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_o16o16(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_o16ag(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_ne20ag(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_mg24ag(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_mg24ap(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_al27pg(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_al27pg_old(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_si28ag(tf, dens_zone, fr, dfrdt, rr, drrdt);
            
            rate_si28ap(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_p31pg(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_s32ag(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_s32ap(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_cl35pg(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_ar36ag(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_ar36ap(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_k39pg(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_ca40ag(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_ca40ap(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_sc43pg(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_ti44ag(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_ti44ap(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_v47pg(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_cr48ag(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_cr48ap(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_mn51pg(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_fe52ag(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_fe52ap(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_co55pg(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_pp(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_png(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_dpg(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_he3ng(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_he3he3(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_he3he4(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_c12pg(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_n14pg(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_n15pg(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_n15pa(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_o16pg(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_n14ag(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_fe52ng(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_fe53ng(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_fe54ng(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_fe54pg(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_fe54ap(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_fe55ng(tf, dens_zone, fr, dfrdt, rr, drrdt);

            rate_fe56pg(tf, dens_zone, fr, dfrdt, rr, drrdt);

            // TODO: langanke and ecapnuc have different inputs/outputs
        });
    }

    // Call the timer again and compute the maximum difference between
    // the start time and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time, IOProc);

    std::string name = "test_aprox_rates.";

    // Write a plotfile
    WriteSingleLevelPlotfile(name + eos_name, state, vars.names, geom, time, 0);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;

}
