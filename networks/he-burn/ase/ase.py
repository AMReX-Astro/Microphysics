import pynucastro as pyna

import he_burn_core

DO_DERIVED_RATES = True


def doit():

    extra_nuclei = ["n", "n14", "f18", "ne21", "na22"]

    subch = he_burn_core.get_core_library(include_n14_approx=False,
                                          include_zn=False,
                                          include_iron_peak=False,
                                          include_low_ye=False,
                                          extra_nuclei=extra_nuclei,
                                          do_detailed_balance=DO_DERIVED_RATES)

    # additional neutron rates to remove
    for r in subch.get_rates():
        if (r == subch.get_rate_by_name("mg24(n,a)ne21") or
            r == subch.get_rate_by_name("ne21(a,n)mg24") or
            r == subch.get_rate_by_name("na22(n,g)na23") or
            r == subch.get_rate_by_name("na23(g,n)na22") or
            r == subch.get_rate_by_name("ne21(g,n)ne20") or
            r == subch.get_rate_by_name("ne20(n,g)ne21")):
            continue

        if pyna.Nucleus("n") in r.reactants or pyna.Nucleus("n") in r.products:
            print("removing neutron rates: ", r)
            subch.remove_rate(r)

    # these are the rates that we are going to allow to be optionally
    # zeroed
    r1 = subch.get_rate_by_name("c12(p,g)n13")
    r2 = subch.get_rate_by_name("n13(he4,p)o16")

    net = pyna.AmrexAstroCxxNetwork(libraries=[subch],
                                    symmetric_screening=False,
                                    disable_rate_params=[r1, r2])

    net.make_ap_pg_approx(intermediate_nuclei=["cl35", "k39", "sc43",
                                               "v47", "mn51", "co55"])
    net.remove_nuclei(["cl35", "k39", "sc43", "v47", "mn51", "co55"])

    net.summary()

    comp = pyna.Composition(net.get_nuclei())
    comp.set_all(0.1)
    comp.set_nuc("he4", 0.95)
    comp.normalize()

    rho = 1.e7
    T = 6.e9

    net.plot(rho, T, comp, outfile="ase.png",
             rotated=True, hide_xalpha=True, curved_edges=True,
             size=(1500, 450),
             node_size=600, node_font_size=11,
             Z_range=(1, 29))

    net.write_network()


if __name__ == "__main__":
    doit()
