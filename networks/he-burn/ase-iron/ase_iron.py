# an approximate network for He/C burning with key rates
# to bypass the C12(a,g)O16 rate.
# This is fully-compatible with NSE.

import pynucastro as pyna
from pynucastro.networks import AmrexAstroCxxNetwork

DO_DERIVED_RATES = True

def get_library():

    reaclib_lib = pyna.ReacLibLibrary()

    all_reactants = ["p",
                     "he4", "c12", "o16", "ne20", "mg24", "si28", "s32",
                     "ar36", "ca40", "ti44", "cr48", "fe52", "ni56",
                     "al27", "p31", "cl35", "k39", "sc43", "v47", "mn51", "co55",
                     "n13", "na23"]

    subch = reaclib_lib.linking_nuclei(all_reactants)

    # in this list, we have the reactants, the actual reactants,
    # and modified products that we will use instead
    other_rates = [("c12(c12,n)mg23", "mg24"),
                   ("o16(o16,n)s31", "s32"),
                   ("o16(c12,n)si27", "si28")]

    for r, mp in other_rates:
        _r = reaclib_lib.get_rate_by_name(r)
        forward_rate = pyna.ModifiedRate(_r, new_products=[mp])
        derived_rate = pyna.DerivedRate(forward_rate, use_pf=True)
        subch += pyna.Library(rates=[forward_rate, derived_rate])

    # C12+Ne20 and reverse
    # (a,g) links between Na23 and Al27
    # (a,g) links between Al27 and P31

    rates_to_remove = ["p31(p,c12)ne20",
                       "si28(a,c12)ne20",
                       "ne20(c12,p)p31",
                       "ne20(c12,a)si28",
                       "na23(a,g)al27",
                       "al27(g,a)na23",
                       "al27(a,g)p31",
                       "p31(g,a)al27"]

    for r in rates_to_remove:
        print("removing: ", r)
        _r = subch.get_rate_by_name(r)
        subch.remove_rate(_r)

    # iron group
    iron_peak = ["n", "p", "he4",
                 "mn51",
                 "fe52", "fe53", "fe54", "fe55", "fe56",
                 "co55", "co56", "co57",
                 "ni56", "ni57", "ni58"]
    subch += reaclib_lib.linking_nuclei(iron_peak)
    weak_lib = pyna.TabularLibrary(ordering=["ffn", "langanke", "oda"])
    iron_weak_lib = weak_lib.linking_nuclei(set(iron_peak + all_reactants))
    subch += iron_weak_lib

    if DO_DERIVED_RATES:
        rates_to_derive = subch.backward().get_rates()

        # now for each of those derived rates, look to see if the pair exists

        for r in rates_to_derive:
            fr = subch.get_rate_by_nuclei(r.products, r.reactants)
            if fr:
                print(f"modifying {r} from {fr}")
                subch.remove_rate(r)
                d = pyna.DerivedRate(fr, use_pf=True)
                subch.add_rate(d)

    subch.eliminate_duplicates(rate_type_preference="tabular")

    return subch


def create_network():

    subch = get_library()

    net = AmrexAstroCxxNetwork(libraries=[subch])
    net.make_ap_pg_approx(intermediate_nuclei=["cl35", "k39", "sc43", "v47"])
    net.remove_nuclei(["cl35", "k39", "sc43", "v47"])

    net.make_nn_g_approx(intermediate_nuclei=["fe53", "fe55", "ni57"])
    net.remove_nuclei(["fe53", "fe55", "ni57"])

    return net


def doit():

    net = create_network()

    net.summary()

    comp = pyna.Composition(net.get_nuclei())
    comp.set_all(0.1)
    comp.set_nuc("he4", 0.95)
    comp.normalize()

    rho = 1.e7
    T = 6.e9

    net.plot(outfile="ase-iron.png",
             rotated=True, hide_xalpha=True,
             size=(1500, 450),
             node_size=600, node_font_size=9,
             Z_range=(1, 29))

    net.write_network()


if __name__ == "__main__":
    doit()
