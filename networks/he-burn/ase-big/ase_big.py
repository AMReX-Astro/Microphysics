# an approximate network for He/C burning with key rates
# to bypass the C12(a,g)O16 rate.
# This is fully-compatible with NSE.

import pynucastro as pyna
from pynucastro.networks import AmrexAstroCxxNetwork

DO_DERIVED_RATES = True

def get_library():

    reaclib_lib = pyna.ReacLibLibrary()

    all_reactants = ["p", "he4",
                     "c12", "c13", "c14",
                     "n13", "n14", "n15",
                     "o16", "o17", "o18",
                     "f18", "f19",
                     "ne20", "ne21", "ne22",
                     "na22", "na23",
                     "mg24", "mg25",
                     "al26", "al27",
                     "si28", "si29",
                     "p30", "p31",
                     "s32", "s33",
                     "cl34", "cl35",
                     "ar36", "ar37",
                     "k38", "k39",
                     "ca40", "ca41",
                     "sc42", "sc43",
                     "ti44", "ti45",
                     "v46", "v47",
                     "cr48", "cr49",
                     "mn50", "mn51", "mn52",
                     "co55",
                     "fe52", "ni56",
                     ]

    subch = reaclib_lib.linking_nuclei(all_reactants)


    # iron group
    iron_peak = ["n", "p", "he4",
                 "mn50", "mn51", "mn52",
                 "fe52", "fe53", "fe54", "fe55", "fe56",
                 "co55", "co56", "co57",
                 "ni56", "ni57", "ni58", "cu59", "zn60"]
    subch += reaclib_lib.linking_nuclei(iron_peak)
    weak_lib = pyna.TabularWeakLibrary(ordering=["ffn", "langanke", "oda"])
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

    #net.make_nn_g_approx(intermediate_nuclei=["fe53", "fe55", "ni57"])
    #net.remove_nuclei(["fe53", "fe55", "ni57"])

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

    net.plot(outfile="ase-big.png",
             rotated=True, hide_xalpha=True,
             size=(1500, 450),
             node_size=600, node_font_size=9,
             Z_range=(1, 31))

    net.write_network()


if __name__ == "__main__":
    doit()
