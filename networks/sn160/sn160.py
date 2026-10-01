# This builds a subch network with all the
# Reaclib rates linking the specified nuclei.

import pynucastro as pyna


def create_network():

    all_nuclei = ["n",
                  "h1,2",
                  "he3,4",
                  "li6,7",
                  "be7,9",
                  "b8,10-11",
                  "c12-14",
                  "n13-15",
                  "o14-18",
                  "f17-19",
                  "ne18-22",
                  "na21-23",
                  "mg23-26",
                  "al25-27",
                  "si28-32",
                  "p29-33",
                  "s32-36",
                  "cl33-37",
                  "ar36-40",
                  "k37-41",
                  "ca40-48",
                  "sc43-49",
                  "ti44-51",
                  "v46-52",
                  "cr48-54",
                  "mn50-55",
                  "fe52-58",
                  "co53-59",
                  "ni56-64",
                  "cu57-65",
                  "zn59-66",
                  "ga62-64",
                  "ge63-64"]

    net = pyna.network_helper(all_nuclei, network_type="amrex",
                              tabular_ordering=["ffn", "langanke"])

    return net


if __name__ == "__main__":

    net = create_network()
    net.write_network()

    net.summary()
