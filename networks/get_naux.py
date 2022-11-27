#!/usr/bin/env python3

import os
import argparse

from general_null import network_param_file


def get_naux(net_file, defines):
    """read through the list network parameter file and output the
    number of aux variables"""

    species = []
    extra_species = []
    aux_vars = []

    # read the species defined in the net_file
    network_param_file.parse(species, extra_species, aux_vars,
                             net_file, defines)

    print(len(aux_vars))


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--defines", type=str, default="",
                        help="preprocessor defines used in building the code")
    parser.add_argument("--microphysics_path", type=str, default="",
                        help="path to Microphysics/")
    parser.add_argument("--net", type=str, default="",
                        help="name of the network")

    args = parser.parse_args()

    micro_path = args.microphysics_path
    net = args.net

    net_file = os.path.join(micro_path, "networks", net, f"{net}.net")
    if not os.path.isfile(net_file):
        net_file = os.path.join(micro_path, "networks", net, "pynucastro.net")

    get_naux(net_file, args.defines)


if __name__ == "__main__":
    main()
