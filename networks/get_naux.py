#!/usr/bin/env python3

import os
import argparse
from pathlib import Path

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
    parser.add_argument("--net-custom-path", type=str, default="",
                        help="full path to the network if it lives outside of Microphysics")

    args = parser.parse_args()

    micro_path = args.microphysics_path
    net = args.net

    if args.net_custom_path is None:
        base_path = Path(micro_path) / "networks" / "net"
    else:
        base_path = Path(args.net_custom_path)

    net_file = base_path / "pynucastro.net"
    if not net_file.exists() and len(net) > 0:
        net_file = base_path / f"{net}.net"

    get_naux(net_file, args.defines)


if __name__ == "__main__":
    main()
