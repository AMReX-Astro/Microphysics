#!/usr/bin/env python3

import os
import argparse

from general_null import write_network


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--microphysics_path", type=str, default="",
                        help="path to Microphysics/")
    parser.add_argument("--net", type=str, default="",
                        help="name of the network")
    parser.add_argument("--odir", type=str, default="",
                        help="output directory")
    parser.add_argument("--defines", type=str, default="",
                        help="any preprocessor defines")

    args = parser.parse_args()

    micro_path = args.microphysics_path
    net = args.net

    net_file = os.path.join(micro_path, "networks", net, f"{net}.net")
    if not os.path.isfile(net_file):
        net_file = os.path.join(micro_path, "networks", net, "pynucastro.net")

    properties_file = os.path.join(micro_path, "networks",
                                   net, "NETWORK_PROPERTIES")

    cxx_template = os.path.join(micro_path, "networks",
                                "general_null/network_header.template")
    cxx_name = os.path.join(args.odir, "network_properties.H")

    try:
        os.makedirs(args.odir)
    except FileExistsError:
        pass

    write_network.write_network(cxx_template,
                                net_file, properties_file,
                                cxx_name, args.defines)


if __name__ == "__main__":
    main()
