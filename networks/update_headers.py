#!/usr/bin/env python3

import os
import argparse
from pathlib import Path

from general_null import write_network


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--microphysics_path", type=str, default="",
                        help="path to Microphysics/")
    parser.add_argument("--net", type=str, default="",
                        help="name of the network")
    parser.add_argument("--net-custom-path", type=str, default="",
                        help="full path to the network if it lives outside of Microphysics")
    parser.add_argument("--odir", type=str, default="",
                        help="output directory")
    parser.add_argument("--defines", type=str, default="",
                        help="any preprocessor defines")

    args = parser.parse_args()

    micro_path = args.microphysics_path
    net = args.net


    if len(args.net_custom_path.strip()) == 0:
        # our network lives in Microphysics/networks
        base_path = Path(micro_path) / "networks" / f"{net}"
    else:
        # our network is standalone (outside of Microphysics/)
        base_path = Path(args.net_custom_path)

    net_file = base_path / "pynucastro.net"
    if not net_file.exists() and len(net) > 0:
        net_file = base_path / f"{net}.net"

    cxx_template = os.path.join(micro_path, "networks",
                                "general_null/network_header.template")
    cxx_name = os.path.join(args.odir, "network_properties.H")

    try:
        os.makedirs(args.odir)
    except FileExistsError:
        pass

    print(f"calling with {cxx_template}, {net_file}")
    write_network.write_network(cxx_template,
                                net_file,
                                cxx_name, args.defines)


if __name__ == "__main__":
    main()
