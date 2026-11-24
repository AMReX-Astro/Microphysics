#!/usr/bin/env python3

import argparse

import network_param_file


def get_naux(net_file, defines):
    """read through the list network parameter file and output the
    number of aux variables"""

    species = []
    extra_species = []
    aux_vars = []

    # read the species defined in the net_file

    network_param_file.parse(species, extra_species,
                             aux_vars, net_file, defines)

    print(len(aux_vars))


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--defines", type=str, default="",
                        help="preprocessor defines used in building the code")
    parser.add_argument("network_file", type=str, nargs=1, default="",
                        help="network file name")

    args = parser.parse_args()

    get_naux(args.network_file[0], args.defines)


if __name__ == "__main__":
    main()
