#!/usr/bin/env python3

import os

from general_null import write_network



def main():

    # find all of the networks to update -- we assume that they will
    # have a .net file with the same name as the network, e.g.,
    # aprox13/aprox13.net
    net_files = []
    for d in os.listdir():
        if os.path.isfile("{}/{}.net".format(d, d)):
            net_files.append(d)

    print(net_files)

    cwd = os.getcwd()

    fortran_template = os.path.join(cwd, "general_null/network_properties.template")
    cxx_template = os.path.join(cwd, "general_null/network_header.template")

    for net in net_files:
        os.chdir(net)

        net_file = "{}.net".format(net)

        write_network.write_network(fortran_template, cxx_template,
                                    net_file,
                                    "network_properties.F90", "network_properties.H")

        os.chdir(cwd)

if __name__ == "__main__":
    main()
