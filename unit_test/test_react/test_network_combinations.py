#!/usr/bin/env python3

from __future__ import print_function

import itertools
import os
import sys

params_file = r"""
&PROBIN
  test_set = "gr0_3d.small"

  xin_file   = "xin.aprox13"
  run_prefix = "react_aprox13_"

  @@test-params@@

/
"""

executable = ""

link_files = ["xin.aprox13",
              "gr0_3d.small",
              "helm_table.dat"]


# this dictionary holds the parameters we want to set.  For each key,
# we use a list to give all the possible values (even if there is only
# a single value). 
params = {
    "dens_min": [1.e4, 1.e6],
    "dens_max": [1.e8],
    "temp_min": [5.e7],
    "temp_max": [5.e9],
    "tmax": [1.e-6, 1.e-3],
    "rtol_spec": [1.e-12, 1.e-8],
    "atol_spec": [1.e-12, 1.e-8],
    "jacobian": [1, 2]
}

# itertools.product() will produce lists with every possible
# combination from the input lists.  Here we use as the input lists
# the values from our dictionary.  This magic comes from:
# http://stackoverflow.com/questions/3873654/combinations-from-dictionary-with-list-values-using-python

combinations = [[{k: v} for (k, v) in zip(params.keys(), values)] 
                for values in itertools.product(*params.values())]


run_no = 0

top_dir = os.getcwd()

for c in combinations:

    cparams = {k: v for d in c for k, v in d.items()}

    # run this combination of test parameters

    # make the directory
    run_no += 1
    odir = "{:02d}".format(run_no)

    try:
        os.mkdir(odir)
    except:
        sys.exit("unable to create directory")

    # copy the executable and suport files
    
    for f in link_files:
        try: 
            os.symlink(os.path.join(top_dir, f), 
                       os.path.join(odir, os.path.basename(f)))
        except:
            sys.exit("unable to link file")

    # write the input file
    infile = "{}/inputs.{}".format(odir, odir)
    with open(infile, "w") as f:
        for line in params_file.splitlines():
            if line.find("@@test-params@@") >= 0:
                for k, v in cparams.items():
                    f.write("  {} = {}\n".format(k, v))
            else:
                f.write("{}\n".format(line))

    # run
    

    # log the result


