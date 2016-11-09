#!/usr/bin/env python3

from __future__ import print_function

import itertools

PARAMS = """
&PROBIN
  test_set = "gr0_3d.small"

  xin_file   = "xin.aprox13"
  run_prefix = "react_aprox13_"

  @@test-params@@

/
"""

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


for c in combinations:

    # run this combination of test parameters

    # make the directory

    # copy the executable and suport files

    # write the input file

    # run

    # log the result


