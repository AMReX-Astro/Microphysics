# C-burning rate module generator

import reaclib

files = ["c12-ag-o16-nac2"]

rc = reaclib.RateCollection(files)

rc.make_network('boxlib')




