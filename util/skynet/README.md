This directory contains:

"compare.py" - loads the output of skynet and starkiller networks
                     for comparison tests
"burn.py"    - runs skynet in parallel on problems stated therein

NOTE: The "presure" option in the skydata class in compare.py can be
      used if the helmoholtz eos executable for finding pressure is
      located in the directory. The helmholtz eos is freely available
      from: http://cococubed.asu.edu/code_pages/eos.shtml   

Quick Start:
  1. Make sure Skynet is installed and working. Be sure to set your
     SKYNET_HOME variable, and to prepend SKYNET_HOME to your
     PYTHON_PATH. The above depend on importing buit in scripts from
     Skynet, so we need to tell python where to find them.

  2. Change desired problem parameters in burn.py as needed.

     >> python burn.py
     
  3. Load results using the sky class in "compare.py"

     >> import compare as com
     >> data = com.skydata(/path/to/h5/file/)
     >> edot = data.edot









