import StarKillerMicrophysics as SKM
from StarKiller.eos import Eos

def starkiller_initialize(probin_filename):
    eos = Eos()
    with eos._initialize_safe():
        skinit = SKM.Starkiller_Initialization_Module()
        skinit.starkiller_initialize(probin_filename)
