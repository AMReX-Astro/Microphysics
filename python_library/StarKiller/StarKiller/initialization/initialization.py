import StarKillerMicrophysics as SKM
from StarKiller.eos import Eos

def starkiller_initialize(probin_filename):
    with Eos():
        skinit = SKM.Starkiller_Initialization_Module()
        skinit.starkiller_initialize(probin_filename)
