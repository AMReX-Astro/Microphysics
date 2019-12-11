import StarKillerMicrophysics as SKM
from StarKiller.eos import Eos
from StarKiller.network import Network

def starkiller_initialize(probin_filename):
    try:
        Eos._initialize_safe()
    except:
        print("EOS cannot be initialized.")
        raise

    try:
        Network._initialize_safe()
    except:
        print("Network cannot be initialized.")
        raise

    skinit = SKM.Starkiller_Initialization_Module()
    skinit.starkiller_initialize(probin_filename)

    eos = Eos()
    net = Network()

    print("\nInitialized StarKiller with ...")
    print("- EOS:     {}".format(eos.name))
    print("- Network: {}".format(net.name))
