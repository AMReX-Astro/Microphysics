import StarKillerMicrophysics as SKM
from StarKiller.interfaces import EosType
import os.path

class Eos(object):
    def __init__(self):
        self.EosModule = SKM.Eos_Module()
        self.name = self.EosModule.get_eos_name().decode("ASCII").strip().lower()

    def evaluate(self, input_mode, input_state, use_raw_inputs=False):
        self.EosModule.eos(input_mode, input_state.state, use_raw_inputs)

    def _initialize_safe(self):
        if (self.name == "helmholtz"):
            # check to see if the Helmholtz table is in the
            # current working directory
            try:
                assert(os.path.isfile("helm_table.dat"))
            except AssertionError:
                print("helm_table.dat file or symlink is missing from the working directory")
                raise
            except:
                raise
        return True
