import StarKillerMicrophysics as SKM
from StarKiller.interfaces import BurnType

class Integrator(object):
    def __init__(self):
        self.IntegratorModule = SKM.Integrator_Module()

    def integrate(self, state_in, dt):
        # state_in is a BurnType object
        # integrate() returns a BurnType object as state_out
        state_out = state_in.copy()

        # always start integration at t=0
        time = 0.0
        
        # call the integrator
        self.IntegratorModule.integrator(state_in.state,
                                         state_out.state,
                                         dt, time)

        return state_out
