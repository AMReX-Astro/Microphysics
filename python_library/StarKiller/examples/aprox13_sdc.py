from StarKiller.initialization import StarKiller_initialize
from StarKiller.integration import SDCOde
from StarKiller.interfaces import BurnType

StarKiller_initialize("probin_aprox13")

sdc = SDCOde()

state_in = BurnType()
state_in.state.rho = 1.0e8
state_in.state.t = 1.0e9
state_in.state.xn[:] = 0.0
state_in.state.xn[0] = 1.0

time = 0.0
dt = 1.0e-6
dt_init = 1.0e-9

state_out = sdc.integrate(state_in, time, dt, dt_init)

print("state in:")
print(state_in.state)

print("state out:")
print(state_out.state)
