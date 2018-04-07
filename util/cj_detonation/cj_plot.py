import matplotlib.pyplot as plt
import numpy as np

hug = np.loadtxt("hugoniot.txt")

# read in the headings
with open("hugoniot.txt", "r") as f:
    line = f.readline()
    rho = float(line.split("=")[-1])

    line = f.readline()
    p = float(line.split("=")[-1])

    line = f.readline()
    rho_det = float(line.split("=")[-1])

    line = f.readline()
    p_det = float(line.split("=")[-1])


fig = plt.figure(1)
ax = fig.add_subplot(111)

ax.plot(1.0/hug[:, 0], hug[:, 2], lw=2, color="C0", label="detonation adiabat")
ax.plot(1.0/hug[:, 0], hug[:, 1], lw=1, color="C1", label="shock adiabat (no combustion)")

# draw the Rayleigh line
ax.scatter(1/rho, p, marker="x", color="r")
ax.scatter(1/rho_det, p_det, marker="x", color="r")

# draw the Rayleigh line
slope = (p_det - p)/(1/rho_det - 1/rho)
v = 1/hug[:, 0]
ps = p + slope * (v - 1/rho)

ax.plot(v, ps, color="0.5", ls=":", label="Rayleigh line")

ax.set_xlim(0.2e-7, 1.1e-7)
ax.set_ylim(5.e23, 4.e25)

ax.legend(frameon=False)

#ax.set_xscale("log")
#ax.set_yscale("log")

ax.set_xlabel(r"$1/\rho$")
ax.set_ylabel(r"$p$")

fig.savefig("cj_det.png")

