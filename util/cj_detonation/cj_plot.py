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

ax.scatter(1/rho, p, marker="x", color="r")

ax.plot([1.0/rho, 1.0/rho_det], [p, p_det])

ax.set_xscale("log")
ax.set_yscale("log")

fig.savefig("cj_det.png")

