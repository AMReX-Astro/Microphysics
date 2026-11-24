import matplotlib.pyplot as plt
import numpy as np

# read in the output from different tolerances
tola = np.loadtxt("state_over_time_tol_a.txt")
tolb = np.loadtxt("state_over_time_tol_b.txt")
tolc = np.loadtxt("state_over_time_tol_c.txt")
ref = np.loadtxt("state_over_time_tol_d.txt")

labela = r"tol = $10^{-3}$"
labelb = r"tol = $10^{-5}$"
labelc = r"tol = $10^{-8}$"

# compute the relative error with respect to the
# tighest tolerance
dta = np.abs(tola - ref) / np.abs(ref)
dtb = np.abs(tolb - ref) / np.abs(ref)
dtc = np.abs(tolc - ref) / np.abs(ref)

# plot the output
fig, ax = plt.subplots()

itemp = 1
ihe4 = 2
ic12 = 3
isi28 = 7
ini56 = 14

fig, axes = plt.subplots(2, 2)

for n, (ax, idx, name) in enumerate(zip(axes.flat, [ihe4, ic12, isi28, ini56], ["He4", "C12", "Si28", "Ni56"])):

    ax.plot(tola[:, 0], dta[:, idx], label=labela)
    ax.plot(tolb[:, 0], dtb[:, idx], label=labelb)
    ax.plot(tolc[:, 0], dtc[:, idx], label=labelc)

    ax.set_xlabel("time (s)")
    ax.set_ylabel(f"relative error {name}")

    ax.set_xscale("log")
    ax.set_yscale("log")

    if n == 0:
        ax.legend(fontsize="small")

fig.tight_layout()
fig.set_size_inches(8, 7)

fig.savefig("tolerance-compare.png")
