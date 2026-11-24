import argparse

import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np


def doit(state_file, xmin, n_plot, outfile):

    data = np.genfromtxt(state_file, names=True)

    # find the keys of the most abundant species
    maxX = {}
    for k in data.dtype.names:
        if k in ["Time", "Temperature"]:
            continue
        maxX[k] = data[k].max()

    abundant = [k for k, v in sorted(maxX.items(), key=lambda item: item[1])][::-1]

    fig , (ax1, ax2) = plt.subplots(2, 1, sharex=True)

    ax1.plot(data["Time"], data["Temperature"])

    # the first time point is t = 0, which doesn't work well
    # on a logscale for x, so we use symlog to transition
    # to log at the second time point
    ax1.set_xscale("symlog", linthresh=data["Time"][1])
    ax1.grid(ls=":")

    ax1.set_ylabel("temperature [K]")
    ax1.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))

    # now the species
    for n in range(n_plot):
        ax2.plot(data["Time"], data[abundant[n]], label=abundant[n])

    ax2.set_xscale("symlog", linthresh=data["Time"][1])
    ax2.set_yscale("log")
    ax2.set_ylim(xmin, 1.5)
    ax2.grid(ls=":")

    ax2.set_xlabel("time (s)")
    ax2.set_ylabel("mass fraction")

    ax2.legend(fontsize="small")

    fig.tight_layout()

    fig.set_size_inches(6, 7)

    fig.savefig(outfile)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--xmin", type=float, default=1.e-10,
                        help="cutoff value for mass fractions")
    parser.add_argument("-n", type=int, default=5,
                        help="number of species to plot (most abundant shown)")
    parser.add_argument("-o", type=str, default="state.png",
                        help="name of plot output file")
    parser.add_argument("state_file", type=str, nargs=1,
                        help="burn cell state history file")

    args = parser.parse_args()

    doit(args.state_file[0], args.xmin, args.n, args.o)
