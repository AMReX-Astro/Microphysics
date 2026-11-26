#!/usr/bin/env python

import argparse

import numpy as np
import matplotlib.pyplot as plt


def doit(datafile):

    # read just the first line of the datafile to get the variable
    # name
    with open(datafile) as df:
        line = df.readline()
        varname = line.split()[-1]

    data = np.loadtxt(datafile)

    # get the unique density and temperatures
    Ts = np.log10(np.unique(data[:, 0]))
    nT = len(Ts)

    rhos = np.log10(np.unique(data[:, 1]))
    nrho = len(rhos)

    # reshape the EOS quantity

    var = data[:, 2].reshape(nT, nrho)

    # if the variable spans more than 3 decades, then log it
    logged = False
    if abs(var.max() / var.min()) > 1000:
        var = np.log10(var)
        logged = True

    fig, ax = plt.subplots(layout="constrained")

    im = ax.imshow(var,
                   origin="lower",
                   extent=[rhos.min(), rhos.max(), Ts.min(), Ts.max()])

    ax.set_xlabel(r"$\rho$ [g/cm$^3$]")
    ax.set_ylabel("T [K]")

    orientation = "horizontal" if nrho > nT else "vertical"

    label = varname
    if logged:
        label = f"log10({label})"

    fig.colorbar(im, ax=ax, label=label, orientation=orientation)

    aspect = nT / nrho

    width = 7
    height = width * aspect

    if nrho > nT:
        # extra for colorbar
        height += 1
    else:
        width += 1

    fig.set_size_inches(width, height)

    fig.savefig(f"eos_{varname}.png")



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="input filename")
    args = parser.parse_args()

    doit(args.filename)
