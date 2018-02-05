# This code is designed to compare the absolute difference between one 
# reference burn_cell test and multiple other burn_cell tests.
# burn_cell_testing.py must be run before running this.

#!/usr/bin/env python
from __future__ import print_function
import argparse
import glob
import numpy as np
from cycler import cycler
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('runprefix', type=str,
                    help='Prefix of the output run files. We look for files named as [prefix]_[0-9]*')
parser.add_argument('--filenum', action='store_true', help='If --filenum, plot vs. file number')
parser.add_argument('--logtime', action='store_true', help='If --logtime, plot Log10(time).')
parser.add_argument('--tlo', type=float, help='Time lower limit')
parser.add_argument('--thi', type=float, help='Time upper limit')
parser.add_argument('--nlo', type=float, help='File num lower limit')
parser.add_argument('--nhi', type=float, help='File num upper limit')
args = parser.parse_args()

# INITIALIZED VARIABLES

runprefix = args.runprefix

file_specs = open('{}_short_spec_names.txt'.format(runprefix), 'r')
short_spec_names = []
for line in file_specs:
    short_spec_names.append(line.strip())
file_specs.close()

nspec = len(short_spec_names)

file_testprefixes = open('{}_testprefixes.txt'.format(runprefix), 'r')
testprefixes = []
for line in file_testprefixes:
    testprefixes.append('{}'.format(line.strip()))
file_testprefixes.close()

inputs = []
for i in range(len(testprefixes)):
    # i corresponds to the index of a test prefix
    inputs.append([])
    file_inputs = open('{}_{}_inputs.txt'.format(runprefix, testprefixes[i]))
    for line in file_inputs:
        inputs[i].append('{}'.format(line.strip()))
    file_inputs.close()

# Init time, temp, ener
xn = []
ydot = []
fnum = []
temp = []
dtime = []
time = []
ener = []
denerdt = []

for prefix in range(len(testprefixes)):
    xn.append([])
    ydot.append([])
    for n in range(nspec):
        xn[prefix].append(np.loadtxt('{}_{}_xn{}.txt'.format(args.runprefix, testprefixes[prefix], n)))
        ydot[prefix].append(np.loadtxt('{}_{}_ydot{}.txt'.format(args.runprefix, testprefixes[prefix], n)))
    temp.append(np.loadtxt('{}_{}_temp.txt'.format(args.runprefix, testprefixes[prefix])))
    ener.append(np.loadtxt('{}_{}_ener.txt'.format(args.runprefix, testprefixes[prefix])))
    denerdt.append(np.loadtxt('{}_{}_denerdt.txt'.format(args.runprefix, testprefixes[prefix])))

dtime = np.loadtxt('{}_{}_dtime.txt'.format(args.runprefix, testprefixes[0]))
time = np.loadtxt('{}_{}_time.txt'.format(args.runprefix, testprefixes[0]))
fnum = np.loadtxt('{}_{}_fnum.txt'.format(args.runprefix, testprefixes[0]))

## Define RGBA to HEX
def rgba_to_hex(rgba):
    r = int(rgba[0]*255.0)
    g = int(rgba[1]*255.0)
    b = int(rgba[2]*255.0)
    return '#{:02X}{:02X}{:02X}'.format(r,g,b)

## PLOTTING

# Figure out time axis limits
if args.tlo and args.thi:
    ltlim = [args.tlo, args.thi]
elif args.tlo:
    ltlim = [args.tlo, time[-1]]
elif args.thi:
    ltlim = [time[0], args.thi]
else:
    ltlim = [time[0], time[-1]]
if args.logtime:
    time = np.log10(time)
    ltlim = np.log10(ltlim)

# Number axis limits
if args.nlo and args.nhi:
    fnlim = [args.nlo, args.nhi]
elif args.tlo:
    fnlim = [args.nlo, fnum[-1]]
elif args.thi:
    fnlim = [fnum[0], args.nhi]
else:
    fnlim = [fnum[0], fnum[-1]]

# Time or file number selection
if args.filenum or args.nlo or args.nhi:
    plot_vs_fnum = True
    xlabel = r'$\mathrm{Output \#}$'
    xvec = fnum
    xlim = fnlim
else:
    xvec = time
    xlim = ltlim
    if args.logtime:
        xlabel = r'$\mathrm{Log_{10}}$'
    else:
        xlabel = r'$\mathrm{Time~(s)}$'
    
# Get set of colors to use for abundances
cm = plt.get_cmap('nipy_spectral')
clist = [cm(1.0*i/4) for i in range(4)]
hexclist = [rgba_to_hex(ci) for ci in clist]

# Plot difference in X vs. time
plt.figure(1, figsize=(5,9))
ax = plt.subplot(211)
errx = plt.subplot(212)
plt.figure(2, figsize=(5,9)) 
ay = plt.subplot(211)
erry = plt.subplot(212)
plt.figure(3, figsize=(5,9))
aT = plt.subplot(211)
errT = plt.subplot(212)
plt.figure(4, figsize=(5,9))
ae = plt.subplot(211)
erre = plt.subplot(212)

diffx = []
diffydot = []
difftemp = []
diffdenerdt = []

line_styles = ['solid', 'dashed', 'dotted', 'dashdot']

# Plotting the reference data
print('Plotting the reference data: {}'.format(testprefixes[0]))
for x in range(len(short_spec_names)):
    # x corresponds to each molecule in the list of species
    plt.figure(1)
    ax.semilogy(xvec, xn[0][x], label='{}-{}'.format(short_spec_names[x], testprefixes[0]), linestyle = line_styles[0])
    plt.figure(2)
    ay.semilogy(xvec, ydot[0][x], label='{}-{}'.format(short_spec_names[x], testprefixes[0]), linestyle = line_styles[0])
plt.figure(3)
aT.semilogy(xvec, temp[0], label=testprefixes[0], linestyle = line_styles[0])
plt.figure(4)
ae.semilogy(xvec, denerdt[0], label=testprefixes[0], linestyle = line_styles[0])

# Plotting the data compared to reference and the error
print('Plotting the compared data and the errors.')
for i in range(1, len(testprefixes)):   
    # In this context i cooresponds to a test prefix to be compared
    #    to the data from a chosen data set
    print('Plotting for {}'.format(testprefixes[i]))
    difftemp.append([])
    diffdenerdt.append([])
    for n in range(len(xvec)):
        # n is for every time step from 0 to tmax
        difftemp[i-1].append(abs(temp[0][n] - temp[i][n]))
        diffdenerdt[i-1].append(abs(denerdt[0][n] - denerdt[i][n]))
    plt.figure(3)
    # Uncomment the following line and the commented ae, ax, and ay
    #   to add additional graphs to the top graph in the output files
    #aT.semilogy(xvec, temp[i], label=testprefixes[i], linestyle = line_styles[i])
    errT.semilogy(xvec, difftemp[i-1], label=testprefixes[i], linestyle = line_styles[i-1])
    plt.figure(4)
    #ae.semilogy(xvec, denerdt[i], label=testprefixes[i], linestyle = line_styles[i])
    erre.semilogy(xvec, diffdenerdt[i-1], label=testprefixes[i], linestyle = line_styles[i-1])
    diffx.append([])
    diffydot.append([])
    for x in range(3):
        # x is for each species involved
        diffx[i-1].append([])
        diffydot[i-1].append([])
        for n in range(len(xvec)):
            # n is for every time step from 0 to tmax
            diffx[i-1][x].append(abs(xn[0][x][n] - xn[i][x][n]))
            diffydot[i-1][x].append(abs(ydot[0][x][n] - ydot[i][x][n]))
        plt.figure(1)
        #ax.semilogy(xvec, xn[i][x], label='{}-{}'.format(short_spec_names[x], testprefixes[i]), linestyle = line_styles[i])
        errx.semilogy(xvec, diffx[i-1][x], label='{}-{}'.format(short_spec_names[x], testprefixes[i]), linestyle = line_styles[i-1])
        plt.figure(2)
        #ay.semilogy(xvec, ydot[i][x], label='{}-{}'.format(short_spec_names[x], testprefixes[i]), linestyle = line_styles[i])
        erry.plot(xvec, diffydot[i-1][x], label='{}-{}'.format(short_spec_names[x], testprefixes[i]), linestyle = line_styles[i-1])

# Mass Fraction Figure
print('Compiling Mass Fraction graph.')
plt.figure(1)
ax.legend(fontsize = 5, loc = 'lower right')
ax.text(0.005, 0.005, '{}    {}'.format(inputs[0][30], inputs[0][31]), fontsize=5, transform=ax.transAxes)
ax.set_prop_cycle(cycler('color', hexclist))
ax.set_xlabel(xlabel, fontsize=10)
ax.set_ylabel('$\\mathrm{Log_{10} X}$', fontsize=10)
ax.set_title('Mass Fraction')
ax.set_xlim(xlim)
ax.tick_params(axis='both', which='major', labelsize=5)
errx.legend(fontsize = 5, loc = 'upper left')
errx.set_prop_cycle(cycler('color', hexclist))
errx.set_xlabel(xlabel, fontsize=10)
errx.set_ylabel('$\\mathrm{Log_{10} X}$', fontsize=10)
errx.set_title('Absolute Errors in Mass Fraction', fontsize=15)
errx.set_xlim(xlim)
errx.tick_params(axis='both', which='major', labelsize=5)
plt.savefig('{}_{}_xn_compare_abs.png'.format(runprefix, testprefixes[0]), dpi=700)

# Moller Fractions
print('Compiling Moller Fraction graph.')
plt.figure(2)
ay.legend(fontsize = 5, loc = 'lower right')
ay.text(0.005, 0.005, '{}    {}'.format(inputs[0][30], inputs[0][31]), fontsize=5, transform=ay.transAxes)
ay.set_prop_cycle(cycler('color', hexclist))
ay.set_xlabel(xlabel, fontsize=10)
ay.set_ylabel('$\\mathrm{Log_{10} \\dot{Y}}$', fontsize=10)
ay.set_title('Moller Fraction')
ay.set_xlim(xlim)
ay.tick_params(axis='both', which='major', labelsize=5)
erry.legend(fontsize = 5, loc = 'lower right')
erry.set_prop_cycle(cycler('color', hexclist))
erry.set_xlabel(xlabel, fontsize=10)
erry.set_ylabel('$\\mathrm{Log_{10} \\dot{Y}}$', fontsize=10)
erry.set_title('Absolute Errors in Moller Fraction', fontsize=15)
erry.set_xlim(xlim)
erry.tick_params(axis='both', which='major', labelsize=5)
plt.savefig('{}_{}_y_compare_abs.png'.format(runprefix, testprefixes[0]), dpi=700)

# Temperature Figure
print('Compiling Temperature graph.')
plt.figure(3)
aT.legend(fontsize = 5, loc = 'lower right')
aT.text(0.005, 0.005, '{}    {}'.format(inputs[0][30], inputs[0][31]), fontsize=5, transform=aT.transAxes)
aT.set_prop_cycle(cycler('color', hexclist))
aT.set_xlabel(xlabel, fontsize=10)
aT.set_ylabel('$\\mathrm{Log_{10} T~(K)}$', fontsize=10)
aT.set_title('Temperature')
aT.set_xlim(xlim)
aT.tick_params(axis='both', which='major', labelsize=5)
errT.legend(fontsize = 5, loc = 'lower right')
errT.set_prop_cycle(cycler('color', hexclist))
errT.set_xlabel(xlabel, fontsize=10)
errT.set_ylabel('$\\mathrm{Log_{10} T~(K)}$', fontsize=10)
errT.set_title('Absolute Error in Temperature', fontsize=15)
errT.set_xlim(xlim)
errT.tick_params(axis='both', which='major', labelsize=5)
plt.savefig('{}_{}_T_compare_abs.png'.format(runprefix, testprefixes[0]), dpi=700)

# Energy Generation Rate
print('Compiling Enerergy Generation Rate graph.')
plt.figure(4)
ae.legend(fontsize = 5, loc = 'lower right')
ae.text(0.005, 0.005, '{}    {}'.format(inputs[0][30], inputs[0][31]), fontsize=5, transform=ae.transAxes)
ae.set_prop_cycle(cycler('color', hexclist))
ae.set_xlabel(xlabel, fontsize=10)
ae.set_ylabel('$\\mathrm{Log_{10} \\dot{e}~(erg/g/s)}$', fontsize=10)
ae.set_title('Energy Generation Rate')
ae.set_xlim(xlim)
ae.tick_params(axis='both', which='major', labelsize=5)
erre.legend(fontsize = 5, loc = 'lower right')
erre.set_prop_cycle(cycler('color', hexclist))
erre.set_xlabel(xlabel, fontsize=10)
erre.set_ylabel('$\\mathrm{Log_{10} \\dot{e}~(erg/g/s)}$', fontsize=10)
erre.set_title('Absolute Error in Energy Generation Rate', fontsize=15)
erre.set_xlim(xlim)
erre.tick_params(axis='both', which='major', labelsize=5)
plt.savefig('{}_{}_edot_compare_abs.png'.format(runprefix, testprefixes[0]), dpi=700)
