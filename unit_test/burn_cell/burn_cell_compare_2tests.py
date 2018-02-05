# This code will compare the mass fraction, moller fraction, temp, and density
# of two different burn_cell tests. 
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
time1 = np.loadtxt('{}_{}_time.txt'.format(args.runprefix, testprefixes[0]))
time2 = np.loadtxt('{}_{}_time.txt'.format(args.runprefix, testprefixes[1]))
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
    ltlim1 = [time1[0], time1[-1]]
    ltlim2 = [time2[0], time2[-1]]
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
    xvec1 = time1
    xvec2 = time2
    xlim1 = ltlim1
    xlim2 = ltlim2
    if args.logtime:
        xlabel = r'$\mathrm{Log_{10}}$'
    else:
        xlabel = r'$\mathrm{Time~(s)}$'
    
# Get set of colors to use for abundances
cm = plt.get_cmap('nipy_spectral')
clist = [cm(1.0*i/4) for i in range(4)]
hexclist = [rgba_to_hex(ci) for ci in clist]

# Initialize plots
# The a plots are plots of the top graph and err plots are the bottom graph
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

# Plotting the density data
print('Plotting the reference data: {}'.format(testprefixes[0]))
for x in range(len(short_spec_names)):
    # x corresponds to each molecule in the list of species
    plt.figure(1)
    ax.semilogy(xvec1, xn[0][x], label='{}-{}'.format(short_spec_names[x], testprefixes[0]))
    errx.semilogy(xvec2, xn[1][x], label='{}-{}'.format(short_spec_names[x], testprefixes[1]))
    plt.figure(2)
    ay.semilogy(xvec1, ydot[0][x], label='{}-{}'.format(short_spec_names[x], testprefixes[0]))
    erry.semilogy(xvec2, ydot[1][x], label='{}-{}'.format(short_spec_names[x], testprefixes[1]))
plt.figure(3)
aT.semilogy(xvec1, temp[0], label=testprefixes[0])
errT.semilogy(xvec2, temp[1], label=testprefixes[1])
plt.figure(4)
ae.semilogy(xvec1, denerdt[0], label=testprefixes[0])
erre.semilogy(xvec2, denerdt[1], label=testprefixes[1])

# Mass Fraction Figure
print('Compiling Mass Fraction graph.')
plt.figure(1)
ax.legend(fontsize = 5, loc = 'lower right')
ax.text(0.005, 0.005, '{}    {}'.format(inputs[0][30], inputs[0][31]), fontsize=5, transform=ax.transAxes)
ax.set_prop_cycle(cycler('color', hexclist))
ax.set_xlabel(xlabel, fontsize=10)
ax.set_ylabel('$\\mathrm{Log_{10} X}$', fontsize=10)
ax.set_title('Mass Fraction for {}'.format(testprefixes[0]), fontsize=15)
ax.set_xlim(xlim1)
ax.tick_params(axis='both', which='major', labelsize=5)
errx.legend(fontsize = 5, loc = 'lower right')
errx.text(0.005, 0.005, '{}    {}'.format(inputs[1][30], inputs[1][31]), fontsize=5, transform=errx.transAxes)
errx.set_prop_cycle(cycler('color', hexclist))
errx.set_xlabel(xlabel, fontsize=10)
errx.set_ylabel('$\\mathrm{Log_{10} X}$', fontsize=10)
errx.set_title('Mass Fraction for {}'.format(testprefixes[1]), fontsize=15)
errx.set_xlim(xlim2)
errx.tick_params(axis='both', which='major', labelsize=5)
plt.savefig('{}_{}_xn_compare_rho.png'.format(runprefix, testprefixes[0]), dpi=700)

# Moller Fractions
print('Compiling Moller Fraction graph.')
plt.figure(2)
ay.legend(fontsize = 5, loc = 'lower right')
ay.text(0.005, 0.005, '{}    {}'.format(inputs[0][30], inputs[0][31]), fontsize=5, transform=ay.transAxes)
ay.set_prop_cycle(cycler('color', hexclist))
ay.set_xlabel(xlabel, fontsize=10)
ay.set_ylabel('$\\mathrm{Log_{10} \\dot{Y}}$', fontsize=10)
ay.set_title('Moller Fraction for {}'.format(testprefixes[0]), fontsize=15)
ay.set_xlim(xlim1)
ay.tick_params(axis='both', which='major', labelsize=5)
erry.legend(fontsize = 5, loc = 'lower right')
erry.text(0.005, 0.005, '{}    {}'.format(inputs[1][30], inputs[1][31]), fontsize=5, transform=erry.transAxes)
erry.set_prop_cycle(cycler('color', hexclist))
erry.set_xlabel(xlabel, fontsize=10)
erry.set_ylabel('$\\mathrm{Log_{10} \\dot{Y}}$', fontsize=10)
erry.set_title('Moller Fraction for {}'.format(testprefixes[1]), fontsize=15)
erry.set_xlim(xlim2)
erry.tick_params(axis='both', which='major', labelsize=5)
plt.savefig('{}_{}_y_compare_rho.png'.format(runprefix, testprefixes[0]), dpi=700)

# Temperature Figure
print('Compiling Temperature graph.')
plt.figure(3)
aT.legend(fontsize = 5, loc = 'lower right')
aT.text(0.005, 0.005, '{}    {}'.format(inputs[0][30], inputs[0][31]), fontsize=5, transform=aT.transAxes)
aT.set_prop_cycle(cycler('color', hexclist))
aT.set_xlabel(xlabel, fontsize=10)
aT.set_ylabel('$\\mathrm{Log_{10} T~(K)}$', fontsize=10)
aT.set_title('Temperature for {}'.format(testprefixes[0]), fontsize=15)
aT.set_xlim(xlim1)
aT.tick_params(axis='both', which='major', labelsize=5)
errT.legend(fontsize = 5, loc = 'lower right')
errT.text(0.005, 0.005, '{}    {}'.format(inputs[1][30], inputs[1][31]), fontsize=5, transform=errT.transAxes)
errT.set_prop_cycle(cycler('color', hexclist))
errT.set_xlabel(xlabel, fontsize=10)
errT.set_ylabel('$\\mathrm{Log_{10} T~(K)}$', fontsize=10)
errT.set_title('Temperature for {}'.format(testprefixes[1]), fontsize=15)
errT.set_xlim(xlim2)
errT.tick_params(axis='both', which='major', labelsize=5)
plt.savefig('{}_{}_T_compare_rho.png'.format(runprefix, testprefixes[0]), dpi=700)

# Energy Generation Rate
print('Compiling Enerergy Generation Rate graph.')
plt.figure(4)
ae.legend(fontsize = 5, loc = 'lower right')
ae.text(0.005, 0.005, '{}    {}'.format(inputs[0][30], inputs[0][31]), fontsize=5, transform=ae.transAxes)
ae.set_prop_cycle(cycler('color', hexclist))
ae.set_xlabel(xlabel, fontsize=10)
ae.set_ylabel('$\\mathrm{Log_{10} \\dot{e}~(erg/g/s)}$', fontsize=10)
ae.set_title('Energy Generation Rate for {}'.format(testprefixes[0]), fontsize=15)
ae.set_xlim(xlim1)
ae.tick_params(axis='both', which='major', labelsize=5)
erre.legend(fontsize = 5, loc = 'lower right')
erre.text(0.005, 0.005, '{}    {}'.format(inputs[1][30], inputs[1][31]), fontsize=5, transform=erre.transAxes)
erre.set_prop_cycle(cycler('color', hexclist))
erre.set_xlabel(xlabel, fontsize=10)
erre.set_ylabel('$\\mathrm{Log_{10} \\dot{e}~(erg/g/s)}$', fontsize=10)
erre.set_title('Energy Generation Rate for {}'.format(testprefixes[1]), fontsize=15)
erre.set_xlim(xlim2)
erre.tick_params(axis='both', which='major', labelsize=5)
plt.savefig('{}_{}_edot_compare_rho.png'.format(runprefix, testprefixes[0]), dpi=700)
