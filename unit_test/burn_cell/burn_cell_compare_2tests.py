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

# Initializing varibales and loading in data

print('Initializing')

runprefix = args.runprefix

file_testprefixes = open('{}_testprefixes.txt'.format(runprefix), 'r')
testprefixes = []
for line in file_testprefixes:
    testprefixes.append('{}'.format(line.strip()))
file_testprefixes.close()

short_spec_names = []
nspec = []
inputs = []
for i in range(len(testprefixes)):
    # i corresponds to the index of a test prefix
    inputs.append([])
    file_inputs = open('{}_{}_inputs.txt'.format(runprefix, testprefixes[i]))
    for line in file_inputs:
        inputs[i].append('{}'.format(line.strip()))
    file_inputs.close()
    short_spec_names.append([])
    file_specs = open('{}_{}_short_spec_names.txt'.format(runprefix, testprefixes[i]), 'r')
    for line in file_specs:
        short_spec_names[i].append(line.strip())
    file_specs.close()
    nspec.append(len(short_spec_names[i]))

# Init time, temp, ener, xn, ydot
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
    if prefix==0:
        for n in range(nspec[0]):
            xn[prefix].append(np.loadtxt('{}_{}_xn{}.txt'.format(args.runprefix, testprefixes[prefix], n)))
            ydot[prefix].append(np.loadtxt('{}_{}_ydot{}.txt'.format(args.runprefix, testprefixes[prefix], n)))
    if prefix==1:
        for n in range(nspec[1]):
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
clist1 = [cm(1.0*i/nspec[0]) for i in range(nspec[0])]
hexclist1 = [rgba_to_hex(ci) for ci in clist1]
clist2 = [cm(1.0*i/nspec[1]) for i in range(nspec[1])]
hexclist2 = [rgba_to_hex(ci) for ci in clist2]

# Initialize plots
plt.figure(1, figsize=(6,9))
ax1 = plt.subplot(211)
ax1.set_prop_cycle(cycler('color', hexclist1))
ax2 = plt.subplot(212)
ax2.set_prop_cycle(cycler('color', hexclist2))
plt.figure(2, figsize=(6,9))
ay1 = plt.subplot(211)
ay1.set_prop_cycle(cycler('color', hexclist1))
ay2 = plt.subplot(212)
ay2.set_prop_cycle(cycler('color', hexclist2))
plt.figure(3, figsize=(5,9))
aT1 = plt.subplot(211)
aT2 = plt.subplot(212)
plt.figure(4, figsize=(5,9))
ae1 = plt.subplot(211)
ae2 = plt.subplot(212)

# Initialize arrays to contain values for plotting
diffx = []
diffydot = []
difftemp = []
diffdenerdt = []

# Plotting the density data
for x in range(nspec[0]):
    # x corresponds to each molecule in the list of species
    plt.figure(1)
    ax1.semilogy(xvec1, xn[0][x], label='{}-{}'.format(short_spec_names[0][x], testprefixes[0]))
    plt.figure(2)
    ay1.semilogy(xvec1, ydot[0][x], label='{}-{}'.format(short_spec_names[0][x], testprefixes[0]))
for x in range(nspec[1]):
    plt.figure(1)
    ax2.semilogy(xvec2, xn[1][x], label='{}-{}'.format(short_spec_names[1][x], testprefixes[1]))
    plt.figure(2)
    ay2.semilogy(xvec2, ydot[1][x], label='{}-{}'.format(short_spec_names[1][x], testprefixes[1]))
plt.figure(3)
aT1.semilogy(xvec1, temp[0], label=testprefixes[0])
aT2.semilogy(xvec2, temp[1], label=testprefixes[1])
plt.figure(4)
ae1.semilogy(xvec1, denerdt[0], label=testprefixes[0])
ae2.semilogy(xvec2, denerdt[1], label=testprefixes[1])

# Mass Fraction Figure
print('Compiling Mass Fraction graph.')
plt.figure(1)
box = ax1.get_position()
ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax1.legend(loc='upper left', bbox_to_anchor=(1,1), fontsize = 5)
ax1.text(0.005, 0.005, '{}    {}'.format(inputs[0][30], inputs[0][31]), fontsize=5, transform=ax1.transAxes)
ax1.set_xlabel(xlabel, fontsize=10)
ax1.set_ylabel('$\\mathrm{Log_{10} X}$', fontsize=10)
ax1.set_title('Mass Fraction for {}'.format(testprefixes[0]), fontsize=15)
ax1.set_xlim(xlim1)
ax1.tick_params(axis='both', which='both', labelsize=5)
box = ax2.get_position()
ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax2.legend(loc='upper left', bbox_to_anchor=(1,1), fontsize = 5)
ax2.text(0.005, 0.005, '{}    {}'.format(inputs[1][30], inputs[1][31]), fontsize=5, transform=ax2.transAxes)
ax2.set_xlabel(xlabel, fontsize=10)
ax2.set_ylabel('$\\mathrm{Log_{10} X}$', fontsize=10)
ax2.set_title('Mass Fraction for {}'.format(testprefixes[1]), fontsize=15)
ax2.set_xlim(xlim2)
ax2.tick_params(axis='both', which='both', labelsize=5)
plt.savefig('{}_{}_xn_compare_2networks.png'.format(runprefix, testprefixes[0]), dpi=700)

# Moller Fractions
print('Compiling Moller Fraction graph.')
plt.figure(2)
box = ay1.get_position()
ay1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ay1.legend(loc='upper left', bbox_to_anchor=(1,1), fontsize = 5)
ay1.text(0.005, 0.005, '{}    {}'.format(inputs[0][30], inputs[0][31]), fontsize=5, transform=ay1.transAxes)
ay1.set_xlabel(xlabel, fontsize=10)
ay1.set_ylabel('$\\mathrm{Log_{10} \\dot{Y}}$', fontsize=10)
ay1.set_title('Moller Fraction for {}'.format(testprefixes[0]), fontsize=15)
ay1.set_xlim(xlim1)
ay1.tick_params(axis='both', which='both', labelsize=5)
box = ay2.get_position()
ay2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ay2.legend(loc='upper left', bbox_to_anchor=(1,1), fontsize = 5)
ay2.text(0.005, 0.005, '{}    {}'.format(inputs[1][30], inputs[1][31]), fontsize=5, transform=ay2.transAxes)
ay2.set_xlabel(xlabel, fontsize=10)
ay2.set_ylabel('$\\mathrm{Log_{10} \\dot{Y}}$', fontsize=10)
ay2.set_title('Moller Fraction for {}'.format(testprefixes[1]), fontsize=15)
ay2.set_xlim(xlim2)
ay2.tick_params(axis='both', which='both', labelsize=5)
plt.savefig('{}_{}_y_compare_2networks.png'.format(runprefix, testprefixes[0]), dpi=700)

# Temperature Figure
print('Compiling Temperature graph.')
plt.figure(3)
aT1.legend(fontsize = 5, loc = 'lower right')
aT1.text(0.005, 0.005, '{}    {}'.format(inputs[0][30], inputs[0][31]), fontsize=5, transform=aT1.transAxes)
aT1.set_xlabel(xlabel, fontsize=10)
aT1.set_ylabel('$\\mathrm{Log_{10} T~(K)}$', fontsize=10)
aT1.set_title('Temperature for {}'.format(testprefixes[0]), fontsize=15)
aT1.set_xlim(xlim1)
aT1.tick_params(axis='both', which='both', labelsize=5)
aT2.legend(fontsize = 5, loc = 'lower right')
aT2.text(0.005, 0.005, '{}    {}'.format(inputs[1][30], inputs[1][31]), fontsize=5, transform=aT2.transAxes)
aT2.set_xlabel(xlabel, fontsize=10)
aT2.set_ylabel('$\\mathrm{Log_{10} T~(K)}$', fontsize=10)
aT2.set_title('Temperature for {}'.format(testprefixes[1]), fontsize=15)
aT2.set_xlim(xlim2)
aT2.tick_params(axis='both', which='both', labelsize=5)
plt.savefig('{}_{}_T_compare_2networks.png'.format(runprefix, testprefixes[0]), dpi=700)

# Energy Generation Rate
print('Compiling Enerergy Generation Rate graph.')
plt.figure(4)
ae1.legend(fontsize = 5, loc = 'lower right')
ae1.text(0.005, 0.005, '{}    {}'.format(inputs[0][30], inputs[0][31]), fontsize=5, transform=ae1.transAxes)
ae1.set_xlabel(xlabel, fontsize=10)
ae1.set_ylabel('$\\mathrm{Log_{10} \\dot{e}~(erg/g/s)}$', fontsize=10)
ae1.set_title('Energy Generation Rate for {}'.format(testprefixes[0]), fontsize=15)
ae1.set_xlim(xlim1)
ae1.tick_params(axis='both', which='both', labelsize=5)
ae2.legend(fontsize = 5, loc = 'lower right')
ae2.text(0.005, 0.005, '{}    {}'.format(inputs[1][30], inputs[1][31]), fontsize=5, transform=ae2.transAxes)
ae2.set_xlabel(xlabel, fontsize=10)
ae2.set_ylabel('$\\mathrm{Log_{10} \\dot{e}~(erg/g/s)}$', fontsize=10)
ae2.set_title('Energy Generation Rate for {}'.format(testprefixes[1]), fontsize=15)
ae2.set_xlim(xlim2)
ae2.tick_params(axis='both', which='both', labelsize=5)
plt.savefig('{}_{}_edot_compare_2networks.png'.format(runprefix, testprefixes[0]), dpi=700)
