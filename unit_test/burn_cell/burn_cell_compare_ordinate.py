#!/usr/bin/env python
from __future__ import print_function
import argparse
import glob
import numpy as np
from cycler import cycler
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('runprefixes', type=str, nargs='+',
                    help='Prefixes of the ordinate and time txt files to compare.')
parser.add_argument('--ordinate', type=str, help='Name of the ordinate variable to plot vs time (in the txt file name)',
                    default='denerdt')
parser.add_argument('--logtime', action='store_true', help='If --logtime, plot Log10(time).')
parser.add_argument('--tlo', type=float, help='Time lower limit')
parser.add_argument('--thi', type=float, help='Time upper limit')
args = parser.parse_args()

ofiles = []
tfiles = []
for prefix in args.runprefixes:
    ofiles += glob.glob('{}_{}.txt'.format(prefix, args.ordinate))
    tfiles += glob.glob('{}_time.txt'.format(prefix))

if len(ofiles) == 0 or len(tfiles) == 0:
    exit()

ordinate_data = [np.genfromtxt(of) for of in ofiles]
time_data = [np.genfromtxt(tf) for tf in tfiles]

## Define RGBA to HEX
def rgba_to_hex(rgba):
    r = int(rgba[0]*255.0)
    g = int(rgba[1]*255.0)
    b = int(rgba[2]*255.0)
    return '#{:02X}{:02X}{:02X}'.format(r,g,b)

# Figure out time axis limits
if args.tlo and args.thi:
    ltlim = [args.tlo, args.thi]
elif args.tlo:
    ltlim = [args.tlo, time_data[0][-1]]
elif args.thi:
    ltlim = [time_data[0][0], args.thi]
else:
    ltlim = [time_data[0][0], time_data[0][-1]]
if args.logtime:
    time_data = [np.log10(td) for td in time_data]
    ltlim = np.log10(ltlim)

xlim = ltlim
if args.logtime:
    xlabel = '$\\mathrm{Log_{10}~Time~(s)}$'
else:
    xlabel = '$\\mathrm{Time~(s)}$'
    
# Get set of colors to use for abundances
cm = plt.get_cmap('nipy_spectral')
clist = [cm(1.0*i/len(time_data)) for i in range(len(time_data))]
hexclist = [rgba_to_hex(ci) for ci in clist]

# Plot ordinate vs. time
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_prop_cycle(cycler('color', hexclist))
for i, (od, td) in enumerate(zip(ordinate_data, time_data)):
    ax.plot(td, np.log10(od), label='$\\mathrm{' + '{}'.format(args.runprefixes[i].replace('_', '\_') + '}$'))
lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
ax.set_xlim(xlim)
plt.xlabel(xlabel)
plt.ylabel('$\\mathrm{' + '{}'.format(args.ordinate) + '}$')
plt.savefig('cf_{}_v_time.eps'.format(args.ordinate),
            bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.savefig('cf_{}_v_time.png'.format(args.ordinate), dpi=300,
            bbox_extra_artists=(lgd,), bbox_inches='tight')

