# This code operates the same way as burn_cell.py except it will save all 
# arrays created to a batch of files named <runprefix>_<testprefix>___.txt
# This can be used with other codes that open these saved files to process
# and analyze the information they hold. 

#!/usr/bin/env python
from __future__ import print_function
import argparse
import glob
import numpy as np
from cycler import cycler
import matplotlib.pyplot as plt
import os

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

files = glob.glob(args.runprefix + r'_output/' + args.runprefix + r'_[0-9]*')
print('Found {} files matching pattern {}'.format(len(files), args.runprefix+'_[0-9]*'))
if len(files) == 0:
    exit()

data = []

for fn in files:
    d = {}
    f = open(fn, 'r')

    # Get file number
    fnsplit = fn.split('_')
    fnum = int(fnsplit[-1])
    d['fnum'] = fnum
    
    # Eat lines
    flines = []
    for l in f:
        flines.append(l.strip())
    f.close()

    # Fill dataset
    n = len(flines)
    for i, xi in enumerate(flines):
        if xi[-1] == ':':
            # This is a field
            key = xi[:-1]
            xd = []
            # Collect following lines of field data
            for j in range(i+1,n):
                xj = flines[j]
                if xj[-1] == ':':
                    break
                else:
                   xd.append(xj.split())
                   
            # Interpret data
            if (key=='nspec' or key=='neqs' or key=='i' or key=='j' or key=='k'
                or key=='n_rhs' or key=='n_jac'):
                # Integer
                d[key] = int(xd[0][0])
            elif key=='short_spec_names':
                # Strings of nuclide abbreviations
                d[key] = [xj[0] for xj in xd]
            elif key=='self_heat':
                # Boolean
                d[key] = (xd[0][0]=='T')
            elif key=='xn' or key=='ydot':
                # 1-D Float Array
                d[key] = np.array([float(y) for y in xd[0]])
            elif key=='jac':
                # 2-D Float Array
                d[key] = np.array([[float(xij) for xij in xdi] for xdi in xd])
            else:
                # Float
                d[key] = float(xd[0][0])

    # Store data
    data.append(d)

## SORT file data by file number
data_sorted = sorted(data, key=lambda d: d['fnum'])
data = data_sorted
    
## INIT VARIABLES
nspec = data[0]['nspec']
neqs  = data[0]['neqs']
short_spec_names = data[0]['short_spec_names']

# Init time, temp, ener
fnum = []
dtime = []
temp = []
ener = []

# Init abundance vectors
xn = [[] for i in range(nspec)]

# Init ydot vectors
ydot = [[] for i in range(neqs)]

## DATA LOOP
# Loop through data and collect
for d in data:
    fnum.append(d['fnum'])
    temp.append(d['T'])
    ener.append(d['e'])
    dtime.append(d['time'])
    for i, xi in enumerate(d['xn']):
        xn[i].append(xi)
    for i, xi in enumerate(d['ydot']):
        ydot[i].append(xi)

# Convert data to numpy arrays
for i in range(nspec):
    xn[i] = np.array(xn[i])
for i in range(neqs):
    ydot[i] = np.array(ydot[i])
    ydot[i][0] = ydot[i][1]
fnum = np.array(fnum)
temp = np.array(temp)
dtime = np.array(dtime)
time = np.cumsum(dtime)
ener = np.array(ener)
denerdt = np.zeros_like(ener)
denerdt[1:] = ener[1:]/dtime[1:]

# Obtain test prefix
testprefix = input('Please input a testing prefix:   ')

#Save Data to a file
for n in range(nspec):
    np.savetxt('{}_{}_xn{}.txt'.format(args.runprefix, testprefix, n), xn[n])
    np.savetxt('{}_{}_ydot{}.txt'.format(args.runprefix, testprefix, n), ydot[n])
np.savetxt('{}_{}_fnum.txt'.format(args.runprefix, testprefix), fnum)
np.savetxt('{}_{}_temp.txt'.format(args.runprefix, testprefix), temp)
np.savetxt('{}_{}_dtime.txt'.format(args.runprefix, testprefix), dtime)
np.savetxt('{}_{}_time.txt'.format(args.runprefix, testprefix), time)
np.savetxt('{}_{}_ener.txt'.format(args.runprefix, testprefix), ener)
np.savetxt('{}_{}_denerdt.txt'.format(args.runprefix, testprefix), denerdt)
np.savetxt('{}_{}_fnum.txt'.format(args.runprefix, testprefix), fnum)

#Save species names to file
specfile = open('{}_{}_short_spec_names.txt'.format(args.runprefix, testprefix), 'w')
for item in short_spec_names:
    specfile.write('{}\n'.format(item))
specfile.close()

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
    xlabel = '$\\mathrm{Output \\#}$'
    xvec = fnum
    xlim = fnlim
else:
    xvec = time
    xlim = ltlim
    if args.logtime:
        xlabel = '$\\mathrm{Log_{10}~Time~(s)}$'
    else:
        xlabel = '$\\mathrm{Time~(s)}$'
    
# Get set of colors to use for abundances
cm = plt.get_cmap('nipy_spectral')
clist = [cm(1.0*i/nspec) for i in range(nspec)]
hexclist = [rgba_to_hex(ci) for ci in clist]

# Plot X vs. time
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_prop_cycle(cycler('color', hexclist))
for i in range(nspec):    
    ax.plot(xvec, np.log10(xn[i]), label=short_spec_names[i])
lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
ax.set_xlim(xlim)
plt.xlabel(xlabel)
plt.ylabel('$\\mathrm{Log_{10} X}$')
plt.savefig(args.runprefix+'_logX.eps',
            bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.savefig(args.runprefix+'_logX.png', dpi=300,
            bbox_extra_artists=(lgd,), bbox_inches='tight')

# Clear figure
plt.clf()

# Plot T, edot vs. time

# Find where edot = 0
def y_where_x_zero(y, x):
    yzero = []
    xiszero = False
    ylo = 0.0
    for yi, xi in zip(y, x):
        if xi == 0.0 and not xiszero:
            xiszero = True
            ylo = yi
        if xi != 0.0 and xiszero:
            xiszero = False
            yzero.append([ylo, yi])
    if xiszero:
        yzero.append([ylo, y[-1]])
    return yzero

edotzero = y_where_x_zero(time, denerdt)
ax = fig.add_subplot(111)
ax.plot(xvec, np.log10(temp), label='Temperature', color='red')
lgd = ax.legend(bbox_to_anchor=(1.1, 1), loc=2, borderaxespad=0.0)
ax.set_xlabel(xlabel)
ax.set_ylabel('$\\mathrm{Log_{10} T~(K)}$')
ax.set_xlim(xlim)
ax2 = ax.twinx()
ax2.plot(xvec, np.log10(denerdt), label='E Gen Rate', color='blue')
ax2.set_ylabel('$\\mathrm{Log_{10} \\dot{e}~(erg/g/s)}$')
ax2.set_xlim(xlim)
lgd2 = ax2.legend(bbox_to_anchor=(1.1, 0.75), loc=2, borderaxespad=0.0)
# hatch where edot=0
for edz in edotzero:
    plt.axvspan(np.log10(edz[0]), np.log10(edz[1]), color='blue', fill=False, linewidth=0, hatch='/', alpha=0.2)
plt.savefig(args.runprefix+'_'+testprefix+'_T-edot.eps',
            bbox_extra_artists=(lgd,lgd2,), bbox_inches='tight')
plt.savefig(args.runprefix+'_'+testprefix+'_T-edot.png', dpi=300,
            bbox_extra_artists=(lgd,lgd2,), bbox_inches='tight')


# Clear figure
plt.clf()

# Plot Nuclide dY/dt vs time
ax = fig.add_subplot(111)

# Separate Ydot tracks into +/- sets
class Track(object):
    def __init__(self):
        self.xvec = []
        self.yval = []
        self.sign = 1
        self.color = ''
        self.linestyle = ''
        self.label = ''
        return

tracks = []
for i in range(neqs-2):
    firsttrack = True
    scratch = Track()
    newtrack = True
    
    for xp, yp in zip(xvec, ydot[i]):
        if (not newtrack) and (np.sign(yp) != scratch.sign):
            # Store track and create new
            tracks.append(scratch)
            scratch = Track()
            newtrack = True
            
        if newtrack:
            # set color, sign, linestyle, etc
            newtrack = False
            scratch.color = hexclist[i]
            scratch.sign = np.sign(yp)
            if firsttrack:
                scratch.label = short_spec_names[i]
                firsttrack = False
            else:
                scratch.label = None
            if scratch.sign == -1:
                scratch.linestyle = 'dashed'
            else:
                scratch.linestyle = 'solid'

        if np.sign(yp) == scratch.sign:
            # append to this track
            scratch.xvec.append(xp)
            scratch.yval.append(np.absolute(yp))
    # Append final track for this species
    tracks.append(scratch)
            
for trc in tracks:
    ax.plot(trc.xvec, np.log10(trc.yval), label=trc.label,
            color=trc.color, linestyle=trc.linestyle)
lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
ax.set_xlim(xlim)
plt.xlabel(xlabel)
plt.ylabel('$\\mathrm{Log_{10} \\dot{Y}}$')
plt.savefig(args.runprefix+'_'+testprefix+'_ydot.eps',
            bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.savefig(args.runprefix+'_'+testprefix+'_ydot.png', dpi=300,
            bbox_extra_artists=(lgd,), bbox_inches='tight')
