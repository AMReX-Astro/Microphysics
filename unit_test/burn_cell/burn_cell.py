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
args = parser.parse_args()

files = glob.glob(args.runprefix + r'_[0-9]*')
print('Found {} files matching pattern {}'.format(len(files), args.runprefix+'_[0-9]*'))

data = []

for fn in files:
    d = {}
    f = open(fn, 'r')
    
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

## INIT VARIABLES
nspec = data[0]['nspec']
neqs  = data[0]['neqs']
short_spec_names = data[0]['short_spec_names']

# Init time
dtime = []

# Init abundance vectors
xn = [[] for i in range(nspec)]

# Init ydot vectors
ydot = [[] for i in range(neqs)]

## DATA LOOP
# Loop through data and collect
for d in data:
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
dtime = np.array(dtime)
time = np.cumsum(dtime)

## Define RGBA to HEX
def rgba_to_hex(rgba):
    r = int(rgba[0]*255.0)
    g = int(rgba[1]*255.0)
    b = int(rgba[2]*255.0)
    return '#{:02X}{:02X}{:02X}'.format(r,g,b)

## PLOTTING
# Get set of colors to use for abundances
cm = plt.get_cmap('nipy_spectral')
clist = [cm(1.0*i/nspec) for i in range(nspec)]
hexclist = [rgba_to_hex(ci) for ci in clist]

# Plot X vs. time
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_prop_cycle(cycler('color', hexclist))
for i in range(nspec):    
    ax.plot(time, np.log10(xn[i]), label=short_spec_names[i])
lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
plt.xlabel('Time (s)')
plt.ylabel('$\\mathrm{Log_{10} X}$')
plt.savefig(args.runprefix+'_abundances.eps',
            bbox_extra_artists=(lgd,), bbox_inches='tight')

# Clear figure
plt.clf()

# Plot Nuclide dY/dt vs time
ax = fig.add_subplot(111)

# Separate Ydot tracks into +/- sets
class Track(object):
    def __init__(self):
        self.time = []
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
    
    for t, yp in zip(time, ydot[i]):
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
            scratch.time.append(t)
            scratch.yval.append(np.absolute(yp))
    # Append final track for this species
    tracks.append(scratch)
            
for trc in tracks:
    ax.plot(trc.time, np.log10(trc.yval), label=trc.label,
            color=trc.color, linestyle=trc.linestyle)
lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
ax.set_xscale('log')
plt.xlabel('Time (s)')
plt.ylabel('$\\mathrm{Log_{10} \\dot{Y}}$')
plt.savefig(args.runprefix+'_ydot.eps',
            bbox_extra_artists=(lgd,), bbox_inches='tight')
