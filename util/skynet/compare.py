#!/usr/bin/env python
from __future__ import print_function
from SkyNet import *
import SkyNet
import sys
import os
import subprocess
import argparse
import glob
import numpy as np
from cycler import cycler
import matplotlib.pyplot as plt

# Get Star Killer output

class starkiller_data():
  # modified from Microphysics/unit_test/burn_cell/burn_cell.py
  def __init__(self, inpath, runprefix='react_urca'):
      files = glob.glob(inpath + runprefix + r'_[0-9]*')
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
	self.short_spec_names = data[0]['short_spec_names']

	# Init time, temp, ener
	fnum = []
	time = []
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
	    time.append(d['time'])
	    for i, xi in enumerate(d['xn']):
		xn[i].append(xi)
	    for i, xi in enumerate(d['ydot']):
		ydot[i].append(xi)

	# Convert data to numpy arrays
	for i in range(nspec):
	    xn[i] = np.array(xn[i])
        self.X = xn
	#for i in range(neqs):
	#    ydot[i] = np.array(ydot[i])
	#    ydot[i][0] = ydot[i][1]
	self.fnum = np.array(fnum)
	self.temp = np.array(temp)
	self.time = np.array(time)
	self.dtime = np.ediff1d(time, to_begin=0)
	self.ener = np.array(ener)
	denerdt = np.zeros_like(ener)
	denerdt[1:] = self.ener[1:]/self.dtime[1:]
	self.edot = denerdt
	# for now ydot is garbage might be useful later?


class skydata():

  def __init__(self, inpath, pressure = False): 

      out = NetworkOutput.ReadFromFile(inpath)
      self.time = np.array(out.Times())
      self.temp = np.array(out.TemperatureVsTime()) * 1.0E9
      self.rho  = np.array(out.DensityVsTime())
      self.edot = np.array(out.HeatingRateVsTime())
      self.Y = np.array(out.YVsTime())
      self.A = np.array(out.A())
      self.Z = np.array(out.Z())
      self.X = np.array([self.A * y for y in self.Y[:]])
      
      if pressure:
	  jmax = len(self.temp)

          # writing intermediate thermo file for EoS
          # calls the Helmholtz EoS for pressure as a.out 
	  f = open("tmp.txt", 'w')
	  f.write(str(jmax) + "\n")

	  for i in range(jmax):
	      wbar = 1.0 / np.sum(self.Y[i,:])
	      abar = wbar * np.sum(np.dot(self.A,self.Y[i,:]))
	  zbar = wbar * np.sum(np.dot(self.Z,self.Y[i,:]))
	  f.write("{:.15e}   {:.15e}   {:.15e}   {:.15e}\n".format(self.temp[i], 
	    self.rho[i], abar, zbar))
	  f.close()

	  process = subprocess.Popen("./a.out", shell=True)
	  process.wait()
	  self.pressure = np.loadtxt("tmp_out.txt")
	  os.system("rm tmp*")
	
  def isoidx(self, a, z):
        '''
        return index of isotope in the self.z, self.z, self.isos, self.misos,
        self.xisos[i,:] arrays
        '''
        wh = np.where((self.A == a) & (self.Z == z))[0]
        if len(wh) == 0:
            print('isotope not found')
        else:
            return wh[0]
      

#sky = skydata(inpath = "./urca_abovethresh/urca_isoc_abovethresh.h5", pressure = False)
#star = starkiller_data("./urca_abovethresh/", runprefix='react_urca')
