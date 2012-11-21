#!/usr/bin/env python

# Example illustrating the application of MBAR to compute a 1D PMF from an umbrella sampling simulation.
#
# number of bins is max([4,20]) such that no bin is empty
#
# usage: ./umbrella-sampling-density.py (n_umbrellas) (max_traj_length) (min_val) (max_val) (centers_file)
#
# files: 
#    centers.dat 
#        [spring center] [spring constant] [path/to/trajfile.dat]
#    data_dir/trajfile.dat
#        [time] [val] [other vals, to be discarded]
#
# REFERENCE
# 
# D. L. Mobley, A. P. Graves, J. D. Chodera, A. C. McReynolds, B. K. Shoichet and K. A. Dill, "Predicting absolute ligand binding free energies to a simple model site," Journal of Molecular Biology 371(4):1118-1134 (2007).
# http://dx.doi.org/10.1016/j.jmb.2007.06.002

import numpy # numerical array library
#from math import *
import sys
import pymbar # multistate Bennett acceptance ratio
import timeseries # timeseries analysis
import math


# Code is meant to be used as somewhat of a black box but you will need to
# modify the following parameters (at the very least) for use with your code:
###########################################

# Parameters
K = 8                       # number of umbrellas
N_max = 100000              # maximum number of snapshots/simulation
val_min = 1                 # min for PMF
val_max = 25                # max for PMF
centersf = "umbrellas"      # file with list of brolly centers and paths to trajectory files
n_skip = 0                  # skip the first n datapoints

###########################################


# Temperatures
kB = 1.0
temperature = 1.0 # assume a single temperature -- can be overridden with data from center.dat 
T_k = numpy.ones(K,float)*temperature # inital temperatures are all equal 
beta = 1.0 / (kB * temperature) # inverse temperature of simulations (in 1/(kJ/mol))

# Allocate storage for simulation data
N_k = numpy.zeros([K], numpy.int32) # N_k[k] is the number of snapshots from umbrella simulation k
K_k = numpy.zeros([K], numpy.float64) # K_k[k] is the spring constant (in kJ/mol/deg**2) for umbrella simulation k
val0_k = numpy.zeros([K], numpy.float64) # val0_k[k] is the spring center location (in deg) for umbrella simulation k
val_kn = numpy.zeros([K,N_max], numpy.float64) # val_kn[k,n] is the torsion angle (in deg) for snapshot n from umbrella simulation k
u_kn = numpy.zeros([K,N_max], numpy.float64) # u_kn[k,n] is the reduced potential energy without umbrella restraints of snapshot n of umbrella simulation k
u_kln = numpy.zeros([K,K,N_max], numpy.float64) # u_kln[k,l,n] is the reduced potential energy of snapshot n from umbrella simulation k evaluated at umbrella l
g_k = numpy.zeros([K],numpy.float32)
trajf_k = []

# Read in umbrella spring constants and centers.
infile = open(centersf, 'r')
lines = infile.readlines()
infile.close()
for k in range(K):
    # Parse line k.
    line = lines[k]
    tokens = line.split()
    val0_k[k] = float(tokens[0]) # spring center location 
    K_k[k] = float(tokens[1]) # spring constant
    #trajf_k.append(str(tokens[2])) # infile name
    trajf_k.append("window_N=" + tokens[0] + ".txt")
        #  if len(tokens) > 3:
        #      T_k[k] = float(tokens[2])  # temperature the kth simulation was run at.

beta_k = 1.0/(kB*T_k)   # beta factor for the different temperatures
DifferentTemperatures = True
if (min(T_k) == max(T_k)):
    DifferentTemperatures = False            # if all the temperatures are the same, then we don't have to read in energies.
# Read the simulation data
for k in range(K):
    # Read torsion angle data.
    filename = trajf_k[k]
    print "Reading %s..." % filename
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
    # Parse data.
    n = 0
    for line in lines:
        if line[0] != '#' and line[0] != '@':
            if n >= n_skip:
                tokens = line.split()
                val = float(tokens[1])
                val_kn[k,n-n_skip] = val
            n += 1
    N_k[k] = n-n_skip

   # if (DifferentTemperatures):  # if different temperatures are specified the metadata file, 
                                 # then we need the energies to compute the PMF
        # Read energies
   #     filename = 'data/prod%d_energies.xvg' % k
   #     print "Reading %s..." % filename
   #     infile = open(filename, 'r')
   #     lines = infile.readlines()
   #     infile.close()
        # Parse data.
   #     n = 0
   #     for line in lines:
   #         if line[0] != '#' and line[0] != '@':
   #             tokens = line.split()            
   #             u_kn[k,n] = beta_k[k] * (float(tokens[2]) - float(tokens[1])) # reduced potential energy without umbrella restraint
   #             n += 1

    # Compute correlation times for potential energy and val
    # timeseries.  If the temperatures differ, use energies to determine samples; otherwise, use the cosine of val
            
    if (DifferentTemperatures):        
        g_k[k] = timeseries.statisticalInefficiency(u_kn[k,:], u_kn[k,:])
        print "Correlation time for set %5d is %10.3f" % (k,g_k[k])
        indices = timeseries.subsampleCorrelatedData(u_kn[k,:])
    else:
        #g_k[k] = timeseries.statisticalInefficiency(val_kn[k,:], val_kn[k,:])
        #print "Correlation time for set %5d is %10.3f" % (k,g_k[k])
        indices = timeseries.subsampleCorrelatedData(val_kn[k,0:N_k[k]], fast=True, verbose=True)
    # Subsample data.
    N_k[k] = len(indices)
    u_kn[k,0:N_k[k]] = u_kn[k,indices]
    val_kn[k,0:N_k[k]] = val_kn[k,indices]
   # print val_kn[k,0:N_k[k]]

# Set zero of u_kn -- this is arbitrary.
u_kn -= u_kn.min()

val_min = numpy.min([numpy.min(val_kn[k,0:N_k[k]]) for k in range(K)])
val_max = numpy.max([numpy.max(val_kn[k,0:N_k[k]]) for k in range(K)])

# Test bin numbers
good_bins = []
for tmpnbins in range(5,50,5):
    print "Testing nbins =", tmpnbins
    delta = (val_max - val_min) / float(tmpnbins)
    # compute bin centers
    bin_center_i = numpy.zeros([tmpnbins], numpy.float64)
    for i in range(tmpnbins):
        bin_center_i[i] = val_min + delta/2 + delta * i
    # Bin data
    bin_kn = numpy.zeros([K,numpy.max(N_k)], numpy.int32)
    for k in range(K):
        for n in range(N_k[k]):
            # Compute bin assignment.
            bin_kn[k,n] = int((val_kn[k,n] - val_min) / delta)
    # Test for empty bins
    isgood = True
    for i in range(tmpnbins):
     #   print bin_center_i[i], numpy.sum(bin_kn==i)
        if numpy.sum(bin_kn==i) < 5:
            isgood = False
    if isgood:
        good_bins.append(tmpnbins)

if len(good_bins) > 0:
    print "good nbins =", good_bins, ", choosing max =", max(good_bins)
    nbins = max(good_bins)
else:
    print "found no good nbins, try resetting min and max?"
    print "real min", numpy.min([numpy.min(val_kn[k,0:N_k[k]]) for k in range(K)]), "vs given min", val_min
    print "real max", numpy.max([numpy.max(val_kn[k,0:N_k[k]]) for k in range(K)]), "vs given max", val_max
    sys.exit()

# Construct bins
print "Binning data for realz..."
delta = (val_max - val_min) / float(nbins)
# compute bin centers
bin_center_i = numpy.zeros([nbins], numpy.float64)
for i in range(nbins):
    bin_center_i[i] = val_min + delta/2 + delta * i
# Bin data
bin_kn = numpy.zeros([K,numpy.max(N_k)], numpy.int32)
for k in range(K):
    for n in range(N_k[k]):
        # Compute bin assignment.
        bin_kn[k,n] = int((val_kn[k,n] - val_min) / delta)
# Test for empty bins
for i in range(nbins):
    print bin_center_i[i], numpy.sum(bin_kn==i)

# Evaluate reduced energies in all umbrellas
print "Evaluating reduced potential energies..."
for k in range(K):
    for n in range(N_k[k]):
        for l in range(K):
            # Compute deviation from umbrella center l
            dval = val_kn[k,n] - val0_k[l]

            # Compute energy of snapshot n from simulation k in umbrella potential l
            u_kln[k,l,n] = u_kn[k,n] + beta_k[k] * (K_k[l]/2.0) * dval**2

# Initialize MBAR.
print "Running MBAR..."
mbar = pymbar.MBAR(u_kln, N_k, verbose = True, method = 'adaptive', initialize='BAR', relative_tolerance=1.0e-6)
#mbar = pymbar.MBAR(u_kln, N_k, verbose = True, method = 'Newton-Raphson')

# Compute PMF in unbiased potential (in units of kT).
(f_i, df_i) = mbar.computePMF(u_kn, bin_kn, nbins)

# Write out PMF
print "PMF (in units of kT)"
print "%8s %8s %8s" % ('bin', 'f', 'df')
for i in range(nbins):
    print "%8.3f %8.3f %8.3f" % (bin_center_i[i], f_i[i], df_i[i])

##############################################

pKfile = open("logP(N).txt", 'w')
for i in range(nbins):
    pKfile.write(str(bin_center_i[i]) + "\t" + str(-f_i[i]) + "\n")
pKfile.close()

