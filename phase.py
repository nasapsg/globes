# ---------------------------------------------------------------
# Script to compute phase curves with PSG/GlobES
# Villanueva, Suissa - NASA Goddard Space Flight Center
# February 2021
# ---------------------------------------------------------------
import sys, os
import numpy as np
import matplotlib.pyplot as plt

# Parameters of the phase
ncfile = 'data/gcm_exocam.nc' # Name of the netCDF file
update = True        # Flag if we want to update the PSG/GlobES simulations or simply read/plot the spectra
phase1 = 0.0         # Initial phase (degrees) for the simulation, 0 is sub-solar point, 180 is night-side
phase2 = 360.0       # Final phase (degrees)
dphase = 90.0        # Phase step (degrees)
binning= 200         # Binning applied to the GCM data for each radiative-transfer (greater is faster, minimum is 1)
lam1   = 0.4         # Initial wavelength of the simulations (um)
lam2   = 30.0        # Final wavelength of the simulations (um)
lamRP  = 500.0       # Resolving power
radunit= 'ppm'       # Desired radiation unit (https://psg.gsfc.nasa.gov/helpmodel.php#units)
#psgurl = 'http://localhost:3000' # URL of the PSG server - For PSG/Docker
psgurl = 'http://localhost' # URL of the PSG server

# Convert netCDF file to PSG/GCM format
from gcm_exocam import convertgcm
convertgcm(ncfile, 'gcm_psg.dat')

# GlobES/API calls can be sequentially, and PSG will remember the previous values
# This means that we can upload parameters step-by-step. To reset your config for GlobES (use type=set), and to simply update (use type=upd)
if update: os.system('curl -s -d app=globes -d type=set --data-urlencode file@gcm_psg.dat %s/api.php' % psgurl)

# Define parameters of this run
fr = open("config.txt", "w")
fr.write('<GENERATOR-RANGE1>%f\n' % lam1)
fr.write('<GENERATOR-RANGE2>%f\n' % lam2)
fr.write('<GENERATOR-RANGEUNIT>um\n')
fr.write('<GENERATOR-RESOLUTION>%f\n' % lamRP)
fr.write('<GENERATOR-RESOLUTIONUNIT>RP\n')
fr.write('<GENERATOR-RADUNITS>%s\n' % radunit)
fr.write('<GENERATOR-RADUNITS>%s\n' % radunit)
fr.write('<GENERATOR-GCM-BINNING>%d' % binning)
fr.close()
if update: os.system('curl -s -d app=globes -d type=upd --data-urlencode file@config.txt %s/api.php' % psgurl)

# Calculate the spectra across the phases
if not os.path.isdir('spectra'): os.system('mkdir spectra')
for phase in np.arange(phase1,phase2+dphase,dphase):
    fr = open("config.txt", "w")
    fr.write('<OBJECT-SEASON>%f\n' % phase)
    fr.write('<OBJECT-OBS-LONGITUDE>%f\n' % phase)
    fr.close()
    if phase>178 and phase<182: phase=182 # Add transit phase
    if update: os.system('curl -s -d app=globes --data-urlencode file@config.txt %s/api.php > spectra/phase%d.txt' % (psgurl,phase))
    data = np.genfromtxt('spectra/phase%d.txt' % phase)
    plt.plot(data[:,0],data[:,1],label=phase)
    print(phase)
#Endfor phase
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.show()
