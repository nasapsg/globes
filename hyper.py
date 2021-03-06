# ---------------------------------------------------------------
# Script to plot hyperspectra as computed by PSG/GlobES
# Villanueva - NASA Goddard Space Flight Center
# March 2021
# ---------------------------------------------------------------
import struct, os
import numpy as np
import matplotlib.pyplot as plt

file = 'data/psg_img.dat'

# Read the binary hyperspectral cube
# The pos array has 5 numbers (0:longitude[deg], 1:latitude[deg], 2:altitude[km], 3:xpos[radius], 4:ypos[radius])
fsz = os.path.getsize(file)
fr = open(file,'rb')
npts = struct.unpack('i', fr.read(4))[0]
xval = struct.unpack('f'*npts, fr.read(4*npts))
nspec = int((fsz-npts*4-4)/(4.0*(npts+5)))
pos = np.zeros([nspec,5])
data = np.zeros([nspec,npts])
for i in range(nspec):
    pos[i,:] = struct.unpack('f'*5, fr.read(4*5))
    data[i,:] = struct.unpack('f'*npts, fr.read(4*npts))
fr.close()

# Plot the spectra and points on the map
pl,ax = plt.subplots(1,3, figsize=(13, 4))
ax[0].plot(pos[:,0],pos[:,1],'o')
ax[0].set_xlabel('Longitude')
ax[0].set_ylabel('Latitude')
ax[1].plot(pos[:,3],pos[:,4],'o')
ax[1].set_xlabel('Position [radius]')
ax[1].set_ylabel('Position [radius]')
cm = plt.get_cmap('inferno')
for i in range(nspec):
    trace = plt.plot(xval, data[i,:])
    cindex = int(230.0*(pos[i,1]+90.0)/180.0)
    trace[0].set_color(cm(cindex))
#Endfor
plt.tight_layout()
plt.show()
