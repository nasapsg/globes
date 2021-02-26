# ---------------------------------------------------------------
# Script to convert netCDF climate files into PSG GCM files
# netCDF files in LMD climate model
# Villanueva, Fauchez - NASA Goddard Space Flight Center
# February 2021
# ---------------------------------------------------------------
import sys
import numpy as np
from netCDF4 import Dataset as ncdf

# Module that performs the conversion
def convertgcm(filein = 'data/gcm_lmd.nc', fileout = 'gcm_psg.dat', itime=0):
	# Read the netCDF file
	nfile = ncdf(filein)
	lat = (nfile.variables['latitude'])[:]
	lon = (nfile.variables['longitude'])[:]
	aps = (nfile.variables['aps'])[:]
	bps = (nfile.variables['bps'])[:]
	ps = (nfile.variables['ps'])[:];       ps = ps[itime,:,:]
	tsurf = (nfile.variables['tsurf'])[:]; tsurf = tsurf[itime,:,:]
	temp = (nfile.variables['temp'])[:];   temp = temp[itime,:,:,:]
	u = (nfile.variables['u'])[:];         u = u[itime,:,:,:]
	v = (nfile.variables['v'])[:];         v = v[itime,:,:,:]
	h2o = (nfile.variables['h2o_vap'])[:]; h2o = h2o[itime,:,:,:] * 28./18.
	ice = (nfile.variables['h2o_ice'])[:]; ice = ice[itime,:,:,:] + 1e-20
	nfile.close()

	# Create the 3D pressures
	sz = np.shape(temp)
	sz = np.flip(sz)
	press3D = np.zeros((sz), dtype=np.float32);
	for i in range(sz[0]):
		for j in range(sz[1]):
			press3D[i,j,:] = (aps + bps * ps[j,i]) * 1e-5
		#Endfor
	#Endfor

	# Save object parameters
	newf = []
	newf.append('<OBJECT>Exoplanet')
	newf.append('<OBJECT-NAME>Exoplanet')
	newf.append('<OBJECT-DIAMETER>12742')
	newf.append('<OBJECT-GRAVITY>9.8')
	newf.append('<OBJECT-GRAVITY-UNIT>g')
	newf.append('<OBJECT-STAR-DISTANCE>1.0')
	newf.append('<OBJECT-STAR-VELOCITY>0.0')
	newf.append('<OBJECT-SOLAR-LONGITUDE>0.0')
	newf.append('<OBJECT-SOLAR-LATITUDE>0.0')
	newf.append('<OBJECT-SEASON>0.0')
	newf.append('<OBJECT-STAR-TYPE>G')
	newf.append('<OBJECT-STAR-TEMPERATURE>5777')
	newf.append('<OBJECT-STAR-RADIUS>1.0')
	newf.append('<OBJECT-STAR-METALLICITY>0.0')
	newf.append('<OBJECT-OBS-LONGITUDE>0.0')
	newf.append('<OBJECT-OBS-LATITUDE>0.0')
	newf.append('<OBJECT-OBS-PERIOD>0.0')
	newf.append('<OBJECT-OBS-VELOCITY>0.0')

	# Save atmosphere parameters
	newf.append('<ATMOSPHERE-DESCRIPTION>Laboratoire de Meteorologie Dynamique (LMD) simulation')
	newf.append('<ATMOSPHERE-STRUCTURE>Equilibrium')
	newf.append('<ATMOSPHERE-PRESSURE>1.0')
	newf.append('<ATMOSPHERE-PUNIT>bar')
	newf.append('<ATMOSPHERE-TEMPERATURE>300')
	newf.append('<ATMOSPHERE-WEIGHT>28.0')
	newf.append('<ATMOSPHERE-LAYERS>0')
	newf.append('<ATMOSPHERE-NGAS>1')
	newf.append('<ATMOSPHERE-GAS>H2O')
	newf.append('<ATMOSPHERE-TYPE>HIT[1]')
	newf.append('<ATMOSPHERE-ABUN>1')
	newf.append('<ATMOSPHERE-UNIT>scl')
	newf.append('<ATMOSPHERE-NAERO>1')
	newf.append('<ATMOSPHERE-AEROS>WaterIce')
	newf.append('<ATMOSPHERE-ATYPE>Warren_ice_HRI')
	newf.append('<ATMOSPHERE-AABUN>1')
	newf.append('<ATMOSPHERE-AUNIT>scl')
	newf.append('<ATMOSPHERE-ASIZE>1')
	newf.append('<ATMOSPHERE-ASUNI>um')

	# Save surface parameters
	newf.append('<SURFACE-TEMPERATURE>300')
	newf.append('<SURFACE-ALBEDO>0.2')
	newf.append('<SURFACE-EMISSIVITY>0.8')
	newf.append('<SURFACE-NSURF>0')

	# Save GCM parameters
	vars = '<ATMOSPHERE-GCM-PARAMETERS>'
	vars = vars + str("%d,%d,%d,%.1f,%.1f,%.2f,%.2f" %(sz[0],sz[1],sz[2],lon[0],lat[0],lon[1]-lon[0],lat[1]-lat[0]))
	vars = vars + ',Winds,Tsurf,Psurf,Temperature,Pressure,H2O,WaterIce'
	newf.append(vars)

	with open(fileout,'w') as fw:
		for i in newf: fw.write(i+'\n')
	with open(fileout,'ab') as fb:
		if sys.version_info>=(3,0,0): bc=fb.write(bytes('<BINARY>',encoding = 'utf-8'))
		else: bc=fb.write('<BINARY>')
		fb.write(np.asarray(u,order='C',dtype=np.float32))
		fb.write(np.asarray(v,order='C',dtype=np.float32))
		fb.write(np.asarray(tsurf,order='C',dtype=np.float32))
		fb.write(np.log10(np.asarray(ps*1e-5,order='C',dtype=np.float32)))
		fb.write(np.asarray(temp,order='C',dtype=np.float32))
		fb.write(np.log10(np.asarray(np.transpose(press3D),order='C')))
		fb.write(np.log10(np.asarray(h2o,order='C',dtype=np.float32)))
		fb.write(np.log10(np.asarray(ice,order='C',dtype=np.float32)))
		if sys.version_info>=(3,0,0): bc=fb.write(bytes('</BINARY>',encoding = 'utf-8'))
		else: bc=fb.write('</BINARY>')
	fb.close()
#End convert

if __name__ == "__main__": convertgcm()
