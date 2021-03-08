# ---------------------------------------------------------------
# Script to convert netCDF climate files into PSG GCM files
# netCDF files in UM general circulation model
# Fauchez, Villanueva - NASA Goddard Space Flight Center
# March 2021
# ---------------------------------------------------------------
import sys
import numpy as np
from netCDF4 import Dataset as ncdf

# Module that performs the conversion
def convertgcm(filein = 'data/gcm_um.nc', fileout = 'gcm_psg.dat', itime=0):
	# Read the netCDF file
	nfile = ncdf(filein)
	lat = (nfile.variables['latitude_t'])[:];
	lon = (nfile.variables['longitude_t'])[:];
	alt = (nfile.variables['thlev_zsea_theta'])[:]
	tsurf = (nfile.variables['STASH_m01s00i024'])[:]; tsurf = tsurf[itime,:,:]
	temp = (nfile.variables['STASH_m01s16i004'])[:];  temp = temp[itime,:,:,:]
	press = (nfile.variables['STASH_m01s00i408'])[:]; press = press[itime,:,:,:]
	uk = (nfile.variables['STASH_m01s00i002'])[:];    uk = uk[itime,:,:,:]
	vk = (nfile.variables['STASH_m01s00i003'])[:];    vk = vk[itime,:,:,:]
	h2o = (nfile.variables['STASH_m01s00i010'])[:];   h2o = h2o[itime,:,:,:] * 28./18.
	ice = (nfile.variables['STASH_m01s00i012'])[:];   ice = ice[itime,:,:,:]
	liq = (nfile.variables['STASH_m01s00i254'])[:];   liq = liq[itime,:,:,:]
	rliq = (nfile.variables['STASH_m01s01i221'])[:];  rliq = rliq[itime,:,:,:]/1e6# Size of ice clouds in m
	nfile.close()

	# Fix variables
	sz = np.shape(temp)
	uw = np.zeros(sz); uw[:-1,:,:] = uk
	vw = np.zeros(sz); vw[:-1,:,:] = vk[:,:-1,:]
	temp = np.where((temp>0) & (np.isfinite(temp)), temp, 300.0)
	h2o = np.where((h2o>0) & (np.isfinite(h2o)), h2o, 1e-30)
	ice = np.where((ice>0) & (np.isfinite(ice)), ice, 1e-30)
	liq = np.where((liq>0) & (np.isfinite(liq)), liq, 1e-30)
	rliq = np.where((rliq>0) & (np.isfinite(rliq)), rliq, 1e-6)

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
	newf.append('<ATMOSPHERE-DESCRIPTION>Met Office Unified Model (UM) simulation')
	newf.append('<ATMOSPHERE-STRUCTURE>Equilibrium')
	newf.append('<ATMOSPHERE-PRESSURE>1.0')
	newf.append('<ATMOSPHERE-PUNIT>bar')
	newf.append('<ATMOSPHERE-WEIGHT>28.0')
	newf.append('<ATMOSPHERE-LAYERS>0')
	newf.append('<ATMOSPHERE-NGAS>3')
	newf.append('<ATMOSPHERE-GAS>N2,CO2,H2O')
	newf.append('<ATMOSPHERE-TYPE>HIT[22],HIT[2],HIT[1]')
	newf.append('<ATMOSPHERE-ABUN>99,400,1')
	newf.append('<ATMOSPHERE-UNIT>pct,ppm,scl')
	newf.append('<ATMOSPHERE-NMAX>2')
	newf.append('<ATMOSPHERE-LMAX>2')
	newf.append('<ATMOSPHERE-NAERO>2')
	newf.append('<ATMOSPHERE-AEROS>Water,WaterIce')
	newf.append('<ATMOSPHERE-ATYPE>AFCRL_Water_HRI,Warren_ice_HRI')
	newf.append('<ATMOSPHERE-AABUN>1,1')
	newf.append('<ATMOSPHERE-AUNIT>scl,scl')
	newf.append('<ATMOSPHERE-ASIZE>1,1')
	newf.append('<ATMOSPHERE-ASUNI>scl,um')

	# Save surface parameters
	newf.append('<SURFACE-TEMPERATURE>300')
	newf.append('<SURFACE-ALBEDO>0.2')
	newf.append('<SURFACE-EMISSIVITY>1.0')
	newf.append('<SURFACE-NSURF>0')

	# Save GCM parameters
	vars = '<ATMOSPHERE-GCM-PARAMETERS>'
	vars = vars + str("%d,%d,%d,%.1f,%.1f,%.2f,%.2f" %(sz[2],sz[1],sz[0],lon[0],lat[0],lon[1]-lon[0],lat[1]-lat[0]))
	vars = vars + ',Winds,Tsurf,Temperature,Pressure,H2O,Water,WaterIce,Water_size'
	newf.append(vars)
	with open(fileout,'w') as fw:
		for i in newf: fw.write(i+'\n')
	with open(fileout,'ab') as fb:
		if sys.version_info>=(3,0,0): bc=fb.write(bytes('<BINARY>',encoding = 'utf-8'))
		else: bc=fb.write('<BINARY>')
		fb.write(np.asarray(uw,order='C',dtype=np.float32))
		fb.write(np.asarray(vw,order='C',dtype=np.float32))
		fb.write(np.asarray(tsurf,order='C',dtype=np.float32))
		fb.write(np.asarray(temp,order='C',dtype=np.float32))
		fb.write(np.log10(np.asarray(press*1e-5,order='C',dtype=np.float32)))
		fb.write(np.log10(np.asarray(h2o,order='C',dtype=np.float32)))
		fb.write(np.log10(np.asarray(liq,order='C',dtype=np.float32)))
		fb.write(np.log10(np.asarray(ice,order='C',dtype=np.float32)))
		fb.write(np.log10(np.asarray(rliq,order='C',dtype=np.float32)))
		if sys.version_info>=(3,0,0): bc=fb.write(bytes('</BINARY>',encoding = 'utf-8'))
		else: bc=fb.write('</BINARY>')
	fb.close()
#End convert

if __name__ == "__main__": convertgcm()
