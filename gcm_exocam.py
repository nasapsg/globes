# ---------------------------------------------------------------
# Script to convert netCDF climate files into PSG GCM files
# netCDF files in ExoCAM (Community Atmospheric Model) format
# Villanueva, Fauchez, Suissa, Wolf - NASA Goddard Space Flight Center
# February 2021
# ---------------------------------------------------------------
import sys
import numpy as np
from netCDF4 import Dataset as ncdf

# Module that performs the conversion
def convertgcm(filein = 'data/gcm_exocam.nc', fileout = 'gcm_psg.dat', itime=0):
	# Read values from netCDF
	nfile = ncdf(filein)
	lat = (nfile.variables['lat'])[:]
	lon = (nfile.variables['lon'])[:];       lon = lon+180.
	hyam = (nfile.variables['hyam'])[:];     hyam = np.flipud(hyam)
	hybm = (nfile.variables['hybm'])[:];     hybm = np.flipud(hybm)
	P0 = (nfile.variables['P0'])[:];         P0 = P0*1e-5
	PS = (nfile.variables['PS'])[:];         PS = PS[itime,:,:]*1e-5
	U = (nfile.variables['U'])[:];           U = U[itime,::-1,:,:]
	V = (nfile.variables['V'])[:];           V = V[itime,::-1,:,:]
	T = (nfile.variables['T'])[:];           T = T[itime,::-1,:,:]
	TS = (nfile.variables['TS'])[:];         TS = TS[itime,:,:]
	#FUS = (nfile.variables['FUS'])[:];   	 FUS = FUS[itime,0,:,:]
	#FDS = (nfile.variables['FDS'])[:];   	 FDS = FDS[itime,0,:,:]; ASDIR = FUS/FDS
	ASDIR = (nfile.variables['ASDIR'])[:];   ASDIR = ASDIR[itime,:,:]
	CLDICE = (nfile.variables['CLDICE'])[:]; CLDICE = CLDICE[itime,::-1,:,:]
	CLDLIQ = (nfile.variables['CLDLIQ'])[:]; CLDLIQ = CLDLIQ[itime,::-1,:,:]
	REI = (nfile.variables['REI'])[:];       REI = REI[itime,::-1,:,:]/1e6
	REL = (nfile.variables['REL'])[:];       REL = REL[itime,::-1,:,:]/1e6
	CH4 = (nfile.variables['ch4vmr'])[:];    CH4 = float(CH4[itime]) + T*0.0
	CO2 = (nfile.variables['co2vmr'])[:];    CO2 = float(CO2[itime]) + T*0.0
	N2O = (nfile.variables['n2ovmr'])[:];    N2O = float(N2O[itime]) + T*0.0
	H2O = (nfile.variables['Q'])[:];         H2O = H2O[itime,::-1,:,:]; H2O = H2O/(1.0-H2O); H2O = H2O*(28.0/18.0)
	nfile.close()

	# Fix variables
	ASDIR = np.where((ASDIR>=0) & (ASDIR<=1.0) & (np.isfinite(ASDIR)), ASDIR, 0.3)
	CLDICE = np.where((CLDICE>0) & (np.isfinite(CLDICE)), CLDICE, 1e-30)
	CLDLIQ = np.where((CLDLIQ>0) & (np.isfinite(CLDLIQ)), CLDLIQ, 1e-30)
	REI = np.where((REI>0) & (np.isfinite(REI)), REI, 1e-6)
	REL = np.where((REL>0) & (np.isfinite(REL)), REL, 1e-6)
	CH4 = np.where((CH4>0) & (np.isfinite(CH4)), CH4, 1e-30)
	CO2 = np.where((CO2>0) & (np.isfinite(CO2)), CO2, 1e-30)
	N2O = np.where((N2O>0) & (np.isfinite(N2O)), N2O, 1e-30)
	H2O = np.where((H2O>0) & (np.isfinite(H2O)), H2O, 1e-30)

	# Compute profiles of pressure and N2 abundances
	sz = np.shape(T)
	sz = np.flip(sz)
	press3D = np.zeros((sz), dtype=np.float32);
	for i in range(sz[0]):
		for j in range(sz[1]):
			press3D[i,j,:] = hyam * P0 + hybm * PS[j,i]
	N2 = 1.0 - (CH4 + CO2 + N2O + H2O)

	# Save object parameters
	# Example with TOI-700d parameters
	newf = []
	newf.append('<OBJECT>Exoplanet')
	newf.append('<OBJECT-NAME>Exoplanet')
	newf.append('<OBJECT-DIAMETER>14207.330')
	newf.append('<OBJECT-GRAVITY>11.2396')
	newf.append('<OBJECT-GRAVITY-UNIT>g')
	newf.append('<OBJECT-STAR-DISTANCE>0.163356')
	newf.append('<OBJECT-STAR-VELOCITY>0.0')
	newf.append('<OBJECT-SOLAR-LONGITUDE>0.0')
	newf.append('<OBJECT-SOLAR-LATITUDE>0.0')
	newf.append('<OBJECT-SEASON>0.0')
	newf.append('<OBJECT-STAR-TYPE>M')
	newf.append('<OBJECT-STAR-TEMPERATURE>3480.29')
	newf.append('<OBJECT-STAR-RADIUS>0.4185')
	newf.append('<OBJECT-STAR-METALLICITY>0.0')
	newf.append('<OBJECT-OBS-LONGITUDE>0.0')
	newf.append('<OBJECT-OBS-LATITUDE>0.0')
	newf.append('<OBJECT-OBS-PERIOD>37.42600445188339')
	newf.append('<GEOMETRY>Observatory')
	newf.append('<GEOMETRY-OBS-ALTITUDE>14.6')
	newf.append('<GEOMETRY-ALTITUDE-UNIT>pc')

	# Save atmosphere parameters
	newf.append('<ATMOSPHERE-DESCRIPTION>Community Atmospheric Model (CAM) simulation')
	newf.append('<ATMOSPHERE-STRUCTURE>Equilibrium')
	newf.append('<ATMOSPHERE-PRESSURE>1.0')
	newf.append('<ATMOSPHERE-PUNIT>bar')
	newf.append('<ATMOSPHERE-WEIGHT>28.0')
	newf.append('<ATMOSPHERE-LAYERS>0')
	newf.append('<ATMOSPHERE-NGAS>3')
	newf.append('<ATMOSPHERE-GAS>N2,CO2,H2O')
	newf.append('<ATMOSPHERE-TYPE>HIT[22],HIT[1],HIT[2]')
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
	newf.append('<ATMOSPHERE-ASUNI>scl,scl')#

	# Save surface parameters
	newf.append('<SURFACE-TEMPERATURE>300')
	newf.append('<SURFACE-ALBEDO>0.2')
	newf.append('<SURFACE-EMISSIVITY>1.0')
	newf.append('<SURFACE-NSURF>0')

	# Simulation parameters
	newf.append('<GENERATOR-RANGE1>0.4')
	newf.append('<GENERATOR-RANGE2>30')
	newf.append('<GENERATOR-RANGEUNIT>um')
	newf.append('<GENERATOR-RESOLUTION>500')
	newf.append('<GENERATOR-RESOLUTIONUNIT>RP')
	newf.append('<GENERATOR-GAS-MODEL>Y')
	newf.append('<GENERATOR-CONT-MODEL>Y')
	newf.append('<GENERATOR-CONT-STELLAR>N')
	newf.append('<GENERATOR-TRANS-APPLY>N')
	newf.append('<GENERATOR-TRANS-SHOW>N')
	newf.append('<GENERATOR-RADUNITS>ppm')
	newf.append('<GENERATOR-LOGRAD>Y')
	newf.append('<GENERATOR-TELESCOPE>SINGLE')
	newf.append('<GENERATOR-BEAM>1')
	newf.append('<GENERATOR-BEAM-UNIT>arcsec')
	newf.append('<GENERATOR-NOISE>NO')

	# Save configuration file
	vars = '<ATMOSPHERE-GCM-PARAMETERS>'
	vars = vars + str("%d,%d,%d,%.1f,%.1f,%.2f,%.2f" %(sz[0],sz[1],sz[2],lon[0],lat[0],lon[1]-lon[0],lat[1]-lat[0]))
	vars = vars + ',Winds,Tsurf,Psurf,Albedo,Temperature,Pressure,H2O,Water,WaterIce,Water_size,WaterIce_size'#
	newf.append(vars)
	with open(fileout,'w') as fw:
		for i in newf: fw.write(i+'\n')
	with open(fileout,'ab') as fb:
		if sys.version_info>=(3,0,0): bc=fb.write(bytes('<BINARY>',encoding = 'utf-8'))
		else: bc=fb.write('<BINARY>')
		fb.write(np.asarray(U,order='C'))
		fb.write(np.asarray(V,order='C'))
		fb.write(np.asarray(TS,order='C'))
		fb.write(np.asarray(PS,order='C'))
		fb.write(np.asarray(ASDIR,order='C'))
		fb.write(np.asarray(T,order='C'))
		fb.write(np.log10(np.asarray(np.transpose(press3D),order='C')))
		fb.write(np.log10(np.asarray(H2O,order='C')))
		fb.write(np.log10(np.asarray(CLDLIQ,order='C')))
		fb.write(np.log10(np.asarray(CLDICE,order='C')))
		fb.write(np.log10(np.asarray(REL,order='C')))
		fb.write(np.log10(np.asarray(REI,order='C')))
		if sys.version_info>=(3,0,0): bc=fb.write(bytes('</BINARY>',encoding = 'utf-8'))
		else: bc=fb.write('</BINARY>')
	fb.close()

	print(U.shape, ASDIR.shape, press3D.shape)
#End convert

if __name__ == "__main__": convertgcm()
