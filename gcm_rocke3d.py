# ---------------------------------------------------------------
# Script to convert netCDF climate files into PSG GCM files
# netCDF files in ROCKE3D climate model
# Villanueva, Way, Fauchez - NASA Goddard Space Flight Center
# This modules requires the three ROCKE3D files: aijk, aijl, aij
# February 2021
# ---------------------------------------------------------------
import sys
import numpy as np
from netCDF4 import Dataset as ncdf

# Module that performs the conversion
def convertgcm(filein = 'data/gcm_rocke3d', fileout = 'gcm_psg.dat', itime=0):
	# Read ROCKE-3D outputs across three different netCDF files
	nfile = ncdf("%s_aijk.nc" % filein)
	uk = np.array((nfile.variables['ub'])[:]);
	vk = np.array((nfile.variables['vb'])[:]);
	dpb = np.array((nfile.variables['dpb'])[:]);
	nfile.close()

	nfile = ncdf("%s_aijl.nc" % filein)
	lat = (nfile.variables['lat'])[:]
	lon = (nfile.variables['lon'])[:]; lon = lon + 180.
	plm = np.array((nfile.variables['plm'])[:]);
	tmp = np.array((nfile.variables['TempL'])[:]); #temperature in K (in level)
	h2o = np.array((nfile.variables['SpHuL'])[:]); #specific humidity in kg/kg
	h2o = h2o/(1.0-h2o); h2o = h2o*(28.0/18.0) + 1E-30; #Water vapor abundance [molecules/molecule]
	cldice = np.array((nfile.variables['icecld'])[:]) + 1E-30; #Water ice clouds in kg/kg
	cldliq = np.array((nfile.variables['wtrcld'])[:]) + 1E-30; #Liquid water clouds in kg/kg
	nfile.close()

	nfile = ncdf("%s_aij.nc" % filein)
	lat = (nfile.variables['lat'])[:]
	lon = (nfile.variables['lon'])[:]; lon = lon + 180.
	albedo = np.array((nfile.variables['grnd_alb'])[:])*1e-2; #ground albedo [0 to 1.0]
	albedo = np.where((~np.isfinite(albedo)) | (albedo<0.) | (albedo>1.0), 0.3, albedo);
	tsurf = np.array((nfile.variables['tgrnd'])[:]) + 273.15; #ground temperature [K]
	nfile.close()

	# Transform variables into a 3D vectors
	sz = np.shape(tmp)
	press3D = np.zeros((sz), dtype=np.float32);
	ul = uk; vl = vk
	for i in range(sz[1]):
		for j in range(sz[2]):
			press3D[:,i,j] = (plm + dpb[:,i,j])*1e-3
			luk = 0.0; lvk = 0.0
			for k in range(sz[0]-1,-1,-1):
				if np.isfinite(uk[k,i,j]) and np.abs(uk[k,i,j])<1e3: luk = uk[k,i,j]
				if np.isfinite(vk[k,i,j]) and np.abs(vk[k,i,j])<1e3: lvk = vk[k,i,j]
				ul[k,i,j] = luk
				vl[k,i,j] = lvk
			#Endfor layers
		#Endfor
	#Endfor

	# Save object parameters
	newf = []
	newf.append('<OBJECT>Exoplanet')
	newf.append('<OBJECT-NAME>Exoplanet')
	newf.append('<OBJECT-DATE>2018/01/01 10:00')
	newf.append('<OBJECT-DIAMETER>12742')
	newf.append('<OBJECT-GRAVITY>9.8')
	newf.append('<OBJECT-GRAVITY-UNIT>g')
	newf.append('<OBJECT-STAR-DISTANCE>1.0')
	newf.append('<OBJECT-STAR-VELOCITY>0.0')
	newf.append('<OBJECT-SOLAR-LONGITUDE>45.0')
	newf.append('<OBJECT-SOLAR-LATITUDE>0.0')
	newf.append('<OBJECT-SEASON>45.0')
	newf.append('<OBJECT-STAR-TYPE>G')
	newf.append('<OBJECT-STAR-TEMPERATURE>5777')
	newf.append('<OBJECT-STAR-RADIUS>1.0')
	newf.append('<OBJECT-STAR-METALLICITY>0.0')
	newf.append('<OBJECT-OBS-LONGITUDE>0.0')
	newf.append('<OBJECT-OBS-LATITUDE>0.0')
	newf.append('<OBJECT-OBS-PERIOD>0.0')
	newf.append('<OBJECT-OBS-VELOCITY>0.0')

	# Save atmosphere parameters
	newf.append('<ATMOSPHERE-DESCRIPTION>ROCKE-3D Simulation')
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
	newf.append('<ATMOSPHERE-ASUNI>um,um')

	# Save surface parameters
	newf.append('<SURFACE-TEMPERATURE>300')
	newf.append('<SURFACE-ALBEDO>0.2')
	newf.append('<SURFACE-EMISSIVITY>0.8')
	newf.append('<SURFACE-NSURF>0')

	# Instrument parameters (defaults)
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
	newf.append('<GENERATOR-GCM-BINNING>200')

	# Save GCM parameters
	vars = '<ATMOSPHERE-GCM-PARAMETERS>'
	vars = vars + str("%d,%d,%d,%.1f,%.1f,%.2f,%.2f" %(sz[2],sz[1],sz[0],lon[0],lat[0],lon[1]-lon[0],lat[1]-lat[0]))
	vars = vars + ',Winds,Temperature,Tsurf,Albedo,Pressure,H2O,Water,WaterIce'
	newf.append(vars)

	with open(fileout,'w') as fw:
		for i in newf: fw.write(i+'\n')
	with open(fileout,'ab') as fb:
		if sys.version_info>=(3,0,0): bc=fb.write(bytes('<BINARY>',encoding = 'utf-8'))
		else: bc=fb.write('<BINARY>')
		fb.write(np.asarray(ul,order='C'))
		fb.write(np.asarray(vl,order='C'))
		fb.write(np.asarray(tmp,order='C'))
		fb.write(np.asarray(tsurf,order='C'))
		fb.write(np.asarray(albedo,order='C'))
		fb.write(np.log10(np.asarray(press3D,order='C')))
		fb.write(np.log10(np.asarray(h2o,order='C')))
		fb.write(np.log10(np.asarray(cldliq,order='C')))
		fb.write(np.log10(np.asarray(cldice,order='C')))
		if sys.version_info>=(3,0,0): bc=fb.write(bytes('</BINARY>',encoding = 'utf-8'))
		else: bc=fb.write('</BINARY>')
	fb.close()
#End convert

if __name__ == "__main__": convertgcm()
