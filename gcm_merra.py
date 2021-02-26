# ---------------------------------------------------------------
# Script to extract MERRA2 climatology and convert into PSG GCM files
# Villanueva - NASA Goddard Space Flight Center
# February 2021
# ---------------------------------------------------------------
import sys, os
import numpy as np
from PIL import Image

# Module that performs the extraction and conversion
def convertgcm(date = '2018/01/01 10:00', fileout = 'gcm_psg.dat', itime=0):

	# Process date/time UT (YYYY/MM/DD HH:MN)
	ss = date.split()
	sd = ss[0].split('/')
	if int(sd[0])>2010: mvr="400"
	elif int(sd[0])>2000: mvr="300"
	elif int(sd[0])>1991: mvr="200"
	else: mvr="100"
	st = ss[1].split(':')
	itme = round((float(st[0])*60 + float(st[1]))/180)
	if itme>=8: itme=7
	if itme<0: itme=0

	# Define EarthData Login info (https://urs.earthdata.nasa.gov/)
	user = 'user'
	pwd = 'password'

	# Parse and run request
	cmd = 'curl -s -u %s:%s -b gcm_merra.ck -c gcm_merra.ck -m 12000 -L --location-trusted' % (user,pwd)
	cmd = '%s -g https://goldsmr5.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2I3NVASM.5.12.4/' % cmd
	cmd = '%s%s/%s/MERRA2_%s.inst3_3d_asm_Nv.%s%s%s.nc4.ascii?' % (cmd,sd[0],sd[1],mvr,sd[0],sd[1],sd[2])
	cmd = '%sPL[%d:%d][0:71][0:360][0:575],T[%d:%d][0:71][0:360][0:575],' % (cmd,itme,itme,itme,itme)
	cmd = '%sU[%d:%d][0:71][0:360][0:575],V[%d:%d][0:71][0:360][0:575],O3[%d:%d][0:71][0:360][0:575],' % (cmd,itme,itme,itme,itme,itme,itme)
	cmd = '%sQV[%d:%d][0:71][0:360][0:575],QI[%d:%d][0:71][0:360][0:575],QL[%d:%d][0:71][0:360][0:575]' % (cmd,itme,itme,itme,itme,itme,itme)
	cmd = '%s > gcm_merra.txt' % cmd
	os.system(cmd)

	# Read file (bin by 4 in lat/lon to reduce the size by 16)
	fp=open('gcm_merra.txt','r'); ll=fp.readline(); lines=fp.readlines(); fp.close()
	if ll[0:15]!="Dataset: MERRA2": exit()
	data = np.zeros([8,72,91,144], dtype=np.float32)
	npts = np.zeros([8,72,91,144], dtype=np.float32)
	for line in lines:
		st = line.split(',')
		vt = st[0].split('.')
		if vt[1]=='lon': continue
		i1 = st[0].find('lev=')+4; i2 = st[0][i1:].find(']')+i1; ialt=71-int(st[0][i1:i2])-1
		i1 = st[0].find('lat=')+4; i2 = st[0][i1:].find(']')+i1; ilat=round((float(st[0][i1:i2])+90.0)/0.5)
		ilat = round(0.25*ilat)                       # Bin the data by 4 to reduce file/memory size
		if vt[0]=='U': ivar=0; kscl=1.0               # Horizontal u-wind [m/s]
		elif vt[0]=='V': ivar=1; kscl=1.0             # Horizontal v-wind [m/s]
		elif vt[0]=='T': ivar=2; kscl=1.0             # Temperature [K]
		elif vt[0]=='PL': ivar=3; kscl=1e-5           # Pressure [bar] (converted from Pa)
		elif vt[0]=='O3': ivar=4; kscl=28.97/48.0     # O3 ozone [mol/mol] (converted from kg/kg)
		elif vt[0]=='QV': ivar=5; kscl=28.97/18.0     # H2O water vapor [mol/mol] (converted from kg/kg)
		elif vt[0]=='QI': ivar=6; kscl=1.0            # H2O cloud water ice [kg/kg]
		elif vt[0]=='QL': ivar=7; kscl=1.0            # H2O cloud liquie water [kg/kg]
		for i in range(0,576):
			ilon = round(0.25*i)
			if ilon>=144: ilon=ilon-144
			val = float(st[i+1])*kscl
			if ivar>2:
				if val<1e-30: val=1e-30
				val=np.log10(val)
			#Endif
			data[ivar,ialt,ilat,ilon] = data[ivar,ialt,ilat,ilon] + val
			npts[ivar,ialt,ilat,ilon] = npts[ivar,ialt,ilat,ilon] + 1
		#Endfor
	#Endfor
	data = data/npts

	# Process albedo data from image
	im = np.array(Image.open('data/earth.png'))
	imsz = im.shape
	tflux = np.sum(im[:,:,0:3], axis=2)/(3.0*255.0)
	albedo = np.zeros([91,144], dtype=np.float32)
	npts = np.zeros([91,144], dtype=np.float32)
	for i in range(imsz[0]):
		il = 90 - round(i*90.0/imsz[0])
		for j in range(imsz[1]):
			jl = round(j*143.0/imsz[1])
			albedo[il,jl] = albedo[il,jl] + tflux[i,j]
			npts[il,jl] = npts[il,jl] + 1.0
		#Endfor
	#Endfor
	albedo = albedo/npts

	# Write file
	newf = []
	newf.append('<OBJECT>Earth')
	newf.append('<OBJECT-NAME>Earth')
	newf.append('<OBJECT-DATE>%s' % date)
	newf.append('<OBJECT-DIAMETER>12742')
	newf.append('<OBJECT-GRAVITY>9.807')
	newf.append('<OBJECT-GRAVITY-UNIT>g')
	newf.append('<OBJECT-STAR-DISTANCE>1.0')
	newf.append('<OBJECT-STAR-VELOCITY>0.0')
	newf.append('<OBJECT-SOLAR-LONGITUDE>0.0')
	newf.append('<OBJECT-SOLAR-LATITUDE>0.0')
	newf.append('<OBJECT-SEASON>0.0')
	newf.append('<OBJECT-STAR-TYPE>G')
	newf.append('<OBJECT-STAR-TEMPERATURE>5777')
	newf.append('<OBJECT-STAR-RADIUS>1.0')
	newf.append('<OBJECT-OBS-LONGITUDE>0.0')
	newf.append('<OBJECT-OBS-LATITUDE>0.0')
	newf.append('<OBJECT-OBS-VELOCITY>0.0')
	newf.append('<OBJECT-PERIOD>1.0')
	newf.append('<OBJECT-STAR-METALLICITY>0.0')
	newf.append('<OBJECT-INCLINATION>90.0')
	newf.append('<OBJECT-ECCENTRICITY>0.0')
	newf.append('<OBJECT-PERIAPSIS>0.00')

	# Save atmosphere parameters
	newf.append('<ATMOSPHERE-DESCRIPTION>NASA/MERRA2 3D data')
	newf.append('<ATMOSPHERE-STRUCTURE>Equilibrium')
	newf.append('<ATMOSPHERE-WEIGHT>28.97')
	newf.append('<ATMOSPHERE-PRESSURE>1')
	newf.append('<ATMOSPHERE-PUNIT>bar')
	newf.append('<ATMOSPHERE-TEMPERATURE>270')
	newf.append('<ATMOSPHERE-NGAS>8')
	newf.append('<ATMOSPHERE-GAS>N2,O2,CO2,CH4,N2O,CO,H2O,O3')
	newf.append('<ATMOSPHERE-TYPE>HIT[22],HIT[7],HIT[2],HIT[6],HIT[4],HIT[5],HIT[1],HIT[3]')
	newf.append('<ATMOSPHERE-ABUN>78,20.9,400,1.6,0.3,0.1,1,1')
	newf.append('<ATMOSPHERE-UNIT>pct,pct,ppm,ppm,ppm,ppm,scl,scl')
	newf.append('<ATMOSPHERE-NAERO>2')
	newf.append('<ATMOSPHERE-NMAX>0')
	newf.append('<ATMOSPHERE-LMAX>0')
	newf.append('<ATMOSPHERE-AEROS>Water,WaterIce')
	newf.append('<ATMOSPHERE-ATYPE>AFCRL_Water_HRI,Warren_ice_HRI')
	newf.append('<ATMOSPHERE-AABUN>1,1')
	newf.append('<ATMOSPHERE-AUNIT>scl,scl')
	newf.append('<ATMOSPHERE-ASIZE>1,1')
	newf.append('<ATMOSPHERE-ASUNI>um,um')

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
	newf.append('<ATMOSPHERE-GCM-PARAMETERS>144,91,72,-180,-90,2.5,2.0,Albedo,Winds,Temperature,Pressure,O3,H2O,WaterIce,Water')
	with open(fileout,'w') as fw:
		for i in newf: fw.write(i+'\n')
	with open(fileout,'ab') as fb:
		if sys.version_info>=(3,0,0): bc=fb.write(bytes('<BINARY>',encoding = 'utf-8'))
		else: bc=fb.write('<BINARY>')
		fb.write(np.asarray(albedo,order='C'))
		fb.write(np.asarray(data,order='C'))
		if sys.version_info>=(3,0,0): bc=fb.write(bytes('</BINARY>',encoding = 'utf-8'))
		else: bc=fb.write('</BINARY>')
	fb.close()
#End convert

if __name__ == "__main__": convertgcm()
