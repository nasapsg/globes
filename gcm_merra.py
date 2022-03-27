# ---------------------------------------------------------------
# Script to extract MERRA2 climatology and convert into PSG GCM files
# Villanueva - NASA Goddard Space Flight Center
# March 2022
# ---------------------------------------------------------------
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import h5py

def rebin(a, *args):
    shape = a.shape
    lenShape = len(shape)
    factor = (np.asarray(shape)/np.asarray(args)).astype(np.int32)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.sum(%d)'%(i+1) for i in range(lenShape)] + \
             ['/factor[%d]'%i for i in range(lenShape)]
    return eval(''.join(evList))
#Enddef

# Module that performs the extraction and conversion
def convertgcm(date = '2019/03/20 00:00', fileout = 'gcm_psg.dat', itime=0):

    # Process albedo data from image
    surf = np.zeros([7,91,144], dtype=np.float32)  # Surface parameters (0:Albedo, 1:Tsurf, 2:Ocean, 3:Snow, 4:Desert, 5:Forest, 6:Grass
    sz = surf.shape[1:]
    im = np.array(Image.open('data/earth.png'))
    dh = np.sum(im[:,:,0:3], axis=2)/(3.0*255.0)
    dt = np.flip(rebin(dh, 135, 270),0)
    sl = dt.shape
    for i in range(91):
        il = round(i*(sl[0]-1.0)/(sz[0]-1.0))
        for j in range(144):
            jl = round((j-0.5)*(sl[1]-1.0)/(sz[1]-1.0))
            if jl<0: jl=sl[1]+jl
            surf[0,i,j] = dt[il,jl]
        #Endfor
    #Endfor

    # Read MODIS (MCD12C1) land cover maps
    hdf5f = h5py.File('data/modis.h5', 'r')
    landh = np.array(list(hdf5f['MOD12C1/Data Fields/Land_Cover_Type_1_Percent']))
    hdf5f.close()
    lands = landh*0.0
    lands[:,:-20,:] = landh[:,20:,:]
    lands[:,20:,:] = landh[:,0:-20,:]
    land = np.flip(rebin(lands, 90, 144, 17),0)
    ocean = (land[:,:,0] + land[:,:,11])/100.0
    snow = land[:,:,15]/100.0
    desert = (land[:,:,16] + land[:,:,13] + land[:,:,7])/100.0
    forest = (land[:,:,1]+land[:,:,2]+land[:,:,3]+land[:,:,4]+land[:,:,5]+land[:,:,6]+land[:,:,8])/100.0
    grass = (land[:,:,9] + land[:,:,10] + land[:,:,12] + land[:,:,14])/100.0
    for i in range(91):
        il = round(i*89.0/90.0)
        surf[2,i,:] = ocean[il,:]
        surf[3,i,:] = snow[il,:]
        surf[4,i,:] = desert[il,:]
        surf[5,i,:] = forest[il,:]
        surf[6,i,:] = grass[il,:]
    #Endfor

    # Process MERRA file
    if not os.path.exists('gcm_merra.txt'):
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
        pwd = 'pwd'

        # Parse and run request
        cmd = 'curl -s -u %s:%s -b gcm_merra.ck -c gcm_merra.ck -m 12000 -L --location-trusted' % (user,pwd)
        cmd = '%s -g https://goldsmr5.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2I3NVASM.5.12.4/' % cmd
        cmd = '%s%s/%s/MERRA2_%s.inst3_3d_asm_Nv.%s%s%s.nc4.ascii?' % (cmd,sd[0],sd[1],mvr,sd[0],sd[1],sd[2])
        cmd = '%sPL[%d:%d][0:71][0:360][0:575],T[%d:%d][0:71][0:360][0:575],' % (cmd,itme,itme,itme,itme)
        cmd = '%sU[%d:%d][0:71][0:360][0:575],V[%d:%d][0:71][0:360][0:575],O3[%d:%d][0:71][0:360][0:575],' % (cmd,itme,itme,itme,itme,itme,itme)
        cmd = '%sQV[%d:%d][0:71][0:360][0:575],QI[%d:%d][0:71][0:360][0:575],QL[%d:%d][0:71][0:360][0:575]' % (cmd,itme,itme,itme,itme,itme,itme)
        cmd = '%s > gcm_merra.txt' % cmd
        os.system(cmd)
    #Endif

    # Read file (bin by 4 in lat/lon to reduce the size by 16)
    fp=open('gcm_merra.txt','r'); line=fp.readline();
    if line[0:15]!="Dataset: MERRA2": exit()
    data = np.zeros([8,72,91,144], dtype=np.float32)
    npts = np.zeros([8,72,91,144], dtype=np.float32)
    while True:
        line = fp.readline()
        if not line: break
        st = line.split(',')
        vt = st[0].split('[')

        if len(vt)==1: continue # Disregard lat/lon/lev arrays
        itime = int(vt[1][0:-1])          # The curl command is for a specific time, so one time is returned only (always 0)
        ialt  = 71 - int(vt[2][0:-1])     # Altitudes are reversed in MERRA (0:top, 71:bottom)
        ilat  = round(int(vt[3][0:-1])/4) # Bin latitude by 4
        if ilat>=91: ilat=90

        if   vt[0]=='U':  ivar=0; kscl=1.0            # Horizontal u-wind [m/s]
        elif vt[0]=='V':  ivar=1; kscl=1.0            # Horizontal v-wind [m/s]
        elif vt[0]=='T':  ivar=2; kscl=1.0            # Temperature [K]
        elif vt[0]=='PL': ivar=3; kscl=1e-5           # Pressure [bar] (converted from Pa)
        elif vt[0]=='O3': ivar=4; kscl=28.97/48.0     # O3 ozone [mol/mol] (converted from kg/kg)
        elif vt[0]=='QV': ivar=5; kscl=28.97/18.0     # H2O water vapor [mol/mol] (converted from kg/kg)
        elif vt[0]=='QI': ivar=6; kscl=1.0            # H2O cloud water ice [kg/kg]
        elif vt[0]=='QL': ivar=7; kscl=1.0            # H2O cloud liquie water [kg/kg]

        for i in range(0,576):
            ilon = round(i/4) # Bin longitude by 4
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
    fp.close()
    data = data/npts
    surf[1,:,:] = data[2,0,:,:]

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
    newf.append('<ATMOSPHERE-NGAS>8')
    newf.append('<ATMOSPHERE-GAS>N2,O2,CO2,CH4,N2O,CO,H2O,O3')
    newf.append('<ATMOSPHERE-TYPE>HIT[22],HIT[7],HIT[2],HIT[6],HIT[4],HIT[5],HIT[1],HIT[3]')
    newf.append('<ATMOSPHERE-ABUN>78,20.9,400,1.6,0.3,0.1,1,1')
    newf.append('<ATMOSPHERE-UNIT>pct,pct,ppm,ppm,ppm,ppm,scl,scl')
    newf.append('<ATMOSPHERE-NAERO>2')
    newf.append('<ATMOSPHERE-NMAX>1')
    newf.append('<ATMOSPHERE-LMAX>2')
    newf.append('<ATMOSPHERE-AEROS>Water,WaterIce')
    newf.append('<ATMOSPHERE-ATYPE>AFCRL_Water_HRI,Warren_ice_HRI')
    newf.append('<ATMOSPHERE-AABUN>1,1')
    newf.append('<ATMOSPHERE-AUNIT>scl,scl')
    newf.append('<ATMOSPHERE-ASIZE>0.1,0.1')
    newf.append('<ATMOSPHERE-ASUNI>um,um')
    newf.append('<SURFACE-TEMPERATURE>290.0')
    newf.append('<SURFACE-ALBEDO>0.30')
    newf.append('<SURFACE-EMISSIVITY>0.70')
    newf.append('<SURFACE-MODEL>Lambert')
    newf.append('<SURFACE-NSURF>5')
    newf.append('<SURFACE-SURF>Ocean,Snow,Grass,Desert,Forest')
    newf.append('<SURFACE-TYPE>McLinden_GSFC,USGS_GSFC,Dry+green_series_30964_USGS,JB-JLB-851-A-BKR1JB851A_RELAB,USGS_GSFC')
    newf.append('<SURFACE-ABUN>0.5,0.1,0.1,0.2,0.1')
    newf.append('<SURFACE-UNIT>scl,scl,scl,scl,scl')

    # Simulation parameters
    newf.append('<GENERATOR-RANGE1>0.1')
    newf.append('<GENERATOR-RANGE2>2.5')
    newf.append('<GENERATOR-RANGEUNIT>um')
    newf.append('<GENERATOR-RESOLUTION>100')
    newf.append('<GENERATOR-RESOLUTIONUNIT>RP')
    newf.append('<GENERATOR-GAS-MODEL>Y')
    newf.append('<GENERATOR-CONT-MODEL>Y')
    newf.append('<GENERATOR-CONT-STELLAR>N')
    newf.append('<GENERATOR-TRANS-APPLY>N')
    newf.append('<GENERATOR-TRANS-SHOW>N')
    newf.append('<GENERATOR-RADUNITS>rif')
    newf.append('<GENERATOR-LOGRAD>N')
    newf.append('<GENERATOR-TELESCOPE>SINGLE')
    newf.append('<GENERATOR-BEAM>1')
    newf.append('<GENERATOR-BEAM-UNIT>diameter')
    newf.append('<GENERATOR-NOISE>NO')
    newf.append('<GEOMETRY>Observatory')
    newf.append('<GEOMETRY-OBS-ALTITUDE>1.00')
    newf.append('<GEOMETRY-ALTITUDE-UNIT>AU')

    # Save configuration file
    newf.append('<ATMOSPHERE-GCM-PARAMETERS>144,91,72,-180,-90,2.5,2.0,Surf_Albedo,Surf_T,Surf_Ocean,Surf_Snow,Surf_Desert,Surf_Forest,Surf_Grass,Winds,Temperature,Pressure,O3,H2O,WaterIce,Water')
    with open(fileout,'w') as fw:
        for i in newf: fw.write(i+'\n')
    with open(fileout,'ab') as fb:
        if sys.version_info>=(3,0,0): bc=fb.write(bytes('<BINARY>',encoding = 'utf-8'))
        else: bc=fb.write('<BINARY>')
        fb.write(np.asarray(surf,order='C'))
        fb.write(np.asarray(data,order='C'))
        if sys.version_info>=(3,0,0): bc=fb.write(bytes('</BINARY>',encoding = 'utf-8'))
        else: bc=fb.write('</BINARY>')
    fb.close()
#End convert

if __name__ == "__main__": convertgcm()
