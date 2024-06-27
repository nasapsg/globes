# -------------------------------------------------------
# PyGlobES - A python interface for PSG-GlobES
# By Kofman & Villanueva  (vincent.kofman@nasa.gov)
# More details in Kofman et al. 2024, Planetary Science Journal
# This code was written for the Kofman+24 work, and is tailored specifically for this work. Feel free to use any aspects and modify as needed.
# if you use any of these codes, please consider citing our work (Kofman+24)
# -------------------------------------------------------
import os
import struct
import requests
import datetime
import numpy as np

# This is the list of phyiscal constants used in PSG (C code from Geronimo)
#define CS        299792458.0  // Speed of light [m/s]
#define HP    6.626070040e-34  // Planck's constant [W s2] or [J s]
#define KB     1.38064852e-23  // Boltzmann' constant [J/K]
#define C2         1.43877736  // Second radiation constant [K/cm-1]
#define STEFBOL   5.670367e-8  // Stefan–Boltzmann constant [W m-2 K-4]
#define AU     149597870700.0  // Astronomical unit [m]
#define PC       3.0857215e16  // Parsec [m]
#define PI      3.14159265359  // PI
#define PK      57.2957795131  // Degree to radians
#define SUNRADIUS   6.95700e8  // Radius of the Sun [m]
#define SUNMASS     1.9884e30  // Sun mass [Kg]
#define GCONST    6.67430e-11  // Gravitational constant [m3 kg-1 s-2]
#define ASR     4.2545166e+10  // Arcseconds2 / sr ((3600.0*180/!PI)^2.0)
#define AVOG   6.022140857e23  // Avogadro number [molecules/mol]
#define ASCL    4.8481368e-06  // Radians per arcsec
#define PTOA             1e-8  // Pressure where we assume Top of the atmosphere [mbar]
#define KDOPS  9.25108423e-14  // Doppler parameter [g/mol/K] - 1e3*(KB*AVOG)/(CS*CS)

psg_link = 'http://localhost:3000'

# constants used by Kofman
h   = 6.62607015e-34     # J s
Na  = 6.02214076e23      # mol-1
kb  = 1.380649e-23       # J/K
c   = 299792458.0        # m/s
amu = 1.660539040e-27    # kg
G   = 6.67430e-11        # m3⋅kg−1⋅s−2
c2  = h*c/kb

def star_sr_pl(Rstar,d_pl):
    d_pl *=  1.496e8 # km/AU
    """
    Returns the size of the star in the sky in steradians at the distance of the planet from the star
    input: rstar in kilometers, d_pl in AU
    Output: single value in sr. 
    """
    return np.pi*(np.rad2deg(np.tan(2*Rstar/d_pl))/2)**2* ((np.pi/180)**2)

def black_body_rad(T,u):
    h   = 6.62607015e-34     # J s
    kb  = 1.380649e-23       # J/K
    c   = 299792458.0        # m/s
    c2  = h*c/kb
    """
    Return black body emission from body at temperature T at wavelength u in meters
    Units are in W/m2/um/sr
    Input parameters: temperature in K and wavelength in meters
    Output: single value in W/m2/um/sr  
    """
    B = 2*h*c**2/(u**5)*1/(np.exp(c2/(u*T))-1)/1e6
    return B

black_body_rad = np.vectorize(black_body_rad)

def find_nearest(array, value):
    """
    Return the closes value in the array
    Input parameters: array, value
    Output: single value from array
    """
    n = (np.abs(np.array(array)-value)).argmin()
    return array[n]

def find_nearest_l(array, value):
    """
    Returns the location of the closest value in the array
    Input parameters: array, value
    Output: integer
    """
    n = (np.abs(np.array(array)-value)).argmin()
    return n

# the number of days per month in a non-leap year
ms = [31,28,31, 30,31,30, 31,31,30, 31,30,31]

def sol_lat_d(day):
    """
    Returns the sub-solar latitude as a function fo the day in the year.
    Input: the day of the year
    Output: the sub-solar latitude.
    """
    
    return -23.44*np.cos(np.linspace(10*2*np.pi/365,2*np.pi+10*2*np.pi/365,365))[day]

def psg_img_to_2d_map(file,nlon,nlat,npts_expected= 100):
    """
    This converts the PSG hyperspectral images into an RGB image. 
    Input parameters: file location, number of longitudinal points, nr. latitudinal points, and the number of spectal points
    Output: 2D RGB image of size nlon by nlat.
    """

    data,pos,xval = open_psg_img(file,npts_expected)
    lat_image = np.linspace(90,-90,nlat)
    lon_image = sorted(np.linspace(180,-180,nlon,endpoint=False))
    image = np.zeros(shape=(nlat,nlon,3))
    for i in range(len(data)):
        lo = find_nearest_l(lon_image,pos[i,0])
        la = find_nearest_l(lat_image,pos[i,1])
        image[la,lo,:] = spec_to_rgb(xval,data[i])
    return image

def rebin_gcm(data,logs=False,sc=1):
    """
    Rebins the size of the MERRA-2 3D files (72,361,576 - vertical,latitudinal,longitudinal) into the size used for this study (72,91,144)
    Input parameters: 3D data file, log setting (required for molecular abundances), and the scaler if values have to be converted (e.g. Pascal to bar)
    Output: 3D array of size 72,91,144.
    """

    dnew = np.zeros([72,91,144], dtype=np.float32)
    npts = np.zeros([72,91,144], dtype=np.float32)
    for j in range(361):
        ilat = round(j/4)
        for h in range(72):
            ialt = 71-h
            for i in range(0,576):
                ilon = round(i/4) # Bin longitude by 4
                if ilon>=144: ilon=ilon-144
                if logs==True:
                    dnew[ialt,ilat,ilon] += np.log10(sc*data[h,j,i]+1e-30)
                else: 
                    dnew[ialt,ilat,ilon] += sc*data[h,j,i]
                npts[ialt,ilat,ilon] += 1
    return dnew/npts

def lamb_2dmap(lon_size,lat_size,sun_lon,sun_lat):
    """
    Provides a 2D map with the solar incidence angle at each of the lat/lon combinations.
    Input parameters: nr. longitudinal points, nr. latitudinal points, the sub-solar longitude, the sub-solar latitude.
    Output: A map of size lon,lat with the solar incidence angle. Negative values are replaced by 0.
    """
    longr,latgr = np.meshgrid(sorted(np.linspace(180,-180,lon_size,endpoint=False)),np.linspace(90,-90,lat_size))
    rmap  = np.cos(haversine_np(sun_lon,sun_lat,longr,latgr,1))
    rmap[rmap<0]=0
    return rmap

def lamb_globe(npix,oblon,oblat,sun_lon,sun_lat):
    """
    Provides a 3D projection of the solar incidence angle at each of the locations on the sphere.
    Input parameters: number of x and y pixels, the sub-observer longitude, the sub-observer latitude, the sub-solar longitue, the sub-solar latitude
    Output: 2D image of the globe with the values corresponding to the solar incidence angle.
    """

    lat_size = 361
    lon_size = 576
    proj = np.zeros(shape=(npix,npix))
    longr,latgr = np.meshgrid(sorted(np.linspace(180,-180,lon_size,endpoint=False)),np.linspace(90,-90,lat_size))
    lat_m,lon_m = pix_to_lat_lon_ds(npix,oblat,oblon)
    lat_image = np.linspace(90,-90,lat_size)
    lon_image = sorted(np.linspace(180,-180,lon_size,endpoint=False))
    rmap  = np.cos(haversine_np(sun_lon,sun_lat,longr,latgr,1))
    rmap[rmap<0]=0
    for x in range(npix):
        for y in range(npix):
            if np.sqrt((x-npix/2)**2+(y-npix/2)**2)>npix/2: continue
            lat = lat_m[y,x]
            lon = lon_m[y,x]
            proj[x,y] = rmap[find_nearest_l(lat_image,lat),find_nearest_l(lon_image,lon)]
    return proj

def pix_to_lat_lon_ds(npix,oblat,oblon):
    """
    Returns the latitude and longitude of the 2D projection of a globe of size npix*npix
    from perspective of the observed latitude (oblat) and longitude (oblon)
    Input parameters: number of x and y pixels, the sub-observer longitude, the sub-observer latitude.
    Output: two 2D arrays: containg the latitude/longitude of the x-y pixels of the projection.
    """
    lat_m = np.zeros(shape=(npix,npix))
    lon_m = np.zeros(shape=(npix,npix))
    ds_m  = np.zeros(shape=(npix,npix))
    k = 180.0/np.pi
    for x in range(npix):
        for y in range(npix):
            px = (x-(npix/2.0))/(npix/2.0)
            py = (y-(npix/2.0))/(npix/2.0)
            ro = np.sqrt(px*px + py*py)
            if ro>1: 
                ds_m[x,y]=np.nan
                continue
            if ro==0.0: ro=1e-9
            c = np.arcsin(ro)
            lon = oblon + k*np.arctan2(px*np.sin(c), (ro*np.cos(oblat/k)*np.cos(c) + py*np.sin(oblat/k)*np.sin(c)))
            while lon<-180.0: lon+=360.0
            while lon>=180.00: lon-=360.0
            # Calculate physical latitude
            lat = k*np.arcsin(np.cos(c)*np.sin(oblat/k) - py*np.sin(c)*np.cos(oblat/k)/ro)
            
            lat_m[x,y]      =  lat
            lon_m[x,y]      =  lon
            
    return lat_m,lon_m

def haversine_np(lon1, lat1, lon2, lat2, R):
    """
    Calculate the great circle distance between two points on a sphere.
    Input parameters: Longitude of point 1, latitude of point 1, longitude of point 2, latitude of point 1, Radius of the sphere
    Output: distance in the same units as R.
    """
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2
    return R*2 * np.arcsin(np.sqrt(a))



def map2d_to_globe(data,lat_l,lon_l,npix,oblat,oblon):
    """
    Projects a 2D map into a globe.
    Input parameters 2D dataframe, list of latitudes, list of longitudes, number of x and y pixels of image, the observed latitude, the observed longitude
    Output: 2D array.
    """

    lat_size = len(lat_l)
    lon_size = len(lon_l)
    #longr,latgr =np.meshgrid(np.linspace(-180,180,lon_size),np.linspace(90,-90,lat_size))
    lat_m,lon_m = pix_to_lat_lon_ds(npix,oblat,oblon)
    lat_image = np.linspace(90,-90,lat_size)
    lon_image = sorted(np.linspace(180,-180,lon_size,endpoint=False))
    globe = np.zeros(shape=(npix,npix))
    for x in range(npix):
        for y in range(npix):
            if np.sqrt((x-npix/2)**2+(y-npix/2)**2)>npix/2: continue
            lat = lat_m[y,x]
            lon = lon_m[y,x]
            globe[x,y] = data[find_nearest_l(lat_image,lat),find_nearest_l(lon_image,lon)] 
    return globe

def map_to_globe(data,lat_l, lon_l, npix,oblat,oblon,xval):
    """
    This function makes a RGB image of a 2D spaxels (data) of dimensions lat_l & lon_l.
    npix refers to the dimensions of the output image (only tested on 250*250)
    oblat, oblon refers to the observed latitude and longitude
    xval is the wavelength information of the spaxels 
    Input Parameters: data,lat_l, lon_l, npix,oblat,oblon,xval
    Output: a RGB image of size npix*npix*3 
    """
    lat_m,lon_m = pix_to_lat_lon_ds(npix,oblat,oblon)
    globe = np.zeros(shape=(npix,npix,3))
    for x in range(npix):
        for y in range(npix):
            if np.sqrt((x-npix/2)**2+(y-npix/2)**2)>npix/2: continue
            lat = lat_m[y,x]
            lon = lon_m[y,x]
            lon_sel = np.where(lon_l==find_nearest(lon_l,lon))
            lat_sub_sel = find_nearest_l(lat_l[lon_sel],lat)
            globe[x,y,:] = spec_to_rgb(xval, data[lon_sel[0][lat_sub_sel],:] )
    return globe

def spec_to_rgb(wvls,spec):
    """
    Converts a spectrum to a percieved color in RGB space
    Input parameters: wvls, wavelength bins of the spectrum (spec) in nanometer.
    Output: 3 values corresponding to the intensity of the RGB challens
    Function from Geronimo Villanueva.
    """

    spec= spec #/(1.4*np.max(spec))
    r = find_nearest_l(wvls,650)
    g = find_nearest_l(wvls,550)
    b = find_nearest_l(wvls,450)
    #return [spec[r]*1.4,spec[g]*1.2,spec[b]]
    return [spec[r],spec[g],spec[b]]


def boost_gamma(image,gamma,scale):
    """
    Boost the gamma of an RGB image by factor scale
    Input parameters: RGB image, gamma value, list of 3 scaling values which the RGB channels should be scaled by
    Output: RGB image.
    """
    image = np.copy(image)
    image[:,:,0] = np.clip(pow(image[:,:,0]/scale[0], gamma),0,255)
    image[:,:,1] = np.clip(pow(image[:,:,1]/scale[1], gamma),0,255)
    image[:,:,2] = np.clip(pow(image[:,:,2]/scale[2], gamma),0,255)
    return image       

def call_api(config_path: str, psg_url: str = 'http://localhost:3000',
    # call_psg function by Ted Johnson #
             api_key: str = None, output_type: str = None, app: str = None,
             outfile: str = None) -> None:
    """
    Call the PSG api
    Build and execute an API query to communicate with PSG.

    Parameters
    ----------
    config_path : str or pathlib.Path
        The path to the `PSG` config file.
    psg_url : str, default='https://psg.gsfc.nasa.gov'
        The URL of the `PSG` API. Use 'http://localhost:3000' if running locally.
    api_key : str, default=None
        The key for the public API. Needed only if not runnning `PSG` locally.
    output_type : str, default=None
        The type of output to retrieve from `PSG`. Options include 'cfg', 'rad',
        'noi', 'lyr', 'all'.
    app : str, default=None
        The PSG app to call. For example: 'globes'
    outfile : str, default=None
        The path to write the PSG output.
    """
    data = {}
    with open(config_path,'rb') as file:
        dat = file.read()
    data['file'] = dat
    if api_key is not None:
        data['key'] = api_key
    if app is not None:
        data['app'] = app
    if output_type is not None:
        data['type'] = output_type
    data['option'] = '-s'

    url = f'{psg_url}/api.php'
    reply = requests.post(url,data=data,timeout=5555555)
    if outfile is not None:
        with open(outfile,'wb') as file:
            file.write(reply.content)
    return reply

def psg_img_to_globe(file,oblat,oblon,npts_expected):
    """
    Input parameters:
    file, oblat, oblon, npts_expected
    Opens PSG-GlobES hyperspectral images and returns a 2D RGB image by using the map_to_globe routine.
    file: link to file,
    oblat, oblon, observed latitude/longitude
    npts_expected: the number of spectral points in the file (in principle this could be read, but if there are errors or warnings this makes it more complicated.)
    """
    fsz = os.path.getsize(file)
    with open (file,'rb') as f: 
        if b'ERROR' in f.read(7): print('It seesm the simulations crashed - there is a PSG ERROR indicating PUMAS crashed'); return 
    fr = open(file,'rb')
    npts = struct.unpack('i', fr.read(4))[0]
    while npts != npts_expected:
        next(fr)
        fsz = fsz - fr.tell()
        npts = struct.unpack('i', fr.read(4))[0]
    xval = struct.unpack('f'*npts, fr.read(4*npts))
    nspec = int((fsz-npts*4-4)/(4.0*(npts+5)))
    pos = np.zeros([nspec,5])
    data = np.zeros([nspec,npts])
    for i in range(nspec):
        pos[i,:] = struct.unpack('f'*5, fr.read(4*5))
        data[i,:] = struct.unpack('f'*npts, fr.read(4*npts))
    fr.close()
    pos2 = np.copy(pos[:,0])
    for i in range(len(pos2)):
        if pos2[i]>180:
            pos2[i]-=360
    return  map_to_globe(data, pos[:,1],pos2,250, oblat,oblon,xval)

def map_to_weight(data,lat_l, lon_l, npix,oblat,oblon):
    """
    Provides the weight of each of the spectra in the hyperspectal image in order to combine these into a single-pixel spectrum
    The spectra are weighed as a function of their visibility on a 2D projection of size 250*250.
    Input parameters: spectral data (only used for the length), list of latitudes, list of longitude, npix, sub-observer latitude, sub-observer longitude
    Output: a list with the weights of each of the lat/lon combinations
    Works together closely with 'open_psg_img':
    Example:
        datan ,pos,xval = open_psg_img(tfile)
        weights,spec_map = map_to_weight(datan,pos[:,1],pos[:,0],250,oblat,oblon)
    """
    lat_m,lon_m = pix_to_lat_lon_ds(npix,oblat,oblon)
    weights = np.zeros(shape=np.shape(data)[0])
    spec_map = -1*np.ones(shape=(npix,npix))
    for x in range(npix):
        for y in range(npix):
            if np.sqrt((x-npix/2)**2+(y-npix/2)**2)>npix/2: continue
            lat = lat_m[y,x]
            lon = lon_m[y,x]
            lon_sel = np.where(lon_l==find_nearest(lon_l,lon))
            lat_sub_sel = find_nearest_l(lat_l[lon_sel],lat)
            spec_map[x,y]= lon_sel[0][lat_sub_sel]
            weights[lon_sel[0][lat_sub_sel]] +=1
    return weights,spec_map

def open_psg_img(file,npts_expected= 100):
    """
    Opens the PSG hyperspectral datafiles.
    Input parameters: PSG hyper spectral file, number of wavelenght points.
    Output: array containing all the spectra, the position vectors, and the range of wavelenght values corresponding to the spectra.
    The position vector contains 5 values per spectra (longitude, latitude, 0 (currenly unused),  x position of the 2D projection, y position of the 2D projection ) 
    """

    fsz = os.path.getsize(file)
    with open (file,'rb') as f: 
        if b'ERROR' in f.read(7): print('More binary errros than just the warning'); return 
    fr = open(file,'rb')
    npts = struct.unpack('i', fr.read(4))[0]
    while npts != npts_expected:
        next(fr)
        fsz = fsz - fr.tell()
        npts = struct.unpack('i', fr.read(4))[0]
    xval = struct.unpack('f'*npts, fr.read(4*npts))
    nspec = int((fsz-npts*4-4)/(4.0*(npts+5)))
    pos = np.zeros([nspec,5])
    data = np.zeros([nspec,npts])
    for i in range(nspec):
        pos[i,:] = struct.unpack('f'*5, fr.read(4*5))
        data[i,:] = struct.unpack('f'*npts, fr.read(4*npts))
    fr.close()
    pos2 = np.copy(pos[:,0])
    for i in range(len(pos2)):
        if pos2[i]>180:
            pos2[i]-=360
    pos[:,0] = pos2
    #print(np.median(pos2))
    if np.shape(data)[0] != np.shape(pos)[0]:
        print('Error, data and positions are not the same')
        return
    return  data, pos, xval


def get_image(band):
    """
    Extract image from an EPIC observation (example taken from EPIC team)
    """

    # Non-Earth pixels are infinity, fill them to black
    img = np.ma.fix_invalid(np.array(band['Image']), fill_value=0)

    # Flip/rotate the image so North is up:
    img = np.fliplr(np.rot90(img, -1))
    return img

def build_rgb(red, green, blue):
    """
    Build RGB image from an EPIC observation (example taken from EPIC team) - not used anymore
    """
    # Apply some simple colour correction
    green *= 0.9
    red *= 1.25
    # Normalise the image to a floating point 0-1 range, and stack the images into a 3D array
    max_value = max(red.max(), green.max(), blue.max())
    rgb = np.dstack((red/max_value, green/max_value, blue/max_value))
    # Exposure correction:
    rgb = np.clip(rgb * 1.8, 0, 1)
    return rgb

def bin_image(image,factor):
    """
    Bin an image down 
    Input: image and the binning factor
    Output: 2D image
    """

    dims_o = np.shape(image)
    dims_n  = (round(dims_o[0]/factor), round(dims_o[1]/factor),dims_o[2])
    print(dims_n)
    newim = np.zeros(shape=dims_n)
    couim = np.zeros(shape=dims_n)
    for i in range(dims_o[0]):
        ib = round(i/factor)
        if ib>=dims_n[0]:continue
        for j in range(dims_o[1]):
            jb = round(j/factor)
            if jb>=dims_n[1]:continue
            newim[ib,jb,:] +=image[i,j,:]
            couim[ib,jb,:] +=1
    return newim/couim

def get_time_nr(time_in):
    """
    Gets the number of the closest observation in the MERRA dataset based on the given time of day.
    Input: time string in format (HH:MM)
    Output: singe integer
    """

    time_in = time_in.split(':')
    tims = []
    for i in range(8):
        tims.append([i,abs(datetime.datetime(2008,3,18,1+3*i,30)-datetime.datetime(2008,3,18,int(time_in[0]),int(time_in[1])))])
    tims=np.array(tims)
    return np.where(tims[:,1]==np.min(tims[:,1]))[0][0]

def get_obs_c(file, sol=False):
    """
    Returns the OBS-latitude and longitude from a configuration file
    if sol=TRUE:
    Returns the OBS and SOLAR longitude and latitude from a configuarion file.
    Input: configuration file
    """

    f =  open(file, "r"); 
    cfg = f.readlines(); 
    f.close()
    obs_lon,obs_lat,sol_lat,sol_lon = np.nan,np.nan,np.nan,np.nan
    for line in cfg: 
        if '<OBJECT-OBS-LONGITUDE>' in line: obs_lon = float(line.split('>')[-1])
        if '<OBJECT-OBS-LATITUDE>' in line: obs_lat = float(line.split('>')[-1])
        if sol is True:
            if '<OBJECT-SOLAR-LATITUDE>' in line: sol_lat = float(line.split('>')[-1])
            if '<OBJECT-SOLAR-LONGITUDE>' in line: sol_lon = float(line.split('>')[-1])
    if sol is True:
        return obs_lon,obs_lat,sol_lon,sol_lat
    else: 
        return obs_lon,obs_lat
    

def psg_img_to_spec_wc(file_clear,file_cloud,cloud_map,obs_lat,obs_lon,nspec=85):
    """
    Returns the disk integrated spectrum, combined from the clear and cloud files using the cloud coverage map.
    Input: clear file location, cloudy file location, cloud coverage projection, obs_latitude, obs_longitude, nspec       
    """


    cm =  np.flip(np.genfromtxt(cloud_map),axis=0)
    datan, pos, xval   = open_psg_img(file_clear,nspec)
    datac, pos, xval   = open_psg_img(file_cloud,nspec)
    weights, spec_map  = map_to_weight(datan,pos[:,1],pos[:,0],250,obs_lat,obs_lon)
    c = 0
    npix=250
    av_spec = np.zeros(shape=(85))
    for x in range(npix):
        for y in range(npix):
            spec = int(spec_map[x,y])
            if spec < 0: continue
            else: 
                av_spec = np.nansum([av_spec,datac[spec,:]*cm[x,y]+datan[spec,:]*(1-cm[x,y])],axis=0)
                c+=1
    return np.array([xval,av_spec/c]).T


def psg_img_to_spec(file_clear,obs_lat,obs_lon,nspec=85):
    """
    Returns the disk integrated spectrum, combined from the clear and cloud files using the cloud coverage map.
    Input: clear file location, cloudy file location, cloud coverage projection, obs_latitude, obs_longitude, nspec       
    """
    #cm =  np.flip(np.genfromtxt(cloud_map),axis=0)
    datan, pos, xval   = open_psg_img(file_clear,nspec)
    #datac, pos, xval   = open_psg_img(file_cloud,nspec)
    weights, spec_map  = map_to_weight(datan,pos[:,1],pos[:,0],250,obs_lat,obs_lon)
    c = 0
    npix=250
    av_spec = np.zeros(shape=(85))
    for x in range(npix):
        for y in range(npix):
            spec = int(spec_map[x,y])
            if spec < 0: continue
            else: 
                av_spec = np.nansum([av_spec,datan[spec,:]],axis=0)
                c+=1
    return np.array([xval,av_spec/c]).T