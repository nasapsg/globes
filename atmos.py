# ---------------------------------------------------------------
# Script to convert ATMOS model outputs to a PSG configuration file
# Fauchez, Villanueva - NASA Goddard Space Flight Center
# February 2021
# ---------------------------------------------------------------
import sys, os
import numpy as np
import matplotlib.pyplot as plt

# Module to convert ATMOS to PSG atmospheric file
def atmosatm(ptfile = 'profile.pt', hcfile = 'hcaer.out', filebase = ''):
    # Read the atmos profiles ------------------------------
    data = np.genfromtxt(ptfile, skip_header=1)
    if data.shape[1]!=16:
        print("Invalid file dimensions: {}.".format(data.shape))
        exit()
    #Endif
    bin = 4
    mols = 'Altitude,H2O,CH4,C2H6,CO2,O2,O3,CO,H2CO,HNO3,NO2,SO2,N2O,N2'; ncols=16
    data[:,0] = data[:,0]/1e5   # Convert altitude from cm to km
    dens = data[::bin,2]        # Density [molecules/cm3]
    nalt = int(data.shape[0]/bin)
    gprof = np.zeros([nalt,19])
    gprof[:,0:15] = data[::bin,[3,1,0,4,5,6,7,8,9,10,11,12,13,14,15]]
    gprof[:,15] = 1.0 - np.sum(gprof[:,3:],axis=1) # Compute N2 abundance
    gprof[:,17] = 1e-6
    mmol = [18.0,16.0,30.0,44.0,32.0,48.0,28.0,30.0,63.0,46.0,64.0,44.0,28.0]
    gatm = np.sum(gprof[1,3:16] * mmol)/np.sum(gprof[1,3:16])

    # Read the haze profile --------------------------------
    if os.path.isfile(hcfile):
        data = np.genfromtxt(hcfile, skip_header=4)
        if data.shape[1]!=9:
            print("Invalid file dimensions: {}.".format(data.shape))
            exit()
        #Endif
        rho = 0.64                                               # Particle density [g/cm3]
        data[:,0] = data[:,0]/1e5                                # Convert altitude from cm to km
        naero = np.interp(gprof[:,2], data[:,0], data[:,1])      # Particle abundance [#/cm3]
        raero = np.interp(gprof[:,2], data[:,0], data[:,2])      # Particle radius [m]
        vaero = np.pi*(4.0/3.0)*raero**3.0                       # Volume of particle [cm3]
        maero = naero*vaero*rho                                  # Haze mass density [g/cm3]
        xaero = maero/(gatm*dens/6.0221409e+23)                  # Haze abundance [g/g]
        gprof[:,17] = xaero
        gprof[:,18] = raero/1e2
        mols = '%s,Haze,Haze_size' % mols; ncols=19
    #Endif

    # Write configuration file --------------------------------
    newf = []
    newf.append('<ATMOSPHERE-DESCRIPTION>ATMOS Photochemistry')  # Description establishing the source/reference for the vertical profile
    newf.append('<ATMOSPHERE-STRUCTURE>Equilibrium')             # The structure of the atmosphere, None / Equilibrium:'Hydrostatic equilibrium' / Coma:'Cometary expanding coma'
    newf.append('<ATMOSPHERE-PRESSURE>%.3f' % gprof[0,0])        # For equilibrium atmospheres, this field defines the surface pressure; while for cometary coma, this field indicates the gas production rate
    newf.append('<ATMOSPHERE-PUNIT>bar')                         # The unit of the ATMOSPHERE-PRESSURE field
    newf.append('<ATMOSPHERE-WEIGHT>%.3f' % gatm)                # Molecular weight of the atmosphere [g/mol] or expansion velocity [m/s] for expanding atmospheres
    newf.append('<ATMOSPHERE-NGAS>11')                           # Number of gas in Atmos profile.pt
    newf.append('<ATMOSPHERE-GAS>N2,H2O,CH4,CO2,O2,O3,CO,H2CO,NO2,SO2,N2O') # Name of the gases to include in the simulation,
    newf.append('<ATMOSPHERE-TYPE>HIT[22],HIT[1],HIT[6],HIT[2],HIT[7],HIT[3],HIT[5],HIT[20],HIT[10],HIT[9],HIT[4]') # Sub-type of the gases, e.g. 'HIT[1], HIT[2]'
    newf.append('<ATMOSPHERE-ABUN>1,1,1,1,1,1,1,1,1,1,1')      # Value of the scaler that would multiply the gas profile
    newf.append('<ATMOSPHERE-UNIT>scl,scl,scl,scl,scl,scl,scl,scl,scl,scl,scl') # Scaler unit

    # Add haze parameters
    if ncols>17:
        newf.append('<ATMOSPHERE-NAERO>1')                       # Number of aerosols to include in the simulation
        newf.append('<ATMOSPHERE-AEROS>Haze')                    # Name of the aerosols to include in the simulation
        newf.append('<ATMOSPHERE-ATYPE>KHARE_TITAN_THOLINS_HRI') # Sub-type of the aerosols
        newf.append('<ATMOSPHERE-AABUN>1')                       # Abundance of aerosols. The values can be assumed to be same across all altitudes/layers [%,ppm,ppb,ppt,Kg/m2], or as a multiplier [scaler] to the provided vertical profile
        newf.append('<ATMOSPHERE-AUNIT>scl')                     # Unit of the ATMOSPHERE-AABUN field, % / ppmv / ppbv / pptv / m2:'molecules/m2' / scl:'scaler of profile
        newf.append('<ATMOSPHERE-ASIZE>1')                       # Effective radius of the aerosol particles. The values can be assumed to be same across all layers [um, m, log(um)], or as a multiplier [scaler] to the provided size vertical profile
        newf.append('<ATMOSPHERE-ASUNI>scl')                     # Unit of the size of the aerosol particles
    else:
        newf.append('<ATMOSPHERE-NAERO>0')                       # Number of aerosols to include in the simulation
    #Endif

    # Add profile
    newf.append('<ATMOSPHERE-LAYERS-MOLECULES>%s' % mols)        # Molecules and aerosols quantified by the vertical profile
    newf.append('<ATMOSPHERE-LAYERS>%d' % nalt)                  # Number of layers of the atmospheric vertical profile
    for i in range(nalt):                                        # Values for that specific layer: Pressure[bar], Temperature[K], gases[mol/mol], aerosols [kg/kg] - Optional fields: Altitude[km], X_size[m, aerosol X particle size]
        str = '%.5e' % gprof[i,0]
        for j in range(1,ncols): str = '%s,%.5e' % (str,gprof[i,j])
        newf.append('<ATMOSPHERE-LAYER-%d>%s' % (i+1,str))
    #Endfor

    # Add surface parameters
    newf.append('<SURFACE-TEMPERATURE>%.3f' % gprof[0,1])        # Temperature of the surface [K]
    newf.append('<SURFACE-ALBEDO>0.25')                          # Albedo the surface [0:non-reflectance, 1:fully-reflective]
    newf.append('<SURFACE-EMISSIVITY>1.0')                       # Emissivity of the surface [0:non-emitting, 1:perfect-emitter]
    newf.append('<SURFACE-NSURF>0')                              # Number of components describing the surface properties [areal mixing]

    # Save atmospheric file
    if len(filebase):
        fw = open("%s_atm.txt" % filebase,'w')
        for line in newf: fw.write('%s\n' % line)
        fw.close()
    #Endif
    return newf
#End module atmosatm

# Module to compute spectra from a config-file and to plot it
def psgspec(config=[], filebase = 'atmos_psg', showplot=True):
    #psgurl = 'http://localhost:3000' # URL of the PSG server - For PSG/Docker
    #psgurl = 'http://localhost' # URL of the test PSG server
    psgurl = 'https://psg.gsfc.nasa.gov' # URL of the PSG server

    # Save configuration and run via the API
    fw = open("%s_cfg.txt" % filebase,'w')
    for line in newf: fw.write('%s\n' % line)
    fw.close()
    os.system('curl -s --data-urlencode file@%s_cfg.txt %s/api.php > %s_rad.txt' % (filebase, psgurl, filebase))

    # Plot results
    if showplot:
        xunit = 'Wavelength [um]'; yunit = 'Contrast [ppm]'; cols = 'Total'
        fw = open("%s_rad.txt" % filebase); lines = fw.readlines(); fw.close()
        for line in lines:
            if line[0]!='#': break
            if line[0:16]=='# Spectral unit:': xunit = line[17:-1]
            if line[0:16]=='# Radiance unit:': yunit = line[17:-1]
            if line[0:11]=='# Wave/freq': cols = line[12:-1]
        #Endfor
        data = np.genfromtxt("%s_rad.txt" % filebase)
        cols = cols.split(); wnoise=0
        plt.figure(figsize=[10,5])
        for i in range(len(cols)):
            if cols[i]=='Noise': wnoise=i; continue
            plt.plot(data[:,0],data[:,i+1],label=cols[i])
        #Endfors
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(xunit)
        plt.ylabel(yunit)
        plt.ylim([max(data[:,1])/100.0,max(data[:,1])*1.2])
        plt.xlim([min(data[:,0]),max(data[:,0])])
        plt.legend()
        plt.title('ATMOS Atmosphere synthesized with NASA/PSG')
        plt.tight_layout()
        plt.savefig('%s_rad.png' % filebase)
        plt.show()
    #Endif
#End

# Main module if the user is calling this as a script
if __name__ == "__main__":
    # Parameters of the simulation
    diameter = 0.920 * 12742.0   # Planet's diameter [km], 12742 km is Earth's diameter
    gravity  = 0.930 * 9.81      # Planet's gravity [m/2], 9.81 m/s2 is Earth's gravity

    # Orbital parameters
    distance = 12.4              # Distance to the star/planet from Earth [pc]
    semiaxis = 28.17e-3          # Distance planet-star [AU]
    period = 6.1                 # Orbital period [days]
    ttransit = 3345.0            # Transit duration [seconds]
    phase = 180                  # Orbital phase, 180:Primary_transit, 0:Secondary, 90/270:Quadrature

    # Stellar parameters
    stype = 'M'                  # Spectral type
    stemp = 2566.0               # TRAPPIST-1 temperature [K] from Agol et al 2020
    sradius = 0.1192             # TRAPPIST-1 radius [Sun radius unit] from Agol et al 2020

    # Integration time
    ntransits = 1.0              # Number of transits

    # Convert ATMOS files to a configuration file
    newf = atmosatm('data/atmos_profile.pt', 'data/atmos_hcaer.out')

    # Add geometrical/observational parameters for the simulation
    newf.append('<OBJECT>Exoplanet')                             # Object type (e.g., Exoplanet) or object name for the main bodies in the Solar System
    newf.append('<OBJECT-NAME>Exoplanet')                        # Object name
    newf.append('<OBJECT-DIAMETER>%.1f' % diameter)              # Diameter of the object [km]
    newf.append('<OBJECT-GRAVITY>%.2f' % gravity)                # Gravity/density/mass of the object
    newf.append('<OBJECT-GRAVITY-UNIT>g')                        # Unit for the OBJECT-GRAVITY field, g:'Surface gravity [m/s2]', rho:'Mean density [g/cm3]', or kg:'Total mass [kg]'
    newf.append('<OBJECT-STAR-DISTANCE>%.5f' % semiaxis)         # Distance of the planet to the Sun [AU], and for exoplanets the semi-major axis [AU]
    newf.append('<OBJECT-STAR-VELOCITY>0.0')                     # Velocity of the planet to the Sun [km/s], and for exoplanets the RV amplitude [km/s]
    newf.append('<OBJECT-STAR-TYPE>%s' % stype)                  # Stellar type of the parent star [O/B/A/F/G/K/M]
    newf.append('<OBJECT-STAR-TEMPERATURE>%.1f' % stemp)         # Temperature of the parent star [K]
    newf.append('<OBJECT-STAR-RADIUS>%.4f' % sradius)            # Radius of the parent star [Rsun]
    newf.append('<OBJECT-SOLAR-LONGITUDE>0.0')                   # Sub-solar east longitude [degrees]
    newf.append('<OBJECT-SOLAR-LATITUDE>0.0')                    # Sub-solar latitude [degrees]
    newf.append('<OBJECT-STAR-METALLICITY>0.0')                  # Metallicity of the parent star and object with respect to the Sun in log [dex]
    newf.append('<OBJECT-OBS-LONGITUDE>%.1f' % phase)            # Sub-observer east longitude [degrees]
    newf.append('<OBJECT-OBS-LATITUDE>0.0')                      # Sub-observer latitude, for exoplanets inclination [degrees]
    newf.append('<OBJECT-INCLINATION>90.00')                     # Orbital inclination [degree], mainly relevant for exoplanets. Zero is phase on, 90 is a transiting orbit
    newf.append('<OBJECT-SEASON>%.1f' % phase)                   # Angular parameter (season/phase) that defines the position of the planet moving along its Keplerian orbit. For exoplanets, 0:Secondary transit, 180:Primary transit, 90/270:Opposition.
    newf.append('<OBJECT-PERIOD>%.2f' % period)                  # Apparent rotational period of the object as seen from the observer [days]

    newf.append('<GEOMETRY>Observatory')                         # Type of observing geometry
    newf.append('<GEOMETRY-OFFSET-UNIT>arcsec')                  # Unit of the GEOMETRY-OFFSET field, arcsec / arcmin / degree / km / diameter
    newf.append('<GEOMETRY-OBS-ALTITUDE>%.3f' % distance)        # Distance between the observer and the surface of the planet
    newf.append('<GEOMETRY-ALTITUDE-UNIT>pc')                    # Unit of the GEOMETRY-OBS-ALTITUDE field, AU / km / diameter and pc:'parsec'
    newf.append('<GEOMETRY-PLANET-FRACTION>1.0')                 # This field is computed by the geometry module - It indicates how much the beam fills the planetary area (1:maximum)
    newf.append('<GEOMETRY-PHASE>%.1f' % phase)                  # This field is computed by the geometry module - It indicates the phase between the Sun and observer
    newf.append('<GEOMETRY-REF>User')                            # Reference geometry (e.g., ExoMars, Maunakea), default is user defined or 'User'

    # Save instrument parameters
    newf.append('<GENERATOR-INSTRUMENT>user')                    # Text describing if an instrument template was used to define the GENERATOR parameters
    newf.append('<GENERATOR-RANGE1>0.6')                         # Lower spectral range for the simulation
    newf.append('<GENERATOR-RANGE2>20')                          # Upper spectral range for the simulation
    newf.append('<GENERATOR-RANGEUNIT>um')                       # Unit of the GENERATOR-RANGE fields, um / nm / mm / An:'Angstrom' / cm:'Wavenumber [cm-1]' / MHz / GHz / kHz
    newf.append('<GENERATOR-RESOLUTION>500')                     # Spectral resolution for the simulation.
    newf.append('<GENERATOR-RESOLUTIONUNIT>RP')                  # Unit of the GENERATOR-RESOLUTION field, RP:'Resolving power' / um / nm / mm / An:'Angstrom' / cm:'Wavenumber [cm-1]' / MHz / GHz / kHz
    newf.append('<GENERATOR-TELESCOPE>SINGLE')                   # Type of telescope, SINGLE:'single dish telescope or instrument', ARRAY:'Interferometric array', CORONA:'Coronagraph', AOTF or LIDAR
    newf.append('<GENERATOR-TELESCOPE2>0.0')                     # This field indicates the zodi-level (1.0:Ecliptic pole/minimum, 2.0:HST/JWST low values, 10.0:Normal values, 100.0:Close to ecliptic/Sun), or order number for the AOTF system. For coronagraphs, this field indicates allows two entries: the exozodi level and the local zodiacal dust level
    newf.append('<GENERATOR-DIAMTELE>5.6')                       # Diameter of the main reflecting surface of the telescope or instrument [m]
    newf.append('<GENERATOR-BEAM>1.0')                           # Full width half-maximum (FWHM) of the instrument's beam or field-of-view (FOV)
    newf.append('<GENERATOR-BEAM-UNIT>diffrac')                  # Unit of the GENERATOR-BEAM field, arcsec / arcmin / degree / km / diameter:'beamsize in terms of the planet's diameter' / diffrac:'defined by the telescope diameter and center wavelength'
    newf.append('<GENERATOR-NOISE>CCD')                          # noise model to consider
    newf.append('<GENERATOR-NOISE1>16.8')                        # First noise model parameter - For RMS, 1-sigma noise; for TRX, the receiver temperature; for BKG, the 1-sigma noise; for NEP, the sensitivity in W/sqrt(Hz); for DET, the sensitivity in cm.sqrt(Hz)/W; for CCD, the read noise [e-]
    newf.append('<GENERATOR-NOISE2>0.005')                       # Second noise model parameter - For RMS, not used; for TRX, the sideband g-factor; for BKG, the not used; for NEP, not used; for DET, the pixel size [um]; for CCD, the dark rate [e-/s]
    newf.append('<GENERATOR-NOISEOTEMP>50')                      # Temperature of the telescope+instrument optics [K]
    newf.append('<GENERATOR-NOISEOEFF>0.4')                      # Total throughput of the telescope+instrument, from photons arriving to the main mirror to photons being quantified by the detector [0:none to 1:perfect].
    newf.append('<GENERATOR-NOISEOEMIS>0.10')                    # Emissivity of the telescope+instrument optics [0 to 1]
    newf.append('<GENERATOR-NOISETIME>0.2')                      # Exposure time per frame [sec]
    newf.append('<GENERATOR-NOISEFRAMES>%.1f' % (ttransit*ntransits/0.2)) # Number of exposures
    newf.append('<GENERATOR-NOISEPIXELS>8')                      # Total number of pixels that encompass the beam (GENERATOR-BEAM) and the spectral unit (GENERATOR-RESOLUTION)
    newf.append('<GENERATOR-TRANS-APPLY>N')                      # Flag indicating whether to show the spectra as observed and multiplied by the telluric transmittance [Y/N]
    newf.append('<GENERATOR-TRANS-SHOW>N')                       # Flag indicating whether we are synthesizing planetary spectra as observed with a ground-based telescope. This flag will ensure that the noise module properly includes telluric signatures [Y/N]
    newf.append('<GENERATOR-GAS-MODEL>Y')                        # Flag indicating whether to include molecular signatures as generated with PUMAS or CEM [Y/N]
    newf.append('<GENERATOR-CONT-MODEL>Y')                       # Flag indicating whether to include continuum signatures as generated by the surface, the star (when in the field) and dust/nucleus (when synthesizing comets) [Y/N]
    newf.append('<GENERATOR-CONT-STELLAR>Y')                     # Flag indicating whether to include stellar absorption signatures in the reflected sunlight / stellar spectra [Y/N]
    newf.append('<GENERATOR-RADUNITS>ppm')                       # Radiation unit for the generated spectra, see full list of permitted keywords in the 'Modeling > Radiation Units' sectio

    # Run simulation
    psgspec(newf)
#End main modules
