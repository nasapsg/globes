# -----------------------------------------
# Nose parameterization code
# More details in Kofman et al. 2024, Planetary Science Journal
# -----------------------------------------
import os
import numpy as np
import scipy as sc
import pandas as pd

rEarth = 6378100    # meters
Rsun = 695700       # kilometers 
AU = 1.496e8        # kilometers
d=AU*1              # d in AU
h = 6.62607015e-34  # Planck's constant [W s2] or [J s]t
kb = 1.380649e-23   # Boltzman constant [J/K]
c = 299792458.0     # speed of light m/s
c2 = h*c/kb         # second radiation constant [K/cm-1]
A = 0.31            # average albedo Earth


def black_body_rad(T,u):
    """
    Return black body emission from body at temperature T at wavelength u in meters
    Units are in W/m2/um/sr
    """
    B = 2*h*c**2/(u**5)*1/(np.exp(c2/(u*T))-1)/1e6
    return B
black_body_rad = np.vectorize(black_body_rad)

def find_nearest_l(array, value):
    n = (np.abs(np.array(array)-value)).argmin()
    return n

def star_sr_pl(Rstar,d_pl):
    d_pl *=  1.496e8 # km/AU
    """
    input: rstar in kilometers, d_pl in AU
    Returns the size of the star in the sky in steradians at the distance of the planet from the star
    """
    return np.pi*(np.rad2deg(np.tan(2*Rstar/d_pl))/2)**2* ((np.pi/180)**2)

def remove_feature(ylist, a,b):
    """Replaces a spectral feature in ylist with a linear interpolation between point a and b"""
    ylist_return = np.copy(ylist)
    ylist_return[a:b] = np.linspace(ylist[a],ylist[b-1],b-a)
    return ylist_return

def detection_strength_feature(ylist1,ylist2,noise):
    """
    Return the signal-to-noise ratio of the difference between ylist1 and ylist2 considering the noise
    """
    return np.sqrt(np.sum( ((ylist1-ylist2)/noise)**2  ))

# Luvoir A/B
noef        = '<GENERATOR-NOISEOEFF>0.0317@0.2000,0.0437@0.2261,0.0589@0.2580,0.0742@0.2986,0.0851@0.3377,0.0917@0.3667,0.0971@0.4029,0.1015@0.4493,0.1004@0.4971,0.1004@0.5140,0.1670@0.5150,0.1659@0.5377,0.1506@0.6304,0.1255@0.7087,0.0939@0.7986,0.0884@0.8435,0.1146@0.9058,0.1419@0.9594,0.1594@0.9942,0.1821@1.2200,0.1958@1.4100,0.2049@1.6200,0.2094@1.8700,0.2140@2.0000'
nod         = '<GENERATOR-NOISE2>3e-5@0.2,3e-5@1,2e-3@1.01,2e-3@2.0'
nor         = '<GENERATOR-NOISE1>0@0.2,0@1,2.5@1.01,2.5@2.0'
trp         = np.genfromtxt(noef.split('>')[1].split(','),delimiter='@')
drn         = np.genfromtxt(nor.split('>')[1].split(','),delimiter='@')
ddc         = np.genfromtxt(nod.split('>')[1].split(','),delimiter='@')
L_tpsf      = sc.interpolate.interp1d(trp[:,1],trp[:,0], kind='linear',bounds_error=False)
L_ddcf      = sc.interpolate.interp1d(ddc[:,1],ddc[:,0], kind='linear',bounds_error=False)
L_drnf      = sc.interpolate.interp1d(drn[:,1],drn[:,0], kind='linear',bounds_error=False)

#HabEx Starshade
noef         = '<GENERATOR-NOISEOEFF>0.2260@0.2000,0.2201@0.2112,0.2182@0.2224,0.2300@0.2587,0.2536@0.3035,0.2673@0.3399,0.2791@0.3874,0.2830@0.4294,0.2850@0.4490,0.1796@0.4500,0.1848@0.5161,0.1749@0.6503,0.1474@0.7734,0.1415@0.8434,0.1592@0.9021,0.1926@0.9608,0.1988@0.9750,0.1988@0.9760,0.2162@1.0587,0.2339@1.2462,0.2437@1.4503,0.2516@1.6657,0.2594@1.8000'
nod          = '<GENERATOR-NOISE2>3e-5@0.2,3e-5@0.975,0.005@0.976,0.005@1.8'
nor          = '<GENERATOR-NOISE1>0.008@0.2,0.008@0.975,0.32@0.976,0.32@1.8'
trp          = np.genfromtxt(noef.split('>')[1].split(','),delimiter='@')
drn          = np.genfromtxt(nor.split('>')[1].split(','),delimiter='@')
ddc          = np.genfromtxt(nod.split('>')[1].split(','),delimiter='@')
HS_tpsf      = sc.interpolate.interp1d(trp[:,1],trp[:,0], kind='linear',bounds_error=False)
HS_ddcf      = sc.interpolate.interp1d(ddc[:,1],ddc[:,0], kind='linear',bounds_error=False)
HS_drnf      = sc.interpolate.interp1d(drn[:,1],drn[:,0], kind='linear',bounds_error=False)

# HabEx coronagraph
noef     = '<GENERATOR-NOISEOEFF>0.0000@0.3459,0.0038@0.3519,0.0301@0.3608,0.0473@0.3638,0.0659@0.3698,0.0710@0.3817,0.0774@0.3996,0.0850@0.4175,0.1023@0.4473,0.1090@0.4500,0.1113@0.4950,0.1093@0.5726,0.1023@0.6561,0.0914@0.7515,0.0831@0.7992,0.0806@0.8410,0.0844@0.8946,0.0895@0.9394,0.0902@0.9602,0.0902@0.9750,0.0902@0.9760,0.0902@0.9930,0.0946@1.1243,0.1036@1.3718,0.1113@1.6044,0.1157@1.8000'
trp       = np.genfromtxt(noef.split('>')[1].split(','),delimiter='@')
HC_tpsf      = HS_tpsf
HC_ddcf      = HS_ddcf
HC_drnf      = sc.interpolate.interp1d(drn[:,1],drn[:,0], kind='linear',bounds_error=False)

# coronagraph throughput values
HC = np.genfromtxt('./telescopes/HabEx_HCG_coronagraph_throughput.txt',delimiter=',')
HS = np.genfromtxt('./telescopes/HabEx_SS_throughput.txt',delimiter=',')
LA = np.genfromtxt('./telescopes/Luvoir_A_coronagraph_throughput.txt',delimiter=',')
LB = np.genfromtxt('./telescopes/Luvoir_B_coronagraph_throughput.txt',delimiter=',')

# made into continuous functions
HCf_ci = sc.interpolate.interp1d(HC[:,0], HC[:,1],kind='linear',bounds_error=False,fill_value='extrapolate')
HSf_ci = sc.interpolate.interp1d(HS[:,0], HS[:,1],kind='linear',bounds_error=False,fill_value='extrapolate')
LAf_ci = sc.interpolate.interp1d(LA[:,0], LA[:,1],kind='linear',bounds_error=False,fill_value='extrapolate')
LBf_ci = sc.interpolate.interp1d(LB[:,0], LB[:,1],kind='linear',bounds_error=False,fill_value='extrapolate')

Luvoir_A = {}
Luvoir_A['n_read']    = 0.008
Luvoir_A['cor_con']   = 1e-10
Luvoir_A['d_npix']    = 10
Luvoir_A['t_tps_f']   = L_tpsf
Luvoir_A['m_ra']      = 7.5
Luvoir_A['T_tel']     = 270
Luvoir_A['E_tel']     = 0.1
Luvoir_A['Starshade'] = False
Luvoir_A['ins_tp_f']  = LAf_ci
Luvoir_A['drn_f']     = L_drnf
Luvoir_A['ddc_f']     = L_ddcf

Luvoir_B = {}
Luvoir_B['n_read']    = 0.008
Luvoir_B['cor_con']   = 1e-10
Luvoir_B['d_npix']    = 10
Luvoir_B['t_tps_f']   = L_tpsf
Luvoir_B['m_ra']      = 4.0
Luvoir_B['T_tel']     = 270
Luvoir_B['E_tel']     = 0.1
Luvoir_B['Starshade'] = False
Luvoir_B['ins_tp_f']  = LBf_ci
Luvoir_B['drn_f']     = L_drnf
Luvoir_B['ddc_f']     = L_ddcf

HWO = {}
HWO['n_read']    = 0
HWO['cor_con']   = 1e-10
HWO['d_npix']    = 10
HWO['t_tps_f']   = L_tpsf
HWO['m_ra']      = 3.0
HWO['T_tel']     = 273
HWO['E_tel']     = 0.1
HWO['Starshade'] = False
HWO['ins_tp_f']  = LBf_ci
HWO['drn_f']     = sc.interpolate.interp1d([0,1], [0,0],kind='linear',bounds_error=False,fill_value='extrapolate')
HWO['ddc_f']     = sc.interpolate.interp1d([0,1], [3e-5,3e-5],kind='linear',bounds_error=False,fill_value='extrapolate')

HabEx_SS = {}
HabEx_SS['n_read']    = 0.008
HabEx_SS['cor_con']   = 1e-10
HabEx_SS['d_npix']    = 10
HabEx_SS['t_tps_f']   = HS_tpsf
HabEx_SS['m_ra']      = 2.0
HabEx_SS['T_tel']     = 270
HabEx_SS['E_tel']     = 0.1
HabEx_SS['ins_tp_f']  = HSf_ci
HabEx_SS['Starshade'] = True
HabEx_SS['drn_f']     = HS_drnf
HabEx_SS['ddc_f']     = HS_ddcf

HabEx_CG = {}
HabEx_CG['n_read']    = 0.008
HabEx_CG['cor_con']   = 2.5e-10
HabEx_CG['d_npix']    = 10
HabEx_CG['t_tps_f']   = HC_tpsf
HabEx_CG['m_ra']      = 2.0
HabEx_CG['T_tel']     = 270
HabEx_CG['E_tel']     = 0.1
HabEx_CG['ins_tp_f']  = HCf_ci
HabEx_CG['Starshade'] = False
HabEx_CG['drn_f']     = HC_drnf
HabEx_CG['ddc_f']     = HC_ddcf

AU_per_parsec = 3.0857E+13/1.4960E+08
Rstar = 695700
star_sr_pl(Rstar,1)
pc = 5


def noise_calculation_t(telescope,system,SN_tar_f,xmin,xmax,fmin,fmax,test=False,ret_noisem = False):
    success=False
    h   = 6.62607015e-34     # J s
    Na  = 6.02214076e23      # mol-1
    kb  = 1.380649e-23       # J/K
    c   = 299792458.0        # m/s
    amu = 1.660539040e-27    # kg
    G   = 6.67430e-11        # m3⋅kg−1⋅s−2
    c2  = h*c/kb

    
    # System parameters   
    xval     = system['xval']         # wavlength array
    IF       = system['I_F']          # I/F
    Rearth   = system['Rplan']        # radius planet km 
    Rstar    = system['Rstar']        # radias star km
    dsp      = system['dsp']          # distance star-planet in AU
    pc       = system['dist']         # distance to system in parsec
    Tstar    = system['Tstar']        # temperature of the star
    deg_arcs = 180*3600/np.pi         # degrees per arc second
    
    # observatory parameters
    #d_dc    = telescope['d_cd']                   # detector dark current
    n_read  = telescope['n_read']                  # electon read noise in e-
    cor_con = telescope['cor_con']                 # coronagraph contrast
    d_npix  = telescope['d_npix']                  # detector pixels
    tpsf    = telescope['t_tps_f']                 # telescope throughput function
    t_tp    = tpsf(xval[:-1]/1000)                 # telescope throughput * det_QE
    m_ra    = telescope['m_ra']                    # mirror radius
    t_ca    = np.pi*m_ra**2                        # telescope collecting area

    # For coronagraph throughput
    sss      = dsp/(pc*AU_per_parsec)*deg_arcs     # spatial separation of the system
    q_l      = xval[:-1]*1e-9                      # wavelengt (m)
    L_D      = 1.22 * q_l/(m_ra*2)*deg_arcs        # L/D in arcseconds

    # background noise
    nzodi    = 1   
    bg_f     = 23.3                                # background flux in V magnitude (550 nm)
    bg_jy    = 3640*2.512**-23.3                   # background flux in Janskys  (1 Jy = 10-26 W m-2 Hz-1)
    bg_phmh  = (bg_jy*1e-26)/(h*c/(550*1e-9))      # photon flux in photons/m2/hz
    bg_phm   = abs(bg_phmh*np.diff(c/(xval*1e-9))) # photon background flux per wvl bin per m2 
    bg_phm  /=  (1/L_D)**2                         # scale zodi flux by pixel size

    # telescope noise (not implemented yet)
    T_tel     = telescope['T_tel']                                # temperature of telescope
    E_tel     = telescope['E_tel']                                # emissivity of the telescope
    xnm       = xval*1e-9                           # xwavelength in nm
    ph_tel    = E_tel * (2e24*h*c**2/xnm)/(np.exp(c2/(T_tel*xnm)) -1)*1e-6  # photon noise telescope


    #print('Planet - star separation: (arcsec, L/D)',str(sss),int(sss/L_D),'\nCoronagraph planet throughput: ', str(LAf(sss/L_D)))
    #print()

    noisem = np.zeros(shape=(len(xval),2000))

    xmin_l = find_nearest_l(xval,xmin)
    xmax_l = find_nearest_l(xval,xmax)
    fmin_l = find_nearest_l(xval,fmin)
    fmax_l = find_nearest_l(xval,fmax)
    
    sr_of_pl = star_sr_pl(Rearth,pc*AU_per_parsec)
    sr_of_st = star_sr_pl(Rstar,pc*AU_per_parsec)

    ph_um_wvl   = (sr_of_pl*black_body_rad(Tstar,[i*1e-9 for i in xval])*star_sr_pl(Rstar,dsp)/np.pi)/[h*c/(i*1e-9) for i in xval]   # note the pi
    ph_um_wvl_s = (black_body_rad(Tstar,[i*1e-9 for i in xval])*star_sr_pl(Rstar,pc*AU_per_parsec)/np.pi)/[h*c/(i*1e-9) for i in xval] # note we do NOT need the factor pi here

    S           = ph_um_wvl[:-1]  * (np.diff(xval)/1000)  
    S_s         = ph_um_wvl_s[:-1]* (np.diff(xval)/1000)
    

    inst_tp_f =  telescope['ins_tp_f']       # instrument throughput function
    drnf      =  telescope['drn_f']          # detector read noise function
    ddcf      =  telescope['ddc_f']          # detector dark current function

    t_min = 60
    t_max = 200 * 3600
    for i,t_obs in enumerate(np.linspace(t_min,t_max,2000)):

        #Observational parameters
        n_exp   = t_obs/60                      # number of exposures
        t       = t_obs/60                      # time in seconds
        # source fluxes 
        if telescope['Starshade'] == False: 
            sig         = S      * t * n_exp * t_ca * t_tp *  inst_tp_f(sss/L_D)  * IF[:-1] 
        if telescope['Starshade'] == True: 
            sig         = S      * t * n_exp * t_ca * t_tp *  inst_tp_f(sss)  * IF[:-1] 
        stel_s      = S_s    * t * n_exp * t_ca * t_tp *  cor_con      
        zo_bg       = bg_phm * t * n_exp * t_ca * t_tp 
        tel_bg      = ph_tel[:-1] * t * n_exp * t_ca * t_tp 

        # noise 
        sourc_n = np.sqrt(sig+stel_s)

        # the factor (0.5*np.sqrt(2)) I cannot explain yet
        detc_n  = np.sqrt(d_npix*n_exp*(drnf(xval[:-1]/1000)**2+ddcf(xval[:-1]/1000)*t))/(0.5*np.sqrt(2))      

        #the factor x2 because it should be counted twice?
        bgzo_n  = np.sqrt(zo_bg)  * 2

        cmb_n   = np.sqrt( detc_n**2 + sourc_n**2 +  bgzo_n**2 )
        noisem[0,i]  = t_obs 
        noisem[1:,i] = sig/cmb_n

        #if test==True:
        #    return sig[xmin:xmax] ,cmb_n[xmin:xmax] , IF[xmin:xmax] ,noise_target[xmin:xmax]
        spec_no_feature = remove_feature(sig,fmin_l,fmax_l)[xmin_l:xmax_l]
        if detection_strength_feature(sig[xmin_l:xmax_l],spec_no_feature,cmb_n[xmin_l:xmax_l])>SN_tar_f: success = True; break
    if ret_noisem  == True:
        return noisem
    if test==True:
        plt.errorbar(xval[xmin_l:xmax_l],sig[xmin_l:xmax_l],cmb_n[xmin_l:xmax_l])
        plt.plot(xval[xmin_l:xmax_l],spec_no_feature)
    if success == True:
        return t_obs
    else:
        return 'SnR target not reached after 100 h'


spectra_d_wc = pd.read_csv('albedos.csv')

IFwc      = np.mean(spectra_d_wc.iloc[:,1:],axis=1)
xval      = spectra_d_wc.iloc[:,0]

# system parameters
system            = {}
system['xval']    = spectra_d_wc.iloc[:,0].values
system['I_F']     = IFwc
system['dist']    = 10
system['Rstar']   = 695700
system['Rplan']   = 12742/2
system['dsp']     = 1   
system['Tstar']   = 5777


IFwc      = spectra_d_wc.iloc[:,1]
xval      = spectra_d_wc.iloc[:,0]

xmin,xmax,fmin,fmax = 700,850,750,790
xmin_l = find_nearest_l(xval,xmin)
xmax_l = find_nearest_l(xval,xmax)
fmin_l = find_nearest_l(xval,fmin)
fmax_l = find_nearest_l(xval,fmax)



print(noise_calculation_t(HWO,system,5,xmin,xmax,fmin,fmax)/3600,' (should be 20.3)')
print(noise_calculation_t(Luvoir_A,system,5,xmin,xmax,fmin,fmax)/3600,' (should be 3.9)')
print(noise_calculation_t(Luvoir_B,system,5,xmin,xmax,fmin,fmax)/3600,' (should be 9.7)')

