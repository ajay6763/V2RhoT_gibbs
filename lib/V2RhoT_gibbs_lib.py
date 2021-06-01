import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from scipy import interpolate
import math


def lookup_P_T(V,P,table):
	index=[]
#	dist=np.array((T[:]-T_LitMod)**2-( P[:]-P_LitMod)**2)
	dist=np.array(((T-table[:,0])**2+(P-table[:,1])**2));
	index=dist.argmin(); 
	#print index, T_LitMod,P_LitMod
	return table[index,0]-273.15,table[index,2],table[index,3],table[index,4]


def lookup_Vp_P(vp,P,table):
	index=[]
#	dist=np.array((T[:]-T_LitMod)**2-( P[:]-P_LitMod)**2)
	dist=np.array(((vp-table[:,3])**2+(P-table[:,1])**2));
	index=dist.argmin(); 
	#print index, T_LitMod,P_LitMod
	return table[index,0]-273.0,table[index,2],table[index,3],table[index,4]

def lookup_vs_P(vs,P,table):
    index=[]
    #dist=np.array((T[:]-T_LitMod)**2-( P[:]-P_LitMod)**2)
    dist=np.array(((vs-table[:,4])**2+(P-table[:,1])**2));
    index=dist.argmin();
    '''
    if table[index,0] <= 1823.0:
        pass
    else:            
        index=0
        #print index, T_LitMod,P_LitMod
    '''
    return table[index,0]-273.0,table[index,2],table[index,3],table[index,4]


def atten_correction(T,P,Vp,Vs,oscill,grain_size):
    ## Parameters from Jackson and Faule 2010, Kumar et al., 2020
    A      = 816       #------------ Pre-exponential factor 
    alfa   = 0.36      #------------ frequency dependence
    energi = 293.0E03  #------------ Activation energy
    volexp = 1.20E-05  #------------ Activation volume
    R      = 8.314472  #------------ Gas constant
    pi     = 3.1415926 #------------ shephard;s pie :)
    
    #################################################
    ## calculating Qp and Qs
    parexp    = math.exp((-(energi+(volexp*P)))/(R*(T)))
    sqatt50   = A*(((oscill*(1.0E0/(grain_size*1000.0E0)))*parexp))**alfa
    Qp        = (1/sqatt50)*(9/4)
    Qs        = 1/sqatt50
    cots50    = ((1.0E0/math.tan((pi*alfa)/2.0E0))*sqatt50)*0.5E0
    cotp50    = ((1.0E0/math.tan((pi*alfa)/2.0E0))*sqatt50)*(2.0E0/9.0E0)
    #################################################
    ## correcting velocities
    Vs_correc = Vs*(1.0E0-cots50)
    Vp_correc = Vp*(1.0E0-cotp50)
    return Vp_correc,Vs_correc



# Defining function to correct for reduction in veloctiy from melts
# Solidus and liquidus temperature for crust
## Temperature is in degree celcius and Pressure is in GPa
''''
# For crust
  T_s_c =  920 + 156*P ---------------------------------------- (1)
  T_l_c = 1120 + 156*P -----------------------------------------(2)
 where T_s_c is the solidus of crustal rocks P is pressure in GPa
 where T_l_c is liquidus of crustal rocks P is pressure in GPa

# For mantle
  T_s_m = 1080 + 134.2*P - 6.581*P*P + 0.1054*P*P*P ----------- (3)
  T_l_m = 1762 + 57-46*P - 3.48*P*P + 0.077*P*P*P   ----------- (4)
  where T_s_m is solidus of mantle rocks and T_l_m is liquidus of mantle rocks

 Solidus temperatures here represents dy granitic rocks and dry peridotites (Hirschmann, 2000 and Winter, 2010),
 hence partial melts predicted here should be taken as indicative only.
# Now, melt fractions can be computed as following using eq. 3 -4
 M_crust = (T-T_s_c)/(T_l_c -T_s_c) ------------------------------------------- (5)
 M_mantle = (T-T_s_m)/(T_l_m -T_s_m) ------------------------------------------- (6)
 where M_crust and M_mantle is melt fraction for the crust and mantle respectively, T is the actual temperature in  the model

# Effects of melts on density and seismic velocitie::
 Text from Afonso et al., 2016 III
 inversion.We neglect the effect of melt on bulk density reduction, as it is always small compared to the uncertainties in data sets 
 constraining this property. The effects of melts on seismic velocities, on the other hand, are significant and we estimate them based 
 on the results of Hammond and Humphreys [2000]. In the absence of more detailed information, we chose the average of the minimum and maximum 
 values proposed by these authors for Vs and Vp velocity reductions for melt fractions <= 1% : (-5.3 = dlnVs/%melt and -2.4 = dlnVp/%melt). 
 Over 99% of our predicted melt fractions are below 1% (see section 5), so these derivatives are considered adequate.
'''
## crust solidus
def crust_T_s_c(T,P):
	T_s_c = 920 + 156*P
	return T_s_c

## crust liquidus
def crust_T_l_c(T,P):
	T_l_c = 1120 + 156*P
	return T_l_c

## mantle solidus
def mantle_T_s_m(T,P):
	T_s_m = 1080 + 134.2*P - 6.581*P*P + 0.1054*P*P*P
	return T_s_m

## mantle liquidus
def mantle_T_l_m(T,P):
	T_l_m = 1762 + 57.46*P - 3.48*P*P + 0.077*P*P*P
	return T_l_m

## melt fraction - general
def melt_frac(T,T_l,T_s):
	melt = (T - T_s)/(T_l - T_s)
	return melt

## melt fraction - crust
def melt_frac_crust(T,P):
	T_s_c = 920 + 156*P
	T_l_c = 1120 + 156*P
	melt = (T - T_s_c)/(T_l_c - T_s_c)
	return melt

## melt fraction - mantle
def melt_frac_mantle(T,P):
	T_s_m = 1080.0 + 134.20*P - 6.581*P**2 + 0.1054*P**3
	T_l_m = 1762.0 + 57.46*P  - 3.487*P**2 + 0.077*P**3
	melt  = (T - T_s_m)/(T_l_m - T_s_m)
	return melt

def velocity_melt_correction(T,P,Vp,Vs):
    ## get the melt fraction at the P and T
    melt_frac = melt_frac_mantle(T,P)
    #########################
    # Vp:   -5.3 = dlnVs/%melt => dlnVs = -5.3 * %melt 
    # Vs:   -2.4 = dlnVp/%melt => dlnVp = -2.3 * %melt
    '''
    dVp = np.exp(-5.3 * melt_frac) 
    dVs = np.exp(-2.3 * melt_frac)
    Vp_corrected = Vp - dVp
    Vs_corrected = Vs - dVs
    '''
    if melt_frac > 0:
        dVp = np.exp(-5.3 * melt_frac)
        dVs = np.exp(-2.3 * melt_frac)
        Vp_corrected = Vp - dVp
        Vs_corrected = Vs - dVs
    else:
        Vp_corrected = Vp
        Vs_corrected = Vs
    
    return Vp_corrected,Vs_corrected, melt_frac


def vel_to_temp(depth,velocity,Table,index):
    Temperature_out = np.zeros_like(depth[:])
    Density_out     = np.zeros_like(depth[:])
    Vp_out          = np.zeros_like(depth[:])
    Vs_out          = np.zeros_like(depth[:])
    for i in range(len(index)):
        P  = pressure_inter(depth[index[i]])
        Vs_in = velocity[index[i]]
        temp,dens,vp,vs=lookup_vs_P(Vs_in,P.tolist(),Table)
        Vp_out[index[i]]=vp
        Vs_out[index[i]]=vs
        Temperature_out[index[i]]=temp
        Density_out[index[i]]=dens
    return Temperature_out,Density_out,Vp_out,Vs_out