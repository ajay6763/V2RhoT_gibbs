import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from scipy import interpolate
import math
import scipy as scipy


########################################
# Defining pressure function based on Ak133
# I choose ak135 becaouse the pressures from the LitMod2D_2.0 ref model
# are practically similar to ak135
########################################
try:
    ak135 = np.loadtxt('./databases/ak135f.txt',skiprows=1)
except:
    ak135 = np.loadtxt('/home/ajay/projects/V2RhoT_gibbs/databases/ak135f.txt',skiprows=1)

ak135_P = 9.8*ak135[:,0]*1e3*ak135[:,1]*1e3*1e-5
#pressure_inter = interpolate.interp1d(ak135[:,0],ak135_P)
#depth_inter = interpolate.interp1d(ak135_P,ak135[:,0])
y2=ak135[:,1]*9.8
gravity=9.8
density=ak135[:,1]*1e3
depth=ak135[:,0]*1e3
pressure_int = scipy.integrate.cumulative_trapezoid(gravity * density, x=depth, initial=0)
pressure_int=pressure_int*1e-5 # pascal to bar
pressure_inter = interpolate.interp1d(ak135[:,0],pressure_int)


## melt fraction - mantle
def melt_frac_mantle_Hirschmann(T,P):
    """
    Hirschmann, 2000 and Winter, 2010
    P in GPa
    T in oC
    """
    T_s_m = 1080.0 + 134.2*P - 6.581*P**2 + 0.1054*P**3
    T_l_m = 1762.0 + 57.46*P  - 3.487*P**2 + 0.077*P**3
    melt  = (T - T_s_m)/(T_l_m - T_s_m)
    return melt
def melt_frac_mantle_Hirschmann(T, P):
    """
    Hirschmann, 2000 and Winter, 2010
    P in GPa, T in °C
    Supports both scalar values and NumPy arrays.
    """
    T_s_m = 1080.0 + 134.2*P - 6.581*P**2 + 0.1054*P**3
    T_l_m = 1762.0 + 57.46*P  - 3.487*P**2 + 0.077*P**3
    
    # Avoid potential division by zero if T_l_m equals T_s_m
    denominator = np.where(T_l_m == T_s_m, 1e-5, T_l_m - T_s_m)
    melt = (T - T_s_m) / denominator
    # Restrict bounds cleanly to [0.0, 1.0]
    return np.clip(melt, 0.0, 1.0)
def velocity_melt_correction_mantle_Hammond_Humphreys(T, P, Vp, Vs):
    """
    Hammond_Humphreys¸
    """
    ## get the melt fraction at the P and T
    melt_frac = melt_frac_mantle_Hirschmann(T,P)
    #########################
    # Vs:   -5.3 = dlnVs/melt_frac => dlnVs = -5.3 * melt_frac
    # Vp:   -2.4 = dlnVp/melt_frac => dlnVp = -2.4 * melt_frac
    if melt_frac > 0:
        dVs = -5.3 * melt_frac
        dVp = -2.4 * melt_frac
        Vs_corrected = Vs*np.exp(dVs)
        Vp_corrected = Vp*np.exp(dVp)
    else:
        Vp_corrected = Vp
        Vs_corrected = Vs
        melt_frac    = 0.
    return Vp_corrected,Vs_corrected, melt_frac*100.

def atten_correction_JF2010(T,P,Vp,Vs,oscill,grain_size,A=68,alpha=0.36,energy=293e3,volexp=1.20e-5):
    '''
    Jackson and Faul, 2010 and Kumar et al., 2020    
    Input:
    T - Kelvin
    P - Pascal
    Vp,Vs - km/s 
    oscill - time period (seconds)
    d - mm
    Output:
    Vp,Vs: km/s    
    '''
    ## Parameters from Jackson and Faule 2010, Kumar et al., 2020
    #A = 68          # 𝑠^(−𝛼) 𝑚m^(−𝛼) Pre-exponential factor (for d in mm)
    #alpha   = 0.34      #------------ frequency dependence
    #energy = 293e3  #------------ Activation energy
    #volexp = 1.20e-5  #------------ Activation volume
    R      = 8.314472  #------------ Gas constant
    pi     = 3.1415926 #------------ shephard;s pie :)

    #################################################
    ## calculating Qp and Qs
    parexp    = math.exp((-(energy+(volexp*P)))/(R*(T)))
    Qs_inv   = A*(((oscill*(1.0/(grain_size)))*parexp))**alpha
    Qp        = (1/Qs_inv)*(9/4)
    Qs        = 1/Qs_inv
    vs_correction    = ((1.0/math.tan((pi*alpha)/2.0))*Qs_inv)*0.5
    # assuming Qp=9/4*Qs, then Qp^-1=4/9*Qs^-1 (i.e., bulk attenuation is negligible)
    vp_correction    = ((1.0/math.tan((pi*alpha)/2.0))*Qs_inv)*(2.0/9.0)
    #################################################
    ## correcting velocities
    Vs_correc = Vs*(1.0-vs_correction)
    Vp_correc = Vp*(1.0-vp_correction)
    return Vp_correc,Vs_correc
def atten_correction_Behn2009(T,P,Vp,Vs,oscill,d,COH,rQ=1.2):
    '''
    Behn et al., 2009 https://doi.org/10.1016/j.epsl.2009.03.014
    Input:
    T - Kelvin
    P - Pascal
    Vp,Vs - km/s 
    oscill - time period (seconds)
    d - meter
    COH - olivine water concentration (in H/10**6Si), e.g., 50 H/10**6Si --> dry; 1000 H/10**6Si --> wet equivalent to 125 +/- 75 ppm weight water cf Behn et al., 2009
    rQ - default (1.2). to include effect of hydration on Qs^-1, rQ=1.2 for wet olivine and rQ=0 for dry olivine
    -- For diffusion creep it is 1.0 and for dislocation creep it is 1.2 cf. Behn et al., 2009 Table S1.
    Note: For anelastic attenuation pq=1 is used i.e., processes are at the grain-boundary so not really diffusion.
    for the dislocation grain size exponent is =0 then rQ=1.2 is apt. But it is assumed that the effect of hydration in anelasticity scales similar to dislocation creep (wistful).
    Parameters are also taken from https://github.com/wshinevar/WISTFUL/blob/main/behn2009Shear.m
    Output:
    '''
    frequency   =   1/oscill # CHECK THIS FOR OMEGA
    R           =   8.314    # gas constant
    pi          =   3.1415926 #------------ shephard;s pie :)
    pq_ref      =   1.09    #reference grain exponent 
    pq          =   1.0     # grain exponent 
    TQ_ref      =   1265+273.15 # refrence temperature in Kelvin
    d_ref       =   1.24e-5 # reference grain size in meters 
    EQ_ref      =   505e3   # referene activation energy in J/mol
    EQ          =   420e3   #activation energy
    VQ_ref      =   1.2e-5  #reference activation volime m**3/mol 
    VQ          =   1.2e-5  #activation volume
    Bo          =   1.28e8  # prefactor for Q for omega=0.122 s^-1
    PQ_ref      =   300e6 # reference pressure in Pa;
    #Hydration effect
    COH_ref     =   50      # H/10^6 Si
    # CR. for rQ: Shinevar et al., 2014, https://doi.org/10.1029/2022GC010329
    alpha       =   0.27
    # Calulate the pre-exponentail factor B based on the given parameters
    B = Bo*d_ref**(pq-pq_ref)*(COH/COH_ref)**rQ*math.exp(((EQ+PQ_ref*VQ)-(EQ_ref+PQ_ref*VQ_ref))/(R*TQ_ref))
    # Inverse of Anelastic factor Qs^-1
    Qs_inv = (B * d**(-1 * pq) / frequency * math.exp(-(EQ + P * VQ) / (R * T)))**alpha 
    Qp        = (1/Qs_inv)*(9/4)
    Qs        = 1/Qs_inv
    vs_correction    = ((1.0/math.tan((pi*alpha)/2.0))*Qs_inv)*0.5
    # assuming Qp=9/4*Qs, then Qp^-1=4/9*Qs^-1 (i.e., bulk attenuation is negligible)
    vp_correction    = ((1.0/math.tan((pi*alpha)/2.0))*Qs_inv)*(2.0/9.0)
    #################################################
    ## correcting velocities
    Vs_correc = Vs*(1.0-vs_correction)
    Vp_correc = Vp*(1.0-vp_correction)
    return Vp_correc,Vs_correc

def lookup_vs_P_accurate(vs_obs, P_obs, table, P_tol=None):
    """
    Lookup/interpolate T and properties from a Perple_X P-T table
    using an isobaric slice and 1D Vs inversion.
    """
    # 1. Isobaric Lock: Find closest pressure level
    P_values = np.unique(table[:, 1])
    P_closest = P_values[np.argmin(np.abs(P_values - P_obs))]

    # Extract and sort slice by Temperature
    mask = np.isclose(table[:, 1], P_closest)
    sub = table[mask]
    sub = sub[np.argsort(sub[:, 0])]

    T_grid = sub[:, 0]
    Vs_grid = sub[:, 4]

    # Evaluate matching conditions
    exact = np.isclose(Vs_grid, vs_obs)
    diff = Vs_grid - vs_obs
    crossing = np.where(diff[:-1] * diff[1:] <= 0)[0]

    # 2. Branch selection to compute property outputs
    if np.any(exact):
        i = np.where(exact)[0][0]
        P_out = sub[i, 1]
        T_celsius = sub[i, 0] - 273.15
        Dens = sub[i, 2]
        Vp = sub[i, 3]
        Vs_out = sub[i, 4]
        melt = sub[i, 5]

    elif len(crossing) > 0:
        i = crossing[0]
        T1, T2 = T_grid[i], T_grid[i+1]
        Vs1, Vs2 = Vs_grid[i], Vs_grid[i+1]

        denom = (Vs2 - Vs1) if Vs2 != Vs1 else 1.0
        w = (vs_obs - Vs1) / denom

        P_out = P_closest
        T_celsius = (T1 + w * (T2 - T1)) - 273.15
        Dens = sub[i, 2] + w * (sub[i+1, 2] - sub[i, 2])
        Vp = sub[i, 3] + w * (sub[i+1, 3] - sub[i, 3])
        Vs_out = vs_obs
        melt = sub[i, 5] + w * (sub[i+1, 5] - sub[i, 5])

    else:
        # Out-of-bounds fallback (clamps to nearest edge)
        i = np.argmin(np.abs(Vs_grid - vs_obs))
        P_out = sub[i, 1]
        T_celsius = sub[i, 0] - 273.15
        Dens = sub[i, 2]
        Vp = sub[i, 3]
        Vs_out = sub[i, 4]
        melt = sub[i, 5]

    return P_out, T_celsius, Dens, Vp, Vs_out, melt

def vel_vs_to_temp_prop_out(depth,Vs,Table):
    '''
    Input:
    depth : depth column in km.
    Vs    : tomography Vs velocity in km/s.
    Table : Perplex lookup table corrected form anelasticity and melt effects.

    Output: Output: [depth,P_out,Temperature_out,Density_out,Vp_out,Vs_out,diff_Vp,melt_out)
    '''
    Temperature_out = []#np.zeros_like(tomo[:,1])
    Density_out     = []#np.zeros_like(tomo[:,1])
    melt_out     = [] #np.zeros_like(tomo[:,1])
    Vp_out     = [] #np.zeros_like(tomo[:,1])
    Vs_out     = [] #np.zeros_like(tomo[:,1])
    diff_Vs         = []
    P_out           = []
    #Vp_out          = []#np.zeros_like(tomo[:,1])
    #Vs_out          = []#np.zeros_like(tomo[:,1])
    for i in range(len(depth)):
        P  = pressure_inter(depth[i])
        Vs_in = Vs[i]
        P_table,temp,dens,vp,vs,m=lookup_vs_P_accurate(Vs_in,P.tolist(),Table)
        #Vp_out.append(vp)
        #Vs_out.append(vs)
        P_out.append(P_table)
        Temperature_out.append(temp)
        Density_out.append(dens)
        Vs_out.append(vs)
        Vp_out.append(vp)
        #diff_Vs.append(((Vs_in-vs)/Vs_in)*100)
        diff_Vs.append(((Vs_in - vs) / Vs_in * 100.0) if not np.isclose(Vs_in, 0.0) else 0.0)
        melt_out.append(m)
    ### pasting the outputs to the input tomo table
    out=depth;
    out=np.column_stack((out,P_out))
    out=np.column_stack((out,Temperature_out))
    out=np.column_stack((out,Density_out))
    out=np.column_stack((out,Vp_out))
    out=np.column_stack((out,Vs_out))
    out=np.column_stack((out,diff_Vs))
    out=np.column_stack((out,melt_out))
    return out


def mantle_melt_atten_correction_JF2010(Table,grain_size,oscillation):
    Table_atten_corrected = np.copy(Table)
    #for i in range(len(Table_atten_corrected)):
    #    Table_atten_corrected[i,3],Table_atten_corrected[i,4] = atten_correction_J_2002(Table_atten_corrected[i,0],Table_atten_corrected[i,1]*1e5,
    #                                                         Table_atten_corrected[i,3],Table_atten_corrected[i,4],oscillation,grain_size)
    for i in range(len(Table_atten_corrected)):
        Table_atten_corrected[i,3],Table_atten_corrected[i,4] = atten_correction_JF2010(Table_atten_corrected[i,0],Table_atten_corrected[i,1]*1e5,
                                                             Table_atten_corrected[i,3],Table_atten_corrected[i,4],oscillation,grain_size)
    Table_atten_melt_corrected = np.copy(Table_atten_corrected)
    melt = np.zeros_like(Table_atten_melt_corrected[:,0])
    for i in range(len(Table_atten_melt_corrected)):
        Table_atten_melt_corrected[i,3],Table_atten_melt_corrected[i,4],melt[i] = velocity_melt_correction_mantle_Hammond_Humphreys(Table_atten_melt_corrected[i,0]-273.15,
                                                                                                 Table_atten_melt_corrected[i,1]/1e4,
                                                                 Table_atten_melt_corrected[i,3],Table_atten_melt_corrected[i,4])
    Table_atten_melt_corrected[:,5]=melt[:]
    return Table_atten_melt_corrected

def mantle_melt_atten_correction_Behn2009(Table,grain_size,oscillation,COH):
    Table_atten_corrected = np.copy(Table)
    for i in range(len(Table_atten_corrected)):
        Table_atten_corrected[i,3],Table_atten_corrected[i,4] = atten_correction_Behn2009(Table_atten_corrected[i,0],Table_atten_corrected[i,1]*1e5,
                                                             Table_atten_corrected[i,3],Table_atten_corrected[i,4],oscillation,grain_size/1e3,COH)
    Table_atten_melt_corrected = np.copy(Table_atten_corrected)
    melt = np.zeros_like(Table_atten_melt_corrected[:,0])
    for i in range(len(Table_atten_melt_corrected)):
        Table_atten_melt_corrected[i,3],Table_atten_melt_corrected[i,4],melt[i] = velocity_melt_correction_mantle_Hammond_Humphreys(Table_atten_melt_corrected[i,0]-273.15,
                                                                                                 Table_atten_melt_corrected[i,1]/1e4,
                                                                 Table_atten_melt_corrected[i,3],Table_atten_melt_corrected[i,4])
    Table_atten_melt_corrected[:,5]=melt[:]
    return Table_atten_melt_corrected