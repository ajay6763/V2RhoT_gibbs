import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from scipy import interpolate
import math

########################################
# Defining pressure function based on Ak133
# I choose ak135 becaouse the pressures from the LitMod2D_2.0 ref model
# are practically similar to ak135
########################################
ak135 = np.loadtxt('./databases/ak135f.txt',skiprows=1)
ak135_P = 9.8*ak135[:,0]*1e3*ak135[:,1]*1e3*1e-5
pressure_inter = interpolate.interp1d(ak135[:,0],ak135_P)
depth_inter = interpolate.interp1d(ak135_P,ak135[:,0])

def lithostatic_pressure(depth,density):
    Pressure_lith = np.zeros_like(depth)
    # Check is the first index is at zero depth or not
    # If it is not then pressures is set equal rho*g*h
    # else is it set to 0
    if depth[0] != 0:
        Pressure_lith[0] = density[0] * depth[0]* 9.8 * 1e3 * 1e-5
    else:
        #thickness = (depth[1]-depth[0]) #*1e3
        Pressure_lith[0] = 0
    # Now looping through depths
    for j in range(len(depth)-2):
        thickness = (depth[j+2] - depth [j+1]) *1e3
        Pressure_lith[j+1]=Pressure_lith[j] + thickness*density[j+1] *9.8*1e-5

    # fixing last index
    Pressure_lith[-1] = Pressure_lith[-2] + (depth[-1] - depth [-2])*1e3 * density[-2]*9.8*1e-5
    #print(depth[:],density[:],Pressure_lith[:])

    p_func = interpolate.interp1d(depth,Pressure_lith)
    return p_func
def lookup_P_T(V,P,table):
	index=[]
    #dist=np.array((T[:]-T_LitMod)**2-( P[:]-P_LitMod)**2)
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

def lookup_vs_P_accurate(vs,P,table):
    """
    Input:

    Output:

    
    This function is a bit accurate than the minimum of the L2 norm.
    So, what I am doing is that first I look for the minimum of the L2 norm, then
    I look for the difference between the observed velocity and node above and below.
    In case if the L2 norm give "bulls eye" hit where observed velocity matches the
    node velocity I pick the properties from that node. If not then I ask which way,
    up of down, difference between the observed and node velocity is minimum and
    take the average of the properties at the minimum L2 norm node and up or down node.
    """
    index=[]
    Vp=[]
    Vs=[]
    Dens=[]
    T=[]
    #dist=np.array((T[:]-T_LitMod)**2-( P[:]-P_LitMod)**2)
    dist=np.array(((vs-table[:,4])**2+(P-table[:,1])**2));
    index=dist.argmin();

    diff_vs=table[index,4] - vs
    diff_vs_up=table[index-1,4] - vs
    diff_vs_down=table[index+1,4] - vs

    if diff_vs==0:
        Vp=table[index,3]
        Vs=table[index,4]
        Dens=table[index,2]
        T=table[index,0]-273.0
    elif diff_vs_up<diff_vs_down:
        Dens=(table[index,2]+table[index-1,2])/2
        Vp=(table[index,3]+table[index-1,3])/2
        Vs=(table[index,4]+table[index-1,4])/2
        T=-273.0+(table[index,0]+table[index-1,0])/2
    else:
        Dens=(table[index,2]+table[index+1,2])/2
        Vp=(table[index,3]+table[index+1,3])/2
        Vs=(table[index,4]+table[index+1,4])/2
        T=-273.0+(table[index,0]+table[index+1,0])/2

        #print index, T_LitMod,P_LitMod
    return table[index,1],T,Dens,Vp,Vs

def atten_correction(T,P,Vp,Vs,oscill,grain_size):
    """
    """
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

def atten_correction_Behn2009(T,P,Vp,Vs,oscill,d,COH):
    """
    Behn et al., 2009 https://doi.org/10.1016/j.epsl.2009.03.014
    Input:
    T - Kelvin
    P - Pascal
    Vp,Vs - km/s 
    oscill - time period (seconds)
    d - meter
    COH
    Output:
    """
    frequency   =   1/oscill # CHECK THIS FOR OMEGA
    R           =   8.314    # gas constant
    pi          =   3.1415926 #------------ shephard;s pie :)
    pq_ref      =   1.09    #reference grain exponent 
    pq          =   1.0     # grain exponent 
    TQ_ref      =   1265    # refrence temperature in oC
    d_ref       =   1.24e-5 # reference grain size in meters 

    EQ_ref      =   505e3   # referene activation energy in J/mol
    EQ          =   420e3   #activation energy
    VQ_ref      =   1.2e-5  #reference activation volime m**3/mol 
    VQ          =   1.2e-5  #activation volume
    
    Bo          =   1.28e8  # prefactor for Q for omega=0.122 s^-1
    COH_ref     =   50      # H/10^6 Si
    PQ_ref      =   300e6 # reference pressure in Pa;
    rQ          =   1.2;
    alpha       =   0.27;

    B = Bo*d_ref**(pq-pq_ref)*(COH/COH_ref)**rQ*math.exp(((EQ+PQ_ref*VQ)-(EQ_ref+PQ_ref*VQ_ref))/(R*TQ_ref))
    
    Qs_inv=(B*d**(-1*pq)*frequency**(-1)*math.exp(-(EQ+P*VQ)/(R*T)))**alpha # Inverse of Anelastic factor
    Qp        = (1/Qs_inv)*(9/4)
    Qs        = 1/Qs_inv
    vs_correction    = ((1.0/math.tan((pi*alpha)/2.0))*Qs_inv)*0.5
    vp_correction    = ((1.0/math.tan((pi*alpha)/2.0))*Qs_inv)*(2.0/9.0)
    #################################################
    ## correcting velocities
    Vs_correc = Vs*(1.0-vs_correction)
    Vp_correc = Vp*(1.0-vp_correction)
    return Vp_correc,Vs_correc


def atten_correction_Behn2009_crust(T,P,Vp,Vs,oscill,d,COH):
    """
    Behn et al., 2009 https://doi.org/10.1016/j.epsl.2009.03.014
    Input:
    T - Kelvin
    P - Pascal
    Vp,Vs - km/s 
    oscill - time period (seconds)
    d - meter
    COH
    Output:
    """
    frequency   =   1/oscill # CHECK THIS FOR OMEGA
    R           =   8.314    # gas constant
    pi          =   3.1415926 #------------ shephard;s pie :)
    pq_ref      =   1.09    #reference grain exponent 
    pq          =   1.0     # grain exponent 
    TQ_ref      =   1265    # refrence temperature in oC
    d_ref       =   1.24e-5 # reference grain size in meters 

    EQ_ref      =   505e3   # referene activation energy in J/mol
    EQ          =   420e3   #activation energy
    VQ_ref      =   1.2e-5  #reference activation volime m**3/mol 
    VQ          =   1.2e-5  #activation volume
    
    Bo          =   1.28e8  # prefactor for Q for omega=0.122 s^-1
    COH_ref     =   50      # H/10^6 Si
    PQ_ref      =   300e6 # reference pressure in Pa;
    rQ          =   1.2;
    alpha       =   0.27;

    B = Bo*d_ref**(pq-pq_ref)*(COH/COH_ref)**rQ*math.exp(((EQ+PQ_ref*VQ)-(EQ_ref+PQ_ref*VQ_ref))/(R*TQ_ref))
    
    Qs_inv=(B*d**(-1*pq)*frequency**(-1)*math.exp(-(EQ+P*VQ)/(R*T)))**alpha # Inverse of Anelastic factor
    Qp        = (1/Qs_inv)*(9/4)
    Qs        = 1/Qs_inv
    if T > 800 +273.15 and P > 1*1e9:
        vs_correction    = ((1.0/math.tan((pi*alpha)/2.0))*Qs_inv)*0.5
        vp_correction    = ((1.0/math.tan((pi*alpha)/2.0))*Qs_inv)*(2.0/9.0)
    else:
        vs_correction    = ((1.0/math.tan((pi*alpha)/2.0))*(1/50))*0.5
        vp_correction    = ((1.0/math.tan((pi*alpha)/2.0))*(1/100))*(2.0/9.0)
    #################################################
    ## correcting velocities
    Vs_correc = Vs*(1.0-vs_correction)
    Vp_correc = Vp*(1.0-vp_correction)
    return Vp_correc,Vs_correc

def atten_correction_J_2002(T,P,Vp,Vs,oscill,grain_size):
    """
    """
    ## Parameters from Jackson et al. 2002
    A      = 750       #------------ Pre-exponential factor
    alfa   = 0.26      #------------ frequency dependence
    energi = 424.0E03  #------------ Activation energy
    volexp = 1.60E-05  #------------ Activation volume
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
    """
    Hirschmann, 2000 and Winter, 2010
    """
    T_s_c = 920 + 156*P
    T_l_c = 1120 + 156*P
    melt = (T - T_s_c)/(T_l_c - T_s_c)
    return melt

## melt fraction - mantle
def melt_frac_mantle_Hirschmann(T,P):
    """
    Hirschmann, 2000 and Winter, 2010
    """
    T_s_m = 1080.0 + 134.2*P - 6.581*P**2 + 0.1054*P**3
    T_l_m = 1762.0 + 57.46*P  - 3.487*P**2 + 0.077*P**3
    melt  = (T - T_s_m)/(T_l_m - T_s_m)
    return melt

## melt fraction - mantle
def melt_frac_mantle_Katz(T,P):
    """
    Katz et al. 2003
    T_solidus = A1 + A2*P + A3*P**2
    T_liquidus_lherzo = B1 + B2*P + B3*P**2
    T_liquidus = C1 + C2*P + C3*P**2, this more general for arguments see Katz et al. 2003. Hence I am using this one.
    This basically means using different coefficients from Table 2 of Katz et al. 2003
    """
    T_s_m = 1085.7 + 132.9*P - 5.1*P**2
    T_l_m = 1780.0 + 45.0*P  - 2.0*P**2
    melt  = (T - T_s_m)/(T_l_m - T_s_m)
    return melt

def velocity_melt_correction_mantle_Hammond_Humphreys(T,P,Vp,Vs):
    """
    Hammond_Humphreys
    """
    ## get the melt fraction at the P and T
    melt_frac = melt_frac_mantle_Hirschmann(T,P)
    #########################
    # Vs:   -5.3 = dlnVs/%melt => dlnVs = -5.3 * %melt
    # Vp:   -2.4 = dlnVp/%melt => dlnVp = -2.3 * %melt
    '''
    dVs = np.exp(-5.3 * melt_frac)
    dVp = np.exp(-2.3 * melt_frac)
    Vp_corrected = Vp - dVp
    Vs_corrected = Vs - dVs
    '''
    if melt_frac > 0:
        dVs = -5.3 * melt_frac #*100
        dVp = -2.3 * melt_frac #*100
        Vp_corrected = Vp + (dVp*Vp)/100
        Vs_corrected = Vs + (dVs*Vs)/100
    else:
        Vp_corrected = Vp
        Vs_corrected = Vs
        melt_frac    = 0
    return Vp_corrected,Vs_corrected, melt_frac

def velocity_melt_correction_mantle_Chantel_2016(T,P,Vp,Vs):
    """
    From Cahntle et al 2016 Sci. Adv.
    Reduction in Vp and Vs as a functoion of melts
    The curves in this paper have absolute value of Vp and Vs
    as intercept and correction are added to it.
    What I will do here is that replace this absolute value
    with the anharmonic velocity. This means that I am simply moving
    correction up or down (scaling and up-scaling?) for all
    pressure ranges i.e. depth levels.
    """
    ## get the melt fraction at the P and T
    melt_frac = melt_frac_mantle_Katz(T,P)
    #########################
    if melt_frac > 0:
        dVs = -5.3 * melt_frac #*100
        dVp = -2.3 * melt_frac #*100
        Vp_corrected = 0.07*melt_frac**2  - 0.5566*melt_frac + Vp
        Vs_corrected = 0.065*melt_frac**2 - 0.5565*melt_frac + Vp
    else:
        Vp_corrected = Vp
        Vs_corrected = Vs
        melt_frac    = 0
    return Vp_corrected,Vs_corrected, melt_frac



def velocity_melt_correction_crust(T,P,Vp,Vs):
    ## get the melt fraction at the P and T
    melt_frac = melt_frac_crust(T,P)
    #########################
    # Vs:   -5.3 = dlnVs/%melt => dlnVs = -5.3 * %melt
    # Vp:   -2.4 = dlnVp/%melt => dlnVp = -2.3 * %melt
    '''
    dVs = np.exp(-5.3 * melt_frac)
    dVp = np.exp(-2.3 * melt_frac)
    Vp_corrected = Vp - dVp
    Vs_corrected = Vs - dVs
    '''
    if melt_frac > 0:
        dVs = -5.3 * melt_frac #*100
        dVp = -2.3 * melt_frac #*100
        Vp_corrected = Vp + (dVp*Vp)/100
        Vs_corrected = Vs + (dVs*Vs)/100
    else:
        Vp_corrected = Vp
        Vs_corrected = Vs
        melt_frac    = 0.
    return Vp_corrected,Vs_corrected, melt_frac

def lookup_vs_P_accurate_prop(vs,P,table):
    """
    
    This function is a bit accurate than the minimum of the L2 norm.
    So, what I am doing is that first I look for the minimum of the L2 norm, then
    I look for the difference between the observed velocity and node above and below.
    In case if the L2 norm give "bulls eye" hit where observed velocity matches the
    node velocity I pick the properties from that node. If not then I ask which way,
    up of down, difference between the observed and node velocity is minimum and
    take the average of the properties at the minimum L2 norm node and up or down node.
    """
    index=[]
    Vp=[]
    Vs=[]
    Dens=[]
    T=[]
    P_out=[]
    melt=[]
    #dist=np.array((T[:]-T_LitMod)**2-( P[:]-P_LitMod)**2)
    dist=np.array(((vs-table[:,4])**2+(P-table[:,1])**2)**0.5);
    index=dist.argmin();

    diff_vs=table[index,4] - vs
    diff_vs_up=table[index-1,4] - vs
    diff_vs_down=table[index+1,4] - vs

    if diff_vs==0:
        T=table[index,0]-273.0
        P_out=table[index,1]
        Dens=table[index,2]
        Vp=table[index,3]
        Vs=table[index,4]
        melt=table[index,5]
    elif diff_vs_up<diff_vs_down:
        T=-273.0+(table[index,0]+table[index-1,0])/2
        P_out=(table[index,1]+table[index-1,1])/2
        Dens=(table[index,2]+table[index-1,2])/2
        Vp=(table[index,3]+table[index-1,3])/2
        Vs=(table[index,4]+table[index-1,4])/2
        melt=table[index,5]
        melt=(table[index,5]+table[index-1,5])/2

    else:
        T=-273.0+(table[index,0]+table[index+1,0])/2
        P_out=(table[index,1]+table[index+1,1])/2
        Dens=(table[index,2]+table[index+1,2])/2
        Vp=(table[index,3]+table[index+1,3])/2
        Vs=(table[index,4]+table[index+1,4])/2
        melt=table[index,5]
        melt=(table[index,5]+table[index+1,5])/2

        #print index, T_LitMod,P_LitMod
    return P_out,T,Dens,Vp,Vs,melt

def mantle_melt_atten_correction(Table,grain_size,oscillation):
    """
    Table : perplex table
    grain_size : grain size in mm.
    oscilation: oscillation period in seconds.
    """
    # correction using grain size = 10 mm and oscillatio period of 75 seconds.
    # Attenuation model of Jackson and Faul 2010
    # Function: lib.atten_correction (T (oC),P (Pascal),Vp (km/s),Vs (km/s),oscilation period (s), grain size (mm))
    Table_atten_corrected = np.copy(Table)
    #for i in range(len(Table_atten_corrected)):
    #    Table_atten_corrected[i,3],Table_atten_corrected[i,4] = atten_correction_J_2002(Table_atten_corrected[i,0],Table_atten_corrected[i,1]*1e5,
    #                                                         Table_atten_corrected[i,3],Table_atten_corrected[i,4],oscillation,grain_size)
    for i in range(len(Table_atten_corrected)):
        Table_atten_corrected[i,3],Table_atten_corrected[i,4] = atten_correction(Table_atten_corrected[i,0],Table_atten_corrected[i,1]*1e5,
                                                             Table_atten_corrected[i,3],Table_atten_corrected[i,4],oscillation,grain_size)

    # correction for melts
    # These are relations from lab experiments. More details in Afonso et al., 2016 III
    # Function: lib.velocity_melt_correction_mantle (T (oC),P (GPa),VP (km/s),Vs (km/s)
    Table_atten_melt_corrected = np.copy(Table_atten_corrected)
    melt = np.zeros_like(Table_atten_melt_corrected[:,0])
    for i in range(len(Table_atten_melt_corrected)):
        Table_atten_melt_corrected[i,3],Table_atten_melt_corrected[i,4],melt[i] = velocity_melt_correction_mantle_Hammond_Humphreys(Table_atten_melt_corrected[i,0]-273.15,
                                                                                                 Table_atten_melt_corrected[i,1]/1e4,
                                                                 Table_atten_melt_corrected[i,3],Table_atten_melt_corrected[i,4])
    # append melt to the table
    #Table_atten_melt_corrected[:,5]=0.0
    Table_atten_melt_corrected[:,5]=melt[:]
    return Table_atten_melt_corrected


def mantle_melt_atten_correction_Behn2009(Table,grain_size,oscillation,COH):
    """
    Table : perplex table
    grain_size : grain size in mm.
    oscilation: oscillation period in seconds.
    """
    # correction using grain size = 10 mm and oscillatio period of 75 seconds.
    # Attenuation model of Jackson and Faul 2010
    # Function: lib.atten_correction (T (oC),P (Pascal),Vp (km/s),Vs (km/s),oscilation period (s), grain size (mm))
    Table_atten_corrected = np.copy(Table)
    """
    atten_correction_Behn2009(T,P,Vp,Vs,oscill,grain_size,COH):
    Behn et al., 2009 https://doi.org/10.1016/j.epsl.2009.03.014
    Input:
    T - Kelvin
    P - Pascal
    Vp,Vs - km/s 
    oscill - time period (seconds) 
    d - meter
    Output:
   """ 
    for i in range(len(Table_atten_corrected)):
        Table_atten_corrected[i,3],Table_atten_corrected[i,4] = atten_correction_Behn2009(Table_atten_corrected[i,0],Table_atten_corrected[i,1]*1e5,
                                                             Table_atten_corrected[i,3],Table_atten_corrected[i,4],oscillation,grain_size/1e3,COH)

    # correction for melts
    # These are relations from lab experiments. More details in Afonso et al., 2016 III
    # Function: lib.velocity_melt_correction_mantle (T (oC),P (GPa),VP (km/s),Vs (km/s)
    Table_atten_melt_corrected = np.copy(Table_atten_corrected)
    melt = np.zeros_like(Table_atten_melt_corrected[:,0])
    for i in range(len(Table_atten_melt_corrected)):
        Table_atten_melt_corrected[i,3],Table_atten_melt_corrected[i,4],melt[i] = velocity_melt_correction_mantle_Hammond_Humphreys(Table_atten_melt_corrected[i,0]-273.15,
                                                                                                 Table_atten_melt_corrected[i,1]/1e4,
                                                                 Table_atten_melt_corrected[i,3],Table_atten_melt_corrected[i,4])
    # append melt to the table
    #Table_atten_melt_corrected[:,5]=0.0
    Table_atten_melt_corrected[:,5]=melt[:]
    return Table_atten_melt_corrected

def mantle_melt_atten_correction_Behn2009_crust(Table,grain_size,oscillation,COH):
    """
    Table : perplex table
    grain_size : grain size in mm.
    oscilation: oscillation period in seconds.
    """
    # correction using grain size = 10 mm and oscillatio period of 75 seconds.
    # Attenuation model of Jackson and Faul 2010
    # Function: lib.atten_correction (T (oC),P (Pascal),Vp (km/s),Vs (km/s),oscilation period (s), grain size (mm))
    Table_atten_corrected = np.copy(Table)
    """
    atten_correction_Behn2009(T,P,Vp,Vs,oscill,grain_size,COH):
    Behn et al., 2009 https://doi.org/10.1016/j.epsl.2009.03.014
    Input:
    T - Kelvin
    P - Pascal
    Vp,Vs - km/s 
    oscill - time period (seconds) 
    d - meter
    Output:
   """ 
    for i in range(len(Table_atten_corrected)):
        Table_atten_corrected[i,3],Table_atten_corrected[i,4] = atten_correction_Behn2009_crust(Table_atten_corrected[i,0],Table_atten_corrected[i,1]*1e5,
                                                             Table_atten_corrected[i,3],Table_atten_corrected[i,4],oscillation,grain_size/1e3,COH)

    # correction for melts
    # These are relations from lab experiments. More details in Afonso et al., 2016 III
    # Function: lib.velocity_melt_correction_mantle (T (oC),P (GPa),VP (km/s),Vs (km/s)
    Table_atten_melt_corrected = np.copy(Table_atten_corrected)
    melt = np.zeros_like(Table_atten_melt_corrected[:,0])
    for i in range(len(Table_atten_melt_corrected)):
        Table_atten_melt_corrected[i,3],Table_atten_melt_corrected[i,4],melt[i] = velocity_melt_correction_mantle_Hammond_Humphreys(Table_atten_melt_corrected[i,0]-273.15,
                                                                                                 Table_atten_melt_corrected[i,1]/1e4,
                                                                 Table_atten_melt_corrected[i,3],Table_atten_melt_corrected[i,4])
    # append melt to the table
    #Table_atten_melt_corrected[:,5]=0.0
    Table_atten_melt_corrected[:,5]=melt[:]
    return Table_atten_melt_corrected
def mantle_melt_atten_correction_J_2002(Table,grain_size,oscillation):
    """
    Table : perplex table
    grain_size : grain size in mm.
    oscilation: oscillation period in seconds.
    """
    # correction using grain size = 10 mm and oscillatio period of 75 seconds.
    # Attenuation model of Jackson and Faul 2010
    # Function: lib.atten_correction (T (oC),P (Pascal),Vp (km/s),Vs (km/s),oscilation period (s), grain size (mm))
    Table_atten_corrected = np.copy(Table)
    #for i in range(len(Table_atten_corrected)):
    #    Table_atten_corrected[i,3],Table_atten_corrected[i,4] = atten_correction_J_2002(Table_atten_corrected[i,0],Table_atten_corrected[i,1]*1e5,
    #                                                         Table_atten_corrected[i,3],Table_atten_corrected[i,4],oscillation,grain_size)
    for i in range(len(Table_atten_corrected)):
        Table_atten_corrected[i,3],Table_atten_corrected[i,4] = atten_correction_J_2002(Table_atten_corrected[i,0],Table_atten_corrected[i,1]*1e5,
                                                             Table_atten_corrected[i,3],Table_atten_corrected[i,4],oscillation,grain_size)

    # correction for melts
    # These are relations from lab experiments. More details in Afonso et al., 2016 III
    # Function: lib.velocity_melt_correction_mantle (T (oC),P (GPa),VP (km/s),Vs (km/s)
    Table_atten_melt_corrected = np.copy(Table_atten_corrected)
    melt = np.zeros_like(Table_atten_melt_corrected[:,0])
    for i in range(len(Table_atten_melt_corrected)):
        Table_atten_melt_corrected[i,3],Table_atten_melt_corrected[i,4],melt[i] = velocity_melt_correction_mantle_Hammond_Humphreys(Table_atten_melt_corrected[i,0]-273.15,
                                                                                                 Table_atten_melt_corrected[i,1]/1e4,
                                                                 Table_atten_melt_corrected[i,3],Table_atten_melt_corrected[i,4])
    # append melt to the table
    #Table_atten_melt_corrected[:,5]=0.0
    Table_atten_melt_corrected[:,5]=melt[:]
    return Table_atten_melt_corrected
def crust_melt_atten_correction(Table,grain_size,oscillation):
    """
    Table : perplex table
    grain_size : grain size in mm.
    oscilation: oscillation period in seconds.
    """
    # correction using grain size = 10 mm and oscillatio period of 75 seconds.
    # Attenuation model of Jackson and Faul 2010
    # Function: lib.atten_correction (T (oC),P (Pascal),Vp (km/s),Vs (km/s),oscilation period (s), grain size (mm))
    Table_atten_corrected = np.copy(Table)
    for i in range(len(Table_atten_corrected)):
        Table_atten_corrected[i,3],Table_atten_corrected[i,4] = atten_correction(Table_atten_corrected[i,0],Table_atten_corrected[i,1]*1e5,
                                                             Table_atten_corrected[i,3],Table_atten_corrected[i,4],oscillation,grain_size)

    # correction for melts
    # These are relations from lab experiments. More details in Afonso et al., 2016 III
    # Function: lib.velocity_melt_correction_mantle (T (oC),P (GPa),VP (km/s),Vs (km/s)
    Table_atten_melt_corrected = np.copy(Table_atten_corrected)
    melt = np.zeros_like(Table_atten_melt_corrected[:,0])
    for i in range(len(Table_atten_melt_corrected)):
        Table_atten_melt_corrected[i,3],Table_atten_melt_corrected[i,4],melt[i] = velocity_melt_correction_crust(Table_atten_melt_corrected[i,0]-273.15,
                                                                                                 Table_atten_melt_corrected[i,1]/1e4,
                                                                 Table_atten_melt_corrected[i,3],Table_atten_melt_corrected[i,4])
    # append melt to the table
    Table_atten_melt_corrected[:,5]=melt[:]
    return Table_atten_melt_corrected
'''
def vel_to_temp(depth,Vs,Table):
    Temperature_out = []#np.zeros_like(tomo[:,1])
    Density_out     = []#np.zeros_like(tomo[:,1])
    diff_Vs         = []
    #Vp_out          = []#np.zeros_like(tomo[:,1])
    #Vs_out          = []#np.zeros_like(tomo[:,1])
    for i in range(len(depth)):
        P  = pressure_inter(depth[i])
        Vs_in = Vs[i]
        temp,dens,vp,vs,ind=lookup_vs_P_accurate(Vs_in,P.tolist(),Table)
        #Vp_out.append(vp)
        #Vs_out.append(vs)
        Temperature_out.append(temp)
        Density_out.append(dens)
        diff_Vs.append(Vs_in-vs)
    ### pasting the outputs to the input tomo table
    out=depth;
    out=np.column_stack((out,Temperature_out))
    out=np.column_stack((out,Density_out))
    out=np.column_stack((out,diff_Vs))

    return out
'''
def vel_to_temp(depth,Vs,Table):
    Temperature_out = []#np.zeros_like(tomo[:,1])
    Density_out     = []#np.zeros_like(tomo[:,1])
    diff_Vs         = []
    P_out           = []
    #Vp_out          = []#np.zeros_like(tomo[:,1])
    #Vs_out          = []#np.zeros_like(tomo[:,1])
    for i in range(len(depth)):
        P  = pressure_inter(depth[i])
        Vs_in = Vs[i]
        P_table,temp,dens,vp,vs=lookup_vs_P_accurate(Vs_in,P.tolist(),Table)
        #Vp_out.append(vp)
        #Vs_out.append(vs)
        P_out.append(P_table)
        Temperature_out.append(temp)
        Density_out.append(dens)
        diff_Vs.append(((Vs_in-vs)/Vs_in)*100)
    ### pasting the outputs to the input tomo table
    out=depth;
    out=np.column_stack((out,P_out))
    out=np.column_stack((out,Temperature_out))
    out=np.column_stack((out,Density_out))
    out=np.column_stack((out,diff_Vs))

    return out

def vel_to_temp_prop_out(depth,Vs,Table):
    """
    Input:
    depth : depth column in km.
    Vs    : tomography Vs velocity in km/s.
    Table : Perplex lookup table corrected form anelasticity and melt effects.

    Output:

    """
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
        P_table,temp,dens,vp,vs,m=lookup_vs_P_accurate_prop(Vs_in,P.tolist(),Table)
        #Vp_out.append(vp)
        #Vs_out.append(vs)
        P_out.append(P_table)
        Temperature_out.append(temp)
        Density_out.append(dens)
        Vs_out.append(vs)
        Vp_out.append(vp)
        diff_Vs.append(((Vs_in-vs)/Vs_in)*100)
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


def vel_to_temp_P_in(depth,Vs,Table,P_func):
    Temperature_out = []#np.zeros_like(tomo[:,1])
    Density_out     = []#np.zeros_like(tomo[:,1])
    diff_Vs         = []
    P_out           = []
    #Vp_out          = []#np.zeros_like(tomo[:,1])
    #Vs_out          = []#np.zeros_like(tomo[:,1])
    for i in range(len(depth)):
        P  = P_func(depth[i])
        Vs_in = Vs[i]
        P_table,temp,dens,vp,vs=lookup_vs_P_accurate(Vs_in,P.tolist(),Table)
        #Vp_out.append(vp)
        #Vs_out.append(vs)
        P_out.append(P_table)
        Temperature_out.append(temp)
        Density_out.append(dens)
        diff_Vs.append(Vs_in-vs)
    ### pasting the outputs to the input tomo table
    out=depth;
    out=np.column_stack((out,P_out))
    out=np.column_stack((out,Temperature_out))
    out=np.column_stack((out,Density_out))
    out=np.column_stack((out,diff_Vs))

    return out
def vel_to_temp_AHCZ(tomo,geo_func,moho_func,mc_func,uc_func,sedi_func,Crust_LC_use,Crust_MC_use,Crust_UC_use,Sediments_use,DMM_use,Archon_use,Proton_use,Tecton_use,Ocean_use):
    Temperature_out = []#np.zeros_like(tomo[:,1])
    Density_out     = []#np.zeros_like(tomo[:,1])
    melt_out     = [] #np.zeros_like(tomo[:,1])
    Vp_out     = [] #np.zeros_like(tomo[:,1])
    Vs_out     = [] #np.zeros_like(tomo[:,1])
    diff_Vs         = []
    P_out           = []
    #Vp_out          = []#np.zeros_like(tomo[:,1])
    #Vs_out          = []#np.zeros_like(tomo[:,1])
    for i in range(len(tomo)):
        moho = moho_func(tomo[i,0],tomo[i,1])
        mc = mc_func(tomo[i,0],tomo[i,1])
        uc = uc_func(tomo[i,0],tomo[i,1])
        sedi = sedi_func(tomo[i,0],tomo[i,1])
        P  = pressure_inter(tomo[i,2])
        Vs_in = tomo[i,3]
        Composition = geo_func(tomo[i,0],tomo[i,1])
        ## asthenosphere
        if tomo[i,2] > moho:
            # First try to DMM composition to find out LAB i.e. 1300 oC
            P_table,temp,dens,vp,vs,m=lookup_vs_P_accurate_prop(Vs_in,P.tolist(),DMM_use)
            # Now check if the converted temperature if above 1300 oC if yes then pass
            # if not then it is lithospheric mantle and used lithospheric mantle composition
            if temp < 1300.0:
                if Composition == 6:

                    P_table,temp,dens,vp,vs,m = lookup_vs_P_accurate_prop(Vs_in,P.tolist(),Archon_use)

                elif Composition == 5 or Composition == 4 or Composition == 3:

                    P_table,temp,dens,vp,vs,m = lookup_vs_P_accurate_prop(Vs_in,P.tolist(),Proton_use)

                elif Composition == 2 or Composition == 1:

                    P_table,temp,dens,vp,vs,m = lookup_vs_P_accurate_prop(Vs_in,P.tolist(),Tecton_use)
                else:
                    P_table,temp,dens,vp,vs,m = lookup_vs_P_accurate_prop(Vs_in,P.tolist(),Ocean_use)
            else:
                pass
        elif tomo[i,2] <= moho and tomo[i,2] > mc:
            P_table,temp,dens,vp,vs,m = lookup_vs_P_accurate_prop(Vs_in,P.tolist(),Crust_LC_use)
        elif tomo[i,2] <= mc and tomo[i,2] > uc:
            P_table,temp,dens,vp,vs,m = lookup_vs_P_accurate_prop(Vs_in,P.tolist(),Crust_MC_use)
        elif tomo[i,2] <= uc and tomo[i,2] > sedi:
            P_table,temp,dens,vp,vs,m = lookup_vs_P_accurate_prop(Vs_in,P.tolist(),Crust_UC_use)
        else:
            P_table,temp,dens,vp,vs,m = lookup_vs_P_accurate_prop(Vs_in,P.tolist(),Sediments_use)
        # gather conversions
        P_out.append(P_table)
        Temperature_out.append(temp)
        Density_out.append(dens)
        Vs_out.append(vs)
        Vp_out.append(vp)
        diff_Vs.append(((Vs_in-vs)/Vs_in)*100)
        melt_out.append(m)
    ### pasting the outputs to the input tomo table
    out=tomo[:,0];
    out=np.column_stack((out,tomo[:,1]))
    out=np.column_stack((out,tomo[:,2]))
    out=np.column_stack((out,P_out))
    out=np.column_stack((out,Temperature_out))
    out=np.column_stack((out,Density_out))
    out=np.column_stack((out,Vp_out))
    out=np.column_stack((out,tomo[:,3]))
    out=np.column_stack((out,Vs_out))
    out=np.column_stack((out,diff_Vs))
    out=np.column_stack((out,melt_out))
    return out
def half_space_ocean(depth,Age,T_mantle,T_surface):
    kappa = 1e-6;
    SecYear = 3600*24*365.25
    T = T_surface + (T_mantle-T_surface)*(math.erf(depth/(2*math.sqrt(kappa*Age*SecYear*1.0e6))))
    return T
