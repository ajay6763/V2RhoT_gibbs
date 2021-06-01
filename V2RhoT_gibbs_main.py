import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from scipy import interpolate
import math
import V2RhoT_gibbs_lib as lib

########################################
# Loading perplex tables
########################################
#DMM_lith = np.loadtxt('./../databases/DMM_lith',skiprows=12)
DMM_asth_atten_corrected = np.loadtxt('./../databases/DMM_asth',skiprows=12)
#Tc_1_no_atten = np.loadtxt('./../databases/Avg_Gnt_Tecton_no_atten')
#DMM_no_atten = np.loadtxt('./../databases/DMM_lith_no_atten')
#Tc_1_atten_perplex = np.loadtxt('./../databases/Avg_Gnt_Tecton',skiprows=12)


########################################
# Correction for melts and anleasticity
########################################
''''
DMM_atten_corrected = np.copy(DMM_asth)
for i in range(len(DMM_atten_corrected)):
    ## for Vp
    DMM_atten_corrected[i,3],DMM_atten_corrected[i,4] = atten_correction(DMM_atten_corrected[i,0],DMM_atten_corrected[i,1]*100000,
''''                                                             DMM_atten_corrected[i,3],DMM_atten_corrected[i,4],75,10)
DMM_atten_melt_corrected = np.copy(DMM_atten_corrected)
melt = np.zeros_like(DMM_atten_melt_corrected[:,0])
for i in range(len(DMM_atten_melt_corrected)):
    ## for Vs
    DMM_atten_melt_corrected[i,3],DMM_atten_melt_corrected[i,4],melt[i] = lib.velocity_melt_correction(DMM_atten_melt_corrected[i,0]-273.15,
                                                                                             DMM_atten_melt_corrected[i,1]/10000,
                                                             DMM_atten_melt_corrected[i,3],DMM_atten_melt_corrected[i,4])

########################################
# Defining pressure function based on Ak133
# I choose ak135 becaouse the pressures from the LitMod2D_2.0 ref model
# are practically similar to ak135
########################################
ak135 = np.loadtxt('./../ak135f.txt',skiprows=1+2)
ak135_P = 9.8*ak135[:,0]*1e3*ak135[:,1]*1e3*1e-5
pressure_inter = interpolate.interp1d(ak135[:,0],ak135_P)
# pressure in the LitMod ref is in pascal so converting it above to bar becaouse
# look-up table e.g. 99 have pressure in bar
#pressure_inter = interpolate.interp1d(LitMod_ref[:,0],LitMod_ref[:,5]/100000)
# pressure in the LitMod ref is in pascal so converting it above to bar becaouse
# look-up table e.g. 99 have pressure in bar


###########################################
# load CSEM tomography
###########################################
tomo_in = np.loadtxt('West_Med_CSEM_UTM.txt',skiprows=1)
Vs = np.zeros_like(tomo_in[:,3])
Vs=(2.0*tomo_in[:,3]+tomo_in[:,4])/3.0

###########################################
# load MeRe2020 tomography
###########################################
tomo_in = np.loadtxt('West_Med_MeRe2020_UTM.txt',skiprows=1)
Vs = np.zeros_like(tomo_in[:,3])
Vs = tomo_in[:,3]

#################################################
# find indices for which to convert
# filter based on depth
###################################################
ind = np.where( (tomo_in[:,2] <= 400.0) & (tomo_in[:,2] >=50) )
for item in ind:
    index=item
index=index.tolist()

################################################
# Doing the conversion
################################################
Temp,D,Vp,Vs=lib.vel_to_temp(tomo_in[index[:],2],Vs[index[:]],Tc_1_atten_melt_corrected,index)

################################################
# Saving
################################################
alps_Tc=np.zeros_like(Temp[:])
alps_Tc=tomo_in[index[:],0]
alps_Tc=np.column_stack((alps_Tc,tomo_in[index[:],1]))
alps_Tc=np.column_stack((alps_Tc,-tomo_in[index[:],2]))
alps_Tc=np.column_stack((alps_Tc,Temp))
alps_Tc=np.column_stack((alps_Tc,D))
alps_Tc=np.column_stack((alps_Tc,Vp))
alps_Tc=np.column_stack((alps_Tc,Vs))

np.savetxt('alps_Tc_1_Temp_MeRe2020.txt',alps_Tc,header="x y z temperature density Vp Vs",comments='')