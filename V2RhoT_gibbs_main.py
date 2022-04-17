import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from scipy import interpolate
import math
import time
import V2RhoT_gibbs_lib as lib
systime = time.time()
########################################
# Loading perplex tables
########################################
#DMM_lith = np.loadtxt('./../databases/DMM_lith',skiprows=12)
DMM_no_atten = np.loadtxt('./databases/DMM_HP',comments='#')
#Tc_1_no_atten = np.loadtxt('./../databases/Avg_Gnt_Tecton_no_atten')
#DMM_no_atten = np.loadtxt('./../databases/DMM_lith_no_atten')
#Tc_1_atten_perplex = np.loadtxt('./../databases/Avg_Gnt_Tecton',skiprows=12)


########################################
# Correction for melts and anleasticity
########################################
# correction using grain size = 10 mm and oscillatio period of 75 seconds.
# Attenuation model of Jackson and Faul 2010
# Function: lib.atten_correction (T (oC),P (Pascal),VP (km/s),Vs (km/s),oscilation period (s), grain size (mm))
DMM_atten_corrected = np.copy(DMM_no_atten)
for i in range(len(DMM_atten_corrected)):
    DMM_atten_corrected[i,3],DMM_atten_corrected[i,4] = lib.atten_correction(DMM_atten_corrected[i,0],DMM_atten_corrected[i,1]*100000,
                                                         DMM_atten_corrected[i,3],DMM_atten_corrected[i,4],75,10)

# correction for melts
# These are relations from lab experiments. More details in Afonso et al., 2016 III
# Function: lib.velocity_melt_correction_mantle (T (oC),P (GPa),VP (km/s),Vs (km/s),oscilation period (s), grain size (mm))
DMM_atten_melt_corrected = np.copy(DMM_atten_corrected)
melt = np.zeros_like(DMM_atten_melt_corrected[:,0])
for i in range(len(DMM_atten_melt_corrected)):
    DMM_atten_melt_corrected[i,3],DMM_atten_melt_corrected[i,4],melt[i] = lib.velocity_melt_correction_mantle(DMM_atten_melt_corrected[i,0]-273.15,
                                                                                             DMM_atten_melt_corrected[i,1]/10000,
                                                             DMM_atten_melt_corrected[i,3],DMM_atten_melt_corrected[i,4])

########################################
# Defining pressure function based on Ak133
# I choose ak135 becaouse the pressures from the LitMod2D_2.0 ref model
# are practically similar to ak135
########################################
ak135 = np.loadtxt('./databases/ak135f.txt',skiprows=1+2)
ak135_P = 9.8*ak135[:,0]*1e3*ak135[:,1]*1e3*1e-5
pressure_inter = interpolate.interp1d(ak135[:,0],ak135_P)
# pressure in the LitMod ref is in pascal so converting it above to bar becaouse
# look-up table e.g. 99 have pressure in bar
#pressure_inter = interpolate.interp1d(LitMod_ref[:,0],LitMod_ref[:,5]/100000)
# pressure in the LitMod ref is in pascal so converting it above to bar becaouse
# look-up table e.g. 99 have pressure in bar



###############
# data format
# age(Ma) depth(km) Vs(km/s)
# 
tomo_NA_stack = np.loadtxt('./data_tomo/V_mean.txt',comments='#')

################################################
# Doing the conversion
################################################
#  Output  = function       (depth(km) , Vs(km/s), Table) 
#   out    = lib.vel_to_temp(tomo_in[index[:],2],Vs[index[:]],Tc_1_atten_melt_corrected)
# The out in sequence is -depth(km), Pressure (bar), Temperature(oC), Density(kg/m3), V_obs-V_model
start_time = time.time()
print("\n Started at:",start_time)
out_gibbs=lib.vel_to_temp(tomo_NA_stack[:,1],tomo_NA_stack[:,2],DMM_atten_melt_corrected)
end_time = time.time()
print("\n Ended at:",end_time)
print("\n Total time taken to run: seconds",end_time-start_time)
################################################
# Saving
################################################
out_save=np.zeros_like(tomo_NA_stack[:,1])
out_save=tomo_NA_stack[:,0]
out_save=np.column_stack((out_save,out_gibbs[:,1]))
out_save=np.column_stack((out_save,out_gibbs[:,0]))
out_save=np.column_stack((out_save,out_gibbs[:,1]))
out_save=np.column_stack((out_save,out_gibbs[:,2]))
out_save=np.column_stack((out_save,out_gibbs[:,3]))
out_save=np.column_stack((out_save,out_gibbs[:,4]))
np.savetxt('NA_Age_vel_converted.txt',out_save,header="#Age(Myr) depth(km) Pressure(bar) Temperature(oC) Density(kg/m3) Vs_diff(km/s)",comments='',fmt='%10.3f')
print("Total execution time: {:.1f} min".format((time.time()-systime)/60.0))