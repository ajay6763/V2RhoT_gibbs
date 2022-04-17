#LITMOD 4.0
#Script for parallel computation of the gravity field
#Eldar Baykiev, 2018-2020

import numpy as np
import time
import multiprocessing
import subprocess
import os

import sys
import glob
import V2RhoT_gibbs_lib as lib

systime = time.time()

#print "**************************************************"
#print "Calculation of Gravity Potentials"
#print ''

tsplit = time.time()

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


###############
# data format
# age(Ma) depth(km) Vs(km/s)
# 
tomo_NA_stack = np.loadtxt('./data_tomo/NA_age_vel_stack.dat',comments='#')

#####################
# find number of processors
no_processors = int(multiprocessing.cpu_count())
print("Number of processors: "+str(no_processors))

####################
# Now what I will do it split the 
# input file into parts which is a multiple 
# of the total no of available processors
#
no_of_parts = int(len(tomo_NA_stack)/no_processors)

#print(filenames)
tsplit = time.time()-tsplit
###########################
## run conversion here
out_gibbs=lib.vel_to_temp(tomo_NA_stack[:,1],tomo_NA_stack[:,2],DMM_atten_melt_corrected)


print('Parallel computation... ', end='')
tcalc = time.time()
print('done')
print("Time spent on computations: {:.2f} sec".format(time.time()-tcalc))

print('Collect output...', end='')
tsum = time.time()
print('done')
print("Time spent on operations with files: {:.2f} sec".format(time.time()-tsum))

print("Total execution time: {:.1f} min".format((time.time()-systime)/60.0))

print('\n')