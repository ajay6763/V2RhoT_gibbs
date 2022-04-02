import multiprocessing
import random
import time
import numpy as np
import V2RhoT_gibbs_lib as lib




def worker(name, data, table: str) -> None:
    print(f'Started worker {name}')
    ## run conversion here
    systime = time.time()
    out_gibbs=lib.vel_to_temp(data[:,2],data[:,3],table)
    out_save=np.zeros_like(data[:,1])
    out_save=data[:,0]
    out_save=np.column_stack((out_save,data[:,1]))
    out_save=np.column_stack((out_save,out_gibbs[:,0]))
    out_save=np.column_stack((out_save,out_gibbs[:,1]))
    out_save=np.column_stack((out_save,out_gibbs[:,2]))
    out_save=np.column_stack((out_save,out_gibbs[:,3]))
    out_save=np.column_stack((out_save,out_gibbs[:,4]))
    np.savetxt(str(name)+'_vel_converted.txt',out_save,header="#x(km) y(km) depth(km) Pressure(bar) Temperature(oC) Density(kg/m3) Vs_diff(km/s)",comments='',fmt='%10.3f')
    worker_time=(time.time()-systime)/60.0
    print(f'{name} worker finished in {worker_time} min.')
    

if __name__ == '__main__':
    ####################################################
    systime = time.time()
    # Load table
    DMM_no_atten = np.loadtxt('./databases/DMM_HP',comments='#')
    ########################################
    #Correction for melts and anleasticity
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
    table=DMM_atten_melt_corrected
    ###############
    # data format
    # age(Ma) depth(km) Vs(km/s)
    tomo_in = np.loadtxt('./data_tomo/Test_parallel.dat',comments='#')
    #tomo_in = np.loadtxt('./data_tomo/V_mean.txt',comments='#')

    #print("Time spent on loading data: {:.2f} sec".format(time.time()-tcalc))
    print("Time spent on loading data: {:.2f} sec".format((time.time()-systime)))

    #####################
    # find number of processors
    no_processors = int(multiprocessing.cpu_count())
    print("Total number of available processors: "+str(no_processors))
    #no_processors = 4
    
    no_of_parts = int(len(tomo_in)/no_processors)
    #no_of_parts = 1000
    print("Length of the input model: "+str(len(tomo_in)))
    print("Length of chunk to run on each processor: "+str(no_of_parts))
    systime = time.time()
    processes = []
    for i in range(no_processors):
        start_index = no_of_parts*i
        # If the input file factorization is not whole numnber
        # then what i am doing here is that i will put the last chunk start
        # and untill the length of the input on the last core
        if i == no_processors-1:
            end_index   = len(tomo_in)
        else:
            end_index   = start_index + no_of_parts  
        process = multiprocessing.Process(target=worker, 
                                          args=(f'Core_{i+1}',tomo_in[start_index:end_index,:],DMM_no_atten))
        processes.append(process)
        process.start()
    for proc in processes:
        proc.join()
    print("Total execution time: {:.1f} min".format((time.time()-systime)/60.0))
    ###############
    # gather all output
    out_save = np.loadtxt('Core_1_vel_converted.txt',comments='#')
    for i in range(1,no_processors):
        out = np.loadtxt('Core_'+str(i+1)+'_vel_converted.txt',comments='#') 
        out_save=np.append(out_save,out,axis=0)
    np.savetxt('Final_vel_converted.txt',out_save,header="#x(km) y(km) depth(km) Pressure(bar) Temperature(oC) Density(kg/m3) Vs_diff(km/s)",comments='',fmt='%10.3f')
    
