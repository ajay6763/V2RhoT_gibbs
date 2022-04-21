from datetime import date
import multiprocessing, os, sys, getopt, time
import numpy as np
import V2RhoT_gibbs_lib as lib


def worker(name, data, table,outdir: str) -> None:
    print(f'Started worker {name}')
    ## run conversion here
    systime = time.time()
    #out_gibbs=lib.vel_to_temp(data[:,2],data[:,3],table)
    out_gibbs=lib.vel_to_temp_prop_out(data[:,2],data[:,3],table)
    out_save=np.zeros_like(data[:,1])
    out_save=data[:,0]
    out_save=np.column_stack((out_save,data[:,1]))
    out_save=np.column_stack((out_save,out_gibbs[:,0]))
    out_save=np.column_stack((out_save,out_gibbs[:,1]))
    out_save=np.column_stack((out_save,out_gibbs[:,2]))
    out_save=np.column_stack((out_save,out_gibbs[:,3]))
    out_save=np.column_stack((out_save,out_gibbs[:,4]))
    out_save=np.column_stack((out_save,out_gibbs[:,5]))
    out_save=np.column_stack((out_save,out_gibbs[:,6]))
    out_save=np.column_stack((out_save,out_gibbs[:,7]))
    np.savetxt(str(outdir)+"/"+str(name)+'_vel_converted.txt',out_save,header="#x(km) y(km) depth(km) Pressure(bar) Temperature(oC) Density(kg/m3) Vp(km/s) Vs(km/s) Vs_diff(%) Pseudo-metls(%)",comments='',fmt='%10.3f')
    #np.savetxt(str(outdir)+"/"+str(name)+'_vel_converted.txt',out_save,header="#x(km) y(km) depth(km) Pressure(bar) Temperature(oC) Density(kg/m3) Vp (km/s) Vs (km/s) Vs_diff(km/s)",comments='',fmt='%10.3f')
    worker_time=(time.time()-systime)/60.0
    print(f'{name} worker finished in {worker_time} min.')
    
def main(argv):
    no_of_processes     = int(1)
    inputfile           = str('')
    outputfile          = str('conversion_out.txt')
    materialfile        = str('DMM_HP')
    grain_size          = float(10)
    oscillation_period  = float(75)
    try :
        opts, args = getopt.getopt(argv,"h:p:i:o:m:g:s",["ifile","ofile","mfile","gsize","operiod"])
    except getopt.GetoptError:
        print('\n###########################################################################################')
        print(' Simple run example: parallel_conversion.py -i <inputfile> -o <outputfile>\n\
             Below are the available options which you can pass as a command line argument:\n\
            -h : help\n\
            -p : no of parts to run in parallel (e.g., no of available cores). Default is 1\n\
            -i : input file name in the data_tomo folder. Format x(*) y(*) depth(km) Vs(km/s)\n\
            -o : output file name which will be saved in output folder\n\
            -m : name of the material file in databases folder. Default is DMM_HP\n\
            -g : grain size in mm. Default is 10mm\n\
            -s : oscillation period in seconds. Default is 75 seconds\nAll of the input options will be written in the output file as comments.')
        print('###########################################################################################\n')

        sys.exit(2)
    if (len(opts)!=0):
        for opt, arg in opts:
            if opt == '-h':
                print('parallel_conversion.py -i <inputfile> -o <outputfile>\n -h : help\n -p : no of parts to run in parallel (e.g., no of available cores)\n\
                 -i : input file name in the data_tomo folder\n -o : output file name which will be saved in output folder. Default output file name is conversion_out.txt')
                sys.exit()
            elif opt in ("-p", "--processes"):
                no_of_processes = int(arg)
            elif opt in ("-i", "--ifile"):
                inputfile = arg
            elif opt in ("-o", "--ofile"):
                outputfile = arg
            elif opt in ("-m", "--mfile"):
                materialfile = arg
            elif opt in ("-g", "--gsize"):
                grain_size = float(arg)
            elif opt in ("-s", "--operiod"):
                oscillation_period = float(arg)
        print('\n###########################################')
        print('Ýour options are:')
        print('Total no of processes',no_of_processes)    
        print('Input file is', inputfile)
        print('Output file is', outputfile)
        print('Material file is', materialfile)
        print('Grain size is', grain_size,'mm')
        print('Oscillation period is', oscillation_period,'seconds')
        print('\n###########################################')
        
        ### get current directory
        path = os.getcwd()
        ### make output directory
        outdir = str(path)+"/output"
        # check if it already exists
        isExist = os.path.exists(outdir)
        if not isExist:
            print('\n###########################################')
            print('Output directory does not exist. Making one for you.')
            print('\n###########################################')
            os.makedirs(outdir)
        else:
            print('\n###########################################') 
            print('Output directory exists. It will be overwritted.')
            '''
            print('Output directory exists. To save it will be renamed with _bu appended to it.\n\
            New output directory will be made for you new results. But re')
            os.mkdirs(outdir,str(path)+"/back_ups")
            os.move(outdir,str(path)+"/back_ups")
            os.makedirs(outdir)
            '''
            print('\n###########################################')
        ####################################################
        systime = time.time()
        ###############
        # data format
        # age(Ma) depth(km) Vs(km/s)
        try:
            tomo_in = np.loadtxt(str(path)+"/data_tomo/"+str(inputfile),comments='#')
        except:
            print('\n###########################################')
            print('Could not find input tomography file',str(path)+"/data_tomo/"+str(inputfile),'\nMake sure you have this file.')
            print('\n###########################################')
            sys.exit()

        #tomo_in = np.loadtxt('./data_tomo/V_mean.txt',comments='#')

        #print("Time spent on loading data: {:.2f} sec".format(time.time()-tcalc))

        # Load table
        try:
            DMM_no_atten = np.loadtxt(str(path)+'./databases/'+str(materialfile),comments='#')
        except:
            print('\n###########################################')
            print('Could not find the material file',str(path)+'./databases/'+str(materialfile),'\nMake sure you have this file.')
            sys.exit(0)
            print('\n###########################################')
        ########################################
        #Correction for melts and anleasticity
        ########################################
        # correction using grain size = 10 mm and oscillatio period of 75 seconds.
        # Attenuation model of Jackson and Faul 2010
        # Function: lib.atten_correction (T (oC),P (Pascal),VP (km/s),Vs (km/s),oscilation period (s), grain size (mm))
        DMM_atten_corrected = np.copy(DMM_no_atten)
        for i in range(len(DMM_atten_corrected)):
            DMM_atten_corrected[i,3],DMM_atten_corrected[i,4] = lib.atten_correction(DMM_atten_corrected[i,0],DMM_atten_corrected[i,1]*100000,
                                                                DMM_atten_corrected[i,3],DMM_atten_corrected[i,4],grain_size,oscillation_period)
        # correction for melts
        # These are relations from lab experiments. More details in Afonso et al., 2016 III
        # Function: lib.velocity_melt_correction_mantle (T (oC),P (GPa),VP (km/s),Vs (km/s),oscilation period (s), grain size (mm))
        DMM_atten_melt_corrected = np.copy(DMM_atten_corrected)
        melt = np.zeros_like(DMM_atten_melt_corrected[:,0])
        for i in range(len(DMM_atten_melt_corrected)):
            DMM_atten_melt_corrected[i,3],DMM_atten_melt_corrected[i,4],melt[i] = lib.velocity_melt_correction_mantle(DMM_atten_melt_corrected[i,0]-273.15,
                                                                                                    DMM_atten_melt_corrected[i,1]/10000,
                                                                    DMM_atten_melt_corrected[i,3],DMM_atten_melt_corrected[i,4])
        #DMM_atten_melt_corrected[:,5]=0
        DMM_atten_melt_corrected[:,5]=melt[:]
        table=DMM_atten_melt_corrected
        print('\n###########################################')
        print("Time spent on loading data: {:.2f} sec".format((time.time()-systime)))
        print('\n###########################################')

        #####################
        # find number of processors
        #no_processors = int(multiprocessing.cpu_count())
        #print("Total number of available processors: "+str(no_processors))
        #no_processors = 4
        
        no_of_parts = int(len(tomo_in)/no_of_processes)
        #no_of_parts = 1000
        print('\n###########################################')
        print("Length of the input model: "+str(len(tomo_in)))
        print("Length of chunk to run on each processor: "+str(no_of_parts))
        print('\n###########################################')
        systime = time.time()
        processes = []
        for i in range(no_of_processes):
            start_index = no_of_parts*i
            # If the input file factorization is not whole numnber
            # then what i am doing here is that i will put the last chunk start
            # and untill the length of the input on the last core
            if i == no_of_processes-1:
                end_index   = len(tomo_in)
            else:
                end_index   = start_index + no_of_parts  
            process = multiprocessing.Process(target=worker, 
                                            args=(f'Process_{i+1}',tomo_in[start_index:end_index,:],table,outdir))
            processes.append(process)
            process.start()
        for proc in processes:
            proc.join()
        print('\n###########################################')
        print("Total execution time: {:.1f} min".format((time.time()-systime)/60.0))
        print('Hope that was fast enough for you. If not then may be try incrasing no of processes with -p option')
        #print('There are lot of for loops in the code and could be improved.\n')
        print('Enjoy your conversion results in the output folder :)')
        print('\n###########################################')
        ###############
        # gather all output
        # make comments that includes all the parameters used for conversion and put in the final converted file
        meta_data = "#Created on: " + str(date.today()) +"\n#Input file is: " + str(inputfile) +"\n#Output file is: " + str(outputfile) +"\n#Material file is: " + str(materialfile) +"\n#Grain size is: " + str(grain_size) +" mm"+"\n#Oscillation period is: " + str(oscillation_period) +" seconds\n"

        out_save = np.loadtxt(str(outdir)+'/Process_1_vel_converted.txt',comments='#')
        for i in range(1,no_of_processes):
            out = np.loadtxt(str(outdir)+'/Process_'+str(i+1)+'_vel_converted.txt',comments='#') 
            out_save=np.append(out_save,out,axis=0)
        np.savetxt(str(outdir)+'/'+str(outputfile),out_save,header="#x(km) y(km) depth(km) Pressure(bar) Temperature(oC) Density(kg/m3) Vp(km/s) Vs(km/s) Vs_diff(%) Pseudo-melts(%)",comments=meta_data,fmt='%10.3f')
    else:
        print('\n###########################################')
        print('You did not provide required input.')
        print('Run the code with -h option for help.')
        print('###########################################\n')

if __name__ == '__main__':
    main(sys.argv[1:])