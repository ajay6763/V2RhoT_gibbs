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
    inputdir            = str('./data_tomo')
    outputdir           = str('./output')
    atten_model         = int(1)
    grain_size          = float(10.0)
    oscillation_period  = float(75.0)
    COH_val 	        = float(50.0)
    #oscillation_period  = 75
    try :
        opts,args = getopt.getopt(argv,"h:p:I:i:O:o:m:A:g:s:COH",["idir","ifile","odir","ofile","mfile","Amodel","gsize","operiod","water"])
        #print(opts)
        #print(args)
    except getopt.GetoptError:
        print('\n###########################################################################################')
        print(' Simple run example: parallel_conversion.py -i <inputfile> -o <outputfile>\n\
             Below are the available options which you can pass as a command line argument:\n\
            -h : help\n\
            -p : no of parts to run in parallel (e.g., no of available cores). Default is 1\n\
            -I : input directory. Default is ./data_tomo\n\
            -i : input file name in the data_tomo folder. Format x(*) y(*) depth(km) Vs(km/s)\n\
            -O : output directory. Default is ./output\n\
            -o : output file name which will be saved in output folder\n\
            -m : name of the material file in databases folder. Default is DMM_HP\n\
            -A : attenuation model flag: 1 for Jackson and Faul 2010, 2 for Behn et al., 2009. If you choose 2 then you will have to supply COH. See COH option.\n\
            -g : grain size in mm. Default is 10mm\n\
            -s : oscillation period in seconds. Default is 75 seconds\n\
            -COH : Water contnent in  H/10**6Si. Defualt is 50 H/10**6Si which almost dry. Note this will be used if you choose attenuation model 2.\n\
            -geo : Geology file in the data_geo folder. Format Format x(*) y(*) geo(codes). Note: This option is not active at the momemnt.\
            \nAll of the input options will be written in the output file as comments.')
        print('###########################################################################################\n')

        sys.exit(2)
    if (len(opts)!=0):
        for opt, arg in opts:
            if opt == '-h':
                print('\n###########################################################################################')
                print(' Simple run example: parallel_conversion.py -i <inputfile> -o <outputfile>\n\
                Below are the available options which you can pass as a command line argument:\n\
                -h : help\n\
                -p : no of parts to run in parallel (e.g., no of available cores). Default is 1\n\
                -i : input file name in the data_tomo folder. Format x(*) y(*) depth(km) Vs(km/s)\n\
                -o : output file name which will be saved in output folder\n\
                -I : input directory. Default is ./data_tomo\n\
                -O : output directory. Default is ./output\n\
                -m : name of the material file in databases folder. Default is DMM_HP\n\
                -A : attenuation model flag: 1 for Jackson and Faul 2010 (default model), 2 for Behn et al., 2009. If you choose 2 then you will have to supply COH. See COH option.\n\
                -g : grain size in mm. Default is 10mm\n\
                -s : oscillation period in seconds. Default is 75 seconds\n\
                -COH : Water contnent in  H/10**6Si. Defualt is 50 H/10**6Si which almost dry. Note this will be used if you choose attenuation model 2.\n\
                -geo : Geology file in the data_geo folder. Format Format x(*) y(*) geo(codes). Note: This option is not active at the momemnt.\
                \nAll of the input options will be written in the output file as comments.')
                print('###########################################################################################\n')               
                sys.exit()
            elif opt in ("-p", "--processes"):
                no_of_processes = int(arg)
            elif opt in ("-i", "--ifile"):
                inputfile = arg
            elif opt in ("-o", "--ofile"):
                outputfile = arg
            elif opt in ("-I", "--idir"):
                inputdir = str(arg)
            elif opt in ("-O", "--odir"):
                outputdir =str(arg)
            elif opt in ("-m", "--mfile"):
                materialfile = arg
            elif opt in ("-A", "--Amodel"):
                atten_model = int(arg)
            elif opt in ("-g", "--gsize"):
                grain_size = float(arg)
            elif opt in ("-s", "--operiod"):
                oscillation_period = float(arg)
            elif opt in ("-COH", "--water"):
                COH_val = float(arg)
            else:
                pass
        print('\n###########################################')
        print('Your options are:')
        print('Total no of processes',no_of_processes)    
        print('Input file is', inputfile)
        print('Output file is', outputfile)
        print('Material file is', materialfile)
        print('Input folder is', inputdir)
        print('Output folder is', outputdir)
        print('Attenuation model is', atten_model)
        print('Grain size is', grain_size,'mm')
        print('Oscillation period is', oscillation_period,'seconds')
        print('COH is (used when Attenuation model is 1)', COH_val,'H/10**6Si')
        print('\n###########################################')
        
        ### get current directory
        path = os.getcwd()
        isExist = os.path.exists(outputdir)
        if not isExist:
            print('\n###########################################')
            print('Output directory does not exist. Making one for you.')
            print('\n###########################################')
            os.makedirs(outputdir)
        else:
            print('\n###########################################')
            print('Output directory exists. It will be overwritten.')
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
            tomo_in = np.loadtxt(inputdir+'/'+inputfile,comments='#')
        except:
            print('\n###########################################')
            print('Could not find input tomography file',inputdir+'/'+inputfile,'\nMake sure you have this file.')
            print('\n###########################################')
            sys.exit()

        #tomo_in = np.loadtxt('./data_tomo/V_mean.txt',comments='#')

        #print("Time spent on loading data: {:.2f} sec".format(time.time()-tcalc))

        # Load table
        try:
            DMM_no_atten = np.loadtxt(str(path)+'/databases/'+str(materialfile),comments='#')
            #print('dwerwerwerwer')
            #print(atten_model)
            #table        = lib.mantle_melt_atten_correction(DMM_no_atten,grain_size,oscillation_period)
            if atten_model == 2:
                print('Chose 2 Atten model.')
                table        = lib.mantle_melt_atten_correction_Behn2009(DMM_no_atten,grain_size,oscillation_period,COH_val)
            elif atten_model == 1:
                print('Chose 1 Atten model.')
                table        =  lib.mantle_melt_atten_correction(DMM_no_atten,grain_size,oscillation_period)
            else:
                print('You choose wronge attenuation model.Available options are 1 and 2 . For help run with -h option')
                sys.exit()
        except:
            print('\n###########################################')
            print('Could not find the material file',str(path)+'/databases/'+str(materialfile),'\nMake sure you have this file.')
            sys.exit(0)
            print('\n###########################################')
        
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
                                            args=(f'Process_{i+1}',tomo_in[start_index:end_index,:],table,outputdir))
            processes.append(process)
            process.start()
        for proc in processes:
            proc.join()
        print('\n###########################################')
        print("Total execution time: {:.1f} min".format((time.time()-systime)/60.0))
        print('Hope that was fast enough for you. If not then maybe try increasing no of processes with -p option')
        #print('There are lot of for loops in the code and could be improved.\n')
        print('Enjoy your conversion results in the output folder :)')
        print('\n###########################################')
        ###############
        # gather all output
        # make comments that includes all the parameters used for conversion and put in the final converted file
        meta_data = "#Created on: " + str(date.today()) +"\n#Input file is: " + str(inputfile) +"\n#Output file is: " + str(outputfile) +"\n#Material file is: " + str(materialfile) +"\n#Attenuation model is: " + str(atten_model) +"\n#Grain size is: " + str(grain_size) +" mm"+"\n#Oscillation period is: " + str(oscillation_period) +" seconds\n#COH is: " + str(COH_val) +"H/10**6Si\n"

        out_save = np.loadtxt(str(outputdir)+'/Process_1_vel_converted.txt',comments='#')
        for i in range(1,no_of_processes):
            out = np.loadtxt(str(outputdir)+'/Process_'+str(i+1)+'_vel_converted.txt',comments='#') 
            out_save=np.append(out_save,out,axis=0)
        np.savetxt(str(outputdir)+'/'+str(outputfile),out_save,delimiter=',',header="#x(km), y(km), depth(km), Pressure(bar), Temperature(oC), Density(kg/m3), Vp(km/s), Vs(km/s), Vs_diff(%), Pseudo-melts(%)",comments=meta_data,fmt='%10.3f')
        
    else:
        print('\n###########################################')
        print('You did not provide required input.')
        print('Run the code with -h option for help.')
        print('###########################################\n')

if __name__ == '__main__':
    main(sys.argv[1:])