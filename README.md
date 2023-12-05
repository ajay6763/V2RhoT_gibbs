### V2RhoT_gibbs

[V2RhoT_gibbs](https://doi.org/10.5281/zenodo.6538257) is a Python tool that allows converting seismic velocities to temperature and density in the mantle using a Gibbs free energy minimization algorithm.

It uses pre-computed look-up tables of density and anharmonic seismic velocities as a function of pressure and temperature for a given bulk mantle composition in terms of major oxides weight % within the Na2O-CaO-FeO-MgO-Al2O3-SiO2 (NCFMAS) system since they account for the ~99% of the Earth’s mantle ([Palme and O’Neill, 2013](#ONeill_Palme_2003)).
These look-up tables are computed using a module of **LitMod2D_2.0** called **Generator** ([Afonso et al. 2008](#Afonso_etal_2008); [Kumar et al. 2020](#Kumar_etal_2020)]) that uses **Perple_X** ([Connolly 2005](#Connolly_2005), [Connolly 2009](#Connolly_2009) for Gibbs free-energy minimization.
Look-up tables for the standard mantle compositions are provided with the distribution in databases folder.
In case you want to generate look-up tables of your own for a given mantle composition then follow the instructions in [Installing Generator](#installing-generator-module). Generator is a module of **LitMod2D_2.0**.
You can choose to do the density and seismic velocities calculations directly in [Perple_X](https://www.perplex.ethz.ch/).
The only requirement is that the look-up tables should have the following column order:
```
T(k) P(bar) Density(kg/m3) Vp(km/s) Vs(km/s)
```
### Citing
If you use this tool please cite: 
1. Kumar, A., Cacace, M., Scheck-Wenderoth, M., Götze, H.-J., & Kaus, B. J. P. (2022). Present-day upper-mantle architecture of the Alps: Insights from data-driven dynamic modeling. Geophysical Research Letters, 49, e2022GL099476. https://doi.org/10.1029/2022GL099476
2. Kumar, A., Fernàndez, M.,Jiménez‐Munt, I., Torne, M., Vergés, J.,& Afonso, J. C. (2020). LitMod2D_2.0:An improved integratedgeophysical‐petrological modeling toolfor the physical interpretation of uppermantle anomalies.Geochemistry,Geophysics, Geosystems,21,e2019GC008777. https://doi.org/10.1029/2019GC008777
#### Installation

**V2RhoT_gibbs** can be downloaded from [Zenodo](https://doi.org/10.5281/zenodo.6538257) or from the [GitHub](https://github.com/ajay6763/V2RhoT_gibbs).
It is compatible with Python 3 and requires the following standard Python libraries:
- multiprocessing
- numpy
- scipy
- matplotlib       

#### Running

The main script to run the conversions is `parallel_conversion.py`.
To get info about the input data format and options run the following in the **V2RhoT_gibbs** folder.
```
python parallel_conversion.py -h
 ```
This should give you the following:
```powershell
PS C:\Users\kumar\work\V2RhoT_Gibbs_projects\V2RhoT_gibbs_dev> python parallel_conversion.py -h

###########################################################################################
 Simple run example: parallel_conversion.py -i <inputfile> -o <outputfile>
             Below are the available options which you can pass as a command line argument:
            -h : help
            -p : no of parts to run in parallel (e.g., no of available cores). Default is 1
            -I : input directory. Default is ./data_tomo
            -i : input file name in the data_tomo folder. Format x(*) y(*) depth(km) Vs(km/s)
            -O : output directory. Default is ./output
            -o : output file name which will be saved in output folder
            -m : name of the material file in databases folder. Default is DMM_HP
            -A : attenuation model flag: 1 for Jackson and Faul 2010, 2 for Behn et al., 2009. If you choose 2 then you will have to supply COH. See COH option.
            -g : grain size in mm. Default is 10mm
            -s : oscillation period in seconds. Default is 75 seconds
            -W : Water contnent in  H/10**6Si. Defualt is 50 H/10**6Si which almost dry. Note this will be used if you choose attenuation model 2.
            -geo : Geology file in the data_geo folder. Format Format x(*) y(*) geo(codes). Note: This option is not active at the momemnt.
All of the input options will be written in the output file as comments.
```

#### Example
Here, we will run a sample conversion for the data provided `Test_parallel.dat` in the `data_tomo` folder.
Input tomography data must be in this folder.
The format of input file is as follows:
```
x(*) y(*) depth(km) Vs(km/s)
```  
The first two columns are location e.g., long, lat and they can have any unit you want (km, m, degrees).
The depth and seismic velocity must have units in km and km/s, respectively.

To run the conversion on this file navigate to the `V2RhoT_gibbs` folder and run the following command:
```powershell
PS C:\Users\kumar\work\V2RhoT_Gibbs_projects\V2RhoT_gibbs_opt\V2RhoT_gibbs> python parallel_conversion.py -p 1 -i Test_parallel.dat -o Test_parallel.csv
```
This should give you following output:
```
PS C:\Users\kumar\work\V2RhoT_Gibbs_projects\V2RhoT_gibbs_opt\V2RhoT_gibbs> python .\parallel_conversion.py -p 1 -i Test_parallel.dat -o Test_parallel.csv

###########################################
Your options are:
Total no of processes 1
Input file is Test_parallel.dat
Output file is Test_parallel.csv
Material file is DMM_HP
Input folder is ./data_tomo
Output folder is ./output
Attenuation model is 1
Grain size is 10.0 mm
Oscillation period is 75.0 seconds
COH is (used when Attenuation model is 1) 50.0 H/10**6Si

###########################################

###########################################
Output directory exists. It will be overwritten.

###########################################
Chose 1 Atten model.

###########################################
Time spent on loading data: 2.35 sec

###########################################

###########################################
Length of the input model: 3120
Length of chunk to run on each processor: 3120

###########################################
Started worker Process_1
Process_1 worker finished in 0.27479236125946044 min.

###########################################
Total execution time: 0.3 min
Hope that was fast enough for you. If not then maybe try increasing no of processes with -p option
Enjoy your conversion results in the output folder :)
###########################################
```
After a successful run, you should have a file named `Test_parallel.csv` in the output folder that has the following header:
```
#Created on: 2023-12-05
#Input file is: Test_parallel.dat
#Output file is: Test_parallel.csv
#Material file is: DMM_HP
#Attenuation model is: 1
#Grain size is: 10.0 mm
#Oscillation period is: 75.0 seconds
#COH is: 50.0H/10**6Si
#x(km), y(km), depth(km), Pressure(bar), Temperature(oC), Density(kg/m3), Vp(km/s), Vs(km/s), Vs_diff(%), Pseudo-melts(%)						
```
In this example, since we did not specify grain size, oscillation period or a material file it took the default values that are printed at the top. Two type of attenuation models are implemented which are input with a option -A where 1 refers to model of Jackson and Faul, 2010 10.1016/j.pepi.2010.09.005 (default case) and 2 refers to model of Behn et al., 2009 10.1016/j.epsl.2009.03.014. 
Now let us say you want to used different grain size, oscillation period and material file then you can do the following:
```powershell
PS C:\Users\kumar\work\V2RhoT_Gibbs_projects\V2RhoT_gibbs_opt\V2RhoT_gibbs> python parallel_conversion.py -p 6 -i Test_parallel.dat -o Test_parallel_user_options.csv -g 5.0 -s 50.0 -m DMM_7km_HP
```
Options that we are using in the above example are:
- `-p 6` : means it will split the input file into 6 parts and run the conversion in parallel. In case you have large tomography file then you can increase it.
- `-i Test_parallel.dat` : this is same as previous example because we wanted to run the conversion with different parameters.
- `-o Test_parallel_user_options.csv` : for this we have added _user_option and it will produce a different file while keeping the output from the previous example. This allows you to have different output files for different input parameters.
- `-g 5.0` : here we are using grain size of 5 mm.
- `-s 50.0` : here we are using oscillation period of 50 seconds
- `-m DMM_7km_HP` : here we are using a material file called `DMM_7km_HP` which is provided in the databases folder.
- `-A 2` : here we are using attenuation model of Behn et al., 2009 10.1016/j.epsl.2009.03.014.
- `-W 1000.0` : here we water content of 1000.0 H/10**6Si which is input for attenuation model 2. The defualt value is 50.0 which equivalant to dry conditions.
This should produce the following output on the screen
```
###########################################
Your options are:
Total no of processes 6
Input file is Test_parallel.dat
Output file is Test_parallel_user_options.csv
Material file is DMM_7km_HP
Input folder is ./data_tomo
Output folder is ./output
Attenuation model is 2
Grain size is 5.0 mm
Oscillation period is 50.0 seconds
COH is (used when Attenuation model is 1) 1000.0 H/10**6Si

###########################################

###########################################
Output directory exists. It will be overwritten.

###########################################
Chose 2 Atten model.

###########################################
Time spent on loading data: 2.53 sec

###########################################

###########################################
Length of the input model: 3120
Length of chunk to run on each processor: 520

###########################################
Started worker Process_1
Started worker Process_2
Process_1 worker finished in 0.04423450231552124 min.
Started worker Process_3
Process_2 worker finished in 0.045877655347188316 min.
Started worker Process_4
Process_3 worker finished in 0.05230891704559326 min.
Started worker Process_5
Process_4 worker finished in 0.049238912264506024 min.
Started worker Process_6
Process_5 worker finished in 0.04650991757710775 min.
Process_6 worker finished in 0.0446491281191508 min.

###########################################
Total execution time: 0.2 min
Hope that was fast enough for you. If not then maybe try increasing no of processes with -p option
Enjoy your conversion results in the output folder :)

###########################################
```

#### Generator Module

#### Provided look-up tables
| **Name** | **SiO2** | **Al2O3** | **FeO** | **MgO** | **CaO** | **Na2O** | **#Mg** | **Reference**                                    |
|----------|----------|-----------|---------|---------|---------|----------|---------|--------------------------------------------------|
| PUM      | 45.0     | 4.50      | 8.10    | 37.80   | 3.60    | 0.360    | 89.3    | ([McDonough and Sun, 1995](#McDonough_Sun_1995)) |
| DMM      | 44.70    | 3.98      | 8.18    | 38.73   | 3.17    | 0.130    | 89.4    | ([Workman and Hart, 2005](#Workman_Hart_2005))   |
| DMM_7km  | 44.43    | 2.97      | 8.23    | 40.78   | 2.70    | 0.045    | 89.8    | ([Kumar et al., 2021](#Kumar_etal_2021))         |
| Pyrolite | 45.1     | 4.6       | 7.6     | 38.1    | 3.1     | 0.40     | 89.9    | ([Ringwood, 1975](#Ringwood_1975))               |
| Tc_1     | 44.5     | 3.5       | 8.0     | 39.8    | 3.1     | 0.24     | 89.8    | ([Griffin et al., 2009](#Griffin_etal_2009))     |

#### Thermodynamic database key:

| Table name    | Description |
|---------------|-------------|
|   _HP      |Augmented‐modified version of [Holland and Powell, 1998](#Holland_Powell_1998) thermodynamic database [Afonso & Zlotnik, 2011](#Afonso_Zlotnik_2011))      |
|   _3       |[Xu et al., 2008](#Xu_etal_2008)| 

##### Installing Generator Module

Generator module is an interface that talks to the **Perple_X** executable by basically producing input files for the gibbs-free energy minimization.

Depending on your OS (Windows, Linux or Mac) you can choose the way of installing the Generator which is explained in [here](https://raw.githubusercontent.com/ajay6763/LitMod2D_2.0_package_dist_users/master/manual/LitMod_GUI_manual.pdf).

**Windows:**\
For this, I recommend using the docker way of installation.
Just briefly, this will have the whole `LitMod2D_2.0` installed and we will only use the Generator part.

Introduction to docker is the following meme and should give you an idea of what it is (from Reddit, do not know how to cite memes )

![image](./figures/V2RhoT_gibbs/meme.png)

Docker can be installed from the official website (https://www.docker.com/products/docker-desktop).
Once you have docker installed and working properly, you need to download the image of LitMod2D 2.0 from the docker hub.
To do this run the following command in a terminal (for Linux or Mac OS) and command prompt (Windows OS):
```
docker pull ajay6763/litmod2d 2.0:final
```
To download it might take a while. Once the downloading finishes check if you have the image downloaded by running following
in powershell to check list of images:
```powershell
docker images
```
for example
```
REPOSITORY              TAG       IMAGE ID       CREATED         SIZE
ajay6763/litmod2d_2.0   final     d97d8869dfee   18 months ago   4.31GB
```
Docker image of the LitMod2D 2.0 behaves as a stand-alone distribution (like a virtual machine). In case you want to see the graphical contents, we need to attach a server to the docker image in order to see the graphical components. This is only important if you decide to use the LitMod2D_2.0 part, not for the Generator part. To be exhaustive server attachment is kept. You can decide to skip the IP and DISPLAY variable part. From this point onward, make sure that you have Docker Desktop program running.

Get IP
```powershell
Get-NetIPAddress -AddressFamily IPv4 -InterfaceIndex $(Get-NetConnectionProfile | Select-Object -ExpandProperty InterfaceIndex) | Select-Object -ExpandProperty IPAddress
```
which results in e.g. `139.17.91.159`.

The use it to set `DISPLAY` variable:
```powershell
set-variable -name DISPLAY 139.17.91.159:0.0
```
Before running the docker image we need to map the local drive (<host directory> e.g., C:\Users\kumar\work) and path inside the container (<container directory> e.g., /home/work). This allows you to save the files from inside the container to your local machine. You do this by moving the files that you want to save to the <container directory> folder.
Finally, run docker with a proper paths ID
```powershell
docker run -ti --rm -v <host directory>:<container directory> -e DISPLAY=$DISPLAY <docker image ID>
```
e.g.
```powershell
docker run -ti --rm -v C:/Users/kumar/work/litmod_docker:/home/work -e DISPLAY=$DISPLAY d97d8869dfee
```
This will take you to the GUI folder. Navigate back to the Generator_Linux folder and follow the steps in Running Generator Module.


**Linux:**\
You do not need to have the whole installation of the LitMod2D_2.0 just the Generator part (written below). Download the [Generator_Linux](https://github.com/ajay6763/LitMod2D_2.0_package_dist_users/tree/master/Generator_Linux). To run this program in Linux based platform you need to install "wine". This allows you to run executable files in Linux.

To install wine run following commands in the terminal
```
sudo apt-get install wine-binfmt "sudo update-binfmts --import /usr/share/binfmts/wine
```
Put all the .exe files from the Generator_Liunx in your path through ~/.bashsrc file and source it. You have to make Generator_LINUX executable by changing the permissions by running the following command:
```bash
chmod 755 Generator_LINUX
```

**Mac OS:**\
*(do not have access to this OS so not sure if the Linux way would work. The docker way should work)*


##### Running Generator Module

To run Generator module run following in the Generator_Linux folder:
```
./Generator_LINUX
```
After typing enter a command prompt will appear and you will be asked certain questions that you answer and thermo-physical property tables will be put in the `Build_dat_(name of the run you gave)`.
Inside this folder there will be a file called `TABLE_LITMOD_no_atten` copy this file into the database folder and rename such that you can identify it e.g. 98,97 etc.
In case you are doing this the docker way make sure to copy this renamed file to `<container directory>` e.g., `/home/work`.


##### Example Generator
Here, we will run an example of generating material file.
We start the docker and open a powershell. In the powershell we start the docker as explained above and you should be in the GUI folder as shown in the image below:

![image.png](./figures/V2RhoT_gibbs/image.png)

Note here that the `"C:/Users/kumar/work/litmod_docker"` is the location on my local computer which I am  mapping to `"/home/work"` inside the docker.
Make sure that when you are done you move file to this location.

Navigate back to the Generator Linux folder and following command:

```
./Generator_Linux
```
This will produce following in the terminal:

![image-1.png](./figures/V2RhoT_gibbs/image-1.png)

Then you choose a thermodynamic database, e.g. 5
[Attach image]

Next, it will ask you is this will the the sublithospheric mantle which simple means a prescribed P-T range to do calculations.
Choosing 0 should suffice for the lithosheric mantle depths.
Then you enter major oxides weight %.

![image-2.png](./figures/V2RhoT_gibbs/image-2.png)

Name of the body is the the calculation name (e.g. test) which will be the folder (Buildfile.dat_test)) where output will be saved.
[Attach image]

The it will ask if you want the full property table with assemblages as well or just the properties.
The difference between the two is that in case of full calculation in addition to properties table it will also output stable phase and mineral assemblages that can be used to visualize the stable minerals. For this example we choose 0  
![image-4.png](./figures/V2RhoT_gibbs/image-4.png)

It will a while to do the calculations.
In the end it will ask for grain size and oscillation period. You can 10 mm and 75 seconds.  
![image-5.png](./figures/V2RhoT_gibbs/image-5.png)

Once the calculations are finished you will see a Buildfile.dat_test folder. Go to this folder and copy the `TABLE_LITMOD_no_atten` file to the `/home/work`.
This is file the final material file that can be used in velocity conversions.

**Recommendation:** it is better to copy the whole folder and have it so that later you have information of the minimization.
Once you copy this to `/home/work` is should be visible at `C:/Users/kumar/work/litmod_docker`.

![image-6.png](./figures/V2RhoT_gibbs/image-6.png)

___
## References
<a name='Afonso_etal_2008'></a>
Afonso, J. C., Fernandez, M., Ranalli, G., Griffin, W., & Connolly, J. (2008). Integrated geophysical‐petrological modeling of the lithosphere and sublithospheric upper mantle: Methodology and applications. Geochemistry, Geophysics, Geosystems, 9(5). https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2007GC001834 

<a name='Afonso_Zlotnik_2011'></a>
Afonso, J. C., & Zlotnik, S. (2011). The subductability of the continental lithosphere: The before and after story. In D. Brown & P. D. Ryan (Eds.), Arc‐Continent Collision (pp. 53–86). Berlin: Frontiers in Earth Sciences, Springer.

<a name='Connolly_2005'></a>
Connolly, J. A. (2005). Computation of phase equilibria by linear programming: a tool for geodynamic modeling and its application to subduction zone decarbonation. Earth and Planetary Science Letters, 236(1-2), 524-541. https://www.sciencedirect.com/science/article/pii/S0012821X05002839 

<a name='Connolly_2009'></a>
Connolly, J. (2009). The geodynamic equation of state: what and how. Geochemistry, Geophysics, Geosystems, 10(10). https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2009GC002540 

<a name='Griffin_etal_2009'></a>
Griffin, W., O'Reilly, S., Afonso, J., & Begg, G. (2009). The composition and evolution of lithospheric mantle: A re-evaluation and its tecton- ic implications. Journal of Petrology, 50(7), 1185–1204. https://doi.org/10.1093/petrology/egn033

<a name='Hollan_Powell_1998'></a>
Holland, T., & Powell, R. (1998). An internally consistent thermodynamic data set for phases of petrological interest. Journal of Metamorphic Geology, 16(3), 309–343.

<a name='Kumar_etal_2021'></a>
Kumar, A., Fernàndez, M., Vergés, J., Torne, M., & Jiménez-Munt, I. (2021). Opposite symmetry in the lithospheric structure of the Alboran and Algerian basins and their margins (Western Mediterranean): Geodynamic implications. Journal of Geophysical Research: Solid Earth, 126, e2020JB021388. https://doi. org/10.1029/2020JB021388

<a name='Kumar_etal_2020'></a>
Kumar, A., Fernàndez, M., Jiménez‐Munt, I., Torne, M., Vergés, J., & Afonso, J. C. (2020). LitMod2D_2. 0: An improved integrated geophysical‐petrological modeling tool for the physical interpretation of upper mantle anomalies. Geochemistry, Geophysics, Geosystems, 21(3), e2019GC008777. https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019GC008777 

<a name='McDonough_Sun_1995'></a>
McDonough, W., & Sun, S.‐S. (1995). The composition of the Earth'. Chemical Geology, 120(120), 223–253. https://www.geol.umd.edu/~mcdonoug/KITP%20Website%20for%20Bill/papers/Mantle-Comp/Workman_Hart_(EPSL_05).pdf 

<a name='ONeill_Palme_2003'></a>
O’Neill, H , and Palme HSC. 2003. '3.1 Cosmochemical Estimates of Mantle Composition', Treatise on geochemistry, 2: 1-38. https://www.geol.umd.edu/~mcdonoug/KITP%20Website%20for%20Bill/papers/Earth_Models/3.1%20Palme%20&%20O'Neill%20Primative%20mantle%20(1).pdf 

<a name='Ringwood_1975'></a>
Ringwood, A.E. (1975) Composition and Petrology of the Earth’s Mantle. McGraw-Hill, London, New York and Sydney, 618 p.

<a name='Workman_Hart_2005'></a>
Workman, R., & Hart, S. (2005). Major and trace element composition of the depleted MORB mantle (DMM). Earth and Planetary Science Letters, 231(1‐2), 53–72. https://doi.org/10.1016/0009-2541(94)00140-4 

<a name='Xu_etal_2008'></a>
Xu, W., Lithgow‐Bertelloni, C., Stixrude, L., & Ritsema, J. (2008). The effect of bulk composition and temperature on mantle seismic structure. Earth and Planetary Science Letters, 275(1‐2), 70–79.
___
[Back to main contents](./README.md)




















# V2RhoT_gibbs
Python scripts to convert seismic velocities to temperature and density based on Gibbs free energy minimization from Perple_X.

Look-up tables are generated using a module of LitMod2D_2.0 (https://github.com/ajay6763/LitMod2D_2.0_package_dist_users) called Generator. If you want to generate look-up tables
for new composition please follow the instruction in LitMod2D_2.0. You do not need to have the whole installation of the LitMod2D_2.0 just the Generator part (written below).
Instructions to install can be found at https://github.com/ajay6763/LitMod2D_2.0_package_dist_users/blob/master/manual/LitMod_GUI_manual.pdf. 

Program to produce thermophysical properties look-table used in LitMod2D_2.0

** To run this program in Linux based platform you need to install "wine" . This allows you to run *.exe files.

To install wine run following commands in the terminal

"sudo apt-get install wine-binfmt "sudo update-binfmts --import /usr/share/binfmts/wine"

Put all the .exe file in your path throught ~/.bashsrc file

After you have done as suggested above, you run "Generator_LINUX" from the same directory. You have to make it executable by changing the permissions by running following command.

"chmod 755 Generator_LINUX"

then "./Generator_LINUX"

After type enter command a command prompt will appear and you will be asked certain question that you answer and thermophysical property tables will be put in the Build_dat_(name of the run you gave). Inside this folder there will be a file called TABLE_LITMOD_no_atten copy this file into the database folder and rename such that you can identify it e.g. 98,97 etc. 


If you do not want to waste time in installling then I highly recommend using Docker way described in the manual (https://github.com/ajay6763/LitMod2D_2.0_package_dist_users/blob/master/manual/LitMod_GUI_manual.pdf). 
)

If you use V2RhoT_gibbs please cite LitMod2D, Afonso et al., 2008 (https://doi.org/10.1029/2007GC001834), Kumar et al., 2020( https://doi.org/10.1029/2019GC008777) and Perple_X Connolly, 2009 (https://doi.org/10.1029/2009GC002540).
