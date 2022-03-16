# V2RhoT_gibbs
Python scripts to convert seismic velocities to temperature and density based on Gibbs free energy minimization from Perple_X

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
