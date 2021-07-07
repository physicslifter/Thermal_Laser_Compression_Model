#A file for writing the .mat input for the FEM in python.
#The subsequent step to this in the process will be to send the input to the MATLAB file
#to run in the python matlab engine
from scipy import io

#===============================
#Writing the Default values
#===============================

#write the matrix variables (a list of the variables can be found in )
#===============================================================================
#Initial data refinement
#===============================================================================
sop_lineouts_data='../Data/s88773_sop_lineouts_TP.xlsx' #.xlsx file with the sop lineout data
reflectivity_data='../Data/Reflectivity_s88773.xlsx' #.xlsx file containing the reflectivity data to use
throughput=1.0 #throughput correction based on slit width and shot number
aeta= 25.5#1/delta_t of SOP data
a0= 481000.0# based on ND filter
t0a= 1.909# based on ND filter

#===============================================================================
#FEM
#===============================================================================

#sample
density=12500.0 #density kg/m^3 previous at 11000
heat_capacity=500.0 #heat capacity J/kg K
t0=2000.0 #initial temp K/counts
width=1.0 #width (um)

#window
km=100.0
dm=4200.0
cm=800.0
t0m=2000.0

#Transient Heated Block variables
peak_temp=28000
start_time=4.0*10.0**(-8.0)

#other variables
a=0.028 #a (in k=a*x+b)
b=30 #b (in k=a*x+b)
#==================================================

#Now write a dictionary for the default values
default_values={
    'sop_lineouts_data':sop_lineouts_data,
    'reflectivity_data':reflectivity_data,
    'xw':throughput,
    'aeta':aeta,
    'a0':a0,
    't0a':t0a,
    'd':density,
    'c':heat_capacity,
    't0':t0,
    'w':width,
    'km':km,
    'dm':dm,
    'cm':cm,
    't0m':t0m,
    'peak_temp':peak_temp,
    'start_time':start_time,
    'a':a,
    'b':b
}

#Finally, save the dictionary as a .mat file
io.savemat('inputs/default_input_matrix.mat',mdict=default_values)


