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

#Constants for defining the thermal conductivity
a=0.028 #a (in k=a*x+b)
b=30.0 #b (in k=a*x+b)

#=============================================================================
#Constants from the BC_external_exp function
#=============================================================================
peak_temp=35000.0
diffusivity=36000000.0
time_shift=1.8*10.0**(-8.0)

#Writing the values to a dictionary
my_values={
    'sop_lineouts_data':sop_lineouts_data,
    'reflectivity_data':reflectivity_data,
    'a':a,
    'b':b,
    'peak_temp':peak_temp,
    'diffusivity':diffusivity,
    'time_shift':time_shift
}

#And saving the dictionary as a .mat file...
#Finally, save the dictionary as a .mat file
io.savemat('inputs/1D_input_matrix.mat',mdict=my_values)
