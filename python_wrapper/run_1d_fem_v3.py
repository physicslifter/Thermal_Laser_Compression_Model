#IMPORTS
from scipy import io
from scipy import interpolate
import matlab
from matplotlib import pyplot as plt
import os
import numpy as np
from matlab import engine

#=============================================================
#FUNCTIONS
#=============================================================
def chi_sq_interp(x0,y0,x1, y1):
    #x0: the independent variable data for the model output
    #y0: the dependent variable data for the model output
    #x1: independent variable data for results
    #y1: dependent variable data for results

    #Step 1 - interpolate the model output to match the dataset we have
    f=interpolate.interp1d(x0,y0, bounds_error=False, fill_value='extrapolate')
    interpolated_data=f(x1)

    avg=np.average(y1)
    variance=(np.sum(np.square(y1-avg)))/(len(y1))

    #Now for the chi^2 assessment
    chi_sq=np.sum( np.divide(np.square(interpolated_data-y1) , variance))/(len(y1))

    return chi_sq

def least_squares_interp(x0,y0,x1,y1):
    f=interpolate.interp1d(x0,y0, bounds_error=False, fill_value='extrapolate')
    interpolated_data=f(x1)

    least_squares=np.sum(np.square(y1-interpolated_data))

    return least_squares

def chi_sq_interp_peak(x0,y0,x1, y1):
    #This function should return the same as chi_sq_interp, except it only evaluates the chi^2 before the peak
    
    #step 1 - find the peak
    peak_val=np.where(y1 == np.amax(y1))[0].item()

    #Chop off the experimental data at the peak
    new_x1=x1[0:peak_val]
    new_y1=y1[0:peak_val]

    #Now return the chi^2
    #Since the chi^2 has a built in interpolation, we can leave the model ranges untouched, and their data will automatically match the domain of the new experimental data
    return chi_sq_interp(x0,y0,new_x1,new_y1)

def get_chi_sq_from_output_mat(filepath):

    #Start by reading in the data
    d=io.loadmat(filepath)

    #write to variables so chi^sq_interp will accept the arguments
    #variables are named as follows:
    #[input name to chi_sq_interp (x0,y0,x1 or y1)]_[face # (1, 2, or 3)]_[run #: 73, 76 or 80]
    #x0 is only defined once, since it will always be tlist
    #the t_data arrays are only defined once per experiment run
    
    x0=d['tlist'][0] #tlist data for all runs @ all faces

    #run s88773
    x1_73=d['t_data_73'] #time data
    y0_1_73=d['T11'][0] #face 1 model output
    y0_2_73=d['T12'][0] # face 2 model output
    y0_3_73=d['T13'][0] # face 3 model output
    y1_1_73=d['temp1_73'] # face 1 experimental data
    y1_2_73=d['temp2_73'] # face 2 experimental data
    y1_3_73=d['temp3_73'] # face 3 experimental data

    #run s88776
    x1_76=d['t_data_76'] #time data
    y0_1_76=d['T21'][0] #face 1 model output
    y0_2_76=d['T22'][0] # face 2 model output
    y0_3_76=d['T23'][0] # face 3 model output
    y1_1_76=d['temp1_76'] # face 1 experimental data
    y1_2_76=d['temp2_76'] # face 2 experimental data
    y1_3_76=d['temp3_76'] # face 3 experimental data

    #run s88780
    x1_80=d['t_data_80'] #time data
    y0_1_80=d['T31'][0] #face 1 model output
    y0_2_80=d['T32'][0] # face 2 model output
    y0_3_80=d['T33'][0] # face 3 model output
    y1_1_80=d['temp1_80'] # face 1 experimental data
    y1_2_80=d['temp2_80'] # face 2 experimental data
    y1_3_80=d['temp3_80'] # face 3 experimental data

    #Get chi^2 values for each run
    cs731=chi_sq_interp(x0,y0_1_73,x1_73,y1_1_73) # s88773 face 1
    cs732=chi_sq_interp(x0,y0_2_73,x1_73,y1_2_73) # s88773 face 2
    cs733=chi_sq_interp(x0,y0_3_73,x1_73,y1_3_73) # s88773 face 3
    cs761=chi_sq_interp(x0,y0_1_76,x1_76,y1_1_76) # s88773 face 1
    cs762=chi_sq_interp(x0,y0_2_76,x1_76,y1_2_76) # s88773 face 2
    cs763=chi_sq_interp(x0,y0_3_76,x1_76,y1_3_76) # s88773 face 3
    cs801=chi_sq_interp(x0,y0_1_80,x1_80,y1_1_80) # s88773 face 1
    cs802=chi_sq_interp(x0,y0_2_80,x1_80,y1_2_80) # s88773 face 2
    cs803=chi_sq_interp(x0,y0_3_80,x1_80,y1_3_80) # s88773 face 3

    #finally, return the values as an array
    return [cs731,cs732,cs733,cs761,cs762,cs763,cs801,cs802,cs803]

def peak_chi_sq_from_output_mat(filepath):

    #Start by reading in the data
    d=io.loadmat(filepath)

    

    #write to variables so chi^sq_interp will accept the arguments
    #variables are named as follows:
    #[input name to chi_sq_interp (x0,y0,x1 or y1)]_[face # (1, 2, or 3)]_[run #: 73, 76 or 80]
    #x0 is only defined once, since it will always be tlist
    #the t_data arrays are only defined once per experiment run
    
    x0=d['tlist'][0] #tlist data for all runs @ all faces

    #run s88773
    x1_73=d['t_data_73'] #time data
    y0_1_73=d['T11'][0] #face 1 model output
    y0_2_73=d['T12'][0] # face 2 model output
    y0_3_73=d['T13'][0] # face 3 model output
    y1_1_73=d['temp1_73'] # face 1 experimental data
    y1_2_73=d['temp2_73'] # face 2 experimental data
    y1_3_73=d['temp3_73'] # face 3 experimental data

    #run s88776
    x1_76=d['t_data_76'] #time data
    y0_1_76=d['T21'][0] #face 1 model output
    y0_2_76=d['T22'][0] # face 2 model output
    y0_3_76=d['T23'][0] # face 3 model output
    y1_1_76=d['temp1_76'] # face 1 experimental data
    y1_2_76=d['temp2_76'] # face 2 experimental data
    y1_3_76=d['temp3_76'] # face 3 experimental data

    #run s88780
    x1_80=d['t_data_80'] #time data
    y0_1_80=d['T31'][0] #face 1 model output
    y0_2_80=d['T32'][0] # face 2 model output
    y0_3_80=d['T33'][0] # face 3 model output
    y1_1_80=d['temp1_80'] # face 1 experimental data
    y1_2_80=d['temp2_80'] # face 2 experimental data
    y1_3_80=d['temp3_80'] # face 3 experimental data

    #Get chi^2 values for each run
    cs731=chi_sq_interp_peak(x0,y0_1_73,x1_73,y1_1_73) # s88773 face 1
    cs732=chi_sq_interp_peak(x0,y0_2_73,x1_73,y1_2_73) # s88773 face 2
    cs733=chi_sq_interp_peak(x0,y0_3_73,x1_73,y1_3_73) # s88773 face 3
    cs761=chi_sq_interp_peak(x0,y0_1_76,x1_76,y1_1_76) # s88773 face 1
    cs762=chi_sq_interp_peak(x0,y0_2_76,x1_76,y1_2_76) # s88773 face 2
    cs763=chi_sq_interp_peak(x0,y0_3_76,x1_76,y1_3_76) # s88773 face 3
    cs801=chi_sq_interp_peak(x0,y0_1_80,x1_80,y1_1_80) # s88773 face 1
    cs802=chi_sq_interp_peak(x0,y0_2_80,x1_80,y1_2_80) # s88773 face 2
    cs803=chi_sq_interp_peak(x0,y0_3_80,x1_80,y1_3_80) # s88773 face 3

    #finally, return the values as an array
    return [cs731,cs732,cs733,cs761,cs762,cs763,cs801,cs802,cs803]

#function for the least_squares fitting test from early August 2021
def s88780test_get_least_sq_from_output_mat(filepath):

    #Start by reading in the data
    d=io.loadmat(filepath)

    #We only care about the data from s88780 for this one
    #This varies from the above functions because we need to read in all of the time data
    x0=d['tlist'][0] #tlist data for all runs @ all faces

    #run s88780
    x1_1_80=d['t_data1_80'] #time data
    x1_2_80=d['t_data2_80'] #time data
    x1_3_80=d['t_data3_80'] #time data
    y0_1_80=d['T31'][0] #face 1 model output
    y0_2_80=d['T32'][0] # face 2 model output
    y0_3_80=d['T33'][0] # face 3 model output
    y1_1_80=d['temp1_80'] # face 1 experimental data
    y1_2_80=d['temp2_80'] # face 2 experimental data
    y1_3_80=d['temp3_80'] # face 3 experimental data

    #Get chi^2 values for each run
    cs801=least_squares_interp(x0,y0_1_80,x1_1_80,y1_1_80) # s88773 face 1
    cs802=least_squares_interp(x0,y0_2_80,x1_2_80,y1_2_80) # s88773 face 2
    cs803=least_squares_interp(x0,y0_3_80,x1_3_80,y1_3_80) # s88773 face 3

    #finally, return the values as an array
    return [cs801,cs802,cs803]

def new_get_least_sq_from_output_mat(filepath):

    #Start by reading in the data
    d=io.loadmat(filepath)

    #We only care about the data from s88780 for this one
    #This varies from the above functions because we need to read in all of the time data
    x0=d['tlist'][0] #tlist data for all runs @ all faces

    #run s88773
    x1_1_73=d['t_data1_73'] #time data
    x1_2_73=d['t_data2_73'] #time data
    x1_3_73=d['t_data3_73'] #time data
    y0_1_73=d['T11'][0] #face 1 model output
    y0_2_73=d['T12'][0] # face 2 model output
    y0_3_73=d['T13'][0] # face 3 model output
    y1_1_73=d['temp1_73'] # face 1 experimental data
    y1_2_73=d['temp2_73'] # face 2 experimental data
    y1_3_73=d['temp3_73'] # face 3 experimental data

    #run s88776
    x1_1_76=d['t_data1_76'] #time data
    x1_2_76=d['t_data2_76'] #time data
    x1_3_76=d['t_data3_76'] #time data
    y0_1_76=d['T21'][0] #face 1 model output
    y0_2_76=d['T22'][0] # face 2 model output
    y0_3_76=d['T23'][0] # face 3 model output
    y1_1_76=d['temp1_76'] # face 1 experimental data
    y1_2_80=d['temp2_76'] # face 2 experimental data
    y1_3_80=d['temp3_76'] # face 3 experimental data

    #run s88780
    x1_1_80=d['t_data1_80'] #time data
    x1_2_80=d['t_data2_80'] #time data
    x1_3_80=d['t_data3_80'] #time data
    y0_1_80=d['T31'][0] #face 1 model output
    y0_2_80=d['T32'][0] # face 2 model output
    y0_3_80=d['T33'][0] # face 3 model output
    y1_1_80=d['temp1_80'] # face 1 experimental data
    y1_2_80=d['temp2_80'] # face 2 experimental data
    y1_3_80=d['temp3_80'] # face 3 experimental data

    #Get chi^2 values for each run
    cs801=least_squares_interp(x0,y0_1_80,x1_1_80,y1_1_80) # s88780 face 1
    cs802=least_squares_interp(x0,y0_2_80,x1_2_80,y1_2_80) # s8878- face 2
    cs803=least_squares_interp(x0,y0_3_80,x1_3_80,y1_3_80) # s88780 face 3

    #finally, return the values as an array
    return [cs801,cs802,cs803]

#function for getting the least_squares value from the output matrix
def get_least_sq_from_output_mat(filepath):

    #Start by reading in the data
    d=io.loadmat(filepath)

    #We only care about the data from s88780 for this one
    #This varies from the above functions because we need to read in all of the time data
    x0=d['tlist'][0] #tlist data for all runs @ all faces

    #run s88773
    x1_1_73=d['t_data1_73'] #time data
    x1_2_73=d['t_data2_73'] #time data
    x1_3_73=d['t_data3_73'] #time data
    y0_1_73=d['T11'][0] #face 1 model output
    y0_2_73=d['T12'][0] # face 2 model output
    y0_3_73=d['T13'][0] # face 3 model output
    y1_1_73=d['temp1_73'] # face 1 experimental data
    y1_2_73=d['temp2_73'] # face 2 experimental data
    y1_3_73=d['temp3_73'] # face 3 experimental data

    #run s88776
    x1_1_76=d['t_data1_76'] #time data
    x1_2_76=d['t_data2_76'] #time data
    x1_3_76=d['t_data3_76'] #time data
    y0_1_76=d['T21'][0] #face 1 model output
    y0_2_76=d['T22'][0] # face 2 model output
    y0_3_76=d['T23'][0] # face 3 model output
    y1_1_76=d['temp1_76'] # face 1 experimental data
    y1_2_76=d['temp2_76'] # face 2 experimental data
    y1_3_76=d['temp3_76'] # face 3 experimental data

    #run s88780
    x1_1_80=d['t_data1_80'] #time data
    x1_2_80=d['t_data2_80'] #time data
    x1_3_80=d['t_data3_80'] #time data
    y0_1_80=d['T31'][0] #face 1 model output
    y0_2_80=d['T32'][0] # face 2 model output
    y0_3_80=d['T33'][0] # face 3 model output
    y1_1_80=d['temp1_80'] # face 1 experimental data
    y1_2_80=d['temp2_80'] # face 2 experimental data
    y1_3_80=d['temp3_80'] # face 3 experimental data

    #Get chi^2 values for each run

    #73
    cs731=least_squares_interp(x0,y0_1_73,x1_1_73,y1_1_73) # s88773 face 1
    cs732=least_squares_interp(x0,y0_2_73,x1_2_73,y1_2_73) # s88773 face 2
    cs733=least_squares_interp(x0,y0_3_73,x1_3_73,y1_3_73) # s88773 face 3

    #76
    cs761=least_squares_interp(x0,y0_1_76,x1_1_76,y1_1_76) # s88776 face 1
    cs762=least_squares_interp(x0,y0_2_76,x1_2_76,y1_2_76) # s88776 face 2
    cs763=least_squares_interp(x0,y0_3_76,x1_3_76,y1_3_76) # s88776 face 3

    #80
    cs801=least_squares_interp(x0,y0_1_80,x1_1_80,y1_1_80) # s88780 face 1
    cs802=least_squares_interp(x0,y0_2_80,x1_2_80,y1_2_80) # s88780 face 2
    cs803=least_squares_interp(x0,y0_3_80,x1_3_80,y1_3_80) # s88780 face 3

    #finally, return the values as an array
    return [cs731,cs732,cs733,cs761,cs762,cs763,cs801,cs802,cs803]

#function for getting the least_squares value from the output matrix
def get_least_sq_from_output_mat_run4(filepath):

    #Start by reading in the data
    d=io.loadmat(filepath)

    #We only care about the data from s88780 for this one
    #This varies from the above functions because we need to read in all of the time data
    x0=d['tlist'][0] #tlist data for all runs @ all faces

    #run s88773
    x1_1_73=d['t_data1_73'] #time data
    x1_2_73=d['t_data2_73'] #time data
    x1_3_73=d['t_data3_73'] #time data
    y0_1_73=d['T11'][0] #face 1 model output
    y0_2_73=d['T12'][0] # face 2 model output
    y0_3_73=d['T13'][0] # face 3 model output
    y1_1_73=d['temp1_73'] # face 1 experimental data
    y1_2_73=d['temp2_73'] # face 2 experimental data
    y1_3_73=d['temp3_73'] # face 3 experimental data

    #run s88776
    x1_1_76=d['t_data1_76'] #time data
    x1_2_76=d['t_data2_76'] #time data
    x1_3_76=d['t_data3_76'] #time data
    y0_1_76=d['T21'][0] #face 1 model output
    y0_2_76=d['T22'][0] # face 2 model output
    y0_3_76=d['T23'][0] # face 3 model output
    y1_1_76=d['temp1_76'] # face 1 experimental data
    y1_2_76=d['temp2_76'] # face 2 experimental data
    y1_3_76=d['temp3_76'] # face 3 experimental data

    #run s88780
    x1_1_80=d['t_data1_80'] #time data
    x1_2_80=d['t_data2_80'] #time data
    x1_3_80=d['t_data3_80'] #time data
    y0_1_80=d['T31'][0] #face 1 model output
    y0_2_80=d['T32'][0] # face 2 model output
    y0_3_80=d['T33'][0] # face 3 model output
    y1_1_80=d['temp1_80'] # face 1 experimental data
    y1_2_80=d['temp2_80'] # face 2 experimental data
    y1_3_80=d['temp3_80'] # face 3 experimental data

    #run s86483
    x1_1_83=d['t_data1_83'] #time data
    x1_2_83=d['t_data2_83'] #time data
    x1_3_83=d['t_data3_83'] #time data
    y0_1_83=d['T41'][0] #face 1 model output
    y0_2_83=d['T42'][0] # face 2 model output
    y0_3_83=d['T43'][0] # face 3 model output
    y1_1_83=d['temp1_83'] # face 1 experimental data
    y1_2_83=d['temp2_83'] # face 2 experimental data
    y1_3_83=d['temp3_83'] # face 3 experimental data

    #Get chi^2 values for each run

    #73
    cs731=least_squares_interp(x0,y0_1_73,x1_1_73,y1_1_73) # s88773 face 1
    cs732=least_squares_interp(x0,y0_2_73,x1_2_73,y1_2_73) # s88773 face 2
    cs733=least_squares_interp(x0,y0_3_73,x1_3_73,y1_3_73) # s88773 face 3

    #76
    cs761=least_squares_interp(x0,y0_1_76,x1_1_76,y1_1_76) # s88776 face 1
    cs762=least_squares_interp(x0,y0_2_76,x1_2_76,y1_2_76) # s88776 face 2
    cs763=least_squares_interp(x0,y0_3_76,x1_3_76,y1_3_76) # s88776 face 3

    #80
    cs801=least_squares_interp(x0,y0_1_80,x1_1_80,y1_1_80) # s88780 face 1
    cs802=least_squares_interp(x0,y0_2_80,x1_2_80,y1_2_80) # s88780 face 2
    cs803=least_squares_interp(x0,y0_3_80,x1_3_80,y1_3_80) # s88780 face 3

    #83
    cs831=least_squares_interp(x0,y0_1_83,x1_1_83,y1_1_83) # s88780 face 1
    cs832=least_squares_interp(x0,y0_2_83,x1_2_83,y1_2_83) # s88780 face 2
    cs833=least_squares_interp(x0,y0_3_83,x1_3_83,y1_3_83) # s88780 face 3

    #finally, return the values as an array
    return [cs731,cs732,cs733,cs761,cs762,cs763,cs801,cs802,cs803,cs831,cs832,cs833]

#function for writing the input parameters (mainly copied from Create_default_mat.py)
def write_input_parameters(peak_temp, a, b, time_shift, diffusivity, filename):
    
    my_values={
        'a':a,
        'b':b,
    }

    #And saving the dictionary as a .mat file...
    #Finally, save the dictionary as a .mat file
    io.savemat('inputs/1D_combined_input_matrix.mat',mdict=my_values)

    #Now for rewriting the BC file...

    #Changing the name of the file to .txt so we can edit
    os.system('cmd /c "ren '+filename+'.m '+filename+'.txt"')

    f=open(filename+'.txt') #open file
    s=f.readlines() #read file contents into a list
    f.close() #close the file

    #change lines accordingly
    s[27]='  '+'peak='+str(peak_temp)+'; \n'#peak_temp
    s[28]='  '+'a='+str(diffusivity)+'; \n'#diffusivity
    s[29]='  '+'b='+str(time_shift)+'; \n'#time_shift

    #Now write the new lines to the .txt BC file
    with open(filename+'.txt','w') as filehandle:
        filehandle.writelines("%s" % line for line in s)

    #And finally, change the file back to .m
    os.system('cmd /c "ren '+filename+'.txt '+filename+'.m"')

def write_sqwv_input_parameters(peak_temp, a, b, time_shift, filename):
    
    my_values={
        'a':a,
        'b':b,
    }

    #And saving the dictionary as a .mat file...
    #Finally, save the dictionary as a .mat file
    io.savemat('inputs/1D_combined_input_matrix.mat',mdict=my_values)

    #Now for rewriting the BC file...

    #Changing the name of the file to .txt so we can edit
    os.system('cmd /c "ren '+filename+'.m '+filename+'.txt"')

    f=open(filename+'.txt') #open file
    s=f.readlines() #read file contents into a list
    f.close() #close the file

    #change lines accordingly
    s[20]='  '+'peak='+str(peak_temp)+'; \n'#peak_temp

    s[21]='  '+'start_time='+str(time_shift)+'; \n'#time_shift

    #Now write the new lines to the .txt BC file
    with open(filename+'.txt','w') as filehandle:
        filehandle.writelines("%s" % line for line in s)

    #And finally, change the file back to .m
    os.system('cmd /c "ren '+filename+'.txt '+filename+'.m"')


def run_model(peak_temp, a, b, time_shift, diffusivity):

    BC_filename='BC_external_exp'

    #write the input parameters for the run file
    write_input_parameters(peak_temp, a, b, time_shift, diffusivity, BC_filename)

    #Set up and run the matlab engine
    eng=matlab.engine.start_matlab()

    #run the file with no output arguments
    eng.heat_equation_data_1D_combo3(nargout=0)

    #once the file runs, go get the output
    mydata=io.loadmat('FEM_output/1D_combo_run_output.mat')

    #save as numpy file .npy
    #np.save(out_filename,allow_pickle=True)

    #And finally, get the chi^2 value
    #chi_2=get_chi_sq_from_output_mat('FEM_output/1D_combo_run_output.mat')
    chi_2=get_least_sq_from_output_mat('FEM_output/1D_combo_run_output.mat')

    #now return the output
    return [mydata, chi_2]

#function for running the model with the square wave
def run_sqwv_model(peak_temp, a, b, time_shift):
    BC_filename='BC_external_exp2'

    #write the input parameters for the run file
    write_sqwv_input_parameters(peak_temp, a, b, time_shift, BC_filename)

    #Set up and run the matlab engine
    eng=matlab.engine.start_matlab()

    #run the file with no output arguments
    eng.heat_equation_data_1D_combo5(nargout=0)

    #once the file runs, go get the output
    mydata=io.loadmat('FEM_output/1D_combo_run_output.mat')

    #save as numpy file .npy
    #np.save(out_filename,allow_pickle=True)

    #And finally, get the chi^2 value
    #chi_2=get_chi_sq_from_output_mat('FEM_output/1D_combo_run_output.mat')
    chi_2=get_least_sq_from_output_mat('FEM_output/1D_combo_run_output.mat')

    #now return the output
    return [mydata, chi_2]

def run_peak_model(peak_temp, a, b, time_shift, diffusivity):

    BC_filename='BC_external_exp'

    #write the input parameters for the run file
    write_input_parameters(peak_temp, a, b, time_shift, diffusivity, BC_filename)

    #Set up and run the matlab engine
    eng=matlab.engine.start_matlab()

    #run the file with no output arguments
    eng.heat_equation_data_1D_combo3(nargout=0)

    #once the file runs, go get the output
    mydata=io.loadmat('FEM_output/1D_combo_run_output.mat')

    #save as numpy file .npy
    #np.save(out_filename,allow_pickle=True)

    #And finally, get the chi^2 value
    chi_2=peak_chi_sq_from_output_mat('FEM_output/1D_combo_run_output.mat')

    #now return the output
    return [mydata, chi_2]

def run_sqwv_peak_model(peak_temp, a, b, time_shift):

    BC_filename='BC_external_exp2'

    #write the input parameters for the run file
    write_sqwv_input_parameters(peak_temp, a, b, time_shift, BC_filename)

    #Set up and run the matlab engine
    eng=matlab.engine.start_matlab()

    #run the file with no output arguments
    eng.heat_equation_data_1D_combo2(nargout=0)

    #once the file runs, go get the output
    mydata=io.loadmat('FEM_output/1D_combo_run_output.mat')

    #save as numpy file .npy
    #np.save(out_filename,allow_pickle=True)

    #And finally, get the chi^2 value
    chi_2=peak_chi_sq_from_output_mat('FEM_output/1D_combo_run_output.mat')

    #now return the output
    return [mydata, chi_2]

def run4(peak_temp, a, b, time_shift):
    BC_filename='BC_external_exp2'

    #write the input parameters for the run file
    write_sqwv_input_parameters(peak_temp, a, b, time_shift, BC_filename)

    #Set up and run the matlab engine
    eng=matlab.engine.start_matlab()

    #run the file with no output arguments
    eng.heat_equation_data_1D_combo7(nargout=0)

    #once the file runs, go get the output
    mydata=io.loadmat('FEM_output/1D_combo_run_output4.mat')

    #save as numpy file .npy
    #np.save(out_filename,allow_pickle=True)

    #And finally, get the chi^2 value
    #chi_2=get_chi_sq_from_output_mat('FEM_output/1D_combo_run_output.mat')
    chi_2=get_least_sq_from_output_mat('FEM_output/1D_combo_run_output4.mat')

    #now return the output
    return [mydata, chi_2]

#==================================================================
# Results plotting
#=================================================================

def plot_results(model_output):
    #corrected temps
    a1=model_output['temp1_corrected']
    a2=model_output['temp2_corrected']
    a3=model_output['temp3_corrected']

    ta=model_output['t_data']


    #model outputs
    b1=model_output['T11'][0] #boundary 1
    b2=model_output['T12'][0] #boundary 2
    b3=model_output['T13'][0] #boundary 3

    t=model_output['tlist'][0]

    #plotting

    #scatter the corrected temperature data
    plt.scatter(ta,a1)
    plt.scatter(ta,a2)
    plt.scatter(ta,a3)

    #plot the model output data
    plt.plot(t,b1,t,b2,t,b3)

    #assign an appropriate title to the plot
    plt.title('Model output')

#===========================================================
#Simple Runs
#===========================================================

def simple_run(parameter_array):

    #Takes a single parameter array as input
    #The parameter array is defined as below
    # parameter_array[0]=peak_temp
    # parameter_array[1]=a
    # parameter_array[2]=b
    # parameter_array[3]=time_shift
    # parameter_array[4]=diffusivity

    #run the model & return single chi^2 value
    return np.sum(run_model(parameter_array[0],parameter_array[1],parameter_array[2],parameter_array[3],parameter_array[4])[1])

def simple_peak_run(parameter_array):

    #Takes a single parameter array as input
    #The parameter array is defined as below
    # parameter_array[0]=peak_temp
    # parameter_array[1]=a
    # parameter_array[2]=b
    # parameter_array[3]=time_shift
    # parameter_array[4]=diffusivity

    #run the model & return single chi^2 value
    run=run_model(parameter_array[0],parameter_array[1],parameter_array[2],parameter_array[3],parameter_array[4])
    chi_2=run[1]
    print(np.sum(chi_2[0:3]))
    #save the value to the optimization data file
    with open('optimization_data.csv', 'a+') as file_object:
        num_iterations=sum(1 for line in open('optimization_data.csv'))-1
        file_object.write('\n')#newline
        file_object.write(str(num_iterations)+', '+str(np.sum(chi_2[0:3]))+', '+str(parameter_array[0])+', '+str(parameter_array[1])+', '+str(parameter_array[2])+', '+str(parameter_array[3])+', '+str(parameter_array[4]))

    return np.sum(chi_2[0:3])

def simple_sqwv_run(parameter_array):

    #Takes a single parameter array as input
    #The parameter array is defined as below
    # parameter_array[0]=peak_temp
    # parameter_array[1]=a
    # parameter_array[2]=b
    # parameter_array[3]=time_shift

    #run the model & return single chi^2 value
    run=run_sqwv_model(parameter_array[0],parameter_array[1],parameter_array[2],parameter_array[3])
    chi_2=run[1]
    print(np.sum(chi_2[0:3]))
    #save the value to the optimization data file
    #sls=chi_2[0]+chi_2[1]+chi_2[6]+chi_2[7]
    sls=np.sum(chi_2)
    with open('optimization_data.csv', 'a+') as file_object:
        num_iterations=sum(1 for line in open('optimization_data.csv'))-1
        file_object.write('\n')#newline
        file_object.write(str(num_iterations)+', '+str(sls)+', '+str(parameter_array[0])+', '+str(parameter_array[1])+', '+str(parameter_array[2])+', '+str(parameter_array[3]))

    return sls

def simple_four_run(parameter_array):

    #Takes a single parameter array as input
    #The parameter array is defined as below
    # parameter_array[0]=peak_temp
    # parameter_array[1]=a
    # parameter_array[2]=b
    # parameter_array[3]=time_shift

    #run the model & return single chi^2 value
    run=run4(parameter_array[0],parameter_array[1],parameter_array[2],parameter_array[3])
    chi_2=run[1]
    print(np.sum(chi_2))
    #save the value to the optimization data file
    #sls=chi_2[0]+chi_2[1]+chi_2[6]+chi_2[7]
    sls=np.sum(chi_2)
    with open('optimization_data.csv', 'a+') as file_object:
        num_iterations=sum(1 for line in open('optimization_data.csv'))-1
        file_object.write('\n')#newline
        file_object.write(str(num_iterations)+', '+str(sls)+', '+str(parameter_array[0])+', '+str(parameter_array[1])+', '+str(parameter_array[2])+', '+str(parameter_array[3]))

    return sls 