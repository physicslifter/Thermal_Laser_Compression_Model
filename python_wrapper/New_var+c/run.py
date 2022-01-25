#IMPORTS
from scipy import io
from scipy import interpolate
import matlab
from matplotlib import pyplot as plt
import os
import numpy as np
from matlab import engine
#==============================================================
#==============================================================
#Matlab stuff
#==============================================================
#==============================================================

#=============================================================
#FUNCTIONS
#=============================================================

def least_squares_interp(x0,y0,x1,y1):
    f=interpolate.interp1d(x0,y0, bounds_error=False, fill_value='extrapolate')
    interpolated_data=f(x1)

    least_squares=np.sum(np.square(y1-interpolated_data))

    return least_squares

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

#Run 1
#================================================================================================
def write_sqwv_input_parameters1(peak_temp, a, b, c, time_shift, filename):
    
    my_values={
        'a':a,
        'b':b,
        'c':c,
    }
    #And saving the dictionary as a .mat file...
    #Finally, save the dictionary as a .mat file
    io.savemat('1D_combined_input_matrix1.mat',mdict=my_values)
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

#Running the optimization
def run41(peak_temp, a, b, c, time_shift):
    BC_filename='BC1'
    #write the input parameters for the run file
    write_sqwv_input_parameters1(peak_temp, a, b, c, time_shift, BC_filename)
    #Set up and run the matlab engine
    eng=matlab.engine.start_matlab() 
    #run the file with no output arguments
    eng.Combo_eq_1(nargout=0)
    #once the file runs, go get the output
    mydata=io.loadmat('combo1.mat')
    #save as numpy file .npy
    #np.save(out_filename,allow_pickle=True)
    #And finally, get the chi^2 value
    #chi_2=get_chi_sq_from_output_mat('FEM_output/1D_combo_run_output.mat')
    chi_2=get_least_sq_from_output_mat_run4('combo1.mat')
    #now return the output
    return [mydata, chi_2]
    
#WRAPPED FUNCTION
#(this is what we optimize)
def sfr1(parameter_array):
    #Takes a single parameter array as input
    #The parameter array is defined as below
    # parameter_array[0]=peak_temp
    # parameter_array[1]=a
    # parameter_array[2]=b
    # parameter_array[3]=c
    # parameter_array[4]=time_shift
    #run the model & return single chi^2 value

    #Run the model
    run=run41(parameter_array[0],parameter_array[1],parameter_array[2],parameter_array[3], parameter_array[4])
    
    #Get Chi^2
    chi_2=run[1]

    #print save the value(s) to the optimization data file
    sls=np.sum(chi_2)
    print(sls)
    with open('data1.csv', 'a+') as file_object:
        num_iterations=sum(1 for line in open('data1.csv'))-1
        file_object.write(str('\n'))#newline
        file_object.write(str(num_iterations)+', '+str(sls)+', '+str(parameter_array[0])+', '+str(parameter_array[1])+', '+str(parameter_array[2])+', '+str(parameter_array[3])+', '+str(parameter_array[4]))
    return sls