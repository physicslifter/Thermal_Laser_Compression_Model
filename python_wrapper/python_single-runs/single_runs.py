# functions for running a single run

#==========================================
#IMPORTS
from scipy import io
from scipy import interpolate
import matlab
from matplotlib import pyplot as plt
import os
import numpy as np
from matlab import engine
#==========================================

#=================================================================================
#_________________________________________________________________________________
#Writing input parameters...
def write_input_parameters(peak_temp, a, b, time_shift, run_ID):
    
    my_values={
        'a':a,
        'b':b,
    }

    filename='BC_'+run_ID #define name of the boundary condition file 

    #And saving the dictionary as a .mat file...
    #Finally, save the dictionary as a .mat file
    save_name=run_ID+'_input.mat'
    io.savemat(save_name,mdict=my_values)

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

#___________________________________________________________________________________
#Getting least_squares

#generic least squares function
def least_squares_interp(x0,y0,x1,y1):
    f=interpolate.interp1d(x0,y0, bounds_error=False, fill_value='extrapolate')
    interpolated_data=f(x1)

    least_squares=np.sum(np.square(y1-interpolated_data))

    return least_squares

def get_least_sq(run_ID):
    short_id=run_ID[-2:] #get last 2 numbers for the run
    filepath=short_id+'_output.mat'

    #now read in the data
    d=io.loadmat(filepath)

    x0=d['tlist'][0] #tlist data for all runs @ all faces

    x1=d['t_data1'] #time data
    x2=d['t_data2']
    x3=d['t_data3']
    y0_1=d['T11'][0] #face 1 model output
    y0_2=d['T12'][0] # face 2 model output
    y0_3=d['T13'][0] # face 3 model output
    y1_1=d['temp1'] # face 1 experimental data
    y1_2=d['temp2'] # face 2 experimental data
    y1_3=d['temp3'] # face 3 experimental data

    least_squares1=least_squares_interp(x0,y0_1,x1,y1_1)
    least_squares2=least_squares_interp(x0,y0_2,x2,y1_2)
    least_squares3=least_squares_interp(x0,y0_3,x3,y1_3)

    return [least_squares1,least_squares2,least_squares3]


#=====================================================================
#Running the model
#=====================================================================
def run_model(peak_temp, a, b, time_shift, run_ID):
    #write the input parameters for the run file
    write_input_parameters(peak_temp, a, b, time_shift, run_ID)

    #Set up and run the matlab engine
    eng=matlab.engine.start_matlab()

    #run the file with no output arguments
    if run_ID=="s88773":
        #run the associated file for s88773
        eng.s88773(nargout=0)
    else:
        print('run ID not found')

    short_id=run_ID[-2:] #get last 2 numbers for the run
    output_filename=short_id+'_output.mat'

    #Finally, get the least_squares value
    ls=get_least_sq(run_ID)

    #now return the output
    return ls

#======================================================================
#Writing a simple version for sending to the optimizer
#======================================================================
def simple_run(parameter_array, run_ID:str, face1:bool, face2:bool, face3:bool):
    #takes single array as input
    # array values are defined as below:

    #parameter_array[0]=peak_temp
    # parameter_array[1]=a
    # parameter_array[2]=b
    # parameter_array[3]=time_shift

    #run the model
    run=run_model(parameter_array[0],parameter_array[1],parameter_array[2],parameter_array[3], run_ID)

    #defining sls
    if face1==True:
        if face2==True:
            if face3==True:
                sls=run
            else:
                sls=run[:2]
        else:
            if face3==True:
                sls=[run[0],run[2]]
            else:
                sls=[run[0]]        
    else:
        if face2==True:
            if face3==True:
                sls=[run[1],run[2]]
            else:
                sls=[run[1]]
        else:
            if face3==True:
                sls=[run[2]]
            else:
                sls=[]      

    sls=np.sum(sls)

    #print values to the terminal
    print(sls, parameter_array)

    #Get filename of the data
    short_id=run_ID[-2:] #get last 2 numbers for the run
    fname=short_id+'_data.csv'

    #write values to the data file
    with open(fname, 'a+') as file_object:
        num_iterations=sum(1 for line in open(fname))-1
        file_object.write('\n')#newline
        file_object.write(str(num_iterations)+', '+str(sls)+', '+str(parameter_array[0])+', '+str(parameter_array[1])+', '+str(parameter_array[2])+', '+str(parameter_array[3]))

    return sls

def simple_run73(parameter_array):
    return simple_run(parameter_array, 's88773', True, True, True)