#IMPORTS
from scipy import io
import matlab
from matplotlib import pyplot as plt
import os

#=============================================================
#FUNCTIONS
#=============================================================
def chi_sq_interp(x0,y0,x1, y1):
    #x0: the independent variable data for the model output
    #y0: the original dependent variable data for the model output
    #x1: independent variable data for results
    #y1: dependent variable data for results

    #Step 1 - interpolate the model output to match the dataset we have
    f=interpolate.interp1d(x0,y0, fill_value="extrapolate")
    interpolated_data=f(x1)

    #Now for the chi^2 assessment
    chi_sq=np.sum( np.divide(np.square(y0-y1) , y1))

    return chi_sq

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

def run_model(peak_temp, a, b, time_shift, diffusivity, BC_filename, out_filename):

    #write the input parameters for the run file
    write_input_parameters(peak_temp, a, b, time_shift, diffusivity, BC_filename)

    #Set up and run the matlab engine
    eng=matlab.engine.start_matlab()

    #run the file with no output arguments
    eng.heat_equation_data_1D_exp(nargout=0)

    #once the file runs, go get the output
    mydata=io.loadmat('FEM_output/1D_run_output.mat')

    #save as numpy file .npy
    np.save(out_filename,allow_pickle=True)

    #And finally, get the chi^2 value
    chi_2=get_chi_sq_from_output_mat('FEM_output/1D_run_output.mat')

    #now return the output
    return [mydata, chi_2]

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