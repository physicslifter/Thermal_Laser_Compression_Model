#IMPORTS
from scipy import io
import matlab
from matplotlib import pyplot as plt
import os

#=============================================================
#FUNCTIONS
#=============================================================

#function for writing the input parameters (mainly copied from Create_default_mat.py)
def write_input_parameters(peak_temp, a, b, time_shift, diffusivity, filename):
    
    my_values={
        'a':a,
        'b':b,
    }

    #And saving the dictionary as a .mat file...
    #Finally, save the dictionary as a .mat file
    io.savemat('inputs/1D_input_matrix.mat',mdict=my_values)

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
    #np.save(out_filename,allow_pickle=True)

    #now return the output
    return mydata

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