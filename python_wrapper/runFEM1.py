#IMPORTS
from scipy import io
import matlab
from matplotlib import pyplot as plt

#=============================================================
#FUNCTIONS
#=============================================================

#function for writing the input parameters (mainly copied from Create_default_mat.py)
def write_input_parameters(peak_temp, a, b, time_shift, diffusivity):
    
    my_values={
        'a':a,
        'b':b,
        'peak_temp':peak_temp,
        'diffusivity':diffusivity,
        'time_shift':time_shift
    }

    #And saving the dictionary as a .mat file...
    #Finally, save the dictionary as a .mat file
    io.savemat('inputs/1D_input_matrix.mat',mdict=my_values)

def run_model(peak_temp, a, b, time_shift, diffusivity):

    #write the input parameters
    write_input_parameters(peak_temp, a, b, time_shift, diffusivity)

    #Set up and run the matlab engine
    eng=matlab.engine.start_matlab()

    #run the file with no output arguments
    eng.runFEM2(nargout=0)

    #once the file runs, go get the output
    mydata=io.loadmat('FEM_output/1D_run_output.mat')

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