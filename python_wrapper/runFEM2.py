#IMPORTS
from scipy import io
import matlab
from matplotlib import pyplot as plt

#=============================================================
#FUNCTIONS
#=============================================================

#function for writing the input parameters (mainly copied from Create_default_mat.py)
def write_input_parameters(peak_temp):
    #Default values...
    sop_lineouts_data='../Data/s88773_sop_lineouts_TP.xlsx' #.xlsx file with the sop lineout data
    reflectivity_data='../Data/Reflectivity_s88773.xlsx' #.xlsx file containing the reflectivity data to use
    throughput=1.0 #throughput correction based on slit width and shot number
    aeta= 25.5#1/delta_t of SOP data
    a0= 481000.0# based on ND filter
    t0a= 1.909# based on ND filter
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

    #Writing to a dictionary
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
        'peak_temp':peak_temp
    }

    #save the matrix in the default location
    io.savemat('inputs/input_matrix.mat',mdict=default_values)

#start the matlab engine
#eng=matlab.engine.start_matlab()

#write an input file with the peak temp as 10000K
#write_input_parameters(10000)

#read back the data

def run_model(input_parameters):

    #write the input parameters
    write_input_parameters(input_parameters)

    #Set up and run the matlab engine
    eng=matlab.engine.start_matlab()

    #run the file with no output arguments
    eng.runFEM2(nargout=0)

    #once the file runs, go get the output
    mydata=io.loadmat('FEM_output/run_output.mat')

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
    plt.title('Model output: Peak temp = '+str(model_output['peak_temp'])+'K')

    #And finally, show the plot