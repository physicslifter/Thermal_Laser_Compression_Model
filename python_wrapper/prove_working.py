#Import the matlab engine for python
import matlab.engine

#Write function for running the model
def fem_run(input_matrix,run_name):
    eng=matlab.engine.start_matlab()

    #Run the matlab function with the associated inputs
    eng.RunFEM(input_matrix,run_name,nargout=0)

    #exit/quit the engine
    eng.quit()

    #show user the run is complete
    print('run complete')

#Point to input data, and the sample run name
input_data_file='inputs/default_input_matrix.mat'
sample_run_name='sample_script_in_python-default_output_matrix'

#Run the model
fem_run(input_data_file,sample_run_name)

