import wrapped_run
import numpy as np
from run_model import main as run
from scipy.optimize import minimize
from scipy.optimize import differential_evolution
import os
import errno

my_shots=(
    's88773',
)

global_vars=wrapped_run.load_json('global_variables.json')
my_shots=global_vars['my_shots']
bounds=tuple(global_vars['bounds'])
num_steps=int(global_vars['num_steps'])
optimization_name=global_vars['optimization_name']
popsize=global_vars['popsize']
root_path=global_vars['optimization_data_path']

#make sure this name is not already in use. If it is, an error will be raised
try:
    os.makedirs(root_path+'/'+optimization_name)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise


num_steps=60
experimental_data=wrapped_run.load_json('data_dict.json')

def optimizable_function(parameters):
    a, b, peak_temp, start_time = parameters
    model_data=run(my_shots, num_steps, a, b, peak_temp, start_time)
    chi_2=wrapped_run.get_chi_sq(model_data, experimental_data)
    output_file=root_path+'/'+optimization_name+'/optimization_output.csv'
    
    #save output for later inspection
    with open(output_file, 'a+') as file_object:
        num_iterations=sum(1 for line in open(output_file))-1
        file_object.write('\n')#newline
        file_object.write(str(num_iterations)+','+str(np.sum(chi_2))+','+str(a)+','+str(b)+','+str(peak_temp)+','+str(start_time))
    
    #return our fitness function values
    return np.sum(chi_2)
    

initial_parameters=[0.1, 60, 42000, 1.5*10**-8] #for minimize method
#bounds=((0.01,0.05),(5,60),(20000,30000),(1.5*10**-8,2.5*10**-8)) #for differential_evolution method
#minimization_optimization=minimize(optimizable_function, initial_parameters, method='Nelder-Mead', options={'maxiter':100, 'disp':True})
DE_optimization=differential_evolution(optimizable_function,bounds, popsize=popsize, disp=True)