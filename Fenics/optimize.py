import wrapped_run
import numpy as np
from run_model import main as run
from plot import plot
from scipy.optimize import minimize
from scipy.optimize import differential_evolution
import json

my_shots=(
    's88773',
)

num_steps=60
experimental_data=wrapped_run.load_json('data_dict.json')

def optimizable_function(parameters):
    a, b, peak_temp, start_time = parameters
    model_data=run(my_shots, num_steps, a, b, peak_temp, start_time)
    chi_2=wrapped_run.get_chi_sq(model_data, experimental_data)
    output_file='optimization_output.csv'
    model_results_file='model_output.csv'
    
    #save output for later inspection
    with open(output_file, 'a+') as file_object:
        num_iterations=sum(1 for line in open(output_file))-1
        file_object.write('\n')#newline
        file_object.write(str(num_iterations)+','+str(np.sum(chi_2))+','+str(a)+','+str(b)+','+str(peak_temp)+','+str(start_time))
    
    '''
    model_results_json=json.dumps(model_data)
    f=open(model_results_file, 'w')
    f.write(model_results_json)
    f.close()
    '''
    
    #return our fitness function values
    return np.sum(chi_2)
    

initial_parameters=[0.1, 60, 42000, 1.5*10**-8] #for minimize method
bounds=((0.01,0.05),(5,60),(20000,30000),(1.5*10**-8,2.5*10**-8)) #for differential_evolution method
#minimization_optimization=minimize(optimizable_function, initial_parameters, method='Nelder-Mead', options={'maxiter':100, 'disp':True})
DE_optimization=differential_evolution(optimizable_function,((0.01,0.05),(5,60),(25000,50000),(1.75*10**-8,2.05*10**-8)), popsize=10, disp=True)