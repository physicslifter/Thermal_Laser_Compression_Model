import wrapped_run
import numpy as np
from run_model import main as run
from plot import plot
from scipy.optimize import differential_evolution
import json

my_shots=(
    's88773',
    's88776',
    's88780',
    's86483'
)


num_steps=60
experimental_data=wrapped_run.load_json('data_dict.json')

def optimizable_function(parameters, shots, num_steps):
    a, b, peak_temp, start_time = parameters
    wrapped_run.load_json('data_dict.json')
    model_data=run(shots, num_steps, a, b, peak_temp, start_time)
    chi_2=wrapped_run.get_chi_sq(model_data, experimental_data)
    output_file='optimization_output.csv'
    model_results_file='model_output.csv'
    
    #save output for later inspection
    with open(output_file, 'a+') as file_object:
        num_iterations=sum(1 for line in open(output_file))-1
        file_object.write('\n')#newline
        file_object.write(str(num_iterations)+','+str(chi_2))
    
    model_results_json=json.dumps(model_data)
    f=open(model_results_file, 'w')
    f.write(model_results_json)
    f.close()
    
    #return our fitness function values
    return np.sum(chi_2)

bounds=((0.01,0.05),(5,100),(20000,30000),(1.5*10**-8,2.5*10**-8))
optimization=differential_evolution(func=optimizable_function, bounds=bounds, args=(my_shots, 60), popsize=5, disp=True)
