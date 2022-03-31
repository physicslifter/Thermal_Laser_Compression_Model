import json
from run_model import main as run
from scipy.optimize import minimize
from scipy.optimize import differential_evolution
import wrapped_run
import errno
import os
import numpy as np


def load_json(fname):
    with open(fname) as json_file:
        data=json.load(json_file)
        
    return data

def rewrite_global_params(num_steps, my_shots, optimization_name, bounds, popsize):
    global_json_as_dict=load_json('global_variables.json')
    global_json_as_dict['num_steps']=num_steps
    global_json_as_dict['my_shots']=my_shots
    global_json_as_dict['optimization_name']=optimization_name
    global_json_as_dict['bounds']=bounds
    global_json_as_dict['popsize']=popsize
    
    global_json=json.dumps(global_json_as_dict)
    f=open('global_variables.json', 'w')
    f.write(global_json)
    f.close()
    
def optimizable_function(parameters):
    #globals
    experimental_data=wrapped_run.load_json('data_dict.json')
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
    a, b, start_time, peak_temp = parameters
    model_data=run(my_shots, num_steps, a, b, peak_temp, start_time)
    chi_2=wrapped_run.get_chi_sq(model_data, experimental_data)
    output_file=root_path+'/'+optimization_name+'/optimization_output.csv'
    
    with open(output_file, 'a+') as file_object:
        num_iterations=sum(1 for line in open(output_file))-1
        file_object.write('\n')#newline
        file_object.write(str(num_iterations)+','+str(np.sum(chi_2))+','+str(a)+','+str(b)+','+str(start_time)+','+str(peak_temp))
    
    #return our fitness function values
    return np.sum(chi_2)
    
def run_optimization(num_steps:int, my_shots:list, optimization_name:str, bounds:tuple, popsize:int):
    rewrite_global_params(num_steps, my_shots, optimization_name, bounds, popsize) #rewrite params to global json
    #import optimize
    #optimize.run_optimization() #Run the optimization
    
    #globals
    experimental_data=wrapped_run.load_json('data_dict.json')
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
    
    diffEV=differential_evolution(optimizable_function,bounds, popsize=popsize, disp=True)
    