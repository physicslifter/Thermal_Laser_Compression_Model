import json
#import optimize

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
    
def run_optimization(num_steps, my_shots, optimization_name, bounds, popsize):
    rewrite_global_params(num_steps, my_shots, optimization_name, bounds, popsize) #rewrite params to global json
    import optimize
    optimize.run_optimization() #Run the optimization
    
    
    
    