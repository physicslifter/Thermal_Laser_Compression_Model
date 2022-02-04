import func1
import json

run_simple_fem=func1.main

def load_json(fname):
    with open(fname) as json_file:
        data=json.load(json_file)
        
    return data

def convert_geometry(length, density=12800, ambient_density=7900):
    ratio=ambient_density/density
    new_length=ratio*length
    
    return new_length

def main(run_IDs, num_steps, a, b, peak_temp):
    num_runs=len(run_IDs)
    geometry_dict=load_json('geometry_dict.json')
    values_dict=load_json('default_values_dict.json')
    results={}
    for x in range(num_runs):
        runID=run_IDs[x] # get the xth name of the run specified in run_IDs input
        solution={
            'face1': run_simple_fem(tf=values_dict['tf'], num_steps=num_steps,a=a, b=b, peak_temp=peak_temp, init_temp=values_dict['init_temp'], k_1=values_dict['k_1'], rho=values_dict['rho'], c=values_dict['c'], Fe_length=convert_geometry(geometry_dict[runID]['face1']) , MgO_length=convert_geometry(values_dict['MgO_length'])),
            'face2': run_simple_fem(tf=values_dict['tf'], num_steps=num_steps,a=a, b=b, peak_temp=peak_temp, init_temp=values_dict['init_temp'], k_1=values_dict['k_1'], rho=values_dict['rho'], c=values_dict['c'], Fe_length=convert_geometry(geometry_dict[runID]['face2']) , MgO_length=convert_geometry(values_dict['MgO_length'])),
            'face3': run_simple_fem(tf=values_dict['tf'], num_steps=num_steps,a=a, b=b, peak_temp=peak_temp, init_temp=values_dict['init_temp'], k_1=values_dict['k_1'], rho=values_dict['rho'], c=values_dict['c'], Fe_length=convert_geometry(geometry_dict[runID]['face3']) , MgO_length=convert_geometry(values_dict['MgO_length']))
        }
        results[runID]=solution
    return results

if __name__ == "__main__":
    main()