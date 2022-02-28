import matplotlib
matplotlib.use('WebAgg')
from matplotlib import pyplot as plt
import json

def load_json(fname):
    with open(fname) as json_file:
        data=json.load(json_file)
        
    return data

#Plotting model output
def plot(dict):
    #Takes a dictionary type object as input (the type
    # output by our run_model function)
    plt.clf()
    keys=dict.keys()
    num_runs=len(keys)
    nrows=int((num_runs+num_runs%2)/2)
    ncols=int((num_runs-num_runs%2)/2)+1
    if ncols==0:
        name=list(dict.keys())[0]
        faces=list(dict[name])
        for face in faces:
            plt.plot(dict[name][face][0],dict[name][face][1], label=face)
        plt.legend(loc='upper left')
    
    else:
        fig, ax = plt.subplots(nrows,ncols)

        for x in range(num_runs):
            run_ID=str(list(keys)[x])
            c=x%ncols
            r=int((x-c)/2)
            faces=list(dict[run_ID])
            for face in faces:
                print(r,c)
                ax[r,c].plot(dict[run_ID][face][0],dict[run_ID][face][1],label=face)
                ax[r,c].set_title(run_ID)

            plt.legend(loc="upper left")
            
    plt.show()
    
def plot_json(json="model_results.json"):
    plot(load_json(json))
    
#Plotting model output with experimental data
def comparePlots():
    model_output=load_json('model_results.json')
    real_data=load_json('data_dict.json')
    plt.clf()
    keys=model_output.keys()
    num_runs=len(keys)
    nrows=int((num_runs+num_runs%2)/2)
    ncols=int((num_runs-num_runs%2)/2)+1

    name=list(model_output.keys())[0]
    faces=list(model_output[name])
    for face in faces:
        plt.scatter(real_data[name][face][0], real_data[name][face][1])
        plt.plot(model_output[name][face][0], model_output[name][face][1], label=face)
        
    plt.xlim(1.5*10**-8,4*10**-8)
    plt.ylim(0, 30000)
    plt.show()
    
def comparePlotsByName(runID):
    model_output=load_json('model_results.json')
    real_data=load_json('data_dict.json')
    plt.clf()
    keys=model_output.keys()
    num_runs=len(keys)
    nrows=int((num_runs+num_runs%2)/2)
    ncols=int((num_runs-num_runs%2)/2)+1

    name=runID
    faces=list(model_output[name])
    for face in faces:
        plt.scatter(real_data[name][face][0], real_data[name][face][1])
        plt.plot(model_output[name][face][0], model_output[name][face][1], label=face)
        
    plt.xlim(1.5*10**-8,4*10**-8)
    plt.ylim(0, 30000)
    plt.show()

    