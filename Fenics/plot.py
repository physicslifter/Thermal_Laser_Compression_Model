import matplotlib
matplotlib.use('WebAgg')
from matplotlib import pyplot as plt

#Plotting model output
def plot(dict):
    #Takes a dictionary type object as input (the type
    # output by our run_model function)
    plt.clf()
    keys=dict.keys()
    num_runs=len(keys)
    nrows=int((num_runs+num_runs%2)/2)
    ncols=int((num_runs-num_runs%2)/2)
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
    
#Plotting model output with experimental data

    