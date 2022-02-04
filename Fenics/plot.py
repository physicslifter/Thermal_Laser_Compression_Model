import matplotlib
matplotlib.use('WebAgg')
from matplotlib import pyplot as plt

def plot(dict):
    #Takes a dictionary type object as input (the type
    # output by our run_model function)
    keys=dict.keys()
    num_runs=len(keys)
    nrows=int((num_runs+num_runs%2)/2)
    ncols=int((num_runs-num_runs%2)/2)
    
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
    