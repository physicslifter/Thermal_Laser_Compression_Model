import matplotlib
matplotlib.use('WebAgg')
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from scipy import io
from plot import plot
import json

def load_json(fname):
    with open(fname) as json_file:
        data=json.load(json_file)
    return data

model_data_filepath='model_results.json'
real_data_filepath='data_dict.json'

#Plot setup
model_data=load_json(model_data_filepath)
keys=model_data.keys()
num_runs=len(keys)
nrows=int((num_runs+num_runs%2)/2)
ncols=int((num_runs-num_runs%2)/2)+1
real_data=load_json(real_data_filepath)

fig1, ax1 = plt.subplots(nrows,ncols)

#Animation
def animate(i):
    pullData=load_json(model_data_filepath)
    for x in range(num_runs):
        run_ID=str(list(keys)[x])
        c=x%ncols
        r=int((x-c)/2)
        faces=list(pullData[run_ID])
        for face in faces:
            print(r,c)
            ax1[r,c].clear()
            ax1[r,c].plot(pullData[run_ID][face][0],pullData[run_ID][face][1],label=face)
            ax1[r,c].plot(real_data[run_ID][face][0],real_data[run_ID][face][1])
            ax1[r,c].set_title(run_ID)
        plt.legend(loc="upper left")
        
ani=animation.FuncAnimation(fig1, animate, interval=10000)
plt.show()