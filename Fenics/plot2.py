import matplotlib
matplotlib.use('WebAgg')
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import json
import math

def load_json(fname):
    with open(fname) as json_file:
        data=json.load(json_file)
        
    return data

class ComparitivePlot:
    def __init__(self, data_filepath):
        self.model_fp=data_filepath
        self.model_output=load_json(data_filepath)
        self.run_output=load_json('data_dict.json')
        self.run_IDs=list(self.model_output.keys())
        
    def update_model_to_current_run(self):
        self.model_output=load_json(self.model_fp)
        
    def show_runs(self):
        mystr=''
        for run in self.run_IDs:
            mystr=mystr+str(run)+', '
        print('Plot runs include: '+ mystr)
        
    def plot(self):
        plt.clf()
        num_runs=len(self.run_IDs)
        nrows=int(math.ceil(num_runs**0.5))
        ncols=nrows
        fig, ax = plt.subplots(nrows, ncols)
        print(nrows, ncols)
        row_count=0
        col_count=0
        fig_count=0
        for run in self.run_IDs:
            for face in list(self.run_output[run].keys()):
                ax[row_count,col_count].plot(self.run_output[run][face][0], self.run_output[run][face][1])
                ax[row_count,col_count].plot(self.model_output[run][face][0], self.model_output[run][face][1])
                ax[row_count,col_count].set_title(run)
                print(run, face, fig_count, row_count, col_count)
            fig_count=fig_count+1
            col_count=fig_count%nrows
            row_count=int(fig_count/nrows)
        
        plt.show()
        
    def animate_plots(self):
        num_runs=len(self.run_IDs)
        ncols=int(math.ceil(num_runs**0.5))
        nrows=int(math.ceil(num_runs/ncols))
        fig, ax = plt.subplots(nrows, ncols)
        
        def animate(i,nrows=nrows):
            if nrows>1:
                self.update_model_to_current_run() #Get current model output
                row_count=0
                col_count=0
                fig_count=0
                for run in self.run_IDs:
                    ax[row_count, col_count].clear() #clear data from previous plot
                    for face in list(self.run_output[run].keys()):
                        print(run, face)
                        ax[row_count,col_count].plot(self.run_output[run][face][0], self.run_output[run][face][1])
                        ax[row_count,col_count].plot(self.model_output[run][face][0], self.model_output[run][face][1])
                        ax[row_count,col_count].set_title(run)
                    fig_count=fig_count+1
                    col_count=fig_count%nrows
                    row_count=int(fig_count/nrows)
            else:
                self.update_model_to_current_run() #Get current model output
                fig_count=0
                for run in self.run_IDs:
                    ax[fig_count].clear() #clear data from previous plot
                    for face in list(self.run_output[run].keys()):
                        print(run, face)
                        ax[fig_count].plot(self.run_output[run][face][0], self.run_output[run][face][1])
                        ax[fig_count].plot(self.model_output[run][face][0], self.model_output[run][face][1])
                        ax[fig_count].set_title(run)
                    fig_count=fig_count+1
                
        ani=animation.FuncAnimation(fig, animate, interval=1000)
        plt.show()
        
m=ComparitivePlot('model_results.json')
m.animate_plots()
            