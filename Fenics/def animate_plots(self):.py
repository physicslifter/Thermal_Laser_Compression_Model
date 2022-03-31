def animate_plots(self):
        self.run_IDs=self.run_dictionary.keys()
        num_runs=len(self.run_IDs)
        ncols=int(math.ceil(num_runs**0.5))
        nrows=int(math.ceil(num_runs/ncols))
        fig, ax = plt.subplots(nrows, ncols)
        
        def animate(i,nrows=nrows, num_runs=num_runs):
            if nrows>1:
                row_count=0
                col_count=0
                fig_count=0
                for run in self.run_IDs:
                    ax[row_count, col_count].clear() #clear data from previous plot
                    for face in list(self.model_data[run].keys()):
                        #print(run, face)
                        ax[row_count,col_count].scatter(self.real_data[run][face][0], self.real_data[run][face][1])
                        ax[row_count,col_count].plot(self.model_data[run][face][0], self.model_data[run][face][1])
                        ax[row_count,col_count].set_title(run)
                    fig_count=fig_count+1
                    col_count=fig_count%nrows
                    row_count=int(fig_count/nrows)
            else:
                if num_runs==1:
                    fig_count=0
                    for run in self.run_IDs:
                        ax.clear() #clear data from previous plot
                        for face in list(self.model_data[run].keys()):
                            print(run, face)
                            ax.scatter(self.real_data[run][face][0], self.real_data[run][face][1])
                            ax.plot(self.model_data[run][face][0], self.model_data[run][face][1])
                            ax.set_title(run)
                        fig_count=fig_count+1
                else:
                    fig_count=0
                    for run in self.run_IDs:
                        ax[fig_count].clear() #clear data from previous plot
                        for face in list(self.model_data[run].keys()):
                            print(run, face)
                            ax[fig_count].scatter(self.real_data[run][face][0], self.real_data[run][face][1])
                            ax[fig_count].plot(self.model_data[run][face][0], self.model_data[run][face][1])
                            ax[fig_count].set_title(run)
                        fig_count=fig_count+1