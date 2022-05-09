from m2 import OptimizationData
from ttkthemes import ThemedTk
import os
from tkinter import *
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg 
from matplotlib.figure import Figure 
import numpy as np

rootdir = '../CleanData'

def find_neighbours(value, df, colname):
    exactmatch = df[df[colname] == value]
    if not exactmatch.empty:
        return exactmatch.index
    else:
        lowerneighbour_ind = df[df[colname] < value][colname].idxmax()
        upperneighbour_ind = df[df[colname] > value][colname].idxmin()
        return [lowerneighbour_ind, upperneighbour_ind]

class DataVisWindow:
    def update_plot(self):
        #self.ax1.clear()
        #print(self.face_choice_var.get())
        #self.ax1.plot(self.data.data.index, self.data.data[self.face_choice_var.get()])
        try:
            min_index = str(find_neighbours(self.minVal.get(), self.data.data,'time2')[0])
        except:
            min_index = str(0)
            
        try:
            max_index = str(find_neighbours(self.maxVal.get(), self.data.data,'time2')[0])
        except:
            max_index = str(0)
                        
        self.min_index_val.set(min_index)
        self.max_index_val.set(max_index)
        print(self.plotter_check_var.get())
        if self.plotter_check_var.get():
            save_path = f'{self.plotter_save_filepath.get()}/{self.data.name}_{self.face_choice_var.get()}.png'
        else:
            save_path = None

        self.plot(self.ax1, savefile = save_path )
        
    def plot(self, ax, savefile = None):
        ax.clear()
        ax.plot(self.data.data['time'], self.data.data[self.face_choice_var.get()])
        ax.axvline(x = self.minVal.get(), color = 'g')
        ax.axvline(x = self.maxVal.get(), color = 'r')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Temp (K)')
        try:
            bounds_string = f'{self.data.name} {self.face_choice_var.get()} || Min: {self.min_index_val.get()}, Max: {self.max_index_val.get()}'
        except:
            bounds_string = f'{self.data.name} {self.face_choice_var.get()}'
        ax.set_title(bounds_string)
        #ax.text(8, 13, bounds_string)
        pct = 0.05
        ax.set_xlim(self.minVal.get()-pct*self.minVal.get(), self.maxVal.get()+pct*self.maxVal.get())
        
        self.graph.draw()
        print(savefile)
        if savefile != None:
            save_path = f'{savefile}'
            self.fig.savefig(save_path)
        
        #ax.show()
        
    def __init__(self, data_path, root = None, theme = None, name = None):
        
        self.data = OptimizationData(data_path)
        self.data.correct_data()
        #max_val = 
        
        if root is None and theme is None:
            self.root = Tk()
        elif root is None and theme is not None:
            self.root = ThemedTk(theme)
        else:
            self.root = Toplevel(root)
        self.root.title(name) if name is not None else 1
        self.root.geometry('600x300')
        
        #==========================================
        #get which steps we have
        #==========================================
        steps = []
        for key in self.data.data.keys():
            if key[-9:] == 'corrected':
                steps.append(key)
        #==========================================
        
        self.face_selection_frm = Frame(self.root)
        self.face_selection_frm.pack(side = TOP)
        self.face_label = Label(master = self.face_selection_frm, text = 'Select Face: ')
        self.face_label.pack(side = TOP)
        self.face_choice_var = StringVar(self.face_selection_frm)
        self.face_choice_var.set(steps[0])

        self.face_list = OptionMenu(self.face_selection_frm, self.face_choice_var, *steps)
        self.face_list.pack(side = TOP)
        
        self.entry_frm = Frame(self.root)
        
        
        self.fig = Figure()
        self.ax1 = self.fig.add_subplot(1,1,1)
        
        
        self.graph = FigureCanvasTkAgg(self.fig, master = self.root)
        self.graph.get_tk_widget().pack(side = 'top', fill = 'both', expand = True)
        self.minVal = DoubleVar()
        self.minVal.set(1e-8)
        self.maxVal = DoubleVar()
        self.maxVal.set(5*10**-8)
        self.plot(self.ax1)
        self.graph.draw()
        #self.ax1.plot(self.data.data.index, self.data.data[self.face_choice_var.get()])
        left_frame = Frame(self.root)
        left_frame.pack(side = LEFT)
        param_frm = Frame(left_frame)
        param_frm.pack(side = TOP)
        minmax_frm = Frame(param_frm)
        minmax_frm.pack(side = LEFT)
        min_frm = Frame(minmax_frm)
        min_frm.pack(side = TOP)
        max_frm = Frame(minmax_frm)
        max_frm.pack(side = TOP)
        minLabel = Label(master = min_frm, text = "Min: ").pack(side = LEFT)

        minEntry = Entry(master = min_frm, textvariable = self.minVal).pack(side = LEFT)
        maxLabel = Label(master = max_frm, text = "Max: ").pack(side = LEFT)
        maxEntry = Entry(master = max_frm, textvariable = self.maxVal).pack(side = LEFT)
        print(self.data.data.keys())
        try:
            self.min_index = str(find_neighbours(self.minVal.get(), self.data.data,'time')[0])
        except:
            self.min_index = str(0)
        self.min_index_val = StringVar()
        self.min_index_val.set(str(self.min_index[0]))
        self.minIndexLabel = Label(master = min_frm, textvariable = self.min_index_val).pack(side = LEFT)
        
        try:
            self.max_index = str(find_neighbours(self.maxVal.get(), self.data.data,'time')[1])
        except:
            self.max_index = str(len(self.data.data.index)-1)
        self.max_index_val = StringVar()
        self.max_index_val.set(str(self.max_index))
        self.maxIndexLabel = Label(master = max_frm, textvariable = self.max_index_val).pack(side = LEFT)
        
        update_Button = Button(master = self.root, text = 'Update Face', command = self.update_plot).pack(side = BOTTOM)
        
        save_frm = Frame(self.root)
        save_frm.pack(side = TOP)
        self.plotter_check_var = BooleanVar()
        self.plotter_check_var.set(False)
        self.plotter_save_filepath = StringVar()
        self.plotter_save_filepath.set('../../winhome/Desktop/Bounds_Saver')
        plotter_check = Checkbutton(self.root, text = 'Save Figure', var = self.plotter_check_var).pack(side = RIGHT)
        plotter_entry = Entry(self.root, textvariable=self.plotter_save_filepath).pack(side = RIGHT)
        
        self.root.mainloop()

class DataViewer:
    def view_shot_data(self):
        data_path = f'{rootdir}/{self.choice_var.get()}.xlsx'
        self.vis_window = DataVisWindow(data_path = data_path, root = self.root, theme = self.gui_theme, name = self.choice_var.get())
    def __init__(self):
        self.choices = []
        for file in os.listdir(rootdir):
            fp = f'{rootdir}/{file}'
            try:
                a = OptimizationData(fp)
                print(file, 'Format Works')
                self.choices.append(a.name)
            except:
                print(file, 'Data Format Wrong!')
        
        self.gui_theme = 'black'

        self.root = ThemedTk(theme = self.gui_theme)
        self.root.title('Data Viewer')
        self.root.geometry('300x100')

        self.choices_frame = Frame(self.root)
        self.choices_frame.pack(side = TOP)
        self.choice_label = Label(master = self.choices_frame, text = "Select Shot:")
        self.choice_label.pack(side = TOP)
        self.choice_var = StringVar(self.choices_frame)
        self.choice_var.set(self.choices[0])

        self.shot_list = OptionMenu(self.choices_frame, self.choice_var, *self.choices)
        self.shot_list.pack(side = TOP)

        View_Shot_Data_Button = Button(master = self.root, text = "visualize Data", bg='orange', fg='black', command = self.view_shot_data)

        View_Shot_Data_Button.pack()
        
        self.root.mainloop()
        
DataViewer()