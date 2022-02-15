import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from tkinter import *
from tkinter import ttk
import json
from tkinter import messagebox
import opt_func
import numpy as np
from threading import *
import matplotlib.animation as animation

#JSOn func
def load_json(fname):
    with open(fname) as json_file:
        data=json.load(json_file)
        
    return data

def donothing():
   filewin = Toplevel(root)
   button = Button(filewin, text="Do nothing button")
   button.pack()
   
def helloCallBack():
    messagebox.showinfo( "INFO", num_steps)
   
def openNewOptimization():   
    
    global_vars=load_json('global_variables.json')
    newWindow=Tk()
    newWindow.title("New optimization")
    newWindow.geometry("300x600")
    
    frm=ttk.Frame(root, padding=10)
    frm.grid()
    
    Label(newWindow, text="Create New Optimization").grid(column=0, row=0)
    Label(newWindow, text="================================").grid(column=0, row=1)
    Label(newWindow, text="Select Shots & Faces to optimize:").grid(column=0, row=2)
    
    #Display checklist for the shots
    default_shots=['s88773','s88776', 's88780', 's86483']
    
    count=0
    for shot in default_shots:
        globals()[shot]=BooleanVar()
        globals()[shot].set(False)
        str=shot+'chk'
        globals()[str]=Checkbutton(newWindow, text=shot, var=globals()[shot])
        globals()[str].grid(column=0, row=count+3)
        count=count+1
        
    current_row=count+3
    folder_path=global_vars['optimization_data_path']
    run_name=StringVar()
    run_name.set("DEFAULT")
    
    Label(newWindow, text="=================================").grid(column=0, row=current_row+3)
    Label(newWindow, text="Enter Run Name").grid(column=0, row=current_row+4)
    rn=Entry(newWindow, textvariable=run_name).grid(column=0, row=current_row+5)
    
    def print_var():
        content=rn.get()
        messagebox.showinfo('VALUES', 'Run Name: '+content) 
    
    current_row=current_row+6
    
    num_steps=StringVar()
    popsize=StringVar()
    Label(newWindow, text="Number of Finite Element steps").grid(column=0, row=current_row+2)
    Label(newWindow, text="Population size").grid(column=0, row=current_row+4)
    Entry(newWindow, textvariable=num_steps).grid(column=0, row=current_row+3)
    Entry(newWindow, textvariable=popsize).grid(column=0, row=current_row+5)
    
    
    current_row=current_row+6
    a_bounds=StringVar()
    b_bounds=StringVar()
    start_time_bounds=StringVar()
    peak_temp_bounds=StringVar()
    Label(newWindow, text="=================================").grid(column=0, row=current_row+3)
    Label(newWindow, text="BOUNDS").grid(column=0, row=current_row+4)
    Label(newWindow, text="enter as (min, max) pair in parentheses").grid(column=0, row=current_row+5)
    
    #a
    Label(newWindow, text="a").grid(column=0, row=current_row+7)
    Entry(newWindow, textvariable=a_bounds).grid(column=0, row=current_row+8)
    
    #b
    Label(newWindow, text="b").grid(column=0, row=current_row+9)
    Entry(newWindow, textvariable=b_bounds).grid(column=0, row=current_row+10) 
    
    #start_time
    Label(newWindow, text="Start Time").grid(column=0, row=current_row+11)
    Entry(newWindow, textvariable=start_time_bounds).grid(column=0, row=current_row+12)
    
    #peak_temp
    Label(newWindow, text="Peak Temperature").grid(column=0, row=current_row+13)
    Entry(newWindow, textvariable=peak_temp_bounds).grid(column=0, row=current_row+14) 
    
    Label(newWindow, text="=================================").grid(column=0, row=current_row+15) 
    
    V=Button(newWindow, text="Run Optimization", command=print_var).grid(column=0, row=current_row+16)
    
    newWindow.mainloop()    
    
        
def openPreviousOptimization():
    newWindow=Toplevel(root)
    newWindow.title("Optimization")
    newWindow.geometry("200x200")
    Label(newWindow, text="Opening an existing optimization").pack()
   
root = Tk()

#======================================================
#Menu
menubar = Menu(root)
filemenu = Menu(menubar, tearoff=0)
filemenu.add_command(label="New", command=openNewOptimization)
filemenu.add_command(label="Open", command=openPreviousOptimization)
filemenu.add_command(label="Save", command=donothing)
filemenu.add_command(label="Save as...", command=donothing)
filemenu.add_command(label="Close", command=donothing)

filemenu.add_separator()

filemenu.add_command(label="Exit", command=root.quit)
menubar.add_cascade(label="File", menu=filemenu)
editmenu = Menu(menubar, tearoff=0)
editmenu.add_command(label="Undo", command=donothing)

editmenu.add_separator()

editmenu.add_command(label="Cut", command=donothing)
editmenu.add_command(label="Copy", command=donothing)
editmenu.add_command(label="Paste", command=donothing)
editmenu.add_command(label="Delete", command=donothing)
editmenu.add_command(label="Select All", command=donothing)

menubar.add_cascade(label="Edit", menu=editmenu)
helpmenu = Menu(menubar, tearoff=0)
helpmenu.add_command(label="Help Index", command=donothing)
helpmenu.add_command(label="About...", command=donothing)
menubar.add_cascade(label="Help", menu=helpmenu)
#=====================================================

#=====================================================
#Body
srs=StringVar()

global_vars=load_json('global_variables.json')

frm=ttk.Frame(root, padding=10)
frm.grid()

Label(root, text="Create New Optimization").grid(column=0, row=0)
Label(root, text="================================").grid(column=0, row=1)
Label(root, text="Select Shots & Faces to optimize:").grid(column=0, row=2)

#Display checklist for the shots
default_shots=['s88773','s88776', 's88780', 's86483']

count=0
for shot in default_shots:
    globals()[shot]=BooleanVar()
    globals()[shot].set(True)
    str=shot+'chk'
    globals()[str]=Checkbutton(root, text=shot, var=globals()[shot])
    globals()[str].grid(column=0, row=count+3)
    count=count+1
    
current_row=count+3
folder_path=StringVar()
folder_path.set(global_vars['optimization_data_path'])

run_name=StringVar()
run_name.set("RUN NAME")

Label(root, text="=================================").grid(column=0, row=current_row+3)
Label(root, text="Root Folder").grid(column=0, row=current_row+4)
Entry(root, textvariable=folder_path).grid(column=0, row=current_row+5)

current_row=current_row+5

Label(root, text="Enter Run Name").grid(column=0, row=current_row+3)
Entry(root, textvariable=run_name).grid(column=0, row=current_row+4)

current_row=current_row+5

num_steps=StringVar()
num_steps.set(60)
popsize=StringVar()
popsize.set(5)
Label(root, text="Number of Finite Element steps").grid(column=0, row=current_row+2)
Label(root, text="Population size").grid(column=0, row=current_row+4)
Entry(root, textvariable=num_steps).grid(column=0, row=current_row+3)
Entry(root, textvariable=popsize).grid(column=0, row=current_row+5)


current_row=current_row+6
a_bounds=StringVar()
a_bounds.set("( 0.01 , 0.05 )")
b_bounds=StringVar()
b_bounds.set("( 10 , 100 )")
start_time_bounds=StringVar()
start_time_bounds.set("( 1.75E-8 , 2.25E-8 )")
peak_temp_bounds=StringVar()
peak_temp_bounds.set("( 25000 , 60000 )")
Label(root, text="=================================").grid(column=0, row=current_row+3)
Label(root, text="BOUNDS").grid(column=0, row=current_row+4)
Label(root, text="enter as (min, max) pair in parentheses").grid(column=0, row=current_row+5)

#a
Label(root, text="a").grid(column=0, row=current_row+7)
Entry(root, textvariable=a_bounds).grid(column=0, row=current_row+8)

#b
Label(root, text="b").grid(column=0, row=current_row+9)
Entry(root, textvariable=b_bounds).grid(column=0, row=current_row+10) 

#start_time
Label(root, text="Start Time").grid(column=0, row=current_row+11)
Entry(root, textvariable=start_time_bounds).grid(column=0, row=current_row+12)

#peak_temp
Label(root, text="Peak Temperature").grid(column=0, row=current_row+13)
Entry(root, textvariable=peak_temp_bounds).grid(column=0, row=current_row+14) 

Label(root, text="=================================").grid(column=0, row=current_row+15) 

def print_var():
    #Get Bounds
    s_a_bounds=a_bounds.get()
    s_b_bounds=b_bounds.get()
    s_start_time_bounds=start_time_bounds.get()
    s_peak_temp_bounds=peak_temp_bounds.get()
    a_min=float(s_a_bounds.split('(')[1].split(')')[0].split(',')[0].replace(" ",""))
    a_max=float(s_a_bounds.split('(')[1].split(')')[0].split(',')[1].replace(" ",""))
    b_min=float(s_b_bounds.split('(')[1].split(')')[0].split(',')[0].replace(" ",""))
    b_max=float(s_b_bounds.split('(')[1].split(')')[0].split(',')[1].replace(" ",""))
    start_time_min=float(s_start_time_bounds.split('(')[1].split(')')[0].split(',')[0].replace(" ",""))
    start_time_max=float(s_start_time_bounds.split('(')[1].split(')')[0].split(',')[1].replace(" ",""))
    peak_temp_min=float(s_peak_temp_bounds.split('(')[1].split(')')[0].split(',')[0].replace(" ",""))
    peak_temp_max=float(s_peak_temp_bounds.split('(')[1].split(')')[0].split(',')[1].replace(" ",""))
    my_bounds=((a_min, a_max),(b_min, b_max),(start_time_min, start_time_max),(peak_temp_min, peak_temp_max))
    
    #Get Shots
    default_shots=['s88773','s88776', 's88780', 's86483']
    shots_run=[]
    for shot in default_shots:
        if globals()[shot].get()==True:
            shots_run.append(shot)
            
    for k in range(len(shots_run)):
        if k==0:
            shots_run_string=shots_run[k]
        else:
            shots_run_string=shots_run_string+','+shots_run[k]
            
    srs.set(shots_run_string)
    
    content=StringVar()
    content.set(my_bounds)
    lines=['PARAMETERS', '=================', ' ','--','RUN NAME: '+run_name.get(), ' ','--','Shots: '+srs.get(),' ','--', 'Parameters:', ' ', 'Bounds: '+content.get(),'Number of FEM steps: '+num_steps.get(), 'Population size: '+popsize.get() ]
    #messagebox.showinfo('VALUES', 'Bounds: '+content.get()) 
    #infobox=messagebox.showinfo('VALUES', "\n".join(lines))
    res=messagebox.askquestion('Run Optimization?', "\n".join(lines))
    run_opt.set(False)
    if res=='yes':
        run_opt.set(True)
        
    elif res=='no':
        messagebox.showinfo('Return','You will now return to the application screen')
    else:
        messagebox.showwarning('error', 'Something went wrong!')
        
    count=0
    while count<10000:
        count=count+1
       

run_opt=BooleanVar()
run_opt.set(False)
 
#Check_Button=Button(root, text="Check Optimization Progress", command=)

def threading():
    t1=Thread(target=run_optimization)
    t1.start()

def run_optimization():
    if run_opt.get()==True:
        #Rewrite the global json to include the new root folder path
        global_json_as_dict=opt_func.load_json('global_variables.json')
        global_json_as_dict['optimization_data_path']=folder_path.get()

        global_json=json.dumps(global_json_as_dict)
        f=open('global_variables.json', 'w')
        f.write(global_json)
        f.close()

        default_shots=['s88773','s88776', 's88780', 's86483']
        shots_run=[]
        for shot in default_shots:
            if globals()[shot].get()==True:
                shots_run.append(shot)

        for k in range(len(shots_run)):
            if k==0:
                shots_run_string=shots_run[k]
            else:
                shots_run_string=shots_run_string+','+shots_run[k]

        s_a_bounds=a_bounds.get()
        s_b_bounds=b_bounds.get()
        s_start_time_bounds=start_time_bounds.get()
        s_peak_temp_bounds=peak_temp_bounds.get()
        a_min=float(s_a_bounds.split('(')[1].split(')')[0].split(',')[0].replace(" ",""))
        a_max=float(s_a_bounds.split('(')[1].split(')')[0].split(',')[1].replace(" ",""))
        b_min=float(s_b_bounds.split('(')[1].split(')')[0].split(',')[0].replace(" ",""))
        b_max=float(s_b_bounds.split('(')[1].split(')')[0].split(',')[1].replace(" ",""))
        start_time_min=float(s_start_time_bounds.split('(')[1].split(')')[0].split(',')[0].replace(" ",""))
        start_time_max=float(s_start_time_bounds.split('(')[1].split(')')[0].split(',')[1].replace(" ",""))
        peak_temp_min=float(s_peak_temp_bounds.split('(')[1].split(')')[0].split(',')[0].replace(" ",""))
        peak_temp_max=float(s_peak_temp_bounds.split('(')[1].split(')')[0].split(',')[1].replace(" ",""))
        my_bounds=((a_min, a_max),(b_min, b_max),(start_time_min, start_time_max),(peak_temp_min, peak_temp_max))

        #Run optimization
        Bounds=my_bounds
        NumSteps=int(num_steps.get())
        MyShots=shots_run
        OptimizationName=run_name.get()
        PopSize=int(popsize.get())
        opt_func.run_optimization(NumSteps, MyShots, OptimizationName, Bounds, PopSize)
        
    elif run_opt.get()==False:
        messagebox.showinfo('Not running', 'optimization not running. Please set up new parameters & run')
        
    else:
        messagebox.showwarning()
        
def combined():
    print_var()
    
    threading()
        
def visualize_plots():
    global_vars=load_json('global_variables.json')
    fp=global_vars['optimization_data_path']+'/'+global_vars['optimization_name']+'/'+'optimization_output.csv'
    newWindow=Tk()
    newWindow.title(global_vars['optimization_name'])
    newWindow.geometry("800x800")
    
    frm=ttk.Frame(root, padding=10)
    frm.grid()
    
    parameters=('a', 'b', 'peak_temp', 'start_time')

    fig=plt.figure()
    ax1=fig.add_subplot(2,3,1) #sls
    ax2=fig.add_subplot(2,3,2) #peak_temp
    ax3=fig.add_subplot(2,3,3) #a
    ax4=fig.add_subplot(2,3,4) #b
    ax5=fig.add_subplot(2,3,5) #start_time
    
    def animate(i):
        
        #========================================
        #Parameters
        #========================================
        pullData = open(fp,"r").read()
        dataArray = pullData.split('\n')

        #defining arrays for each relevant value...
        itar = [] #iterations
        slsar = [] #sum of least squares
        ptar = [] #peak temp
        aar = [] #a
        bar = [] #b
        tsar = [] #time shift

        for eachLine in dataArray:
            if len(eachLine)>1:
                it, sls, a, b, pt, ts = eachLine.split(',')
                itar.append(float(it))
                slsar.append(float(sls))
                ptar.append(float(pt))
                aar.append(float(a))
                bar.append(float(b))
                tsar.append(float(ts))

        #==========================================
        #plotting
        #==========================================
        ax1.clear()
        ax2.clear()
        ax3.clear()
        ax4.clear()
        ax5.clear()

        ax1.set_title('sls')
        ax2.set_title('peak_temp')
        ax3.set_title('a')
        ax4.set_title('b')
        ax5.set_title('start_time')

        ax1.plot(itar,slsar)
        ax2.plot(itar,ptar)
        ax3.plot(itar,aar)
        ax4.plot(itar,bar)
        ax5.plot(itar,tsar)
        
    ani=animation.FuncAnimation(fig,animate, interval=20000)
    
    newWindow.mainloop()
        
Run_Button=Button(root, text="Run Optimization", command=combined).grid(column=0, row=current_row+16) 
Visualization_Button=Button(root, text="Visualize Optimization", command=visualize_plots).grid(column=0, row=current_row+17)

root.config(menu=menubar)
root.mainloop()