from tkinter import *
from tkinter import ttk
import json
from tkinter import messagebox

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
folder_path=global_vars['optimization_data_path']
run_name=StringVar()
run_name.set("RUN NAME")

Label(root, text="=================================").grid(column=0, row=current_row+3)
Label(root, text="Enter Run Name").grid(column=0, row=current_row+4)
Entry(root, textvariable=run_name).grid(column=0, row=current_row+5)

def print_var():
    content=run_name.get()
    messagebox.showinfo('VALUES', 'Run Name: '+content) 

current_row=current_row+6

num_steps=IntVar()
num_steps.set(60)
popsize=IntVar()
popsize.set(5)
Label(root, text="Number of Finite Element steps").grid(column=0, row=current_row+2)
Label(root, text="Population size").grid(column=0, row=current_row+4)
Entry(root, textvariable=num_steps).grid(column=0, row=current_row+3)
Entry(root, textvariable=popsize).grid(column=0, row=current_row+5)


current_row=current_row+6
a_bounds=StringVar()
a_bounds.set("( , )")
b_bounds=StringVar()
b_bounds.set("( , )")
start_time_bounds=StringVar()
start_time_bounds.set("( , )")
peak_temp_bounds=StringVar()
peak_temp_bounds.set("( , )")
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

V=Button(root, text="Run Optimization", command=print_var).grid(column=0, row=current_row+16)  

root.config(menu=menubar)
root.mainloop()