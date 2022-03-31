from tkinter import *
from tkinter import ttk
import json
import checklistcombobox as c
from tkinter import messagebox
import main

def load_json(fname):
    with open(fname) as json_file:
        data=json.load(json_file)
        
    return data

all_shots=['s88773', 's88776', 's88780', 's86483']

def donothing():
    filewin = Toplevel(root)
    button = Button(filewin, text="Do nothing button")
    button.pack()
    
def run_opt(run_dict, bounds):
    pass
    
def show_run_dict():
    s_a_bounds=a_bounds.get()
    s_b_bounds=b_bounds.get()
    s_start_time_bounds=st_bounds.get()
    s_peak_temp_bounds=pt_bounds.get()
    a_min=float(s_a_bounds.split('(')[1].split(')')[0].split(',')[0].replace(" ",""))
    a_max=float(s_a_bounds.split('(')[1].split(')')[0].split(',')[1].replace(" ",""))
    b_min=float(s_b_bounds.split('(')[1].split(')')[0].split(',')[0].replace(" ",""))
    b_max=float(s_b_bounds.split('(')[1].split(')')[0].split(',')[1].replace(" ",""))
    start_time_min=float(s_start_time_bounds.split('(')[1].split(')')[0].split(',')[0].replace(" ",""))
    start_time_max=float(s_start_time_bounds.split('(')[1].split(')')[0].split(',')[1].replace(" ",""))
    peak_temp_min=float(s_peak_temp_bounds.split('(')[1].split(')')[0].split(',')[0].replace(" ",""))
    peak_temp_max=float(s_peak_temp_bounds.split('(')[1].split(')')[0].split(',')[1].replace(" ",""))
    my_bounds=((a_min, a_max),(b_min, b_max),(start_time_min, start_time_max),(peak_temp_min, peak_temp_max))
    
    lines=[]
    shots_run={}
    for shot in all_shots:
        print(shot, globals()[shot].get())
        for i in range(3):
            face_string=shot+'_face'+str(i+1)
            face_num=i+1
            print('shot: '+shot+', face: '+face_string+', val: ', globals()[face_string].get())
        if globals()[shot].get()==True:
            lines.append(shot)
            shot_faces_list=[]
            for i in range(3):
                face_string=shot+'_face'+str(i+1)
                face_num=i+1
                if globals()[face_string].get()==True:
                    shot_faces_list.append(face_num)
                    lines.append(face_string)
            shots_run[shot]=shot_faces_list
            
    res=messagebox.askquestion('Run Optimization?', shots_run)
    
    if res=='yes':
        run_opt(shots_run, my_bounds)

        
def visualize_optimization():
    pass
        
        
                
   
   
root=Tk()
root.geometry('250x500')
root.title('TLCM Optimizer')

#==========================================
#Load Global variables
global_vars=load_json('global_variables.json')
#==========================================

#===========================================
#Header
header_frame=ttk.Frame(root, padding=10)
header_frame.pack(side=TOP)
Label(header_frame, text="Create New Optimization").pack()
Label(header_frame, text="========================").pack()
#===========================================

#===========================================
#Checkboxes for shots and frames
#===========================================
shot_and_face_frame=ttk.Frame(root)
shot_and_face_frame.pack(side=TOP)
Label(shot_and_face_frame, text="Select Shots & faces to Optimize:").pack()
Label(shot_and_face_frame, text=' ').pack()

shot_frames={}
for shot in all_shots:
    shot_frames[shot]=ttk.Frame(shot_and_face_frame)
    shot_frames[shot].pack(side=TOP)
    globals()[shot]=BooleanVar()
    globals()[shot].set(False)
    str1=shot+'chk'
    globals()[str1]=Checkbutton(shot_frames[shot], text=shot, var=globals()[shot], font='helvetica 11 bold')
    globals()[str1].pack(side=LEFT)
    for i in range(3):
        num=str(i+1)
        var_string=shot+'_'+'face'+num+'_chk'
        face_text=shot+'_face'+str(i+1)
        disp_face_text='face'+str(i+1)
        globals()[face_text]=BooleanVar()
        globals()[face_text].set(False)
        globals()[var_string]=Checkbutton(shot_frames[shot], text=disp_face_text, var=globals()[face_text])
        globals()[var_string].pack(side=LEFT)
#================================================

#================================================
#Optimization Parameters Section
opt_params_frm=ttk.Frame(root, padding=10)
opt_params_frm.pack(side=TOP)

Label(opt_params_frm, text="===========================").pack()
Label(opt_params_frm, text="Optimization Parameters").pack()
Label(opt_params_frm, text=" ").pack()

#Folder Path
fp_frm=ttk.Frame(opt_params_frm)
fp_frm.pack(side=TOP)
folder_path=StringVar()
folder_path.set(global_vars['optimization_data_path'])
Label(fp_frm, text="Root Folder").pack(side=LEFT)
Entry(fp_frm, textvariable=folder_path).pack(side=LEFT)

#Optimization run name
rn_frm=ttk.Frame(opt_params_frm)
rn_frm.pack(side=TOP)
run_name=StringVar()
run_name.set("Default")
Label(rn_frm, text="Run Name ").pack(side=LEFT)
Entry(rn_frm, textvariable=run_name).pack(side=LEFT)

#Number of steps
ns_frm=ttk.Frame(opt_params_frm)
ns_frm.pack(side=TOP)
num_steps=IntVar()
num_steps.set(60)
Label(ns_frm, text="Num Steps ").pack(side=LEFT)
Entry(ns_frm, textvariable=num_steps).pack(side=LEFT)

#Population size
ps_frm=ttk.Frame(opt_params_frm)
ps_frm.pack(side=TOP)
popsize=IntVar()
popsize.set(5)
Label(ps_frm, text="Pop Size     ").pack(side=LEFT)
Entry(ps_frm, textvariable=popsize).pack(side=LEFT)
#======================================

#======================================
#Bounds
#======================================
bounds_frm=ttk.Frame(root, padding=10)

Label(opt_params_frm, text="===========================").pack()
Label(opt_params_frm, text="Bounds").pack()

#a
a_frm=ttk.Frame(root)
a_frm.pack(side=TOP)
a_bounds=StringVar()
a_bounds.set('( 0.01 , 0.05 )')
Label(a_frm, text="a:                ").pack(side=LEFT)
Entry(a_frm, textvariable=a_bounds).pack(side=LEFT)

#b
b_frm=ttk.Frame(root)
b_frm.pack(side=TOP)
b_bounds=StringVar()
b_bounds.set('( 0 , 100 )')
Label(b_frm, text="b:                ").pack(side=LEFT)
Entry(b_frm, textvariable=b_bounds).pack(side=LEFT)

#start_time
st_frm=ttk.Frame(root)
st_frm.pack(side=TOP)
st_bounds=StringVar()
st_bounds.set('( 1.75E-8 , 2.25E-8 )')
Label(st_frm, text="Start time: ").pack(side=LEFT)
Entry(st_frm, textvariable=st_bounds).pack(side=LEFT)

#peak_temp
pt_frm=ttk.Frame(root)
pt_frm.pack(side=TOP)
pt_bounds=StringVar()
pt_bounds.set('( 30000 , 100000 )')
Label(pt_frm, text="Peak temp: ").pack(side=LEFT)
Entry(pt_frm, textvariable=pt_bounds).pack(side=LEFT)
#===============================================

#===============================================
#Buttons
#============================================
# ===
Label(root, text="===========================").pack(side=TOP)
Vis_Button=Button(root, text="Visualize Optimization", command=visualize_optimization).pack(side=BOTTOM)
Run_Button=Button(root, text="Run Optimization", command=show_run_dict).pack(side=BOTTOM)

#===============================================
root.mainloop()