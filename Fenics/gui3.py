from tkinter import *
from tkinter import ttk
import json
import checklistcombobox as c
from tkinter import messagebox
import main
from pprint import pprint
from multiprocessing import Process
import pprint
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg 
from matplotlib.figure import Figure 
import numpy as np
import matplotlib.animation as animation
from tkinter.scrolledtext import ScrolledText

def load_json(fname):
    with open(fname) as json_file:
        data=json.load(json_file)
        
    return data

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    mov_avg = ret[n - 1:] / n
    
    return mov_avg


all_shots=['s88773', 's88776', 's88780', 's86483', 's86480', 's86484', 's88774']

class DataChoppingWindow:
    def __init__(self, root):
        pass
    
class NewOptimizationWindow:
    def donothing(self):
        filewin = Toplevel(self.root)
        button = Button(filewin, text="Do nothing button")
        button.pack()
        
    def get_start_stop_dict(self):
        new_ss_dict={}
        disp_dict={}
        for shot in self.default_start_stop_dict.keys():
            #Check if the shot is in the starts and stops listed
            ss_has_shot=False
            for ss_shot in self.starts_and_stops.keys():
                if ss_shot == shot:
                    ss_has_shot=True
                    
            if ss_has_shot==True:
                new_ss_dict[shot]={}
                disp_dict[shot]={}
                for face in self.default_start_stop_dict[shot].keys():
                    ss_has_face=False
                    for ss_face in self.starts_and_stops[shot].keys():
                        if face==ss_face:
                            ss_has_face=True
                    if ss_has_face==True:
                        disp_dict[shot][face]=[self.starts_and_stops[shot][face][0].get(), self.starts_and_stops[shot][face][1].get()]
                        new_ss_dict[shot][face]=[self.starts_and_stops[shot][face][0].get(), self.starts_and_stops[shot][face][1].get()]
                    else:
                        new_ss_dict[shot][face]=self.default_start_stop_dict[shot][face]
            else:
                new_ss_dict[shot]=self.default_start_stop_dict[shot]
                
        lines=['Start & stop points for each face', '================================']
        for shot in disp_dict.keys():
            lines.append(shot)
            for face in disp_dict[shot].keys():
                fancy_name='face'+str(face)
                fancy_string=fancy_name+': '+str(new_ss_dict[shot][face])
                lines.append(fancy_string)
                
                
        res=messagebox.askquestion('Chop Data?',"\n".join(lines))
        
        if res=='yes':
            print(new_ss_dict)
            self.real_data.get_data(new_ss_dict)
        else:
            messagebox.showinfo('Return', 'Please Chop data to your preference')
        
    def chop_data(self):
        datawin=Toplevel(self.root)
        datawin.title(self.run_name.get()+'Data Chopping')
        #========================================
        #Entry Windows for the bounds for each run
        #========================================
        
        #========================================
        #Default Data Chopping
        '''self.default_start_stop_dict={
            's88773':{
                    1:[179, 240],
                    2:[259, 390],       
                    3:[310, 700]
            },
            's88776':{
                    1:[88, 155],        
                    2:[93, 400],        
                    3:[246, 468]
            },
            's88780':{
                    1:[145, 190],       
                    2:[240, 350],       
                    3:[310, 502]
            },
            's86483':{
                    1:[355, 400],       
                    2:[582, 732],       
                    3:[377, 800]
            }
        }'''
        
        self.default_start_stop_dict={
        's88773':{
                1:[179, 240],
                2:[259, 390],       
                3:[310, 700]
        },
        's88776':{
                1:[88, 155],        
                2:[93, 400],        
                3:[246, 468]
        },
        's88780':{
                1:[145, 190],       
                2:[240, 350],       
                3:[310, 502]
        },
        's86483':{
                1:[355, 400],       
                2:[582, 732],       
                3:[377, 800]
        },
        's86480':{
                1:[202, 234],
                2:[293, 345],
                3:[310, 473]
        },
        's86484':{
                1:[187, 256],
                2:[195, 322],
                3:[670, 772]
        },        
        's88774':{
                1:[199, 284],
                2:[190, 458],
                3:[333, 725]
        },        
}
        #========================================
        
        
        shots_run={}
        for shot in all_shots:
            print(shot, self.shot_var_dict[shot].get())
            for i in range(3):
                face_string=shot+'_face'+str(i+1)
                face_num=i+1
                print('shot: '+shot+', face: '+face_string+', val: ', self.shot_face_var_dict[shot][i].get())
            if self.shot_var_dict[shot].get()==True:
                shot_faces_list=[]
                for i in range(3):
                    face_string=shot+'_face'+str(i)
                    face_num=i+1
                    if self.shot_face_var_dict[shot][i].get()==True:
                        shot_faces_list.append(face_num)
                shots_run[shot]=shot_faces_list

        self.shots_and_faces_run=shots_run
        
        data_dict=self.shots_and_faces_run
        starts_and_stops={}
        for shot in data_dict.keys(): #iterate through the shots
            shot_frame=ttk.Frame(datawin)
            shot_frame.pack(side=TOP)
            Label(shot_frame, text='==========================').pack(side=TOP)
            Label(shot_frame, text='Shot '+shot+': ').pack(side=TOP)
            face_start_stop={}
            for face in data_dict[shot]:
                face_name='face'+str(face)
                face_frame=ttk.Frame(shot_frame) #Create a frame for the shot within the face frame
                face_frame.pack(side=LEFT)
                start=self.default_start_stop_dict[shot][face][0]
                stop=self.default_start_stop_dict[shot][face][1]
                Label(face_frame, text=face_name+' Bounds:').pack(side=TOP)
                start_frame=ttk.Frame(face_frame)
                stop_frame=ttk.Frame(face_frame)
                start_frame.pack(side=LEFT)
                stop_frame.pack(side=LEFT)
                Label(start_frame, text="Start: ").pack(side=LEFT)
                start_var=IntVar()
                start_var.set(int(start))
                Entry(start_frame, textvariable=start_var).pack(side=LEFT)
                Label(stop_frame, text="Stop: ").pack(side=LEFT)
                stop_var=IntVar()
                stop_var.set(int(stop))
                Entry(stop_frame, textvariable=stop_var).pack(side=LEFT)
                
                start_and_stop=[start_var, stop_var]
                face_start_stop[face]=start_and_stop
            starts_and_stops[shot]=face_start_stop
        self.starts_and_stops=starts_and_stops

        updateParametersButton=Button(master=datawin, text="Chop data", command=self.get_start_stop_dict).pack(side=BOTTOM)
       
        datawin.mainloop()
        
        #--
        #a_frm=ttk.Frame(bounds_frm)
        #a_frm.pack(side=TOP)
        #self.a_bounds=StringVar()
        #self.a_bounds.set('( 0.01 , 0.05 )')
        
        #--
        #fp_frm=ttk.Frame(opt_params_frm)
        #fp_frm.pack(side=TOP)
        #self.folder_path=StringVar()
        #self.folder_path.set(self.global_vars['optimization_data_path'])
        #@Label(fp_frm, text="Root Folder").pack(side=LEFT)
        #Entry(fp_frm, textvariable=self.folder_path).pack(side=LEFT)
        
    def save_plots(self):
        fp=self.folder_path.get()     #Get Filepath
        rn=self.run_name.get()
        folder_path='../../winhome/Desktop/optimization_data/20220315/'+fp+'/'+rn
        #Saving the raw data plot
        
    
    def opt_func(self):
        main.optimize(self.group, self.opt_data, bounds=self.my_bounds, popsize=self.popsize.get())
    
    def run_threaded_opt(self):
        path_for_optimization=self.folder_path.get()+'/'+self.run_name.get()
        self.group=main.OptimizationGroup(self.meshpoints.get(), self.shots_and_faces_run,path_for_optimization, self.real_data)
        global t1 
        t1=Process(target=self.opt_func)
        t1.start()
        #while t1.is_alive():
        #    self.opt_running_flag.set(True)
        #self.save_plots()
        
    
    def threaded_run(self):
        pass

    def show_run_dict(self):
        s_a_bounds=self.a_bounds.get()
        s_b_bounds=self.b_bounds.get()
        s_start_time_bounds=self.st_bounds.get()
        s_peak_temp_bounds=self.pt_bounds.get()
        a_min=float(s_a_bounds.split('(')[1].split(')')[0].split(',')[0].replace(" ",""))
        a_max=float(s_a_bounds.split('(')[1].split(')')[0].split(',')[1].replace(" ",""))
        b_min=float(s_b_bounds.split('(')[1].split(')')[0].split(',')[0].replace(" ",""))
        b_max=float(s_b_bounds.split('(')[1].split(')')[0].split(',')[1].replace(" ",""))
        start_time_min=float(s_start_time_bounds.split('(')[1].split(')')[0].split(',')[0].replace(" ",""))
        start_time_max=float(s_start_time_bounds.split('(')[1].split(')')[0].split(',')[1].replace(" ",""))
        peak_temp_min=float(s_peak_temp_bounds.split('(')[1].split(')')[0].split(',')[0].replace(" ",""))
        peak_temp_max=float(s_peak_temp_bounds.split('(')[1].split(')')[0].split(',')[1].replace(" ",""))
        self.my_bounds=((a_min, a_max),(b_min, b_max),(start_time_min, start_time_max),(peak_temp_min, peak_temp_max))

        lines=[]
        shots_run={}
        for shot in all_shots:
            print(shot, self.shot_var_dict[shot].get())
            for i in range(3):
                face_string=shot+'_face'+str(i+1)
                face_num=i+1
                print('shot: '+shot+', face: '+face_string+', val: ', self.shot_face_var_dict[shot][i].get())
            if self.shot_var_dict[shot].get()==True:
                lines.append(shot+':')
                shot_faces_list=[]
                for i in range(3):
                    face_string=shot+'_face'+str(i)
                    face_num=i+1
                    if self.shot_face_var_dict[shot][i].get()==True:
                        shot_faces_list.append(face_num)
                        lines.append(face_string)
                shots_run[shot]=shot_faces_list

        lines=['Parameters', 
               '================', 
               ' ', 
               '--', 
               'Run Name: '+str(self.run_name.get()), 
               ' ', 
               '--', 
               ' ', 
               'Bounds: '+str(self.my_bounds), 
               ' ',
               'No. of FEM meshpoints: '+str(self.meshpoints.get()), 
               ' ',
               'Population size: '+str(self.popsize.get()), 
               ' ',
               'Shots & Faces: ',
               ' ']
        for shot in shots_run.keys(): #For each shot we have that is being run...
            lines.append('  '+shot)
            for face in shots_run[shot]:
                key='       face'+str(face+1)
                lines.append(key)

        self.shots_and_faces_run=shots_run

        res=messagebox.askquestion('Run Optimization?', "\n".join(lines))

        if res=='yes':
            self.run_threaded_opt()
            
        else:
            messagebox.showinfo('Return', 'Please change parameters to be of your liking')

    def vis_current_fit(self):
        window=Toplevel(self.root)
        window.title('Current Fit - LIVE PLOT')
        
    def save_log(self):
        fp=self.folder_path.get()
        rn=self.run_name.get()
        folder_path='../../winhome/Desktop/optimization_data/20220315/'+fp+'/'+rn
        cur_inp=self.text_box.get("1.0", END)
        log_filepath=folder_path+'/log.txt'
        f1=open(log_filepath, "w")
        f1.write(cur_inp)
        f1.close()
        
    def log_comments(self):
        window=Toplevel(self.root)
        window.title('Optimization Log: '+self.run_name.get())
        window.geometry('600x275')
        window.config(bg="#84BF04")
        frame=Frame(window)
        self.text_box=Text(
            frame,
            height=15,
            width=75,
            wrap='word'
        )
        message='Enter Log text here: '
        self.text_box.insert('end', message)
        self.text_box.pack(side=LEFT, expand=True)
        
        sb=Scrollbar(frame)
        sb.pack(side=RIGHT, fill=BOTH)
        
        self.text_box.config(yscrollcommand=sb.set)
        sb.config(command=self.text_box.yview)
        
        frame.pack(expand=True)
        
        log_button=Button(window, text="Save Log", command=self.save_log).pack(side=TOP)
        
        window.mainloop()

    def visualize_optimization(self):
        fp=self.folder_path.get()     #Get Filepath
        rn=self.run_name.get()
        optimization_filepath='../../winhome/Desktop/optimization_data/20220315/'+fp+'/'+rn+'/optimization_output.csv'
        
        newWindow=Toplevel(self.root)
        newWindow.title('Optimization LIVE PLOT')
        newWindow.geometry("650x500")
        lab=Label(newWindow, text="Live Plotting", bg="White").pack()

        fig=Figure()
        ax1=fig.add_subplot(2,3,1) #sls
        ax2=fig.add_subplot(2,3,2) #peak_temp
        ax3=fig.add_subplot(2,3,3) #a
        ax4=fig.add_subplot(2,3,4) #b
        ax5=fig.add_subplot(2,3,5) #start_time
        
        graph=FigureCanvasTkAgg(fig, master=newWindow)
    
        graph.get_tk_widget().pack(side="top",fill='both',expand=True)
        
        def animate(i):
            pullData=open(optimization_filepath, "r").read()
            dataArray=pullData.split('\n')
            
            #Defining list in which to store info
            itar, slsar, aar, bar, star, ptar = [], [], [], [], [], []
            
            for eachLine in dataArray:
                line=eachLine.replace('[','').replace(']','') #Get rid of brackets
                if len(line)>0:
                    it, sls, a, b, st, pt = line.split(',')    
                    itar.append(float(it))
                    slsar.append(float(sls))
                    aar.append(float(a))
                    bar.append(float(b))
                    star.append(float(st))
                    ptar.append(float(pt))
            
            #===========================
            # Get Moving Averages of values
            #===========================
            
            mov_avg=15
            
            if int(it)>mov_avg:
                ma_sls=moving_average(slsar, n=mov_avg)
                ma_a=moving_average(aar, n=mov_avg)
                ma_b=moving_average(bar, n=mov_avg)
                ma_st=moving_average(star, n=mov_avg)
                ma_pt=moving_average(ptar, n=mov_avg)
            
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

                ax1.scatter(itar,slsar, s=10, alpha=0.7, edgecolors='k')
                ax1.plot(itar[(mov_avg-1):],ma_sls, linewidth=0.5, color='r')
                ax2.scatter(itar,ptar, s=10, alpha=0.7, edgecolors='k')
                ax2.plot(itar[(mov_avg-1):], ma_pt, linewidth=0.5, color='r')
                ax3.scatter(itar,aar,  s=10, alpha=0.7, edgecolors='k')
                ax3.plot(itar[(mov_avg-1):], ma_a, linewidth=0.5, color='r')
                ax4.scatter(itar,bar,  s=10, alpha=0.7, edgecolors='k')
                ax4.plot(itar[(mov_avg-1):], ma_b, linewidth=0.5, color='r')
                ax5.scatter(itar,star,  s=10, alpha=0.7, edgecolors='k')
                ax5.plot(itar[(mov_avg-1):], ma_st, linewidth=0.5, color='r')
        
        ani=animation.FuncAnimation(fig,animate, interval=1000)
        
                         
        #sample=np.arange(10)
        #ax1.plot(sample, sample)
        #ax2.plot(sample, sample)
        #ax3.plot()
        
        newWindow.mainloop()
        
    
    def __init__(self, root):
        
        
        #Boolean flag for whether the optimization is running
        self.opt_running_flag=BooleanVar()
        self.opt_running_flag.set(True)
        
        self.global_vars=load_json('global_variables.json')
        
        self.opt_data=main.OptimizationData('')
        
        self.real_data=main.Data()
        self.real_data.get_data()
        
        #Setting up the root
        self.root=Toplevel(root)
        self.root.geometry('250x550')
        self.root.title('TLCM Optimizer')

        #===========================================
        #Header
        self.header_frame=ttk.Frame(self.root, padding=10)
        self.header_frame.pack(side=TOP)
        Label(self.header_frame, text="Create New Optimization").pack()
        Label(self.header_frame, text="========================").pack()
        #===========================================
        
        #===========================================
        #Checkboxes for shots and frames
        #===========================================
        self.shot_and_face_frame=ttk.Frame(self.root)
        self.shot_and_face_frame.pack(side=TOP)
        Label(self.shot_and_face_frame, text="Select Shots & faces to Optimize:").pack()
        Label(self.shot_and_face_frame, text=' ').pack()

        self.shot_frames={}
        self.shot_var_dict={}
        self.shot_chk_dict={}
        self.shot_face_var_dict={}
        self.shot_face_chk_dict={}
        for shot in all_shots:
            self.shot_frames[shot]=ttk.Frame(self.shot_and_face_frame)
            self.shot_frames[shot].pack(side=TOP)
            self.shot_var_dict[shot]=BooleanVar()
            self.shot_var_dict[shot].set(False)
            
            #globals()[shot]=BooleanVar()
            #globals()[shot].set(False)
            
            str1=shot+'chk'
            print(self.shot_frames[shot])
            print(self.shot_var_dict[shot])
            self.shot_chk_dict[str1]=Checkbutton(self.shot_frames[shot], text=shot, var=self.shot_var_dict[shot], font='helvetica 11 bold')
            self.shot_chk_dict[str1].pack(side=LEFT)
            
            #globals()[str1]=Checkbutton(self.shot_frames[shot], text=shot, var=self.check_shot_dict, font='helvetica 11 bold')
            #globals()[str1].pack(side=LEFT)
            
            self.face_var_dict={}
            self.face_chk_dict={}
            for i in range(3):
                num=str(i+1)
                var_string=shot+'_'+'face'+num+'_chk'
                face_text=shot+'_face'+str(i+1)
                disp_face_text='face'+str(i+1)
                
                self.face_var_dict[i]=BooleanVar()
                self.face_var_dict[i].set(False)
                
                self.face_chk_dict[i]=Checkbutton(self.shot_frames[shot], text=disp_face_text, var=self.face_var_dict[i])
                self.face_chk_dict[i].pack(side=LEFT)
                
            self.shot_face_var_dict[shot]=self.face_var_dict
            self.shot_face_chk_dict[shot]=self.face_chk_dict
                
                #globals()[face_text]=BooleanVar()
                #globals()[face_text].set(False)
                #globals()[var_string]=Checkbutton(self.shot_frames[shot], text=disp_face_text, var=globals()[face_text])
                #globals()[var_string].pack(side=LEFT)
                
        Data_Chopping_Button=Button(self.root, text="Chop Data", command=self.chop_data).pack(side=TOP)
        
        #================================================
        
        #================================================
        #Optimization Parameters Section
        opt_params_frm=ttk.Frame(self.root, padding=10)
        opt_params_frm.pack(side=TOP)

        Label(opt_params_frm, text="===========================").pack()
        Label(opt_params_frm, text="Optimization Parameters").pack()
        Label(opt_params_frm, text=" ").pack()

        #Folder Path
        fp_frm=ttk.Frame(opt_params_frm)
        fp_frm.pack(side=TOP)
        self.folder_path=StringVar()
        self.folder_path.set(self.global_vars['optimization_data_path'])
        Label(fp_frm, text="Root Folder").pack(side=LEFT)
        Entry(fp_frm, textvariable=self.folder_path).pack(side=LEFT)

        #Optimization run name
        rn_frm=ttk.Frame(opt_params_frm)
        rn_frm.pack(side=TOP)
        self.run_name=StringVar()
        self.run_name.set("Default")
        Label(rn_frm, text="Run Name ").pack(side=LEFT)
        Entry(rn_frm, textvariable=self.run_name).pack(side=LEFT)

        #Number of steps
        ns_frm=ttk.Frame(opt_params_frm)
        ns_frm.pack(side=TOP)
        self.meshpoints=IntVar()
        self.meshpoints.set(60)
        Label(ns_frm, text="Meshpoints ").pack(side=LEFT)
        Entry(ns_frm, textvariable=self.meshpoints).pack(side=LEFT)

        #Population size
        ps_frm=ttk.Frame(opt_params_frm)
        ps_frm.pack(side=TOP)
        self.popsize=IntVar()
        self.popsize.set(5)
        Label(ps_frm, text="Pop Size     ").pack(side=LEFT)
        Entry(ps_frm, textvariable=self.popsize).pack(side=LEFT)
        #======================================
        
        #======================================
        #Bounds
        #======================================
        bounds_frm=ttk.Frame(self.root, padding=10)
        bounds_frm.pack(side=TOP)

        Label(bounds_frm, text="===========================").pack()
        Label(bounds_frm, text="Bounds").pack()

        #a
        a_frm=ttk.Frame(bounds_frm)
        a_frm.pack(side=TOP)
        self.a_bounds=StringVar()
        self.a_bounds.set('( 0.01 , 0.05 )')
        Label(a_frm, text="a:                ").pack(side=LEFT)
        Entry(a_frm, textvariable=self.a_bounds).pack(side=LEFT)

        #b
        b_frm=ttk.Frame(bounds_frm)
        b_frm.pack(side=TOP)
        self.b_bounds=StringVar()
        self.b_bounds.set('( 0 , 100 )')
        Label(b_frm, text="b:                ").pack(side=LEFT)
        Entry(b_frm, textvariable=self.b_bounds).pack(side=LEFT)

        #start_time
        st_frm=ttk.Frame(bounds_frm)
        st_frm.pack(side=TOP)
        self.st_bounds=StringVar()
        self.st_bounds.set('( 1.75E-8 , 2.25E-8 )')
        Label(st_frm, text="Start time: ").pack(side=LEFT)
        Entry(st_frm, textvariable=self.st_bounds).pack(side=LEFT)

        #peak_temp
        pt_frm=ttk.Frame(bounds_frm)
        pt_frm.pack(side=TOP)
        self.pt_bounds=StringVar()
        self.pt_bounds.set('( 30000 , 100000 )')
        Label(pt_frm, text="Peak temp: ").pack(side=LEFT)
        Entry(pt_frm, textvariable=self.pt_bounds).pack(side=LEFT)
        #===============================================
        
        #===============================================
        #Buttons
        #===============================================
        Label(self.root, text="===========================").pack(side=TOP)
        #Data_Chopping_Button=Button(self.root, text="Chop Data", command=self.chop_data).pack(side=BOTTOM)
        Vis_Button=Button(self.root, text="Visualize Optimization", command=self.visualize_optimization).pack(side=BOTTOM)
        Run_Button=Button(self.root, text="Run Optimization", command=self.show_run_dict).pack(side=BOTTOM)
        Log_Button=Button(self.root, text="Log comments", command=self.log_comments).pack(side=BOTTOM)
        #===============================================
        
        self.root.mainloop()
    
class PrimaryWindow:
    
    def donothing(self):
        filewin = Toplevel(self.root)
        button = Button(filewin, text="Do nothing button")
        button.pack()
    
    def open_new_optimization_window(self):
        new_window=NewOptimizationWindow(self.root)
    
    def __init__(self):
        self.root=Tk()
        self.root.title('Optimization Model - Main Window')
        self.root.geometry('200x200')
        #=============================
        #Menu
        #=============================
        menubar = Menu(self.root)
        filemenu = Menu(menubar, tearoff=0)
        editmenu = Menu(menubar, tearoff=0)
        filemenu.add_command(label="New", command=self.open_new_optimization_window)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=self.root.quit)
        menubar.add_cascade(label="File", menu=filemenu)
        menubar.add_cascade(label="Edit", menu=editmenu)
        helpmenu = Menu(menubar, tearoff=0)
        helpmenu.add_command(label="Help Index", command=self.donothing)
        helpmenu.add_command(label="About...", command=self.donothing)
        menubar.add_cascade(label="Help", menu=helpmenu)
        #==============================
        self.root.config(menu=menubar)
        self.root.mainloop()
        
PrimaryWindow()