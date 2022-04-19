'''
Main file for Thermal Laser Compression Model
Pat LaChapelle
Created: February 10, 2022
Last updated: March 9, 2022
'''
from selectors import EpollSelector
from tkinter import E
import fenics
fenics.set_log_level(30) #surpress fenics log except for warnings
import numpy as np
import matplotlib
matplotlib.use('WebAgg')
from matplotlib import pyplot as plt
import json
import pandas as pd
import math
from scipy import interpolate
from matplotlib import animation
from scipy.optimize import differential_evolution
from threading import Thread
import os
import errno
import datetime
from pprint import pprint
from matplotlib.figure import Figure
import pickle
import dill
from matplotlib.pyplot import figure
from datetime import datetime

default_best_SLS=10e+20 #Set initial best sls very large

default_parameters=[0.03, 30, 1.8*10**-8, 30000] #Define some default parameters for passing to function files when convenient
start_stop_list=[               #s88773
                 [[179, 240],       # face 1
                  [259, 390],       # face 2
                  [310, 700]],      # face 3 
                                #s88776
                 [[88, 155],        # face 1
                  [93, 400],        # face 2
                  [246, 468]],      # face 3
                                #s88780
                 [[145, 190],       # face 1
                  [240, 350],       # face 2
                  [310, 502]],      # face 3
                                #s86483
                 [[355, 400],       # face 1
                  [582, 732],       # face 2
                  [377, 800]]]      # face 3

start_stop_dict={
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

def least_squares_interp(x0,y0,x1,y1):
    #x0: the independent variable data for the model output
    #y0: the dependent variable data for the model output
    #x1: independent variable data for results
    #y1: dependent variable data for results
    f=interpolate.interp1d(x0,y0, bounds_error=False, fill_value='extrapolate')
    interpolated_data=f(x1)

    least_squares=np.sum(np.square(y1-interpolated_data))

    return least_squares


def view_optimization(fp, save_plot=False):
    
    fig=Figure()
    ax1=fig.add_subplot(2,3,1) #sls
    ax2=fig.add_subplot(2,3,2) #peak_temp
    ax3=fig.add_subplot(2,3,3) #a
    ax4=fig.add_subplot(2,3,4) #b
    ax5=fig.add_subplot(2,3,5) #start_time
    
    it=[]
    sls=[]
    a=[]
    b=[]
    peak_temp=[]
    start_time=[]
    pullData=open(fp, "r").read()
    data=pullData.split('\n')
    clean_data=[]
    for i in data:
        line=i.replace('[',"").replace(']',"")
        if len(line)>0:
            clean_data.append(line)
            _it, _sls, _a, _b, _start_time, _peak_temp = line.split(',')
            print(clean_data)
            it.append(float(_it))
            sls.append(float(_sls))
            a.append(float(_a))
            b.append(float(_b))
            peak_temp.append(float(_peak_temp))
            start_time.append(float(_start_time))
            
    plt.show()
    
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
    ax1.plot(it,sls)
    ax2.plot(it, a)
    ax3.plot(it, b)
    ax4.plot(it, start_time)
    ax5.plot(it, peak_temp)

        
def load_json(fname): #function for loading in json files
    with open(fname) as json_file: #assign opened file to object json_file
        data=json.load(json_file) #load data
    return data

def boundary ( x ): #Boundary function used for finite element model
    value = x[0] < fenics.DOLFIN_EPS #Dirichlet boundary condition 
    return value

class RawPlot:
    def __init__(self, data_dictionary):
        self.run_output=data_dictionary
        self.run_IDs=list(self.run_output.keys())
        
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
                ax[row_count,col_count].set_title(run)
            fig_count=fig_count+1
            col_count=fig_count%nrows
            row_count=int(fig_count/nrows)
        
        plt.show()

class Data: # Class for data objects. Should hold all data that we are fitting
    def __init__(self): # initialization function
        self.data=load_json('data_dict.json') #load in the data from the data dictionary
        self.metadata={'s88773': {'filepath':'../Data/s88773_sop_lineouts_TP.xlsx',
                              'sheet_name':'June2021',
                              'faces':[1,2,3]}, 
                   's88776': {'filepath':'../Data/s88776_SOP_TP.xlsx',
                              'sheet_name':'Sheet1',
                              'faces':[1,2,3]},
                   's88780': {'filepath':'../Data/s88780_SOP_TP.xlsx',
                              'sheet_name':'Sheet1',
                              'faces':[1,2,3]},
                   's86483': {'filepath':'../Data/s86483_SOP_TP.xlsx',
                              'sheet_name':'Sheet1',
                              'faces':[1,2,3]},
                   's86480': {'filepath':'../Data/s86480_SOP_2022.xlsx',
                              'sheet_name':'Sheet1',
                              'faces':[1,2,3]
                              },
                   's86484': {'filepath':'../Data/s86484_SOP.xlsx',
                              'sheet_name':'Sheet1',
                              'faces':[1,2,3]
                              },
                   's88774': {'filepath':'../Data/s88774_SOP.xlsx',
                              'sheet_name':'Sheet1',
                              'faces':[1,2,3]
                              }                                      
                   }
        
    def get_data(self, start_stop:dict=start_stop_dict):
        self.start_stop=start_stop
        SOP_data={}
        count=0
        for shot in self.metadata.keys():
            fp=self.metadata[shot]['filepath']
            sheet=self.metadata[shot]['sheet_name']
            shot_data=pd.read_excel(fp, sheet_name=sheet)
            
            for key in shot_data.keys():
                if key=='Time':
                    my_key='Time'
                elif key=='time':
                    my_key='time'

            face_data={}
            for face in self.metadata[shot]['faces']:
                face_key='face'+str(face)
                key='step'+str(face)+'_corrected'
                start=start_stop[shot][face][0]
                stop=start_stop[shot][face][1]
                #face_data[face_key]=[(shot_data[my_key][start_stop[count][face-1][0]:start_stop[count][face-1][1]].to_numpy()*10**-9).astype(float), (shot_data[key][start_stop[count][face-1][0]:start_stop[count][face-1][1]].to_numpy()).astype(float)]
                face_data[face_key]=[(shot_data[my_key][start:stop].to_numpy()*10**-9).astype(float), (shot_data[key][start:stop].to_numpy()).astype(float)]
            SOP_data[shot]=face_data
            count=count+1
        self.data=SOP_data
            
    def correct_SOP(self, shot, xw=1, aeta=25.5, a0=481000, t0a=1.909, reference=0.55 ):
        #for shot in self.data.keys():
        
        #correctedSOP={}
        
        shot_dict=self.data[shot]
        face_data={}
        for face in shot_dict.keys():
            SOP_data=shot_dict[face][1]
            newSOP=np.real(11605*t0a)/(np.log(1+((1-reference)*a0/(aeta*np.abs(SOP_data)))))
            face_data[face]=shot_dict[face]
            face_data[face][1]=newSOP
        return face_data
        #self.data=correctedSOP
        
    def correct_all_SOP(self):
        correctedSOP={}
        for shot in self.data.keys():
            if shot=='s88774':
                corrected_data=self.correct_SOP('s88774', t0a=1.909, a0=481000)
            else:
                corrected_data=self.correct_SOP(shot)
            correctedSOP[shot]=corrected_data
        self.data=correctedSOP
                
                
    def plot_data(self, fname):
        plot_obj=RawPlot(self.data)
        plot_obj.plot()
        
class ComparitivePlot:
    def __init__(self, model_data, real_data, savefile=None):
        self.savefile=savefile
        self.run_output=real_data
        self.model_output=model_data
        self.run_IDs=list(self.model_output.keys())
        print(self.run_IDs)
        
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
        fig, ax = plt.subplots(nrows, ncols, squeeze=False)
        print(nrows, ncols)
        row_count=0
        col_count=0
        fig_count=0
        print(nrows, ncols)
        print(self.model_output[self.run_IDs[0]])
        for run in self.run_IDs:
            for face in list(self.model_output[run].keys()):
                print(self.model_output[run].keys())
                ax[row_count,col_count].plot(self.run_output[run][face][0], self.run_output[run][face][1])
                ax[row_count,col_count].plot(self.model_output[run][face][0], self.model_output[run][face][1])
                ax[row_count,col_count].set_title(run)
                print(run, face, fig_count, row_count, col_count)
            fig_count=fig_count+1
            col_count=fig_count%nrows
            row_count=int(fig_count/nrows)
        
        if self.savefile==None:
            plt.show()
        else:
            plt.savefig(self.savefile)
            
    def plot2(self):
        for shot in self.model_output.keys():
            shot_dict=self.model_output[shot]
            for face in shot_dict.keys():
                #Add each face to the plot
                plt.plot

class SingleFaceModel:
    def __init__(self, data, tf=30*10**-9, num_steps=120, init_temp=2000, k_1=100, rho=12000, c=450, MgO_length=2*10**-6, Fe_length=1.07*10**-6, meshpoints=120, equation=1):
        self.tf, self.data, self.num_steps, self.init_temp, self.k_1, self.rho, self.c1, self.MgO_length, self.Fe_length, self.meshpoints = tf, data, num_steps, init_temp, k_1, rho, c, MgO_length, Fe_length, meshpoints #assign values to self so we can reference them later as needed
        self.dt=self.tf/self.num_steps #time step (dt) is total time divided by the number of time steps
        #Handle the different equations
        if equation==1:
            self.expression='x[0] <= l + tol ? b+a*u : k_1'
        elif equation==2:
            self.expression='x[0] <= l + tol ? b+a*u+d*pow((u+1),-0.5) : k_1'
        elif equation==3:
            self.expression='x[0] <= l + tol ? b+a*u+d*pow((u+1),2) : k_1'
        elif equation==4:
            self.expression='x[0] <= 1 + tol ? b+a*u+d*pow((u+1), 0.5) : k_1'
        self.equation=equation
    
        if equation==1:
            self.default_parameters=[0.03, 30, 1.8*10**-8, 30000]
        elif equation==2 or equation==3 or equation==4:
            self.default_parameters=[0.03, 30, 1, 1.8*10**-8, 30000]
    
    def run(self, parameters): # Function for running the single face model
        if self.equation==1:
            a, b, start_time, peak_temp = parameters #get values from parameters
            self.a, self.b, self.start_time, self.peak_temp=a,b,start_time,peak_temp #assign values to self
            self.my_mesh = fenics.IntervalMesh(self.meshpoints, 0, self.Fe_length+self.MgO_length) #Use Fenics to create a mesh for the finite element model
            self.V = fenics.FunctionSpace(self.my_mesh, 'P', 1) #create a function space for the finite element model
            self.u_D = fenics.Expression(str(self.peak_temp),degree=1) #Define the peak temperature of the square wave boundary condition
            self.bc = fenics.DirichletBC(self.V, self.u_D, boundary) #Define the Dirichlet boundary condition as variable bc
            self.u_0 = fenics.Expression(str(self.init_temp), degree=1) #Define the initial value
            self.u_n = fenics.interpolate(self.u_0, self.V)  # ||
            self.u = fenics.Function(self.V) #Define variational problem
            self.v = fenics.TestFunction(self.V) # ||
            self.kappa = fenics.Expression(self.expression, degree=1, tol=fenics.DOLFIN_EPS, k_1=self.k_1, u=self.u, a=a, b=b, l=self.Fe_length) # ||
            inc = 0.01*10**-6 #Define some small distance to measure off of the boundary, so that we are not pulling the data from exactly on the boundary
            loc = self.Fe_length - inc # measuring location
            self.times = np.linspace(0,self.tf, self.num_steps) # Get list of times for the model. This will be the independent variable of the model output
            self.times = (np.asarray(self.times)+start_time).tolist() # shift times to account for the start time of the model
            self.time_line = [] # Empty list for getting the model output
            t=0 #start time @ 0 (will be shifted by start time from line 46)
            self.F = self.u*self.v*fenics.dx + self.kappa/(self.rho*self.c1)*self.dt*fenics.dot(fenics.grad(self.u), fenics.grad(self.v))*fenics.dx - (self.u_n)*self.v*fenics.dx # Defining variational problem
            for n in range(self.num_steps): # Loop through model to get output
                t+=self.dt # time step through model
                self.u_D.t=t # update current time
                fenics.solve(self.F==0, self.u, self.bc) # Solve the variational problem for the current time to get the Temperature in terms of the position
                self.time_line.append(self.u(loc)) # Get the Temperature @ the measurement location & add it to the output data
                self.u_n.assign(self.u) # update previous solution
                
        elif self.equation==2 or self.equation==3 or self.equation==4:
            a, b, c, start_time, peak_temp = parameters
            self.a, self.b, self.c, self.start_time, self.peak_temp = a, b, c, start_time, peak_temp
            self.my_mesh=fenics.IntervalMesh(self.meshpoints, 0, self.Fe_length+self.MgO_length)
            self.V = fenics.FunctionSpace(self.my_mesh, 'P', 1)
            self.u_D = fenics.Expression(str(self.peak_temp),degree=1)
            self.bc = fenics.DirichletBC(self.V, self.u_D, boundary)
            self.u_0 = fenics.Expression(str(self.init_temp), degree=1)
            self.u_n = fenics.interpolate(self.u_0, self.V)
            self.u = fenics.Function(self.V)
            self.v = fenics.TestFunction(self.V)
            self.kappa = fenics.Expression(self.expression, degree=1, tol=fenics.DOLFIN_EPS, k_1=self.k_1, u=self.u, a=a, b=b, d=c, l=self.Fe_length)
            inc = 0.01*10**-6
            loc = self.Fe_length - inc
            self.times = np.linspace(0,self.tf, self.num_steps)
            self.times = (np.asarray(self.times)+start_time).tolist() # shift times to account for the start time of the model
            self.time_line = []
            t=0
            self.F = self.u*self.v*fenics.dx + self.kappa/(self.rho*self.c1)*self.dt*fenics.dot(fenics.grad(self.u), fenics.grad(self.v))*fenics.dx - (self.u_n)*self.v*fenics.dx
            for n in range(self.num_steps):
                t+=self.dt 
                self.u_D.t=t 
                fenics.solve(self.F==0, self.u, self.bc)
                self.time_line.append(self.u(loc))
                self.u_n.assign(self.u)
            
            
            
    def get_chi_squared(self):
        chi_2=least_squares_interp(self.times,self.time_line,self.data[0], self.data[1])
        return chi_2
            
    def plot(self): # Function for getting the plot of the model's output
        plt.clf() # clear previous data from the plot
        plt.plot(self.times, self.time_line) # plot data output for the model
        plt.show() # show the plot
        
class Geometry: # Geometry class for storing geometry data of the shot
    
    def __init__(self, runID:str, facelist, dictionary_filepath='geometry_dict.json'): # Initialization function requires specification of a single shot name as input runID (type:string), and a list of the faces to get, along with the filename to go look at for the geometry. Filename is defaulted as geometry_dict.json. Eventually, I will add a function for reading in the data and saving it to self
        data=load_json(dictionary_filepath)[runID] # Load the geometry dictionary for the run that has been specified
        new_dict={} # assign a new blank dictionary to write the geometry data to
        for face in data: # iterate through each face
            for id in facelist: # look through which faces are specified in the facelist
                name='face'+str(id) # Get what the name is for the corresponding face
                if name==face: # If the face is the one we are looking for, pull the geometry data & store in new_dict accordingly
                    new_dict[face]=data[face] # ||
        self.geo=new_dict # Assign geometry  dictionary to property geo

class FiniteElementModel: #Class for the finite element model of a single shot
    
    def __init__(self, meshpoints, data, geometry_dict:Geometry, equation:int=1): # Initialize an object of class FiniteElementModel, defining the number of meshpoints desired
        vals=load_json('default_values_dict.json') # load in default values from the dictionary of default values
        self.tf, self.num_steps, self.MgO_length, self.init_temp, self.k_1, self.rho, self.c1, self.meshpoints = vals['tf'], vals['num_steps'], vals['MgO_length'], vals['init_temp'], vals['k_1'], vals['rho'], vals['c'], meshpoints #assign values
        self.dt=self.tf/self.num_steps # time step
        self.data=data #assign data to self
        self.faces=geometry_dict.geo.keys() # Get the names of the faces that we need to build finite element models for
        face_count=1 # first face is always 1, so write face_count to 1 before iterating through the faces
        self.geometry=geometry_dict.geo # assign the geometry dictionary to a self object
        models={} # Write new blank dictionary named models for storing the model that has been created for each face
        for face in self.faces: # iterate through the faces
            key='face'+str(face_count) # Get the key for the face in the geometry dictionary
            Fe_length=self.geometry[face] # Get the finite element length for the face from the geometry dictionary
            models[key]=(SingleFaceModel(self.data[face], self.tf, self.num_steps, self.init_temp, self.k_1, self.rho, self.c1, self.MgO_length, Fe_length, self.meshpoints, equation=equation)) # Write the model into the model dictionary as the values associated with the key for the face
            face_count=face_count+1 # increment face_count by 1
        self.face_models=models # Now assign models dictionary as a property of self so that we can reference it later 
        
    def run_model(self, parameters): # Method for running the model that has been defined by the initialization of the object of class FiniteElementModel
        results={} # Create empty dictionary for storing the results of the model run
        for face in self.face_models.keys(): # Iterate through the faces
            model=self.face_models[face] # Get the model that is associated with this face
            model.run(parameters) #run the model
            results[face]=model.times, model.time_line #store the model results
        self.model_results=results
            
    def plot_model(self): # Define function for plotting the current output of the model
        plt.clf() # Clear plot
        for face in self.face_models.keys(): # iterate through the faces
            model=self.face_models[face] # Get the model for the face
            plt.plot(model.times, model.time_line) # Plot the model output data from the model
        plt.show() # Display the plot
        
class OptimizationGroup: # A class for an optimization group. An optimization group can contain several runs and several faces, which are specified in the run_dictionary
    
    def __init__(self, meshpoints, run_dictionary:dict, name:str, data:Data, equation:int=1): # Initialize an OptimizationGroup object. run_dictionary should have keys for each runID (ie 's88776'). Associated value for each runID key is a list of faces (ie [1,2,3]). Also takes number of meshpoints to use for the finite element model of each run and face
        #================================================
        #Assigning input variables to self & getting the temperature data + geometry for the runs
        self.run_dictionary=run_dictionary
        self.models={}
        #self.full_shot_data=Data() # get data
        #self.full_shot_data.get_data()
        self.meshpoints=meshpoints
        self.full_shot_data=data
         
        self.full_shot_data.correct_all_SOP() # correct SOP data
        self.real_data=self.full_shot_data.data
        self.name=name
        for run in run_dictionary.keys(): # iterate through the different runs
            run_geometry=Geometry(run, run_dictionary[run]) # Get geometry for the run and faces
            f=FiniteElementModel(meshpoints, self.real_data[run], run_geometry, equation=equation) # Create a Finite element model for each run in run_dictionary
            self.models[run]=f

        #==============================================
        
        #==============================================
        #Now for Setting up the optimization environment
        global_vars=load_json('global_variables.json')
        root_folder_path=global_vars['optimization_data_path']
        try:
            os.makedirs(root_folder_path+'/'+name)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        self.folder_path=root_folder_path+'/'+name
        
        #============================================
        #Saving models as pickle object
        '''with open(self.folder_path+'/test_pick.pickle', 'wb') as saveg:
            pickle.dump([self.models], saveg, protocol=pickle.HIGHEST_PROTOCOL)
        '''#============================================
        
        self.optimization_data_path=self.folder_path+'/optimization_output.csv'
        
        #==============================================
        self.equation=equation
        
            
    def run(self, params): #Write function for running the model for each shot
        model_data={}
        for run in self.models.keys(): # iterate through each shot
            model=self.models[run] # Get the model for the run
            model.run_model(params) # run the model for the given parameters
            model_data[run]=model.model_results
        self.model_data=model_data
            
    def get_chi_sq(self):
        chi_sq=0
        for shot in self.run_dictionary.keys(): #iterate through the shots specified in the run dictionary
            model=self.models[shot]
            for face in model.face_models.keys(): #iterate through each face being fit for this shot
                face_model=model.face_models[face] # Get the model for this face
                face_chi_sq=face_model.get_chi_squared() # Get the chi^2 value for this face
                chi_sq=chi_sq+face_chi_sq
        return chi_sq
    
    def plot_optimization_group(self, fname=None):
        if fname==None:
            self.my_plot=ComparitivePlot(self.model_data, self.real_data)  
        else:
            self.my_plot=ComparitivePlot(self.model_data, self.real_data, savefile=fname)
        self.my_plot.plot()      
    
    def optimizable_function(self, params):
        self.run(params)
        chi_sq = self.get_chi_sq()
        return chi_sq
    
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
            
                        
            ani=animation.FuncAnimation(fig, animate, interval=1000)
            #ani.save(gif_file, writer='imagemagick',fps=1)
            plt.show()
    
    def optimize(self, bounds, popsize):
        def optimizable_function(params):
            self.run(params)
            chi_sq = self.get_chi_sq()
            return chi_sq
        diffEV=differential_evolution(optimizable_function, bounds=bounds, popsize=popsize, disp=True)



'''        self.run(default_parameters) #run default parameters to ensure model_results attribute exists
        self.my_plot=ComparitivePlot(self.model_data, self.real_data)  
        thread1=Thread(target=differential_evolution, args= (optimizable_function, bounds, popsize))
        thread2=Thread(target=self.animate_plots())
        #threading.start_new_thread(self.my_plot.animate_plots())
        thread1.start()
        thread2.start()
        thread1.join()
        thread2.join()'''
        
class OptimizationData:
    def __init__(self, filepath):
        self.optimization_data=[]
        self.filepath=filepath
        
    def get_current_vals(self, current_vals):
        #Current vals should be represented as:
        # iteration no., sls, a, b, start_time, peak_temp
        self.optimization_data.append(current_vals)
        
    def show_vals(self):
        print(self.optimization_data)
        

def optimize(group:OptimizationGroup, data_obj:OptimizationData, bounds, popsize, equation:int=1):
    #==================================
    #Write metadata file
    #==================================
    def bounds_to_bounds_dict(bounds, equation=equation):
        if equation==1:
            keys=['a', 'b', 'start_time', 'peak_temp']
            bounds_dict={}
            count=0
            for item in bounds:
                bounds_dict[keys[count]]=item
                count=count+1
        elif equation==2 or equation==3 or equation==4:
            keys=['a', 'b', 'c', 'start_time', 'peak_temp']
            bounds_dict={}
            count=0
            for item in bounds:
                print(count)
                print(keys[count])
                bounds_dict[keys[count]]=item
                count=count+1
        return bounds_dict
                
    metadata_fp=group.folder_path+'/'+'metadata.txt'
    if equation==1:
        equation_string='k = a * T + b'
    elif equation==2:
        equation_string='k = a*T + b + c/(sqrt(T))'
    elif equation==3:
        equation_string='k = a*T + b + c*T^2'
    elif equation==4:
        equation_string='k = a*T + b + c*sqrt(T)'    
        
    bounds_dict=bounds_to_bounds_dict(bounds, equation)
        
    with open(metadata_fp, 'wt') as out:
        out.write(group.name+': Metadata')
        out.write('\n')
        out.write('\n')
        out.write('===========================')
        out.write('\n')
        out.write('Optimization Parameters: ')
        out.write('\n')
        out.write('Equation: '+equation_string)
        out.write('\n')
        out.write('Optimization initialized @: ')
        out.write('\n')
        out.write(str(datetime.now()))
        out.write('\n')
        out.write('popsize: '+str(popsize))
        out.write('\n')
        out.write('Meshpoints: '+str(group.meshpoints))
        out.write('\n')
        out.write('Optimization type: Differential Evolution')
        out.write(' ')
        out.write('\n')
        out.write('===========================')
        out.write('\n')
        out.write('Runs & Faces: ')
        out.write('\n')
        pprint(group.run_dictionary, stream=out, width=1)
        out.write(' ')
        out.write('\n')
        out.write('==========================')
        out.write('\n')
        out.write('Bounds: ')
        out.write('\n')
        count=0
        bounds_names=list(bounds_dict.keys()) #['a', 'b', 'start_time', 'peak_temp']
        for item in bounds:
                out.write(bounds_names[count]+': '+str(item))
                out.write('\n')
                count=count+1
        out.write('\n')
        out.write('===========================')
        out.write('\n')
        out.write('Data Chopping: ')
        out.write('\n')
        out.write('\n')
        for real_data_shot_name in group.real_data.keys(): # iterate through the shot names
                for run_dict_shot_name in group.run_dictionary.keys(): # iterate through the shots specified in the run dictionary
                        if real_data_shot_name==run_dict_shot_name: #If the shot from the real data is in the run dictionary
                                out.write(str(run_dict_shot_name)+': ')
                                out.write('\n')
                                data=group.full_shot_data.start_stop[run_dict_shot_name]
                                for face in group.run_dictionary[run_dict_shot_name]: # Get the faces for each shot
                                        pretty_name='face'+str(face)
                                        out.write('     '+pretty_name+': '+str(data[face]))
                                        out.write('\n')
                                out.write('\n')    
    
    #==================================
    
    #==========================================
    #Run Optimization
    #==========================================
    def optimizable_function(params:list):
        has_best_fit=False
        print(params)
        group.run(params) #Run the model with the given parameters
        least_squares=group.get_chi_sq() #Get least squares from the run with these parameters
        iteration=len(data_obj.optimization_data) #Get the iteration. Since it starts @ 0, iteration number should just be the length of the optimization data
        current_data=[iteration, least_squares]
        current_data.extend(params) #Write list for the current data
        if iteration==0:
            best_fit=current_data
            has_best_fit=True
            group.plot_optimization_group(fname=group.folder_path+'/Initial_Fit.png')
        else:
            #If the iteration is >0, start by getting the pickled best fit
            best_fit_filepath=group.folder_path+'/Best_Fit.pickle'
            best_fit_file=open(best_fit_filepath, 'rb')
            best_fit=pickle.load(best_fit_file)
            
            #Now check if current sls is an improvement on the current best
            current_best_sls=best_fit[1]
            if least_squares<current_best_sls:
                has_best_fit=True
                best_fit=current_data
            else:
                has_best_fit=False
                
        data_obj.get_current_vals(current_data) #Assign data from this iteration to the data object 
        print(iteration, current_data)
        print(group.optimization_data_path)
        with open(group.optimization_data_path, 'a+') as file_object:
            file_object.write(str(current_data))
            file_object.write('\n')
        
        #Save the best fit as an obj
        best_fit_filepath=group.folder_path+'/Best_Fit.pickle'
        best_fit_file=open(best_fit_filepath, "wb")
        pickle.dump(best_fit, best_fit_file)
        best_fit_file.close()
        
        #Plot & save if this is the best fit of the run
        if has_best_fit: #If this iteration is the current best fit...
            group.plot_optimization_group(fname=group.folder_path+'/Best_Fit.png')
        return least_squares    
    
    #=============================
    #Using Pickle to save our group to a file so that we may re-access it again at a time of our choosing
    
    print(dir(group))
    for i in dir(group):
        if not i.startswith('__'):
            print(i)
            i
    for v in vars(group):
        print(v, type(v))
        
    group_s=group
        
    group_obj_file_path=group.folder_path+'GROUP_OBJ.pickle'
    
    '''    with open(group_obj_file_path, 'wb') as save_group:
        print(type(group))
        pickle.dump(list(group), save_group)
        '''
    
    #==============================
    
    group.DiffEV=differential_evolution(func=optimizable_function, bounds=bounds, popsize=popsize)
    
    return group
    #===========================================
    
class PostOptOperations: #A class for writing and defining operations to be done on the optimization after it has completed
    def __init__(self, folder_path):
        self.has_opt_group=False #Variable for whether the optimization group exists for these post optimization operations
        f=open(folder_path+'/metadata.txt', 'r')
        equation_string=f.readlines()[4][10:][:-1]
        print(equation_string)
        f=open(folder_path+'/metadata.txt', 'r')
        self.metadata=f.readlines()
        if equation_string=='k = a * T + b':
            equation=1
        elif equation_string=='k = a*T + b + c/(sqrt(T))':
            equation=2
        elif equation_string=='k = a*T + b + c*T^2':
            equation=3
        elif equation_string=='k = a*T + b + c*sqrt(T)':
            equation=4
        else:
            print('------------------')
            print('Equation String: ')
            print(equation_string)
            print('Equation not FOUND')
            print('------------------')
   
        self.equation=equation
        
        self.fp=folder_path
        self.data_path=folder_path+'/optimization_output.csv'
        group_obj_path=folder_path+'/GROUP_OBJ.pickle'
        files_in_path=os.listdir(folder_path)
        self.hasGroup=False
        self.hasCleanData=False
        for file in files_in_path:
            if file=='GROUP_OBJ.pickle':
                self.hasGroup=True
                group_f=open(group_obj_path, 'rb')
                self.group=pickle.load(group_f)
                group_f.close() 
                print('==============================')
                print(str(self.group.equation))
                print('==============================')
            if file=='Clean_Opt_Data.csv':
                self.hasCleanData=True
                
    def clean_data(self):
        print(self.hasCleanData)
        if self.hasCleanData==False:
            pullData=open(self.data_path, 'r').read()
            dataArray=pullData.split('\n')
            print(dataArray[0].split(','), len(dataArray[0].split(',')))
            if len(dataArray[0].split(','))==6:
                clean_lines=['iteration,sls,a,b,start_time,peak_temp']
            elif len(dataArray[0].split(','))==7:
                clean_lines=['iteration,sls,a,b,c,start_time,peak_temp']
            else:
                if len(dataArray[1])==6:
                    clean_lines=['iteration,sls,a,b,start_time,peak_temp']
                elif len(dataArray[1])==7:
                    clean_lines=['iteration,sls,a,b,c,start_time,peak_temp']
                else:
                    print('LINES 1 & 2 ARE BLANK')
            for eachLine in dataArray:
                line=eachLine.replace('[','').replace(']','')
                if len(line)>0:
                    clean_lines.append(line)
                    print(clean_lines)
            clean_data_path=self.fp+'/Clean_Opt_Data.csv'       
            with open(clean_data_path, 'a+') as file_object:
                for line in clean_lines:
                    file_object.write(str(line))
                    file_object.write('\n')        
            self.hasCleanData=True
            return clean_lines
        else:
            print('DATA ALREADY EXISTS')
        
    def get_run_dict(self):
        should_print=False
        run_dict_lines=[]
        for l in self.metadata:
            if should_print==True:
                #print(l)
                run_dict_lines.append(l)
            if l == 'Runs & Faces: \n':
                should_print=True
            if l == ' \n':
                should_print=False
    
        sl=''
        for line in run_dict_lines:
            sl=sl+line
    
        in_dict=False #Boolean val for whether we are inside a dictionary
        is_key=False #Boolean val for whether the current char is part of a dictionary key
        looking_for_val=False #Boolean val for whether we are currently looking for the val in a key/val pair
        is_val=False #Boolean for telling us whether the current char is part of a string associated w/ a value
        is_num=False #Boolean for telling whether 
        has_key=False #Boolean indicator for whether or not we currently have a dictionary key
        current_vals=[] #Empty list for storing initial vals to first key
        my_dict={}
        for char in sl: #Iterate through each char in the single line

            #=====================================================================================
            #Detecting whether we are in a dictionary

            if char=="{": # { is the character denoting the beginning of a dictionary
                in_dict=True
            elif char=="}": # } is the char denoting the closing/ending of a dictionary
                in_dict=False
            #=====================================================================================

            #=====================================================================================
            #Detecting whether current char is part of a dictionary key

            if in_dict==True: #We only care to find keys if we are currently in a dictionary
                if char=="'": #Quotation marks signal the end/beginning of a key
                    if is_key==True: #If we're in a key, this signals the end of the key
                        is_key=False
                        looking_for_val=True #Now that we've found the key, we're looking for associated values
                        has_key=True
                    else: #If we are not in a key, this signals the beginning of a key
                        is_key=True
                        key_string='' #Blank string for saving the key once we have each char associated w it
                else:
                    if is_key==True: #If we are currently getting a key, add this character to the key string
                        key_string=key_string+char

            #=====================================================================================

            #=====================================================================================
            #Getting Values for the current key

            if looking_for_val==True: #If we're looking for a value
                if is_val==True: #If we've detected a val
                    if is_num==False: #If we have detected a val, but not yet a number
                        if char.isnumeric()==True: #If our character is a number
                            is_num=True #We are currently at a number, which is part of the val, so set is_num to True
                            string_num=char #Set a string character of our current number so we can access it when we loop over the next character
                    else: #otherwise, if we are currently looking at a number
                        if char.isnumeric()==True: #There is another char in the number, so append this value
                            string_num=string_num+char
                        else: #If the character is not numeric
                            if char=='.': #If the non-numeric character is a decimal, then we are still in a #
                                string_num=string_num+char 
                            else: #If the string is not part of a number...
                                is_num=False
                                current_vals.append(int(string_num)) #write integer form of string_num to current_vals
                if char=='[': #If the character is an opening bracket, this is the beginning of our values
                    is_val=True #Change is_val to True so that we know we are in a value string
                if char==']': #This signals the end of a set of values
                    is_val=False
                    looking_for_val=False
                    my_dict[key_string]=current_vals #Write the new entry into the dictionary
                    current_vals=[] #Reset current_vals to be an empty list
                    has_key=False
                    
        return my_dict
    
    def get_opt_group(self):
        if not self.has_opt_group: #If we do not already have the opt group
            run_dict=self.get_run_dict()
            run_data=Data()
            run_data.get_data()
            self.group=OptimizationGroup(120, run_dict, 'PostOptTester', run_data, equation=self.equation)
            self.has_opt_group=True    
    
    def get_initial_fit(self):
        pass
    
    def get_best_fit(self):
        self.get_opt_group()
        self.clean_data()
        q=pd.read_csv(self.fp+'/Clean_Opt_Data.csv')
        if self.equation==1:
            min_a=q.a.values[q.idxmin().sls]
            min_b=q.b.values[q.idxmin().sls]
            min_pt=q.peak_temp.values[q.idxmin().sls]
            min_st=q.start_time.values[q.idxmin().sls]
            params=(min_a, min_b, min_pt, min_st)
        #test_group=OptimizationGroup(120, )
        elif self.equation==2 or self.equation==3 or self.equation==4:
            min_a=q.a.values[q.idxmin().sls]
            min_b=q.b.values[q.idxmin().sls]
            min_c=q.c.values[q.idxmin().sls]
            min_pt=q.peak_temp.values[q.idxmin().sls]
            min_st=q.start_time.values[q.idxmin().sls]      
            params=(min_a, min_b, min_c, min_pt, min_st)
            
        # & Getting the best fit output
        self.group.run(params)     

        #return params
    
    def plot_fit(self):
        pass
        
    def get_initial_fit(self):
        self.get_opt_group()
        self.clean_data()
        q=pd.read_csv(self.fp+'/Clean_Opt_Data.csv')
        initial_params=q.iloc[0]
        if self.equation==1:
            params=q.iloc[0][1], q.iloc[0][2], q.iloc[0][3], q.iloc[0][4]
        elif self.equation==2 or self.equation==3 or self.equation==4:
            params=q.iloc[0][1], q.iloc[0][2], q.iloc[0][3], q.iloc[0][4], q.iloc[0][5]
        self.group.run(params)
        
    def get_optimization_plot(self):
        #Script for getting the final optimization plot
        if self.equation==1:
            fig=Figure(figsize=(15,10))
            ax_sls=fig.add_subplot(2,3,1)
            ax_a=fig.add_subplot(2,3,2)
            ax_b=fig.add_subplot(2,3,3)
            ax_start_time=fig.add_subplot(2,3,4)
            ax_peak_temp=fig.add_subplot(2,3,5)
            
            dataset=pd.read_csv(self.fp+'/Clean_Opt_Data.csv')
            ax_sls.scatter(dataset.iteration.values, dataset.sls.values, color='b', label='Real Values')
            ax_sls.plot(dataset.iteration.values, dataset.sls.rolling(30).mean().values, color='r', label='Moving avg')
            ax_a.scatter(dataset.iteration.values, dataset.a.values, s=10, color='b', label='Real Values')
            ax_a.plot(dataset.iteration.values, dataset.a.rolling(30).mean().values, color='r', label='Moving avg')
            ax_b.scatter(dataset.iteration.values, dataset.b.values, s=10, color='b', label='Real Values')
            ax_b.plot(dataset.iteration.values, dataset.b.rolling(30).mean().values, color='r', label='Moving avg')
            ax_start_time.scatter(dataset.iteration.values, dataset.start_time.values, s=10, color='b', label='Real Values')
            ax_start_time.plot(dataset.iteration.values, dataset.start_time.rolling(30).mean().values, color='r', label='Moving avg')
            ax_peak_temp.scatter(dataset.iteration.values, dataset.peak_temp.values, s=10, color='b', label='Real Values')
            ax_peak_temp.plot(dataset.iteration.values, dataset.peak_temp.rolling(30).mean().values, color='r', label='Moving avg')

            #Titles
            ax_sls.set_title('Sum of Least Squares', pad=13)
            ax_a.set_title('a', pad=13)
            ax_b.set_title('b', pad=13)
            ax_start_time.set_title('Start Time', pad=13)
            ax_peak_temp.set_title('Peak Temp', pad=13)

            ax_sls.legend()
            ax_a.legend()
            ax_b.legend()
            ax_start_time.legend()
            ax_peak_temp.legend()

            run_name=self.fp.split('/')[-1]
            fig.suptitle(run_name)
            fig.tight_layout()
            #plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.75)
            fig.savefig(self.fp+'/Optimization_Results.png')
        else:
            fig=Figure(figsize=(15,10))
            ax_sls=fig.add_subplot(2,3,1)
            ax_a=fig.add_subplot(2,3,2)
            ax_b=fig.add_subplot(2,3,3)
            ax_c=fig.add_subplot(2,3,4)
            ax_start_time=fig.add_subplot(2,3,5)
            ax_peak_temp=fig.add_subplot(2,3,6)
            
            dataset=pd.read_csv(self.fp+'/Clean_Opt_Data.csv')
            ax_sls.scatter(dataset.iteration.values, dataset.sls.values, color='b', label='Real Values')
            ax_sls.plot(dataset.iteration.values, dataset.sls.rolling(30).mean().values, color='r', label='Moving avg')
            ax_a.scatter(dataset.iteration.values, dataset.a.values, s=10, color='b', label='Real Values')
            ax_a.plot(dataset.iteration.values, dataset.a.rolling(30).mean().values, color='r', label='Moving avg')
            ax_b.scatter(dataset.iteration.values, dataset.b.values, s=10, color='b', label='Real Values')
            ax_b.plot(dataset.iteration.values, dataset.b.rolling(30).mean().values, color='r', label='Moving avg')
            ax_c.scatter(dataset.iteration.values, dataset.c.values, s=10, color='b', label='Real Values')
            ax_c.plot(dataset.iteration.values, dataset.c.rolling(30).mean().values, color='r', label='Moving avg')            
            ax_start_time.scatter(dataset.iteration.values, dataset.start_time.values, s=10, color='b', label='Real Values')
            ax_start_time.plot(dataset.iteration.values, dataset.start_time.rolling(30).mean().values, color='r', label='Moving avg')
            ax_peak_temp.scatter(dataset.iteration.values, dataset.peak_temp.values, s=10, color='b', label='Real Values')
            ax_peak_temp.plot(dataset.iteration.values, dataset.peak_temp.rolling(30).mean().values, color='r', label='Moving avg')
            
            #Titles
            ax_sls.set_title('Sum of Least Squares', pad=13)
            ax_a.set_title('a', pad=13)
            ax_b.set_title('b', pad=13)
            ax_c.set_title('c', pad=13)
            ax_start_time.set_title('Start Time', pad=13)
            ax_peak_temp.set_title('Peak Temp', pad=13)
            
            ax_sls.legend()
            ax_a.legend()
            ax_b.legend()
            ax_c.legend()
            ax_start_time.legend()
            ax_peak_temp.legend()
            
            run_name=self.fp.split('/')[-1]
            fig.suptitle(run_name)
            fig.tight_layout()
            #plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.75)
            fig.savefig(self.fp+'/Optimization_Results.png')
        
    def write_summary_file(self):
        pass
        
    def run_end_operations(self):
        self.clean_data()
        self.get_optimization_plot()
        
class EndOperationsBot:
    '''
    Class for automating and getting the end operations for optimizations that have been run
    '''
    def __init__(self):
        date=datetime.now()
        self.year, self.month, self.day = str(date.year), str(date.month), str(date.day)
        if len(self.month) == 1:
            self.month='0'+self.month
        if len(self.day) == 1:
            self.day = '0'+self.day
        self.general_opt_path = "../../winhome/Desktop/optimization_data/silly/../../winhome/Desktop/optimization_data/"
        self.days_data_path = "../../winhome/Desktop/optimization_data/silly/../../winhome/Desktop/optimization_data/"+str(self.year)+str(self.month)+str(self.day)
        self.yesterdays_data_path = "../../winhome/Desktop/optimization_data/silly/../../winhome/Desktop/optimization_data/"+str(self.year)+str(self.month)+str(int(self.day)-1)
        
        
    def get_todays_plots(self):
        if os.path.isdir(self.days_data_path):
            print(self.days_data_path)
            for file in os.listdir(self.days_data_path): #iterate through the files in the directory
                new_path = self.days_data_path+'/'+file
                if os.path.isdir(new_path):
                    has_end_plot=False #Assume this optimization folder does not have a plot at the end
                    has_opt_data=False #Assume there is not optimization data for this run
                    has_metadata = False #Assume no metadata file is present
                    for file1 in os.listdir(new_path):
                        if file1 == 'Optimization_Results.png':
                            has_end_plot = True
                        if file1 == 'optimization_output.csv':
                            has_opt_data = True
                        if file1 == 'metadata.txt':
                            #print(file1)
                            has_metadata = True
                    if has_end_plot == False and has_opt_data == True and has_metadata == True:
                        print(new_path)
                        post_opt_operator = PostOptOperations(new_path)
                        post_opt_operator.run_end_operations()
            #post_opt_operator = PostOptOperations(self.days_data_path)
            #post_opt_operator.run_end_operations()
        
    def get_yesterdays_plots(self):
        print('Running')
        if os.path.isdir(self.yesterdays_data_path):
            print(self.yesterdays_data_path)
            for file in os.listdir(self.yesterdays_data_path): #iterate through the files in the directory
                new_path = self.yesterdays_data_path+'/'+file
                if os.path.isdir(new_path):
                    has_end_plot=False #Assume this optimization folder does not have a plot at the end
                    has_opt_data=False #Assume there is not optimization data for this run
                    has_metadata = False #Assume no metadata file is present
                    for file1 in os.listdir(new_path):
                        if file1 == 'Optimization_Results.png':
                            has_end_plot = True
                        if file1 == 'optimization_output.csv':
                            has_opt_data = True
                        if file1 == 'metadata.txt':
                            #print(file1)
                            has_metadata = True
                    if has_end_plot == False and has_opt_data == True and has_metadata == True:
                        print(new_path)
                        post_opt_operator = PostOptOperations(new_path)
                        post_opt_operator.run_end_operations()
        
    def update_all_plots(self):
        for file in os.listdir(self.general_opt_path):
            #print(file)
            if len(file) == 8 and str(file).isnumeric()==True:
                #print('---------------')
                #print(file)
                #print('---------------')
                year=int(file[0:4])
                month=int(file[-4:-2])
                day=int(file[-2:])
                get_plots=False #Assume we don't want to get plots for this day
                if year == 2022:
                    if month > 4:
                        get_plots=True
                    elif month == 4:
                        if day >= 11:
                            get_plots=True
                else:
                    get_plots=True
                print('---------------')
                print(file, get_plots, year, month, day)
                print('---------------')
                    
                if get_plots == True:
                    d = os.path.join(self.general_opt_path, file)
                    if os.path.isdir(d): #If this is a folder, not a file
                        #print(d)

                        for opt_run in os.listdir(d): #Get each optimization run for the day
                            opt_fold_path = os.path.join(d, opt_run)
                            has_end_plot=False #Assume this optimization folder does not have a plot at the end
                            has_opt_data=False #Assume there is not optimization data for this run
                            has_metadata = False #Assume no metadata file is present
                            for file1 in os.listdir(opt_fold_path):
                                print(file1)
                                fname = os.path.join(d, file1)
                                #print(fname[-23:])
                                #print(file1)
                                if file1 == 'Optimization_Results.png':
                                    has_end_plot = True
                                if file1 == 'optimization_output.csv':
                                    has_opt_data = True
                                if file1 == 'metadata.txt':
                                    #print(file1)
                                    has_metadata = True

                        if has_opt_data == True and has_end_plot == False and has_metadata == True: #If there is optimization data but no end plot, add an end plot
                            o = PostOptOperations(opt_fold_path)
                            o.run_end_operations()