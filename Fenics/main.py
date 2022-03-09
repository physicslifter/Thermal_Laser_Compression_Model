'''
Main file for Thermal Laser Compression Model
Pat LaChapelle
Created: February 10, 2022
Last updated: March 9, 2022
'''
import fenics
import numpy as np
import matplotlib
matplotlib.use('WebAgg')
from matplotlib import pyplot as plt
import json

default_parameters=[0.03, 30, 1.8*10**-8, 30000] #Define some default parameters for passing to function files when convenient

def load_json(fname): #function for loading in json files
    with open(fname) as json_file: #assign opened file to object json_file
        data=json.load(json_file) #load data
    return data

def boundary ( x ): #Boundary function used for finite element model
    value = x[0] < fenics.DOLFIN_EPS #Dirichlet boundary condition 
    return value

class SingleFaceModel:
    def __init__(self, tf=30*10**-9, num_steps=60, init_temp=2000, k_1=100, rho=12000, c=450, MgO_length=2*10**-6, Fe_length=1.07*10**-6, meshpoints=120):
        self.tf, self.num_steps, self.init_temp, self.k_1, self.rho, self.c, self.MgO_length, self.Fe_length, self.meshpoints = tf, num_steps, init_temp, k_1, rho, c, MgO_length, Fe_length, meshpoints #assign values to self so we can reference them later as needed
        self.dt=self.tf/self.num_steps #time step (dt) is total time divided by the number of time steps
    
    def run(self, parameters=default_parameters): # Function for running the single face model
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
        self.kappa = fenics.Expression('x[0] <= l + tol ? b+a*u : k_1', degree=1, tol=fenics.DOLFIN_EPS, k_1=self.k_1, u=self.u, a=a, b=b, l=self.Fe_length) # ||
        inc = 0.01*10**-6 #Define some small distance to measure off of the boundary, so that we are not pulling the data from exactly on the boundary
        loc = self.Fe_length - inc # measuring location
        self.times = np.linspace(0,self.tf, self.num_steps) # Get list of times for the model. This will be the independent variable of the model output
        self.times = (np.asarray(self.times)+start_time).tolist() # shift times to account for the start time of the model
        self.time_line = [] # Empty list for getting the model output
        t=0 #start time @ 0 (will be shifted by start time from line 46)
        self.F = self.u*self.v*fenics.dx + self.kappa/(self.rho*self.c)*self.dt*fenics.dot(fenics.grad(self.u), fenics.grad(self.v))*fenics.dx - (self.u_n)*self.v*fenics.dx # Defining variational problem
        for n in range(self.num_steps): # Loop through model to get output
            t+=self.dt # time step through model
            self.u_D.t=t # update current time
            fenics.solve(self.F==0, self.u, self.bc) # Solve the variational problem for the current time to get the Temperature in terms of the position
            self.time_line.append(self.u(loc)) # Get the Temperature @ the measurement location & add it to the output data
            self.u_n.assign(self.u) # update previous solution
            
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
    
    def __init__(self, meshpoints, geometry_dict:Geometry): # Initialize an object of class FiniteElementModel, defining the number of meshpoints desired
        vals=load_json('default_values_dict.json') # load in default values from the dictionary of default values
        self.tf, self.num_steps, self.MgO_length, self.init_temp, self.k_1, self.rho, self.c, self.meshpoints = vals['tf'], vals['num_steps'], vals['MgO_length'], vals['init_temp'], vals['k_1'], vals['rho'], vals['c'], meshpoints #assign values
        self.dt=self.tf/self.num_steps # time step
        self.faces=geometry_dict.geo.keys() # Get the names of the faces that we need to build finite element models for
        face_count=1 # first face is always 1, so write face_count to 1 before iterating through the faces
        self.geometry=geometry_dict.geo # assign the geometry dictionary to a self object
        models={} # Write new blank dictionary named models for storing the model that has been created for each face
        for face in self.faces: # iterate through the faces
            key='face'+str(face_count) # Get the key for the face in the geometry dictionary
            Fe_length=self.geometry[face] # Get the finite element length for the face from the geometry dictionary
            models[key]=(SingleFaceModel(self.tf, self.num_steps, self.init_temp, self.k_1, self.rho, self.c, self.MgO_length, Fe_length, self.meshpoints)) # Write the model into the model dictionary as the values associated with the key for the face
            face_count=face_count+1 # increment face_count by 1
        self.face_models=models # Now assign models dictionary as a property of self so that we can reference it later 
        
    def run_model(self, parameters): # Method for running the model that has been defined by the initialization of the object of class FiniteElementModel
        results={} # Create empty dictionary for storing the results of the model run
        for face in self.face_models.keys(): # Iterate through the faces
            model=self.face_models[face] # Get the model that is associated with this face
            model.run(parameters) #run the model
            results[face]=model.times, model.time_line #store the model results
            
    def plot_model(self): # Define function for plotting the current output of the model
        plt.clf() # Clear plot
        for face in self.face_models.keys(): # iterate through the faces
            model=self.face_models[face] # Get the model for the face
            plt.plot(model.times, model.time_line) # Plot the model output data from the model
        plt.show() # Display the plot
        
class OptimizationGroup: # A class for an optimization group. An optimization group can contain several runs and several faces, which are specified in the run_dictionary
    
    def __init__(self, meshpoints, run_dictionary): # Initialize an OptimizationGroup object. run_dictionary should have keys for each runID (ie 's88776'). Associated value for each runID key is a list of faces (ie [1,2,3]). Also takes number of meshpoints to use for the finite element model of each run and face
        self.run_dictionary=run_dictionary
        self.models={}
        for run in run_dictionary.keys(): # iterate through the different runs
            run_geometry=Geometry(run, run_dictionary[run]) # Get geometry for the run and faces
            f=FiniteElementModel(meshpoints, run_geometry) # Create a Finite element model for each run in run_dictionary
            self.models[run]=f
            
    def run(self, params): #Write function for running the model for each shot
        for run in self.models.keys(): # iterate through each shot
            model=self.models[run] # Get the model for the run
            model.run_model(params) # run the model for the given parameters