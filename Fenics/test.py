from xml.sax.handler import DTDHandler
from fenics import *
import numpy as np
import matplotlib
matplotlib.use('WebAgg')
from matplotlib import pyplot as plt
import json

default_parameters=[0.03, 30, 1.8*10**-8, 30000]

def load_json(fname):
    with open(fname) as json_file:
        data=json.load(json_file)   
    return data

def boundary ( x ):
    value = x[0] < DOLFIN_EPS
    return value

class SingleFaceModel:
    def __init__(self, tf=30*10**-9, num_steps=60, init_temp=2000, k_1=100, rho=12000, c=450, MgO_length=2*10**-6, Fe_length=1.07*10**-6, meshpoints=120):
        self.tf, self.num_steps, self.init_temp, self.k_1, self.rho, self.c, self.MgO_length, self.Fe_length, self.meshpoints = tf, num_steps, init_temp, k_1, rho, c, MgO_length, Fe_length, meshpoints
        self.dt=self.tf/self.num_steps
    
    def run(self, parameters=default_parameters):
        #self.meshpoints, self.a, self.b, self.start_time, self.num_steps, self.Fe_length, self.MgO_length, self.peak_temp, self.rho, self.c, self.k_1, self.tf, self.init_temp = meshpoints, a, b, start_time, num_steps, Fe_length, MgO_length, peak_temp, rho, c, k_1, tf, init_temp
        a, b, start_time, peak_temp = parameters
        self.a, self.b, self.start_time, self.peak_temp=a,b,start_time,peak_temp
        self.dt=self.tf/self.num_steps
        self.my_mesh = IntervalMesh(self.meshpoints, 0, self.Fe_length+self.MgO_length)
        self.V = FunctionSpace(self.my_mesh, 'P', 1)
        self.u_D = Expression(str(self.peak_temp),degree=1)
        self.bc = DirichletBC(self.V, self.u_D, boundary)
        self.u_0 = Expression(str(self.init_temp), degree=1)
        self.u_n = interpolate(self.u_0, self.V)
        self.u = Function(self.V)
        self.v = TestFunction(self.V)
        self.kappa = Expression('x[0] <= l + tol ? b+a*u : k_1', degree=1, tol=DOLFIN_EPS, k_1=self.k_1, u=self.u, a=a, b=b, l=self.Fe_length)
        
        inc = 0.01*10**-6
        loc = self.Fe_length - inc
        points = np.linspace(0,self.Fe_length+self.MgO_length,self.meshpoints)
        self.times = np.linspace(0,self.tf, self.num_steps)
        self.time_line = []
        t=0
        
        self.u_n=interpolate(self.u_0, self.V)
        
        self.F = self.u*self.v*dx + self.kappa/(self.rho*self.c)*self.dt*dot(grad(self.u), grad(self.v))*dx - (self.u_n)*self.v*dx
        self.time=[]
        for n in range(self.num_steps):
            t+=self.dt 
            self.u_D.t=t
            solve(self.F==0, self.u, self.bc)
            self.time_line.append(self.u(loc))
            self.u_n.assign(self.u)
            
    def plot(self):
        plt.clf()
        plt.plot(self.times, self.time_line)
        plt.show()
        
class Geometry:
    def __init__(self, runID:str, facelist, dictionary_filepath='geometry_dict.json'):
        with open(dictionary_filepath) as json_file:
            data=json.load(json_file)
            data=data[runID]
        new_dict={}
        for face in data:
            for id in facelist:
                name='face'+str(id)
                if name==face:
                    new_dict[face]=data[face]
        self.geo=new_dict

class FiniteElementModel:
    
    def __init__(self, meshpoints, geometry_dict:Geometry):
        vals=load_json('default_values_dict.json')
        self.tf=vals['tf']
        self.num_steps=vals['num_steps']
        self.MgO_length=vals['MgO_length']
        self.init_temp=vals['init_temp']
        self.k_1=vals['k_1']
        self.rho=vals['rho']
        self.c=vals['c']
        self.meshpoints=meshpoints
        self.dt=self.tf/self.num_steps
        self.faces=geometry_dict.geo.keys()
        face_count=1
        self.geometry=geometry_dict.geo
        models={}
        for face in self.faces:
            key='face'+str(face_count)
            Fe_length=self.geometry[key]
            models[key]=(SingleFaceModel(self.tf, self.num_steps, self.init_temp, self.k_1, self.rho, self.c, self.MgO_length, Fe_length, self.meshpoints))
            #setattr(self,key,SingleFaceModel(self.meshpoints, self.tf, self.num_steps, Fe_length, self.MgO_length))
            face_count=face_count+1
        self.face_models=models