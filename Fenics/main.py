import matplotlib
matplotlib.use('WebAgg')
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from json import tool
from xmlrpc.server import DocXMLRPCRequestHandler
from fenics import *
import numpy as np
from sympy import *
import json

def load_json(fname):
    with open(fname) as json_file:
        data=json.load(json_file)   
    return data

def boundary(x):
    value=x[0] < DOLFIN_EPS
    return value

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
        
        
class SingleFaceModel:
    def __init__(self, tf, num_steps, init_temp, k_1, rho, c, MgO_length, Fe_length, meshpoints):
        self.tf, self.num_steps, self.init_temp, self.k_1, self.rho, self.c, self.MgO_length, self.Fe_length = tf, num_steps, init_temp, k_1, rho, c, MgO_length, Fe_length
        self.dt=self.tf/self.num_steps
        self.meshpoints=meshpoints

    def run(self, parameters):
        a, b, start_time, peak_temp = parameters
        self.dt=self.tf/self.num_steps
        self.my_mesh = IntervalMesh(self.meshpoints, 0, self.Fe_length+self.MgO_length)
        self.V = FunctionSpace(self.my_mesh, 'P', 1)
        self.u_D = Expression(str(peak_temp),degree=1)
        self.bc = DirichletBC(self.V, self.u_D, boundary)
        self.u_0 = Expression(str(self.init_temp), degree=1)
        self.u_n = interpolate(self.u_0, self.V)
        
        '''my_mesh = IntervalMesh(self.meshpoints, 0, self.Fe_length+self.MgO_length)
        V = FunctionSpace(my_mesh, 'P', 1)
        u_D = Expression(str(peak_temp),degree=1)
        bc = DirichletBC(V, u_D, SingleFaceModel.boundary)
        u_0 = Expression(str(self.init_temp), degree=1)
        u_n = interpolate(u_0, V)
        u = Function(V)
        v = TestFunction(V)
        kappa = Expression('x[0] <= l + tol ? b+a*u : k_1', degree=1, tol=DOLFIN_EPS, k_1=self.k_1, u=u, a=a, b=b, l=self.Fe_length)
        
        inc = 0.01*10**-6
        loc = Fe_length - inc
        points = np.linspace(0,Fe_length+MgO_length,meshpoints)
        times = np.linspace(0,tf, num_steps)
        time_line = []
        num_loops = 1
        t=0
        F=self.u*self.v*dx + (self.kappa/(self.rho*self.c))*self.dt*dot(grad(self.u), grad(self.v))*dx - (self.u_n)*self.v*DocXMLRPCRequestHandler
        for n in range(num_steps):
            t += dt
            self.u_D.t=t
            solve(F==0, self.u, self.bc)
            time_line.append(self.u(loc))
            self.u_n.assign(u)
        self.times=(np.asarray(times)+start_time).tolist() #+start_time
        self.time_line=time_line'''
        
    def plot_run(self):
        plt.plot(self.times, self.time_line)
        
    