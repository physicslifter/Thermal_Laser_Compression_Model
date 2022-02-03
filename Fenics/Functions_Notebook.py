#%%
from json import tool
from matplotlib import pyplot as plt
from fenics import *
import numpy as np
import json

'''
File for getting a simple FEM running that is like a single FEM of our model in Matlab
'''

#Defining Boundary Condition square wave
def boundary ( x ):
    value = x[0] < DOLFIN_EPS
    return value

def run_simple_fem(tf, num_steps, a, b, peak_temp, init_temp, k_1, rho, c, Fe_length, MgO_length):
    # | T  - Final time
    # | num_steps - number of time steps
    # | dt - time step size
    
    dt = tf/num_steps
    meshpoints=120
    
    my_mesh = IntervalMesh(meshpoints, 0, Fe_length+MgO_length)
    V = FunctionSpace(my_mesh, 'P', 1)
    
    #Define Boundary Condition
    u_D = Expression(str(peak_temp),degree=1)

    bc = DirichletBC(V, u_D, boundary)
    
    # Define initial value
    u_0 = Expression(str(init_temp), degree=1)
    u_n = interpolate(u_0, V)
    #u_n = project(u_D, V)
    
    # Define variational problem
    u = Function(V)
    v = TestFunction(V)
    kappa = Expression('x[0] <= l + tol ? b+a*u : k_1', degree=1, tol=DOLFIN_EPS, k_1=k_1, u=u, a=a, b=b, l=Fe_length)
    
    # Time-stepping
    inc = 0.01*10**-6
    loc = Fe_length - inc
    points = np.linspace(0,Fe_length+MgO_length,meshpoints)
    times = np.linspace(0,tf, num_steps)
    time_line = []
    num_loops = 1
    t=0
    
    for loop in range(num_loops):
        t=0
        u_n = interpolate(u_0, V)
        F = u*v*dx + kappa/(rho*c)*dt*dot(grad(u), grad(v))*dx - (u_n)*v*dx
        time_line = []
        for n in range(num_steps):

            # Update current time
            t += dt
            u_D.t = t
            # Compute solution
            solve(F==0, u, bc)

            # Plot solution
            #u_line = [u(point) for point in points]
            time_line.append(u(loc))

            # Update previous solution
            u_n.assign(u)
        
    return times, time_line

#standard values for each
tf = 30*10**-9     # final time
num_steps = 60     # number of time steps
dt = tf / num_steps # time step size
init_temp = 2000
peak_temp = 30000
a = 0.03
b = 30
k_1 = 100
rho=12000
c=450
Fe_length=2.126*10**-6
MgO_length=2*10**-6

test_run=run_simple_fem(tf=tf, num_steps=num_steps, a=a, b=b, peak_temp=peak_temp, init_temp=init_temp, k_1=k_1, rho=rho, c=c, Fe_length=Fe_length, MgO_length=MgO_length)

plt.plot(test_run[0], test_run[1])
plt.title('Test run of simple FEM function')
plt.show()
# %%

#Constants

#Fe
Fe_density=12800 #kg/m^3 - determined from Smith et al. 2018 isentrope (previous at 11000) 190 GPa for MgO
Fe_specific_heat=500 #J/kg K
Fe_inc = 0.1 #amount measuring into window

#MgO
MgO_density=5800 #kg/m^3 (4200) from Jin et al. 2010
MgO_specific_heat=800 #J/kg K (800)
MgO_inc = 0.1 #measurement into iron

#%%
#Geometry
5

#%%

#Function for opening our JSON file
def load_json(fname):
    with open(fname) as json_file:
        data=json.load(json_file)
        
    return data

#Now we can define simple functions for running each of the models independently
def run_model(run_IDs, num_steps, a, b, peak_temp):
    num_runs=len(run_IDs)
    geometry_dict=load_json('geometry_dict.json')
    values_dict=load_json('default_values_dict.json')
    results={}
    for x in range(num_runs):
        runID=run_IDs[x] # get the xth name of the run specified in run_IDs input
        solution={
            'face1':run_simple_fem(tf=values_dict['tf'], num_steps=num_steps,a=a, b=b, peak_temp=peak_temp, init_temp=values_dict['init_temp'], k_1=values_dict['k_1'], rho=values_dict['rho'], c=values_dict['c'], Fe_length=geometry_dict[runID]['face1'] , MgO_length=values_dict['MgO_length']),
            'face2':run_simple_fem(tf=values_dict['tf'], num_steps=num_steps,a=a, b=b, peak_temp=peak_temp, init_temp=values_dict['init_temp'], k_1=values_dict['k_1'], rho=values_dict['rho'], c=values_dict['c'], Fe_length=geometry_dict[runID]['face2'] , MgO_length=values_dict['MgO_length']),
            'face3':run_simple_fem(tf=values_dict['tf'], num_steps=num_steps,a=a, b=b, peak_temp=peak_temp, init_temp=values_dict['init_temp'], k_1=values_dict['k_1'], rho=values_dict['rho'], c=values_dict['c'], Fe_length=geometry_dict[runID]['face3'] , MgO_length=values_dict['MgO_length'])
        }
        results[runID]=solution
    return results
#%%
#Now time to show the function works
test_run=run_model(['s88773'], 60, 0.03, 30, 30000)
# %%
runID='s88773'
geometry_dict=load_json('geometry_dict.json')
values_dict=load_json('default_values_dict.json')
test1=run_simple_fem(tf=values_dict['tf'], num_steps=60,a=0.03, b=30, peak_temp=30000, init_temp=values_dict['init_temp'], k_1=values_dict['k_1'], rho=values_dict['rho'], c=values_dict['c'], Fe_length=geometry_dict[runID]['face1'] , MgO_length=values_dict['MgO_length'])

#%%
run_simple_fem(tf=30*10**-9, num_steps=60, a=0.03, b=30, peak_temp=30000, init_temp=2000, k_1=100, rho=12000, c=450, Fe_length=convert_geometry(1.07), MgO_length=2)
# %%
