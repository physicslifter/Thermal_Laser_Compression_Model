from json import tool
from matplotlib import pyplot as plt
from fenics import *
import numpy as np

'''
File for getting a simple FEM running that is like a single FEM of our model in Matlab
'''

#Defining Boundary Condition square wave
def boundary ( x ):
    value = x[0] < DOLFIN_EPS
    return value

def main(tf, num_steps, a, b, peak_temp, init_temp, k_1, rho, c, Fe_length, MgO_length, start_time):
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
            
    del dt, meshpoints, my_mesh, V, u_D, bc, u_0, u_n, u, v, kappa,inc, loc, points, num_loops, t, F 
        
    times=np.asarray(times)+start_time #+start_time
    
    return times, time_line

if __name__ == "__main__":
    main()