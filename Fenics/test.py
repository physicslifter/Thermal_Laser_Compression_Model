from fenics import *

def boundary (self, x ):
        value = x[0] < DOLFIN_EPS
        return value

class Test:
    def __init(self):
        pass
    
    def main(self, meshpoints, a=0.03, b=30, start_time=1.8*10**-8, k_1=100, init_temp=2000, tf=30*10**-9, num_steps=60, Fe_length=1.07*10**-6, MgO_length=2*10**-6, peak_temp=30000):
        self.meshpoints, self.a, self.b, self.start_time, self.num_steps, self.Fe_length, self.MgO_length, self.peak_temp = meshpoints, a, b, start_time, num_steps, Fe_length, MgO_length, peak_temp
        self.dt=tf/num_steps
        self.my_mesh = IntervalMesh(meshpoints, 0, Fe_length+MgO_length)
        self.V = FunctionSpace(self.my_mesh, 'P', 1)
        self.u_D = Expression(str(peak_temp),degree=1)
        self.bc = DirichletBC(self.V, self.u_D, boundary)
        self.u_0 = Expression(str(init_temp), degree=1)
        self.u_n = interpolate(self.u_0, self.V)
        self.u = Function(self.V)
        self.v = TestFunction(self.V)
        self.kappa = Expression('x[0] <= l + tol ? b+a*u : k_1', degree=1, tol=DOLFIN_EPS, k_1=k_1, u=self.u, a=a, b=b, l=Fe_length)
        
Test()