#A script for running a single iteration of the FEM model from the python wrapper, and printing the results
#===========================================================

#importing packages
from scipy import io 
from matlab import engine
from matplotlib import pyplot as plt

#import functions from runFEM2.py as f
import runFEM2 as f

#==========================================================
#Defining values
#CHANGE THIS:
#values have to be defines as float (add a decimal to round numbers; ex: ".0")

#Boundary conditions
peak_temp=28000.0 #peak temperature parameter from the boundary condition file
start_time=1.8*10**(-8) #start time parameter from boundary condition file

#thermal conductivity parameters (for equation k=a*x+b)
a=0.035
b=30.0

#=========================================================
#Running the function
mydata=f.run_model(peak_temp,start_time, a, b)

#printing the results
f.plot_results(mydata)
plt.show()
