#functions for testing the chi^2 interpolation function to ensure that it works
import numpy as np
from run_1d_fem_v3 import chi_sq_interp
import matplotlib.pyplot as plt
import numpy as np
import datetime


def sine(x,a,b,c,d):
    return a*np.sin(b*x+c)+d

#function for returning the chi^2 value between the modeled and real sine curve
def chi_2(params):

    a=params[0]
    b=params[1]
    c=params[2]
    d=params[3]

    #space the x data for the model
    x0=np.arange(1,10)

    #Get model output
    model_data=sine(x0,a,b,c,d)

    #Now the "actual" data, for the function f(x)=sin(x)
    x1=np.arange(1,10,0.1)
    y1=np.sin(x1)

    #Now let's add some noise to the data
    y1noise=np.random.normal(0,0.1,y1.shape)
    y1=y1+y1noise

    #Now time to assess the chi^2 between the model value and the output value
    chi_2=chi_sq_interp(x0,model_data,x1,y1)
    print(chi_2)

    #save a plot of the data
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    ax1.scatter(x1,y1,s=10,marker='o',label='data')
    ax1.plot(x0,model_data, label="model")
    plt.legend()

    #get timestamp for saving name and save plot to file
    time_stamp=str(np.datetime64(datetime.datetime.now()).item()).replace(' ','').replace(':','').replace('.','')
    plt.savefig('chi_sq_test_plots/test2_'+time_stamp+'.png')
    


    #save to output file
    with open('chi_sq_test_output.txt', 'a+') as file_object:
        file_object.write('\n')#newlinw
        file_object.write(str(chi_2))

    #Now wait for some time to ensure the next iteration has a unique filename
    count=0
    while(count<200):
        count=count+1

    #Now, return the value
    return chi_2


