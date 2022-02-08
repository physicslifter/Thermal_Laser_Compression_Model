import matplotlib
matplotlib.use('WebAgg')
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib import pyplot as plt
from scipy import io
from plot import plot
import json

def load_json(fname):
    with open(fname) as json_file:
        data=json.load(json_file)
    return data

real_data_filepath='data_dict.json'
model_data_filepath='model_results.json'

parameters=('a', 'b', 'peak_temp', 'start_time')

fig=plt.figure()
ax1=fig.add_subplot(2,3,1) #sls
ax2=fig.add_subplot(2,3,2) #peak_temp
ax3=fig.add_subplot(2,3,3) #a
ax4=fig.add_subplot(2,3,4) #b
ax5=fig.add_subplot(2,3,5) #start_time

def animate(i):
    
    #========================================
    #Parameters
    #========================================
    pullData = open("optimization_output.csv","r").read()
    dataArray = pullData.split('\n')
    
    #defining arrays for each relevant value...
    itar = [] #iterations
    slsar = [] #sum of least squares
    ptar = [] #peak temp
    aar = [] #a
    bar = [] #b
    tsar = [] #time shift
    
    for eachLine in dataArray:
        if len(eachLine)>1:
            it, sls, a, b, pt, ts = eachLine.split(',')
            itar.append(float(it))
            slsar.append(float(sls))
            ptar.append(float(pt))
            aar.append(float(a))
            bar.append(float(b))
            tsar.append(float(ts))
            
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
    
    ax1.plot(itar,slsar)
    ax2.plot(itar,ptar)
    ax3.plot(itar,aar)
    ax4.plot(itar,bar)
    ax5.plot(itar,tsar)
    
ani=animation.FuncAnimation(fig,animate, interval=20000)
plt.show()