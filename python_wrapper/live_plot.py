#import necessary packages & our GUI function...
import matplotlib.animation as animation
import pandas as pd 
import time
import GUI_functions as gf 
from matplotlib import pyplot as plt

#set up plot
fig=plt.figure()
ax1=fig.add_subplot(1,1,1)

def animate(i):
    pullData = open("optimization_data.csv","r").read()
    dataArray = pullData.split('\n')

    #defining arrays for each relevant value...
    itar = [] #iterations
    slsar = [] #sum of least squares
    ptar = [] #peak temp
    aar = [] #a
    bar = [] #b
    tsar = [] #time shift
    diffar = [] #diffusivity

    for eachLine in dataArray:
        if len(eachLine)>1:
            it, sls, pt, a, b, ts, diff = eachLine.split(',')
            itar.append(float(it))
            slsar.append(float(sls))
            ptar.append(float(pt))
            aar.append(float(a))
            bar.append(float(b))
            tsar.append(float(ts))
            diffar.append(float(diff))

    ax1.clear()
    ax1.plot(itar,slsar)


ani=animation.FuncAnimation(fig,animate, interval=1000)
plt.show()
