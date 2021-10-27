#import necessary packages & our GUI function...
import matplotlib.animation as animation
import pandas as pd 
import time
import GUI_functions as gf 
from matplotlib import pyplot as plt

#set up plot
fig=plt.figure()
ax1=fig.add_subplot(2,3,1)
ax2=fig.add_subplot(2,3,2)
ax3=fig.add_subplot(2,3,3)
ax4=fig.add_subplot(2,3,4)
ax5=fig.add_subplot(2,3,5)
#ax5=fig.add_subplot(2,3,1)

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


ani=animation.FuncAnimation(fig,animate, interval=1000)
plt.show()
