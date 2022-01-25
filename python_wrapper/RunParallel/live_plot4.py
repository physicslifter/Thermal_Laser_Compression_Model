#import necessary packages & our GUI function...
import matplotlib.animation as animation
import pandas as pd 
import time
import GUI_functions as gf 
from matplotlib import pyplot as plt
from scipy import io

#set up plot
fig1=plt.figure()
fig1.suptitle('Run 4')
ax1=fig1.add_subplot(3,3,1) #sls
ax2=fig1.add_subplot(3,3,2) #peak_temp
ax3=fig1.add_subplot(3,3,3) #a
ax4=fig1.add_subplot(3,3,4) #b
ax5=fig1.add_subplot(3,3,5) #start_time
ax6=fig1.add_subplot(3,3,6) #s88773
ax7=fig1.add_subplot(3,3,7) #s88776
ax8=fig1.add_subplot(3,3,8) #s88780
ax9=fig1.add_subplot(3,3,9) #s86483

def animate(i):
    pullData = open("data4.csv","r").read()
    dataArray = pullData.split('\n')

    #========================================
    #Parameters
    #========================================
    #defining arrays for each relevant value...
    itar = [] #iterations
    slsar = [] #sum of least squares
    ptar = [] #peak temp
    aar = [] #a
    bar = [] #b
    tsar = [] #time shift
    #diffar = [] #diffusivity

    for eachLine in dataArray:
        if len(eachLine)>1:
            it, sls, pt, a, b, ts = eachLine.split(',')
            itar.append(float(it))
            slsar.append(float(sls))
            ptar.append(float(pt))
            aar.append(float(a))
            bar.append(float(b))
            tsar.append(float(ts))
            #diffar.append(float(diff))

    #===========================================
    #Fits
    #===========================================
    d=io.loadmat('combo4.mat')

    x0=d['tlist'][0] #tlist model data for all runs @ all faces

    #run s88773
    x1_1_73=d['t_data1_73'] #time data
    x1_2_73=d['t_data2_73'] #time data
    x1_3_73=d['t_data3_73'] #time data
    y0_1_73=d['T11'][0] #face 1 model output
    y0_2_73=d['T12'][0] # face 2 model output
    y0_3_73=d['T13'][0] # face 3 model output
    y1_1_73=d['temp1_73'] # face 1 experimental data
    y1_2_73=d['temp2_73'] # face 2 experimental data
    y1_3_73=d['temp3_73'] # face 3 experimental data

    #run s88776
    x1_1_76=d['t_data1_76'] #time data
    x1_2_76=d['t_data2_76'] #time data
    x1_3_76=d['t_data3_76'] #time data
    y0_1_76=d['T21'][0] #face 1 model output
    y0_2_76=d['T22'][0] # face 2 model output
    y0_3_76=d['T23'][0] # face 3 model output
    y1_1_76=d['temp1_76'] # face 1 experimental data
    y1_2_76=d['temp2_76'] # face 2 experimental data
    y1_3_76=d['temp3_76'] # face 3 experimental data

    #run s88780
    x1_1_80=d['t_data1_80'] #time data
    x1_2_80=d['t_data2_80'] #time data
    x1_3_80=d['t_data3_80'] #time data
    y0_1_80=d['T31'][0] #face 1 model output
    y0_2_80=d['T32'][0] # face 2 model output
    y0_3_80=d['T33'][0] # face 3 model output
    y1_1_80=d['temp1_80'] # face 1 experimental data
    y1_2_80=d['temp2_80'] # face 2 experimental data
    y1_3_80=d['temp3_80'] # face 3 experimental data

    #run s86483
    x1_1_83=d['t_data1_83'] #time data
    x1_2_83=d['t_data2_83'] #time data
    x1_3_83=d['t_data3_83'] #time data
    y0_1_83=d['T41'][0] #face 1 model output
    y0_2_83=d['T42'][0] # face 2 model output
    y0_3_83=d['T43'][0] # face 3 model output
    y1_1_83=d['temp1_83'] # face 1 experimental data
    y1_2_83=d['temp2_83'] # face 2 experimental data
    y1_3_83=d['temp3_83'] # face 3 experimental data


    #==========================================
    #plotting
    #==========================================
    ax1.clear()
    ax2.clear()
    ax3.clear()
    ax4.clear()
    ax5.clear()
    ax6.clear()
    ax6.clear()
    ax7.clear()
    ax8.clear()
    ax9.clear()
    ax1.set_title('sls')
    ax2.set_title('peak_temp')
    ax3.set_title('a')
    ax4.set_title('b')
    ax5.set_title('start_time')
    ax6.set_title('s88773')
    ax7.set_title('s88776')
    ax8.set_title('s88780')
    ax9.set_title('s86483')

    #parameters
    ax1.plot(itar,slsar)
    ax2.plot(itar,ptar)
    ax3.plot(itar,aar)
    ax4.plot(itar,bar)
    ax5.plot(itar,tsar)

    #s88773
    ax6.plot(x0, y0_1_73, color='mediumorchid')
    ax6.scatter(x1_1_73,y1_1_73, color='blue')
    ax6.plot(x0, y0_2_73, color='green')
    ax6.scatter(x1_2_73, y1_2_73, color='darkorange')
    #ax6.scatter(x1_3_73, y1_3_73, color='gold')
    #ax6.plot(x0, y0_3_73, color='cyan')

    #s88776
    ax7.plot(x0, y0_1_76, color='mediumorchid')
    ax7.scatter(x1_1_76,y1_1_76, color='blue')
    #ax7.plot(x0, y0_2_76, color='green')
    #ax7.scatter(x1_2_76, y1_2_76, color='darkorange')
    #ax7.scatter(x1_3_76, y1_3_76, color='gold')
    #ax7.plot(x0, y0_3_76, color='cyan')

    #s88780
    ax8.plot(x0, y0_1_80, color='mediumorchid')
    ax8.scatter(x1_1_80,y1_1_80, color='blue')
    ax8.plot(x0, y0_2_80, color='green')
    ax8.scatter(x1_2_80, y1_2_80, color='darkorange')
    ax8.scatter(x1_3_80, y1_3_80, color='gold')
    ax8.plot(x0, y0_3_80, color='cyan')

    #s88773
    ax9.plot(x0, y0_1_83, color='mediumorchid')
    ax9.scatter(x1_1_83,y1_1_83, color='blue')
    ax9.plot(x0, y0_2_83, color='green')
    ax9.scatter(x1_2_83, y1_2_83, color='darkorange')
    ax9.scatter(x1_3_83, y1_3_83, color='gold')
    ax9.plot(x0, y0_3_83, color='cyan')


ani=animation.FuncAnimation(fig1,animate, interval=20000)
plt.show()