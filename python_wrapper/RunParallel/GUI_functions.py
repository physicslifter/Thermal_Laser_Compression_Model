#Functions to be utilized by the live data interface
import pandas as pd
import matplotlib.animation as animation
import time

#Reading from the CSV file
def read_file():
    pd.read_csv('optimization_data.csv')

def get_file_length():
    return len(read_file())

def animate(axis):
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

    axis.clear()
    axis.plot(itar,slsar)

