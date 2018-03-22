
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
import os
import math
import sys

############################################################################################
### For a single run, plot:         ########################################################
### Order parameter over the run    ########################################################
### System orientation over the run ########################################################
### Mean-squared displacement       ########################################################
### Log(1-MSD)                      ########################################################
### Density fluctuation scaling     ########################################################
### Autocorrelation function        ########################################################
### Correlation function            ########################################################
### Correlation function log-log    ########################################################
### Correlation function log-lin    ########################################################
### Pair correlation function       ########################################################
### Velocity distribution           ########################################################
############################################################################################

runID = str(sys.argv[1])
fullRunID = str(sys.argv[2])
NDIM = int(sys.argv[4])

### Print a maximum of 1000 points per line to avoid the graphs being too crowded

npoints = 1000
inv = int(int(sys.argv[3])/npoints)
if inv==0:
    inv = 1

#directory = '/Users/Daniel1/Desktop/ActiveMatterResearch/clusterOutput/'.'$2/'
#directory = '/Users/Daniel1/Desktop/ActiveMatterResearch/clusterOutput/'+fullRunID+'/'+fullRunID+'/'+runID
#directory = '/home/dmccusker/remote/jamming-dynamics/output/'.'$2/'
directory = "/Users/Daniel1/Desktop/ActiveMatterResearch/jamming-dynamics/output/"+fullRunID+"/"+runID
#directory = "/Users/Daniel1/Desktop/ActiveMatterResearch/jamming-dynamics/output/3Dtest/test4"

titles = [ "Order parameter",
           "Orientation",
           "Mean-squared displacement",
           "Log of mean-squared displacement",
           "Density fluctuations",
           "Autocorrelation function",
           "Correlation function",
           "Pair correlation function",
           "Velocity distribution"]

axislabels = [ [ "Steps", "Phi" ],
               [ "x", "y", "z", "t" ],
               [ "steps", "MSD" ],
               [ "steps", "Log(1-MSD)" ],
               [ "<V>", "rms fluctuation" ],
               [ "lag time steps", "C" ],
               [ "distance", "C" ],
               [ "distance", "C" ],
               [ "speed", "frequency" ] ]

inputs = [ "order.dat",
           "orientation.dat",
           "MSD.dat",
           "MSD.dat",
           "fluct.dat",
           "autoCorr.dat",
           "corr.dat",
           "pairCorr.dat",
           "velDist.dat"]

outputs = [ "order.eps",
            "orientation.eps",
            "MSD.eps",
            "MSDLog.eps",
            "fluct.eps",
            "autoCorr.eps",
            "corr.eps",
            "pairCorr.eps",
            "velDist.eps"]

def linear(x,m,b):
    return m*x+b
def exponential(x,m,b,c):
    return b*np.power( np.e, -1.0*x/m ) + c
def powerlaw(x,m,b,c):
    return b*x**m + c

os.system("rm "+directory+"/dat/summary3.dat")
os.system("rm "+directory+"/dat/summary4.dat")
summary3 = open(directory+"/dat/summary3.dat","a")
summary4 = open(directory+"/dat/summary4.dat","a")

nplots = len(inputs)
for i in range(0,nplots,1):
    inputs[i] = directory+"/dat/"+inputs[i]
    outputs[i] = directory+"/eps/"+outputs[i]
    
    plt.figure(figsize=(8,6))
    axes = plt.gca()
    
    temp = pd.read_csv(inputs[i], delimiter = "\t", usecols=[0,1], header=None)
    x = np.array(temp[0].tolist(), dtype=float)
    y = np.array(temp[1].tolist(), dtype=float)
    
### Order parameter ###
    if i==0:
        plt.ylim(0,1)
        plt.plot(x,y,'ro',markevery=inv)
        plt.plot(x,np.repeat(np.average(y),len(x)), label='Average = %.4f' % np.average(y))

### Orientation ###
    elif i==1:
        fig = plt.figure()
        
        if NDIM==2:
            ax = fig.add_subplot(111)
            temp = pd.read_csv(inputs[i], delimiter = "\t", usecols=[0,1,2], names=["t", "x", "y"], header=None)
            pnt = ax.scatter(temp.x,temp.y,c=temp.t)
            cbar = plt.colorbar(pnt)
            ax.set_xlim(-1,1)
            ax.set_ylim(-1,1)
        
        if NDIM==3:
            ax = fig.add_subplot(111, projection='3d')
            temp = pd.read_csv(inputs[i], delimiter = "\t", usecols=[0,1,2,3], names=["t", "x", "y", "z"], header=None)
            pnt= ax.scatter(temp.x,temp.y,temp.z,marker='o',c=temp.t)
            cbar = plt.colorbar(pnt)
            ax.set_xlim(-1,1)
            ax.set_ylim(-1,1)
            ax.set_zlim(-1,1)
            fig.show()

### MSD ###
    elif i==2:
        plt.plot(x,y,'ro', label='data',markevery=inv)
        
        fitParams, fitCovariances = curve_fit(linear, x, y)
        m,b = fitParams
        
        if NDIM==2:
            D = m/4.0
            np.savetxt(summary3, ["MSD linear fit 4*D:\t %s" % m], fmt='%s')
            plt.plot(x,linear(x,m,b), label='4Dt \n D = %.3g'%D)
        if NDIM==3:
            D = m/6.0
            np.savetxt(summary3, ["MSD linear fit 6*D:\t %s" % m], fmt='%s')
            plt.plot(x,linear(x,m,b), label='6Dt \n D = %.3g'%D)
                
        np.savetxt(summary4, ["%s" % m], fmt='%s')
    
### Log MSD ###
    elif i==3:
        try:
            y = np.log(np.array(1.0)-y)
            plt.plot(x,y,'ro',markevery=inv)
        except Exception:
            sys.exc_clear()

### Density fluctuations: linear fit on log(x) and log(y) ###
    elif i==4:
        plt.loglog(x,y,'ro', label='data')
       
        xl = np.log(x)
        yl = np.log(y)
        
        fitParams, fitCovariances = curve_fit(linear, xl, yl)
        m,b = fitParams
        plt.loglog(x,powerlaw(x,m,np.power(np.e,b),0), label='power law m = %.3g'%m)
        
        np.savetxt(summary3, ["Density fluctuation scaling: \t %s" % m], fmt='%s')
        np.savetxt(summary4, ["%s" % m], fmt='%s')

### Autocorrelation function ###
    elif i==5:
        plt.plot(x,y,'ro', label='data')
        
        fitParams, fitCovariances = curve_fit(exponential, x, y)
        m,b,c = fitParams
        plt.plot(x,exponential(x,m,b,c), label = str(b)+'*$e^{-t/%.3g}$' % m)

### Correlation function ###
### Ignore the first data point in the fit, which deviates from the fit ###
### (neighboring particles repel each other). ###
    elif i==6:
        fig = plt.figure(figsize=(8,18))
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)
        
        # The first element in the list does not fit well because neighboring particles repel.
        x1 = x[1:]
        y1 = y[1:]
        xmax = x1[len(x1)-1]
        xrange = np.linspace(2,xmax,num=100)
        
        p=[-2, 1, 0] # Initial guess for power law fit
        fitParamsP, fitCovariancesP = curve_fit(powerlaw, x1, y1, p0=p)
        m,b,c = fitParamsP
        powerlawerr = np.sqrt(np.diag(fitCovariancesP))
        ax1.plot(xrange,powerlaw(xrange,m,b,c),     label='power law m = %.3g'%m)
        ax2.semilogy(xrange,powerlaw(xrange,m,b,0), label='power law m = %.3g'%m)
        ax3.loglog(xrange,powerlaw(xrange,m,b,0),   label='power law m = %.3g'%m)
        np.savetxt(summary3, ["Power law coefficient: \t %s" % b], fmt='%s')
        np.savetxt(summary4, ["%s" % b], fmt='%s')
        np.savetxt(summary3, ["Power law exponent: \t %s" % m], fmt='%s')
        np.savetxt(summary4, ["%s" % m], fmt='%s')
        np.savetxt(summary3, ["Error of power law fit: \t %s" % powerlawerr], fmt='%s')
        np.savetxt(summary4, ["%s" % powerlawerr], fmt='%s')
        
        ax1.plot(x,y,'ro', label='data')
        ax2.semilogy(x,y-c,'ro', label='data')
        ax3.loglog(x,y-c,'ro', label='data')
        
        p=[10, 1, 0]
        fitParamsE, fitCovariancesE = curve_fit(exponential, x1, y1, p0=p)
        m,b,c = fitParamsE
        experr = np.sqrt(np.diag(fitCovariances))
        ax1.plot(xrange,exponential(xrange,m,b,c),     label='exponential lc = %.3g'%m)
        ax2.semilogy(xrange,exponential(xrange,m,b,0), label='exponential lc = %.3g'%m)
        ax3.loglog(xrange,exponential(xrange,m,b,0),   label='exponential lc = %.3g'%m)
        np.savetxt(summary3, ["Exponential coefficient: \t %s" % b], fmt='%s')
        np.savetxt(summary4, ["%s" % b], fmt='%s')
        np.savetxt(summary3, ["Exponential correlation length: \t %s" % m], fmt='%s')
        np.savetxt(summary4, ["%s" % m], fmt='%s')
        np.savetxt(summary3, ["Error of exponential fit: \t %s" % experr], fmt='%s')
        np.savetxt(summary4, ["%s" % experr], fmt='%s')
        
        ax1.set_xlim(1,xmax)
        ax2.set_xlim(1,xmax)
        ax3.set_xlim(1,xmax)
        ax1.legend()
        ax2.legend()
        ax3.legend()
        ax1.grid(True)
        ax2.grid(True)
        ax3.grid(True)

### Pair correlation function ###
    elif i==7:
        plt.xlim(0,20)
        plt.plot(x,y,'r-')

### Velocity distribution function ###
    elif i==8:
        plt.plot(x,y,'ro')

    plt.suptitle(titles[i])
    plt.xlabel(axislabels[i][0])
    plt.ylabel(axislabels[i][1])
    plt.grid(True)
    plt.legend()
    plt.savefig(outputs[i], format="eps")
