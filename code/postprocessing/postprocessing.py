
### Makes plots for data from all runs in the phase diagram.
### 1D scaling exponent plot
### MSD...

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import pandas as pd
import numpy as np
import csv
import os
import math

# Still want to add order parameter
# 1D GNF plots
# double check linear fit for MSD

### Get the total number of folders in this directory, which includes "." ###

n = sum(os.path.isdir(i) for i in os.listdir("."))
directory = "./plots/"
if os.path.exists(directory):
    n = n-1
else:
    os.makedirs(directory)

runs = range(1, n+1, 1)     # "range" does not include its stopping value, so we increase n by 1.
nrow = 7                    # Number of self-propulsion values in the phase space
ncol = 10                   # Number of noise values in the phase space
ndensities = 3              # Number of densities in this set of runs

titles = [ "Density fluctuation scaling exponent m", "Logarithm of effective diffusion constant log(D)",
          "Binder cumulant G", "Order parameter $\phi$", "Order parameter variance $\Delta\phi$" ]

plot_titles = ["Pair_Correlation","Velocity_Distribution",
               "Mean_squared_Displacement","Velocity_correlation_function"]

my_dpi = 227                # Dots per inch (for making nice eps plots), this depends on your monitor.
colors = plt.cm.spectral(np.linspace(0,1,31))
lines = ['--', '-.', '-', ':', '--', '-.', '-', ':']

### We want to read, for each run, its ID, noise, self-propulsion force, density, effective diffusion
### constant, density fluctuation scaling exponent, average order parameter, order parameter variance,
### and binder cumulant. We also want the x- and y- values for its pair correlation function, velocity
### distribution, MSD curve, and velocity correlation function.

### "summary2.dat" contains the values output from the main code as follows:
### RunID: index 0
### Number of cells: 1
### Linear grid size: 2
### Number of time steps: 3
### Steps per unit time: 4
### Self-propulsion force: 5
### Noise: 6
### Density: 7
### Total simulation time in seconds: 8
### Number of Verlet list refreshes: 9
### Nearest neighbor correlation: 10 (ignore)
### Binder cumulant: 11
### Average order parameter: 12
### Order parameter variance: 13
### Average pressure: 14
### Pressure variance: 15

### "summary4.dat" contains the fit parameters output from gnuplot as follows
### (once saved in the list "values", each of these indices is increased by 16):
### Pair correlation peak: 0
### Pair correlation peak radius: 1
### Effective diffusion constant: 2
### Density fluctuation scaling exponent: 3
### Mode (most frequent value) of the velocity distribution: 4
### Velocity distribution width: 5
### Velocity distribution peak value: 6
### Correlation power law coefficient: 7
### Power of the correlation function decay: 8
### Correlation exponential coefficient: 9
### Correlation length of the exponential decay: 10

####################
### Read in data ###
####################

ln = []
ls = []
values = []
xdata  = [[] for i in range(len(plot_titles))]
ydata  = [[] for i in range(len(plot_titles))]

for i in runs:
    path = "run"+str(i)+"/dat/summary2.dat"
    infile = open(path, "rb")
    reader = csv.reader(infile, delimiter = "\t")
    temp = []
    for row in reader:
        temp.append(row[0])
    infile.close()

    path = "run"+str(i)+"/dat/summary4.dat"
    infile = open(path, "rb")
    reader = csv.reader(infile, delimiter = "\t")
    temp2 = []
    for row in reader:
        temp2.append(row[0])
    infile.close()

    if float(temp2[2]) < 0:   # In case the linear fit for the MSD curve happens to be negative
        temp2[2] = float(1e-12)

    add = temp+temp2
    values.append(add)

    plots = ["run"+str(i)+"/dat/pairCorr.dat", "run"+str(i)+"/dat/velDist.dat",
             "run"+str(i)+"/dat/MSD.dat", "run"+str(i)+"/dat/angleCorr.dat"]

    for j in range(0,4,1):
        temp = pd.read_csv(plots[j], delimiter = "\t",usecols=[0,1], names=["xtemp", "ytemp"])
        xdata[j].append(temp.xtemp)
        ydata[j].append(temp.ytemp)

values = np.array(values)

#values = values[values[:,7].argsort()]  # Sort data by increasing density
#values = values[values[:,5].argsort()]  # Sort data by increasing self-propulsion
#values = values[values[:,6].argsort()]  # Sort data by increasing noise
values = np.split(values,ndensities,0)  # Split data into arrays, each with a given density

for m in range(ncol):
    ln.append(float(values[0][m*nrow][6]))
for n in range(nrow):
    ls.append(float(values[0][n][5]))

#####################
### End read data ###
#####################

################
### 3D plots ###
################

### Plots for scaling exponent, diffusion constant, binder cumulant, order parameter and its fluctuations

threeDPlots = [ "density_fluctuations.eps", "diffusion.eps", "binder.eps", "order.eps",
               "order_fluctuations.eps" ]
plotIndices = [19,18,11,12,13]
xlabels = ["Noise parameter"]
ylabels = ["Self-propulsion force"]
zlabels = [ "m", "log(D)", "G", r'$\phi$', r'$\Delta\phi$' ]
zlim = [ [0,1], [-7,-2], [0,1], [0,1], [0,0.005] ]

# Define color ranges:
zcol = [[0,1],             # Density fluctuations range from 0 to 1, > 0.5 indicates activity
        [-7,-3.6],         # log(D): for 200,000 time steps, log10(1/200000) = -5.3 is jammed
        [0,1],             # Binder cumulant goes from 1/3 to 2/3 in 2D
        [0,1],             # Order parameter (0,1)
        [0,0.005]]         # Order parameter fluctuations

zmap = [cm.bwr_r, cm.bwr_r, cm.Purples, cm.Purples, cm.Purples]

for p in range(0,ndensities,1):
    k=0
    for j in plotIndices:
        fig = plt.figure(figsize=(16, 6))
        ax = fig.add_subplot(121, projection='3d')
        z = np.zeros((ncol,nrow))
        for m in range(ncol):
            for n in range(nrow):
                if j==18:   # Logarithm of diffusion constant
                    z[m][n] = math.log10(float(values[p][(m*nrow)+n][j]))
                else:
                    z[m][n] = float(values[p][(m*nrow)+n][j])
        z = z.transpose()
        X, Y = np.meshgrid(ln, ls)
        surf = ax.plot_surface(X, Y, z, vmin=zcol[k][0], vmax= zcol[k][1], cmap=zmap[k])
        ax.set_title(titles[k])
        ax.set_zlim(zlim[k][0], zlim[k][1])
        ax.set_xlabel(xlabels[0])
        ax.set_ylabel(ylabels[0])
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
        ax.set_zlabel(zlabels[k])
        fig.colorbar(surf)

        ax = fig.add_subplot(122)
        scatter = ax.scatter(X, Y, c=z, vmin=zcol[k][0], vmax= zcol[k][1], cmap=zmap[k])
        ax.set_title(titles[k])
        ax.set_xlabel(xlabels[0])
        ax.set_ylabel(ylabels[0])

        fig.tight_layout()
        plt.savefig(directory+str(p)+threeDPlots[k], format="eps", dpi=my_dpi/4)
        k=k+1

####################
### End 3D Plots ###
####################

##############################################
### Plot calculated functions for each run ###
##############################################

# Plt positions subplots starting from the top left. Rearrange the indices to place graphs
# in the correct location in the phase diagram. This index depends on how the run numbers are sequenced
# in the initial script which generated the run names.

subplot_arrange = []
for i in range(0, ncol, 1):
    subplot_arrange.append([])
    for j in range(0, nrow, 1):
        subplot_arrange[i].append((j+1)+(i*nrow))
subplot_arrange_temp = list(map(list, zip(*subplot_arrange))) # Invert rows and columns
subplot_arrange_temp = list(reversed(subplot_arrange_temp))   # Reverse index numbering from bottom to top
subplot_list = []
for i in range(0, nrow, 1):
    for j in range(0, ncol, 1):
        subplot_list.append(subplot_arrange_temp[i][j])

# Arrange pair correlation, velocity distribution, and MSD curves in the phase diagram.

for j in range(0,3,1):
    plt.figure(figsize=(8*ncol, 6*nrow))
    plt.suptitle(plot_titles[j])
    for i in runs:
        k = subplot_list[i-1]
        plt.subplot(nrow,ncol,i)
        axes = plt.gca()
        if j==0:
            axes.set_xlim([0,15]) # Show only local pair correlation structure
            axes.set_ylim([0,5])
            # In case the pair correlation function isn't properly normalized:
            finalvalue = ydata[j][k-1][len(ydata[j][k-1]) - 1]
            ydata[j][k-1] = [ m/finalvalue for m in ydata[j][k-1] ]
        if j==1:
            axes.set_ylim([0,0.04])
        plt.grid(True)
        plt.title("run"+str(k))
        plt.plot(xdata[j][k-1],ydata[j][k-1],color='red')
    plt.subplots_adjust(wspace=0.25, hspace=0.25)
    plt.savefig(directory+plot_titles[j]+".eps", format="eps", dpi=my_dpi/4)

### For the correlation function, we want log-log and log-lin plots to distinguish exponential/power law
### decay. For each noise (x) value, we put all self-propulsions on the same subplot.

xlim = [ [1,120] ]
ylim = [ [1e-2,0.5] ]

plt.figure(figsize=(8*ncol, 3*nrow))
for m in range(0,3):
    for i in range(0,ncol,1):
        plt.subplot(3,ncol,i+m*ncol+1)
        axes = plt.gca()
        axes.set_xlim([xlim[0][0], xlim[0][1]])
        axes.set_ylim([ylim[0][0], ylim[0][1]])
        plt.xlabel('r')
        plt.ylabel('g')
        plt.grid(True)
        for j in range(0,nrow,1):
            k = subplot_arrange[i][j]
            x = xdata[3][k-1]
            y = ydata[3][k-1]
            if m==0:
                plt.plot(x,y)
            if m==1:
                plt.semilogy(x,y)
            if m==2:
                plt.loglog(x,y)

plt.suptitle(plot_titles[3])
plt.savefig(directory+"CorrelationFunction.eps", format="eps", dpi=my_dpi/4)

### Plot a few curves at the lowest self-propulsion speed, for clarity
### Change index for xdata,ydata to pick a different speed.

plt.figure(figsize=(9,6), dpi=my_dpi/4)
plt.suptitle(plot_titles[3])

for i in range(0,3):
    plt.subplot(2,2,i+1)
    axes = plt.gca()
    axes.set_xlim([xlim[0][0], xlim[0][1]])
    axes.set_ylim([ylim[0][0], ylim[0][1]])
    plt.ylabel('g')
    plt.grid(True)
    for k in range(1,ncol):
        x = xdata[3][k*nrow-1]
        y = ydata[3][k*nrow-1]
        if i==0:
            plt.plot(x,y,color=colors[k])
        if i==1:
            plt.semilogy(x,y,color=colors[k])
        if i==2:
            plt.loglog(x,y,color=colors[k])

runtitles = np.arange(0.4, 0.57, 0.005)
plt.legend( runtitles, loc='upper left', bbox_to_anchor=(1.2, 1.05), ncol=3)
plt.savefig(directory+"fewcorrelations.eps", format="eps", dpi=my_dpi/4)

##########################
### End plot functions ###
##########################

##########################
###    End script      ###
##########################

