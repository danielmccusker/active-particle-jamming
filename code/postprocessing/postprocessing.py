#######################################################################
###  First, loop through all runs to create plots and fit data.		###
###																	###
###  Then, build a phase diagram.									###
#######################################################################

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np
import csv
import os
import math
import sys

###################################### ******************** ######################################
###################################### 	    User Input		######################################
###################################### ******************** ######################################

### In part 1, we plot and fit data for each run. The fitting takes a little while.
### In part 2, we take all runs and compile a phase diagram.

### Enter 1 for each part you'd like to run, 0 if not.

fit =
phase =

### Check that the local directory is correct:

#location = '/Users/Daniel1/Desktop/ActiveMatterResearch/jamming-dynamics/remote_output/'
#location = '/Users/Daniel1/Desktop/ActiveMatterResearch/jamming-dynamics/local_output/'

### Enter the name of the set of runs (with a / at the end).

fullRunID =

### Enter the number of dimensions in the simulations.

NDIM =

### Enter the number of time steps in the simulation.

nsteps =

### Enter the number of self-propulsion values (nrow), noise values (ncol), and densities.
### Also enter the values of the density as a list of strings, for plotting.

nrow =                    # Number of self-propulsion values in the phase space
ncol =                    # Number of noise values in the phase space
ndensities =

densities = 

###################################### ******************** ######################################
######################################    End User Input    ######################################
###################################### ******************** ######################################

### Make a list of runs [ 1, 2, ... , N ].

runs = range(1, nrow*ncol*ndensities+1, 1)

### We print a maximum of 1000 points per line to avoid the graphs being too crowded.

npoints = 1000
inv = int(nsteps/npoints)
if inv==0:
	inv = 1

#######################################
### Part 1: Curve fitting and plots ###
#######################################
### For a single run, plot:         ###
### Order parameter over the run    ###
### System orientation over the run ###
### Mean-squared displacement       ###
### Log(1-MSD)                      ###
### Density fluctuation scaling     ###
### Autocorrelation function        ###
### Correlation function            ###
### Correlation function log-log    ###
### Correlation function log-lin    ###
### Pair correlation function       ###
### Velocity distribution           ###
#######################################

if fit:

### Define the fitting functions.

	def linear(x,m,b):
		return m*x+b
	def exponential(x,m,b,c):
		return b*np.power( np.e, -1.0*x/m ) + c
	def exponentialZero(x,m,b):
		return b*np.power( np.e, -1.0*x/m )
	def powerlaw(x,m,b,c):
		return b*x**m + c

### Title, axis, I/O information for the plots. ###

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

	nplots = len(inputs)

### Make a loop over run names, and first open each directory.

	for r in runs:

		print "Fitting run "+str(r)+" of "+str(nrow*ncol*ndensities)+"..."
		
		directory = location+fullRunID+'run'+str(r)
		
		if os.path.exists(directory+"/dat/summary3.dat"):
			os.system("rm "+directory+"/dat/summary3.dat")
		if os.path.exists(directory+"/dat/summary4.dat"):
			os.system("rm "+directory+"/dat/summary4.dat")

		summary3 = open(directory+"/dat/summary3.dat","w")
		summary4 = open(directory+"/dat/summary4.dat","w")

### For each run, loop over the desired plots.

		for i in range(0,nplots,1):
			readin  = directory+"/dat/"+inputs[i]
			readout = directory+"/eps/"+outputs[i]
			
			plt.figure(figsize=(8,6))
			axes = plt.gca()
			
			temp = pd.read_csv(readin, delimiter = "\t", usecols=[0,1], header=None)
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
					temp = pd.read_csv(readin, delimiter = "\t", usecols=[0,1,2], names=["t", "x", "y"], header=None)
					pnt = ax.scatter(temp.x,temp.y,c=temp.t)
					cbar = plt.colorbar(pnt)
					ax.set_xlim(-1,1)
					ax.set_ylim(-1,1)
				
				if NDIM==3:
					ax = fig.add_subplot(111, projection='3d')
					temp = pd.read_csv(readin, delimiter = "\t", usecols=[0,1,2,3], names=["t", "x", "y", "z"], header=None)
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
				
				ax1.plot(x,y,'ro', label='data')
				p=[-1, 1, 0] # Initial guess for power law fit
				try:
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
					
					ax2.semilogy(x,y-c,'ro', label='data')
					ax3.loglog(x,y-c,'ro', label='data')

				except RuntimeError:
					np.savetxt(summary3, ["Power law coefficient: \t %s" % 0], fmt='%s')
					np.savetxt(summary4, ["%s" % 0], fmt='%s')
					np.savetxt(summary3, ["Power law exponent: \t %s" % 0], fmt='%s')
					np.savetxt(summary4, ["%s" % 0], fmt='%s')
					np.savetxt(summary3, ["Error of power law fit: \t %s" % 0], fmt='%s')
					np.savetxt(summary4, ["%s" % 0], fmt='%s')

				p=[10, 1]
				try:
					fitParamsE, fitCovariancesE = curve_fit(exponentialZero, x1, y1, p0=p)
					m,b = fitParamsE
					experr = np.sqrt(np.diag(fitCovariances))
					ax1.plot(xrange,exponentialZero(xrange,m,b),     label='exponential lc = %.3g'%m)
					ax2.semilogy(xrange,exponentialZero(xrange,m,b), label='exponential lc = %.3g'%m)
					ax3.loglog(xrange,exponentialZero(xrange,m,b),   label='exponential lc = %.3g'%m)
					np.savetxt(summary3, ["Exponential coefficient: \t %s" % b], fmt='%s')
					np.savetxt(summary4, ["%s" % b], fmt='%s')
					np.savetxt(summary3, ["Exponential correlation length: \t %s" % m], fmt='%s')
					np.savetxt(summary4, ["%s" % m], fmt='%s')
					np.savetxt(summary3, ["Error of exponential fit: \t %s" % experr], fmt='%s')
					np.savetxt(summary4, ["%s" % experr], fmt='%s')

				except RuntimeError:
					np.savetxt(summary3, ["Exponential coefficient: \t %s" % 0], fmt='%s')
					np.savetxt(summary4, ["%s" % 0], fmt='%s')
					np.savetxt(summary3, ["Exponential correlation length: \t %s" % 0], fmt='%s')
					np.savetxt(summary4, ["%s" % 0], fmt='%s')
					np.savetxt(summary3, ["Error of exponential fit: \t %s" % 0], fmt='%s')
					np.savetxt(summary4, ["%s" % 0], fmt='%s')
					
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
			plt.savefig(readout, format="eps")
			plt.close('all')

###############################################
### Part 2: Phase diagram           		###
###############################################
### Plot in the phase space (3D plots):		###
### Density fluctuation scaling exponent 	###
### Log of effective diffusion constant	 	###
### Binder cumulant						 	###
### Order parameter						 	###
### Order parameter fluctuations         	###
### Fitted exponential correlation length	###
### Fitted power law                    	###
###											###
### Also plot these values for a fixed      ###
### 	self-propulsion force.				###
### Make some nice correlation				###
###		function plots.						###
###############################################

if phase:

### Open a folder for the phase space plots to go.

	directory = "./plots/"
	exist = os.path.exists(directory)
	if (not exist):
		os.makedirs(directory)

#### Dots per inch (for making nice eps plots), this depends on your monitor.

	my_dpi = 227
	
### We want to read, for each run, its ID, noise, self-propulsion force, density, effective diffusion
### constant, density fluctuation scaling exponent, average order parameter, order parameter variance,
### and binder cumulant. We also want the x- and y- values for its pair correlation function, velocity
### distribution, MSD curve, and velocity correlation function.

	titles = [ "Density fluctuation scaling exponent m", "Logarithm of effective diffusion constant log(D)",
			  "Binder cumulant G", "Order parameter $\phi$", "Order parameter variance $\Delta\phi$",
			  "Correlation length", "Power law" ]

	plot_titles = ["Pair_Correlation","Velocity_Distribution",
				   "Mean_squared_Displacement","Velocity_correlation_function","Order_parameter"]

### Read in data.

	ln = []
	ls = []
	values = []
	xdata  = [[] for i in range(len(plot_titles))]
	ydata  = [[] for i in range(len(plot_titles))]

	print "Loading data..."
	
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

		add = temp+temp2
		values.append(add)

		datafiles = ["run"+str(i)+"/dat/pairCorr.dat", "run"+str(i)+"/dat/velDist.dat",
				 "run"+str(i)+"/dat/MSD.dat", "run"+str(i)+"/dat/corr.dat",
				 "run"+str(i)+"/dat/order.dat" ]

		for j in range(0,len(datafiles),1):
			temp = pd.read_csv(datafiles[j], delimiter = "\t",usecols=[0,1], names=["xtemp", "ytemp"])
			xdata[j].append(temp.xtemp)
			ydata[j].append(temp.ytemp)

	# Split values into arrays, each with a given density.
	
	values = np.array(values)
	values = np.split(values,ndensities,0)

	for m in range(ncol):
		ln.append(float(values[0][m*nrow][6]))
	for l in range(nrow):
		ls.append(float(values[0][l][5]))

### End read data.

### "values" now stores points for phase diagram.
### "values" is indexed [which density][which (ls,ln) point][which quantity].
### xdata and ydata store the output generated by the main code.
### They are indexed [which run][which data set]

### The key to the indices is as follows:

	# RunID: index 0
	# Number of cells: 1
	# Linear grid size: 2
	# Number of time steps: 3
	# Steps per unit time: 4
	# Self-propulsion force: 5
	# Noise: 6
	# Density: 7
	# Total simulation time in seconds: 8
	# Number of Verlet list refreshes: 9
	# Nearest neighbor correlation: 10 (ignore this one)
	# Binder cumulant: 11
	# Average order parameter: 12
	# Order parameter variance: 13
	# Effective diffusion constant: 14
	# Density fluctuation scaling exponent: 15
	# Correlation power law coefficient: 16
	# Power of the correlation function decay: 17
	# "goodness" of power law fit: 18
	# Correlation exponential coefficient: 19
	# Correlation length of the exponential decay: 20
	# "goodness" of exponential fit: 21

### Generate 3D plots.

	print "Generating 3D plots for "+str(len(values))+" densities, with "+str(len(values[0]))+ \
	" mesh points each. I imported "+str(len(values[0][0]))+ " quantities for each run."

	plotIndices = [15,14,11,12,13,20,17]

	output_names = [ "density_fluctuations.eps", "diffusion.eps", "binder.eps", "order.eps",
				   "order_fluctuations.eps", "correlationLength.eps", "powerLaw.eps" ]

	xlabels = ["Noise parameter"]
	ylabels = ["Self-propulsion force"]
	zlabels = [ "m", "log10(D)", "G", r'$\phi$', r'$\Delta\phi$', "lc", "m" ]

	zlim = [ [0,1], [-10,-1], [0,1], [0,1], [0,0.005], [1,40], [-3,0] ]

	# Define color ranges:
	# For log(D), define the midpoint which separates jammed from unjammed
	if NDIM==2:
		mid = -np.log10(4*nsteps)
	if NDIM==3:
		mid = -np.log10(6*nsteps)

	zcolor = [ [0,1], [mid-6,mid+6], [0,1], [0,1], [0,0.005], [1,40], [-3,0] ]
	zmap = [cm.bwr_r, cm.bwr_r, cm.Purples, cm.Purples, cm.Purples, cm.Purples, cm.Purples]

	k=0
	for j in plotIndices:
		fig = plt.figure(figsize=(16, ndensities*6))
		for p in range(0,ndensities,1):
			ax = fig.add_subplot(3,2,2*p+1, projection='3d')
			z = np.zeros((ncol,nrow))
			for m in range(ncol):
				for l in range(nrow):
					r = (m*nrow)+l
					if j==14:
						if NDIM==2:
							D = float(values[p][r][j])/4.0
						if NDIM==3:
							D = float(values[p][r][j])/6.0
						
						# In case the linear fit for the MSD curve happens to be negative:
						if D < 0:
							D = 1e-8
						z[m][l] = math.log10(D)			# Logarithm of diffusion constant
						
					else:
						z[m][l] = float(values[p][r][j])
		
			z = z.transpose()
			X, Y = np.meshgrid(ln, ls)
			surf = ax.plot_surface(X, Y, z, vmin=zcolor[k][0], vmax= zcolor[k][1], cmap=zmap[k])
			ax.set_title(titles[k])
			ax.set_zlim(zlim[k][0], zlim[k][1])
			ax.set_xlabel(xlabels[0])
			ax.set_ylabel(ylabels[0])
			ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
			ax.set_zlabel(zlabels[k])
			fig.colorbar(surf)

			ax = fig.add_subplot(3,2,2*p+2)
			scatter = ax.scatter(X, Y, c=z, vmin=zcolor[k][0], vmax= zcolor[k][1], cmap=zmap[k])
			ax.set_title(titles[k])
			ax.set_xlabel(xlabels[0])
			ax.set_ylabel(ylabels[0])

		fig.tight_layout()
		plt.savefig(directory+output_names[k], format="eps", dpi=my_dpi/4)
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
			subplot_arrange[i].append(j+(i*nrow))
	subplot_arrange_temp = list(map(list, zip(*subplot_arrange))) # Invert rows and columns
	subplot_arrange_temp = list(reversed(subplot_arrange_temp))   # Reverse index numbering from bottom to top
	subplot_list = []
	for i in range(0, nrow, 1):
		for j in range(0, ncol, 1):
			subplot_list.append(subplot_arrange_temp[i][j])

	# Arrange pair correlation, velocity distribution, MSD curves, and order parameter in the phase diagram.

	print "Arranging graphs in the phase diagram..."
	print "(There are many graphs, this may take a few minutes.)"

	nperplot = nrow*ncol
	for p in range(0,ndensities,1):
		for j in range(0,5,1):
			fig = plt.figure(figsize=(8*ncol, 6*nrow))
			plt.suptitle(plot_titles[j])
			for i in range(0,nperplot,1):
				# print subplot_list[i]
				k = subplot_list[i]
				#print k
				plt.subplot(nrow,ncol,i+1)
				axes = plt.gca()
				if j==0:
					axes.set_xlim([0,15]) # Show only local pair correlation structure
				if j==1:
					axes.set_ylim([0,0.04])
				if j==4:
					axes.set_ylim([0,1])
				plt.grid(True)
				plt.title("run"+str(k+1))
				#print xdata[j][k+p*nperplot]
				#print ydata[j][k+p*nperplot]
				plt.plot(xdata[j][k+p*nperplot],ydata[j][k+p*nperplot],color='red')
			plt.subplots_adjust(wspace=0.25, hspace=0.25)
			plt.savefig(directory+densities[p]+plot_titles[j]+".eps", format="eps", dpi=my_dpi/4)
			plt.close(fig)

	### For the correlation function, we want log-log and log-lin plots to distinguish exponential/power law
	### decay. For each noise (x) value, we put all self-propulsions on the same subplot.

	xlim = [ [1,120] ]
	ylim = [ [1e-2,0.5] ]

	colors = plt.cm.spectral(np.linspace(0,1,31))
	lines = ['--', '-.', '-', ':', '--', '-.', '-', ':']

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

	runtitles = ln
	plt.legend( runtitles, loc='upper left', bbox_to_anchor=(1.2, 1.05), ncol=3)
	plt.savefig(directory+"fewcorrelations.eps", format="eps", dpi=my_dpi/4)

	##########################
	### End plot functions ###
	##########################

##########################
###    End script      ###
##########################

