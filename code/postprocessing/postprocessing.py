import csv
import os
import numpy as np
import math

n = sum(os.path.isdir(i) for i in os.listdir(".")) #get the total number of folders in this directory, which includes "."
outfile = open("temp.dat","a")

for i in range(0, n-1): #does not include the stopping value
	path = "run"+str(i)+"/dat/summary.dat"
	infile = open(path, "rb")
	reader = csv.reader(infile, delimiter = "\t")
	names = []
	values = []
	for row in reader:
		names.append(row[0])
		values.append(row[1])
	infile.close()
	outfile.write(values[5]+"\t"+values[4]+"\t"+values[6]+"\t"+str(math.log10(float(values[7])))+"\t"+values[8]+"\n")
	## noise ## self-propulsion force ## density ## log(diffusion constant) ## scaling exponent ##
outfile.close()

data = np.loadtxt("temp.dat", delimiter='\t')
os.system("rm -f temp.dat")
data = data[data[:,2].argsort()] #sort rows by increasing density
data = np.split(data,4,0) #split into four arrays, each with a given density
np.savetxt("dens0.dat", data[0])
np.savetxt("dens1.dat", data[1])
np.savetxt("dens2.dat", data[2])
np.savetxt("dens3.dat", data[3])
os.system("./GNFvis.gnu")
os.system("./MSDvis.gnu")
