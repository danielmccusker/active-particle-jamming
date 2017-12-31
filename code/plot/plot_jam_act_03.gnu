#!/bin/bash

# 01: Added MSD
# 01: Added position xavg and three individual cells
# 01: Fit discards first 10000 points
# 01: Plot for order plots every 10th point. This saves 90% space
# 02: Plot for MSD
# 02: Plot for MSD with errorbars
# 02: Plot auto correlation function
# 02: Plot means squared difference of radii of neighbouring particles
# 02: Rearranged columns
# 03: fixed directories


gnuplot << EOF

N = 1000            # Number of data points on a graph
int = ($3)/(N*10)            # Step size between data points
#int = 1

set terminal postscript eps color enhanced "Helvetica" 20
#directory = '/Users/Daniel1/Desktop/ActiveMatterResearch/clusterOutput/'.'$2/'
directory = '/Users/Daniel1/Desktop/ActiveMatterResearch/jamming-dynamics/output/'.'$2/'
#directory = '/home/dmccusker/remote/jamming-dynamics/output/'.'$2/'

# ------------------------------------------------------------------------------------------ #
# Plot order parameter

filein  = directory.'$1/dat/data.dat'
fileout = directory.'$1/eps/order.eps'

set output fileout

set title "Order"
set xlabel "time [steps]"
set ylabel "{/Symbol f}"
set yrange [0:1]

set key horizontal below height 3
set key box lt 2 lc -1 lw 3

#set fit errorvariables
#set fit logfile "/dev/null"
set fit quiet

f0(x) = A0
A0 = 0.5

fit f0(x) filein using 1:4 every int:int via A0

plot filein u 1:4 every int w points t '$1', f0(x) t sprintf("mean: A = %.4g, error = %.4g",A0, sqrt(FIT_WSSR / (FIT_NDF + 1 )))

# ------------------------------------------------------------------------------------------ #
# Plot pressure

fileout = directory.'$1/eps/pressure.eps'

set output fileout

set title "Pressure"
set xlabel "time [steps]"
set ylabel "P"
set yrange [0:0.5]

plot filein u 1:12 every int w points t '$1'

# ------------------------------------------------------------------------------------------ #
# Plot system orientation

fileout = directory.'$1/eps/orientation.eps'

set output fileout

set title "Global orientation"
set xlabel "time [steps]"
set ylabel "{/Symbol Y}"
set yrange [-pi:pi]

plot filein u 1:5 every int w points t '$1'

# ------------------------------------------------------------------------------------------ #
# Plot trajectory of the system's center of mass

fileout = directory.'$1/eps/COM.eps'

set output fileout

set title "Position with time of the original system"
set xlabel "x"
set ylabel "y"
set xrange [-750:750]
set yrange [-750:750]

plot filein u 2:3:1 every int with points palette t "CoM"

# ------------------------------------------------------------------------------------------ #
# Plot trajectories of three cells

fileout = directory.'$1/eps/x123.eps'

set output fileout

stats filein u 6:7 nooutput
L = (STATS_max_x > STATS_max_y) ? STATS_max_x : STATS_max_y
L = (-STATS_min_x > L) ? -STATS_min_x : L
L = (-STATS_min_y > L) ? -STATS_min_y : L
stats filein u 8:9 nooutput
L = (STATS_max_x > L) ? STATS_max_x : L
L = (STATS_max_y > L) ? STATS_max_y : L
L = (-STATS_min_x > L) ? -STATS_min_x : L
L = (-STATS_min_y > L) ? -STATS_min_y : L
stats filein u 10:11 nooutput
L = (STATS_max_x > L) ? STATS_max_x : L
L = (STATS_max_y > L) ? STATS_max_y : L
L = (-STATS_min_x > L) ? -STATS_min_x : L
L = (-STATS_min_y > L) ? -STATS_min_y : L

f(x) = x - 2*L*floor((x+L)/(2*L))

set title "Position with time of the original system"
set xlabel "x"
set ylabel "y"
set xrange [-1000:1000]
set yrange [-1000:1000]

plot filein u (f(\$6-\$2)):(f(\$7-\$3)):1 every int with points palette t "Cell 1", filein u (f(\$8-\$2)):(f(\$9-\$3)):1 every int with points palette t "Cell 2", filein u (f(\$10-\$2)):(f(\$11-\$3)):1 every int with points palette t "Cell 3"

# ------------------------------------------------------------------------------------------ #
# Pair correlation function

filein  = directory.'$1/dat/pairCorr.dat'
fileout = directory.'$1/eps/pairCorr.eps'

set output fileout

set title "Pair correlation function"
set xlabel "distance [@^{/=18--}a]"
set ylabel "g(r)"
set autoscale x
set autoscale y

stats filein using 2 every ::1 name "Y" nooutput
stats filein using 1 every ::Y_index_max::Y_index_max name "X" nooutput

plot filein using 1:2 every ::1::2000 t sprintf("$1: peak of g1 = %.4g at r = %.4g", Y_max, X_max) with linespoints

set print directory.'$1/dat/summary3.dat' append
print 'Pair correlation function peak:  ', Y_max
print 'Peaks at radius: ', X_max

set print directory.'$1/dat/summary4.dat' append
print Y_max
print X_max

# ------------------------------------------------------------------------------------------ #
# Long pair correlation function

fileout = directory.'$1/eps/pairCorrLong.eps'

set output fileout

set title "Pair correlation function"
set xlabel "distance [@^{/=18--}a]"
set ylabel "g(r)"
set autoscale x
set autoscale y

plot filein using 1:2 with linespoints

# ------------------------------------------------------------------------------------------ #
# Mean squared displacement

filein  = directory.'$1/dat/MSD.dat'
fileout = directory.'$1/eps/MSD.eps'

set output fileout

set title "MSD"
set xlabel "lagtime [steps]"
set ylabel "<(x({/Symbol t}) - x(0))^2>"

f(x) = 4*D*x+d
D=0.1
d=0.5
FIT_LIMIT = 1e-08
fit f(x) filein u 1:2 every int via d,D

plot filein u 1:2 every int t '$1', f(x) t sprintf("mean: D = %.4g",D)

set print directory.'$1/dat/summary3.dat' append
print 'effective diffusion constant:	', D

set print directory.'$1/dat/summary4.dat' append
print D

# ------------------------------------------------------------------------------------------ #
# Logarithm of MSD

fileout = directory.'$1/eps/logMSD.eps'
set output fileout

set title "log-lin MSD"
set xlabel "lagtime [steps]"
set ylabel "log <(x({/Symbol t}) - x(0))^2>"
unset logscale xy

plot filein u 1:4 every int t '$1', log(f(x)) t sprintf("mean: D = %.4g, error = %.4g",D, sqrt(FIT_WSSR / (FIT_NDF + 1 )))

# ------------------------------------------------------------------------------------------ #
# Plot density fluctuations

filein  = directory.'$1/dat/GNF.dat'
fileout = directory.'$1/eps/GNF.eps'

set logscale xy

f(x) = m*x+b
fit f(x) filein using 4:5 every ::0::9 via m, b
f2(x) = exp(b)*x**m

set output fileout

set title "Density fluctuations"
set xlabel "<A>"
set ylabel "rms fluctuation"

plot filein using 2:3 every ::0::9 with linespoints t '$1', x with lines t 'slope = 1', sqrt(x) with lines t 'slope = 0.5', f2(x) with lines lt rgb "#ff00ff" t sprintf("slope of fit = %.4g", m)

set print directory.'$1/dat/summary3.dat' append
print 'density fluctuations scale like: ', m
set print directory.'$1/dat/summary4.dat' append
print m

# ------------------------------------------------------------------------------------------ #
# Plot order autocorrelation function

filein  = directory.'$1/dat/oaf.dat'
fileout = directory.'$1/eps/oaf.eps'

set output fileout

set title "order autocorrelation function"
set xlabel "lagtime [steps]"
set ylabel "<(phi({/Symbol t}).phi(0))^2>"
unset logscale xy

plot filein using 1:2 t '$1' with linespoints

# ------------------------------------------------------------------------------------------ #
# Plot velocity autocorrelation function

filein  = directory.'$1/dat/vaf.dat'
fileout = directory.'$1/eps/vaf.eps'

set output fileout

set title "velocity autocorrelation function"
set xlabel "lagtime [steps]"
set ylabel "<(v({/Symbol t}).v(0))^2>"
unset logscale xy

plot filein using 1:2 t '$1' with linespoints

# ------------------------------------------------------------------------------------------ #
# Velocity distribution

filein  = directory.'$1/dat/velDist.dat'
fileout = directory.'$1/eps/velDist.eps'

set output fileout

gauss(x)=a/(sqrt(2*pi)*sigma)*exp(-(x-mean)**2/(2*sigma**2))
fit gauss(x) filein via a,sigma,mean

set title "Velocity distribution"
set xlabel "velocity [@^{/=18-}a/{/Symbol t}]"
set ylabel "count"

plot filein u 1:2 t '$1', gauss(x) t sprintf("$1: a = %.4g, mean = %.4g, sigma = %.4g",a,mean,sigma) with linespoints

set print directory.'$1/dat/summary3.dat' append
print 'average velocity:    ', mean
print 'velocity distribution width: ', sigma
print 'velocity distribution height:    ', a

set print directory.'$1/dat/summary4.dat' append
print mean
print sigma
print a

# ------------------------------------------------------------------------------------------ #
# Plot velocity correlation function

filein  = directory.'$1/dat/velCorr.dat'
fileout = directory.'$1/eps/velCorr.eps'

unset logscale xy
set output fileout

f(x) = A*x**b + C*exp(-x/d)
A=0.1
b=-0.5
C=0.1
d=10
FIT_MAXITER = 10
fit f(x) filein via A, b, C, d

set title "Velocity correlation function"
set xlabel "r"
set ylabel "g"

plot filein using 1:2 with linespoints t '$1', f(x) with lines lt rgb "#ff00ff" t sprintf("exponential plus power fit \n %.4g x^{%.4g} + %.4g e^{(-x/%.4g)}", A, b, C, d)

set print directory.'$1/dat/summary3.dat' append
print 'correlation power law strength:  ', A
print 'power:   ', b
print 'exponential strengh: ', C
print 'exponential correlation length:  ', d

set print directory.'$1/dat/summary4.dat' append
print A
print b
print C
print d

# ------------------------------------------------------------------------------------------ #
# Plot velocity correlation function on log-lin scale

filein  = directory.'$1/dat/velCorr.dat'
fileout = directory.'$1/eps/velCorrLogLin.eps'

set logscale y
set output fileout

set title "Velocity correlation function"
set xlabel "r"
set ylabel "C"

plot filein using 1:2 with linespoints t '$1'

# ------------------------------------------------------------------------------------------ #
# Plot velocity correlation function on log-log scale

filein  = directory.'$1/dat/velCorr.dat'
fileout = directory.'$1/eps/velCorrLogLog.eps'

set logscale xy
set output fileout

set title "Velocity correlation function"
set xlabel "r"
set ylabel "C"

plot filein using 1:2 with linespoints t '$1'

# ------------------------------------------------------------------------------------------ #
# Plot velocity correlation corrected function

filein  = directory.'$1/dat/velCorrC.dat'
fileout = directory.'$1/eps/velCorrC.eps'

unset logscale xy
set output fileout

f(x) = A*x**b + C*exp(-x/d)
A=0.1
b=-0.5
C=0.1
d=10
FIT_MAXITER = 10
fit f(x) filein via A, b, C, d

set title "Velocity correlation function"
set xlabel "r"
set ylabel "g"

plot filein using 1:2 with linespoints t '$1', f(x) with lines lt rgb "#ff00ff" t sprintf("exponential plus power fit \n %.4g x^{%.4g} + %.4g e^{(-x/%.4g)}", A, b, C, d)

# ------------------------------------------------------------------------------------------ #
# Plot velocity correlation function on log-lin scale

filein  = directory.'$1/dat/velCorrC.dat'
fileout = directory.'$1/eps/velCorrLogLinC.eps'

set logscale y
set output fileout

set title "Velocity correlation function"
set xlabel "r"
set ylabel "C"

plot filein using 1:2 with linespoints t '$1'

# ------------------------------------------------------------------------------------------ #
# Plot velocity correlation function on log-log scale

filein  = directory.'$1/dat/velCorrC.dat'
fileout = directory.'$1/eps/velCorrLogLogC.eps'

set logscale xy
set output fileout

set title "Velocity correlation function"
set xlabel "r"
set ylabel "C"

plot filein using 1:2 with linespoints t '$1'

# ------------------------------------------------------------------------------------------ #
# Plot angle correlation function

filein  = directory.'$1/dat/angleCorr.dat'
fileout = directory.'$1/eps/angleCorr.eps'

unset logscale xy
set output fileout

set title "Orientation correlation function"
set xlabel "r"
set ylabel "g"

plot filein using 1:2 with linespoints t '$1'

# ------------------------------------------------------------------------------------------ #
# Plot angle correlation function on log-lin scale

filein  = directory.'$1/dat/angleCorr.dat'
fileout = directory.'$1/eps/angleCorrLogLin.eps'

set logscale y
set output fileout

set title "Velocity correlation function"
set xlabel "r"
set ylabel "C"

plot filein using 1:2 with linespoints t '$1'

# ------------------------------------------------------------------------------------------ #
# Plot angle correlation function on log-log scale

filein  = directory.'$1/dat/angleCorr.dat'
fileout = directory.'$1/eps/angleCorrLogLog.eps'

set logscale xy
set output fileout

set title "Velocity correlation function"
set xlabel "r"
set ylabel "C"

plot filein using 1:2 with linespoints t '$1'

# ------------------------------------------------------------------------------------------ #
# Plot angle corrected correlation function

filein  = directory.'$1/dat/angleCorrC.dat'
fileout = directory.'$1/eps/angleCorrC.eps'

unset logscale xy
set output fileout

set title "Orientation correlation function"
set xlabel "r"
set ylabel "g"

plot filein using 1:2 with linespoints t '$1'

# ------------------------------------------------------------------------------------------ #
# Plot angle corrected correlation function on log-lin scale

filein  = directory.'$1/dat/angleCorrC.dat'
fileout = directory.'$1/eps/angleCorrLogLinC.eps'

set logscale y
set output fileout

set title "Velocity correlation function"
set xlabel "r"
set ylabel "C"

plot filein using 1:2 with linespoints t '$1'

# ------------------------------------------------------------------------------------------ #
# Plot angle corrected correlation function on log-log scale

filein  = directory.'$1/dat/angleCorrC.dat'
fileout = directory.'$1/eps/angleCorrLogLogC.eps'

set logscale xy
set output fileout

set title "Velocity correlation function"
set xlabel "r"
set ylabel "C"

plot filein using 1:2 with linespoints t '$1'

# ------------------------------------------------------------------------------------------ #

EOF

