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
#int = ($3)/(N*100)            # Step size between data points
int = 1

set terminal postscript eps color enhanced "Helvetica" 20
#directory = '/Users/Daniel1/Desktop/ActiveMatterResearch/jamming-dynamics/output/'.'$2/'
directory = '/home/dmccusker/remote/jamming-dynamics/output/'.'$2/'

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

plot filein u 1:5 every int w points t '$1', f0(x) t sprintf("mean: A = %.4g, error = %.4g",A0, sqrt(FIT_WSSR / (FIT_NDF + 1 )))

# ------------------------------------------------------------------------------------------ #
# Plot system orientation

fileout = directory.'$1/eps/orientation.eps'

set output fileout

set title "Global orientation"
set xlabel "time [steps]"
set ylabel "{/Symbol Y}"
set yrange [-pi:pi]

plot filein u 1:5 every int w points t '$1'

## ------------------------------------------------------------------------------------------ #
## Plot average velocity
#
#fileout = directory.'$1/eps/vavg.eps'
#
#set output fileout
#
#set title "Average velocity"
#set xlabel "time [steps]"
#set ylabel "v_{gem}"
#set autoscale y
#
#f1(x) = A1
#A1 = 0.5
#
#fit f1(x) filein using 1:4 every int:int via A1
#
#plot filein u 1:4 every int::10 t '$1', f1(x) t sprintf("mean: A = %.4g, error = %.4g",A1, sqrt(FIT_WSSR / (FIT_NDF + 1 )))

# ------------------------------------------------------------------------------------------ #
# Plot trajectory of the system's center of mass

fileout = directory.'$1/eps/COM.eps'

set output fileout

set title "Position with time of the original system"
set xlabel "x"
set ylabel "y"
set xrange [-500:500]
set yrange [-500:500]

plot filein u 2:3:1 every int with points palette t "CoM"

# ------------------------------------------------------------------------------------------ #
# Plot trajectories of three cells

fileout = directory.'$1/eps/x123.eps'

set output fileout

#stats filein u 6:7 nooutput
#L = (STATS_max_x > STATS_max_y) ? STATS_max_x : STATS_max_y
#L = (-STATS_min_x > L) ? -STATS_min_x : L
#L = (-STATS_min_y > L) ? -STATS_min_y : L
#stats filein u 8:9 nooutput
#L = (STATS_max_x > L) ? STATS_max_x : L
#L = (STATS_max_y > L) ? STATS_max_y : L
#L = (-STATS_min_x > L) ? -STATS_min_x : L
#L = (-STATS_min_y > L) ? -STATS_min_y : L
#stats filein u 10:11 nooutput
#L = (STATS_max_x > L) ? STATS_max_x : L
#L = (STATS_max_y > L) ? STATS_max_y : L
#L = (-STATS_min_x > L) ? -STATS_min_x : L
#L = (-STATS_min_y > L) ? -STATS_min_y : L

#f(x) = x - 2*L*floor((x+L)/(2*L))

#set title "Position with time of the original system"
#set xlabel "x"
#set ylabel "y"
#set xrange [-30:30]
#set yrange [-30:30]

#plot filein u (f(\$6-\$2)):(f(\$7-\$3)):1 every int with points palette t "Cell 1", filein u (f(\$8-\$2)):(f(\$9-\$3)):1 every int with points palette t "Cell 2", filein u (f(\$10-\$2)):(f(\$11-\$3)):1 every int with points palette t "Cell 3"

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

plot filein using 1:2 t sprintf("$1: peak of g1 = %.4g at r = %.4g", Y_max, X_max) with linespoints

# ------------------------------------------------------------------------------------------ #
# Velocity distribution

filein  = directory.'$1/dat/velDist.dat'
fileout = directory.'$1/eps/velDist.eps'

set output fileout

set title "Velocity distribution"
set xlabel "velocity [@^{/=18-}a/{/Symbol t}]"
set ylabel "count"

plot filein u 1:2 t '$1'

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
#unset fit quiet
FIT_LIMIT = 1e-08
fit f(x) filein u 1:2 every int via d,D

plot filein u 1:2 every int t '$1', f(x) t sprintf("mean: D = %.4g",D)# , error = %.4g",D, sqrt(FIT_WSSR / (FIT_NDF + 1 )))

set print directory.'$1/dat/summary.dat' append
print 'effective diffusion constant:	', D

#fileout = directory.'$1/eps/MSDlogerror.eps'

#set output fileout
#set logscale xy

#plot filein u 1:2:3 every int w errorbars t '$1', f(x) t sprintf("mean: D = %.4g, error = %.4g",D, sqrt(FIT_WSSR / (FIT_NDF + 1 )))

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
# Plot velocity autocorrelation function

filein  = directory.'$1/dat/vaf.dat'
fileout = directory.'$1/eps/vaf.eps'

set output fileout

set title "velocity autocorrelation function"
set xlabel "lagtime [steps]"
set ylabel "<(v({/Symbol t}).v(0))^2>"
unset logscale xy

plot filein using 1:2 every int t '$1' with linespoints

# ------------------------------------------------------------------------------------------ #
# Plot density fluctuations

filein  = directory.'$1/dat/GNF.dat'
fileout = directory.'$1/eps/GNF.eps'

set logscale xy

f(x) = m*x+b
fit f(x) filein using 4:5 via m,b
f2(x) = exp(b)*x**m

set output fileout

set title "Density fluctuations"
set xlabel "<A>"
set ylabel "rms fluctuation"

plot filein using 2:3 with linespoints t '$1', x with lines t 'slope = 1', sqrt(x) with lines t 'slope = 0.5', f2(x) with lines lt rgb "#ff00ff" t sprintf("slope of fit = %.4g", m)

set print directory.'$1/dat/summary.dat' append
print 'density fluctuations scale like:	', m

# ------------------------------------------------------------------------------------------ #
EOF
