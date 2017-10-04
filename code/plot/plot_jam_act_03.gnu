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

set terminal postscript eps color enhanced "Helvetica" 20
directory = '/home/dmccusker/remote/jamming_cluster/output/'.'$2/'
filein  = directory.'$1/dat/data.dat'
fileout = directory.'$1/eps/order.eps'

set output fileout

set title "Order"
set xlabel "time [steps]"
set ylabel "{/Symbol f}"
set yrange [0:1]

set key horizontal below height 3
set key box lt 2 lc -1 lw 3

set fit errorvariables
set fit logfile "/dev/null"
set fit quiet

f0(x) = A0
A0 = 0.5

fit f0(x) filein using 1:5 every ::1000 via A0

plot filein u 1:5 every 10 t '$1', f0(x) t sprintf("mean: A = %.4g, error = %.4g",A0, sqrt(FIT_WSSR / (FIT_NDF + 1 )))

# ------------------------------------------------------------------------------------------ #

fileout = directory.'$1/eps/orientation.eps'

set output fileout

set title "Global orientation"
set xlabel "time [steps]"
set ylabel "{/Symbol Y}"
set yrange [-pi:pi]

plot filein u 1:6 every 10 t '$1'

# ------------------------------------------------------------------------------------------ #

fileout = directory.'$1/eps/vgem.eps'

set output fileout

set title "Average velocity"
set xlabel "time [steps]"
set ylabel "v_{gem}"
set autoscale y

f1(x) = A1
A1 = 0.5

fit f1(x) filein using 1:4 every ::1000 via A1

plot filein u 1:4 every 10::10 t '$1', f1(x) t sprintf("mean: A = %.4g, error = %.4g",A1, sqrt(FIT_WSSR / (FIT_NDF + 1 )))

# ------------------------------------------------------------------------------------------ #

fileout = directory.'$1/eps/adif2.eps'

set output fileout

set title "Mean squared difference of radii"
set xlabel "time [steps]"
set ylabel "<(a_i - a_j)^2>"

plot filein u 1:7 every 10 with lines

# ------------------------------------------------------------------------------------------ #

fileout = directory.'$1/eps/xavg.eps'

set output fileout

set title "Position with time of the original system"
set xlabel "x"
set ylabel "y"

plot filein u 2:3:1 every 100 with points palette t "CoM"

# ------------------------------------------------------------------------------------------ #

stats filein u 8:9 nooutput
L = (STATS_max_x > STATS_max_y) ? STATS_max_x : STATS_max_y
L = (-STATS_min_x > L) ? -STATS_min_x : L
L = (-STATS_min_y > L) ? -STATS_min_y : L
stats filein u 10:11 nooutput
L = (STATS_max_x > L) ? STATS_max_x : L
L = (STATS_max_y > L) ? STATS_max_y : L
L = (-STATS_min_x > L) ? -STATS_min_x : L
L = (-STATS_min_y > L) ? -STATS_min_y : L
stats filein u 12:13 nooutput
L = (STATS_max_x > L) ? STATS_max_x : L
L = (STATS_max_y > L) ? STATS_max_y : L
L = (-STATS_min_x > L) ? -STATS_min_x : L
L = (-STATS_min_y > L) ? -STATS_min_y : L

f(x) = x - 2*L*floor((x+L)/(2*L))
fileout = directory.'$1/eps/x123.eps'

set output fileout

set title "Position with time of the original system"
set xlabel "x"
set ylabel "y"
set xrange [-35:35]
set yrange [-35:35]

plot filein u (f(\$8-\$2)):(f(\$9-\$3)):1 every 100 with points palette t "Cell 1", filein u (f(\$10-\$2)):(f(\$11-\$3)):1 every 100 with points palette t "Cell 2", filein u (f(\$12-\$2)):(f(\$13-\$3)):1 every 100 with points palette t "Cell 3"

# ------------------------------------------------------------------------------------------ #

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

plot filein using 1:2 t sprintf("$1: peak of g1 = %.4g at r = %.4g", Y_max, X_max) with lines

# ------------------------------------------------------------------------------------------ #

filein  = directory.'$1/dat/velDist.dat'
fileout = directory.'$1/eps/velDist.eps'

set output fileout

set title "Velocity distribution"
set xlabel "velocity [@^{/=18-}a/{/Symbol t}]"
set ylabel "count"

plot filein u 1:2 t '$1'

# ------------------------------------------------------------------------------------------ #

filein  = directory.'$1/dat/MSD.dat'
fileout = directory.'$1/eps/MSD.eps'

set output fileout

set title "MSD"
set xlabel "lagtime [steps]"
set ylabel "<(x({/Symbol t}) - x(0))^2>"

f(x) = 4*D*x
D=0.1

fit f(x) filein u 1:2:3 via D

plot filein u 1:2 t '$1', f(x) t sprintf("mean: D = %.4g, error = %.4g",D, sqrt(FIT_WSSR / (FIT_NDF + 1 )))

set print directory.'$1/dat/summary.dat' append
print 'effective diffusion constant:	', D

fileout = directory.'$1/eps/MSDlogerror.eps'

set output fileout
set logscale xy

plot filein u 1:2:3 w errorbars t '$1', f(x) t sprintf("mean: D = %.4g, error = %.4g",D, sqrt(FIT_WSSR / (FIT_NDF + 1 )))

# ------------------------------------------------------------------------------------------ #

fileout = directory.'$1/eps/vaf.eps'

set output fileout

set title "velocity autocorrelation function"
set xlabel "lagtime [steps]"
set ylabel "<(v({/Symbol t}).v(0))^2>"
unset logscale xy

plot filein using 1:4 t '$1' with lines

# ------------------------------------------------------------------------------------------ #

filein  = directory.'$1/dat/GNF.dat'
fileout = directory.'$1/eps/GNF.eps'

set logscale xy

f(x) = m*x+b
fit f(x) filein using 3:4 via m,b
f2(x) = exp(b)*x**m

set output fileout

set title "Density fluctuations"
set xlabel "<A>"
set ylabel "rms fluctuation"

plot filein with linespoints t '$1', x with lines t 'slope = 1', sqrt(x) with lines t 'slope = 0.5', f2(x) with lines lt rgb "#ff00ff" t sprintf("slope of fit = %.4g", m)

set print directory.'$1/dat/summary.dat' append
print 'density fluctuations scale like:	', m

# ------------------------------------------------------------------------------------------ #
EOF
