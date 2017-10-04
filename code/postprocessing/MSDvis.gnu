gnuplot << EOF

set terminal postscript eps color enhanced "Helvetica" 300 size 150,150
set ticslevel 0
#set logscale cb
set cbrange [-8:-4]
set palette defined (-8 "red",  -7 "white", -4 "blue")
set pm3d interpolate 10,10
set yrange [0:0.31]
set xrange [0:1.1]
#set logscale z
set zrange [-10:-3]
set view map
set dgrid3d

density = "0.840 0.845 0.860 0.900"
fileout = 'MSD.eps'
filein = 'dens0.dat'
set output fileout
set key off

set multiplot layout 2,2 rowsfirst title "MSD"

set size 0.5,0.5
set title 'rho = '.word(density, 1)
set xlabel "Noise"
set ylabel "Self-propulsion"
set zlabel "D"
splot filein u 1:2:4 with pm3d

set size 0.5,0.5
set title 'rho = '.word(density, 2)
filein = 'dens1.dat'
splot filein u 1:2:4 with pm3d

set size 0.5,0.5
set title 'rho = '.word(density, 3)
filein = 'dens2.dat'
splot filein u 1:2:4 with pm3d

set size 0.5,0.5
set title 'rho = '.word(density, 4)
filein = 'dens3.dat'
splot filein u 1:2:4 with pm3d

EOF
