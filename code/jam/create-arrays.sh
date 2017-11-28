#!/bin/bash

rho=( 1.0 )
l_s=( 0.010 0.015 0.020 0.025 0.030 0.035 0.040 0.045 0.050 0.055 )
l_n=( 0.450 0.455 0.460 0.465 0.470 0.475 0.480 0.485 0.490 0.495 0.500 0.505 0.510 0.515 0.520 0.525 0.530 0.535 0.540 0.545 0.550 0.555 0.560 0.565 0.570 )

noCells=100000
noSteps=1100000
run=1
ID="100000TriplePoint1ZoomOut"

rm -f input.txt
g++ active_jam_nbr_17.cpp -I boost_1_64_0/ -O3 -o a.out -std=c++11
for i in ${rho[@]}
do
	for j in ${l_n[@]}
	do
		for k in ${l_s[@]}
		do
                        printf "$ID run$run $noCells $noSteps $k $j $i\n" >> input.txt
			run=$((run+1))
		done
	done
done
exit

